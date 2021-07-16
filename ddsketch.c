/*
 * ddsketch - implementation of ddsketch for PostgreSQL.
 *
 * Copyright (C) Tomas Vondra, 2021
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#include <limits.h>

#include "postgres.h"
#include "libpq/pqformat.h"
#include "utils/array.h"
#include "utils/lsyscache.h"
#include "catalog/pg_type.h"

PG_MODULE_MAGIC;

/*
 * On-disk representation of the ddsketch.
 */
typedef struct ddsketch_t {
	int32		vl_len_;		/* varlena header (do not touch directly!) */
	int32		flags;			/* reserved for future use (versioning, ...) */
	int64		count;			/* number of items added to the ddsketch */
	float4		alpha;			/* alpha used to size the buckets */
	int32		maxbuckets;		/* maximum number of buckets sketch */
	int32		nbuckets;		/* current number of buckets */
	int64		buckets[FLEXIBLE_ARRAY_MEMBER];
} ddsketch_t;

/*
 * An aggregate state, representing the sketch and some additional info
 * (requested percentiles, ...).
 *
 * XXX We only ever use one of values/percentiles, never both at the same
 * time. In the future the values may use a different data types than double
 * (e.g. numeric), so we keep both fields.
 */
typedef struct ddsketch_aggstate_t {
	/* basic sketch fields */
	int64		count;			/* number of items added to the ddsketch */
	float4		alpha;			/* alpha used to size the buckets */
	int32		maxbuckets;		/* maximum number of buckets sketch */
	int32		nbuckets;		/* current number of buckets */

	/* array of requested percentiles and values */
	int			npercentiles;	/* number of percentiles */
	int			nvalues;		/* number of values */

	/* variable-length fields at the end */
	double	   *percentiles;	/* array of percentiles (if any) */
	double	   *values;			/* array of values (if any) */
	int64	   *buckets;		/* sketch buckets */
} ddsketch_aggstate_t;

#define PG_GETARG_DDSKETCH(x)	(ddsketch_t *) PG_DETOAST_DATUM(PG_GETARG_DATUM(x))

/* prototypes */
PG_FUNCTION_INFO_V1(ddsketch_add_double_array);
PG_FUNCTION_INFO_V1(ddsketch_add_double_array_count);
PG_FUNCTION_INFO_V1(ddsketch_add_double_array_values);
PG_FUNCTION_INFO_V1(ddsketch_add_double_array_values_count);
PG_FUNCTION_INFO_V1(ddsketch_add_double);
PG_FUNCTION_INFO_V1(ddsketch_add_double_count);
PG_FUNCTION_INFO_V1(ddsketch_add_double_values);
PG_FUNCTION_INFO_V1(ddsketch_add_double_values_count);

PG_FUNCTION_INFO_V1(ddsketch_add_sketch_array);
PG_FUNCTION_INFO_V1(ddsketch_add_sketch_array_values);
PG_FUNCTION_INFO_V1(ddsketch_add_sketch);
PG_FUNCTION_INFO_V1(ddsketch_add_sketch_values);

PG_FUNCTION_INFO_V1(ddsketch_array_percentiles);
PG_FUNCTION_INFO_V1(ddsketch_array_percentiles_of);
PG_FUNCTION_INFO_V1(ddsketch_percentiles);
PG_FUNCTION_INFO_V1(ddsketch_percentiles_of);
PG_FUNCTION_INFO_V1(ddsketch_sketch);

PG_FUNCTION_INFO_V1(ddsketch_serial);
PG_FUNCTION_INFO_V1(ddsketch_deserial);
PG_FUNCTION_INFO_V1(ddsketch_combine);

PG_FUNCTION_INFO_V1(ddsketch_in);
PG_FUNCTION_INFO_V1(ddsketch_out);
PG_FUNCTION_INFO_V1(ddsketch_send);
PG_FUNCTION_INFO_V1(ddsketch_recv);

PG_FUNCTION_INFO_V1(ddsketch_count);

PG_FUNCTION_INFO_V1(ddsketch_add_double_increment);
PG_FUNCTION_INFO_V1(ddsketch_add_double_array_increment);
PG_FUNCTION_INFO_V1(ddsketch_union_double_increment);

Datum ddsketch_add_double_array(PG_FUNCTION_ARGS);
Datum ddsketch_add_double_array_count(PG_FUNCTION_ARGS);
Datum ddsketch_add_double_array_values(PG_FUNCTION_ARGS);
Datum ddsketch_add_double_array_values_count(PG_FUNCTION_ARGS);
Datum ddsketch_add_double(PG_FUNCTION_ARGS);
Datum ddsketch_add_double_count(PG_FUNCTION_ARGS);
Datum ddsketch_add_double_values(PG_FUNCTION_ARGS);
Datum ddsketch_add_double_values_count(PG_FUNCTION_ARGS);

Datum ddsketch_add_sketch_array(PG_FUNCTION_ARGS);
Datum ddsketch_add_sketch_array_values(PG_FUNCTION_ARGS);
Datum ddsketch_add_sketch(PG_FUNCTION_ARGS);
Datum ddsketch_add_sketch_values(PG_FUNCTION_ARGS);

Datum ddsketch_array_percentiles(PG_FUNCTION_ARGS);
Datum ddsketch_array_percentiles_of(PG_FUNCTION_ARGS);
Datum ddsketch_percentiles(PG_FUNCTION_ARGS);
Datum ddsketch_percentiles_of(PG_FUNCTION_ARGS);

Datum ddsketch_sketch(PG_FUNCTION_ARGS);

Datum ddsketch_serial(PG_FUNCTION_ARGS);
Datum ddsketch_deserial(PG_FUNCTION_ARGS);
Datum ddsketch_combine(PG_FUNCTION_ARGS);

Datum ddsketch_in(PG_FUNCTION_ARGS);
Datum ddsketch_out(PG_FUNCTION_ARGS);
Datum ddsketch_send(PG_FUNCTION_ARGS);
Datum ddsketch_recv(PG_FUNCTION_ARGS);

Datum ddsketch_count(PG_FUNCTION_ARGS);

Datum ddsketch_add_double_increment(PG_FUNCTION_ARGS);
Datum ddsketch_add_double_array_increment(PG_FUNCTION_ARGS);
Datum ddsketch_union_double_increment(PG_FUNCTION_ARGS);

static Datum double_to_array(FunctionCallInfo fcinfo, double * d, int len);
static double *array_to_double(FunctionCallInfo fcinfo, ArrayType *v, int * len);

/* boundaries for relative error */
#define	MIN_SKETCH_ALPHA	0.0001
#define MAX_SKETCH_ALPHA	0.1

#define	MIN_SKETCH_BUCKETS	16
#define MAX_SKETCH_BUCKETS	32768


/* basic checks on the ddsketch (proper sum of counts, ...) */
static void
AssertCheckDDSketch(ddsketch_t *sketch)
{
#ifdef USE_ASSERT_CHECKING
	int 	i;
	int64	count;

	Assert(sketch->alpha >= MIN_SKETCH_ALPHA);
	Assert(sketch->alpha <= MAX_SKETCH_ALPHA);

	Assert(sketch->maxbuckets >= MIN_SKETCH_BUCKETS);
	Assert(sketch->maxbuckets <= MAX_SKETCH_BUCKETS);

	Assert(sketch->nbuckets >= MIN_SKETCH_BUCKETS);
	Assert(sketch->nbuckets <= MAX_SKETCH_BUCKETS);

	Assert(sketch->maxbuckets >= sketch->nbuckets);

	count = 0;
	for (i = 0; i < sketch->nbuckets; i++)
	{
		Assert(sketch->buckets[i] >= 0);
		count += sketch->buckets[i];
	}

	Assert(count == sketch->count);
#endif
}

static void
AssertCheckDDSketchAggState(ddsketch_aggstate_t *state)
{
#ifdef USE_ASSERT_CHECKING
	int 	i;
	int64	count;

	Assert(state->alpha >= MIN_SKETCH_ALPHA);
	Assert(state->alpha <= MAX_SKETCH_ALPHA);

	Assert(state->maxbuckets >= MIN_SKETCH_BUCKETS);
	Assert(state->maxbuckets <= MAX_SKETCH_BUCKETS);

	Assert(state->nbuckets >= MIN_SKETCH_BUCKETS);
	Assert(state->nbuckets <= MAX_SKETCH_BUCKETS);

	Assert(state->maxbuckets >= state->nbuckets);

	Assert(state->npercentiles >= 0);
	Assert(state->nvalues >= 0);

	Assert(!((state->npercentiles > 0) && (state->nvalues > 0)));

	count = 0;
	for (i = 0; i < state->nbuckets; i++)
	{
		Assert(state->buckets[i] >= 0);
		count += state->buckets[i];
	}

	Assert(count == state->count);
#endif
}

/*
 * Estimate requested quantiles from the sketch agg state.
 */
static void
ddsketch_compute_quantiles(ddsketch_aggstate_t *state, double *result)
{
	int			i;
	double		gamma = (1 + state->alpha) / (1 - state->alpha);

	AssertCheckDDSketchAggState(state);

	for (i = 0; i < state->npercentiles; i++)
	{
		int		j;
		int		index = 0;
		int64	count = 0;
		double	goal = (state->percentiles[i] * state->count);

		/* FIXME This should probably look at the "next" non-zero bucket. */
		for (j = 0; j < state->nbuckets; j++)
		{
			if (state->buckets[j] == 0)
				continue;

			count += state->buckets[j];
			index = j;

			if (count > goal)
				break;
		}

		result[i] = 2 * pow(gamma, index) / (gamma + 1);
	}
}

/*
 * Estimate inverse of quantile given a value from the sketch agg state.
 *
 * Essentially an inverse to ddsketch_compute_quantiles.
 */
static void
ddsketch_compute_quantiles_of(ddsketch_aggstate_t *state, double *result)
{
	int		i;
	double	alpha = state->alpha;
	double	gamma = (1 + alpha) / (1 - alpha);

	AssertCheckDDSketchAggState(state);

	for (i = 0; i < state->nvalues; i++)
	{
		int		j;
		int		index;
		int64	count;
		double	value = state->values[i];

		/* FIXME offset the value by 1.0 to always produce positive value from log() */
		index = ceil(log(1.0 + value) / log(gamma));

		count = state->buckets[index] / 2;
		for (j = 0; j < index; j++)
			count += state->buckets[j];

		result[i] = count / (double) state->count;
	}
}


/* add a value to the sketch */
static void
ddsketch_add(ddsketch_aggstate_t *state, double v, int64 c)
{
	int		index;
	double	alpha = state->alpha;
	double	gamma = (1 + alpha) / (1 - alpha);

	AssertCheckDDSketchAggState(state);

	/* FIXME offset the value by 1.0 to always produce positive value from log() */
	index = ceil(log(1.0 + v) / log(gamma));

	/* FIXME maybe resize the sketch */
	if (index >= state->nbuckets)
		elog(ERROR, "index too high %d > %d", index, state->nbuckets);

	/* for a single point, the value is both sum and mean */
	state->buckets[index] += c;
	state->count += c;
}

/* allocate sketch with enough space for a requested number of buckets */
static ddsketch_t *
ddsketch_allocate(double alpha, int maxbuckets, int nbuckets)
{
	Size		len;
	ddsketch_t *sketch;
	char	   *ptr;

	len = offsetof(ddsketch_t, buckets) + nbuckets * sizeof(int64);

	/* we pre-allocate the array for all buckets */
	ptr = palloc0(len);
	SET_VARSIZE(ptr, len);

	sketch = (ddsketch_t *) ptr;

	sketch->flags = 0;
	sketch->count = 0;

	sketch->maxbuckets = maxbuckets;
	sketch->nbuckets = nbuckets;
	sketch->alpha = alpha;

	return sketch;
}

/*
 * allocate a ddsketch aggregate state, along with space for percentile(s)
 * and value(s) requested when calling the aggregate function
 */
static ddsketch_aggstate_t *
ddsketch_aggstate_allocate(int npercentiles, int nvalues, double alpha, int nbuckets)
{
	Size				len;
	ddsketch_aggstate_t *state;
	char			   *ptr;

	/* at least one of those values is 0 */
	Assert(nvalues == 0 || npercentiles == 0);

	/*
	 * We allocate a single chunk for the struct including percentiles and
	 * buckets.
	 */
	len = MAXALIGN(sizeof(ddsketch_aggstate_t)) +
		  MAXALIGN(sizeof(double) * npercentiles) +
		  MAXALIGN(sizeof(double) * nvalues);

	ptr = palloc0(len);

	state = (ddsketch_aggstate_t *) ptr;
	ptr += MAXALIGN(sizeof(ddsketch_aggstate_t));

	state->nvalues = nvalues;
	state->npercentiles = npercentiles;
	state->nbuckets = nbuckets;
	state->maxbuckets = nbuckets;
	state->alpha = alpha;

	if (npercentiles > 0)
	{
		state->percentiles = (double *) ptr;
		ptr += MAXALIGN(sizeof(double) * npercentiles);
	}

	if (nvalues > 0)
	{
		state->values = (double *) ptr;
		ptr += MAXALIGN(sizeof(double) * nvalues);
	}

	Assert(ptr == (char *) state + len);

	/* we may need to repalloc this later */
	state->buckets = palloc0(nbuckets * sizeof(int64));

	AssertCheckDDSketchAggState(state);

	return state;
}

static ddsketch_t *
ddsketch_aggstate_to_ddsketch(ddsketch_aggstate_t *state)
{
	ddsketch_t  *sketch;
	int			nbuckets;

	nbuckets = state->nbuckets;

	sketch = ddsketch_allocate(state->alpha, state->maxbuckets, nbuckets);

	sketch->count = state->count;
	sketch->nbuckets = nbuckets;
	sketch->maxbuckets = state->maxbuckets;
	sketch->alpha = state->alpha;

	memcpy(sketch->buckets, state->buckets, state->nbuckets * sizeof(int64));

	Assert(sketch->nbuckets == nbuckets);

	return sketch;
}

/* check that the requested percentiles are valid */
static void
check_percentiles(double *percentiles, int npercentiles)
{
	int i;

	for (i = 0; i < npercentiles; i++)
	{
		if ((percentiles[i] < 0.0) || (percentiles[i] > 1.0))
			elog(ERROR, "invalid percentile value %f, should be in [0.0, 1.0]",
				 percentiles[i]);
	}
}

static void
check_sketch_parameters(double alpha, int nbuckets)
{
	if (alpha < MIN_SKETCH_ALPHA || alpha > MAX_SKETCH_ALPHA)
		elog(ERROR, "invalid alpha value %f", alpha);

	if (nbuckets < MIN_SKETCH_BUCKETS || nbuckets > MAX_SKETCH_BUCKETS)
		elog(ERROR, "invalid number of buckets %d", nbuckets);
}

/*
 * Add a value to the sketch (create one if needed). Transition function
 * for ddsketch aggregate with a single percentile.
 */
Datum
ddsketch_add_double(PG_FUNCTION_ARGS)
{
	ddsketch_aggstate_t *state;

	MemoryContext aggcontext;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_add_double called in non-aggregate context");

	/*
	 * We want to skip NULL values altogether - we return either the existing
	 * ddsketch (if it already exists) or NULL.
	 */
	if (PG_ARGISNULL(1))
	{
		if (PG_ARGISNULL(0))
			PG_RETURN_NULL();

		/* if there already is a state accumulated, don't forget it */
		PG_RETURN_DATUM(PG_GETARG_DATUM(0));
	}

	/* if there's no ddsketch aggstate allocated, create it now */
	if (PG_ARGISNULL(0))
	{
		double	alpha = PG_GETARG_FLOAT8(2);
		int32	nbuckets = PG_GETARG_INT32(3);

		double *percentiles = NULL;
		int		npercentiles = 0;
		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, nbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		if (PG_NARGS() >= 5)
		{
			percentiles = (double *) palloc(sizeof(double));
			percentiles[0] = PG_GETARG_FLOAT8(4);
			npercentiles = 1;

			check_percentiles(percentiles, npercentiles);
		}

		state = ddsketch_aggstate_allocate(npercentiles, 0, alpha, nbuckets);

		if (percentiles)
		{
			memcpy(state->percentiles, percentiles, sizeof(double) * npercentiles);
			pfree(percentiles);
		}

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	ddsketch_add(state, PG_GETARG_FLOAT8(1), 1);

	PG_RETURN_POINTER(state);
}

/*
 * Add a value with count to the ddsketch (create one if needed). Transition
 * function for ddsketch aggregate with a single percentile.
 */
Datum
ddsketch_add_double_count(PG_FUNCTION_ARGS)
{
	int64				count;
	ddsketch_aggstate_t *state;
	MemoryContext		aggcontext;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_add_double_count called in non-aggregate context");

	/*
	 * We want to skip NULL values altogether - we return either the existing
	 * ddsketch (if it already exists) or NULL.
	 */
	if (PG_ARGISNULL(1))
	{
		if (PG_ARGISNULL(0))
			PG_RETURN_NULL();

		/* if there already is a state accumulated, don't forget it */
		PG_RETURN_DATUM(PG_GETARG_DATUM(0));
	}

	/* if there's no ddsketch aggstate allocated, create it now */
	if (PG_ARGISNULL(0))
	{
		double	alpha = PG_GETARG_FLOAT8(3);
		int32	nbuckets = PG_GETARG_INT32(4);

		double *percentiles = NULL;
		int		npercentiles = 0;
		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, nbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		if (PG_NARGS() >= 6)
		{
			percentiles = (double *) palloc(sizeof(double));
			percentiles[0] = PG_GETARG_FLOAT8(5);
			npercentiles = 1;

			check_percentiles(percentiles, npercentiles);
		}

		state = ddsketch_aggstate_allocate(npercentiles, 0, alpha, nbuckets);

		if (percentiles)
		{
			memcpy(state->percentiles, percentiles, sizeof(double) * npercentiles);
			pfree(percentiles);
		}

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	if (PG_ARGISNULL(2))
		count = 1;
	else
		count = PG_GETARG_INT64(2);

	Assert(count > 0);

	ddsketch_add(state, PG_GETARG_FLOAT8(1), count);

	PG_RETURN_POINTER(state);
}

/*
 * Add a value to the ddsketch (create one if needed). Transition function
 * for ddsketch aggregate with a single value.
 */
Datum
ddsketch_add_double_values(PG_FUNCTION_ARGS)
{
	ddsketch_aggstate_t *state;

	MemoryContext aggcontext;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_add_double called in non-aggregate context");

	/*
	 * We want to skip NULL values altogether - we return either the existing
	 * ddsketch (if it already exists) or NULL.
	 */
	if (PG_ARGISNULL(1))
	{
		if (PG_ARGISNULL(0))
			PG_RETURN_NULL();

		/* if there already is a state accumulated, don't forget it */
		PG_RETURN_DATUM(PG_GETARG_DATUM(0));
	}

	/* if there's no ddsketch aggstate allocated, create it now */
	if (PG_ARGISNULL(0))
	{
		double	alpha = PG_GETARG_FLOAT8(2);
		int32	nbuckets = PG_GETARG_INT32(3);

		double *values = NULL;
		int		nvalues = 0;
		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, nbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		if (PG_NARGS() >= 5)
		{
			values = (double *) palloc(sizeof(double));
			values[0] = PG_GETARG_FLOAT8(4);
			nvalues = 1;
		}

		state = ddsketch_aggstate_allocate(0, nvalues, alpha, nbuckets);

		if (values)
		{
			memcpy(state->values, values, sizeof(double) * nvalues);
			pfree(values);
		}

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	ddsketch_add(state, PG_GETARG_FLOAT8(1), 1);

	PG_RETURN_POINTER(state);
}

/*
 * Add a value to the ddsketch (create one if needed). Transition function
 * for ddsketch aggregate with a single value.
 */
Datum
ddsketch_add_double_values_count(PG_FUNCTION_ARGS)
{
	int64				count;
	ddsketch_aggstate_t *state;

	MemoryContext aggcontext;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_add_double called in non-aggregate context");

	/*
	 * We want to skip NULL values altogether - we return either the existing
	 * ddsketch (if it already exists) or NULL.
	 */
	if (PG_ARGISNULL(1))
	{
		if (PG_ARGISNULL(0))
			PG_RETURN_NULL();

		/* if there already is a state accumulated, don't forget it */
		PG_RETURN_DATUM(PG_GETARG_DATUM(0));
	}

	/* if there's no ddsketch aggstate allocated, create it now */
	if (PG_ARGISNULL(0))
	{
		double	alpha = PG_GETARG_FLOAT8(2);
		int32	nbuckets = PG_GETARG_INT32(3);

		double *values = NULL;
		int		nvalues = 0;
		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, nbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		if (PG_NARGS() >= 5)
		{
			values = (double *) palloc(sizeof(double));
			values[0] = PG_GETARG_FLOAT8(4);
			nvalues = 1;
		}

		state = ddsketch_aggstate_allocate(0, nvalues, alpha, nbuckets);

		if (values)
		{
			memcpy(state->values, values, sizeof(double) * nvalues);
			pfree(values);
		}

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	if (PG_ARGISNULL(2))
		count = 1;
	else
		count = PG_GETARG_INT64(2);

	Assert(count > 0);

	ddsketch_add(state, PG_GETARG_FLOAT8(1), count);

	PG_RETURN_POINTER(state);
}

/*
 * Add a value to the ddsketch (create one if needed). Transition function
 * for ddsketch aggregate with a single percentile.
 */
Datum
ddsketch_add_sketch(PG_FUNCTION_ARGS)
{
	int					i;
	ddsketch_aggstate_t *state;
	ddsketch_t		   *sketch;

	MemoryContext aggcontext;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_add_sketch called in non-aggregate context");

	/*
	 * We want to skip NULL values altogether - we return either the existing
	 * ddsketch (if it already exists) or NULL.
	 */
	if (PG_ARGISNULL(1))
	{
		if (PG_ARGISNULL(0))
			PG_RETURN_NULL();

		/* if there already is a state accumulated, don't forget it */
		PG_RETURN_DATUM(PG_GETARG_DATUM(0));
	}

	sketch = (ddsketch_t *) PG_DETOAST_DATUM(PG_GETARG_DATUM(1));

	/* if there's no aggregate state allocated, create it now */
	if (PG_ARGISNULL(0))
	{
		double *percentiles = NULL;
		int		npercentiles = 0;

		MemoryContext	oldcontext;

		oldcontext = MemoryContextSwitchTo(aggcontext);

		if (PG_NARGS() >= 3)
		{
			percentiles = (double *) palloc(sizeof(double));
			percentiles[0] = PG_GETARG_FLOAT8(2);
			npercentiles = 1;

			check_percentiles(percentiles, npercentiles);
		}

		state = ddsketch_aggstate_allocate(npercentiles, 0, sketch->alpha, sketch->maxbuckets);

		if (percentiles)
		{
			memcpy(state->percentiles, percentiles, sizeof(double) * npercentiles);
			pfree(percentiles);
		}

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	/* check that the sketch and aggstate are compatible */

	if (state->alpha != sketch->alpha)
		elog(ERROR, "state and sketch are not compatible: alpha %lf != %lf",
			 state->alpha, sketch->alpha);

	if (state->maxbuckets != sketch->maxbuckets)
		elog(ERROR, "state and sketch are not compatible: nbuckets %d != %d",
			 state->nbuckets, sketch->nbuckets);

	/* copy data from the sketch into the aggstate */
	for (i = 0; i < sketch->nbuckets; i++)
		state->buckets[i] += sketch->buckets[i];

	state->count += sketch->count;

	PG_RETURN_POINTER(state);
}

/*
 * Add a value to the ddsketch (create one if needed). Transition function
 * for ddsketch aggregate with a single value.
 */
Datum
ddsketch_add_sketch_values(PG_FUNCTION_ARGS)
{
	int					i;
	ddsketch_aggstate_t *state;
	ddsketch_t		   *sketch;

	MemoryContext aggcontext;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_add_sketch_values called in non-aggregate context");

	/*
	 * We want to skip NULL values altogether - we return either the existing
	 * ddsketch (if it already exists) or NULL.
	 */
	if (PG_ARGISNULL(1))
	{
		if (PG_ARGISNULL(0))
			PG_RETURN_NULL();

		/* if there already is a state accumulated, don't forget it */
		PG_RETURN_DATUM(PG_GETARG_DATUM(0));
	}

	sketch = (ddsketch_t *) PG_DETOAST_DATUM(PG_GETARG_DATUM(1));

	/* if there's no aggregate state allocated, create it now */
	if (PG_ARGISNULL(0))
	{
		double *values = NULL;
		int		nvalues = 0;

		MemoryContext	oldcontext;

		oldcontext = MemoryContextSwitchTo(aggcontext);

		if (PG_NARGS() >= 3)
		{
			values = (double *) palloc(sizeof(double));
			values[0] = PG_GETARG_FLOAT8(2);
			nvalues = 1;
		}

		state = ddsketch_aggstate_allocate(0, nvalues, sketch->alpha, sketch->maxbuckets);

		if (values)
		{
			memcpy(state->values, values, sizeof(double) * nvalues);
			pfree(values);
		}

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	/* copy data from sketch to aggstate */
	for (i = 0; i < state->nbuckets; i++)
		state->buckets[i] += sketch->buckets[i];

	state->count += sketch->count;

	PG_RETURN_POINTER(state);
}

/*
 * Add a value to the ddsketch (create one if needed). Transition function
 * for ddsketch aggregate with an array of percentiles.
 */
Datum
ddsketch_add_double_array(PG_FUNCTION_ARGS)
{
	ddsketch_aggstate_t *state;

	MemoryContext aggcontext;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_add_double_array called in non-aggregate context");

	/*
	 * We want to skip NULL values altogether - we return either the existing
	 * ddsketch or NULL.
	 */
	if (PG_ARGISNULL(1))
	{
		if (PG_ARGISNULL(0))
			PG_RETURN_NULL();

		/* if there already is a state accumulated, don't forget it */
		PG_RETURN_DATUM(PG_GETARG_DATUM(0));
	}

	/* if there's no ddsketch aggstate allocated, create it now */
	if (PG_ARGISNULL(0))
	{
		double	alpha = PG_GETARG_FLOAT8(2);
		int32	nbuckets = PG_GETARG_INT32(3);

		double *percentiles;
		int		npercentiles;
		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, nbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		percentiles = array_to_double(fcinfo,
									  PG_GETARG_ARRAYTYPE_P(4),
									  &npercentiles);

		check_percentiles(percentiles, npercentiles);

		state = ddsketch_aggstate_allocate(npercentiles, 0, alpha, nbuckets);

		memcpy(state->percentiles, percentiles, sizeof(double) * npercentiles);

		pfree(percentiles);

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	ddsketch_add(state, PG_GETARG_FLOAT8(1), 1);

	PG_RETURN_POINTER(state);
}

/*
 * Add a value to the ddsketch (create one if needed). Transition function
 * for ddsketch aggregate with an array of percentiles.
 */
Datum
ddsketch_add_double_array_count(PG_FUNCTION_ARGS)
{
	int64				count;
	ddsketch_aggstate_t *state;

	MemoryContext aggcontext;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_add_double_array called in non-aggregate context");

	/*
	 * We want to skip NULL values altogether - we return either the existing
	 * ddsketch or NULL.
	 */
	if (PG_ARGISNULL(1))
	{
		if (PG_ARGISNULL(0))
			PG_RETURN_NULL();

		/* if there already is a state accumulated, don't forget it */
		PG_RETURN_DATUM(PG_GETARG_DATUM(0));
	}

	/* if there's no ddsketch aggstate allocated, create it now */
	if (PG_ARGISNULL(0))
	{
		double	alpha = PG_GETARG_FLOAT8(3);
		int32	nbuckets = PG_GETARG_INT32(4);

		double *percentiles;
		int		npercentiles;
		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, nbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		percentiles = array_to_double(fcinfo,
									  PG_GETARG_ARRAYTYPE_P(5),
									  &npercentiles);

		check_percentiles(percentiles, npercentiles);

		state = ddsketch_aggstate_allocate(npercentiles, 0, alpha, nbuckets);

		memcpy(state->percentiles, percentiles, sizeof(double) * npercentiles);

		pfree(percentiles);

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	if (PG_ARGISNULL(2))
		count = 1;
	else
		count = PG_GETARG_INT64(2);

	ddsketch_add(state, PG_GETARG_FLOAT8(1), count);

	PG_RETURN_POINTER(state);
}

/*
 * Add a value to the ddsketch (create one if needed). Transition function
 * for ddsketch aggregate with an array of values.
 */
Datum
ddsketch_add_double_array_values(PG_FUNCTION_ARGS)
{
	ddsketch_aggstate_t *state;

	MemoryContext aggcontext;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_add_double_array called in non-aggregate context");

	/*
	 * We want to skip NULL values altogether - we return either the existing
	 * ddsketch or NULL.
	 */
	if (PG_ARGISNULL(1))
	{
		if (PG_ARGISNULL(0))
			PG_RETURN_NULL();

		/* if there already is a state accumulated, don't forget it */
		PG_RETURN_DATUM(PG_GETARG_DATUM(0));
	}

	/* if there's no ddsketch aggstate allocated, create it now */
	if (PG_ARGISNULL(0))
	{
		double	alpha = PG_GETARG_FLOAT8(2);
		int32	nbuckets = PG_GETARG_INT32(3);

		double *values;
		int		nvalues;
		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, nbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		values = array_to_double(fcinfo,
								 PG_GETARG_ARRAYTYPE_P(4),
								 &nvalues);

		state = ddsketch_aggstate_allocate(0, nvalues, alpha, nbuckets);

		memcpy(state->values, values, sizeof(double) * nvalues);

		pfree(values);

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	ddsketch_add(state, PG_GETARG_FLOAT8(1), 1);

	PG_RETURN_POINTER(state);
}

/*
 * Add a value to the ddsketch (create one if needed). Transition function
 * for ddsketch aggregate with an array of values.
 */
Datum
ddsketch_add_double_array_values_count(PG_FUNCTION_ARGS)
{
	int64				count;
	ddsketch_aggstate_t *state;

	MemoryContext aggcontext;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_add_double_array called in non-aggregate context");

	/*
	 * We want to skip NULL values altogether - we return either the existing
	 * ddsketch or NULL.
	 */
	if (PG_ARGISNULL(1))
	{
		if (PG_ARGISNULL(0))
			PG_RETURN_NULL();

		/* if there already is a state accumulated, don't forget it */
		PG_RETURN_DATUM(PG_GETARG_DATUM(0));
	}

	/* if there's no ddsketch aggstate allocated, create it now */
	if (PG_ARGISNULL(0))
	{
		double	alpha = PG_GETARG_FLOAT8(3);
		int32	nbuckets = PG_GETARG_INT32(4);

		double *values;
		int		nvalues;
		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, nbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		values = array_to_double(fcinfo,
								 PG_GETARG_ARRAYTYPE_P(5),
								 &nvalues);

		state = ddsketch_aggstate_allocate(0, nvalues, alpha, nbuckets);

		memcpy(state->values, values, sizeof(double) * nvalues);

		pfree(values);

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	if (PG_ARGISNULL(2))
		count = 1;
	else
		count = PG_GETARG_INT64(2);

	Assert(count > 0);

	ddsketch_add(state, PG_GETARG_FLOAT8(1), count);

	PG_RETURN_POINTER(state);
}

/*
 * Add a ddsketch to the ddsketch (create one if needed). Transition function
 * for ddsketch aggregate with an array of percentiles.
 */
Datum
ddsketch_add_sketch_array(PG_FUNCTION_ARGS)
{
	int					i;
	ddsketch_aggstate_t *state;
	ddsketch_t		   *sketch;

	MemoryContext aggcontext;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_add_sketch_array called in non-aggregate context");

	/*
	 * We want to skip NULL values altogether - we return either the existing
	 * ddsketch (if it already exists) or NULL.
	 */
	if (PG_ARGISNULL(1))
	{
		if (PG_ARGISNULL(0))
			PG_RETURN_NULL();

		/* if there already is a state accumulated, don't forget it */
		PG_RETURN_DATUM(PG_GETARG_DATUM(0));
	}

	sketch = (ddsketch_t *) PG_DETOAST_DATUM(PG_GETARG_DATUM(1));

	/* if there's no aggregate state allocated, create it now */
	if (PG_ARGISNULL(0))
	{
		double *percentiles;
		int		npercentiles;
		MemoryContext	oldcontext;

		oldcontext = MemoryContextSwitchTo(aggcontext);

		percentiles = array_to_double(fcinfo,
									  PG_GETARG_ARRAYTYPE_P(2),
									  &npercentiles);

		check_percentiles(percentiles, npercentiles);

		state = ddsketch_aggstate_allocate(npercentiles, 0, sketch->alpha, sketch->maxbuckets);

		memcpy(state->percentiles, percentiles, sizeof(double) * npercentiles);

		pfree(percentiles);

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	/* copy data from sketch to aggstate */
	for (i = 0; i < sketch->nbuckets; i++)
		state->buckets[i] += sketch->buckets[i];

	state->count += sketch->count;

	PG_RETURN_POINTER(state);
}

/*
 * Add a sketch to the ddsketch (create one if needed). Transition function
 * for ddsketch aggregate with an array of values.
 */
Datum
ddsketch_add_sketch_array_values(PG_FUNCTION_ARGS)
{
	int					i;
	ddsketch_aggstate_t *state;
	ddsketch_t		   *sketch;

	MemoryContext aggcontext;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_add_sketch_array called in non-aggregate context");

	/*
	 * We want to skip NULL values altogether - we return either the existing
	 * ddsketch (if it already exists) or NULL.
	 */
	if (PG_ARGISNULL(1))
	{
		if (PG_ARGISNULL(0))
			PG_RETURN_NULL();

		/* if there already is a state accumulated, don't forget it */
		PG_RETURN_DATUM(PG_GETARG_DATUM(0));
	}

	sketch = (ddsketch_t *) PG_DETOAST_DATUM(PG_GETARG_DATUM(1));

	/* if there's no aggregate state allocated, create it now */
	if (PG_ARGISNULL(0))
	{
		double *values;
		int		nvalues;
		MemoryContext	oldcontext;

		oldcontext = MemoryContextSwitchTo(aggcontext);

		values = array_to_double(fcinfo,
								 PG_GETARG_ARRAYTYPE_P(2),
								 &nvalues);

		state = ddsketch_aggstate_allocate(0, nvalues, sketch->alpha, sketch->maxbuckets);

		memcpy(state->values, values, sizeof(double) * nvalues);

		pfree(values);

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	/* copy data from sketch to aggstate */
	for (i = 0; i < sketch->nbuckets; i++)
		state->buckets[i] += sketch->buckets[i];

	state->count += sketch->count;

	PG_RETURN_POINTER(state);
}

/*
 * Compute percentile from a ddsketch. Final function for ddsketch aggregate
 * with a single percentile.
 */
Datum
ddsketch_percentiles(PG_FUNCTION_ARGS)
{
	ddsketch_aggstate_t	   *state;
	MemoryContext	aggcontext;
	double			ret;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_percentiles called in non-aggregate context");

	/* if there's no ddsketch, return NULL */
	if (PG_ARGISNULL(0))
		PG_RETURN_NULL();

	state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	ddsketch_compute_quantiles(state, &ret);

	PG_RETURN_FLOAT8(ret);
}

/*
 * Compute percentile from a ddsketch. Final function for ddsketch aggregate
 * with a single percentile.
 */
Datum
ddsketch_percentiles_of(PG_FUNCTION_ARGS)
{
	ddsketch_aggstate_t	   *state;
	MemoryContext	aggcontext;
	double			ret;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_percentiles_of called in non-aggregate context");

	/* if there's no ddsketch, return NULL */
	if (PG_ARGISNULL(0))
		PG_RETURN_NULL();

	state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	ddsketch_compute_quantiles_of(state, &ret);

	PG_RETURN_FLOAT8(ret);
}

/*
 * Build a ddsketch varlena value from the aggegate state.
 */
Datum
ddsketch_sketch(PG_FUNCTION_ARGS)
{
	ddsketch_t			   *sketch;
	ddsketch_aggstate_t	   *state;
	MemoryContext	aggcontext;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_sketch called in non-aggregate context");

	/* if there's no ddsketch, return NULL */
	if (PG_ARGISNULL(0))
		PG_RETURN_NULL();

	state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	sketch = ddsketch_aggstate_to_ddsketch(state);

	PG_RETURN_POINTER(sketch);
}

/*
 * Compute percentiles from a ddsketch. Final function for ddsketch aggregate
 * with an array of percentiles.
 */
Datum
ddsketch_array_percentiles(PG_FUNCTION_ARGS)
{
	double	*result;
	MemoryContext aggcontext;

	ddsketch_aggstate_t *state;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_array_percentiles called in non-aggregate context");

	if (PG_ARGISNULL(0))
		PG_RETURN_NULL();

	state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	result = palloc(state->npercentiles * sizeof(double));

	ddsketch_compute_quantiles(state, result);

	return double_to_array(fcinfo, result, state->npercentiles);
}

/*
 * Compute percentiles from a ddsketch. Final function for ddsketch aggregate
 * with an array of values.
 */
Datum
ddsketch_array_percentiles_of(PG_FUNCTION_ARGS)
{
	double	*result;
	MemoryContext aggcontext;

	ddsketch_aggstate_t *state;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_array_percentiles_of called in non-aggregate context");

	if (PG_ARGISNULL(0))
		PG_RETURN_NULL();

	state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	result = palloc(state->nvalues * sizeof(double));

	ddsketch_compute_quantiles_of(state, result);

	return double_to_array(fcinfo, result, state->nvalues);
}

Datum
ddsketch_serial(PG_FUNCTION_ARGS)
{
	bytea	   *v;
	ddsketch_aggstate_t  *state;
	Size		len;
	char	   *ptr;

	state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	len = offsetof(ddsketch_aggstate_t, percentiles) +
		  state->npercentiles * sizeof(double) +
		  state->nvalues * sizeof(double) +
		  state->nbuckets * sizeof(int64);

	v = palloc(len + VARHDRSZ);

	SET_VARSIZE(v, len + VARHDRSZ);
	ptr = VARDATA(v);

	memcpy(ptr, state, offsetof(ddsketch_aggstate_t, percentiles));
	ptr += offsetof(ddsketch_aggstate_t, percentiles);

	if (state->npercentiles > 0)
	{
		memcpy(ptr, state->percentiles, sizeof(double) * state->npercentiles);
		ptr += sizeof(double) * state->npercentiles;
	}

	if (state->nvalues > 0)
	{
		memcpy(ptr, state->values, sizeof(double) * state->nvalues);
		ptr += sizeof(double) * state->nvalues;
	}

	/* FIXME maybe don't serialize full buckets, but just the count */
	memcpy(ptr, state->buckets,
		   sizeof(int64) * state->nbuckets);
	ptr += sizeof(int64) * state->nbuckets;

	Assert(VARDATA(v) + len == ptr);

	PG_RETURN_POINTER(v);
}

Datum
ddsketch_deserial(PG_FUNCTION_ARGS)
{
	bytea  *v = (bytea *) PG_GETARG_POINTER(0);
	char   *ptr = VARDATA_ANY(v);
	ddsketch_aggstate_t	tmp;
	ddsketch_aggstate_t *state;
	double			   *percentiles = NULL;
	double			   *values = NULL;

	/* copy aggstate header into a local variable */
	memcpy(&tmp, ptr, offsetof(ddsketch_aggstate_t, percentiles));
	ptr += offsetof(ddsketch_aggstate_t, percentiles);

	/* allocate and copy percentiles */
	if (tmp.npercentiles > 0)
	{
		percentiles = palloc(tmp.npercentiles * sizeof(double));
		memcpy(percentiles, ptr, tmp.npercentiles * sizeof(double));
		ptr += tmp.npercentiles * sizeof(double);
	}

	/* allocate and copy values */
	if (tmp.nvalues > 0)
	{
		values = palloc(tmp.nvalues * sizeof(double));
		memcpy(values, ptr, tmp.nvalues * sizeof(double));
		ptr += tmp.nvalues * sizeof(double);
	}

	state = ddsketch_aggstate_allocate(tmp.npercentiles, tmp.nvalues,
									   tmp.alpha, tmp.nbuckets);

	if (tmp.npercentiles > 0)
	{
		memcpy(state->percentiles, percentiles, tmp.npercentiles * sizeof(double));
		pfree(percentiles);
	}

	if (tmp.nvalues > 0)
	{
		memcpy(state->values, values, tmp.nvalues * sizeof(double));
		pfree(values);
	}

	/* copy the data into the newly-allocated state */
	memcpy(state, &tmp, offsetof(ddsketch_aggstate_t, percentiles));
	/* we don't need to move the pointer */

	/* copy the buckets back */
	memcpy(state->buckets, ptr,
		   sizeof(int64) * state->nbuckets);
	ptr += sizeof(int64) * state->nbuckets;

	PG_RETURN_POINTER(state);
}

static ddsketch_aggstate_t *
ddsketch_copy(ddsketch_aggstate_t *state)
{
	ddsketch_aggstate_t *copy;

	copy = ddsketch_aggstate_allocate(state->npercentiles, state->nvalues,
									  state->alpha, state->nbuckets);

	memcpy(copy, state, offsetof(ddsketch_aggstate_t, percentiles));

	if (state->nvalues > 0)
		memcpy(copy->values, state->values,
			   sizeof(double) * state->nvalues);

	if (state->npercentiles > 0)
		memcpy(copy->percentiles, state->percentiles,
			   sizeof(double) * state->npercentiles);

	memcpy(copy->buckets, state->buckets, state->nbuckets * sizeof(int64));

	return copy;
}

Datum
ddsketch_combine(PG_FUNCTION_ARGS)
{
	int					i;
	ddsketch_aggstate_t	 *src;
	ddsketch_aggstate_t	 *dst;

	MemoryContext aggcontext;
	MemoryContext oldcontext;

	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_combine called in non-aggregate context");

	/* the second parameter must not be NULL */
	Assert(!PG_ARGISNULL(1));

	/* so just grab it */
	src = (ddsketch_aggstate_t *) PG_GETARG_POINTER(1);

	/* when NULL in the first parameter, just return a copy of the second one */
	if (PG_ARGISNULL(0))
	{
		/* copy the ddsketch into the right long-lived memory context */
		oldcontext = MemoryContextSwitchTo(aggcontext);
		src = ddsketch_copy(src);
		MemoryContextSwitchTo(oldcontext);

		PG_RETURN_POINTER(src);
	}

	dst = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	AssertCheckDDSketchAggState(dst);
	AssertCheckDDSketchAggState(src);

	for (i = 0; i < dst->nbuckets; i++)
		dst->buckets[i] += src->buckets[i];

	dst->count += src->count;

	AssertCheckDDSketchAggState(dst);

	PG_RETURN_POINTER(dst);
}

/* API for incremental updates */

/*
 * expand the ddsketch into an in-memory aggregate state
 */
static ddsketch_aggstate_t *
ddsketch_sketch_to_aggstate(ddsketch_t *sketch)
{
	ddsketch_aggstate_t *state;

	state = ddsketch_aggstate_allocate(0, 0, sketch->alpha, sketch->maxbuckets);

	state->count = sketch->count;

	/* copy data from the ddsketch into the aggstate */
	memcpy(state->buckets, sketch->buckets, sketch->nbuckets * sizeof(int64));

	return state;
}

/*
 * Add a single value to the ddsketch. This is not very efficient, as it has
 * to deserialize the ddsketch into the in-memory aggstate representation
 * and serialize it back for each call, but it's convenient and acceptable
 * for some use cases.
 *
 * When efficiency is important, it may be possible to use the batch variant
 * with first aggregating the updates into a ddsketch, and then merge that
 * into an existing ddsketch in one step using ddsketch_union_double_increment
 *
 * This is similar to hll_add, while the "union" is more like hll_union.
 */
Datum
ddsketch_add_double_increment(PG_FUNCTION_ARGS)
{
	ddsketch_aggstate_t *state;

	/*
	 * We want to skip NULL values altogether - we return either the existing
	 * ddsketch (if it already exists) or NULL.
	 */
	if (PG_ARGISNULL(1))
	{
		if (PG_ARGISNULL(0))
			PG_RETURN_NULL();

		/* if there already is a state accumulated, don't forget it */
		PG_RETURN_DATUM(PG_GETARG_DATUM(0));
	}

	/* if there's no aggstate allocated, create it now */
	if (PG_ARGISNULL(0))
	{
		double	alpha;
		int32	nbuckets;

		/*
		 * We don't require compression, but only when there is an existing
		 * ddsketch value. Make sure the value was supplied.
		 */
		if (PG_ARGISNULL(2))
			elog(ERROR, "alpha value not supplied, but ddsketch is NULL");

		if (PG_ARGISNULL(3))
			elog(ERROR, "nbuckets value not supplied, but ddsketch is NULL");

		alpha = PG_GETARG_FLOAT8(2);
		nbuckets = PG_GETARG_INT32(3);

		check_sketch_parameters(alpha, nbuckets);

		state = ddsketch_aggstate_allocate(0, 0, alpha, nbuckets);
	}
	else
		state = ddsketch_sketch_to_aggstate(PG_GETARG_DDSKETCH(0));

	ddsketch_add(state, PG_GETARG_FLOAT8(1), 1);

	PG_RETURN_POINTER(ddsketch_aggstate_to_ddsketch(state));
}

/*
 * Add an array of values to the ddsketch. This amortizes the overhead of
 * deserializing and serializing the ddsketch, compared to the per-value
 * version.
 *
 * When efficiency is important, it may be possible to use the batch variant
 * with first aggregating the updates into a ddsketch, and then merge that
 * into an existing ddsketch in one step using ddsketch_union_double_increment
 *
 * This is similar to hll_add, while the "union" is more like hll_union.
 */
Datum
ddsketch_add_double_array_increment(PG_FUNCTION_ARGS)
{
	ddsketch_aggstate_t *state;
	double			   *values;
	int					nvalues;
	int					i;

	/*
	 * We want to skip NULL values altogether - we return either the existing
	 * ddsketch (if it already exists) or NULL.
	 */
	if (PG_ARGISNULL(1))
	{
		if (PG_ARGISNULL(0))
			PG_RETURN_NULL();

		/* if there already is a state accumulated, don't forget it */
		PG_RETURN_DATUM(PG_GETARG_DATUM(0));
	}

	/* if there's no ddsketch allocated, create it now */
	if (PG_ARGISNULL(0))
	{
		double	alpha;
		int		nbuckets;

		/*
		 * We don't require compression, but only when there is an existing
		 * ddsketch value. Make sure the value was supplied.
		 */
		if (PG_ARGISNULL(2))
			elog(ERROR, "alpha value not supplied, but ddsketch is NULL");

		if (PG_ARGISNULL(3))
			elog(ERROR, "nbuckets value not supplied, but ddsketch is NULL");

		alpha = PG_GETARG_FLOAT8(2);
		nbuckets = PG_GETARG_INT32(3);

		check_sketch_parameters(alpha, nbuckets);

		state = ddsketch_aggstate_allocate(0, 0, alpha, nbuckets);
	}
	else
		state = ddsketch_sketch_to_aggstate(PG_GETARG_DDSKETCH(0));

	values = array_to_double(fcinfo,
							 PG_GETARG_ARRAYTYPE_P(1),
							 &nvalues);

	for (i = 0; i < nvalues; i++)
		ddsketch_add(state, values[i], 1);

	PG_RETURN_POINTER(ddsketch_aggstate_to_ddsketch(state));
}

/*
 * Merge a ddsketch into another ddsketch. This is somewaht inefficient, as
 * it has to deserialize the sketches into the in-memory aggstate values,
 * and serialize it back for each call, but it's better than doing it for
 * each individual value (like ddsketch_union_double_increment).
 *
 * This is similar to hll_union.
 */
Datum
ddsketch_union_double_increment(PG_FUNCTION_ARGS)
{
	int					i;
	ddsketch_aggstate_t *state;
	ddsketch_t		   *sketch;

	if (PG_ARGISNULL(0) && PG_ARGISNULL(1))
		PG_RETURN_NULL();
	else if (PG_ARGISNULL(0))
		PG_RETURN_POINTER(PG_GETARG_POINTER(1));
	else if (PG_ARGISNULL(1))
		PG_RETURN_POINTER(PG_GETARG_POINTER(0));

	/* now we know both arguments are non-null */

	/* parse the first ddsketch (we'll merge the other one into this) */
	state = ddsketch_sketch_to_aggstate(PG_GETARG_DDSKETCH(0));
	AssertCheckDDSketchAggState(state);

	/* parse the second ddsketch */
	sketch = PG_GETARG_DDSKETCH(1);
	AssertCheckDDSketch(sketch);

	/* copy data from sketch to aggstate */
	for (i = 0; i < sketch->nbuckets; i++)
		state->buckets[i] += sketch->buckets[i];

	state->count += sketch->count;

	AssertCheckDDSketchAggState(state);

	PG_RETURN_POINTER(ddsketch_aggstate_to_ddsketch(state));
}


Datum
ddsketch_in(PG_FUNCTION_ARGS)
{
	int			i, r;
	char	   *str = PG_GETARG_CSTRING(0);
	ddsketch_t  *sketch = NULL;

	/* ddsketch header fields */
	int32       flags;
	int64		count;
	double		alpha;
	int			maxbuckets;
	int			nbuckets;
	int			header_length;
	char	   *ptr;

	r = sscanf(str, "flags %d count " INT64_FORMAT " alpha %lf maxbuckets %d buckets %d%n",
			   &flags, &count, &alpha, &maxbuckets, &nbuckets, &header_length);

	if (r != 5)
		elog(ERROR, "failed to parse ddsketch value");

	if ((alpha < MIN_SKETCH_ALPHA) || (alpha > MAX_SKETCH_ALPHA))
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("alpha for ddsketch (%f) must be in [%f, %f]",
						alpha, MIN_SKETCH_ALPHA, MAX_SKETCH_ALPHA)));

	if ((maxbuckets < MIN_SKETCH_BUCKETS) || (maxbuckets > MAX_SKETCH_BUCKETS))
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("number of buckets (%d) for ddsketch must be in [%d, %d]",
						nbuckets, MIN_SKETCH_BUCKETS, MAX_SKETCH_BUCKETS)));

	if ((nbuckets < MIN_SKETCH_BUCKETS) || (nbuckets > MAX_SKETCH_BUCKETS))
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("number of buckets (%d) for ddsketch must be in [%d, %d]",
						nbuckets, MIN_SKETCH_BUCKETS, MAX_SKETCH_BUCKETS)));

	if (nbuckets > maxbuckets)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("number of buckets (%d) for ddsketch must not exceed maxbuckets (%d)",
						nbuckets, maxbuckets)));

	if (count <= 0)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("count value for the ddsketch must be positive")));

	sketch = ddsketch_allocate(alpha, maxbuckets, nbuckets);

	sketch->flags = flags;
	sketch->count = count;
	sketch->maxbuckets = maxbuckets;
	sketch->nbuckets = nbuckets;
	sketch->alpha = alpha;
	sketch->count = 0;

	ptr = str + header_length;

	i = 0;
	while (true)
	{
		int		index;

		if (sscanf(ptr, " (%d, " INT64_FORMAT ")", &index, &count) != 2)
			break;

		/* check that index is between 0 and (nbuckets-1) */
		if (index < 0)
			ereport(ERROR,
					(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
					 errmsg("bucket index (%d) in a ddsketch must not be negative", index)));
		else if (index >= sketch->nbuckets)
			ereport(ERROR,
					(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
					 errmsg("bucket index (%d) in a ddsketch overflows number of buckets (%d)",
							index, sketch->nbuckets)));

		/* we don't include empty buckets */
		if (count <= 0)
			ereport(ERROR,
					(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
					 errmsg("count value for all indexes in a ddsketch must be positive")));

		sketch->buckets[index] = count;
		sketch->count += count;

		/* skip to the end of the centroid */
		ptr = strchr(ptr, ')') + 1;
		i++;
	}

	Assert(ptr == str + strlen(str));

	AssertCheckDDSketch(sketch);

	PG_RETURN_POINTER(sketch);
}

Datum
ddsketch_out(PG_FUNCTION_ARGS)
{
	int			i;
	ddsketch_t  *sketch = (ddsketch_t *) PG_DETOAST_DATUM(PG_GETARG_DATUM(0));
	StringInfoData	str;

	AssertCheckDDSketch(sketch);

	initStringInfo(&str);

	appendStringInfo(&str, "flags %d count " INT64_FORMAT " alpha %lf maxbuckets %d buckets %d",
					 sketch->flags, sketch->count, sketch->alpha, sketch->maxbuckets, sketch->nbuckets);

	for (i = 0; i < sketch->nbuckets; i++)
	{
		if (sketch->buckets[i] == 0)
			continue;

		appendStringInfo(&str, " (%d, " INT64_FORMAT ")", i, sketch->buckets[i]);
	}

	PG_RETURN_CSTRING(str.data);
}

Datum
ddsketch_recv(PG_FUNCTION_ARGS)
{
	StringInfo	buf = (StringInfo) PG_GETARG_POINTER(0);
	ddsketch_t  *sketch;
	int			i;
	int64		count;
	int32		flags;
	int32		maxbuckets;
	int32		nbuckets;
	double		alpha;

	flags = pq_getmsgint(buf, sizeof(int32));

	count = pq_getmsgint64(buf);
	alpha = pq_getmsgfloat8(buf);
	maxbuckets = pq_getmsgint(buf, sizeof(int32));
	nbuckets = pq_getmsgint(buf, sizeof(int32));

	sketch = ddsketch_allocate(alpha, maxbuckets, nbuckets);

	sketch->flags = flags;
	sketch->count = count;
	sketch->nbuckets = nbuckets;
	sketch->maxbuckets = maxbuckets;

	for (i = 0; i < sketch->nbuckets; i++)
		sketch->buckets[i] = pq_getmsgint64(buf);

	PG_RETURN_POINTER(sketch);
}

Datum
ddsketch_send(PG_FUNCTION_ARGS)
{
	ddsketch_t  *sketch = (ddsketch_t *) PG_DETOAST_DATUM(PG_GETARG_DATUM(0));
	StringInfoData buf;
	int			i;

	pq_begintypsend(&buf);

	pq_sendint(&buf, sketch->flags, 4);
	pq_sendint64(&buf, sketch->count);
	pq_sendfloat8(&buf, sketch->alpha);
	pq_sendint(&buf, sketch->maxbuckets, 4);
	pq_sendint(&buf, sketch->nbuckets, 4);

	for (i = 0; i < sketch->nbuckets; i++)
		pq_sendint64(&buf, sketch->buckets[i]);

	PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

Datum
ddsketch_count(PG_FUNCTION_ARGS)
{
	ddsketch_t  *sketch = (ddsketch_t *) PG_DETOAST_DATUM(PG_GETARG_DATUM(0));

	PG_RETURN_INT64(sketch->count);
}

/*
 * Transform an input FLOAT8 SQL array to a plain double C array.
 *
 * This expects a single-dimensional float8 array, fails otherwise.
 */
static double *
array_to_double(FunctionCallInfo fcinfo, ArrayType *v, int *len)
{
	double *result;
	int		nitems,
		   *dims,
			ndims;
	Oid		element_type;
	int16	typlen;
	bool	typbyval;
	char	typalign;
	int		i;

	/* deconstruct_array */
	Datum	   *elements;
	bool	   *nulls;
	int			nelements;

	ndims = ARR_NDIM(v);
	dims = ARR_DIMS(v);
	nitems = ArrayGetNItems(ndims, dims);

	/* this is a special-purpose function for single-dimensional arrays */
	if (ndims != 1)
		elog(ERROR, "expected a single-dimensional array (dims = %d)", ndims);

	/*
	 * if there are no elements, set the length to 0 and return NULL
	 *
	 * XXX Can this actually happen? for empty arrays we seem to error out
	 * on the preceding check, i.e. ndims = 0.
	 */
	if (nitems == 0)
	{
		(*len) = 0;
		return NULL;
	}

	element_type = ARR_ELEMTYPE(v);

	/* XXX not sure if really needed (can it actually happen?) */
	if (element_type != FLOAT8OID)
		elog(ERROR, "array_to_double expects FLOAT8 array");

	/* allocate space for enough elements */
	result = (double*) palloc(nitems * sizeof(double));

	get_typlenbyvalalign(element_type, &typlen, &typbyval, &typalign);

	deconstruct_array(v, element_type, typlen, typbyval, typalign,
					  &elements, &nulls, &nelements);

	/* we should get the same counts here */
	Assert(nelements == nitems);

	for (i = 0; i < nelements; i++)
	{
		if (nulls[i])
			elog(ERROR, "NULL not allowed as a percentile value");

		result[i] = DatumGetFloat8(elements[i]);
	}

	(*len) = nelements;

	return result;
}

/*
 * construct an SQL array from a simple C double array
 */
static Datum
double_to_array(FunctionCallInfo fcinfo, double *d, int len)
{
	ArrayBuildState *astate = NULL;
	int		 i;

	for (i = 0; i < len; i++)
	{
		/* stash away this field */
		astate = accumArrayResult(astate,
								  Float8GetDatum(d[i]),
								  false,
								  FLOAT8OID,
								  CurrentMemoryContext);
	}

	PG_RETURN_ARRAYTYPE_P(makeArrayResult(astate,
										  CurrentMemoryContext));
}
