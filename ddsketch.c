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
#include <float.h>

#include "postgres.h"
#include "libpq/pqformat.h"
#include "utils/array.h"
#include "utils/lsyscache.h"
#include "catalog/pg_type.h"

PG_MODULE_MAGIC;

/*
 * On-disk representation of the ddsketch.
 *
 * XXX We store just non-empty buckets. Most sketches tend to be sparse, i.e.
 * only a small fraction of buckets is non-empty. Consider for example a sketch
 * of latency of an API call - the values are probably from a very narrow
 * interval most of the time, hence only very few buckets will be non-empty.
 * In particular, the lower buckets tend to be empty (e.g. for latencies),
 * because the latency is usually non-zero and when something breaks the
 * latency increases. So by not storing the empty buckets we probably save
 * quite a bit of space, although it depends on how many non-empty ranges
 * are there, actually. The other option would be to use varint instead of
 * int64, which would make storing buckets are almost empty more efficient.
 *
 * XXX Consider using varint when seralizing the on-disk data.
 */
typedef struct bucket_t {
	int32	index;
	int64	count;
} bucket_t;

typedef struct ddsketch_t {
	int32		vl_len_;		/* varlena header (do not touch directly!) */
	int32		flags;			/* reserved for future use (versioning, ...) */
	int64		count;			/* number of items added to the ddsketch */
	float8		alpha;			/* alpha used to size the buckets */
	int32		maxbuckets;		/* maximum number of buckets sketch */
	int32		nbuckets;		/* number of buckets used */
	bucket_t	buckets[FLEXIBLE_ARRAY_MEMBER];
} ddsketch_t;

/*
 * An aggregate state, representing the sketch and some additional info
 * (requested percentiles, ...).
 *
 * XXX We only ever use one of values/percentiles, never both at the same
 * time. In the future the values may use a different data types than double
 * (e.g. numeric), so we keep both fields.
 *
 * XXX Most of the considerations about representing sparse sketches applies
 * here too, except that the varlena compression won't help us with in-memory
 * representation. So it's a bit more pressing issue here.
 *
 * XXX Currently the code simply errors-out if the value being added would
 * require a bucket outside the maxbuckets range. But we could also combine
 * buckets on either upper or lower end of the sketch - we're probably
 * interested in percentiles on one of the tails, so we may sacrifice
 * precision on the other end (or in the middle). That's pretty much the
 * same idea as t-digest, although with a guarantee on relative error.
 */
typedef struct ddsketch_aggstate_t {
	/* basic sketch fields */
	int64		count;			/* number of items added to the ddsketch */
	float8		alpha;			/* alpha used to size the buckets */
	int32		maxbuckets;		/* maximum number of buckets sketch */
	int32		nbuckets_used;	/* number of buckets (used) */
	int32		nbuckets_allocated;	/* number of buckets (allocated) */

	/* pre-calculated parameters for mapping etc. */
	int32		offset;
	double		min_indexable_value;
	double		max_indexable_value;
	double		multiplier;
	double		gamma;

	/* array of requested percentiles and values */
	int			npercentiles;	/* number of percentiles */
	int			nvalues;		/* number of values */

	/* variable-length fields at the end */
	double	   *percentiles;	/* array of percentiles (if any) */
	double	   *values;			/* array of values (if any) */
	bucket_t   *buckets;		/* sketch buckets */
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

/* mapping to bucket indexes etc. */
static double ddsketch_log_gamma(ddsketch_aggstate_t *state, double value);
static int    ddsketch_map_index(ddsketch_aggstate_t *state, double value);

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

	Assert(sketch->nbuckets <= MAX_SKETCH_BUCKETS);

	Assert(sketch->maxbuckets >= sketch->nbuckets);

	count = 0;
	for (i = 0; i < sketch->nbuckets; i++)
	{
		Assert((i == 0) || (sketch->buckets[i-1].index < sketch->buckets[i].index));
		Assert(sketch->buckets[i].count > 0);
		count += sketch->buckets[i].count;
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

	Assert(state->nbuckets_used <= MAX_SKETCH_BUCKETS);

	Assert(state->nbuckets_allocated <= MAX_SKETCH_BUCKETS);

	Assert(state->maxbuckets >= state->nbuckets_used);
	Assert(state->maxbuckets >= state->nbuckets_allocated);

	Assert(state->npercentiles >= 0);
	Assert(state->nvalues >= 0);

	Assert(!((state->npercentiles > 0) && (state->nvalues > 0)));

	count = 0;
	for (i = 0; i < state->nbuckets_used; i++)
	{
		Assert((i == 0) || (state->buckets[i-1].index < state->buckets[i].index));
		Assert(state->buckets[i].count > 0);
		count += state->buckets[i].count;
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
		for (j = 0; j < state->nbuckets_used; j++)
		{
			count += state->buckets[j].count;
			index = state->buckets[j].index;

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

		if (value < 1.0)
			elog(ERROR, "ddsketch_add: the value has to be at least 1.0");

		index = ceil(log(value) / log(gamma));

		count = state->buckets[index].count / 2;
		for (j = 0; j < index; j++)
			count += state->buckets[j].count;

		result[i] = count / (double) state->count;
	}
}

static int
bucket_comparator(const void *a, const void *b)
{
	bucket_t *ba = (bucket_t *) a;
	bucket_t *bb = (bucket_t *) b;

	if (ba->index < bb->index)
		return -1;
	else if (ba->index > bb->index)
		return 1;

	return 0;
}

/* add a value to the sketch */
static void
ddsketch_add(ddsketch_aggstate_t *state, double v, int64 c)
{
	int		i,
			index;

	AssertCheckDDSketchAggState(state);

	if (v > state->min_indexable_value)
		index = ddsketch_map_index(state, v);
	else if (v < -state->min_indexable_value)
		index = ddsketch_map_index(state, -v);
	else
	{
		/* FIXME increment the zero bucket */
		return;
	}

	/*
	 * See if we already have a bucket with the calculated index.
	 *
	 * XXX Linear search - should be optimized to use bsearch or some kind
	 * of simple hash table.
	 */
	for (i = 0; i < state->nbuckets_used; i++)
	{
		/* If we found an existing bucket, we're done. */
		if (state->buckets[i].index == index)
		{
			state->buckets[i].count += c;
			state->count += c;
			return;
		}
	}

	/*
	 * Bucket does not exist yet. so we need to add it. If we already have
	 * enough space pre-allocated, we just add it and then sort the buckets.
	 * Otherwise allocate more space, using the usual doubling approach.
	 *
	 * XXX If we reach the maximum number of buckets, we error-out. We could
	 * also combine some of the buckets, in the less interesting part of the
	 * sketch (middle, lower buckets).
	 */
	if (state->nbuckets_used == state->nbuckets_allocated)
	{
		/* double the space for buckets, but cap by maxbuckets */
		state->nbuckets_allocated *= 2;

		state->nbuckets_allocated = Min(state->nbuckets_allocated,
										state->maxbuckets);

		/* if still equal, we've reached the maximum */
		if (state->nbuckets_used == state->nbuckets_allocated)
			elog(ERROR, "bucket overflow (used %d, allocated %d, max %d)",
				 state->nbuckets_used,
				 state->nbuckets_allocated,
				 state->maxbuckets);

		/* otherwise reallocate the space to add space */
		state->buckets = repalloc(state->buckets,
								  state->nbuckets_allocated * sizeof(bucket_t));
	}

	/* at this point there has to be space for at least one more bucket */
	Assert(state->nbuckets_used < state->nbuckets_allocated);

	state->buckets[state->nbuckets_used].index = index;
	state->buckets[state->nbuckets_used].count = c;
	state->count += c;

	state->nbuckets_used++;

	/* sort the buckets by index */
	pg_qsort(state->buckets, state->nbuckets_used, sizeof(bucket_t), bucket_comparator);
}

/*
 * ddsketch_allocate
 *		allocate sketch with enough space for a requested number of buckets
 *
 * We allocate space only for nbuckets buckets, but we also store the maximum
 * allowed number of buckets the sketch is allowed to use.
 */
static ddsketch_t *
ddsketch_allocate(int64 count, double alpha, int maxbuckets, int nbuckets)
{
	Size		len;
	ddsketch_t *sketch;
	char	   *ptr;

	Assert(nbuckets >= 0);
	Assert(nbuckets <= maxbuckets);

	len = offsetof(ddsketch_t, buckets) + nbuckets * sizeof(bucket_t);

	/* we pre-allocate the array for all buckets */
	ptr = palloc0(len);
	SET_VARSIZE(ptr, len);

	sketch = (ddsketch_t *) ptr;

	sketch->flags = 0;
	sketch->count = count;
	sketch->maxbuckets = maxbuckets;
	sketch->nbuckets = nbuckets;
	sketch->alpha = alpha;

	return sketch;
}

/*
 * ddsketch_aggstate_allocate
 *		allocate a ddsketch aggregate state, along with space for percentile(s)
 * and value(s) requested when calling the aggregate function
 */
static ddsketch_aggstate_t *
ddsketch_aggstate_allocate(int npercentiles, int nvalues, double alpha,
						   int maxbuckets, int nbuckets)
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
	state->maxbuckets = maxbuckets;
	state->nbuckets_allocated = nbuckets;
	state->nbuckets_used = 0;
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
	state->buckets = palloc0(nbuckets * sizeof(bucket_t));

	/* precalculate various parameters used for mapping */
	state->gamma = (1 + alpha) / (1 - alpha);
	state->multiplier = log(2.0) / log1p(2 * alpha / (1 - alpha));
	state->min_indexable_value = DBL_MIN * state->gamma;
	state->max_indexable_value = DBL_MAX / state->gamma;

	AssertCheckDDSketchAggState(state);

	return state;
}

static ddsketch_t *
ddsketch_aggstate_to_ddsketch(ddsketch_aggstate_t *state)
{
	ddsketch_t *sketch;

	sketch = ddsketch_allocate(state->count, state->alpha,
							   state->maxbuckets, state->nbuckets_used);

	memcpy(sketch->buckets, state->buckets, sizeof(bucket_t) * state->nbuckets_used);

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
		int32	maxbuckets = PG_GETARG_INT32(3);

		double *percentiles = NULL;
		int		npercentiles = 0;
		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, maxbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		if (PG_NARGS() >= 5)
		{
			percentiles = (double *) palloc(sizeof(double));
			percentiles[0] = PG_GETARG_FLOAT8(4);
			npercentiles = 1;

			check_percentiles(percentiles, npercentiles);
		}

		state = ddsketch_aggstate_allocate(npercentiles, 0, alpha,
										   maxbuckets, MIN_SKETCH_BUCKETS);

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
		int32	maxbuckets = PG_GETARG_INT32(4);

		double *percentiles = NULL;
		int		npercentiles = 0;
		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, maxbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		if (PG_NARGS() >= 6)
		{
			percentiles = (double *) palloc(sizeof(double));
			percentiles[0] = PG_GETARG_FLOAT8(5);
			npercentiles = 1;

			check_percentiles(percentiles, npercentiles);
		}

		state = ddsketch_aggstate_allocate(npercentiles, 0, alpha,
										   maxbuckets, MIN_SKETCH_BUCKETS);

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

	/* can't add values with non-positive counts */
	if (count <= 0)
		elog(ERROR, "invalid count value %ld, must be a positive value", count);

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
		int32	maxbuckets = PG_GETARG_INT32(3);

		double *values = NULL;
		int		nvalues = 0;
		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, maxbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		if (PG_NARGS() >= 5)
		{
			values = (double *) palloc(sizeof(double));
			values[0] = PG_GETARG_FLOAT8(4);
			nvalues = 1;
		}

		state = ddsketch_aggstate_allocate(0, nvalues, alpha,
										   maxbuckets, MIN_SKETCH_BUCKETS);

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
		int32	maxbuckets = PG_GETARG_INT32(3);

		double *values = NULL;
		int		nvalues = 0;
		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, maxbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		if (PG_NARGS() >= 5)
		{
			values = (double *) palloc(sizeof(double));
			values[0] = PG_GETARG_FLOAT8(4);
			nvalues = 1;
		}

		state = ddsketch_aggstate_allocate(0, nvalues, alpha,
										   maxbuckets, MIN_SKETCH_BUCKETS);

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

	/* can't add values with non-positive counts */
	if (count <= 0)
		elog(ERROR, "invalid count value %ld, must be a positive value", count);

	ddsketch_add(state, PG_GETARG_FLOAT8(1), count);

	PG_RETURN_POINTER(state);
}

/*
 * Add a value to the ddsketch (create one if needed). Transition function
 * for ddsketch aggregate with a single percentile.
 *
 * FIXME this is wrong, we need to add the buckets one by one, and check
 * if the bucket already exists.
 */
Datum
ddsketch_add_sketch(PG_FUNCTION_ARGS)
{
	int					i;
	ddsketch_aggstate_t *state;
	ddsketch_t		   *sketch;
	int					nbuckets_used;

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

		state = ddsketch_aggstate_allocate(npercentiles, 0, sketch->alpha,
										   sketch->maxbuckets, sketch->nbuckets);

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

	/*
	 * XXX Maybe this check is too strict and we should either ignore it or use
	 * max of the two values? Because why fail if both pieces have nbucket well
	 * within those limits?
	 */
	if (state->maxbuckets != sketch->maxbuckets)
		elog(ERROR, "state and sketch are not compatible: nbuckets %d != %d",
			 state->maxbuckets, sketch->maxbuckets);

	/*
	 * Walk through the sketch buckets, and see if there's a matching bucket
	 * in the aggregate state. If yes, just increment the counters. If not,
	 * add the bucket at the end.
	 *
	 * XXX We remember the number of used buckets, so that we don't search
	 * in the newly added ones.
	 */
	nbuckets_used = state->nbuckets_used;

	for (i = 0; i < sketch->nbuckets; i++)
	{
		int		j;
		bool	found = false;

		/* FIXME Linear search, should be improved. */
		for (j = 0; j < nbuckets_used; j++)
		{
			if (state->buckets[j].index != sketch->buckets[i].index)
				continue;

			state->buckets[j].count += sketch->buckets[i].index;
			found = true;

			break;
		}

		/* if added to an existing bucket, we're done */
		if (found)
			continue;

		/* otherwise add the bucket at the very end */
		if (state->nbuckets_used == state->nbuckets_allocated)
		{
			/* double the space for buckets, but cap by maxbuckets */
			state->nbuckets_allocated *= 2;

			state->nbuckets_allocated = Min(state->nbuckets_allocated,
											state->maxbuckets);

			/* if still equal, we've reached the maximum */
			if (state->nbuckets_used == state->nbuckets_allocated)
				elog(ERROR, "bucket overflow (used %d, allocated %d, max %d)",
					 state->nbuckets_used,
					 state->nbuckets_allocated,
					 state->maxbuckets);

			/* otherwise reallocate the space to add space */
			state->buckets = repalloc(state->buckets,
									  state->nbuckets_allocated * sizeof(bucket_t));

		}

		/* at this point there has to be space for at least one more bucket */
		Assert(state->nbuckets_used < state->nbuckets_allocated);

		state->buckets[state->nbuckets_used] = sketch->buckets[i];
		state->count += sketch->buckets[i].count;
		state->nbuckets_used++;
	}

	/* sort the buckets by index */
	pg_qsort(state->buckets, state->nbuckets_used, sizeof(bucket_t), bucket_comparator);

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
	int					nbuckets_used;
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

		state = ddsketch_aggstate_allocate(0, nvalues, sketch->alpha,
										   sketch->maxbuckets, sketch->nbuckets);

		if (values)
		{
			memcpy(state->values, values, sizeof(double) * nvalues);
			pfree(values);
		}

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	/* check that the sketch and aggstate are compatible */

	if (state->alpha != sketch->alpha)
		elog(ERROR, "state and sketch are not compatible: alpha %lf != %lf",
			 state->alpha, sketch->alpha);

	/*
	 * XXX Maybe this check is too strict and we should either ignore it or use
	 * max of the two values? Because why fail if both pieces have nbucket well
	 * within those limits?
	 */
	if (state->maxbuckets != sketch->maxbuckets)
		elog(ERROR, "state and sketch are not compatible: nbuckets %d != %d",
			 state->maxbuckets, sketch->maxbuckets);

	/*
	 * Walk through the sketch buckets, and see if there's a matching bucket
	 * in the aggregate state. If yes, just increment the counters. If not,
	 * add the bucket at the end.
	 *
	 * XXX We remember the number of used buckets, so that we don't search
	 * in the newly added ones.
	 */
	nbuckets_used = state->nbuckets_used;

	for (i = 0; i < sketch->nbuckets; i++)
	{
		int		j;
		bool	found = false;

		/* FIXME Linear search, should be improved. */
		for (j = 0; j < nbuckets_used; j++)
		{
			if (state->buckets[j].index != sketch->buckets[i].index)
				continue;

			state->buckets[j].count += sketch->buckets[i].index;
			found = true;

			break;
		}

		/* if added to an existing bucket, we're done */
		if (found)
			continue;

		/* otherwise add the bucket at the very end */
		if (state->nbuckets_used == state->nbuckets_allocated)
		{
			/* double the space for buckets, but cap by maxbuckets */
			state->nbuckets_allocated *= 2;

			state->nbuckets_allocated = Min(state->nbuckets_allocated,
											state->maxbuckets);

			/* if still equal, we've reached the maximum */
			if (state->nbuckets_used == state->nbuckets_allocated)
				elog(ERROR, "bucket overflow (used %d, allocated %d, max %d)",
					 state->nbuckets_used,
					 state->nbuckets_allocated,
					 state->maxbuckets);

			/* otherwise reallocate the space to add space */
			state->buckets = repalloc(state->buckets,
									  state->nbuckets_allocated * sizeof(bucket_t));

		}

		/* at this point there has to be space for at least one more bucket */
		Assert(state->nbuckets_used < state->nbuckets_allocated);

		state->buckets[state->nbuckets_used] = sketch->buckets[i];
		state->count += sketch->buckets[i].count;
		state->nbuckets_used++;
	}

	/* sort the buckets by index */
	pg_qsort(state->buckets, state->nbuckets_used, sizeof(bucket_t), bucket_comparator);

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
		int32	maxbuckets = PG_GETARG_INT32(3);

		double *percentiles;
		int		npercentiles;
		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, maxbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		percentiles = array_to_double(fcinfo,
									  PG_GETARG_ARRAYTYPE_P(4),
									  &npercentiles);

		check_percentiles(percentiles, npercentiles);

		state = ddsketch_aggstate_allocate(npercentiles, 0, alpha,
										   maxbuckets, MIN_SKETCH_BUCKETS);

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
		int32	maxbuckets = PG_GETARG_INT32(4);

		double *percentiles;
		int		npercentiles;
		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, maxbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		percentiles = array_to_double(fcinfo,
									  PG_GETARG_ARRAYTYPE_P(5),
									  &npercentiles);

		check_percentiles(percentiles, npercentiles);

		state = ddsketch_aggstate_allocate(npercentiles, 0, alpha,
										   maxbuckets, MIN_SKETCH_BUCKETS);

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
		int32	maxbuckets = PG_GETARG_INT32(3);

		double *values;
		int		nvalues;
		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, maxbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		values = array_to_double(fcinfo,
								 PG_GETARG_ARRAYTYPE_P(4),
								 &nvalues);

		state = ddsketch_aggstate_allocate(0, nvalues, alpha,
										   maxbuckets, MIN_SKETCH_BUCKETS);

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
		int32	maxbuckets = PG_GETARG_INT32(4);

		double *values;
		int		nvalues;
		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, maxbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		values = array_to_double(fcinfo,
								 PG_GETARG_ARRAYTYPE_P(5),
								 &nvalues);

		state = ddsketch_aggstate_allocate(0, nvalues, alpha,
										   maxbuckets, MIN_SKETCH_BUCKETS);

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

	/* can't add values with non-positive counts */
	if (count <= 0)
		elog(ERROR, "invalid count value %ld, must be a positive value", count);

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

		state = ddsketch_aggstate_allocate(npercentiles, 0, sketch->alpha,
										   sketch->maxbuckets, sketch->nbuckets);

		memcpy(state->percentiles, percentiles, sizeof(double) * npercentiles);

		pfree(percentiles);

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	/* copy data from sketch to aggstate */
	for (i = 0; i < sketch->nbuckets; i++)
		state->buckets[i] = sketch->buckets[i];

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

		state = ddsketch_aggstate_allocate(0, nvalues, sketch->alpha,
										   sketch->maxbuckets, sketch->nbuckets);

		memcpy(state->values, values, sizeof(double) * nvalues);

		pfree(values);

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	/* copy data from sketch to aggstate */
	for (i = 0; i < sketch->nbuckets; i++)
		state->buckets[i] = sketch->buckets[i];

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
		  state->nbuckets_used * sizeof(bucket_t);

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
		   sizeof(bucket_t) * state->nbuckets_used);
	ptr += sizeof(bucket_t) * state->nbuckets_used;

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

	state = ddsketch_aggstate_allocate(tmp.npercentiles, tmp.nvalues, tmp.alpha,
									   tmp.maxbuckets, tmp.nbuckets_used);

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
		   sizeof(bucket_t) * state->nbuckets_used);
	ptr += sizeof(bucket_t) * state->nbuckets_used;

	PG_RETURN_POINTER(state);
}

static ddsketch_aggstate_t *
ddsketch_copy(ddsketch_aggstate_t *state)
{
	ddsketch_aggstate_t *copy;

	copy = ddsketch_aggstate_allocate(state->npercentiles, state->nvalues,
									  state->alpha, state->maxbuckets,
									  state->nbuckets_used);

	memcpy(copy, state, offsetof(ddsketch_aggstate_t, percentiles));

	if (state->nvalues > 0)
		memcpy(copy->values, state->values,
			   sizeof(double) * state->nvalues);

	if (state->npercentiles > 0)
		memcpy(copy->percentiles, state->percentiles,
			   sizeof(double) * state->npercentiles);

	memcpy(copy->buckets, state->buckets, state->nbuckets_used * sizeof(bucket_t));

	return copy;
}

Datum
ddsketch_combine(PG_FUNCTION_ARGS)
{
	int					i,
						nbuckets_used;
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

	/*
	 * Walk through the source buckets, and see if there's a matching bucket
	 * in the destination state. If yes, just increment the counters. If not,
	 * add the bucket at the end.
	 *
	 * XXX We remember the number of destination buckets, so that we don't
	 * search in the new ones (which are only sorted at the very end).
	 */
	nbuckets_used = dst->nbuckets_used;

	for (i = 0; i < src->nbuckets_used; i++)
	{
		int		j;
		bool	found = false;

		/* FIXME Linear search, should be improved. */
		for (j = 0; j < nbuckets_used; j++)
		{
			if (dst->buckets[j].index != src->buckets[i].index)
				continue;

			dst->buckets[j].count += src->buckets[i].index;
			found = true;

			break;
		}

		/* if added to an existing bucket, we're done */
		if (found)
			continue;

		/* otherwise add the bucket at the very end */
		if (dst->nbuckets_used == src->nbuckets_allocated)
		{
			/* double the space for buckets, but cap by maxbuckets */
			dst->nbuckets_allocated *= 2;

			dst->nbuckets_allocated = Min(dst->nbuckets_allocated,
										  dst->maxbuckets);

			/* if still equal, we've reached the maximum */
			if (dst->nbuckets_used == dst->nbuckets_allocated)
				elog(ERROR, "bucket overflow (used %d, allocated %d, max %d)",
					 dst->nbuckets_used,
					 dst->nbuckets_allocated,
					 dst->maxbuckets);

			/* otherwise reallocate the space to add space */
			dst->buckets = repalloc(dst->buckets,
									dst->nbuckets_allocated * sizeof(bucket_t));

		}

		/* at this point there has to be space for at least one more bucket */
		Assert(dst->nbuckets_used < dst->nbuckets_allocated);

		dst->buckets[dst->nbuckets_used] = src->buckets[i];
		dst->count += src->buckets[i].count;
		dst->nbuckets_used++;
	}

	/* sort the buckets by index */
	pg_qsort(dst->buckets, dst->nbuckets_used, sizeof(bucket_t), bucket_comparator);

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

	state = ddsketch_aggstate_allocate(0, 0, sketch->alpha,
									   sketch->maxbuckets, sketch->nbuckets);

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
		int32	maxbuckets;

		/*
		 * We don't require compression, but only when there is an existing
		 * ddsketch value. Make sure the value was supplied.
		 */
		if (PG_ARGISNULL(2))
			elog(ERROR, "alpha value not supplied, but ddsketch is NULL");

		if (PG_ARGISNULL(3))
			elog(ERROR, "nbuckets value not supplied, but ddsketch is NULL");

		alpha = PG_GETARG_FLOAT8(2);
		maxbuckets = PG_GETARG_INT32(3);

		check_sketch_parameters(alpha, maxbuckets);

		state = ddsketch_aggstate_allocate(0, 0, alpha, maxbuckets, MIN_SKETCH_BUCKETS);
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
		int		maxbuckets;

		/*
		 * We don't require compression, but only when there is an existing
		 * ddsketch value. Make sure the value was supplied.
		 */
		if (PG_ARGISNULL(2))
			elog(ERROR, "alpha value not supplied, but ddsketch is NULL");

		if (PG_ARGISNULL(3))
			elog(ERROR, "nbuckets value not supplied, but ddsketch is NULL");

		alpha = PG_GETARG_FLOAT8(2);
		maxbuckets = PG_GETARG_INT32(3);

		check_sketch_parameters(alpha, maxbuckets);

		state = ddsketch_aggstate_allocate(0, 0, alpha,
										   maxbuckets, MIN_SKETCH_BUCKETS);
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
		state->buckets[i] = sketch->buckets[i];

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

	sketch = ddsketch_allocate(count, alpha, maxbuckets, nbuckets);

	sketch->flags = flags;
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

		sketch->buckets[i].index = index;
		sketch->buckets[i].count = count;
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
		if (sketch->buckets[i].count == 0)
			continue;

		appendStringInfo(&str, " (%d, " INT64_FORMAT ")", sketch->buckets[i].index, sketch->buckets[i].count);
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

	sketch = ddsketch_allocate(count, alpha, maxbuckets, nbuckets);

	sketch->flags = flags;

	for (i = 0; i < sketch->nbuckets; i++)
	{
		sketch->buckets[i].index = pq_getmsgint(buf, sizeof(int32));
		sketch->buckets[i].count = pq_getmsgint64(buf);
	}

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
	{
		pq_sendint(&buf, sketch->buckets[i].index, 4);
		pq_sendint64(&buf, sketch->buckets[i].count);
	}

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

static double
ddsketch_log_gamma(ddsketch_aggstate_t *state, double value)
{
	return log(value) / log(2.0) * state->multiplier;
}

static int
ddsketch_map_index(ddsketch_aggstate_t *state, double value)
{
	return (int)(ceil(ddsketch_log_gamma(state, value)) + state->offset);
}
