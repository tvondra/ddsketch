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
 * We store just non-empty buckets. Most sketches tend to be sparse, i.e.
 * only a small fraction of buckets is non-empty. Consider e.g. a sketch
 * of API latencies. The values are probably from a fairly narrow range,
 * hence only very few buckets will be non-empty. In particular, buckets
 * close to 0 tend to be empty, because the latency is usually non-zero,
 * and when something breaks the latency increases significantly, leaving
 * a significant gap of empty buckets. So by not storing the empty buckets
 * we usually save quite a bit of space. We could store all buckets and
 * rely on TOAST compression to fix this, but this seems more reliable.
 *
 * XXX We could use varint instead of uint16/int64 to store the buckets,
 * which would help if most buckets are almost empty.
 */
typedef struct bucket_t {
	int32	index;
	int64	count;
} bucket_t;

#define		BUCKETS_BYTES(cnt)		((cnt) * sizeof(bucket_t))

typedef struct ddsketch_t {
	int32		vl_len_;		/* varlena header (do not touch directly!) */
	int32		flags;			/* reserved for future use (versioning, ...) */
	int64		count;			/* number of items added to the ddsketch */
	float8		alpha;			/* alpha used to size the buckets */
	int32		maxbuckets;		/* maximum number of buckets sketch */
	int64		zero_count;		/* zero buckets */
	int32		nbuckets;		/* number of buckets / total */
	int32		nbuckets_negative;	/* number of buckets / negative part */
	bucket_t	buckets[FLEXIBLE_ARRAY_MEMBER];
} ddsketch_t;

#define	SKETCH_BUCKETS(sketch) ((sketch)->buckets)
#define	SKETCH_BUCKETS_BYTES(sketch) (BUCKETS_BYTES((sketch)->nbuckets))

#define	SKETCH_BUCKETS_NEGATIVE(sketch) ((sketch)->buckets)
#define	SKETCH_BUCKETS_NEGATIVE_COUNT(sketch) ((sketch)->nbuckets_negative)

#define	SKETCH_BUCKETS_POSITIVE(sketch) ((sketch)->buckets + (sketch)->nbuckets_negative)
#define	SKETCH_BUCKETS_POSITIVE_COUNT(sketch) ((sketch)->nbuckets - (sketch)->nbuckets_negative)

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

	/* pre-calculated parameters for mapping etc. */
	int32		offset;
	double		min_indexable_value;
	double		max_indexable_value;
	double		multiplier;
	double		gamma;

	/* store with buckets (positive and negative) */
	int64		zero_count;		/* values close to zero */
	int32		maxbuckets;		/* maximum number of buckets */
	int32		nbuckets;		/* number of buckets (used) */
	int32		nbuckets_negative;	/* number of buckets in negative part */
	int32		nbuckets_allocated;	/* number of buckets (allocated) */
	bucket_t   *buckets;		/* buckets (negative and positive) */

	/* array of requested percentiles and values */
	int			npercentiles;	/* number of percentiles */
	int			nvalues;		/* number of values */

	/* variable-length fields at the end */
	double	   *percentiles;	/* array of percentiles (if any) */
	double	   *values;			/* array of values (if any) */
} ddsketch_aggstate_t;

#define STATE_BUCKETS_FULL(state)	((state)->nbuckets == (state)->nbuckets_allocated)

#define	STATE_BUCKETS_USED(state)			((state)->nbuckets)
#define	STATE_BUCKETS_BYTES(state)			BUCKETS_BYTES(STATE_BUCKETS_USED(state))

#define	STATE_BUCKETS_NEGATIVE_COUNT(state)	((state)->nbuckets_negative)
#define	STATE_BUCKETS_POSITIVE_COUNT(state)	((state)->nbuckets - (state)->nbuckets_negative)
#define	STATE_BUCKETS_POSITIVE_BYTES(state)	BUCKETS_BYTES(STATE_BUCKETS_POSITIVE_COUNT(state))

#define	STATE_BUCKETS(state)			((state)->buckets)
#define	STATE_BUCKETS_NEGATIVE(state)	((state)->buckets)
#define	STATE_BUCKETS_POSITIVE(state)	((state)->buckets + (state)->nbuckets_negative)


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
	Assert(sketch->nbuckets_negative <= MAX_SKETCH_BUCKETS);

	Assert(sketch->maxbuckets >= sketch->nbuckets);
	Assert(sketch->nbuckets >= sketch->nbuckets_negative);
	Assert(sketch->nbuckets_negative >= 0);

	count = sketch->zero_count;

	for (i = 0; i < sketch->nbuckets_negative; i++)
	{
		Assert((i == 0) || (sketch->buckets[i-1].index > sketch->buckets[i].index));
		Assert(sketch->buckets[i].count > 0);
		count += sketch->buckets[i].count;
	}

	for (i = sketch->nbuckets_negative; i < sketch->nbuckets; i++)
	{
		Assert((i == sketch->nbuckets_negative) || (sketch->buckets[i-1].index < sketch->buckets[i].index));
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
	int		i;
	int64	count;

	Assert(state->alpha >= MIN_SKETCH_ALPHA);
	Assert(state->alpha <= MAX_SKETCH_ALPHA);

	Assert(state->maxbuckets >= MIN_SKETCH_BUCKETS);
	Assert(state->maxbuckets <= MAX_SKETCH_BUCKETS);

	Assert(state->nbuckets <= MAX_SKETCH_BUCKETS);

	Assert(state->nbuckets_allocated <= MAX_SKETCH_BUCKETS);

	Assert(state->maxbuckets >= state->nbuckets_allocated);
	Assert(state->nbuckets_allocated >= state->nbuckets);
	Assert(state->nbuckets >= state->nbuckets_negative);
	Assert(state->nbuckets_negative >= 0);

	count = 0;
	for (i = 0; i < state->nbuckets_negative; i++)
	{
		/* negative part has indexes in decreasing order */
		Assert((i == 0) || (state->buckets[i-1].index > state->buckets[i].index));
		Assert(state->buckets[i].count > 0);
		count += state->buckets[i].count;
	}

	for (i = state->nbuckets_negative; i < state->nbuckets; i++)
	{
		/* positive part has indexes in increasing order */
		Assert((i == state->nbuckets_negative) || (state->buckets[i-1].index < state->buckets[i].index));
		Assert(state->buckets[i].count > 0);
		count += state->buckets[i].count;
	}

	Assert(state->npercentiles >= 0);
	Assert(state->nvalues >= 0);

	Assert(!((state->npercentiles > 0) && (state->nvalues > 0)));
#endif
}

/*
 * Estimate requested quantiles from the sketch agg state.
 */
static void
ddsketch_compute_quantiles(ddsketch_aggstate_t *state, double *result)
{
	int			i;

	AssertCheckDDSketchAggState(state);

	for (i = 0; i < state->npercentiles; i++)
	{
		int		j;
		int		index = 0;
		int64	count = 0;
		double	goal = (state->percentiles[i] * state->count);
		bucket_t *buckets;

		/*
		 * Process the negative, zero and positive stores, in this order.
		 */
		buckets = STATE_BUCKETS_NEGATIVE(state);
		for (j = 0; j < STATE_BUCKETS_NEGATIVE_COUNT(state); j++)
		{
			/* accumulate the count, remember the last bucket index */
			count += buckets[j].count;
			index = buckets[j].index;

			if (count > goal)
				break;
		}

		/* are we done after processing the negative store? */
		if (count > goal)
		{
			result[i] = - 2 * pow(state->gamma, index) / (state->gamma + 1);
			continue;
		}

		/* now the zero bucket */
		count += state->zero_count;

		/* are we done after processing the zero bucket? */
		if (count > goal)
		{
			result[i] = 0; /* FIXME is it correct to just use 0? */
			continue;
		}

		/* and finally the positive store */
		buckets = STATE_BUCKETS_POSITIVE(state);
		for (j = 0; j < STATE_BUCKETS_POSITIVE_COUNT(state); j++)
		{
			count += buckets[j].count;
			index = buckets[j].index;

			if (count > goal)
				break;
		}

		Assert(count >= goal);

		result[i] = 2 * pow(state->gamma, index) / (state->gamma + 1);
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

	AssertCheckDDSketchAggState(state);

	for (i = 0; i < state->nvalues; i++)
	{
		int64	count = 0;
		double	value = state->values[i];

		if (value > state->min_indexable_value)	/* value in positive part */
		{
			int		j;
			int		index = ddsketch_map_index(state, value);
			bucket_t *buckets;

			/* add the whole negative part */
			buckets = STATE_BUCKETS_NEGATIVE(state);
			for (j = 0; j < STATE_BUCKETS_NEGATIVE_COUNT(state); j++)
				count += buckets[j].count;

			/* add the zero bucket */
			count += state->zero_count;

			/* and now add the positive part, up to the index */
			buckets = STATE_BUCKETS_POSITIVE(state);
			for (j = 0; j < STATE_BUCKETS_POSITIVE_COUNT(state); j++)
			{
				if (buckets[j].index > index)
					break;

				if (buckets[j].index < index)
					count += buckets[j].count;
				else
					/* FIXME should this add just half the bucket? */
					count += buckets[j].count / 2;
			}
		}
		else if (value < -state->min_indexable_value)	/* value in negative part */
		{
			int		j;
			int		index = ddsketch_map_index(state, -value);
			bucket_t *buckets;

			buckets = STATE_BUCKETS_NEGATIVE(state);
			for (j = 0; j < STATE_BUCKETS_NEGATIVE_COUNT(state); j++)
			{
				/* negative part is sorted in reverse order */
				if (buckets[j].index < index)
					break;

				if (buckets[j].index > index)
					count += buckets[j].count;
				else
					/* FIXME should this add just half the bucket? */
					count += buckets[j].count / 2;
			}
		}
		else
		{
			int		j;
			bucket_t *buckets;

			/* add the whole negative part */
			buckets = STATE_BUCKETS_NEGATIVE(state);
			for (j = 0; j < STATE_BUCKETS_NEGATIVE_COUNT(state); j++)
				count += buckets[j].count;

			/* add the zero bucket */
			count += state->zero_count;
		}

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

static int
bucket_comparator_reverse(const void *a, const void *b)
{
	bucket_t *ba = (bucket_t *) a;
	bucket_t *bb = (bucket_t *) b;

	if (ba->index < bb->index)
		return 1;
	else if (ba->index > bb->index)
		return -1;

	return 0;
}

static void
ddsketch_store_add(ddsketch_aggstate_t *state, bool positive, int index, int64 count)
{
	int			i;
	bucket_t   *buckets;
	int			nbuckets;

	/*
	 * See if we already have a bucket with the calculated index. We search
	 * either in the negative or positive part of the array.
	 */
	if (positive)	/* positive part */
	{
		buckets = STATE_BUCKETS_POSITIVE(state);
		nbuckets = STATE_BUCKETS_POSITIVE_COUNT(state);
	}
	else	/* negative part */
	{
		buckets = STATE_BUCKETS_NEGATIVE(state);
		nbuckets = STATE_BUCKETS_NEGATIVE_COUNT(state);
	}

	/*
	 * XXX Linear search - should be optimized to use bsearch or some kind
	 * of simple hash table.
	 */
	for (i = 0; i < nbuckets; i++)
	{
		/* If we found an existing bucket, we're done. */
		if (buckets[i].index == index)
		{
			buckets[i].count += count;
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
	if (STATE_BUCKETS_FULL(state))
	{
		/* double the space for buckets, but cap by maxbuckets */
		state->nbuckets_allocated *= 2;

		/* cap it by the maximum allowed number of buckets */
		state->nbuckets_allocated = Min(state->nbuckets_allocated,
										state->maxbuckets);

		/* if still full, we've reached the maximum */
		if (STATE_BUCKETS_FULL(state))
			elog(ERROR, "bucket overflow (used %d, allocated %d, max %d)",
				 state->nbuckets,
				 state->nbuckets_allocated,
				 state->maxbuckets);

		/* otherwise reallocate the space to add space */
		state->buckets = repalloc(state->buckets,
								  BUCKETS_BYTES(state->nbuckets_allocated));
	}

	/* at this point there has to be space for at least one more bucket */
	Assert(state->nbuckets <= state->maxbuckets);

	/*
	 * At this point we know there's space for a new bucket at the end of
	 * the buckets array. For the positive part we can simply append the
	 * bucket at the end, while for the negative part we have to move the
	 * positive array by one bucket.
	 *
	 * XXX We could also move everything and add it at the beginning, in a
	 * symmetric way to the positive part.
	 */
	if (positive)
	{
		bucket_t   *buckets = STATE_BUCKETS_POSITIVE(state);

		buckets[STATE_BUCKETS_POSITIVE_COUNT(state)].index = index;
		buckets[STATE_BUCKETS_POSITIVE_COUNT(state)].count = count;

		STATE_BUCKETS_USED(state)++;

		/* now sort the positive buckets by index */
		pg_qsort(STATE_BUCKETS_POSITIVE(state),
				 STATE_BUCKETS_POSITIVE_COUNT(state),
				 sizeof(bucket_t), bucket_comparator);
	}
	else
	{
		bucket_t   *buckets = STATE_BUCKETS_NEGATIVE(state);

		/* move the positive buckets to allow adding negative bucket */
		memmove(STATE_BUCKETS_POSITIVE(state) + 1,
				STATE_BUCKETS_POSITIVE(state),
				STATE_BUCKETS_POSITIVE_BYTES(state));

		buckets[STATE_BUCKETS_NEGATIVE_COUNT(state)].index = index;
		buckets[STATE_BUCKETS_NEGATIVE_COUNT(state)].count = count;

		STATE_BUCKETS_USED(state)++;
		STATE_BUCKETS_NEGATIVE_COUNT(state)++;

		/* sort the negative buckets by index in reverse */
		pg_qsort(STATE_BUCKETS_NEGATIVE(state),
				 STATE_BUCKETS_NEGATIVE_COUNT(state),
				 sizeof(bucket_t), bucket_comparator_reverse);
	}
}

/* add a value to the sketch */
static void
ddsketch_add(ddsketch_aggstate_t *state, double v, int64 count)
{
	int		index;

	AssertCheckDDSketchAggState(state);

	state->count += count;

	if (v > state->min_indexable_value)
	{
		index = ddsketch_map_index(state, v);
		ddsketch_store_add(state, true, index, count);
	}
	else if (v < -state->min_indexable_value)
	{
		index = ddsketch_map_index(state, -v);
		ddsketch_store_add(state, false, index, count);
	}
	else
	{
		/* FIXME increment the zero bucket */
		state->zero_count += count;
	}

	AssertCheckDDSketchAggState(state);
}

/*
 * ddsketch_allocate
 *		allocate sketch with enough space for a requested number of buckets
 *
 * We allocate space only for nbuckets buckets, but we also store the maximum
 * allowed number of buckets the sketch is allowed to use.
 */
static ddsketch_t *
ddsketch_allocate(int64 count, double alpha, int64 zero_count, int maxbuckets,
				  int nbuckets, int nbuckets_negative)
{
	Size		len;
	ddsketch_t *sketch;
	char	   *ptr;

	Assert(nbuckets_negative >= 0);
	Assert(nbuckets_negative <= maxbuckets);

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
	sketch->nbuckets_negative = nbuckets_negative;
	sketch->alpha = alpha;
	sketch->zero_count = zero_count;

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
	int					nbuckets_allocated;

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

	nbuckets_allocated = 1;
	while (nbuckets_allocated < nbuckets)
		nbuckets_allocated *= 2;

	Assert(nbuckets_allocated <= maxbuckets);

	/* initialize the bucket store */
	state->maxbuckets = maxbuckets;
	state->nbuckets_allocated = nbuckets_allocated;
	state->nbuckets = 0;

	/* we may need to repalloc this later */
	state->buckets = palloc(nbuckets_allocated * sizeof(bucket_t));

	state->zero_count = 0;

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

	sketch = ddsketch_allocate(state->count,
							   state->alpha,
							   state->zero_count,
							   state->maxbuckets,
							   state->nbuckets,
							   STATE_BUCKETS_NEGATIVE_COUNT(state));

	memcpy(sketch->buckets, STATE_BUCKETS(state), STATE_BUCKETS_BYTES(state));

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
										   maxbuckets,
										   MIN_SKETCH_BUCKETS);

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
										   maxbuckets,
										   MIN_SKETCH_BUCKETS);

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
										   maxbuckets,
										   MIN_SKETCH_BUCKETS);

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
										   maxbuckets,
										   MIN_SKETCH_BUCKETS);

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

/* merge buckets into the aggregate state */
static void
ddsketch_merge_buckets(ddsketch_aggstate_t *state,
					   bool positive, bucket_t *buckets, int nbuckets)
{
	int			i,
				j;
	int			n;
	bucket_t   *b;

	/*
	 * We simply copy buckets from both sources into a single array, sort
	 * it and then combine buckets with the same index. Then we copy the
	 * buckets back into the state (might require enlarging).
	 */
	if (positive)
	{
		n = STATE_BUCKETS_POSITIVE_COUNT(state) + nbuckets;
		b = (bucket_t *) palloc(BUCKETS_BYTES(n));

		/* copy the new buckets */
		memcpy(b, buckets, nbuckets);

		/* copy the existing positive buckets */
		memcpy(b + nbuckets, STATE_BUCKETS_POSITIVE(state),
			   STATE_BUCKETS_POSITIVE_BYTES(state));

		/* sort the combined array */
		pg_qsort(b, n, sizeof(bucket_t), bucket_comparator);
	}
	else
	{
		n = STATE_BUCKETS_NEGATIVE_COUNT(state) + nbuckets;
		b = (bucket_t *) palloc(BUCKETS_BYTES(n));

		/* copy the new buckets */
		memcpy(b, buckets, nbuckets);

		/* copy the existing negative buckets */
		memcpy(b + nbuckets, STATE_BUCKETS_NEGATIVE(state),
			   STATE_BUCKETS_NEGATIVE_COUNT(state));

		/* sort the combined array (in reverse, as it's negative) */
		pg_qsort(b, n, sizeof(bucket_t), bucket_comparator_reverse);
	}

	/* walk through the sorted array and combine buckets with equal index */
	j = 0;
	for (i = 1; i < n; i++)
	{
		/* same as preceding bucket, so merge it into */
		if (b[i].index == b[j].index)
		{
			b[j].count += b[i].count;
			continue;
		}

		/* we'll need to store it in new slot */
		b[++j] = b[i];
	}

	/* remember the number of combined buckets */
	n = j;

	/* how many total buckets we'll need in the aggstate */
	if (positive)
		nbuckets = STATE_BUCKETS_NEGATIVE_COUNT(state) + n;
	else
		nbuckets = STATE_BUCKETS_POSITIVE_COUNT(state) + n;

	/* check if we exceed the allowd number of buckets */
	if (nbuckets > state->maxbuckets)
		elog(ERROR, "too many buckets needed %d > %d",
			 nbuckets, state->maxbuckets);

	/* grow the number of buckets to allocate */
	while (state->nbuckets_allocated < nbuckets)
		state->nbuckets_allocated *= 2;

	/* don't exceed the allowed number (maxbuckets may not be power of two) */
	state->nbuckets_allocated = Min(state->nbuckets_allocated,
									state->maxbuckets);

	state->buckets = repalloc(state->buckets,
							  BUCKETS_BYTES(state->nbuckets_allocated));

	/*
	 * Copy the sorted array back into the state (for the negative case we
	 * need to shift the positive part).
	 *
	 * XXX The array may be exactly the same size as the old one, in which
	 * case we only need to do the memcpy (and the rest is mostly no-op).
	 */
	if (positive)
	{
		/* copy the merged positive array into the state */
		memcpy(STATE_BUCKETS_POSITIVE(state), b, BUCKETS_BYTES(n));

		/* update the number of total buckets (determines positive) */
		STATE_BUCKETS_USED(state) = STATE_BUCKETS_NEGATIVE_COUNT(state) + n;
	}
	else
	{
		/* make sure there's space for the new negative array */
		memmove(STATE_BUCKETS_NEGATIVE(state) + n,
				STATE_BUCKETS_POSITIVE(state),
				STATE_BUCKETS_POSITIVE_BYTES(state));

		/* copy the negative array */
		memcpy(STATE_BUCKETS_NEGATIVE(state), b, BUCKETS_BYTES(n));

		/* update number of negative buckets (and total) */
		STATE_BUCKETS_USED(state) = STATE_BUCKETS_POSITIVE_COUNT(state) + n;
		STATE_BUCKETS_NEGATIVE_COUNT(state) = n;
	}
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

		state = ddsketch_aggstate_allocate(npercentiles, 0, sketch->alpha,
										   sketch->maxbuckets,
										   sketch->nbuckets);

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

	ddsketch_merge_buckets(state, false,
						   SKETCH_BUCKETS_NEGATIVE(sketch),
						   SKETCH_BUCKETS_NEGATIVE_COUNT(sketch));

	ddsketch_merge_buckets(state, true,
						   SKETCH_BUCKETS_POSITIVE(sketch),
						   SKETCH_BUCKETS_POSITIVE_COUNT(sketch));

	state->zero_count += sketch->zero_count;
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
										   sketch->maxbuckets,
										   sketch->nbuckets);

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

	ddsketch_merge_buckets(state, false,
						   SKETCH_BUCKETS_NEGATIVE(sketch),
						   SKETCH_BUCKETS_NEGATIVE_COUNT(sketch));

	ddsketch_merge_buckets(state, true,
						   SKETCH_BUCKETS_POSITIVE(sketch),
						   SKETCH_BUCKETS_POSITIVE_COUNT(sketch));

	state->zero_count += sketch->zero_count;
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
										   maxbuckets,
										   MIN_SKETCH_BUCKETS);

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
										   maxbuckets,
										   MIN_SKETCH_BUCKETS);

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
										   maxbuckets,
										   MIN_SKETCH_BUCKETS);

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
										   maxbuckets,
										   MIN_SKETCH_BUCKETS);

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
										   sketch->maxbuckets,
										   sketch->nbuckets);

		memcpy(state->percentiles, percentiles, sizeof(double) * npercentiles);

		pfree(percentiles);

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	/* check that the sketch and aggstate are compatible */

	if (state->alpha != sketch->alpha)
		elog(ERROR, "state and sketch are not compatible: alpha %lf != %lf",
			 state->alpha, sketch->alpha);

	ddsketch_merge_buckets(state, false,
						   SKETCH_BUCKETS_NEGATIVE(sketch),
						   SKETCH_BUCKETS_NEGATIVE_COUNT(sketch));

	ddsketch_merge_buckets(state, true,
						   SKETCH_BUCKETS_POSITIVE(sketch),
						   SKETCH_BUCKETS_POSITIVE_COUNT(sketch));

	state->zero_count += sketch->zero_count;
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
										   sketch->maxbuckets,
										   sketch->nbuckets);

		memcpy(state->values, values, sizeof(double) * nvalues);

		pfree(values);

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	/* check that the sketch and aggstate are compatible */

	if (state->alpha != sketch->alpha)
		elog(ERROR, "state and sketch are not compatible: alpha %lf != %lf",
			 state->alpha, sketch->alpha);

	ddsketch_merge_buckets(state, false,
						   SKETCH_BUCKETS_NEGATIVE(sketch),
						   SKETCH_BUCKETS_NEGATIVE_COUNT(sketch));

	ddsketch_merge_buckets(state, true,
						   SKETCH_BUCKETS_POSITIVE(sketch),
						   SKETCH_BUCKETS_POSITIVE_COUNT(sketch));

	state->zero_count += sketch->zero_count;
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
		  STATE_BUCKETS_BYTES(state);

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
	memcpy(ptr, STATE_BUCKETS(state), STATE_BUCKETS_BYTES(state));
	ptr += STATE_BUCKETS_BYTES(state);

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
									   tmp.maxbuckets,
									   tmp.nbuckets);

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
	memcpy(STATE_BUCKETS(state), ptr, STATE_BUCKETS_BYTES(state));
	ptr += STATE_BUCKETS_BYTES(state);

	PG_RETURN_POINTER(state);
}

static ddsketch_aggstate_t *
ddsketch_copy(ddsketch_aggstate_t *state)
{
	ddsketch_aggstate_t *copy;

	copy = ddsketch_aggstate_allocate(state->npercentiles, state->nvalues,
									  state->alpha,
									  state->maxbuckets,
									  state->nbuckets);

	memcpy(copy, state, offsetof(ddsketch_aggstate_t, percentiles));

	if (state->nvalues > 0)
		memcpy(copy->values, state->values,
			   sizeof(double) * state->nvalues);

	if (state->npercentiles > 0)
		memcpy(copy->percentiles, state->percentiles,
			   sizeof(double) * state->npercentiles);

	memcpy(STATE_BUCKETS(copy), STATE_BUCKETS(state), STATE_BUCKETS_BYTES(state));

	return copy;
}

Datum
ddsketch_combine(PG_FUNCTION_ARGS)
{
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

	ddsketch_merge_buckets(dst, false,
						   STATE_BUCKETS_NEGATIVE(src),
						   STATE_BUCKETS_NEGATIVE_COUNT(src));

	ddsketch_merge_buckets(dst, true,
						   STATE_BUCKETS_POSITIVE(src),
						   STATE_BUCKETS_POSITIVE_COUNT(src));

	dst->zero_count += src->zero_count;
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

	state = ddsketch_aggstate_allocate(0, 0, sketch->alpha,
									   sketch->maxbuckets,
									   sketch->nbuckets);

	state->count = sketch->count;

	/* copy data from the ddsketch into the aggstate */
	memcpy(STATE_BUCKETS(state), SKETCH_BUCKETS_NEGATIVE(sketch),
		   SKETCH_BUCKETS_BYTES(sketch));

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

		state = ddsketch_aggstate_allocate(0, 0, alpha,
										   maxbuckets,
										   MIN_SKETCH_BUCKETS);
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
										   maxbuckets,
										   MIN_SKETCH_BUCKETS);
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
	ddsketch_merge_buckets(state, false,
						   SKETCH_BUCKETS_NEGATIVE(sketch),
						   SKETCH_BUCKETS_NEGATIVE_COUNT(sketch));

	ddsketch_merge_buckets(state, true,
						   SKETCH_BUCKETS_POSITIVE(sketch),
						   SKETCH_BUCKETS_POSITIVE_COUNT(sketch));

	state->zero_count += sketch->zero_count;
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
	int64		zero_count;
	double		alpha;
	int			maxbuckets;
	int			nbuckets;
	int			nbuckets_negative;
	int			header_length;
	char	   *ptr;

	r = sscanf(str, "flags %d count " INT64_FORMAT " alpha %lf zero_count " INT64_FORMAT " maxbuckets %d buckets %d %d%n",
			   &flags, &count, &alpha, &zero_count, &maxbuckets,
			   &nbuckets, &nbuckets_negative, &header_length);

	if (r != 7)
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
						maxbuckets, MIN_SKETCH_BUCKETS, MAX_SKETCH_BUCKETS)));

	if (nbuckets_negative > nbuckets)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("number of negative buckets (%d) for ddsketch must not exceed nbuckets (%d)",
						nbuckets_negative, nbuckets)));

	if (nbuckets > maxbuckets)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("number of buckets (%d) for ddsketch must not exceed maxbuckets (%d)",
						nbuckets, maxbuckets)));

	if (count <= 0)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("count value for the ddsketch must be positive")));

	if (zero_count <= 0)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("zero_count value for the ddsketch must be positive")));

	if (count < zero_count)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("zero_count value for the ddsketch must not exceed count")));

	sketch = ddsketch_allocate(count, alpha, zero_count,
							   maxbuckets, nbuckets, nbuckets_negative);

	sketch->flags = flags;
	sketch->count = 0;

	ptr = str + header_length;

	i = 0;
	while (true)
	{
		int		index;

		if (sscanf(ptr, " (%d, " INT64_FORMAT ")", &index, &count) != 2)
			break;

		if (i >= nbuckets)
			elog(ERROR, "too many buckets parsed");

		/* XXX can we do some check of the index value? */

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

	appendStringInfo(&str, "flags %d count " INT64_FORMAT " alpha %lf maxbuckets %d buckets %d %d",
					 sketch->flags, sketch->count, sketch->alpha, sketch->maxbuckets,
					 sketch->nbuckets, sketch->nbuckets_negative);

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
	int64		zero_count;
	int32		flags;
	int32		maxbuckets;
	int32		nbuckets_negative;
	int32		nbuckets_positive;
	double		alpha;

	flags = pq_getmsgint(buf, sizeof(int32));

	count = pq_getmsgint64(buf);
	alpha = pq_getmsgfloat8(buf);
	zero_count = pq_getmsgint64(buf);
	maxbuckets = pq_getmsgint(buf, sizeof(int32));
	nbuckets_negative = pq_getmsgint(buf, sizeof(int32));
	nbuckets_positive = pq_getmsgint(buf, sizeof(int32));

	sketch = ddsketch_allocate(count, alpha, zero_count, maxbuckets,
							   nbuckets_negative, nbuckets_positive);

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
	pq_sendint64(&buf, sketch->zero_count);
	pq_sendfloat8(&buf, sketch->alpha);
	pq_sendint(&buf, sketch->maxbuckets, 4);
	pq_sendint(&buf, sketch->nbuckets, 4);
	pq_sendint(&buf, sketch->nbuckets_negative, 4);

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
