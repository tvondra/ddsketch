/*
 * ddsketch - implementation of ddsketch for PostgreSQL
 *
 * DDSketch: A fast and fully-mergeable quantile sketch with relative-error
 * guarantees; Charles Masson, Jee E. Rim, Homin K. Lee;
 * PVLDB, 12(12): 2195-2205, 2019; DOI 10.14778/3352063.3352135
 * https://arxiv.org/abs/1908.10693
 *
 * UDDSketch: Accurate Tracking of Quantiles in Data Streams; Italo
 * Epicoco, Catiuscia Melle, Massimo Cafaro, Marco Pulimeno, Giuseppe
 * Morleo; https://arxiv.org/abs/2004.08604
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
#include "funcapi.h"

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
 * Each bucket stores an index (calculated by the mapping function) and
 * the number of values stored in the bucket. The index is calculated by
 * the mapping function (ddsketch_map_index) and is not the same as index
 * in the bucket array. It may be negative for values close to 0.
 *
 * The array stores buckets for both the positive and negative values, as
 * it's easier to manage (enforcing the maxbuckets limit etc.) than two
 * separate arrays. The array stores negative and positive buckets, in this
 * order. The "nbuckets" tracks the total number of buckets (both parts)
 * and "nbuckets_negative" tracks the negative part.
 *
 * The two parts are sorted by index (separately) - negative buckets in
 * descending order, while positive in ascending order.
 *
 * This struct represents the serialized (on-disk) sketch, and we never
 * add buckets to it. So the array is allocated with exactly the right
 * number of buckets.
 *
 * Values too close to zero can't be represented by the negative/positive
 * buckets, and are stored in a separate "zero bucket" (zero_count).
 *
 * XXX We could use varint instead of uint16/int64 to store the buckets,
 * which would help if most buckets are almost empty.
 *
 * XXX Maybe a bloom filter tracking existing buckets wold be helpful?
 */
typedef struct bucket_t {
	int32	index;	/* mapping index */
	int64	count;	/* bucket counter */
} bucket_t;

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

#define	SKETCH_DEFAULT_FLAGS	0

#define	BUCKETS_BYTES(cnt)	\
	((cnt) * sizeof(bucket_t))

#define	SKETCH_BUCKETS(sketch)	\
	((sketch)->buckets)

#define	SKETCH_BUCKETS_BYTES(sketch)	\
	(BUCKETS_BYTES((sketch)->nbuckets))

#define	SKETCH_BUCKETS_NEGATIVE(sketch)	\
	((sketch)->buckets)

#define	SKETCH_BUCKETS_NEGATIVE_COUNT(sketch)	\
	((sketch)->nbuckets_negative)

#define	SKETCH_BUCKETS_POSITIVE(sketch)	\
	((sketch)->buckets + (sketch)->nbuckets_negative)

#define	SKETCH_BUCKETS_POSITIVE_COUNT(sketch)	\
	((sketch)->nbuckets - (sketch)->nbuckets_negative)

#define	SKETCH_BYTES(sketch)	\
	(offsetof(ddsketch_t, buckets) + BUCKETS_BYTES((sketch)->nbuckets))

#define PG_GETARG_DDSKETCH(x)	\
	(ddsketch_t *) PG_DETOAST_DATUM(PG_GETARG_DATUM(x))

/*
 * An aggregate state, representing the sketch and some additional info
 * (requested percentiles, ...).
 *
 * This is similar to the ddsketch_t struct, but includes various values
 * necessary for mapping values to buckets, determining the indexable
 * range and so on.
 *
 * The array of buckets is allocated using the usual doubling strategy.
 * The allocated space is tracked in nbuckets_allocated, while nbuckets
 * stores the number of buckets actually used. The array is split into
 * negative/positive buckets, just like for ddsketch_t.
 *
 * Values too close to zero can't be represented by the negative/positive
 * buckets, and are stored in a separate "zero bucket" (zero_count).
 *
 * XXX We only ever use one of values/percentiles, never both at the same
 * time. In the future the values may use a different data types than double
 * (e.g. numeric), so we keep both fields.
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

	/* trimmed aggregates */
	double		trim_low;		/* low threshold (for trimmed aggs) */
	double		trim_high;		/* high threshold (for trimmed aggs) */

	/* store with buckets (positive and negative) */
	int64		zero_count;		/* values close to zero */
	int32		maxbuckets;		/* maximum number of buckets */
	int32		nbuckets;		/* number of buckets (used) */
	int32		nbuckets_negative;	/* number of buckets in negative part */
	int32		nbuckets_allocated;	/* number of buckets (allocated) */

	/* array of requested percentiles and values */
	int			npercentiles;	/* number of percentiles */
	int			nvalues;		/* number of values */

	/* variable-length fields at the end */
	double	   *percentiles;	/* array of percentiles (if any) */
	double	   *values;			/* array of values (if any) */
	bucket_t   *buckets;		/* buckets (negative and positive) */
} ddsketch_aggstate_t;

#define STATE_BUCKETS_FULL(state)	\
	((state)->nbuckets == (state)->nbuckets_allocated)

#define	STATE_BUCKETS_USED(state)	\
	((state)->nbuckets)

#define	STATE_BUCKETS_BYTES(state)	\
	BUCKETS_BYTES(STATE_BUCKETS_USED(state))

#define	STATE_BUCKETS_NEGATIVE_COUNT(state)	\
	((state)->nbuckets_negative)

#define	STATE_BUCKETS_NEGATIVE_BYTES(state)	\
	BUCKETS_BYTES(STATE_BUCKETS_NEGATIVE_COUNT(state))

#define	STATE_BUCKETS_POSITIVE_COUNT(state)	\
	((state)->nbuckets - (state)->nbuckets_negative)

#define	STATE_BUCKETS_POSITIVE_BYTES(state)	\
	BUCKETS_BYTES(STATE_BUCKETS_POSITIVE_COUNT(state))

#define	STATE_BUCKETS(state)	\
	((state)->buckets)

#define	STATE_BUCKETS_NEGATIVE(state)	\
	((state)->buckets)

#define	STATE_BUCKETS_POSITIVE(state)	\
	((state)->buckets + (state)->nbuckets_negative)

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
PG_FUNCTION_INFO_V1(ddsketch_add_double_count_increment);
PG_FUNCTION_INFO_V1(ddsketch_add_double_array_increment);
PG_FUNCTION_INFO_V1(ddsketch_union_double_increment);

PG_FUNCTION_INFO_V1(ddsketch_sketch_info);
PG_FUNCTION_INFO_V1(ddsketch_sketch_buckets);
PG_FUNCTION_INFO_V1(ddsketch_param_info);
PG_FUNCTION_INFO_V1(ddsketch_param_buckets);

PG_FUNCTION_INFO_V1(ddsketch_add_double_trimmed);
PG_FUNCTION_INFO_V1(ddsketch_add_double_count_trimmed);
PG_FUNCTION_INFO_V1(ddsketch_add_sketch_trimmed);
PG_FUNCTION_INFO_V1(ddsketch_trimmed_avg);
PG_FUNCTION_INFO_V1(ddsketch_trimmed_sum);

PG_FUNCTION_INFO_V1(ddsketch_sketch_sum);
PG_FUNCTION_INFO_V1(ddsketch_sketch_avg);

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
Datum ddsketch_add_double_count_increment(PG_FUNCTION_ARGS);
Datum ddsketch_add_double_array_increment(PG_FUNCTION_ARGS);
Datum ddsketch_union_double_increment(PG_FUNCTION_ARGS);

Datum ddsketch_sketch_info(PG_FUNCTION_ARGS);
Datum ddsketch_sketch_buckets(PG_FUNCTION_ARGS);
Datum ddsketch_param_info(PG_FUNCTION_ARGS);
Datum ddsketch_param_buckets(PG_FUNCTION_ARGS);

Datum ddsketch_add_double_trimmed(PG_FUNCTION_ARGS);
Datum ddsketch_add_double_count_trimmed(PG_FUNCTION_ARGS);
Datum ddsketch_add_sketch_trimmed(PG_FUNCTION_ARGS);
Datum ddsketch_trimmed_avg(PG_FUNCTION_ARGS);
Datum ddsketch_trimmed_sum(PG_FUNCTION_ARGS);

Datum ddsketch_sketch_sum(PG_FUNCTION_ARGS);
Datum ddsketch_sketch_avg(PG_FUNCTION_ARGS);

static Datum double_to_array(FunctionCallInfo fcinfo, double * d, int len);
static double *array_to_double(FunctionCallInfo fcinfo, ArrayType *v, int * len);

/* mapping to bucket indexes etc. */
static double ddsketch_log_gamma(ddsketch_aggstate_t *state, double value);
static double ddsketch_pow_gamma(ddsketch_aggstate_t *state, double value);
static int    ddsketch_map_index(ddsketch_aggstate_t *state, double value);
static double ddsketch_map_value(ddsketch_aggstate_t *state, double index);

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
	bucket_t *buckets;

	Assert(sketch->flags == SKETCH_DEFAULT_FLAGS);

	Assert(sketch->alpha >= MIN_SKETCH_ALPHA);
	Assert(sketch->alpha <= MAX_SKETCH_ALPHA);

	Assert(sketch->maxbuckets >= MIN_SKETCH_BUCKETS);
	Assert(sketch->maxbuckets <= MAX_SKETCH_BUCKETS);

	Assert(sketch->maxbuckets >= sketch->nbuckets);
	Assert(sketch->nbuckets >= sketch->nbuckets_negative);
	Assert(sketch->nbuckets_negative >= 0);

	count = sketch->zero_count;

	/* negative part */
	buckets = SKETCH_BUCKETS_NEGATIVE(sketch);
	for (i = 0; i < SKETCH_BUCKETS_NEGATIVE_COUNT(sketch); i++)
	{
		/* negative part sorted by index in desdending order */
		Assert((i == 0) || (buckets[i-1].index > buckets[i].index));
		Assert(buckets[i].count > 0);
		count += buckets[i].count;
	}

	/* positive part */
	buckets = SKETCH_BUCKETS_POSITIVE(sketch);
	for (i = 0; i < SKETCH_BUCKETS_POSITIVE_COUNT(sketch); i++)
	{
		/* positive part sorted by index in ascending order */
		Assert((i == 0) || (buckets[i-1].index < buckets[i].index));
		Assert(buckets[i].count > 0);
		count += buckets[i].count;
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
	bucket_t *buckets;

	Assert(state->alpha >= MIN_SKETCH_ALPHA);
	Assert(state->alpha <= MAX_SKETCH_ALPHA);

	Assert(state->maxbuckets >= MIN_SKETCH_BUCKETS);
	Assert(state->maxbuckets <= MAX_SKETCH_BUCKETS);

	Assert(state->maxbuckets >= state->nbuckets_allocated);
	Assert(state->nbuckets_allocated >= state->nbuckets);
	Assert(state->nbuckets >= state->nbuckets_negative);
	Assert(state->nbuckets_negative >= 0);

	count = state->zero_count;

	/* negative part */
	buckets = STATE_BUCKETS_NEGATIVE(state);
	for (i = 0; i < STATE_BUCKETS_NEGATIVE_COUNT(state); i++)
	{
		/* negative part sorted by index in desdending order */
		Assert((i == 0) || (buckets[i-1].index > buckets[i].index));
		Assert(buckets[i].count > 0);
		count += buckets[i].count;
	}

	/* positive part */
	buckets = STATE_BUCKETS_POSITIVE(state);
	for (i = 0; i < STATE_BUCKETS_POSITIVE_COUNT(state); i++)
	{
		/* positive part sorted by index in ascending order */
		Assert((i == 0) || (buckets[i-1].index < buckets[i].index));
		Assert(buckets[i].count > 0);
		count += buckets[i].count;
	}

	// Assert(count == state->count);

	Assert(state->npercentiles >= 0);
	Assert(state->nvalues >= 0);

	/* both can't be set at the same time */
	Assert(!((state->npercentiles > 0) && (state->nvalues > 0)));
#endif
}

/*
 * Estimate requested quantiles from the sketch agggregate state.
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
		double	goal = (state->percentiles[i] * (state->count - 1));
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
			result[i] = -ddsketch_map_value(state, index);;
			continue;
		}

		/* now the zero bucket */
		count += state->zero_count;

		/* are we done after processing the zero bucket? */
		if (count > goal)
		{
			result[i] = 0;
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

		result[i] = ddsketch_map_value(state, index);
	}
}

/*
 * Estimate inverse of quantile given a value from the sketch agg state.
 *
 * Essentially an inverse to ddsketch_compute_quantiles.
 *
 * XXX Unlike ddsketch_compute_quantiles, there's no guarantee regarding
 * errors guarantees - the relative error guarantees are due to sizing
 * the bucket ranges [min,max] in a smart way, so that (max-min)/min is
 * less than the desired error. But we have no control over how many
 * values fall into the bucket, which is what matter for quantiles_of.
 * In extreme case all the values may be in a single bucket, and we don't
 * know if all are below/above the parameter, or what. The best thing
 * we can do is assuming it's in the middle of the bucket.
 *
 * XXX We might also calculate the min/max percentiles, and return a range
 * of possible quantiles (a bit like confidence interval).
 *
 * XXX Maybe instead of using half the bucket, we could use linear
 * approximation between the bucket min/max.
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

		result[i] = count / ((double) state->count - 1);
	}
}

/* comparator by index in ascending order (positive buckets) */
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

/* comparator by index in descending order (negative buckets) */
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

/*
 * Add the number of values to the bucket with the given index in either
 * the positive or negative part of the buckets array.
 *
 * The buckets are sorted by index (each part separately), so we can simply
 * do a binary search. If the bucket already exists, we simply increment the
 * counter and we're done.
 *
 * If the bucket does not exist, we make sure there's enough space for a
 * new bucket (we may enlarge the array) and append it at the end of the
 * appropriate part (negative or positive buckets). And then we sort the
 * buckets, so that it's sorted again.
 *
 * XXX It may seem sorting the whole negative/positive array is expensive,
 * but we expect doing that only very rarely - we should create all the
 * necessary buckets fairly quickly, and then sorting should not be needed.
 */
static void
ddsketch_store_add(ddsketch_aggstate_t *state, bool positive, int index, int64 count)
{
	bucket_t   *bucket;

	/*
	 * See if we already have a bucket with the calculated index. We search
	 * either in the negative or positive part of the array.
	 */
	if (positive)	/* positive part */
	{
		bucket_t	key;

		key.index = index;
		bucket = bsearch(&key,
						 STATE_BUCKETS_POSITIVE(state),
						 STATE_BUCKETS_POSITIVE_COUNT(state),
						 sizeof(bucket_t), bucket_comparator);
	}
	else	/* negative part */
	{
		bucket_t	key;

		key.index = index;
		bucket = bsearch(&key,
						 STATE_BUCKETS_NEGATIVE(state),
						 STATE_BUCKETS_NEGATIVE_COUNT(state),
						 sizeof(bucket_t), bucket_comparator_reverse);
	}

	/* If we found a matching bucket, we're done. */
	if (bucket)
	{
		bucket->count += count;
		return;
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

		/* add a new bucket at the end of the positive part */
		buckets[STATE_BUCKETS_POSITIVE_COUNT(state)].index = index;
		buckets[STATE_BUCKETS_POSITIVE_COUNT(state)].count = count;

		STATE_BUCKETS_USED(state)++;

		/* sort the positive buckets by index (ascending) */
		pg_qsort(STATE_BUCKETS_POSITIVE(state),
				 STATE_BUCKETS_POSITIVE_COUNT(state),
				 sizeof(bucket_t), bucket_comparator);
	}
	else
	{
		bucket_t   *buckets = STATE_BUCKETS_NEGATIVE(state);

		/* move the positive buckets to make space for a negative bucket */
		memmove(STATE_BUCKETS_POSITIVE(state) + 1,
				STATE_BUCKETS_POSITIVE(state),
				STATE_BUCKETS_POSITIVE_BYTES(state));

		/* add a new bucket at the end of the negative part */
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

/*
 * Add a double value to the sketch aggstate. Check if the value belongs
 * to the indexable range or zero bucket. If it can be indexed, add it to
 * the negative or positive part.
 *
 * XXX What about values exceeding the maximum indexable values?
 */
static void
ddsketch_add(ddsketch_aggstate_t *state, double value, int64 count)
{
	int		index;

	AssertCheckDDSketchAggState(state);

	state->count += count;

	if (value > state->min_indexable_value)
	{
		index = ddsketch_map_index(state, value);
		ddsketch_store_add(state, true, index, count);
	}
	else if (value < -state->min_indexable_value)
	{
		index = ddsketch_map_index(state, -value);
		ddsketch_store_add(state, false, index, count);
	}
	else
	{
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
ddsketch_allocate(int32 flags, int64 count, double alpha, int64 zero_count,
				  int maxbuckets, int nbuckets, int nbuckets_negative)
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

	sketch->flags = flags;
	sketch->count = count;
	sketch->maxbuckets = maxbuckets;
	sketch->nbuckets = nbuckets;
	sketch->nbuckets_negative = nbuckets_negative;
	sketch->alpha = alpha;
	sketch->zero_count = zero_count;

	/*
	 * FIXME At this point the number of buckets is set but buckets are not
	 * copied yet, so it's somewhat broken.
	 */

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

	/* we may need to repalloc this later */
	state->buckets = palloc(nbuckets_allocated * sizeof(bucket_t));

	state->nbuckets = 0;
	state->nbuckets_negative = 0;
	state->count = 0;
	state->zero_count = 0;

	/* precalculate various parameters used for mapping */
	state->offset = 0;
	state->gamma = (1 + alpha) / (1 - alpha);
	state->multiplier = log(2.0) / log1p(2 * alpha / (1 - alpha));
	state->min_indexable_value = DBL_MIN * state->gamma;
	state->max_indexable_value = DBL_MAX / state->gamma;

	AssertCheckDDSketchAggState(state);

	return state;
}

/*
 * Serialize the aggregate state into the compact ddsketch representation.
 */
static ddsketch_t *
ddsketch_aggstate_to_ddsketch(ddsketch_aggstate_t *state)
{
	ddsketch_t *sketch;

	AssertCheckDDSketchAggState(state);

	sketch = ddsketch_allocate(SKETCH_DEFAULT_FLAGS,
							   state->count,
							   state->alpha,
							   state->zero_count,
							   state->maxbuckets,
							   state->nbuckets,
							   STATE_BUCKETS_NEGATIVE_COUNT(state));

	memcpy(sketch->buckets, STATE_BUCKETS(state), STATE_BUCKETS_BYTES(state));

	AssertCheckDDSketch(sketch);

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

/* check that the user-specified sketch parameters are valid */
static void
check_sketch_parameters(double alpha, int nbuckets)
{
	if (alpha < MIN_SKETCH_ALPHA || alpha > MAX_SKETCH_ALPHA)
		elog(ERROR, "invalid alpha value %f", alpha);

	if (nbuckets < MIN_SKETCH_BUCKETS || nbuckets > MAX_SKETCH_BUCKETS)
		elog(ERROR, "invalid number of buckets %d", nbuckets);
}

static void
check_trim_values(double low, double high)
{
	if (low < 0.0)
		elog(ERROR, "invalid low percentile value %f, should be in [0.0, 1.0]",
			 low);

	if (high > 1.0)
		elog(ERROR, "invalid high percentile value %f, should be in [0.0, 1.0]",
			 high);

	if (low >= high)
		elog(ERROR, "invalid low/high percentile values %f/%f, should be low < high",
			 low, high);
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
		double	alpha = PG_GETARG_FLOAT8(3);
		int32	maxbuckets = PG_GETARG_INT32(4);

		double *values = NULL;
		int		nvalues = 0;
		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, maxbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		if (PG_NARGS() >= 6)
		{
			values = (double *) palloc(sizeof(double));
			values[0] = PG_GETARG_FLOAT8(5);
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

/* merge buckets into the aggregate state */
static void
ddsketch_merge_buckets(ddsketch_aggstate_t *state,
					   bool positive, bucket_t *buckets, int nbuckets)
{
	int			i,
				j;
	int			n;
	bucket_t   *b;

	if (nbuckets == 0)
		return;

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
		memcpy(b, buckets, BUCKETS_BYTES(nbuckets));

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
		memcpy(b, buckets, BUCKETS_BYTES(nbuckets));

		/* copy the existing negative buckets */
		memcpy(b + nbuckets, STATE_BUCKETS_NEGATIVE(state),
			   STATE_BUCKETS_NEGATIVE_BYTES(state));

		/* sort the combined array (in reverse, as it's negative) */
		pg_qsort(b, n, sizeof(bucket_t), bucket_comparator_reverse);
	}

	/* walk through the sorted array and combine buckets with equal index */
	j = 0;
	for (i = 1; i < n; i++)
	{
		/* not the same as preceding bucket, so a new one */
		if (b[i].index != b[i-1].index)
		{
			b[++j] = b[i];
			continue;
		}

		/* just add it to the "current" bucket  */
		b[j].count += b[i].count;
	}

	/* number of combined buckets (index of last bucket plus one) */
	n = (j+1);

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

	AssertCheckDDSketch(sketch);
	AssertCheckDDSketchAggState(state);

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

	AssertCheckDDSketchAggState(state);

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

	AssertCheckDDSketch(sketch);
	AssertCheckDDSketchAggState(state);

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

	AssertCheckDDSketchAggState(state);

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

	AssertCheckDDSketch(sketch);

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

	AssertCheckDDSketch(sketch);
	AssertCheckDDSketchAggState(state);

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

	AssertCheckDDSketchAggState(state);

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

	AssertCheckDDSketch(sketch);

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

	AssertCheckDDSketch(sketch);
	AssertCheckDDSketchAggState(state);

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

	AssertCheckDDSketchAggState(state);

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
	char   *endptr PG_USED_FOR_ASSERTS_ONLY;
	ddsketch_aggstate_t	tmp;
	ddsketch_aggstate_t *state;
	double			   *percentiles = NULL;
	double			   *values = NULL;

	endptr = ptr + VARSIZE_ANY_EXHDR(v);

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
									   tmp.maxbuckets, tmp.nbuckets);

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

	Assert(ptr == endptr);

	PG_RETURN_POINTER(state);
}

static ddsketch_aggstate_t *
ddsketch_copy(ddsketch_aggstate_t *state)
{
	ddsketch_aggstate_t *copy;

	AssertCheckDDSketchAggState(state);

	copy = ddsketch_aggstate_allocate(state->npercentiles, state->nvalues,
									  state->alpha, state->maxbuckets,
									  state->nbuckets);

	memcpy(copy, state, offsetof(ddsketch_aggstate_t, percentiles));

	if (state->nvalues > 0)
		memcpy(copy->values, state->values,
			   sizeof(double) * state->nvalues);

	if (state->npercentiles > 0)
		memcpy(copy->percentiles, state->percentiles,
			   sizeof(double) * state->npercentiles);

	memcpy(STATE_BUCKETS(copy), STATE_BUCKETS(state), STATE_BUCKETS_BYTES(state));

	AssertCheckDDSketchAggState(copy);

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

		AssertCheckDDSketchAggState(src);

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

	AssertCheckDDSketch(sketch);

	state = ddsketch_aggstate_allocate(0, 0, sketch->alpha,
									   sketch->maxbuckets, sketch->nbuckets);

	state->count = sketch->count;
	state->nbuckets = sketch->nbuckets;

	/* copy data from the ddsketch into the aggstate */
	memcpy(STATE_BUCKETS(state), SKETCH_BUCKETS_NEGATIVE(sketch),
		   SKETCH_BUCKETS_BYTES(sketch));

	AssertCheckDDSketchAggState(state);

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

	AssertCheckDDSketchAggState(state);

	ddsketch_add(state, PG_GETARG_FLOAT8(1), 1);

	AssertCheckDDSketchAggState(state);

	PG_RETURN_POINTER(ddsketch_aggstate_to_ddsketch(state));
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
ddsketch_add_double_count_increment(PG_FUNCTION_ARGS)
{
	int64				count;
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
		if (PG_ARGISNULL(3))
			elog(ERROR, "alpha value not supplied, but ddsketch is NULL");

		if (PG_ARGISNULL(4))
			elog(ERROR, "nbuckets value not supplied, but ddsketch is NULL");

		alpha = PG_GETARG_FLOAT8(3);
		maxbuckets = PG_GETARG_INT32(4);

		check_sketch_parameters(alpha, maxbuckets);

		state = ddsketch_aggstate_allocate(0, 0, alpha,
										   maxbuckets,
										   MIN_SKETCH_BUCKETS);
	}
	else
		state = ddsketch_sketch_to_aggstate(PG_GETARG_DDSKETCH(0));

	if (PG_ARGISNULL(2))
		count = 1;
	else
		count = PG_GETARG_INT64(2);

	AssertCheckDDSketchAggState(state);

	ddsketch_add(state, PG_GETARG_FLOAT8(1), count);

	AssertCheckDDSketchAggState(state);

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

	/* parse the second ddsketch */
	sketch = PG_GETARG_DDSKETCH(1);

	AssertCheckDDSketch(sketch);
	AssertCheckDDSketchAggState(state);

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

	if (flags != SKETCH_DEFAULT_FLAGS)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("invalid sketch flags %d", flags)));

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

	if (nbuckets <= 0)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("number of buckets (%d) for ddsketch must be positive",
						nbuckets)));

	if (nbuckets_negative < 0)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("number of negative buckets (%d) for ddsketch must not be negative",
						nbuckets_negative)));

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

	if (zero_count < 0)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("zero_count value for the ddsketch must be positive")));

	if (count < zero_count)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("zero_count value for the ddsketch must not exceed count")));

	sketch = ddsketch_allocate(flags, count, alpha, zero_count,
							   maxbuckets, nbuckets, nbuckets_negative);

	ptr = str + header_length;

	count = zero_count;

	i = 0;
	while (true)
	{
		int		index;
		int64	bucket_count;

		if (sscanf(ptr, " (%d, " INT64_FORMAT ")", &index, &bucket_count) != 2)
			break;

		if (i >= nbuckets)
			elog(ERROR, "too many buckets parsed");

		/*
		 * Basic checks that the indexes are decreasing in the negative part
		 * and increasing in the positive part.
		 *
		 * XXX Can we check the index value is valid (not too low/high)?
		 */
		if ((i != 0) && (i < nbuckets_negative))
		{
			/* negative store - descending index values */
			if (sketch->buckets[i-1].index <= index)
				ereport(ERROR,
						(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						 errmsg("invalid sketch - ascending indexes in the negative part")));
		}
		else if (i  > nbuckets_negative)
		{
			/* positive store - ascending index values */
			if (sketch->buckets[i-1].index >= index)
				ereport(ERROR,
						(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						 errmsg("invalid sketch - descending indexes in the positive part")));
		}

		/* we don't include empty buckets */
		if (bucket_count <= 0)
			ereport(ERROR,
					(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
					 errmsg("count value for all indexes in a ddsketch must be positive")));

		sketch->buckets[i].index = index;
		sketch->buckets[i].count = bucket_count;

		count += bucket_count;

		/* skip to the end of the centroid */
		ptr = strchr(ptr, ')') + 1;
		i++;
	}

	/* Did we parse exactly the expected number of buckets? */
	if (i != nbuckets)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("parsed invalid number of buckets (%d != %d)",
						i, nbuckets)));

	/* Did we get buckets matching the header? */
	if (count != sketch->count)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("total count (%ld) does not match buckets (%ld)",
						sketch->count, count)));

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

	appendStringInfo(&str, "flags %d count " INT64_FORMAT " alpha %lf zero_count " INT64_FORMAT " maxbuckets %d buckets %d %d",
					 sketch->flags, sketch->count, sketch->alpha, sketch->zero_count,
					 sketch->maxbuckets, sketch->nbuckets, sketch->nbuckets_negative);

	for (i = 0; i < sketch->nbuckets; i++)
		appendStringInfo(&str, " (%d, " INT64_FORMAT ")", sketch->buckets[i].index, sketch->buckets[i].count);

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

	sketch = ddsketch_allocate(flags, count, alpha, zero_count, maxbuckets,
							   nbuckets_negative, nbuckets_positive);

	for (i = 0; i < sketch->nbuckets; i++)
	{
		if (i >= sketch->nbuckets)
			elog(ERROR, "too many buckets parsed");

		sketch->buckets[i].index = pq_getmsgint(buf, sizeof(int32));
		sketch->buckets[i].count = pq_getmsgint64(buf);
	}

	AssertCheckDDSketch(sketch);

	PG_RETURN_POINTER(sketch);
}

Datum
ddsketch_send(PG_FUNCTION_ARGS)
{
	ddsketch_t  *sketch = (ddsketch_t *) PG_DETOAST_DATUM(PG_GETARG_DATUM(0));
	StringInfoData buf;
	int			i;

	AssertCheckDDSketch(sketch);

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

static double
ddsketch_pow_gamma(ddsketch_aggstate_t *state, double value)
{
	return pow(2.0, (value / state->multiplier));
}

static double
ddsketch_map_lower_bound(double alpha, int index)
{
	int		offset = 0;
	double	multiplier = log(2.0) / log1p(2 * alpha / (1 - alpha));

	/* XXX not sure about the ceil() inverse */
	return exp(log(2.0) * ((double) index - offset - 1) / multiplier);
}

static double
ddsketch_map_upper_bound(double alpha, int index)
{
	/* lower bound of the next bucket */
	return ddsketch_map_lower_bound(alpha, index + 1);
}

static int
ddsketch_map_index(ddsketch_aggstate_t *state, double value)
{
	return (int)(ceil(ddsketch_log_gamma(state, value)) + state->offset);
}

static int
ddsketch_map_index2(double alpha, double value)
{
	double	multiplier = log(2.0) / log1p(2 * alpha / (1 - alpha));
	double	log_gamma = log(value) / log(2.0) * multiplier;
	double	offset = 0;

	return (int)(ceil(log_gamma) + offset);
}

static double
ddsketch_map_value(ddsketch_aggstate_t *state, double index)
{
	return ddsketch_pow_gamma(state, index - state->offset) * (2.0 / (1 + state->gamma));
}

Datum
ddsketch_sketch_info(PG_FUNCTION_ARGS)
{
	ddsketch_t *sketch = PG_GETARG_DDSKETCH(0);
	TupleDesc	tupdesc;

	Datum		result;
	HeapTuple	tuple;
	Datum		values[10];
	bool		nulls[10];

	double		gamma;
	double		min_indexable_value;
	double		max_indexable_value;

	/* Build a tuple descriptor for our result type */
	if (get_call_result_type(fcinfo, NULL, &tupdesc) != TYPEFUNC_COMPOSITE)
		elog(ERROR, "return type must be a row type");

	gamma = (1 + sketch->alpha) / (1 - sketch->alpha);
	min_indexable_value = DBL_MIN * gamma;
	max_indexable_value = DBL_MAX / gamma;

	values[0] = UInt64GetDatum(SKETCH_BYTES(sketch));
	values[1] = UInt32GetDatum(sketch->flags);
	values[2] = Float8GetDatum(sketch->alpha);
	values[3] = Int64GetDatum(sketch->count);
	values[4] = Int64GetDatum(sketch->zero_count);
	values[5] = Int32GetDatum(sketch->maxbuckets);
	values[6] = Int32GetDatum(sketch->nbuckets_negative);
	values[7] = Int32GetDatum(sketch->nbuckets - sketch->nbuckets_negative);
	values[8] = Float8GetDatum(min_indexable_value);
	values[9] = Float8GetDatum(max_indexable_value);

	/* Build and return the tuple. */

	memset(nulls, 0, sizeof(nulls));

	tuple = heap_form_tuple(tupdesc, values, nulls);
	result = HeapTupleGetDatum(tuple);

	PG_RETURN_DATUM(result);
}

Datum
ddsketch_sketch_buckets(PG_FUNCTION_ARGS)
{
	ddsketch_t *sketch = PG_GETARG_DDSKETCH(0);
	FuncCallContext *fctx;
	TupleDesc		tupdesc;

	if (SRF_IS_FIRSTCALL())
	{
		MemoryContext mctx;

		fctx = SRF_FIRSTCALL_INIT();

		mctx = MemoryContextSwitchTo(fctx->multi_call_memory_ctx);

		/* Build a tuple descriptor for our result type */
		if (get_call_result_type(fcinfo, NULL, &tupdesc) != TYPEFUNC_COMPOSITE)
			elog(ERROR, "return type must be a row type");

		fctx->user_fctx = tupdesc;
		fctx->max_calls = sketch->nbuckets;

		MemoryContextSwitchTo(mctx);
	}

	fctx = SRF_PERCALL_SETUP();

	if (fctx->call_cntr < fctx->max_calls)
	{
		bucket_t   *bucket = &sketch->buckets[fctx->call_cntr];
		HeapTuple	resultTuple;
		Datum		result;
		Datum		values[6];
		bool		nulls[6];

		double		lower_bound = ddsketch_map_lower_bound(sketch->alpha, bucket->index);
		double		upper_bound = ddsketch_map_upper_bound(sketch->alpha, bucket->index);

		tupdesc = fctx->user_fctx;

		memset(nulls, 0, sizeof(nulls));

		/* Extract information from the line pointer */
		values[0] = Int32GetDatum(fctx->call_cntr);
		values[1] = Int32GetDatum(bucket->index);

		if (fctx->call_cntr > sketch->nbuckets_negative)
		{
			values[2] = Float8GetDatum(lower_bound);
			values[3] = Float8GetDatum(upper_bound);
		}
		else
		{
			values[2] = Float8GetDatum(-upper_bound);
			values[3] = Float8GetDatum(-lower_bound);
		}

		values[4] = Float8GetDatum(fabs(upper_bound - lower_bound));
		values[5] = Int32GetDatum(bucket->count);

		/* Build and return the result tuple. */
		resultTuple = heap_form_tuple(tupdesc, values, nulls);
		result = HeapTupleGetDatum(resultTuple);

		SRF_RETURN_NEXT(fctx, result);
	}
	else
		SRF_RETURN_DONE(fctx);
}

Datum
ddsketch_param_info(PG_FUNCTION_ARGS)
{
	double		alpha = PG_GETARG_FLOAT8(0);
	TupleDesc	tupdesc;

	Datum		result;
	HeapTuple	tuple;
	Datum		values[2];
	bool		nulls[2];

	double		gamma;
	double		min_indexable_value;
	double		max_indexable_value;

	/* Build a tuple descriptor for our result type */
	if (get_call_result_type(fcinfo, NULL, &tupdesc) != TYPEFUNC_COMPOSITE)
		elog(ERROR, "return type must be a row type");

	gamma = (1 + alpha) / (1 - alpha);
	min_indexable_value = DBL_MIN * gamma;
	max_indexable_value = DBL_MAX / gamma;

	values[0] = Float8GetDatum(min_indexable_value);
	values[1] = Float8GetDatum(max_indexable_value);

	/* Build and return the tuple. */

	memset(nulls, 0, sizeof(nulls));

	tuple = heap_form_tuple(tupdesc, values, nulls);
	result = HeapTupleGetDatum(tuple);

	PG_RETURN_DATUM(result);
}

typedef struct ddsketch_buckets_state_t {
	TupleDesc	tupdesc;
	bool		negative;
	int			index;
	int			switch_index;
} ddsketch_buckets_state_t;

Datum
ddsketch_param_buckets(PG_FUNCTION_ARGS)
{
	double		alpha = PG_GETARG_FLOAT8(0);
	double		min_value = PG_GETARG_FLOAT8(1);
	double		max_value = PG_GETARG_FLOAT8(2);

	FuncCallContext *fctx;
	TupleDesc		tupdesc;

	if (SRF_IS_FIRSTCALL())
	{
		MemoryContext mctx;
		ddsketch_buckets_state_t *state;

		double	gamma = (1 + alpha) / (1 - alpha);

		double	min_indexable_value = (DBL_MIN * gamma),
				max_indexable_value = (DBL_MAX / gamma);

		fctx = SRF_FIRSTCALL_INIT();

		mctx = MemoryContextSwitchTo(fctx->multi_call_memory_ctx);

		/* Build a tuple descriptor for our result type */
		if (get_call_result_type(fcinfo, NULL, &tupdesc) != TYPEFUNC_COMPOSITE)
			elog(ERROR, "return type must be a row type");

		state = palloc(sizeof(ddsketch_buckets_state_t));
		fctx->user_fctx = state;
		fctx->max_calls = 0;

		state->tupdesc = tupdesc;

		/* Did we get a sensible range? */
		if (min_value > max_value)
			elog(ERROR, "invalid range (%e > %e)", min_value, max_value);

		/*
		 * Consider the indexable range. For the upper bound, we can't do much
		 * about those values - the ddsketch will fail anyway, so just report
		 * the issue here.
		 */
		if (fabs(min_value) > max_indexable_value)
			elog(ERROR, "maximum value is outside indexable range (%e > %e)",
				 max_value, max_indexable_value);

		if (fabs(max_value) > max_indexable_value)
			elog(ERROR, "minimum value is outside indexable range (%e > %e)",
				 max_value, max_indexable_value);

		/*
		 * For the other end of the indexable range (values close to 0), we can
		 * track such values in the zero bucket. So we just replace the value
		 * with min_indexable_value, if needed.
		 */
		if (fabs(min_value) < min_indexable_value)
			min_value = (max_value > 0) ? (min_indexable_value) : (-min_indexable_value);

		if (fabs(max_value) < min_indexable_value)
			max_value = (min_value > 0) ? (-min_indexable_value) : min_indexable_value;

		/*
		 * Now calculate the number of buckets to generate - we need to be
		 * careful about the case containing 0.
		 */
		if (((min_value > 0) && (max_value > 0)) ||
			((min_value < 0) && (max_value < 0)))
		{
			int	min_index = ddsketch_map_index2(alpha, fabs(min_value));
			int	max_index = ddsketch_map_index2(alpha, fabs(max_value));

			fctx->max_calls = (fabs(max_index - min_index) + 1);
			state->index = min_index;
			state->switch_index = (max_value < 0) ? (min_index + 1) : (min_index - 1);
			state->negative = (max_value < 0);
		}
		else
		{
			int	min_index = ddsketch_map_index2(alpha, fabs(min_value));
			int	max_index = ddsketch_map_index2(alpha, fabs(max_value));

			int	switch_index = ddsketch_map_index2(alpha, min_indexable_value);

			fctx->max_calls = (fabs(max_index - switch_index) + fabs(switch_index - min_index) + 2);
			state->index = min_index;
			state->switch_index = switch_index;
			state->negative = (min_value < 0);
		}

		MemoryContextSwitchTo(mctx);
	}

	fctx = SRF_PERCALL_SETUP();

	if (fctx->call_cntr < fctx->max_calls)
	{
		HeapTuple	resultTuple;
		Datum		result;
		Datum		values[4];
		bool		nulls[4];

		ddsketch_buckets_state_t *state = fctx->user_fctx;

		double		lower_bound = ddsketch_map_lower_bound(alpha, state->index);
		double		upper_bound = ddsketch_map_upper_bound(alpha, state->index);

		tupdesc = state->tupdesc;

		memset(nulls, 0, sizeof(nulls));

		/* Extract information from the line pointer */
		values[0] = Int32GetDatum(fctx->call_cntr);
		values[1] = Int32GetDatum(state->index);

		if (state->negative)
		{
			values[2] = Float8GetDatum(-upper_bound);
			values[3] = Float8GetDatum(-lower_bound);
		}
		else
		{
			values[2] = Float8GetDatum(lower_bound);
			values[3] = Float8GetDatum(upper_bound);
		}

		/* Build and return the result tuple. */
		resultTuple = heap_form_tuple(tupdesc, values, nulls);
		result = HeapTupleGetDatum(resultTuple);

		/* proceed */
		if (state->negative && state->index == state->switch_index)
			state->negative = false;
		else if (state->negative)
			state->index--;
		else
			state->index++;

		SRF_RETURN_NEXT(fctx, result);
	}
	else
		SRF_RETURN_DONE(fctx);
}

Datum
ddsketch_add_double_trimmed(PG_FUNCTION_ARGS)
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
		double	low = PG_GETARG_FLOAT8(4);
		double	high = PG_GETARG_FLOAT8(5);

		MemoryContext	oldcontext;

		check_trim_values(low, high);

		check_sketch_parameters(alpha, maxbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		state = ddsketch_aggstate_allocate(0, 0, alpha,
										   maxbuckets, MIN_SKETCH_BUCKETS);

		state->trim_low = low;
		state->trim_high = high;

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	ddsketch_add(state, PG_GETARG_FLOAT8(1), 1);

	PG_RETURN_POINTER(state);
}

Datum
ddsketch_add_double_count_trimmed(PG_FUNCTION_ARGS)
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
		double	low = PG_GETARG_FLOAT8(5);
		double	high = PG_GETARG_FLOAT8(6);

		MemoryContext	oldcontext;

		check_trim_values(low, high);

		check_sketch_parameters(alpha, maxbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		state = ddsketch_aggstate_allocate(0, 0, alpha,
										   maxbuckets, MIN_SKETCH_BUCKETS);

		state->trim_low = low;
		state->trim_high = high;

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
 * Add a value to the sketch (create one if needed). Transition function
 * for trimmed sketch aggregate with a single value.
 */
Datum
ddsketch_add_sketch_trimmed(PG_FUNCTION_ARGS)
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
		double	low = PG_GETARG_FLOAT8(2);
		double	high = PG_GETARG_FLOAT8(3);

		MemoryContext	oldcontext;

		check_trim_values(low, high);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		state = ddsketch_aggstate_allocate(0, 0, sketch->alpha,
										   sketch->maxbuckets, sketch->nbuckets);
		state->trim_low = low;
		state->trim_high = high;

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	AssertCheckDDSketch(sketch);
	AssertCheckDDSketchAggState(state);

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

	AssertCheckDDSketchAggState(state);

	PG_RETURN_POINTER(state);
}

/*
 * Calculate trimmed aggregates from buckets.
 */
static void
ddsketch_trimmed_agg(bucket_t *buckets, int nbuckets, int nbuckets_negative,
					double alpha, int64 count, double low, double high,
					double *sump, int64 *countp)
{
	int		i;
	double	sum = 0;
	int64	count_done = 0,
			count_low,
			count_high;

	/* translate the percentiles to counts */
	count_low = floor(count * low);
	count_high = ceil(count * high);

	count = 0;
	for (i = 0; i < nbuckets; i++)
	{
		int64	count_add = 0;
		int64	count_skip;

		double	bucket_from,
				bucket_to;

		double	start, end;

		bucket_from = ddsketch_map_lower_bound(alpha, buckets[i].index);
		bucket_to = ddsketch_map_upper_bound(alpha, buckets[i].index);

		/* How many items to skip in order to cross the lower threshold? */
		count_skip = Max(0, (count_low - count_done - 1));
		count_skip = Min(count_skip, buckets[i].count);

		/* How many items to consider including in the sum? */
		count_add = buckets[i].count - count_skip;

		Assert((count_skip >= 0) && (count_skip <= buckets[i].count));
		Assert((count_add >= 0) && (count_add <= buckets[i].count));
		Assert(count_add + count_skip == buckets[i].count);

		/*
		 * We might cross the upper threshold, ignore those too, so remove
		 * those items from the count.
		 */
		count_add -= Max(0, count_done + buckets[i].count - count_high);

		Assert((count_add >= 0) && (count_add <= buckets[i].count));
		Assert(count_add + count_skip <= buckets[i].count);

		/*
		 * Assume the values in the bucket are distributed uniformly, so
		 * make sure we include just the appropriate part of the bucket.
		 */
		start = bucket_from + (count_skip * (bucket_to - bucket_from)) / buckets[i].count;
		end = bucket_from + ((count_skip + count_add) * (bucket_to - bucket_from)) / buckets[i].count;

		/* increment the sum / count */
		sum += (start + end) / 2.0 * count_add;
		count += count_add;

		/* consider the whole bucket processed */
		count_done += buckets[i].count;

		/* break once we cross the high threshold */
		if (count_done >= count_high)
			break;
	}

	*sump = sum;
	*countp = count;
}


/*
 * Compute trimmed average from a sketch. Final function for ddsketch
 * aggregate with a low/high thresholds.
 */
Datum
ddsketch_trimmed_avg(PG_FUNCTION_ARGS)
{
	ddsketch_aggstate_t	   *state;
	MemoryContext	aggcontext;
	double			sum;
	int64			count;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_percentiles called in non-aggregate context");

	/* if there's no sketch, return NULL */
	if (PG_ARGISNULL(0))
		PG_RETURN_NULL();

	state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	ddsketch_trimmed_agg(state->buckets, state->nbuckets, state->nbuckets_negative,
						 state->alpha, state->count, state->trim_low, state->trim_high,
						 &sum, &count);

	if (count > 0)
		PG_RETURN_FLOAT8(sum / count);

	PG_RETURN_NULL();
}

/*
 * Compute trimmed sum from a sketch. Final function for ddsketch aggregate
 * with a low/high threshold.
 */
Datum
ddsketch_trimmed_sum(PG_FUNCTION_ARGS)
{
	ddsketch_aggstate_t	   *state;
	MemoryContext	aggcontext;
	double			sum;
	int64			count;

	/* cannot be called directly because of internal-type argument */
	if (!AggCheckCallContext(fcinfo, &aggcontext))
		elog(ERROR, "ddsketch_percentiles called in non-aggregate context");

	/* if there's no sketch, return NULL */
	if (PG_ARGISNULL(0))
		PG_RETURN_NULL();

	state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	ddsketch_trimmed_agg(state->buckets, state->nbuckets, state->nbuckets_negative,
						 state->alpha, state->count, state->trim_low, state->trim_high,
						 &sum, &count);

	if (count > 0)
		PG_RETURN_FLOAT8(sum);

	PG_RETURN_NULL();
}

/*
 * Trimmed sum of a single sketch (non-aggregate function).
 */
Datum
ddsketch_sketch_sum(PG_FUNCTION_ARGS)
{
	ddsketch_t *sketch = PG_GETARG_DDSKETCH(0);
	double		low = PG_GETARG_FLOAT8(1);
	double		high = PG_GETARG_FLOAT8(2);

	double		sum;
	int64		count;

	AssertCheckDDSketch(sketch);

	ddsketch_trimmed_agg(sketch->buckets, sketch->nbuckets, sketch->nbuckets_negative,
						 sketch->alpha, sketch->count, low, high, &sum, &count);

	if (count > 0)
		PG_RETURN_FLOAT8(sum);

	PG_RETURN_NULL();
}

/*
 * Trimmed average of a single sketch (non-aggregate function)
 */
Datum
ddsketch_sketch_avg(PG_FUNCTION_ARGS)
{
	ddsketch_t *sketch = PG_GETARG_DDSKETCH(0);
	double		low = PG_GETARG_FLOAT8(1);
	double		high = PG_GETARG_FLOAT8(2);

	double		sum;
	int64		count;

	AssertCheckDDSketch(sketch);

	ddsketch_trimmed_agg(sketch->buckets, sketch->nbuckets, sketch->nbuckets_negative,
						 sketch->alpha, sketch->count, low, high, &sum, &count);

	if (count > 0)
		PG_RETURN_FLOAT8(sum / count);

	PG_RETURN_NULL();
}
