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
#include "access/htup_details.h"
#include "common/int.h"
#include "libpq/pqformat.h"
#include "miscadmin.h"
#include "utils/array.h"
#include "utils/builtins.h"
#include "utils/lsyscache.h"
#include "catalog/pg_type.h"
#include "funcapi.h"

#if PG_VERSION_NUM >= 120000
#include "utils/float.h"	/* float8out_internal */
#endif

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

	/* store with buckets (positive and negative) */
	int64		zero_count;		/* values close to zero */
	int32		maxbuckets;		/* maximum number of buckets */
	int32		nbuckets;		/* number of buckets (used) */
	int32		nbuckets_negative;	/* number of buckets in negative part */
	int32		nbuckets_allocated;	/* number of buckets (allocated) */

	/* buckets (negative and positive) */
	bucket_t  *buckets;
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
PG_FUNCTION_INFO_V1(ddsketch_add_double);
PG_FUNCTION_INFO_V1(ddsketch_add_double_count);

PG_FUNCTION_INFO_V1(ddsketch_add_sketch);

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

PG_FUNCTION_INFO_V1(ddsketch_sketch_sum);
PG_FUNCTION_INFO_V1(ddsketch_sketch_avg);

Datum ddsketch_add_double(PG_FUNCTION_ARGS);
Datum ddsketch_add_double_count(PG_FUNCTION_ARGS);

Datum ddsketch_add_sketch(PG_FUNCTION_ARGS);

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

Datum ddsketch_sketch_sum(PG_FUNCTION_ARGS);
Datum ddsketch_sketch_avg(PG_FUNCTION_ARGS);

static ArrayType *double_array_allocate(int nitems);
static const double *array_to_double(ArrayType *v, const char *what, int *len);

static ddsketch_aggstate_t *ddsketch_copy(ddsketch_aggstate_t *state);

/* mapping to bucket indexes etc. */
static double ddsketch_log_gamma(double multiplier, double value);
static double ddsketch_pow_gamma(double multiplier, double value);
static int    ddsketch_map_index(int offset, double multiplier, double value);
static double ddsketch_map_value(int offset, double multiplier, double gamma, double index);

#if PG_VERSION_NUM < 150000
/*
 * Thin wrappers that convert strings to exactly 64-bit integers, matching our
 * definition of int64.  (For the naming, compare that POSIX has
 * strtoimax()/strtoumax() which return intmax_t/uintmax_t.)
 *
 * XXX Backward compatibility
 */
#ifdef HAVE_LONG_INT_64
/* int64 is "long int", so strtol() returns exactly the right width */
#define strtoi64(str, endptr, base) ((int64) strtol(str, endptr, base))
#else
/* int64 is "long long int" (C99 guarantees it is at least 64 bits) */
#define strtoi64(str, endptr, base) ((int64) strtoll(str, endptr, base))
#endif

#endif	/* PG_VERSION_NUM < 150000 */

/* boundaries for relative error */
#define	MIN_SKETCH_ALPHA	0.0001
#define MAX_SKETCH_ALPHA	0.1

#define	MIN_SKETCH_BUCKETS	16
#define MAX_SKETCH_BUCKETS	32768

#if (PG_VERSION_NUM < 160000)
/*
 * repalloc0
 *		Adjust the size of a previously allocated chunk and zero out the added
 *		space.
 */
static void *
repalloc0(void *pointer, Size oldsize, Size size)
{
	void       *ret;

	/* catch wrong argument order */
	if (unlikely(oldsize > size))
		elog(ERROR, "invalid repalloc0 call: oldsize %zu, new size %zu",
			 oldsize, size);

	ret = repalloc(pointer, size);
	memset((char *) ret + oldsize, 0, (size - oldsize));

	return ret;
}
#endif

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

	Assert(count == state->count);
#endif
}

/*
 * Estimate requested quantiles from the sketch agggregate state.
 */
static double *
ddsketch_compute_quantiles(ddsketch_t *sketch,
						   int npercentiles, const double *percentiles)
{
	int			i;
	double	   *result = palloc(sizeof(double) * npercentiles);

	/* parameters used for mapping values to buckets */
	int32		offset = 0;
	double		gamma = (1 + sketch->alpha) / (1 - sketch->alpha);
	double		multiplier = log(2.0) / log1p(2 * sketch->alpha / (1 - sketch->alpha));

	AssertCheckDDSketch(sketch);

	for (i = 0; i < npercentiles; i++)
	{
		int		j;
		int		index = 0;
		int64	count = 0;
		double	goal = (percentiles[i] * (sketch->count - 1));
		bucket_t *buckets;

		/*
		 * Process the negative, zero and positive stores, in this order.
		 */
		buckets = SKETCH_BUCKETS_NEGATIVE(sketch);
		for (j = 0; j < SKETCH_BUCKETS_NEGATIVE_COUNT(sketch); j++)
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
			result[i] = -ddsketch_map_value(offset, multiplier, gamma, index);
			continue;
		}

		/* now the zero bucket */
		count += sketch->zero_count;

		/* are we done after processing the zero bucket? */
		if (count > goal)
		{
			result[i] = 0;
			continue;
		}

		/* and finally the positive store */
		buckets = SKETCH_BUCKETS_POSITIVE(sketch);
		for (j = 0; j < SKETCH_BUCKETS_POSITIVE_COUNT(sketch); j++)
		{
			count += buckets[j].count;
			index = buckets[j].index;

			if (count > goal)
				break;
		}

		Assert(count >= goal);

		result[i] = ddsketch_map_value(offset, multiplier, gamma, index);
	}

	return result;
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
static double *
ddsketch_compute_quantiles_of(ddsketch_t *sketch,
							  int nvalues, const double *values)
{
	int		i;
	double	   *result = palloc(sizeof(double) * nvalues);

	/* parameters used for mapping values to buckets */
	int32		offset = 0;
	double		gamma = (1 + sketch->alpha) / (1 - sketch->alpha);
	double		min_indexable_value = DBL_MIN * gamma;
	double		multiplier = log(2.0) / log1p(2 * sketch->alpha / (1 - sketch->alpha));

	AssertCheckDDSketch(sketch);

	for (i = 0; i < nvalues; i++)
	{
		int64	count = 0;
		double	value = values[i];

		if (value > min_indexable_value)	/* value in positive part */
		{
			int		j;
			int		index = ddsketch_map_index(offset, multiplier, value);
			bucket_t *buckets;

			/* add the whole negative part */
			buckets = SKETCH_BUCKETS_NEGATIVE(sketch);
			for (j = 0; j < SKETCH_BUCKETS_NEGATIVE_COUNT(sketch); j++)
				count += buckets[j].count;

			/* add the zero bucket */
			count += sketch->zero_count;

			/* and now add the positive part, up to the index */
			buckets = SKETCH_BUCKETS_POSITIVE(sketch);
			for (j = 0; j < SKETCH_BUCKETS_POSITIVE_COUNT(sketch); j++)
			{
				if (buckets[j].index > index)
					break;

				if (buckets[j].index < index)
					count += buckets[j].count;
				else
					count += buckets[j].count / 2;
			}
		}
		else if (value < -min_indexable_value)	/* value in negative part */
		{
			int		j;
			int		index = ddsketch_map_index(offset, multiplier, -value);
			bucket_t *buckets;

			buckets = SKETCH_BUCKETS_NEGATIVE(sketch);
			for (j = 0; j < SKETCH_BUCKETS_NEGATIVE_COUNT(sketch); j++)
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
			buckets = SKETCH_BUCKETS_NEGATIVE(sketch);
			for (j = 0; j < SKETCH_BUCKETS_NEGATIVE_COUNT(sketch); j++)
				count += buckets[j].count;

			/* add the zero bucket */
			count += sketch->zero_count;
		}

		result[i] = count / ((double) sketch->count - 1);
	}

	return result;
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
		int32		nbuckets_old = state->nbuckets_allocated;

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
		state->buckets = repalloc0(state->buckets,
								   BUCKETS_BYTES(nbuckets_old),
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

	/* make sure we're not adding bogus NaN/infinity values as centroids */
	if (!isfinite(value))
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("all values added to ddsketch must be finite")));

	/* checking the total also bounds every individual bucket count */
	if (pg_add_s64_overflow(state->count, count, &state->count))
		ereport(ERROR,
				(errcode(ERRCODE_NUMERIC_VALUE_OUT_OF_RANGE),
				 errmsg("ddsketch count overflow")));

	if (value > state->min_indexable_value)
	{
		index = ddsketch_map_index(state->offset, state->multiplier, value);
		ddsketch_store_add(state, true, index, count);
	}
	else if (value < -state->min_indexable_value)
	{
		index = ddsketch_map_index(state->offset, state->multiplier, -value);
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
ddsketch_aggstate_allocate(double alpha, int maxbuckets, int nbuckets)
{
	Size				len;
	ddsketch_aggstate_t *state;

	/*
	 * We allocate a single chunk for the struct including percentiles and
	 * buckets.
	 */
	len = MAXALIGN(sizeof(ddsketch_aggstate_t));

	state = (ddsketch_aggstate_t *) palloc0(len);
	state->alpha = alpha;

	/* initialize the bucket store */
	state->maxbuckets = maxbuckets;
	state->nbuckets_allocated = 1;

	while (state->nbuckets_allocated < nbuckets)
		state->nbuckets_allocated *= 2;

	/* don't exceed the allowed number (maxbuckets may not be power of two) */
	state->nbuckets_allocated = Min(state->nbuckets_allocated,
									state->maxbuckets);

	Assert(state->nbuckets_allocated <= state->maxbuckets);

	/* we may need to repalloc this later */
	state->buckets = palloc0(state->nbuckets_allocated * sizeof(bucket_t));

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

static void
ddsketch_aggstate_free(ddsketch_aggstate_t *state)
{
	pfree(state->buckets);
	pfree(state);
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
check_percentiles(const double *percentiles, int npercentiles)
{
	int i;

	for (i = 0; i < npercentiles; i++)
	{
		if (!((percentiles[i] >= 0.0) && (percentiles[i] <= 1.0)))
			elog(ERROR, "invalid percentile value %f, should be in [0.0, 1.0]",
				 percentiles[i]);
	}
}

/* check that the user-specified sketch parameters are valid */

static void
check_alpha(double alpha)
{
	if (!((alpha >= MIN_SKETCH_ALPHA) && (alpha <= MAX_SKETCH_ALPHA)))
			elog(ERROR, "invalid alpha value %f", alpha);
}

static void
check_sketch_parameters(double alpha, int nbuckets)
{
	check_alpha(alpha);

	if (nbuckets < MIN_SKETCH_BUCKETS || nbuckets > MAX_SKETCH_BUCKETS)
		elog(ERROR, "invalid number of buckets %d", nbuckets);
}

static void
check_trim_values(double low, double high)
{
	if (!((low >= 0.0) && (low <= 1.0)))
		elog(ERROR, "invalid low percentile value %f, should be in [0.0, 1.0]",
			 low);

	if (!((high >= 0.0) && (high <= 1.0)))
		elog(ERROR, "invalid high percentile value %f, should be in [0.0, 1.0]",
			 high);

	if (low > high)
		elog(ERROR, "invalid low/high percentile values %f/%f, should be low <= high",
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

		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, maxbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		state = ddsketch_aggstate_allocate(alpha,
										   maxbuckets, MIN_SKETCH_BUCKETS);

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

		MemoryContext	oldcontext;

		check_sketch_parameters(alpha, maxbuckets);

		oldcontext = MemoryContextSwitchTo(aggcontext);

		state = ddsketch_aggstate_allocate(alpha,
										   maxbuckets, MIN_SKETCH_BUCKETS);

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
		elog(ERROR, "invalid count value " INT64_FORMAT ", must be a positive value", count);

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
	int32		nbuckets_old;

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

	nbuckets_old = state->nbuckets_allocated;

	/* grow the number of buckets to allocate */
	while (state->nbuckets_allocated < nbuckets)
		state->nbuckets_allocated *= 2;

	/* don't exceed the allowed number (maxbuckets may not be power of two) */
	state->nbuckets_allocated = Min(state->nbuckets_allocated,
									state->maxbuckets);

	state->buckets = repalloc0(state->buckets,
							   BUCKETS_BYTES(nbuckets_old),
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

	/* free the working array */
	pfree(b);
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

	sketch = PG_GETARG_DDSKETCH(1);

	/* if there's no aggregate state allocated, create it now */
	if (PG_ARGISNULL(0))
	{
		MemoryContext	oldcontext;

		oldcontext = MemoryContextSwitchTo(aggcontext);

		state = ddsketch_aggstate_allocate(sketch->alpha,
										   sketch->maxbuckets, sketch->nbuckets);

		MemoryContextSwitchTo(oldcontext);
	}
	else
		state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	AssertCheckDDSketch(sketch);
	AssertCheckDDSketchAggState(state);

	/* check that the sketch and aggstate are compatible */
	if (state->alpha != sketch->alpha)
		elog(ERROR, "can't merge sketches with different alpha values");

	/*
	 * XXX Should we compare the maxbuckets too? We reject values that don't
	 * fit into the sketch, so one sketch might contain values the other would
	 * have rejected - that doesn't seem great. Or we should at least pick the
	 * maxbuckets in some consistent way, instead of picking the value from
	 * the first sketch.
	 */

	/* checking the total also bounds every individual bucket count */
	if (pg_add_s64_overflow(state->count, sketch->count, &state->count))
		ereport(ERROR,
				(errcode(ERRCODE_NUMERIC_VALUE_OUT_OF_RANGE),
				 errmsg("ddsketch count overflow")));

	state->zero_count += sketch->zero_count;

	ddsketch_merge_buckets(state, false,
						   SKETCH_BUCKETS_NEGATIVE(sketch),
						   SKETCH_BUCKETS_NEGATIVE_COUNT(sketch));

	ddsketch_merge_buckets(state, true,
						   SKETCH_BUCKETS_POSITIVE(sketch),
						   SKETCH_BUCKETS_POSITIVE_COUNT(sketch));

	AssertCheckDDSketchAggState(state);

	PG_FREE_IF_COPY(sketch, 1);

	PG_RETURN_POINTER(state);
}

/*
 * Compute percentile from a ddsketch. Final function for ddsketch aggregate
 * with a single percentile.
 */
Datum
ddsketch_percentiles(PG_FUNCTION_ARGS)
{
	ddsketch_t	   *sketch = PG_GETARG_DDSKETCH(0);
	double		   *ret;
	double			percentile = PG_GETARG_FLOAT8(1);

	check_percentiles(&percentile, 1);

	ret = ddsketch_compute_quantiles(sketch, 1, &percentile);

	PG_RETURN_FLOAT8(*ret);
}

/*
 * Compute percentile from a ddsketch. Final function for ddsketch aggregate
 * with a single percentile.
 */
Datum
ddsketch_percentiles_of(PG_FUNCTION_ARGS)
{
	ddsketch_t	   *sketch = PG_GETARG_DDSKETCH(0);
	double		   *ret;
	double			value = PG_GETARG_FLOAT8(1);

	ret = ddsketch_compute_quantiles_of(sketch, 1, &value);

	PG_RETURN_FLOAT8(*ret);
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
	ddsketch_t *sketch;
	const double	   *percentiles;
	int			npercentiles;
	ArrayType	   *array;

	sketch = PG_GETARG_DDSKETCH(0);

	array = PG_GETARG_ARRAYTYPE_P(1);
	percentiles = array_to_double(array,
								  "a percentile value", &npercentiles);

	check_percentiles(percentiles, npercentiles);

	result = ddsketch_compute_quantiles(sketch, npercentiles, percentiles);

	PG_FREE_IF_COPY(sketch, 0);
	PG_FREE_IF_COPY(array, 1);

	/* copy the results into the array */
	array = double_array_allocate(npercentiles);
	memcpy((double *) ARR_DATA_PTR(array), result, sizeof(double) * npercentiles);

	/* free the result */
	pfree(result);

	PG_RETURN_ARRAYTYPE_P(array);
}

/*
 * Compute percentiles from a ddsketch. Final function for ddsketch aggregate
 * with an array of values.
 */
Datum
ddsketch_array_percentiles_of(PG_FUNCTION_ARGS)
{
	double	*result;
	const double	   *values;
	int			nvalues;
	ArrayType	   *array;

	ddsketch_t *sketch;

	array = PG_GETARG_ARRAYTYPE_P(1);
	values = array_to_double(array,
							 "a value", &nvalues);

	sketch = PG_GETARG_DDSKETCH(0);

	result = ddsketch_compute_quantiles_of(sketch, nvalues, values);

	PG_FREE_IF_COPY(sketch, 0);
	PG_FREE_IF_COPY(array, 1);

	/* copy the results into the array */
	array = double_array_allocate(nvalues);
	memcpy((double *) ARR_DATA_PTR(array), result, sizeof(double) * nvalues);

	/* free the result */
	pfree(result);

	PG_RETURN_ARRAYTYPE_P(array);
}

Datum
ddsketch_serial(PG_FUNCTION_ARGS)
{
	bytea	   *v;
	ddsketch_aggstate_t  *state;
	Size		len;
	char	   *ptr;

	state = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	len = offsetof(ddsketch_aggstate_t, buckets) +
		  STATE_BUCKETS_BYTES(state);

	v = palloc(len + VARHDRSZ);

	SET_VARSIZE(v, len + VARHDRSZ);
	ptr = VARDATA(v);

	memcpy(ptr, state, offsetof(ddsketch_aggstate_t, buckets));
	ptr += offsetof(ddsketch_aggstate_t, buckets);

	/* FIXME maybe don't serialize full buckets, but just the count */
	memcpy(ptr, STATE_BUCKETS(state), STATE_BUCKETS_BYTES(state));
	ptr += STATE_BUCKETS_BYTES(state);

	Assert(VARDATA(v) + len == ptr);

	PG_RETURN_POINTER(v);
}

/*
 * XXX Unlike the other "input" functions (ddsketch_in/ddsketch_recv), this
 * does not validate the sketch at all. We assume this function is used only
 * on data we created in the same query (possibly in a parallel worker),
 * and not on untrusted values controlled by the user (which is why the other
 * input functions need the validation).
 */
Datum
ddsketch_deserial(PG_FUNCTION_ARGS)
{
	bytea  *v = (bytea *) PG_GETARG_POINTER(0);
	char   *ptr = VARDATA_ANY(v);
	char   *endptr PG_USED_FOR_ASSERTS_ONLY;
	ddsketch_aggstate_t	tmp;
	ddsketch_aggstate_t *state;
	int		nbuckets_allocated;

	endptr = ptr + VARSIZE_ANY_EXHDR(v);

	/* copy aggstate header into a local variable */
	memcpy(&tmp, ptr, offsetof(ddsketch_aggstate_t, buckets));
	ptr += offsetof(ddsketch_aggstate_t, buckets);

	state = ddsketch_aggstate_allocate(tmp.alpha,
									   tmp.maxbuckets, tmp.nbuckets);

	 /*
	 * Copy the header, but keep the number of buckets we actually allocated
	 * above - the serialized value describes the allocation of the state it
	 * was produced from, which may well have been larger.
	*/
	nbuckets_allocated = state->nbuckets_allocated;
	memcpy(state, &tmp, offsetof(ddsketch_aggstate_t, buckets));
	state->nbuckets_allocated = nbuckets_allocated;

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
	int					nbuckets_allocated;

	AssertCheckDDSketchAggState(state);

	copy = ddsketch_aggstate_allocate(state->alpha, state->maxbuckets,
									  state->nbuckets);

	 /*
	 * Copy the header, but keep the number of buckets we actually allocated
	 * above - the serialized value describes the allocation of the state it
	 * was produced from, which may well have been larger.
	*/
	nbuckets_allocated = copy->nbuckets_allocated;
	memcpy(copy, state, offsetof(ddsketch_aggstate_t, buckets));
	copy->nbuckets_allocated = nbuckets_allocated;

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

	/* if no "merged" state yet, try creating it */
	if (PG_ARGISNULL(0))
	{
		/* nope, the second argument is NULL too, so return NULL */
		if (PG_ARGISNULL(1))
			PG_RETURN_NULL();

		/* the second argument is not NULL, so copy it */
		src = (ddsketch_aggstate_t *) PG_GETARG_POINTER(1);

		/* copy the sketch into the right long-lived memory context */
		oldcontext = MemoryContextSwitchTo(aggcontext);
		src = ddsketch_copy(src);
		MemoryContextSwitchTo(oldcontext);

		PG_RETURN_POINTER(src);
	}

	/*
	 * If the second argument is NULL, just return the first one (we know
	 * it's not NULL at this point).
	 */
	if (PG_ARGISNULL(1))
		PG_RETURN_DATUM(PG_GETARG_DATUM(0));

	src = (ddsketch_aggstate_t *) PG_GETARG_POINTER(1);
	dst = (ddsketch_aggstate_t *) PG_GETARG_POINTER(0);

	AssertCheckDDSketchAggState(src);
	AssertCheckDDSketchAggState(dst);

	/* check that the two sketches are compatible */
	if (src->alpha != dst->alpha)
		elog(ERROR, "can't merge sketches with different alpha values");

	/*
	 * XXX Should we compare the maxbuckets too? We reject values that don't
	 * fit into the sketch, so one sketch might contain values the other would
	 * have rejected - that doesn't seem great. Or we should at least pick the
	 * maxbuckets in some consistent way, instead of picking the value from
	 * the first sketch.
	 */

	/* checking the total also bounds every individual bucket count */
	if (pg_add_s64_overflow(dst->count, src->count, &dst->count))
		ereport(ERROR,
				(errcode(ERRCODE_NUMERIC_VALUE_OUT_OF_RANGE),
				 errmsg("ddsketch count overflow")));

	dst->zero_count += dst->zero_count;

	ddsketch_merge_buckets(dst, false,
						   STATE_BUCKETS_NEGATIVE(src),
						   STATE_BUCKETS_NEGATIVE_COUNT(src));

	ddsketch_merge_buckets(dst, true,
						   STATE_BUCKETS_POSITIVE(src),
						   STATE_BUCKETS_POSITIVE_COUNT(src));

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

	state = ddsketch_aggstate_allocate(sketch->alpha,
									   sketch->maxbuckets, sketch->nbuckets);

	state->count = sketch->count;
	state->zero_count = sketch->zero_count;

	state->nbuckets = sketch->nbuckets;
	state->nbuckets_negative = sketch->nbuckets_negative;

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
	ddsketch_t *sketch;
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

		state = ddsketch_aggstate_allocate(alpha, maxbuckets, MIN_SKETCH_BUCKETS);
	}
	else
	{
		sketch = PG_GETARG_DDSKETCH(0);
		state = ddsketch_sketch_to_aggstate(sketch);

		PG_FREE_IF_COPY(sketch, 0);
	}

	AssertCheckDDSketchAggState(state);

	ddsketch_add(state, PG_GETARG_FLOAT8(1), 1);

	AssertCheckDDSketchAggState(state);

	sketch = ddsketch_aggstate_to_ddsketch(state);
	ddsketch_aggstate_free(state);

	PG_RETURN_POINTER(sketch);
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
	ddsketch_t *sketch;

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

		state = ddsketch_aggstate_allocate(alpha,
										   maxbuckets,
										   MIN_SKETCH_BUCKETS);
	}
	else
	{
		sketch = PG_GETARG_DDSKETCH(0);
		state = ddsketch_sketch_to_aggstate(sketch);

		PG_FREE_IF_COPY(sketch, 0);
	}

	if (PG_ARGISNULL(2))
		count = 1;
	else
		count = PG_GETARG_INT64(2);

	AssertCheckDDSketchAggState(state);

	/* can't add values with non-positive counts */
	if (count <= 0)
		elog(ERROR, "invalid count value " INT64_FORMAT ", must be a positive value", count);

	ddsketch_add(state, PG_GETARG_FLOAT8(1), count);

	AssertCheckDDSketchAggState(state);

	sketch = ddsketch_aggstate_to_ddsketch(state);
	ddsketch_aggstate_free(state);

	PG_RETURN_POINTER(sketch);
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
	ddsketch_t		  *sketch;
	ddsketch_aggstate_t *state;
	const double	   *values;
	int					nvalues;
	int					i;
	ArrayType		   *array;

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

		state = ddsketch_aggstate_allocate(alpha,
										   maxbuckets, MIN_SKETCH_BUCKETS);
	}
	else
	{
		sketch = PG_GETARG_DDSKETCH(0);
		state = ddsketch_sketch_to_aggstate(sketch);

		PG_FREE_IF_COPY(sketch, 0);
	}

	array = PG_GETARG_ARRAYTYPE_P(1);
	values = array_to_double(array,
							 "an element", &nvalues);

	for (i = 0; i < nvalues; i++)
		ddsketch_add(state, values[i], 1);

	sketch = ddsketch_aggstate_to_ddsketch(state);
	ddsketch_aggstate_free(state);

	PG_RETURN_POINTER(sketch);
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
	sketch = PG_GETARG_DDSKETCH(0);
	state = ddsketch_sketch_to_aggstate(sketch);

	PG_FREE_IF_COPY(sketch, 0);

	/* parse the second ddsketch */
	sketch = PG_GETARG_DDSKETCH(1);

	AssertCheckDDSketch(sketch);
	AssertCheckDDSketchAggState(state);

	/* check that the two sketches are compatible */
	if (sketch->alpha != state->alpha)
		elog(ERROR, "can't merge sketches with different alpha values");

	/*
	 * XXX Should we compare the maxbuckets too? We reject values that don't
	 * fit into the sketch, so one sketch might contain values the other would
	 * have rejected - that doesn't seem great. Or we should at least pick the
	 * maxbuckets in some consistent way, instead of picking the value from
	 * the first sketch.
	 */

	/* checking the total also bounds every individual bucket count */
	if (pg_add_s64_overflow(state->count, sketch->count, &state->count))
		ereport(ERROR,
				(errcode(ERRCODE_NUMERIC_VALUE_OUT_OF_RANGE),
				 errmsg("ddsketch count overflow")));

	state->zero_count += sketch->zero_count;

	/* copy data from sketch to aggstate */
	ddsketch_merge_buckets(state, false,
						   SKETCH_BUCKETS_NEGATIVE(sketch),
						   SKETCH_BUCKETS_NEGATIVE_COUNT(sketch));

	ddsketch_merge_buckets(state, true,
						   SKETCH_BUCKETS_POSITIVE(sketch),
						   SKETCH_BUCKETS_POSITIVE_COUNT(sketch));

	AssertCheckDDSketchAggState(state);

	PG_FREE_IF_COPY(sketch, 1);

	sketch = ddsketch_aggstate_to_ddsketch(state);
	ddsketch_aggstate_free(state);

	PG_RETURN_POINTER(sketch);
}

/*
 * Parsing of the textual ddsketch representation.
 *
 * We can't use sscanf, because it does not report overflows in any way - the
 * value simply saturates to the maximum for the data type. That's a problem
 * for the count, where the saturated value is a perfectly valid count, so we
 * can't detect it after the fact. Use strtoll/strtod, which do set errno.
 *
 * All of these advance the pointer past the parsed part on success, and never
 * return on failure.
 */

/*
 * Match a literal string, after skipping (optional) leading space.
 */
static void
parse_str(char **ptr, const char *value, bool space)
{
	char   *str = *ptr;
	size_t	len = strlen(value);

	/* if requested, skip the one initial space character */
	if (space)
	{
		if (isspace((unsigned char) *str))
			str++;
		else
			ereport(ERROR,
					(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
					 errmsg("failed to parse ddsketch value, missing space")));
	}

	/* at this point there must be no whitespace */
	if (isspace((unsigned char) *str))
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("failed to parse ddsketch value, unexpected space")));

	/* the prefix should match our string */
	if (strncmp(str, value, len) != 0)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("failed to parse ddsketch value, expected \"%s\"",
						value)));

	*ptr = str + len;
}

/*
 * Parse an int64 value, and make sure it's in range.
 */
static int64
parse_int64(char **ptr, const char *field)
{
	char   *endptr;
	int64	value;

	errno = 0;
	value = strtoi64(*ptr, &endptr, 10);

	if (endptr == *ptr)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("failed to parse %s of a ddsketch", field)));

	if (errno == ERANGE)
		ereport(ERROR,
				(errcode(ERRCODE_NUMERIC_VALUE_OUT_OF_RANGE),
				 errmsg("%s of a ddsketch is out of range for bigint", field)));

	*ptr = endptr;

	return value;
}

/*
 * Parse an int32 value, and make sure it's in range.
 *
 * Parse it as int64 first, so that we can range check it before narrowing it
 * down, instead of relying on an implementation-defined narrowing conversion.
 */
static int32
parse_int32(char **ptr, const char *field)
{
	int64	value = parse_int64(ptr, field);

	if (value < PG_INT32_MIN || value > PG_INT32_MAX)
		ereport(ERROR,
				(errcode(ERRCODE_NUMERIC_VALUE_OUT_OF_RANGE),
				 errmsg("%s of a ddsketch is out of range for integer", field)));

	return (int32) value;
}

/*
 * Parse a double value, and make sure it's in range.
 */
static double
parse_double(char **ptr, const char *field)
{
	char   *endptr;
	double	value;

	errno = 0;
	value = strtod(*ptr, &endptr);

	if (endptr == *ptr)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("failed to parse %s of a ddsketch", field)));

	/* ERANGE may also signal a nonzero subnormal, which we can store. */
	if ((errno == ERANGE) && ((value == 0.0) || !isfinite(value)))
		ereport(ERROR,
				(errcode(ERRCODE_NUMERIC_VALUE_OUT_OF_RANGE),
				 errmsg("%s of a ddsketch is out of range for double precision",
						field)));

	*ptr = endptr;

	return value;
}

Datum
ddsketch_in(PG_FUNCTION_ARGS)
{
	char	   *str = PG_GETARG_CSTRING(0);
	ddsketch_t  *sketch = NULL;
	size_t		slen;

	/* ddsketch header fields */
	int32       flags;
	int64		count;
	int64		zero_count;
	double		alpha;
	int			maxbuckets;
	int			nbuckets;
	int			nbuckets_negative;
	char	   *ptr;

	slen = strlen(str);
	ptr = str;

	parse_str(&ptr, "flags", false);
	flags = parse_int32(&ptr, "flags");

	parse_str(&ptr, "count", true);
	count = parse_int64(&ptr, "count");

	parse_str(&ptr, "alpha", true);
	alpha = parse_double(&ptr, "alpha");

	parse_str(&ptr, "zero_count", true);
	zero_count = parse_int64(&ptr, "zero_count");

	parse_str(&ptr, "maxbuckets", true);
	maxbuckets = parse_int32(&ptr, "maxbuckets");

	parse_str(&ptr, "buckets", true);
	nbuckets = parse_int32(&ptr, "buckets");
	nbuckets_negative = parse_int32(&ptr, "negative buckets");

	if (flags != SKETCH_DEFAULT_FLAGS)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("invalid sketch flags %d", flags)));

	if (!((alpha >= MIN_SKETCH_ALPHA) && (alpha <= MAX_SKETCH_ALPHA)))
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("alpha for ddsketch (%f) must be in [%f, %f]",
						alpha, MIN_SKETCH_ALPHA, MAX_SKETCH_ALPHA)));

	if ((maxbuckets < MIN_SKETCH_BUCKETS) || (maxbuckets > MAX_SKETCH_BUCKETS))
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("number of buckets (%d) for ddsketch must be in [%d, %d]",
						maxbuckets, MIN_SKETCH_BUCKETS, MAX_SKETCH_BUCKETS)));

	if (nbuckets < 0)
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

	count = zero_count;

	nbuckets = 0;
	for (int i = 0; i < sketch->nbuckets; i++)
	{
		int		index;
		int64	bucket_count;

		CHECK_FOR_INTERRUPTS();

		parse_str(&ptr, "(", true);
		index = parse_int32(&ptr, "index");
		parse_str(&ptr, ", ", false);
		bucket_count = parse_int64(&ptr, "bucket count");
		parse_str(&ptr, ")", false);

		/* we've parsed a bucket, but we have too many already */
		if (nbuckets >= sketch->nbuckets)
			elog(ERROR, "too many buckets parsed");

		/*
		 * Basic checks that the indexes are decreasing in the negative part
		 * and increasing in the positive part.
		 *
		 * XXX Can we check the index value is valid (not too low/high)?
		 */
		if ((nbuckets != 0) && (nbuckets < nbuckets_negative))
		{
			/* negative store - descending index values */
			if (sketch->buckets[nbuckets - 1].index <= index)
				ereport(ERROR,
						(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						 errmsg("invalid sketch - ascending indexes in the negative part")));
		}
		else if (nbuckets  > nbuckets_negative)
		{
			/* positive store - ascending index values */
			if (sketch->buckets[nbuckets - 1].index >= index)
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
		nbuckets++;

		/*
		 * track the total count so that we can check later
		 *
		 * Make sure the count does not overflow at any point. It could
		 * overflow and then wrap around to the expected total, but it would
		 * still cause an issue.
		 */
		if (pg_add_s64_overflow(count, bucket_count, &count))
			ereport(ERROR,
					(errcode(ERRCODE_NUMERIC_VALUE_OUT_OF_RANGE),
					 errmsg("ddsketch count overflow")));

		/* must not scan past the end of the input string */
		Assert(ptr <= str + slen);
	}

	/*
	 * Malformed inputs may have the wrong number of buckets, in which case
	 * we either don't consume the whole input (nbuckets too high), or we
	 * don't get all the expected buckets (nbuckets too high).
	 */
	if (ptr < str + slen)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("input ddsketch value too long")));

	/* Did we parse exactly the expected number of buckets? */
	if (nbuckets != sketch->nbuckets)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("parsed invalid number of buckets (%d != %d)",
						nbuckets, sketch->nbuckets)));

	/*
	 * If we consumed just the right number of buckets, we must have read
	 * the whole input value exactly.
	 */
	Assert(ptr == str + strlen(str));

	/* Did we get buckets matching the header? */
	if (count != sketch->count)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("total count (" INT64_FORMAT ") does not match buckets (" INT64_FORMAT ")",
						sketch->count, count)));

	AssertCheckDDSketch(sketch);

	PG_RETURN_POINTER(sketch);
}

Datum
ddsketch_out(PG_FUNCTION_ARGS)
{
	int			i;
	ddsketch_t  *sketch = PG_GETARG_DDSKETCH(0);
	StringInfoData	str;
	char	    *alpha = float8out_internal(sketch->alpha);

	AssertCheckDDSketch(sketch);

	initStringInfo(&str);

	appendStringInfo(&str, "flags %d count " INT64_FORMAT " alpha %s zero_count " INT64_FORMAT " maxbuckets %d buckets %d %d",
					 sketch->flags, sketch->count, alpha, sketch->zero_count,
					 sketch->maxbuckets, sketch->nbuckets, sketch->nbuckets_negative);

	for (i = 0; i < sketch->nbuckets; i++)
		appendStringInfo(&str, " (%d, " INT64_FORMAT ")", sketch->buckets[i].index, sketch->buckets[i].count);

	PG_FREE_IF_COPY(sketch, 0);
	pfree(alpha);

	PG_RETURN_CSTRING(str.data);
}

Datum
ddsketch_recv(PG_FUNCTION_ARGS)
{
	StringInfo	buf = (StringInfo) PG_GETARG_POINTER(0);
	ddsketch_t  *sketch;
	int			i;
	int64		count;
	int64		total_count;
	int64		zero_count;
	int32		flags;
	int32		maxbuckets;
	int32		nbuckets;
	int32		nbuckets_negative;
	double		alpha;

	flags = pq_getmsgint(buf, sizeof(int32));

	count = pq_getmsgint64(buf);
	zero_count = pq_getmsgint64(buf);
	alpha = pq_getmsgfloat8(buf);
	maxbuckets = pq_getmsgint(buf, sizeof(int32));
	nbuckets = pq_getmsgint(buf, sizeof(int32));
	nbuckets_negative = pq_getmsgint(buf, sizeof(int32));

	if (flags != SKETCH_DEFAULT_FLAGS)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("invalid sketch flags %d", flags)));

	if (!((alpha >= MIN_SKETCH_ALPHA) && (alpha <= MAX_SKETCH_ALPHA)))
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("alpha for ddsketch (%f) must be in [%f, %f]",
						alpha, MIN_SKETCH_ALPHA, MAX_SKETCH_ALPHA)));

	if ((maxbuckets < MIN_SKETCH_BUCKETS) || (maxbuckets > MAX_SKETCH_BUCKETS))
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("number of buckets (%d) for ddsketch must be in [%d, %d]",
						maxbuckets, MIN_SKETCH_BUCKETS, MAX_SKETCH_BUCKETS)));

	if (nbuckets < 0)
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

	total_count = zero_count;
	for (i = 0; i < sketch->nbuckets; i++)
	{
		CHECK_FOR_INTERRUPTS();

		if (i >= sketch->nbuckets)
			ereport(ERROR,
					(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
					 errmsg("too many buckets parsed")));

		sketch->buckets[i].index = pq_getmsgint(buf, sizeof(int32));
		sketch->buckets[i].count = pq_getmsgint64(buf);

		/*
		 * track the total count so that we can check later
		 *
		 * Make sure the count does not overflow at any point. It could
		 * overflow and then wrap around to the expected total, but it would
		 * still cause an issue.
		 */
		if (pg_add_s64_overflow(total_count, sketch->buckets[i].count,
								&total_count))
			ereport(ERROR,
					(errcode(ERRCODE_NUMERIC_VALUE_OUT_OF_RANGE),
					 errmsg("ddsketch count overflow")));

		if (sketch->buckets[i].count <= 0)
			ereport(ERROR,
					(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
					 errmsg("count value for all buckets in a ddsketch must be positive")));
		else if (sketch->buckets[i].count > sketch->count)
			ereport(ERROR,
					(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
					 errmsg("count value of a bucket exceeds total count")));

		/* negative part, sorted by index in descending order */
		if ((i > 0) && (i < sketch->nbuckets_negative))
		{
			if (sketch->buckets[i-1].index <= sketch->buckets[i].index)
				ereport(ERROR,
						(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						 errmsg("negative buckets in incorrect order")));
		}
		else if (i > sketch->nbuckets_negative)	/* positive part */
		{
			if (sketch->buckets[i-1].index >= sketch->buckets[i].index)
				ereport(ERROR,
						(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
						 errmsg("positive buckets in incorrect order")));
		}
	}

	/* check that the total matches */
	if (total_count != sketch->count)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("total count does not match the data (%lld != %lld)",
						(long long) total_count, (long long) sketch->count)));

	AssertCheckDDSketch(sketch);

	PG_RETURN_POINTER(sketch);
}

Datum
ddsketch_send(PG_FUNCTION_ARGS)
{
	ddsketch_t  *sketch = PG_GETARG_DDSKETCH(0);
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

	PG_FREE_IF_COPY(sketch, 0);

	PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

Datum
ddsketch_count(PG_FUNCTION_ARGS)
{
	ddsketch_t  *sketch = PG_GETARG_DDSKETCH(0);
	int64	count = sketch->count;

	PG_FREE_IF_COPY(sketch, 0);

	PG_RETURN_INT64(count);
}

/*
 * Return a read-only view of an input FLOAT8 SQL array as C doubles.
 *
 * This expects a single-dimensional float8 array, fails otherwise.
 *
 * "what" names a single element of the array, with the article, and is used
 * when reporting a NULL element. The callers pass arrays of different things
 * (percentiles, hypothetical values, values to add to a digest), and a message
 * naming the wrong one points at the wrong argument.
 *
 * The caller must keep the detoasted array alive while using the view.
 */
static const double *
array_to_double(ArrayType *v, const char *what, int *len)
{
	int		nitems,
		   *dims,
			ndims;
	Oid		element_type;

	ndims = ARR_NDIM(v);
	dims = ARR_DIMS(v);
	nitems = ArrayGetNItems(ndims, dims);

	/*
	 * Reject empty arrays explicitly. An empty array has ndims = 0, so
	 * without this it would be caught by the single-dimension check below
	 * and reported as a dimensionality problem, which is misleading - the
	 * array is well-formed, it just has nothing in it.
	 */
	if (nitems == 0)
		elog(ERROR, "the array must not be empty");

	/* this is a special-purpose function for single-dimensional arrays */
	if (ndims != 1)
		elog(ERROR, "expected a single-dimensional array (dims = %d)", ndims);

	element_type = ARR_ELEMTYPE(v);

	/* XXX not sure if really needed (can it actually happen?) */
	if (element_type != FLOAT8OID)
		elog(ERROR, "array_to_double expects FLOAT8 array");

	if (array_contains_nulls(v))
		elog(ERROR, "NULL not allowed as %s", what);

	(*len) = nitems;

	/* Non-NULL float8 elements have the same layout as a C double array. */
	return (const double *) ARR_DATA_PTR(v);
}

/*
 * Allocate a one-dimensional, non-NULL float8 array for direct result writes.
 * Array storage uses native doubles even on platforms with pass-by-reference
 * float8 Datums, so no per-element Datum allocations are needed.
 */
static ArrayType *
double_array_allocate(int nitems)
{
	ArrayType  *array;
	Size		size;

	/* should not happen */
	if (nitems <= 0)
		ereport(ERROR,
				(errcode(ERRCODE_INVALID_PARAMETER_VALUE),
				 errmsg("invalid array size (%d)", nitems)));

	if (nitems > MaxArraySize)
		ereport(ERROR,
				(errcode(ERRCODE_PROGRAM_LIMIT_EXCEEDED),
				 errmsg("array size exceeds the maximum allowed (%d)",
						(int) MaxArraySize)));

	size = ARR_OVERHEAD_NONULLS(1) + nitems * sizeof(double);

	array = (ArrayType *) palloc(size);
	SET_VARSIZE(array, size);
	ARR_NDIM(array) = 1;
	array->dataoffset = 0;
	ARR_ELEMTYPE(array) = FLOAT8OID;
	ARR_DIMS(array)[0] = nitems;
	ARR_LBOUND(array)[0] = 1;

	return array;
}

static double
ddsketch_log_gamma(double multiplier, double value)
{
	return log(value) / log(2.0) * multiplier;
}

static double
ddsketch_pow_gamma(double multiplier, double value)
{
	return pow(2.0, (value / multiplier));
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
ddsketch_map_index(int offset, double multiplier, double value)
{
	return (int) (ceil(ddsketch_log_gamma(multiplier, value)) + offset);
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
ddsketch_map_value(int offset, double multiplier, double gamma, double index)
{
	return ddsketch_pow_gamma(multiplier, index - offset) * (2.0 / (1 + gamma));
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
	values[1] = Int32GetDatum(sketch->flags);
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
		values[5] = Int64GetDatum(bucket->count);

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

	check_alpha(alpha);

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

		double	gamma,
				min_indexable_value,
				max_indexable_value;

		check_alpha(alpha);

		/* now that we know alpha is OK, calculate the other parameters */
		gamma = (1 + alpha) / (1 - alpha);
		min_indexable_value = (DBL_MIN * gamma);
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

			fctx->max_calls = (abs(max_index - min_index) + 1);
			state->index = min_index;
			state->switch_index = (max_value < 0) ? (min_index + 1) : (min_index - 1);
			state->negative = (max_value < 0);
		}
		else
		{
			int	min_index = ddsketch_map_index2(alpha, fabs(min_value));
			int	max_index = ddsketch_map_index2(alpha, fabs(max_value));

			int	switch_index = ddsketch_map_index2(alpha, min_indexable_value);

			fctx->max_calls = (abs(max_index - switch_index) + abs(switch_index - min_index) + 2);
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

	check_trim_values(low, high);

	ddsketch_trimmed_agg(sketch->buckets, sketch->nbuckets, sketch->nbuckets_negative,
						 sketch->alpha, sketch->count, low, high, &sum, &count);

	PG_FREE_IF_COPY(sketch, 0);

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

	check_trim_values(low, high);

	ddsketch_trimmed_agg(sketch->buckets, sketch->nbuckets, sketch->nbuckets_negative,
						 sketch->alpha, sketch->count, low, high, &sum, &count);

	PG_FREE_IF_COPY(sketch, 0);

	if (count > 0)
		PG_RETURN_FLOAT8(sum / count);

	PG_RETURN_NULL();
}
