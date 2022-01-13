# ddsketch extension

[![make installcheck](https://github.com/tvondra/ddsketch/actions/workflows/ci.yml/badge.svg)](https://github.com/tvondra/ddsketch/actions/workflows/ci.yml)

> :warning: **Warning**: This extension is still an early WIP version, not
> suitable for production use. The on-disk format, function signatures etc.
> may change and so on.

This PostgreSQL extension implements ddsketch, a data structure for on-line
accumulation of quantiles, as described in a paper

    DDSketch: A Fast and Fully-Mergeable Quantile Sketch with
    Relative-Error Guarantees, Charles Masson, Jee E. Rim, Homin K. Lee,
    Proceedings of the VLDB Endowment, Vol. 12, No. 12, ISSN 2150-8097
    DOI: https://doi.org/10.14778/3352063.3352135

The algorithm is very friendly to parallel programs, fully mergeable, etc.
A second paper published in 2020 introduces a variant of the sketch, using
a more elaborate procedure when collapsing buckets

    UDDSketch: Accurate Tracking of Quantiles in Data Streams; Italo
    Epicoco, Catiuscia Melle, Massimo Cafaro, Marco Pulimeno, Giuseppe
    Morleo; https://arxiv.org/abs/2004.08604

This allows providing formal accuracy guarantees even for sketches with
collapsed buckets, which is not possible for ddsketch.


## Basic usage

The extension provides two functions, which you can see as a replacement of
`percentile_cont` aggregate:

* `ddsketch_percentile(value double precision,
                       alpha double precision, nbuckets int,
                       quantile double precision)`

* `ddsketch_percentile(value double precision,
                       alpha double precision, nbuckets int,
                       quantiles double precision[])`

* `ddsketch_percentile_of(value double precision,
                          alpha double precision, nbuckets int,
                          value double precision)`

* `ddsketch_percentile_of(value double precision,
                          alpha double precision, nbuckets int,
                          values double precision[])`

That is, instead of running

```
SELECT percentile_cont(0.95) WITHIN GROUP (ORDER BY a) FROM t
```

you might now run

```
SELECT ddsketch_percentile(a, 0.05, 1024, 0.95) FROM t
```

and similarly for the variants with array of percentiles. This should run
much faster, as the ddsketch does not require any sorting of the data and
can be parallelized. Also, the memory usage is very limited, depending on
the `alpha` and `nbuckets` parameters.


## Accuracy

All functions building the ddsketch summaries accept `alpha` parameter
that determines how closely is the CDF approximated. The value limits
the size of the "buckets" in the ddsketch, so the lower the value the
larger the sketch.

Each bucket is represented by a single 8B counter, so 1000 buckets means
the ddsketch is ~8kB. That is however before the transparent compression
all varlena types go through, so the on-disk size may be much smaller.


## Advanced usage

The extension also provides a `ddsketch` data type, which makes it possible
to precompute sketches for subsets of data, and then quickly combine those
"partial" sketched into a sketch representing the whole data set. Those
prebuilt sketches should be much smaller compared to the original data set,
allowing significantly faster response times.

To compute the `ddsketch` use `ddsketch` aggregate function. The sketches can
then be stored on disk and later summarized using the `ddsketch_percentile`
functions (with `ddsketch` as the first argument).

* `ddsketch(value double precision,
            alpha double precision, nbuckets int)`

* `ddsketch_percentile(sketch ddsketch,
					   alpha double precision, nbuckets int,
                       quantile double precision)`

* `ddsketch_percentile(sketch ddsketch,
                       alpha double precision, nbuckets int,
                       quantiles double precision[])`

* `ddsketch_percentile_of(sketch ddsketch,
						  alpha double precision, nbuckets int,
                          value double precision)`

* `ddsketch_percentile_of(sketch ddsketch,
						  alpha double precision, nbuckets int,
                          values double precision[])`

So for example you may do this:

```
-- table with some random source data
CREATE TABLE t (a int, b int, c double precision);

INSERT INTO t SELECT 1000 * random(), 1000 * random(), 1000 * random()
                FROM generate_series(1,10000000);

-- table with pre-aggregated sketches into table "p"
CREATE TABLE p AS SELECT a, b, ddsketch(c, 0.05, 1024) AS d FROM t GROUP BY a, b;

-- summarize the data from "p" (compute the 95-th percentile)
SELECT a, ddsketch_percentile(d, 0.95) FROM p GROUP BY a ORDER BY a;
```

The pre-aggregated table is indeed much smaller:

~~~
db=# \d+
                                  List of relations
 Schema | Name | Type  | Owner | Persistence | Access method |  Size   | Description 
--------+------+-------+-------+-------------+---------------+---------+-------------
 public | p    | table | user  | permanent   | heap          | 2264 kB | 
 public | t    | table | user  | permanent   | heap          | 422 MB  | 
(2 rows)
~~~

And on my machine the last query takes ~1.5ms. Compare that to queries on
the source data:

~~~
\timing on

-- exact results
SELECT a, percentile_cont(0.95) WITHIN GROUP (ORDER BY c)
  FROM t GROUP BY a ORDER BY a;
  ...
Time: 6956.566 ms (00:06.957)

-- ddsketch estimate (no parallelism)
SET max_parallel_workers_per_gather = 0;
SELECT a, ddsketch_percentile(c, 0.05, 1024, 0.95) FROM t GROUP BY a ORDER BY a;
  ...
Time: 2873.116 ms (00:02.873)

-- ddsketch estimate (4 workers)
SET max_parallel_workers_per_gather = 4;
SELECT a, ddsketch_percentile(c, 0.05, 1024, 0.95) FROM t GROUP BY a ORDER BY a;
  ...
Time: 893.538 ms
~~~

This shows how much more efficient the ddsketch estimate is compared to the
exact query with `percentile_cont` (the difference would increase for larger
data sets, due to increased overhead for spilling to disk).

It also shows how effective the pre-aggregation can be. There are 121 rows
in table `p` so with 2264kB disk space that's ~20kB per row, each representing
about 80k values. With 8B per value, that's ~640kB, i.e. a compression ratio
of 30:1. As the sketch size is not tied to the number of items, this will
only improve for larger data set.


## Pre-aggregated data

When dealing with data sets with a lot of redundancy (values repeating
many times), it may be more efficient to partially pre-aggregate the data
and use functions that allow specifying the number of occurrences for each
value. This reduces the number of SQL-function calls.

There are five such aggregate functions:

* `ddsketch_percentile(value double precision, count bigint, compression int,
                      quantile double precision)`

* `ddsketch_percentile(value double precision, count bigint, compression int,
                      quantiles double precision[])`

* `ddsketch_percentile_of(value double precision, count bigint, compression int,
                         value double precision)`

* `ddsketch_percentile_of(value double precision, count bigint, compression int,
                         values double precision[])`

* `ddsketch(value double precision, count bigint, compression int)`


## Incremental updates

An existing ddsketch may be updated incrementally, either by adding a single
value, or by merging-in a whole ddsketch. For example, it's possible to add
1000 random values to the ddsketch like this:

```
DO LANGUAGE plpgsql $$
DECLARE
  r record;
BEGIN
  FOR r IN (SELECT random() AS v FROM generate_series(1,1000)) LOOP
    UPDATE t SET d = ddsketch_add(d, r.v);
  END LOOP;
END $$;
```

The overhead of doing this is fairly high, though - the ddsketch has to be
deserialized and serialized over and over, for each value we're adding.
That overhead may be reduced by pre-aggregating data, either into an array
or a ddsketch.

```
DO LANGUAGE plpgsql $$
DECLARE
  a double precision[];
BEGIN
  SELECT array_agg(random()) INTO a FROM generate_series(1,1000);
  UPDATE t SET d = ddsketch_add(d, a);
END $$;
```

Alternatively, it's possible to use pre-aggregated sketch values instead
of the arrays:

```
DO LANGUAGE plpgsql $$
DECLARE
  r record;
BEGIN
  FOR r IN (SELECT mod(i,3) AS a, ddsketch(random(), 0.05, 1024) AS d FROM generate_series(1,1000) s(i) GROUP BY mod(i,3)) LOOP
    UPDATE t SET d = ddsketch_union(d, r.d);
  END LOOP;
END $$;
```


## Trimmed aggregates

The extension provides several variants of trimmed (truncated) average and
sum aggregates, both for individual values and pre-aggregated sketches.

* `ddsketch_sum(value double precision, alpha double precision, nbuckets int, low double precision, high double precision)`

* `ddsketch_sum(value double precision, count bigint, alpha double precision, nbuckets int, low double precision, high double precision)`

* `ddsketch_sum(sketch ddsketch, low double precision, high double precision)`

* `ddsketch_avg(value double precision, alpha double precision, nbuckets int, low double precision, high double precision)`

* `ddsketch_avg(value double precision, count bigint, alpha double precision, nbuckets int, low double precision, high double precision)`

* `ddsketch_avg(sketch ddsketch, low double precision, high double precision)`

The extension also provides two regular functions allowing to calculate
trimmed (truncated) sum and average from an existing sketch.

* `ddsketch_sketch_sum(digest ddsketch, low double precision, high double precision)`

* `ddsketch_sketch_avg(digest ddsketch, low double precision, high double precision)`

The `low` and `high` parameters specify where to truncate the data. This
is useful when processing individual sketches without having to add an
unnecessary aggregation.


## Functions

### `ddsketch_percentile(value, alpha, nbuckets, percentile)`

Computes a requested percentile from the data, using a sketch with the
specified accuracy.

#### Synopsis

```
SELECT ddsketch_percentile(t.c, 0.05, 1024, 0.95) FROM t
```

#### Parameters

- `value` - values to aggregate
- `alpha` - accuracy of the sketch
- `nbuckets` - number of buckets in the sketch
- `percentile` - value in [0, 1] specifying the percentile


### `ddsketch_percentile(value, count, alpha, nbuckets, percentile)`

Computes a requested percentile from the data, using a sketch with the
specified accuracy.

#### Synopsis

```
SELECT ddsketch_percentile(t.c, t.a, 0.05, 1024, 0.95) FROM t
```

#### Parameters

- `value` - values to aggregate
- `count` - number of occurrences of the value
- `alpha` - accuracy of the sketch
- `nbuckets` - number of buckets in the sketch
- `percentile` - value in [0, 1] specifying the percentile


### `ddsketch_percentile(value, alpha, nbuckets, percentile[])`

Computes requested percentiles from the data, using a sketch with the
specified accuracy.

#### Synopsis

```
SELECT ddsketch_percentile(t.c, 0.05, 1024, ARRAY[0.95, 0.99]) FROM t
```

#### Parameters

- `value` - values to aggregate
- `alpha` - accuracy of the sketch
- `nbuckets` - number of buckets in the sketch
- `percentile[]` - array of values in [0, 1] specifying the percentiles


### `ddsketch_percentile(value, count, alpha, nbuckets, percentile[])`

Computes requested percentiles from the data, using a sketch with the
specified accuracy.

#### Synopsis

```
SELECT ddsketch_percentile(t.c, t.a, 0.05, 1024, ARRAY[0.95, 0.99]) FROM t
```

#### Parameters

- `value` - values to aggregate
- `count` - number of occurrences of the value
- `alpha` - accuracy of the sketch
- `nbuckets` - number of buckets in the sketch
- `percentile[]` - array of values in [0, 1] specifying the percentiles


### `ddsketch_percentile_of(value, alpha, nbuckets, hypothetical_value)`

Computes relative rank of a hypothetical value, using a sketch with the
specified accuracy.

#### Synopsis

```
SELECT ddsketch_percentile_of(t.c, 0.05, 1024, 139832.3) FROM t
```

#### Parameters

- `value` - values to aggregate
- `alpha` - accuracy of the sketch
- `nbuckets` - number of buckets in the sketch
- `hypothetical_value` - hypothetical value


### `ddsketch_percentile_of(value, count, alpha, nbuckets, hypothetical_value)`

Computes relative rank of a hypothetical value, using a sketch with the
specified accuracy.

#### Synopsis

```
SELECT ddsketch_percentile_of(t.c, t.a, 0.05, 1024, 139832.3) FROM t
```

#### Parameters

- `value` - values to aggregate
- `count` - number of occurrences of the value
- `alpha` - accuracy of the sketch
- `nbuckets` - number of buckets in the sketch
- `hypothetical_value` - hypothetical value


### `ddsketch_percentile_of(value, alpha, nbuckets, hypothetical_value[])`

Computes relative ranks of a hypothetical values, using a sketch with
the specified accuracy.

#### Synopsis

```
SELECT ddsketch_percentile_of(t.c, 0.05, 1024, ARRAY[6343.43, 139832.3]) FROM t
```

#### Parameters

- `value` - values to aggregate
- `alpha` - accuracy of the sketch
- `nbuckets` - number of buckets in the sketch
- `hypothetical_value` - hypothetical values


### `ddsketch_percentile_of(value, count, alpha, nbuckets, hypothetical_value[])`

Computes relative ranks of a hypothetical values, using a sketch with
the specified accuracy.

#### Synopsis

```
SELECT ddsketch_percentile_of(t.c, t.a, 0.05, 1024, ARRAY[6343.43, 139832.3]) FROM t
```

#### Parameters

- `value` - values to aggregate
- `count` - number of occurrences of the value
- `alpha` - accuracy of the sketch
- `nbuckets` - number of buckets in the sketch
- `hypothetical_value` - hypothetical values


### `ddsketch(value, alpha, nbuckets)`

Computes sketch with the specified accuracy.

#### Synopsis

```
SELECT ddsketch(t.c, 0.05, 1024) FROM t
```

#### Parameters

- `value` - values to aggregate
- `alpha` - accuracy of the sketch
- `nbuckets` - number of buckets in the sketch


### `ddsketch(value, count, alpha, nbuckets)`

Computes sketch with the specified accuracy. The values are added with
as many occurrences as determined by the count parameter.

#### Synopsis

```
SELECT ddsketch(t.c, t.a, 0.05, 1024) FROM t
```

#### Parameters

- `value` - values to aggregate
- `count` - number of occurrences for each value
- `alpha` - accuracy of the sketch
- `nbuckets` - number of buckets in the sketch


### `ddsketch_count(ddsketch)`

Returns number of items represented by the sketch.

#### Synopsis

```
SELECT ddsketch_count(d) FROM (
    SELECT ddsketch(t.c, 0.05, 1024) FROM t
) foo
```


### `ddsketch_percentile(ddsketch, percentile)`

Computes requested percentile from the pre-computed ddsketch.

#### Synopsis

```
SELECT ddsketch_percentile(d, 0.99) FROM (
    SELECT ddsketch(t.c, 0.05, 1024) FROM t
) foo
```

#### Parameters

- `sketch` - ddsketch to aggregate and process
- `percentile` - value in [0, 1] specifying the percentile


### `ddsketch_percentile(sketch, percentile[])`

Computes requested percentiles from the pre-computed ddsketch.

#### Synopsis

```
SELECT ddsketch_percentile(d, ARRAY[0.95, 0.99]) FROM (
    SELECT ddsketch(t.c, 0.05, 1024) FROM t
) foo
```

#### Parameters

- `sketch` - sketch to aggregate and process
- `percentile` - values in [0, 1] specifying the percentiles


### `ddsketch_percentile_of(sketch, hypothetical_value)`

Computes relative rank of a hypothetical value, using a pre-computed sketch.

#### Synopsis

```
SELECT ddsketch_percentile_of(d, 349834.1) FROM (
    SELECT ddsketch(t.c, 0.05, 1024) FROM t
) foo
```

#### Parameters

- `sketch` - ddsketch to aggregate and process
- `hypothetical_value` - hypothetical value


### `ddsketch_percentile_of(sketch, hypothetical_value[])`

Computes relative ranks of hypothetical values, using a pre-computed sketch.

#### Synopsis

```
SELECT ddsketch_percentile_of(d, ARRAY[438.256, 349834.1]) FROM (
    SELECT ddsketch(t.c, 0.05, 1024) FROM t
) foo
```

#### Parameters

- `sketch` - ddsketch to aggregate and process
- `hypothetical_value` - hypothetical values


### `ddsketch_add(ddsketch, double precision)`

Performs incremental update of the sketch by adding a single value.

#### Synopsis

```
UPDATE t SET d = ddsketch_add(d, random());
```

#### Parameters

- `sketch` - ddsketch to update
- `element` - value to add to the sketch
- `alpha` - accuracy of the sketch
- `nbuckets` - number of buckets in the sketch


### `ddsketch_add(ddsketch, double precision[])`

Performs incremental update of the sketch by adding values from an array.

#### Synopsis

```
UPDATE t SET d = ddsketch_add(d, ARRAY[random(), random(), random()]);
```

#### Parameters

- `sketch` - ddsketch to update
- `elements` - array of values to add to the sketch
- `alpha` - accuracy of the sketch
- `nbuckets` - number of buckets in the sketch


### `ddsketch_union(ddsketch, ddsketch)`

Performs incremental update of the sketch by merging-in another sketch.

#### Synopsis

```
WITH x AS (SELECT ddsketch(random(), 0.05, 1024) AS d FROM generate_series(1,1000))
UPDATE t SET d = ddsketch_union(t.d, x.d) FROM x;
```

#### Parameters

- `ddsketch` - ddsketch to update
- `ddsketch_add` - sketch to merge into `sketch`
- `alpha` - accuracy of the sketch
- `nbuckets` - number of buckets in the sketch


### `ddsketch_sum(value, alpha, nbuckets, low, high)`

Computes trimmed sum of values, discarding values at the low and high end.
The `low` and `high` values specify which part of the sample should be
included in the result, so e.g. `low = 0.1` and `high = 0.9` means 10% low
and high values will be discarded.

#### Synopsis

```
SELECT ddsketch_sum(t.v, 0.05, 1024, 0.1, 0.9) FROM t
```

#### Parameters

- `value` - values to aggregate
- `alpha` - accuracy of the t-digest
- `nbuckets` - number of buckets in the sketch
- `low` - low threshold percentile (values below are discarded)
- `high` - high threshold percentile (values above are discarded)


### `ddsketch_sum(value, count, alpha, nbuckets, low, high)`

Computes trimmed sum of values, discarding values at the low and high end.
The `low` and `high` values specify which part of the sample should be
included in the result, so e.g. `low = 0.1` and `high = 0.9` means 10% low
and high values will be discarded.

#### Synopsis

```
SELECT ddsketch_sum(t.v, t.c, 0.05, 1024, 0.1, 0.9) FROM t
```

#### Parameters

- `value` - values to aggregate
- `count` - number of occurrences of the value
- `alpha` - accuracy of the t-digest
- `nbuckets` - number of buckets in the sketch
- `low` - low threshold percentile (values below are discarded)
- `high` - high threshold percentile (values above are discarded)


### `ddsketch_sum(sketch, low, high)`

Computes trimmed sum of values, discarding values at the low and high end.
The `low` and `high` values specify which part of the sample should be
included in the result, so e.g. `low = 0.1` and `high = 0.9` means 10% low
and high values will be discarded.

#### Synopsis

```
SELECT ddsketch_sum(d, 0.1, 0.9) FROM (
    SELECT ddsketch(t.c, 0.05, 1024) FROM t
) foo
```

#### Parameters

- `count` - number of occurrences of the value
- `low` - low threshold percentile (values below are discarded)
- `high` - high threshold percentile (values above are discarded)


### `ddsketch_avg(value, alpha, nbuckets, low, high)`

Computes trimmed average of values, discarding values at the low and high end.
The `low` and `high` values specify which part of the sample should be
included in the result, so e.g. `low = 0.1` and `high = 0.9` means 10% low
and high values will be discarded.

#### Synopsis

```
SELECT ddsketch_avg(t.v, 0.05, 1024, 0.1, 0.9) FROM t
```

#### Parameters

- `value` - values to aggregate
- `alpha` - accuracy of the t-digest
- `nbuckets` - number of buckets in the sketch
- `low` - low threshold percentile (values below are discarded)
- `high` - high threshold percentile (values above are discarded)


### `ddsketch_avg(value, count, alpha, nbuckets, low, high)`

Computes trimmed average of values, discarding values at the low and high end.
The `low` and `high` values specify which part of the sample should be
included in the result, so e.g. `low = 0.1` and `high = 0.9` means 10% low
and high values will be discarded.

#### Synopsis

```
SELECT ddsketch_avg(t.v, 0.05, 1024, 0.1, 0.9) FROM t
```

#### Parameters

- `value` - values to aggregate
- `count` - number of occurrences of the value
- `alpha` - accuracy of the t-digest
- `nbuckets` - number of buckets in the sketch
- `low` - low threshold percentile (values below are discarded)
- `high` - high threshold percentile (values above are discarded)


### `ddsketch_avg(sketch, low, high)`

Computes trimmed average of values, discarding values at the low and high end.
The `low` and `high` values specify which part of the sample should be
included in the result, so e.g. `low = 0.1` and `high = 0.9` means 10% low
and high values will be discarded.

#### Synopsis

```
SELECT ddsketch_avg(d, 0.1, 0.9) FROM (
    SELECT ddsketch(t.c, 0.05, 1024) FROM t
) foo
```

#### Parameters

- `sketch` - ddsketch to aggregate and process
- `low` - low threshold percentile (values below are discarded)
- `high` - high threshold percentile (values above are discarded)


### `ddsketch_sketch_sum(sketch, low, high)`

Calculates trimmed sum from a single sketch, without aggregation. The `low`
and `high` values specify which part of the sample should be included in the
result, so e.g. `low = 0.1` and `high = 0.9` means 10% low and high values
will be discarded.

#### Synopsis

```
SELECT ddsketch_sketch_sum(
    (SELECT ddsketch(t.c, 0.05, 1024) FROM t),
    0.1, 0.9)
```

#### Parameters

- `sketch` - ddsketch to calculate trimmed sum for
- `low` - low threshold percentile (values below are discarded)
- `high` - high threshold percentile (values above are discarded)


### `ddsketch_sketch_avg(sketch, low, high)`

Calculates trimmed average from a single sketch, without aggregation. The
`low` and `high` values specify which part of the sample should be included
in the result, so e.g. `low = 0.1` and `high = 0.9` means 10% low and high
values will be discarded.

#### Synopsis

```
SELECT ddsketch_sketch_avg(
    (SELECT ddsketch(t.c, 0.05, 1024) FROM t),
    0.1, 0.9)
```

#### Parameters

- `sketch` - ddsketch to calculate trimmed average for
- `low` - low threshold percentile (values below are discarded)
- `high` - high threshold percentile (values above are discarded)


Notes
-----

At the moment, the extension only supports `double precision` values, but
it should not be very difficult to extend it to other numeric types (both
integer and/or floating point, including `numeric`). Ultimately, it could
support any data type with a concept of ordering and mean.

The estimates do depend on the order of incoming data, and so may differ
between runs. This applies especially to parallel queries, for which the
workers generally see different subsets of data for each run (and build
different sketches, which are then combined together).


License
-------
This software is distributed under the terms of PostgreSQL license.
See LICENSE or http://www.opensource.org/licenses/bsd-license.php for
more details.


[1] http://www.vldb.org/pvldb/vol12/p2195-masson.pdf
