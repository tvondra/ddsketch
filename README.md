# ddsketch extension

[![make installcheck](https://github.com/tvondra/ddsketch/actions/workflows/ci.yml/badge.svg)](https://github.com/tvondra/ddsketch/actions/workflows/ci.yml)

> :warning: **Warning**: This extension is still an early WIP version, not
> suitable for production use. The on-disk format, function signatures etc.
> may change and so on.

This PostgreSQL extension implements ddsketch, a data structure for on-line
accumulation of quantiles, as described in a paper [1]

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
collapsed buckets, which is not possible for ddsketch. This extension
implements neither procedure: it never collapses buckets (see Notes).


## Contents

* [Basic usage](#basic-usage)
* [Accuracy](#accuracy)
* [Advanced usage](#advanced-usage)
* [Pre-aggregated data](#pre-aggregated-data)
* [Incremental updates](#incremental-updates)
* [Trimmed statistics](#trimmed-statistics)
* [Installation](#installation)
* [Functions](#functions)
* [Notes](#notes)
* [Security](#security)
* [License](#license)


## Basic usage

The extension provides several aggregate functions, building the `ddsketch`
sketch from source data:

* `ddsketch(value double precision, alpha double precision, max_buckets int)`

And then a couple functions calculating percentiles (or inverse values)
for a given `ddsketch` sketch:

* `ddsketch_percentile(sketch ddsketch, percentile double precision)`

* `ddsketch_percentile(sketch ddsketch, percentiles double precision[])`

* `ddsketch_percentile_of(sketch ddsketch, value double precision)`

* `ddsketch_percentile_of(sketch ddsketch, values double precision[])`

That is, instead of running

```sql
SELECT percentile_cont(0.95) WITHIN GROUP (ORDER BY a) FROM t
```

you might now run

```sql
SELECT ddsketch_percentile(ddsketch(a, 0.05, 1024), 0.95) FROM t
```

and similarly for the variants with array of percentiles. This should run
much faster, as the ddsketch does not require any sorting of the data and
can be parallelized. Also, the memory usage is very limited, depending on
the `alpha` and `max_buckets` parameters.


## Accuracy

All functions building ddsketch summaries accept an `alpha` parameter
that bounds relative error in quantile values, not error in the CDF.
It limits logarithmic bucket width, so lower values generally require
more buckets. Inverse-rank estimates have no corresponding error bound.
The `alpha` value must be between 0.0001 and 0.1.

The estimated quantile is the "lower" quantile, i.e. the value at rank
`floor(q * (n - 1))` (counting from 0) in the sorted data, and the error
bound is relative to that value. There is no interpolation between values,
so for small data sets or data with large gaps between values the estimate
may differ from `percentile_cont` or `percentile_disc` by much more than
`alpha`.

The `max_buckets` parameter is the capacity of the sketch, i.e. the maximum
number of non-empty buckets (not counting the zero bucket). Only non-empty
buckets are stored, and adding a value that would need more buckets fails
with a "bucket overflow" error. It must be between 16 and 32768.

Inputs must be finite and their magnitude must not exceed `max_indexable`
reported by `ddsketch_info(alpha)`. Values with magnitude at or below
`min_indexable` are represented by the zero bucket; the relative-error
guarantee does not apply to these values.

Each bucket stores a 32-bit index and a 64-bit counter. Including
alignment, this is typically 16 bytes per bucket, plus a fixed header
per sketch. Thus 1000 buckets use about 16kB before PostgreSQL's
transparent compression, which may reduce the on-disk size.


## Advanced usage

The extension also provides a `ddsketch` data type, which makes it possible
to precompute sketches for subsets of data, and then quickly combine those
"partial" sketches into a sketch representing the whole data set. Those
prebuilt sketches should be much smaller compared to the original data set,
allowing significantly faster response times.

To compute the `ddsketch` use `ddsketch` aggregate function. The sketches can
then be stored on disk and later summarized using the `ddsketch_percentile`
functions (with `ddsketch` as the first argument).

* `ddsketch(sketch ddsketch) -> ddsketch`

So for example you may do this:

```sql
-- table with some random source data
CREATE TABLE t (a int, b int, c double precision);

INSERT INTO t SELECT 50 * random(), 50 * random(), 1000 * random()
                FROM generate_series(1,10000000);

-- table with pre-aggregated sketches into table "p"
CREATE TABLE p AS SELECT a, b, ddsketch(c, 0.05, 1024) AS d FROM t GROUP BY a, b;

-- summarize the data from "p" (compute the 95-th percentile)
SELECT a, ddsketch_percentile(ddsketch(d), 0.95) FROM p GROUP BY a ORDER BY a;
```

The outer `ddsketch(d)` merges all `(a, b)` sketches for each `a` before
calculating its percentile. Calling `ddsketch_percentile(d, 0.95)` directly
would return a separate result for every `(a, b)` group, not a percentile
over all the original values with the same `a`.

The pre-aggregated table is indeed much smaller:

~~~
db=# \d+
                                    List of relations
 Schema |   Name   | Type  | Owner | Persistence | Access method |  Size   | Description 
--------+----------+-------+-------+-------------+---------------+---------+-------------
 public | p        | table | user  | permanent   | heap          | 3112 kB | 
 public | t        | table | user  | permanent   | heap          | 422 MB  | 
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
SELECT a, ddsketch_percentile(ddsketch(c, 0.05, 1024), 0.95) FROM t GROUP BY a ORDER BY a;
  ...
Time: 2873.116 ms (00:02.873)

-- ddsketch estimate (4 workers)
SET max_parallel_workers_per_gather = 4;
SELECT a, ddsketch_percentile(ddsketch(c, 0.05, 1024), 0.95) FROM t GROUP BY a ORDER BY a;
  ...
Time: 893.538 ms
~~~

This shows how much more efficient the ddsketch estimate is compared to the
exact query with `percentile_cont` (the difference would increase for larger
data sets, due to increased overhead for spilling to disk).

It also shows how effective the pre-aggregation can be. There are ~2600 rows
in table `p` so with 3112kB disk space that's ~1.2kB per row, each representing
about 4000 values. With 8B per value, that's ~32kB, i.e. a compression ratio
of 25:1. As the sketch size is not tied to the number of items, this will
only improve for larger data set.

Pre-aggregation trades the cost of maintaining sketches for less data to
scan at query time. Storage and runtime savings depend on the number of
groups and occupied buckets, not just the number of input rows. Small
groups may use more space as sketches than as raw values.


## Pre-aggregated data

When dealing with data sets with a lot of redundancy (values repeating
many times), it may be more efficient to partially pre-aggregate the data
and use an aggregate function that allows specifying the number of
occurrences for each value. This reduces the number of SQL-function calls.

The weighted aggregate accepts the count explicitly:

* `ddsketch(value double precision, count bigint, alpha double precision, max_buckets int)`

The count must be positive; a NULL count is treated as a single occurrence.


## Incremental updates

An existing ddsketch may be updated incrementally, either by adding a single
value, or by merging-in a whole ddsketch. For example, it's possible to add
1000 random values to the sketches in table `p` like this:

```sql
DO LANGUAGE plpgsql $$
DECLARE
  r record;
BEGIN
  FOR r IN (SELECT random() AS v FROM generate_series(1,1000)) LOOP
    UPDATE p SET d = ddsketch_add(d, r.v);
  END LOOP;
END $$;
```

The overhead of doing this is fairly high, though - the ddsketch has to be
deserialized and serialized over and over, for each value we're adding.
That overhead may be reduced by pre-aggregating data, either into an array
or a ddsketch.

```sql
DO LANGUAGE plpgsql $$
DECLARE
  vals double precision[];
BEGIN
  SELECT array_agg(random()) INTO vals FROM generate_series(1,1000);
  UPDATE p SET d = ddsketch_add(d, vals);
END $$;
```

Alternatively, it's possible to use pre-aggregated sketch values instead
of the arrays:

```sql
DO LANGUAGE plpgsql $$
DECLARE
  r record;
BEGIN
  FOR r IN (SELECT mod(i,3) AS a, ddsketch(random(), 0.05, 1024) AS d FROM generate_series(1,1000) s(i) GROUP BY mod(i,3)) LOOP
    UPDATE p SET d = ddsketch_union(d, r.d);
  END LOOP;
END $$;
```


## Trimmed statistics

The extension provides trimmed (truncated) average and sum functions for
precomputed sketches. Build a sketch with `ddsketch` before calling them.

* `ddsketch_sum(sketch ddsketch, low double precision = 0.0, high double precision = 1.0)`

* `ddsketch_avg(sketch ddsketch, low double precision = 0.0, high double precision = 1.0)`

The `low` and `high` parameters specify where to truncate the data. Both
must be in [0, 1], with `low` not greater than `high`. With the default
bounds, the functions estimate the sum and average of all values. In SQL,
the parameters are named `p_low` and `p_high`, so a single bound may be
specified using named notation, e.g. `ddsketch_avg(d, p_high => 0.9)`.

Equal bounds select no observations, so both functions return `NULL`.
For a nonempty interval, the lower rank is rounded down and the upper
rank is rounded up to whole observations.


## Installation

The extension supports PostgreSQL 11 and newer, and is built using PGXS.
That requires the server development files (e.g. the
`postgresql-server-dev-NN` package on Debian and Ubuntu). To build and
install the extension, run

```sh
make
make install
```

This uses the PostgreSQL installation with `pg_config` found in `PATH`.
A different installation may be specified using `PG_CONFIG`, e.g.
`make PG_CONFIG=/path/to/pg_config install`. Installing usually requires
root privileges (e.g. `sudo make install`).

The extension then has to be created in each database by a superuser:

```sql
CREATE EXTENSION ddsketch;
```

The extension is also available on [PGXN](https://pgxn.org/dist/ddsketch/),
and may be installed using the
[PGXN client](https://pgxn.github.io/pgxnclient/) (e.g. the `pgxnclient`
package on Debian and Ubuntu), which downloads, builds and installs it:

```sh
pgxn install [--unstable] ddsketch
```
By default, the client installs only stable releases, and the extension
is not released as stable yet (see the warning at the top), hence the
`--unstable` option. The build has the same requirements as above. The
client also uses the `pg_config` found in `PATH`, and a different
installation may be specified using `--pg_config /path/to/pg_config`. If
installing requires root privileges, add `--sudo` after the extension
name (e.g. `pgxn install --unstable ddsketch --sudo`), to run just the
installation step using `sudo`.

The regression tests are executed by `make installcheck`, against a
running server with the extension installed (the server is specified by
the usual libpq environment variables, e.g. `PGHOST` and `PGPORT`). The
tests have to connect as a superuser, and require the `lower_quantile`
extension (e.g. `pgxn install lower_quantile`), which calculates exact
results for comparison.


## Functions

### `ddsketch(value, alpha, max_buckets)`

Computes a ddsketch with the specified accuracy. NULL values are skipped,
and the result is NULL if there are no non-NULL values.

#### Synopsis

```sql
SELECT ddsketch(t.c, 0.05, 1024) FROM t
```

#### Parameters

- `value` (`double precision`) - values to aggregate
- `alpha` (`double precision`) - accuracy of the sketch
- `max_buckets` (`int`) - capacity of the sketch (maximum number of non-empty buckets)

#### Returns

- sketch (`ddsketch` type) built on the input values


### `ddsketch(value, count, alpha, max_buckets)`

Computes a ddsketch with the specified accuracy. The values are added with
as many occurrences as determined by the count parameter. As with the
previous aggregate, `alpha` and `max_buckets` should be the same for all
rows.

#### Synopsis

```sql
SELECT ddsketch(s.a, s.n, 0.05, 1024) FROM (
    SELECT t.a, count(*) AS n FROM t GROUP BY t.a
) s
```

#### Parameters

- `value` (`double precision`) - values to aggregate
- `count` (`bigint`) - positive number of occurrences of the value; NULL is treated
  as one
- `alpha` (`double precision`) - accuracy of the sketch
- `max_buckets` (`int`) - capacity of the sketch (maximum number of non-empty buckets)

#### Returns

- sketch (`ddsketch` type) built on the input values


### `ddsketch(sketch)`

Computes ddsketch by combining the input sketches. All the sketches must
use the same `alpha`, and the result uses the largest capacity (see
`ddsketch_union`).

#### Synopsis

```sql
WITH tmp AS (SELECT ddsketch(t.c, 0.05, 1024) AS d FROM t GROUP BY t.a)
SELECT ddsketch(d) FROM tmp
```

#### Parameters

- `sketch` (`ddsketch`) - sketch to merge into the result

#### Returns

- sketch (`ddsketch` type) built from the input sketches


### `ddsketch_percentile(sketch, percentile)`

Computes requested percentile from the pre-computed ddsketch.

#### Synopsis

```sql
SELECT ddsketch_percentile(d, 0.99) FROM (
    SELECT ddsketch(t.c, 0.05, 1024) AS d FROM t
) foo
```

#### Parameters

- `sketch` (`ddsketch`) - sketch to process
- `percentile` (`double precision`) - value in [0, 1] specifying the percentile

#### Return value

- percentile estimate (`double precision` type)


### `ddsketch_percentile(sketch, percentiles)`

Computes requested percentiles from the pre-computed ddsketch.

#### Synopsis

```sql
SELECT ddsketch_percentile(d, ARRAY[0.95, 0.99]) FROM (
    SELECT ddsketch(t.c, 0.05, 1024) AS d FROM t
) foo
```

#### Parameters

- `sketch` (`ddsketch`) - sketch to process
- `percentiles` (`double precision[]`) - values in [0, 1] specifying the percentiles (a non-empty
  one-dimensional array without NULLs)

#### Return value

- array of percentile estimates (`double precision[]` type)


### `ddsketch_percentile_of(sketch, value)`

Computes relative rank of a hypothetical value, using a pre-computed sketch.
The estimate is the fraction of values in lower buckets, plus half the values
in the bucket containing the hypothetical value (all of them for the zero
bucket). `-Infinity` produces 0, `Infinity` produces 1 and `NaN` produces
`NaN`. When passing these as string literals, cast them explicitly (e.g.
`'Infinity'::double precision`), as an untyped literal matches both this
function and the array variant.

#### Synopsis

```sql
SELECT ddsketch_percentile_of(d, 349834.1) FROM (
    SELECT ddsketch(t.c, 0.05, 1024) AS d FROM t
) foo
```

#### Parameters

- `sketch` (`ddsketch`) - sketch to process
- `value` (`double precision`) - hypothetical value

#### Return value

- relative rank of a value (`double precision` type)


### `ddsketch_percentile_of(sketch, values)`

Computes relative ranks of hypothetical values, using a pre-computed sketch.

#### Synopsis

```sql
SELECT ddsketch_percentile_of(d, ARRAY[438.256, 349834.1]) FROM (
    SELECT ddsketch(t.c, 0.05, 1024) AS d FROM t
) foo
```

#### Parameters

- `sketch` (`ddsketch`) - sketch to process
- `values` (`double precision[]`) - hypothetical values (a non-empty one-dimensional
  array without NULLs); in named notation, the name has to be quoted
  (`"values" => ...`), as `VALUES` is an SQL keyword

#### Return value

- array of relative rank estimates (`double precision[]` type)


### `ddsketch_count(sketch)`

Returns number of items represented by the sketch.

#### Synopsis

```sql
SELECT ddsketch_count(d) FROM (
    SELECT ddsketch(t.c, 0.05, 1024) AS d FROM t
) foo
```

#### Parameters

- `sketch` (`ddsketch`) - sketch to inspect

#### Return value

- number of elements (`bigint` type) added to the sketch


### `ddsketch_avg(sketch, low, high)`

Computes trimmed average of values, discarding values at the low and high end.
The `low` and `high` values specify which part of the sample should be
included in the result, so e.g. `low = 0.1` and `high = 0.9` means 10% low
and high values will be discarded.

#### Synopsis

```sql
SELECT ddsketch_avg(d, 0.1, 0.9) FROM (
    SELECT ddsketch(t.c, 0.05, 1024) AS d FROM t
) foo
```

#### Parameters

- `sketch` (`ddsketch`) - sketch to process
- `low` (`double precision`) - low threshold percentile (values below are discarded), 0.0 by default
- `high` (`double precision`) - high threshold percentile (values above are discarded), 1.0 by default

#### Return value

- the average value (`double precision`) estimated from the sketch


### `ddsketch_sum(sketch, low, high)`

Calculates trimmed sum from a single sketch, without aggregation. The `low`
and `high` values specify which part of the sample should be included in the
result, so e.g. `low = 0.1` and `high = 0.9` means 10% low and high values
will be discarded.

#### Synopsis

```sql
SELECT ddsketch_sum(
    (SELECT ddsketch(t.c, 0.05, 1024) FROM t),
    0.1, 0.9)
```

#### Parameters

- `sketch` (`ddsketch`) - sketch to process
- `low` (`double precision`) - low threshold percentile (values below are discarded), 0.0 by default
- `high` (`double precision`) - high threshold percentile (values above are discarded), 1.0 by default

#### Return value

- the sum value (`double precision`) estimated from the sketch


### `ddsketch_add(sketch, value, alpha, max_buckets)`

Performs incremental update of the sketch by adding a single value.

If `sketch` is NULL, a new sketch is created using `alpha` and
`max_buckets`. Otherwise both parameters are ignored (and not validated),
and the sketch keeps its own accuracy and capacity. The same applies to
the other `ddsketch_add` variants.

#### Synopsis

```sql
UPDATE p SET d = ddsketch_add(d, random());
```

#### Parameters

- `sketch` (`ddsketch`) - ddsketch to update
- `value` (`double precision`) - value to add to the sketch
- `alpha` (`double precision`) - accuracy, required only when creating a sketch from NULL
- `max_buckets` (`int`) - capacity, required only when creating a sketch from NULL

#### Return value

- a sketch (`ddsketch`) with the value added


### `ddsketch_add(sketch, value, count, alpha, max_buckets)`

Adds a value with an explicit number of occurrences.

The count has to be a `bigint`. For an `integer` argument (e.g. a plain
literal like `100`), PostgreSQL picks the overload with `alpha` instead,
so `ddsketch_add(d, 2.0, 100)` adds the value only once (or fails, if `d`
is NULL). Use a cast, as in the synopsis, or named notation
(`count => 100`).

#### Synopsis

```sql
UPDATE p SET d = ddsketch_add(d, 2.0, 100::bigint);
```

#### Parameters

- `sketch` (`ddsketch`) - sketch to update
- `value` (`double precision`) - value to add to the sketch
- `count` (`bigint`) - positive number of occurrences; NULL is treated as one
- `alpha` (`double precision`) - accuracy, required only when creating a sketch from NULL
- `max_buckets` (`int`) - capacity, required only when creating a sketch from NULL

#### Return value

- a sketch (`ddsketch`) with the value added


### `ddsketch_add(sketch, values, alpha, max_buckets)`

Performs incremental update of the sketch by adding values from an array.

#### Synopsis

```sql
UPDATE p SET d = ddsketch_add(d, ARRAY[random(), random(), random()]);
```

#### Parameters

- `sketch` (`ddsketch`) - ddsketch to update
- `values` (`double precision[]`) - array of values to add to the sketch (non-empty, one-dimensional,
  without NULLs); in named notation, the name has to be quoted
  (`"values" => ...`), as `VALUES` is an SQL keyword
- `alpha` (`double precision`) - accuracy, required only when creating a sketch from NULL
- `max_buckets` (`int`) - capacity, required only when creating a sketch from NULL

#### Return value

- a sketch (`ddsketch`) with the values added


### `ddsketch_union(sketch1, sketch2)`

Performs incremental update of the sketch by merging-in another sketch.
If one of the sketches is NULL, the other one is returned.

Both sketches must use the same `alpha`. The result uses the larger input
capacity, including sketch aggregation and parallel combination. It may
error out if the combined buckets exceed that capacity.

A later high-capacity input cannot rescue an earlier intermediate merge
that already exceeded its limit; use consistent, sufficiently large
capacities when aggregating a sequence of sketches.

#### Synopsis

```sql
WITH x AS (SELECT ddsketch(random(), 0.05, 1024) AS d FROM generate_series(1,1000))
UPDATE p SET d = ddsketch_union(p.d, x.d) FROM x;
```

#### Parameters

- `sketch1` (`ddsketch`) - first sketch to merge
- `sketch2` (`ddsketch`) - second sketch to merge

#### Return value

- a sketch (`ddsketch`) representing the two input sketches merged


### `ddsketch_info(sketch)`

Returns information about the sketch: `bytes` (size of the sketch, before
compression), `flags` (reserved, currently 0), `alpha`, `count` (number of
values), `zero_count` (number of values in the zero bucket), `max_buckets`
(capacity), `negative_buckets` and `positive_buckets` (number of non-empty
buckets), and `min_indexable` and `max_indexable` (see Accuracy).

#### Synopsis

```sql
SELECT * FROM ddsketch_info((SELECT ddsketch(t.c, 0.05, 1024) FROM t));
```

#### Parameters

- `sketch` (`ddsketch`) - sketch to inspect

#### Return value

- a row (`record`) with various information about the sketch (see above)


### `ddsketch_info(alpha)`

Returns `min_indexable` and `max_indexable` for a given `alpha`, i.e. the
range of magnitudes a sketch can index (outside the zero bucket).

#### Synopsis

```sql
SELECT * FROM ddsketch_info(0.05);
```

#### Parameters

- `alpha` (`double precision`) - accuracy for a sketch

#### Return value

- a row (`record`) with `min_indexable` and `max_indexable` fields


### `ddsketch_buckets(sketch)`

Returns one row for each non-empty bucket of the sketch, ordered by values:
`index` (row number, starting from 0), `bucket_index`, `bucket_lower` and
`bucket_upper` (bucket boundaries, negative for buckets with negative
values), `bucket_length` (width of the bucket) and `bucket_count` (number
of values in the bucket). The zero bucket is not included, its count is
reported as `zero_count` by `ddsketch_info`.

#### Synopsis

```sql
SELECT * FROM ddsketch_buckets((SELECT ddsketch(t.c, 0.05, 1024) FROM t));
```

#### Parameters

- `sketch` (`ddsketch`) - sketch to inspect

#### Return value

- a set of rows, with one row for each non-empty bucket (see above for the
  attributes of each row)


### `ddsketch_buckets(alpha, min_value, max_value)`

Returns boundaries of the buckets covering values from `min_value` to
`max_value` for a given `alpha`: `index` (row number, starting from 0),
`bucket_index`, `bucket_min` and `bucket_max`. The bounds must be finite,
within the indexable range, and `min_value` must not be greater than
`max_value`. The zero bucket is not included, so a range spanning zero
returns all buckets down to `min_indexable` on both sides (thousands or
even millions of rows, depending on `alpha`).

#### Synopsis

```sql
SELECT * FROM ddsketch_buckets(0.05, 1.0, 1000.0);
```

#### Parameters

- `alpha` (`double precision`) - accuracy for a sketch
- `min_value` (`double precision`) - minimum value added to a sketch
- `max_value` (`double precision`) - maximum value added to a sketch

#### Return value

- a set of rows, with one row for each possible bucket (see above for the
  attributes of each row)


Notes
-----

At the moment, the extension only supports `double precision` values, but
it should not be very difficult to extend it to other numeric types (both
integer and/or floating point, including `numeric`). Ultimately, it could
support any data type with a concept of ordering and mean.

With constant configuration, bucket counts and quantile estimates do not
depend on input order or parallel worker assignment. This implementation
does not collapse buckets: it reports an error if the capacity is
exceeded, subject to the intermediate merge limits described above.


Security
--------

If you believe you have found a security vulnerability in this repository,
please report it using [this form](https://github.com/tvondra/ddsketch/security/advisories/new)
[5] of this GitHub project. This creates a private communication channel
between the reporter and the maintainers.

If you are absolutely unable to or have strong reasons not to use GitHub's
vulnerability reporting workflow, please reach out to the maintainer at
[tomas@vondra.me](mailto:tomas@vondra.me).

Notes:

* The code assumes digests stored on-disk are valid and not corrupted.
  If the suspected vulnerability requires a corrupted digest, without a way
  to create such digests (using the current version), it's not a security
  issue. This is in line with general assumptions in the Postgres code.

* A valid vulnerability must not require superuser privileges. A superuser
  can do almost anything (ultimately can read/write memory) and does not
  need to bother with vulnerabilities.


License
-------
This software is distributed under the terms of the PostgreSQL license.
See LICENSE or https://www.postgresql.org/about/licence/ for
more details.


[1] http://www.vldb.org/pvldb/vol12/p2195-masson.pdf
