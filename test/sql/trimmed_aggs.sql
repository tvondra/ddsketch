CREATE TABLE random_data (v double precision, i int);
INSERT INTO random_data SELECT prng(10000, 45547, 34471541, 3, 1000000), generate_series(1,10000);

-- trimmed aggregates
SELECT ddsketch_sum(ddsketch(1000 * v, 0.05, 1024), 0.0, 1.0)   BETWEEN 4750000 AND 5250000 FROM random_data;
SELECT ddsketch_sum(ddsketch(1000 * v, 0.05, 1024), 0.0, 0.5)   BETWEEN 1200000 AND 1300000 FROM random_data;
SELECT ddsketch_sum(ddsketch(1000 * v, 0.05, 1024), 0.5, 1.0)   BETWEEN 3700000 AND 3800000 FROM random_data;
SELECT ddsketch_sum(ddsketch(1000 * v, 0.05, 1024), 0.25, 0.75) BETWEEN 2450000 AND 2550000 FROM random_data;

SELECT ddsketch_avg(ddsketch(1000 * v, 0.05, 1024), 0.0, 1.0)   BETWEEN 490 AND 510 FROM random_data;
SELECT ddsketch_avg(ddsketch(1000 * v, 0.05, 1024), 0.0, 0.5)   BETWEEN 240 AND 260 FROM random_data;
SELECT ddsketch_avg(ddsketch(1000 * v, 0.05, 1024), 0.5, 1.0)   BETWEEN 740 AND 760 FROM random_data;
SELECT ddsketch_avg(ddsketch(1000 * v, 0.05, 1024), 0.25, 0.75) BETWEEN 490 AND 510 FROM random_data;

-- trimmed aggregates, <value, count> API
SELECT ddsketch_sum(ddsketch(1000 * v, 1 + mod(i,5), 0.05, 1024), 0.0, 1.0)   BETWEEN 3 * 4750000 AND 3 * 5250000 FROM random_data;
SELECT ddsketch_sum(ddsketch(1000 * v, 1 + mod(i,5), 0.05, 1024), 0.0, 0.5)   BETWEEN 3 * 1200000 AND 3 * 1300000 FROM random_data;
SELECT ddsketch_sum(ddsketch(1000 * v, 1 + mod(i,5), 0.05, 1024), 0.5, 1.0)   BETWEEN 3 * 3700000 AND 3 * 3800000 FROM random_data;
SELECT ddsketch_sum(ddsketch(1000 * v, 1 + mod(i,5), 0.05, 1024), 0.25, 0.75) BETWEEN 3 * 2450000 AND 3 * 2550000 FROM random_data;

SELECT ddsketch_avg(ddsketch(1000 * v, 1 + mod(i,5), 0.05, 1024), 0.0, 1.0)   BETWEEN 490 AND 510 FROM random_data;
SELECT ddsketch_avg(ddsketch(1000 * v, 1 + mod(i,5), 0.05, 1024), 0.0, 0.5)   BETWEEN 240 AND 260 FROM random_data;
SELECT ddsketch_avg(ddsketch(1000 * v, 1 + mod(i,5), 0.05, 1024), 0.5, 1.0)   BETWEEN 740 AND 760 FROM random_data;
SELECT ddsketch_avg(ddsketch(1000 * v, 1 + mod(i,5), 0.05, 1024), 0.25, 0.75) BETWEEN 490 AND 510 FROM random_data;

-- trimmed aggregates when aggregating sketches
WITH sketches AS (SELECT mod(i,5) AS c, ddsketch(1000 * v, 0.05, 1024) AS s FROM random_data GROUP BY 1)
SELECT ddsketch_sum(ddsketch(s), 0.0, 1.0)   BETWEEN 4750000 AND 5250000 FROM sketches;

WITH sketches AS (SELECT mod(i,5) AS c, ddsketch(1000 * v, 0.05, 1024) AS s FROM random_data GROUP BY 1)
SELECT ddsketch_sum(ddsketch(s), 0.0, 0.5)   BETWEEN 1200000 AND 1300000 FROM sketches;

WITH sketches AS (SELECT mod(i,5) AS c, ddsketch(1000 * v, 0.05, 1024) AS s FROM random_data GROUP BY 1)
SELECT ddsketch_sum(ddsketch(s), 0.5, 1.0)   BETWEEN 3700000 AND 3800000 FROM sketches;

WITH sketches AS (SELECT mod(i,5) AS c, ddsketch(1000 * v, 0.05, 1024) AS s FROM random_data GROUP BY 1)
SELECT ddsketch_sum(ddsketch(s), 0.25, 0.75) BETWEEN 2450000 AND 2550000 FROM sketches;

WITH sketches AS (SELECT mod(i,5) AS c, ddsketch(1000 * v, 0.05, 1024) AS s FROM random_data GROUP BY 1)
SELECT ddsketch_avg(ddsketch(s), 0.0, 1.0)   BETWEEN 490 AND 510 FROM sketches;

WITH sketches AS (SELECT mod(i,5) AS c, ddsketch(1000 * v, 0.05, 1024) AS s FROM random_data GROUP BY 1)
SELECT ddsketch_avg(ddsketch(s), 0.0, 0.5)   BETWEEN 240 AND 260 FROM sketches;

WITH sketches AS (SELECT mod(i,5) AS c, ddsketch(1000 * v, 0.05, 1024) AS s FROM random_data GROUP BY 1)
SELECT ddsketch_avg(ddsketch(s), 0.5, 1.0)   BETWEEN 740 AND 760 FROM sketches;

WITH sketches AS (SELECT mod(i,5) AS c, ddsketch(1000 * v, 0.05, 1024) AS s FROM random_data GROUP BY 1)
SELECT ddsketch_avg(ddsketch(s), 0.25, 0.75) BETWEEN 490 AND 510 FROM sketches;

-- trimmed aggregates for a single sketch
SELECT ddsketch_sum(ddsketch(1000 * v, 0.05, 1024), 0.0, 1.0)   BETWEEN 4750000 AND 5250000 FROM random_data;
SELECT ddsketch_sum(ddsketch(1000 * v, 0.05, 1024), 0.0, 0.5)   BETWEEN 1200000 AND 1300000 FROM random_data;
SELECT ddsketch_sum(ddsketch(1000 * v, 0.05, 1024), 0.5, 1.0)   BETWEEN 3700000 AND 3800000 FROM random_data;
SELECT ddsketch_sum(ddsketch(1000 * v, 0.05, 1024), 0.25, 0.75) BETWEEN 2450000 AND 2550000 FROM random_data;

SELECT ddsketch_avg(ddsketch(1000 * v, 0.05, 1024), 0.0, 1.0)   BETWEEN 490 AND 510 FROM random_data;
SELECT ddsketch_avg(ddsketch(1000 * v, 0.05, 1024), 0.0, 0.5)   BETWEEN 240 AND 260 FROM random_data;
SELECT ddsketch_avg(ddsketch(1000 * v, 0.05, 1024), 0.5, 1.0)   BETWEEN 740 AND 760 FROM random_data;
SELECT ddsketch_avg(ddsketch(1000 * v, 0.05, 1024), 0.25, 0.75) BETWEEN 490 AND 510 FROM random_data;

-- check trim parameters
SELECT ddsketch_sum(ddsketch(1000 * v, 0.05, 1024), -0.1, 1.0)   BETWEEN 4750000 AND 5250000 FROM random_data;
SELECT ddsketch_sum(ddsketch(1000 * v, 0.05, 1024), 0.1, 1.1)   BETWEEN 4750000 AND 5250000 FROM random_data;
SELECT ddsketch_sum(ddsketch(1000 * v, 0.05, 1024), 0.9, 0.1)   BETWEEN 4750000 AND 5250000 FROM random_data;
SELECT ddsketch_sum(ddsketch(1000 * v, 0.05, 1024), 0.5, 0.5) IS NULL FROM random_data;

-- Negative buckets must contribute negative sums and averages.
WITH sketch AS (
    SELECT ddsketch(v, 0.05, 16) AS s FROM (VALUES (-2.0), (-1.0)) t(v)
)
SELECT ddsketch_sum(s) BETWEEN -3.15 AND -2.85 AS negative_sum,
       ddsketch_avg(s) BETWEEN -1.575 AND -1.425 AS negative_avg
FROM sketch;

WITH sketch AS (
    SELECT ddsketch(v, 0.05, 16) AS s
    FROM (VALUES (-2.0), (-1.0), (1.0), (2.0)) t(v)
)
SELECT abs(ddsketch_sum(s)) < 1e-12 AS symmetric_sum,
       ddsketch_sum(s, 0.0, 0.5) BETWEEN -3.15 AND -2.85 AS negative_half
FROM sketch;

-- An all-zero sketch has non-NULL zero sum and average.
SELECT ddsketch_sum(ddsketch(0.0, 3::bigint, 0.05, 16)) AS zero_sum,
       ddsketch_avg(ddsketch(0.0, 3::bigint, 0.05, 16)) AS zero_avg;

-- Zeros count toward both the average divisor and trimmed rank positions.
WITH sketch AS (
    SELECT ddsketch(v, 0.05, 16) AS s FROM (VALUES (0.0), (0.0), (0.0), (2.0)) t(v)
)
SELECT ddsketch_sum(s) / ddsketch_avg(s) = ddsketch_count(s) AS counts_zeros,
       ddsketch_sum(s, 0.0, 0.5) AS zero_prefix_sum,
       ddsketch_avg(s, 0.0, 0.5) AS zero_prefix_avg
FROM sketch;

-- overflows close to INT64_MAX - The results should be about the same as for a
-- sketch with much lower counts, but we can't expect the results to match exactly
-- due to rounding errors, so we use check_relative_erro.

-- two positive buckets, counts adding up to INT64_MAX-1
SELECT check_relative_error(
    (SELECT ddsketch_avg(ddsketch(v, n, 0.05, 16), 0.0, 1.0)
       FROM (VALUES (1.0::double precision, 4611686018427387903::bigint),
             (2.0::double precision, 4611686018427387903::bigint)) t(v, n)),
    (SELECT ddsketch_avg(ddsketch(v, n, 0.05, 16), 0.0, 1.0)
       FROM (VALUES (1.0::double precision, 1::bigint),
             (2.0::double precision, 1::bigint)) t(v, n)),
    1e-9) AS avg_matches;

-- negative, zero and positive buckets, counts adding up to INT64_MAX-1
SELECT check_relative_error(
    (SELECT ddsketch_avg(ddsketch(v, n, 0.05, 16), 0.0, 1.0)
       FROM (VALUES (-2.0::double precision, 3074457345618258602::bigint),
             ( 0.0::double precision, 3074457345618258602::bigint),
             ( 1.0::double precision, 3074457345618258602::bigint)) t(v, n)),
    (SELECT ddsketch_avg(ddsketch(v, n, 0.05, 16), 0.0, 1.0)
       FROM (VALUES (-2.0::double precision, 1::bigint),
             ( 0.0::double precision, 1::bigint),
             ( 1.0::double precision, 1::bigint)) t(v, n)),
    1e-9) AS avg_matches;

-- a single bucket with exactly INT64_MAX items
SELECT check_relative_error(
    (SELECT ddsketch_avg(ddsketch(1.0, 9223372036854775807::bigint, 0.05, 16), 0.0, 1.0)),
    (SELECT ddsketch_avg(ddsketch(1.0, 1::bigint, 0.05, 16), 0.0, 1.0)),
    1e-9) AS avg_matches;

-- partially trimmed ranges have to be well behaved too - with the same number
-- of items in each of the three buckets, trimming the bottom and the top third
-- has to leave just the middle (zero) bucket
CREATE TEMP TABLE trimmed_overflow_sketch AS
SELECT ddsketch(v, n, 0.05, 16) AS s
  FROM (VALUES (-2.0::double precision, 3074457345618258602::bigint),
        ( 0.0::double precision, 3074457345618258602::bigint),
        ( 1.0::double precision, 3074457345618258602::bigint)) t(v, n);

SELECT ddsketch_count(s) AS total_count FROM trimmed_overflow_sketch;

-- With the counts this high, the 1/3 and 2/3 will not align perfectly with
-- the buckets, but will be rounded to get a couple items from the first and
-- third bucket, so don't expect exact result.
--
-- Note: The 1e3 tolerance is needed, because double arithmetics has ULP > 1.0,
-- and in this case (3074457345618258602 * 3) / 3.0 = 3074457345618258432.
SELECT abs(ddsketch_sum(s, 1.0/3.0, 2.0/3.0)) < 1e3 AS middle_is_zero,
       abs(ddsketch_avg(s, 1.0/3.0, 2.0/3.0)) < 1e-9 AS middle_avg_is_zero
  FROM trimmed_overflow_sketch;

-- the whole range must not be empty, and the bottom/top parts must have the
-- expected signs
SELECT ddsketch_sum(s, 0.0, 1.0) < 0 AS total_negative,
       ddsketch_sum(s, 0.0, 1.0/3.0) < 0 AS bottom_negative,
       ddsketch_sum(s, 2.0/3.0, 1.0) > 0 AS top_positive
  FROM trimmed_overflow_sketch;

-- empty ranges are still empty
SELECT ddsketch_sum(s, 0.0, 0.0) AS empty_low,
       ddsketch_avg(s, 1.0, 1.0) AS empty_high
  FROM trimmed_overflow_sketch;

DROP TABLE trimmed_overflow_sketch;
