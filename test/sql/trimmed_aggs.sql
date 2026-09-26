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
