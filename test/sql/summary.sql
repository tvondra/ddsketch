-- verify we can store ddsketch in a summary table
CREATE TABLE intermediate_ddsketch (grouping int, summary ddsketch);

WITH data AS (SELECT row_number() OVER () AS i, 1.0 + pow(z, 4) AS x FROM random_normal(10000) s(z))
INSERT INTO intermediate_ddsketch
SELECT
    i % 10 AS grouping,
    ddsketch(x, 0.05, 1024) AS summary
FROM data
GROUP BY i % 10;

WITH data AS (SELECT 1.0 + pow(z, 4) AS x FROM random_normal(10000) s(z)),
     intermediate AS (SELECT ddsketch_percentile(ddsketch(summary), ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS a FROM intermediate_ddsketch),
     pg_percentile AS (SELECT lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS b FROM data)
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(a) AS a,
        unnest(b) AS b
    FROM intermediate,
         pg_percentile
) foo;

-- verify we can store ddsketch in a summary table (individual percentiles)
TRUNCATE intermediate_ddsketch;

WITH data AS (SELECT row_number() OVER () AS i, 100 * pow(z, 4) AS x FROM random_normal(10000) s(z))
INSERT INTO intermediate_ddsketch
SELECT
    i % 10 AS grouping,
    ddsketch(x, 0.05, 1024) AS summary
FROM data
GROUP BY i % 10;

WITH data AS (SELECT 100 * pow(z, 4) AS x FROM random_normal(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        foo.p AS p,
        (SELECT ddsketch_percentile(ddsketch(summary), foo.p) FROM intermediate_ddsketch) AS a,
        (SELECT lower_quantile(x, p) AS b FROM data) AS b
    FROM
         (SELECT unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) p) foo
) bar;

-- verify 'extreme' percentiles for the dataset would not read out of bounds buckets
WITH data AS (SELECT x FROM generate_series(1,10) AS x)
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.99]) AS p,
        unnest(ddsketch_percentile(ddsketch(x, 0.05, 1024), ARRAY[0.01, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.99])) AS b
    FROM data
) foo;

-- check that the computed percentiles are perfectly correlated (don't decrease for higher p values)
WITH
-- percentiles to compute
perc AS (SELECT array_agg((i / 100.0)::double precision) AS percentiles FROM generate_series(1,99) s(i)),
-- input data (just 15 points)
input_data AS (select i::double precision AS val FROM generate_series(1,15) s(i))
SELECT * FROM (
    SELECT p, v AS v1, lag(v, 1) OVER (ORDER BY p) v2 FROM (
        SELECT
            unnest(perc.percentiles) p,
            unnest(ddsketch_percentile(ddsketch(input_data.val, 0.05, 1024), perc.percentiles)) v
        FROM perc, input_data
        GROUP BY perc.percentiles
    ) foo
) bar where v2 > v1;
