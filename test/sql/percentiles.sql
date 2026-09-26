-- percentile_of with an array of values
WITH data AS (SELECT i / 50.0 AS v FROM generate_series(-5000, 5000) s(i))
SELECT
    f,
    abs(a - b) < 0.05
FROM (
    SELECT
      unnest(ARRAY[-100.0, -75.0, -50.0, -25.0, 0.0, 25.0, 50.0, 75.0, 100.0]) f,
      unnest((SELECT ddsketch_percentile_of(ddsketch(data.v, 0.05, 1024), ARRAY[-100.0, -75.0, -50.0, -25.0, 0.0, 25.0, 50.0, 75.0, 100.0]) FROM data)) a,
      unnest((SELECT array_agg((SELECT percent_rank(f) WITHIN GROUP (ORDER BY v) FROM data)) foo FROM unnest(ARRAY[-100.0, -75.0, -50.0, -25.0, 0.0, 25.0, 50.0, 75.0, 100.0]) f)) AS b
) foo;

-- percentile_of and individual values
WITH
    data AS (SELECT i / 50.0 AS v FROM generate_series(-5000, 5000) s(i))
SELECT
    f,
    abs(a - b) < 0.05
FROM (
    SELECT
        f,
        (SELECT ddsketch_percentile_of(ddsketch(data.v, 0.05, 1024), f) FROM data) a,
        (SELECT percent_rank(f) WITHIN GROUP (ORDER BY v) FROM data) AS b
    FROM unnest(ARRAY[-100.0, -75.0, -50.0, -25.0, 0, 25.0, 50.0, 75.0, 100.0]) AS f
) foo;

-- <value,count> API with percentile_of and individual values
WITH
    data AS (SELECT i / 50.0 AS v, 1 + abs(mod(i,13)) AS c FROM generate_series(-500, 500) s(i)),
    data_expanded AS (SELECT foo.v FROM (SELECT data.c, data.v FROM data) foo, LATERAL generate_series(1, c))
SELECT
    f,
    abs(a - b) < 0.05
FROM (
    SELECT
      f,
      (SELECT ddsketch_percentile_of(ddsketch(data.v, data.c, 0.05, 1024), f) FROM data) a,
      (SELECT percent_rank(f) WITHIN GROUP (ORDER BY v) FROM data_expanded) AS b
    FROM unnest(ARRAY[-100.0, -75.0, -50.0, -25.0, 0, 25.0, 50.0, 75.0, 100.0]) AS f
) foo;

-- <value,count> API with percentile_of and an array
WITH
    data AS (SELECT i / 50.0 AS v, 1 + abs(mod(i,13)) AS c FROM generate_series(-500, 500) s(i)),
    data_expanded AS (SELECT foo.v FROM (SELECT data.c, data.v FROM data) foo, LATERAL generate_series(1, c))
SELECT
    f,
    abs(a - b) < 0.05
FROM (
    SELECT
      unnest(ARRAY[-100.0, -75.0, -50.0, -25.0, 0.0, 25.0, 50.0, 75.0, 100.0]) f,
      unnest((SELECT ddsketch_percentile_of(ddsketch(data.v, data.c, 0.05, 1024), ARRAY[-100.0, -75.0, -50.0, -25.0, 0.0, 25.0, 50.0, 75.0, 100.0]) FROM data)) a,
      unnest((SELECT array_agg((SELECT percent_rank(f) WITHIN GROUP (ORDER BY v) FROM data_expanded)) foo FROM unnest(ARRAY[-100.0, -75.0, -50.0, -25.0, 0.0, 25.0, 50.0, 75.0, 100.0]) f)) AS b
) foo;

-- Special probes have defined results independent of bucket mapping.
WITH sketch AS (
    SELECT ddsketch(v, 0.05, 16) AS s
    FROM (VALUES (-1.0), (0.0), (1.0)) t(v)
)
SELECT ddsketch_percentile_of(s, 'NaN'::float8) AS nan,
       ddsketch_percentile_of(s, '-Infinity'::float8) AS below,
       ddsketch_percentile_of(s, 'Infinity'::float8) AS above,
       ddsketch_percentile_of(s, ARRAY['NaN', '-Infinity', 'Infinity']::float8[]) AS probes
FROM sketch;
