-- hypothetical-set aggregates
--
-- there's no guarantee for relative-errors of hypothetical-aggregates,
-- but for uniform distribution it's fairly close to the relative error
-- the ddsketch was defined with
WITH
  data AS (SELECT i / 10.0 AS v FROM generate_series(1,2500) s(i)), 
  sketch AS (SELECT ddsketch(data.v, 0.05, 1024) AS s FROM data)
SELECT
  v,
  check_relative_error(a, b, 0.06) AS check_error,
  print_relative_error(a, b, 0.06) AS error_info
FROM (
  SELECT
    foo.v,
    (SELECT ddsketch_percentile_of(sketch.s, foo.v) FROM sketch) a,
    (SELECT percent_rank(foo.v) WITHIN GROUP (ORDER BY v) FROM data) b
  FROM
    (SELECT i AS v FROM generate_series(0,1000,25) s(i)) foo
) bar;

WITH
  data AS (SELECT mod(i,10) as x, i / 10.0 AS v FROM generate_series(1,10000) s(i)), 
  sketches AS (SELECT ddsketch(data.v, 0.05, 1024) AS s FROM data GROUP BY data.x)
SELECT
  v,
  check_relative_error(a, b, 0.06) AS check_error,
  print_relative_error(a, b, 0.06) AS error_info
FROM (
  SELECT
    foo.v,
    (SELECT ddsketch_percentile_of(ddsketch(sketches.s), foo.v) FROM sketches) a,
    (SELECT percent_rank(foo.v) WITHIN GROUP (ORDER BY v) FROM data) b
  FROM
    (SELECT i AS v FROM generate_series(0,1000,25) s(i)) foo
) bar;

WITH
  data AS (SELECT i / 10.0 AS v FROM generate_series(1,10000) s(i)), 
  sketch AS (SELECT ddsketch(data.v, 0.05, 1024) AS s FROM data),
  vals AS (SELECT array_agg(i) AS v FROM generate_series(0,1000,25) s(i)),
  sketch_values AS (SELECT ddsketch_percentile_of(sketch.s, vals.v) AS v FROM sketch, vals),
  percent_ranks AS (SELECT array_agg((SELECT percent_rank(i) WITHIN GROUP (ORDER BY v) FROM data)) AS v FROM generate_series(0,1000,25) s(i))
SELECT
  v,
  check_relative_error(a, b, 0.06) AS check_error,
  print_relative_error(a, b, 0.06) AS error_info
FROM (
  SELECT
    unnest(vals.v) as v,
    unnest(sketch_values.v) as a,
    unnest(percent_ranks.v) as b
  FROM
    vals, sketch_values, percent_ranks
) bar;

WITH
  data AS (SELECT mod(i,10) AS x, i / 10.0 AS v FROM generate_series(1,10000) s(i)), 
  sketches AS (SELECT ddsketch(data.v, 0.05, 1024) AS s FROM data GROUP BY x),
  vals AS (SELECT array_agg(i) AS v FROM generate_series(0,1000,25) s(i)),
  sketch_values AS (SELECT ddsketch_percentile_of(ddsketch(sketches.s), vals.v) AS v FROM sketches, vals GROUP BY vals.v),
  percent_ranks AS (SELECT array_agg((SELECT percent_rank(i) WITHIN GROUP (ORDER BY v) FROM data)) AS v FROM generate_series(0,1000,25) s(i))
SELECT
  v,
  check_relative_error(a, b, 0.06) AS check_error,
  print_relative_error(a, b, 0.06) AS error_info
FROM (
  SELECT
    unnest(vals.v) as v,
    unnest(sketch_values.v) as a,
    unnest(percent_ranks.v) as b
  FROM
    vals, sketch_values, percent_ranks
) bar;

-- now the same thing, but pass the values as a single array
WITH
  data AS (SELECT i / 10.0 AS v FROM generate_series(1,10000) s(i)), 
  vals AS (SELECT array_agg(i) AS v FROM generate_series(0,1000,25) s(i)),
  sketch_values AS (SELECT ddsketch_percentile_of(ddsketch(data.v, 0.05, 1024), vals.v) AS v FROM data, vals GROUP BY vals.v),
  percent_ranks AS (SELECT array_agg((SELECT percent_rank(i) WITHIN GROUP (ORDER BY v) FROM data)) AS v FROM generate_series(0,1000,25) s(i))
SELECT
  v,
  check_relative_error(a, b, 0.06) AS check_error,
  print_relative_error(a, b, 0.06) AS error_info
FROM (
  SELECT
    unnest(vals.v) as v,
    unnest(sketch_values.v) as a,
    unnest(percent_ranks.v) as b
  FROM
    vals, sketch_values, percent_ranks
) bar;

WITH
  data AS (SELECT i / 10.0 AS v FROM generate_series(1,10000) s(i)), 
  vals AS (SELECT array_agg(i) AS v FROM generate_series(0,1000,25) s(i)),
  sketch_values AS (SELECT ddsketch_percentile_of(ddsketch(data.v, 0.05, 1024), vals.v) AS v FROM data, vals GROUP BY vals.v),
  percent_ranks AS (SELECT array_agg((SELECT percent_rank(i) WITHIN GROUP (ORDER BY v) FROM data)) AS v FROM generate_series(0,1000,25) s(i))
SELECT
  v,
  check_relative_error(a, b, 0.06) AS check_error,
  print_relative_error(a, b, 0.06) AS error_info
FROM (
  SELECT
    unnest(vals.v) as v,
    unnest(sketch_values.v) as a,
    unnest(percent_ranks.v) as b
  FROM
    vals, sketch_values, percent_ranks
) bar;
