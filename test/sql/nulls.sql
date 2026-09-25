-- NULL handling

-- individual values, individual percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text))
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
    p,
    (SELECT ddsketch_percentile(ddsketch(data.v, 0.01, 1024), p) FROM data) AS a,
    (SELECT lower_quantile(v, p) FROM data) AS b
  FROM
    unnest((SELECT p FROM perc)) p
) foo;

-- <value,count> API, individual percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text)),
  data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c))
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
    p,
    (SELECT ddsketch_percentile(ddsketch(data.v, data.c, 0.01, 1024), p) FROM data) AS a,
    (SELECT lower_quantile(v, p) FROM data_exp) AS b
  FROM
    unnest((SELECT p FROM perc)) p
) foo;

-- <value,count> API, but count is NULL (should be treated as 1), individual percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text))
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
    p,
    (SELECT ddsketch_percentile(ddsketch(data.v, NULL, 0.01, 1024), p) FROM data) AS a,
    (SELECT lower_quantile(v, p) FROM data) AS b
  FROM
    unnest((SELECT p FROM perc)) p
) foo;

-- individual values, array of percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text))
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
  SELECT
    unnest((SELECT p FROM perc)) p,
    unnest((SELECT ddsketch_percentile(ddsketch(data.v, 0.01, 1024), perc.p) FROM data, perc GROUP BY perc.p)) AS a,
    unnest((SELECT lower_quantile(v, p) FROM data, perc)) AS b
) foo;

-- <value,count> API, array of percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text)),
  data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c))
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
  SELECT
    unnest((SELECT p FROM perc)) p,
    unnest((SELECT ddsketch_percentile(ddsketch(data.v, data.c, 0.01, 1024), perc.p) FROM data, perc GROUP BY perc.p)) AS a,
    unnest((SELECT lower_quantile(v, p) FROM data_exp, perc)) AS b
) foo;

-- <value,count> API, but count is NULL (should be treated as 1), array of percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text))
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
  SELECT
    unnest((SELECT p FROM perc)) p,
    unnest((SELECT ddsketch_percentile(ddsketch(data.v, NULL, 0.01, 1024), perc.p) FROM data, perc GROUP BY perc.p)) AS a,
    unnest((SELECT lower_quantile(v, p) FROM data, perc)) AS b
) foo;

-- NULL handling, but this time make sure the first value is NULL

-- individual values, individual percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT NULL AS v, 1 AS c UNION ALL SELECT * FROM (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text)) foo)
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
    p,
    (SELECT ddsketch_percentile(ddsketch(data.v, 0.01, 1024), p) FROM data) AS a,
    (SELECT lower_quantile(v, p) FROM data) AS b
  FROM
    unnest((SELECT p FROM perc)) p
) foo;

-- <value,count> API, individual percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT NULL AS v, 1 AS c UNION ALL SELECT * FROM (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text)) foo),
  data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c))
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
    p,
    (SELECT ddsketch_percentile(ddsketch(data.v, data.c, 0.01, 1024), p) FROM data) AS a,
    (SELECT lower_quantile(v, p) FROM data_exp) AS b
  FROM
    unnest((SELECT p FROM perc)) p
) foo;

-- <value,count> API, but count is NULL (should be treated as 1), individual percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT NULL AS v, 1 AS c UNION ALL SELECT * FROM (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text)) foo)
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
    p,
    (SELECT ddsketch_percentile(ddsketch(data.v, NULL, 0.01, 1024), p) FROM data) AS a,
    (SELECT lower_quantile(v, p) FROM data) AS b
  FROM
    unnest((SELECT p FROM perc)) p
) foo;

-- individual values, array of percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT NULL AS v, 1 AS c UNION ALL SELECT * FROM (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text)) foo)
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
  SELECT
    unnest((SELECT p FROM perc)) p,
    unnest((SELECT ddsketch_percentile(ddsketch(data.v, 0.01, 1024), perc.p) FROM data, perc GROUP BY perc.p)) AS a,
    unnest((SELECT lower_quantile(v, p) FROM data, perc)) AS b
) foo;

-- <value,count> API, array of percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT NULL AS v, 1 AS c UNION ALL SELECT * FROM (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text)) foo),
  data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c))
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
  SELECT
    unnest((SELECT p FROM perc)) p,
    unnest((SELECT ddsketch_percentile(ddsketch(data.v, data.c, 0.01, 1024), perc.p) FROM data, perc GROUP BY perc.p)) AS a,
    unnest((SELECT lower_quantile(v, p) FROM data_exp, perc)) AS b
) foo;

-- <value,count> API, but count is NULL (should be treated as 1), array of percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT NULL AS v, 1 AS c UNION ALL SELECT * FROM (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text)) foo)
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
  SELECT
    unnest((SELECT p FROM perc)) p,
    unnest((SELECT ddsketch_percentile(ddsketch(data.v, NULL, 0.01, 1024), perc.p) FROM data, perc GROUP BY perc.p)) AS a,
    unnest((SELECT lower_quantile(v, p) FROM data, perc)) AS b
) foo;


-- percentile-of values

-- individual values
WITH
    vals AS (SELECT (i*10.0) AS f FROM generate_series(-10,10) s(i)),
    data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text))
SELECT
    f,
    abs(a - b) < 0.05
FROM (
    SELECT
       f,
       (SELECT ddsketch_percentile_of(ddsketch(data.v, 0.01, 1024), f) FROM data) AS a,
       (SELECT percent_rank(f) WITHIN GROUP (ORDER BY v) FROM data WHERE data.v IS NOT NULL) AS b
    FROM vals
) foo;

WITH
    vals AS (SELECT (i*10.0) AS f FROM generate_series(-10,10) s(i)),
    data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text)),
    data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c))
SELECT
    f,
    abs(a - b) < 0.05
FROM (
    SELECT
       f,
       (SELECT ddsketch_percentile_of(ddsketch(data.v, data.c, 0.01, 1024), f) FROM data) AS a,
       (SELECT percent_rank(f) WITHIN GROUP (ORDER BY v) FROM data_exp WHERE data_exp.v IS NOT NULL) AS b
    FROM vals
) foo;

WITH
    vals AS (SELECT (i*10.0) AS f FROM generate_series(-10,10) s(i)),
    data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text))
SELECT
    f,
    abs(a - b) < 0.05
FROM (
    SELECT
       f,
       (SELECT ddsketch_percentile_of(ddsketch(data.v, NULL, 0.01, 1024), f) FROM data) AS a,
       (SELECT percent_rank(f) WITHIN GROUP (ORDER BY v) FROM data WHERE data.v IS NOT NULL) AS b
    FROM vals
) foo;

-- array of values

WITH data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text)),
     vals AS (SELECT array_agg(i) AS v FROM generate_series(-100,100,10) AS s(i)),
     sketch AS (SELECT ddsketch_percentile_of(ddsketch(data.v, data.c, 0.01, 1024), vals.v) s FROM data, vals GROUP BY vals.v),
     ranks AS (SELECT array_agg((SELECT percent_rank(vals) WITHIN GROUP (ORDER BY data.v) FROM data WHERE data.v IS NOT NULL)) AS r FROM unnest((SELECT v FROM vals)) vals)
SELECT
    v,
    abs(a - b) < 0.05
FROM (
  SELECT
    unnest(vals.v) AS v,
    unnest(sketch.s) AS a,
    unnest(ranks.r) AS b
  FROM vals, sketch, ranks
) foo;

WITH data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text)),
     data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c)),
     vals AS (SELECT array_agg(i) AS v FROM generate_series(-100,100,10) AS s(i)),
     sketch AS (SELECT ddsketch_percentile_of(ddsketch(data.v, data.c, 0.01, 1024), vals.v) s FROM data, vals GROUP BY vals.v),
     ranks AS (SELECT array_agg((SELECT percent_rank(vals) WITHIN GROUP (ORDER BY data_exp.v) FROM data_exp WHERE data_exp.v IS NOT NULL)) AS r FROM unnest((SELECT v FROM vals)) vals)
SELECT
    v,
    abs(a - b) < 0.05
FROM (
  SELECT
    unnest(vals.v) AS v,
    unnest(sketch.s) AS a,
    unnest(ranks.r) AS b
  FROM vals, sketch, ranks
) foo;

WITH data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text)),
     vals AS (SELECT array_agg(i) AS v FROM generate_series(-100,100,10) AS s(i)),
     sketch AS (SELECT ddsketch_percentile_of(ddsketch(data.v, NULL, 0.01, 1024), vals.v) s FROM data, vals GROUP BY vals.v),
     ranks AS (SELECT array_agg((SELECT percent_rank(vals) WITHIN GROUP (ORDER BY data.v) FROM data WHERE data.v IS NOT NULL)) AS r FROM unnest((SELECT v FROM vals)) vals)
SELECT
    v,
    abs(a - b) < 0.05
FROM (
  SELECT
    unnest(vals.v) AS v,
    unnest(sketch.s) AS a,
    unnest(ranks.r) AS b
  FROM vals, sketch, ranks
) foo;


-- NULL at the beginning of the data

WITH data AS (SELECT NULL AS v, 1 AS c UNION ALL SELECT * FROM (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,1000) s(i) ORDER BY md5(i::text)) foo),
     vals AS (SELECT array_agg(i) AS v FROM generate_series(-100,100,10) AS s(i)),
     sketch AS (SELECT ddsketch_percentile_of(ddsketch(data.v, data.c, 0.01, 1024), vals.v) s FROM data, vals GROUP BY vals.v),
     ranks AS (SELECT array_agg((SELECT percent_rank(vals) WITHIN GROUP (ORDER BY data.v) FROM data WHERE data.v IS NOT NULL)) AS r FROM unnest((SELECT v FROM vals)) vals)
SELECT
    v,
    abs(a - b) < 0.05
FROM (
  SELECT
    unnest(vals.v) AS v,
    unnest(sketch.s) AS a,
    unnest(ranks.r) AS b
  FROM vals, sketch, ranks
) foo;

WITH data AS (SELECT NULL AS v, 1 AS c UNION ALL SELECT * FROM (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text)) foo),
     data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c)),
     vals AS (SELECT array_agg(i) AS v FROM generate_series(-100,100,10) AS s(i)),
     sketch AS (SELECT ddsketch_percentile_of(ddsketch(data.v, data.c, 0.01, 1024), vals.v) s FROM data, vals GROUP BY vals.v),
     ranks AS (SELECT array_agg((SELECT percent_rank(vals) WITHIN GROUP (ORDER BY data_exp.v) FROM data_exp WHERE data_exp.v IS NOT NULL)) AS r FROM unnest((SELECT v FROM vals)) vals)
SELECT
    v,
    abs(a - b) < 0.05
FROM (
  SELECT
    unnest(vals.v) AS v,
    unnest(sketch.s) AS a,
    unnest(ranks.r) AS b
  FROM vals, sketch, ranks
) foo;

WITH data AS (SELECT NULL AS v, 1 AS c UNION ALL SELECT * FROM (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,2500) s(i) ORDER BY md5(i::text)) foo),
     vals AS (SELECT array_agg(i) AS v FROM generate_series(-100,100,10) AS s(i)),
     sketch AS (SELECT ddsketch_percentile_of(ddsketch(data.v, NULL, 0.01, 1024), vals.v) s FROM data, vals GROUP BY vals.v),
     ranks AS (SELECT array_agg((SELECT percent_rank(vals) WITHIN GROUP (ORDER BY data.v) FROM data WHERE data.v IS NOT NULL)) AS r FROM unnest((SELECT v FROM vals)) vals)
SELECT
    v,
    abs(a - b) < 0.05
FROM (
  SELECT
    unnest(vals.v) AS v,
    unnest(sketch.s) AS a,
    unnest(ranks.r) AS b
  FROM vals, sketch, ranks
) foo;


-- sketches
WITH data AS (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,2500) s(i)),
     sketches AS (SELECT ddsketch(data.v, 0.05, 1024) AS d FROM data GROUP BY c),
     perc AS (SELECT array_agg(i / 100.0) AS p FROM generate_series(0,100,10) s(i))
SELECT
  round(p, 2) as p,
  check_relative_error(a, b, 0.05) AS check_error,
  print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile(ddsketch(sketches.d), p) FROM sketches) AS a,
     (SELECT lower_quantile(v, p) FROM data) AS b
  FROM unnest((SELECT p FROM perc)) p) foo;

WITH data AS (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,2500) s(i)),
     data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c)),
     sketches AS (SELECT ddsketch(data.v, data.c, 0.05, 1024) AS d FROM data GROUP BY c),
     perc AS (SELECT array_agg(i / 100.0) AS p FROM generate_series(0,100,10) s(i))
SELECT
  round(p, 2) as p,
  check_relative_error(a, b, 0.05) AS check_error,
  print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile(ddsketch(sketches.d), p) FROM sketches) AS a,
     (SELECT lower_quantile(v, p) FROM data_exp) AS b
  FROM unnest((SELECT p FROM perc)) p) foo;

-- NULL is equal to count=1
WITH data AS (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,2500) s(i)),
     sketches AS (SELECT ddsketch(data.v, NULL, 0.05, 1024) AS d FROM data GROUP BY c),
     perc AS (SELECT array_agg(i / 100.0) AS p FROM generate_series(0,100,10) s(i))
SELECT
  round(p, 2) as p,
  check_relative_error(a, b, 0.05) AS check_error,
  print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile(ddsketch(sketches.d), p) FROM sketches) AS a,
     (SELECT lower_quantile(v, p) FROM data) AS b
  FROM unnest((SELECT p FROM perc)) p) foo;

-- percentile_of
WITH data AS (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,5000) s(i)),
     sketches AS (SELECT ddsketch(data.v, 0.05, 1024) AS d FROM data GROUP BY c),
     vals AS (SELECT array_agg(i) AS p FROM generate_series(-100,100,10) s(i))
SELECT
  p,
  abs(a - b) < 0.05
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile_of(ddsketch(sketches.d), p) FROM sketches) AS a,
     (SELECT percent_rank(p) WITHIN GROUP (ORDER BY v) FROM data WHERE data.v IS NOT NULL) AS b
  FROM unnest((SELECT p FROM vals)) p) foo;

WITH data AS (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,5000) s(i)),
     data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c)),
     sketches AS (SELECT ddsketch(data.v, data.c, 0.05, 1024) AS d FROM data GROUP BY c),
     vals AS (SELECT array_agg(i) AS p FROM generate_series(-100,100,10) s(i))
SELECT
  p,
  abs(a - b) < 0.05
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile_of(ddsketch(sketches.d), p) FROM sketches) AS a,
     (SELECT percent_rank(p) WITHIN GROUP (ORDER BY v) FROM data_exp WHERE data_exp.v IS NOT NULL) AS b
  FROM unnest((SELECT p FROM vals)) p) foo;

WITH data AS (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,5000) s(i)),
     sketches AS (SELECT ddsketch(data.v, 0.05, 1024) AS d FROM data GROUP BY c),
     vals AS (SELECT array_agg(i) AS p FROM generate_series(-100,100,10) s(i))
SELECT
  p,
  abs(a - b) < 0.05
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile_of(ddsketch(sketches.d), p) FROM sketches) AS a,
     (SELECT percent_rank(p) WITHIN GROUP (ORDER BY v) FROM data WHERE data.v IS NOT NULL) AS b
  FROM unnest((SELECT p FROM vals)) p) foo;

-- sketches (first value is NULL)
WITH data AS (SELECT 1 AS c, NULL AS v UNION ALL SELECT * FROM (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,2500) s(i)) foo),
     sketches AS (SELECT ddsketch(data.v, 0.05, 1024) AS d FROM data GROUP BY c),
     perc AS (SELECT array_agg(i / 100.0) AS p FROM generate_series(0,100,10) s(i))
SELECT
  round(p, 2) as p,
  check_relative_error(a, b, 0.05) AS check_error,
  print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile(ddsketch(sketches.d), p) FROM sketches) AS a,
     (SELECT lower_quantile(v, p) FROM data) AS b
  FROM unnest((SELECT p FROM perc)) p) foo;

WITH data AS (SELECT 1 AS c, NULL AS v UNION ALL SELECT * FROM (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,2500) s(i)) foo),
     data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c)),
     sketches AS (SELECT ddsketch(data.v, data.c, 0.05, 1024) AS d FROM data GROUP BY c),
     perc AS (SELECT array_agg(i / 100.0) AS p FROM generate_series(0,100,10) s(i))
SELECT
  round(p, 2) as p,
  check_relative_error(a, b, 0.05) AS check_error,
  print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile(ddsketch(sketches.d), p) FROM sketches) AS a,
     (SELECT lower_quantile(v, p) FROM data_exp) AS b
  FROM unnest((SELECT p FROM perc)) p) foo;

WITH data AS (SELECT 1 AS c, NULL AS v UNION ALL SELECT * FROM (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,2500) s(i)) foo),
     sketches AS (SELECT ddsketch(data.v, NULL, 0.05, 1024) AS d FROM data GROUP BY c),
     perc AS (SELECT array_agg(i / 100.0) AS p FROM generate_series(0,100,10) s(i))
SELECT
  round(p, 2) as p,
  check_relative_error(a, b, 0.05) AS check_error,
  print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile(ddsketch(sketches.d), p) FROM sketches) AS a,
     (SELECT lower_quantile(v, p) FROM data) AS b
  FROM unnest((SELECT p FROM perc)) p) foo;

-- NULL is equal to count=1
WITH data AS (SELECT 1 AS c, NULL AS v UNION ALL SELECT * FROM (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,2500) s(i)) foo),
     sketches AS (SELECT ddsketch(data.v, NULL, 0.05, 1024) AS d FROM data GROUP BY c),
     perc AS (SELECT array_agg(i / 100.0) AS p FROM generate_series(0,100,10) s(i))
SELECT
  round(p, 2) as p,
  check_relative_error(a, b, 0.05) AS check_error,
  print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile(ddsketch(sketches.d), p) FROM sketches) AS a,
     (SELECT lower_quantile(v, p) FROM data) AS b
  FROM unnest((SELECT p FROM perc)) p) foo;

-- percentile_of
WITH data AS (SELECT 1 AS c, NULL AS v UNION ALL SELECT * FROM (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,5000) s(i)) foo),
     sketches AS (SELECT ddsketch(data.v, 0.05, 1024) AS d FROM data GROUP BY c),
     vals AS (SELECT array_agg(i) AS p FROM generate_series(-100,100,10) s(i))
SELECT
  p,
  abs(a - b) < 0.05
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile_of(ddsketch(sketches.d), p) FROM sketches) AS a,
     (SELECT percent_rank(p) WITHIN GROUP (ORDER BY v) FROM data WHERE data.v IS NOT NULL) AS b
  FROM unnest((SELECT p FROM vals)) p) foo;

WITH data AS (SELECT 1 AS c, NULL AS v UNION ALL SELECT * FROM (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,5000) s(i)) foo),
     data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c)),
     sketches AS (SELECT ddsketch(data.v, data.c, 0.05, 1024) AS d FROM data GROUP BY c),
     vals AS (SELECT array_agg(i) AS p FROM generate_series(-100,100,10) s(i))
SELECT
  p,
  abs(a - b) < 0.05
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile_of(ddsketch(sketches.d), p) FROM sketches) AS a,
     (SELECT percent_rank(p) WITHIN GROUP (ORDER BY v) FROM data_exp WHERE data_exp.v IS NOT NULL) AS b
  FROM unnest((SELECT p FROM vals)) p) foo;

WITH data AS (SELECT 1 AS c, NULL AS v UNION ALL SELECT * FROM (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,5000) s(i)) foo),
     sketches AS (SELECT ddsketch(data.v, 0.05, 1024) AS d FROM data GROUP BY c),
     vals AS (SELECT array_agg(i) AS p FROM generate_series(-100,100,10) s(i))
SELECT
  p,
  abs(a - b) < 0.05
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile_of(ddsketch(sketches.d), p) FROM sketches) AS a,
     (SELECT percent_rank(p) WITHIN GROUP (ORDER BY v) FROM data WHERE data.v IS NOT NULL) AS b
  FROM unnest((SELECT p FROM vals)) p) foo;



SELECT ddsketch_percentile(NULL::ddsketch, 0.9);
SELECT ddsketch_percentile(NULL::ddsketch, ARRAY[0.5, 0.9]);
SELECT ddsketch_percentile_of(NULL::ddsketch, 0.9);
SELECT ddsketch_percentile_of(NULL::ddsketch, ARRAY[0.1, 0.9]);

SELECT ddsketch(NULL::ddsketch) FROM generate_series(1,10);
