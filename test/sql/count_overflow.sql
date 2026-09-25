-- INT64_MAX is valid, but one more observation must not wrap the count.
SELECT ddsketch_count(ddsketch(v, n, 0.05, 16))
FROM (VALUES (0.0, 9223372036854775806::bigint), (1.0, 1)) t(v, n);
SELECT ddsketch(v, n, 0.05, 16)
FROM (VALUES (0.0, 9223372036854775807::bigint), (1.0, 1)) t(v, n);

CREATE TEMP TABLE full_sketch AS SELECT ddsketch(1.0, 9223372036854775807::bigint, 0.05, 16) AS s;
SELECT ddsketch_add(s, 1.0, 0.05, 16) FROM full_sketch;
SELECT ddsketch_add(s, 0.0, 1::bigint, 0.05, 16) FROM full_sketch;
SELECT ddsketch_add(s, ARRAY[-1.0, 1.0]::float8[], 0.05, 16) FROM full_sketch;
SELECT ddsketch_union(s, ddsketch_add(NULL::ddsketch, 1.0, 0.05, 16)) FROM full_sketch;
SELECT ddsketch(s) FROM (
    SELECT s FROM full_sketch
    UNION ALL
    SELECT ddsketch_add(NULL::ddsketch, 1.0, 0.05, 16)
) t;
DROP TABLE full_sketch;

-- Partitionwise aggregation exercises the combine function without workers.
SET enable_partitionwise_aggregate = on;
SET max_parallel_workers_per_gather = 0;
CREATE TEMP TABLE overflow_parts (part integer, v double precision, n bigint) PARTITION BY LIST (part);
CREATE TEMP TABLE overflow_parts_a PARTITION OF overflow_parts FOR VALUES IN (0);
CREATE TEMP TABLE overflow_parts_b PARTITION OF overflow_parts FOR VALUES IN (1);
INSERT INTO overflow_parts VALUES (0, 0.0, 9223372036854775807), (1, 1.0, 1);
INSERT INTO overflow_parts SELECT i % 2, NULL, 1 FROM generate_series(1, 100) s(i);
ANALYZE overflow_parts;
EXPLAIN (COSTS OFF) SELECT ddsketch(v, n, 0.05, 16) FROM overflow_parts;
SELECT ddsketch(v, n, 0.05, 16) FROM overflow_parts;
DROP TABLE overflow_parts;
RESET enable_partitionwise_aggregate;
RESET max_parallel_workers_per_gather;
