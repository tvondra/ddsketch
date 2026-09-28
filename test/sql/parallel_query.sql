-----------------------------------------------------------
-- initialize configuration validation
-----------------------------------------------------------

DO $$
DECLARE
    v_version numeric;
BEGIN

    SELECT substring(setting from '\d+')::numeric INTO v_version FROM pg_settings WHERE name = 'server_version';

    -- GUCs common for all versions
    PERFORM set_config('parallel_setup_cost', '0', false);
    PERFORM set_config('parallel_tuple_cost', '0', false);
    PERFORM set_config('max_parallel_workers_per_gather', '2', false);

    -- 9.6 used somewhat different GUC name for relation size
    IF v_version < 10 THEN
        PERFORM set_config('min_parallel_relation_size', '1kB', false);
    ELSE
        PERFORM set_config('min_parallel_table_scan_size', '1kB', false);
    END IF;

    -- in 14 disable Memoize nodes, to make explain more consistent
    IF v_version >= 14 THEN
        PERFORM set_config('enable_memoize', 'off', false);
    END IF;

END;
$$ LANGUAGE plpgsql;

CREATE TABLE src_data (v double precision);

-- small amount of data, to allow weird cases when combining per-worker states
INSERT INTO src_data SELECT z FROM random_normal(1000, mean := 0.0, stddev := 0.1, minval := -1.0, maxval := 1.0) s(z);
ANALYZE src_data;

-- with parallelism
EXPLAIN (COSTS OFF) SELECT ddsketch(v, 0.05, 1024) FROM src_data;
SELECT ddsketch(v, 0.05, 1024) FROM src_data;

EXPLAIN (COSTS OFF) SELECT ddsketch(v, 0.05, 1024) FROM src_data;
SELECT trunc_value(ddsketch_percentile(ddsketch(v, 0.05, 1024), 0.9)) FROM src_data;

EXPLAIN (COSTS OFF) SELECT ddsketch(v, 0.05, 1024) FROM src_data;
SELECT ddsketch_percentile_of(ddsketch(v, 0.05, 1024), 0.9) FROM src_data;

-- without  parallelism
SET max_parallel_workers_per_gather = 0;
EXPLAIN (COSTS OFF) SELECT ddsketch(v, 0.05, 1024) FROM src_data;
SELECT ddsketch(v, 0.05, 1024) FROM src_data;

DROP TABLE src_data;


-- can't merge sketches with different alpha values
WITH sketches AS (
    SELECT ddsketch(i/100.0, 0.01, 1024) AS s FROM generate_series(1,10000) s(i)
    UNION ALL
    SELECT ddsketch(i/100.0, 0.05, 1024) AS s FROM generate_series(1,10000) s(i)
)
SELECT ddsketch(sketches.s) FROM sketches;

WITH sketches AS (
    SELECT ddsketch(i/100.0, 0.01, 1024) AS s FROM generate_series(1,10000) s(i)
    UNION ALL
    SELECT ddsketch(i/100.0, 0.05, 1024) AS s FROM generate_series(1,10000) s(i)
)
SELECT ddsketch_percentile(ddsketch(sketches.s), 0.9) FROM sketches;

WITH sketches AS (
    SELECT ddsketch(i/100.0, 0.01, 1024) AS s FROM generate_series(1,10000) s(i)
    UNION ALL
    SELECT ddsketch(i/100.0, 0.05, 1024) AS s FROM generate_series(1,10000) s(i)
)
SELECT ddsketch_percentile(ddsketch(sketches.s), ARRAY[0.95, 0.99]) FROM sketches;

WITH sketches AS (
    SELECT ddsketch(i/100.0, 0.01, 1024) AS s FROM generate_series(1,10000) s(i)
    UNION ALL
    SELECT ddsketch(i/100.0, 0.05, 1024) AS s FROM generate_series(1,10000) s(i)
)
SELECT ddsketch_percentile_of(ddsketch(sketches.s), 95) FROM sketches;

WITH sketches AS (
    SELECT ddsketch(i/100.0, 0.01, 1024) AS s FROM generate_series(1,10000) s(i)
    UNION ALL
    SELECT ddsketch(i/100.0, 0.05, 1024) AS s FROM generate_series(1,10000) s(i)
)
SELECT ddsketch_percentile_of(ddsketch(sketches.s), ARRAY[95, 99]) FROM sketches;


WITH sketches AS (
    SELECT ddsketch(i/100.0, 0.01, 1024) AS s FROM generate_series(1,10000) s(i)
)
SELECT ddsketch_union(sketches.s, NULL::ddsketch) FROM sketches;

WITH sketches AS (
    SELECT ddsketch(i/100.0, 0.05, 1024) AS s FROM generate_series(1,10000) s(i)
)
SELECT ddsketch_union(NULL::ddsketch, sketches.s) FROM sketches;

SELECT ddsketch_union(NULL::ddsketch, NULL::ddsketch);

-- Unequal partial zero counts detect doubling the destination instead of adding the source.
SET enable_partitionwise_aggregate = on;
SET max_parallel_workers_per_gather = 0;
CREATE TEMP TABLE zero_parts (part integer, v double precision, n bigint) PARTITION BY LIST (part);
CREATE TEMP TABLE zero_parts_a PARTITION OF zero_parts FOR VALUES IN (0);
CREATE TEMP TABLE zero_parts_b PARTITION OF zero_parts FOR VALUES IN (1);
INSERT INTO zero_parts SELECT 0, 0.0, 1 FROM generate_series(1, 10);
INSERT INTO zero_parts SELECT 1, 0.0, 3 FROM generate_series(1, 10);
ANALYZE zero_parts;
EXPLAIN (COSTS OFF) SELECT ddsketch(v, n, 0.05, 16) FROM zero_parts;
SELECT count, zero_count
FROM ddsketch_info((SELECT ddsketch(v, n, 0.05, 16) FROM zero_parts));
DROP TABLE zero_parts;
RESET enable_partitionwise_aggregate;
RESET max_parallel_workers_per_gather;

-- verify parallel safety of scalar functions
CREATE TABLE parallel_safety_data AS
SELECT (i % 10 - 5)::float8 AS v FROM generate_series(1, 10000) s(i);
ANALYZE parallel_safety_data;
SET max_parallel_workers_per_gather = 2;
SET min_parallel_table_scan_size = 0;
SET parallel_setup_cost = 0;
SET parallel_tuple_cost = 0;

-- transition functions allow partial aggregation despite missing labels
EXPLAIN (COSTS OFF)
SELECT ddsketch_count(ddsketch(v, 0.05, 16)),
       ddsketch_count(ddsketch(v, 2::bigint, 0.05, 16))
FROM parallel_safety_data;
SELECT ddsketch_count(ddsketch(v, 0.05, 16)) AS plain_count,
       ddsketch_count(ddsketch(v, 2::bigint, 0.05, 16)) AS weighted_count
FROM parallel_safety_data;

-- incrementally constructed inputs can also be built in parallel workers
EXPLAIN (COSTS OFF)
SELECT ddsketch(ddsketch_add(NULL::ddsketch, v, 0.05, 16)) FROM parallel_safety_data;
SELECT ddsketch_count(ddsketch(ddsketch_add(NULL::ddsketch, v, 0.05, 16))) AS count
FROM parallel_safety_data;

-- text casts and direct send calls are planner-visible, unlike support calls
EXPLAIN (COSTS OFF)
SELECT ddsketch(ddsketch_add(NULL::ddsketch, v, 0.05, 16)::text::ddsketch),
       count(ddsketch_send(ddsketch_add(NULL::ddsketch, v, 0.05, 16)))
FROM parallel_safety_data;
SELECT ddsketch_count(ddsketch(ddsketch_add(NULL::ddsketch, v, 0.05, 16)::text::ddsketch)) AS restored_count,
       count(ddsketch_send(ddsketch_add(NULL::ddsketch, v, 0.05, 16))) AS sent_count
FROM parallel_safety_data;

RESET max_parallel_workers_per_gather;
RESET min_parallel_table_scan_size;
RESET parallel_setup_cost;
RESET parallel_tuple_cost;
DROP TABLE parallel_safety_data;
