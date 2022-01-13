\set ECHO none

-- disable the notices for the create script (shell types etc.)
SET client_min_messages = 'WARNING';
\i ddsketch--1.0.0.sql
SET client_min_messages = 'NOTICE';

create extension lower_quantile;

\set ECHO all

-- SRF function implementing a simple deterministict PRNG

CREATE OR REPLACE FUNCTION prng(nrows int, seed int = 23982, p1 bigint = 16807, p2 bigint = 0, n bigint = 2147483647) RETURNS SETOF double precision AS $$
DECLARE
    val INT := seed;
BEGIN
    FOR i IN 1..nrows LOOP
        val := (val * p1 + p2) % n;

        RETURN NEXT (val::double precision / n);
    END LOOP;

    RETURN;
END;
$$ LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION random_normal(nrows int, mean double precision = 0.5, stddev double precision = 0.1, minval double precision = 0.0, maxval double precision = 1.0, seed int = 23982, p1 bigint = 16807, p2 bigint = 0, n bigint = 2147483647) RETURNS SETOF double precision AS $$
DECLARE
    v BIGINT := seed;
    x DOUBLE PRECISION;
    y DOUBLE PRECISION;
    s DOUBLE PRECISION;
    r INT := nrows;
BEGIN

    WHILE true LOOP

        -- random x
        v := (v * p1 + p2) % n;
        x := 2 * v / n::double precision - 1.0;

        -- random y
        v := (v * p1 + p2) % n;
        y := 2 * v / n::double precision - 1.0;

        s := x^2 + y^2;

        IF s != 0.0 AND s < 1.0 THEN

            s = sqrt(-2 * ln(s) / s);

            x := mean + stddev * s * x;

            IF x >= minval AND x <= maxval THEN
                RETURN NEXT x;
                r := r - 1;
            END IF;

            EXIT WHEN r = 0;

            y := mean + stddev * s * y;

            IF y >= minval AND y <= maxval THEN
                RETURN NEXT y;
                r := r - 1;
            END IF;

            EXIT WHEN r = 0;

        END IF;

    END LOOP;

END;
$$ LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION relative_error(estimated_value double precision, actual_value double precision) RETURNS double precision AS $$
DECLARE
    divisor double precision;
BEGIN
    -- if both values are 0, then the relative error is 0
    if (abs(actual_value) = 0) and (abs(estimated_value) = 0) then
        return 0.0;
    end if;

    -- use actual value az divisor, but if it's 0 then use the estimate
    -- this allows us to calculate the error without division by zero
    -- (we already know both can't be zero)
    divisor := abs(actual_value);
    if divisor = 0 then
      divisor := abs(estimated_value);
    end if;

    return abs(estimated_value - actual_value) / divisor;
END;
$$ LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION check_relative_error(estimated_value double precision, actual_value double precision, allowed_error double precision) RETURNS bool AS $$
DECLARE
    err double precision;
BEGIN

    IF ((estimated_value < 0) AND (actual_value > 0)) OR ((estimated_value > 0) AND (actual_value < 0)) THEN
        RETURN NULL;
    END IF;

    -- add 1% fuzz factor, to account for rounding errors etc.
    err := relative_error(estimated_value, actual_value);

    RETURN (err < allowed_error * 1.01);

END;
$$ LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION print_relative_error(estimated_value double precision, actual_value double precision, allowed_error double precision) RETURNS text AS $$
DECLARE
    err double precision;
BEGIN

    IF ((estimated_value < 0) AND (actual_value > 0)) OR ((estimated_value > 0) AND (actual_value < 0)) THEN
        RETURN format('estimate = %s, actual = %s', estimated_value, actual_value);
    END IF;

    err := relative_error(estimated_value, actual_value);

    -- add 1% fuzz factor, to account for rounding errors etc.
    IF err < (allowed_error * 1.01) THEN
        RETURN NULL;
    END IF;

    RETURN format('estimate = %s, actual = %s, error = %s', estimated_value, actual_value, err);

END;
$$ LANGUAGE plpgsql;

-- functions to round double precision values (and arrays of)
CREATE OR REPLACE FUNCTION trunc_value(v double precision, s integer = 12) RETURNS text AS $$
SELECT substring(v::text, 1, s);
$$ LANGUAGE sql;

CREATE OR REPLACE FUNCTION trunc_value(v double precision[], s integer = 12) RETURNS text[] AS $$
SELECT array_agg(trunc_value(v, s)) FROM unnest(v) v;
$$ LANGUAGE sql;

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

-----------------------------------------------------------
-- parameter validation
-----------------------------------------------------------

-- invalid percentile value
SELECT ddsketch_percentile(i / 1.0, 0.01, 1024, ARRAY[0.1, -0.1]) FROM generate_series(1,10000) s(i);

-- alpha too low
SELECT ddsketch_percentile(i / 1.0, 0.00009, 1024, 0.5) FROM generate_series(1,10000) s(i);

-- alpha too high
SELECT ddsketch_percentile(i / 1.0, 0.11, 1024, 0.5) FROM generate_series(1,10000) s(i);

-- fewer than minimum number of buckets
SELECT ddsketch_percentile(i / 1.0, 0.01, 15, 0.5) FROM generate_series(1,10000) s(i);

-- more than maximum number of buckets
SELECT ddsketch_percentile(i / 1.0, 0.01, 32769, 0.5) FROM generate_series(1,10000) s(i);

-- too many buckets needed
SELECT ddsketch_percentile(i / 1.0, 0.01, 32, 0.5) FROM generate_series(1,10000) s(i);

-----------------------------------------------------------
-- nice data set with ordered (asc) / evenly-spaced data --
-----------------------------------------------------------

-- alpha = 0.05
WITH data AS (SELECT i AS x FROM generate_series(1,10000) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]::double precision[])) AS b
    FROM data
) foo;

WITH data AS (SELECT i - 5000 AS x FROM generate_series(1,10000) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]::double precision[])) AS b
    FROM data
) foo;

WITH data AS (SELECT i - 10000 AS x FROM generate_series(1,10000) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]::double precision[])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM generate_series(1,10000) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT i - 5000  AS x FROM generate_series(1,10000) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT i - 10000 AS x FROM generate_series(1,10000) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- alpha 0.01
WITH data AS (SELECT i AS x FROM generate_series(1,10000) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT i - 5000 AS x FROM generate_series(1,10000) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT i - 10000 AS x FROM generate_series(1,10000) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM generate_series(1,10000) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT i - 5000 AS x FROM generate_series(1,10000) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT i - 10000 AS x FROM generate_series(1,10000) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- 0.001 alpha
WITH data AS (SELECT i AS x FROM generate_series(1,10000) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 2048, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT i - 5000 AS x FROM generate_series(1,10000) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 4096, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT i - 10000 AS x FROM generate_series(1,10000) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 2048, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM generate_series(1,10000) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 2048, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT i - 5000 AS x FROM generate_series(1,10000) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 4096, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT i - 10000 AS x FROM generate_series(1,10000) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 2048, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

------------------------------------------------------------
-- nice data set with ordered (desc) / evenly-spaced data --
------------------------------------------------------------

-- 0.05 alpha
WITH data AS (SELECT i AS x FROM generate_series(10000,1,-1) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;


WITH data AS (SELECT i - 5000 AS x FROM generate_series(10000,1,-1) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;


WITH data AS (SELECT i - 10000 AS x FROM generate_series(10000,1,-1) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM generate_series(10000,1,-1) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT i - 5000 AS x FROM generate_series(10000,1,-1) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT i - 10000 AS x FROM generate_series(10000,1,-1) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- alpha 0.01
WITH data AS (SELECT i AS x FROM generate_series(10000,1,-1) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT i - 5000 AS x FROM generate_series(10000,1,-1) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT i - 10000 AS x FROM generate_series(10000,1,-1) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM generate_series(10000,1,-1) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT i - 5000 AS x FROM generate_series(10000,1,-1) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT i - 10000 AS x FROM generate_series(10000,1,-1) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- alpha 0.001
WITH data AS (SELECT i AS x FROM generate_series(10000,1,-1) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 2048, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT i - 5000 AS x FROM generate_series(10000,1,-1) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 4096, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT i - 10000 AS x FROM generate_series(10000,1,-1) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 2048, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM generate_series(10000,1,-1) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 2048, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT i - 5000 AS x FROM generate_series(10000,1,-1) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 4096, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT i - 10000 AS x FROM generate_series(10000,1,-1) s(i)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 2048, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

----------------------------------------------------
-- nice data set with random / evenly-spaced data --
----------------------------------------------------

-- alpha 0.05
WITH data AS (SELECT i AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo)
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT i - 5000 AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo)
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT i - 10000 AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo)
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;


WITH data AS (SELECT i - 5000 AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT i - 10000 AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- alpha 0.01
WITH data AS (SELECT i AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo)
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT i - 5000 AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo)
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT i - 10000 AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo)
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT i - 5000 AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT i - 10000 AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- alpha 0.001
WITH data AS (SELECT i AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo)
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT i - 5000 AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo)
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT i - 10000 AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo)
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 8192, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT i - 5000 AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 8192, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT i - 10000 AS x FROM (SELECT generate_series(1,10000) AS i, prng(10000, 49979693) AS x ORDER BY x) foo),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 8192, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

----------------------------------------------
-- nice data set with random data (uniform) --
----------------------------------------------

-- alpha 0.05
WITH data AS (SELECT 10000 * x AS x FROM prng(10000) s(x))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * x - 5000 AS x FROM prng(10000) s(x))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * x - 10000 AS x FROM prng(10000) s(x))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * x AS x FROM prng(10000) s(x)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * x - 5000 AS x FROM prng(10000) s(x)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * x - 10000 AS x FROM prng(10000) s(x)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- alpha 0.01
WITH data AS (SELECT 10000 * x AS x FROM prng(10000) s(x))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * x - 5000 AS x FROM prng(10000) s(x))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * x - 10000 AS x FROM prng(10000) s(x))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * x AS x FROM prng(10000) s(x)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * x - 5000 AS x FROM prng(10000) s(x)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * x - 10000 AS x FROM prng(10000) s(x)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- alpha 0.001
WITH data AS (SELECT 10000 * x AS x FROM prng(10000) s(x))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * x - 5000 AS x FROM prng(10000) s(x))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * x - 10000 AS x FROM prng(10000) s(x))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * x AS x FROM prng(10000) s(x)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 8192, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * x - 5000 AS x FROM prng(10000) s(x)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 8192, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * x - 10000 AS x FROM prng(10000) s(x)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 8192, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

--------------------------------------------------
-- nice data set with random data (skewed sqrt) --
--------------------------------------------------

-- alpha 0.05
WITH data AS (SELECT 10000 * sqrt(z) AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * sqrt(z) - 5000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * sqrt(z) - 10000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * sqrt(z) AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * sqrt(z) - 5000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * sqrt(z) - 10000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- alpha 0.01
WITH data AS (SELECT 10000 * sqrt(z) AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * sqrt(z) - 5000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * sqrt(z) - 10000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * sqrt(z) AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * sqrt(z) - 5000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * sqrt(z) - 10000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- 0.001 alpha
WITH data AS (SELECT 10000 * sqrt(z) AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * sqrt(z) - 5000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * sqrt(z) - 5000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * sqrt(z) AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 8192, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * sqrt(z) - 5000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 8192, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * sqrt(z) - 10000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 8192, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-------------------------------------------------------
-- nice data set with random data (skewed sqrt+sqrt) --
-------------------------------------------------------

-- alpha 0.05
WITH data AS (SELECT 10000 * sqrt(sqrt(z)) AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * sqrt(sqrt(z)) - 5000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * sqrt(sqrt(z)) - 10000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * sqrt(sqrt(z)) AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * sqrt(sqrt(z)) - 5000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * sqrt(sqrt(z)) - 10000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- alpha 0.01
WITH data AS (SELECT 10000 * sqrt(sqrt(z)) AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * sqrt(sqrt(z)) - 5000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * sqrt(sqrt(z)) - 10000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * sqrt(sqrt(z)) AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * sqrt(sqrt(z)) - 5000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * sqrt(sqrt(z)) - 10000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- alpha 0.001
WITH data AS (SELECT 10000 * sqrt(sqrt(z)) AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * sqrt(sqrt(z)) - 5000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * sqrt(sqrt(z)) - 10000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * sqrt(sqrt(z)) AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 8192, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * sqrt(sqrt(z)) - 5000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 8192, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * sqrt(sqrt(z)) - 10000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 8192, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-------------------------------------------------
-- nice data set with random data (skewed pow) --
-------------------------------------------------

-- alpha 0.05
WITH data AS (SELECT 10000 * pow(z, 2) AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * pow(z, 2) - 5000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * pow(z, 2) AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * pow(z, 2) - 5000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- alpha 0.01
WITH data AS (SELECT 10000 * pow(z, 2) AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * pow(z, 2) - 5000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * pow(z, 2) AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * pow(z, 2) - 5000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- 0.001 alpha
WITH data AS (SELECT 10000 * pow(z, 2) AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * pow(z, 2) - 5000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * pow(z, 2) AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 8192, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * pow(z, 2) - 5000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 8192, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-----------------------------------------------------
-- nice data set with random data (skewed pow+pow) --
-----------------------------------------------------

-- 0.05 alpha
WITH data AS (SELECT 10000 * pow(z, 4) AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * pow(z, 4) - 5000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * pow(z, 4) AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * pow(z, 4) - 5000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- alpha 0.01
WITH data AS (SELECT 10000 * pow(z, 4) AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 2048, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * pow(z, 4) - 5000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * pow(z, 4) AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 2048, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * pow(z, 4) - 5000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- 0.001 alpha
WITH data AS (SELECT 10000 * pow(z, 4) AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 16384, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

WITH data AS (SELECT 10000 * pow(z, 4) - 5000 AS x FROM prng(10000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * pow(z, 4) AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 16384, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

WITH data AS (SELECT 10000 * pow(z, 4) - 5000 AS x FROM prng(10000) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 8192, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

----------------------------------------------------------
-- nice data set with random data (normal distribution) --
----------------------------------------------------------

-- 0.05 alpha
WITH data AS (SELECT 10000 * pow(z, 3) AS x FROM random_normal(10000, mean := 0.0, stddev := 0.1, minval := -1.0, maxval := 1.0) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * pow(z, 3) AS x FROM random_normal(10000, mean := 0.0, stddev := 0.1, minval := -1.0, maxval := 1.0) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.05, 1024, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- 0.01 alpha
WITH data AS (SELECT 10000 * pow(z, 3) AS x FROM random_normal(10000, mean := 0.0, stddev := 0.1, minval := -1.0, maxval := 1.0) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 4096, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * pow(z, 3) AS x FROM random_normal(10000, mean := 0.0, stddev := 0.1, minval := -1.0, maxval := 1.0) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.01, 4096, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- 0.001 alpha
WITH data AS (SELECT 10000 * pow(z, 3) AS x FROM random_normal(10000, mean := 0.0, stddev := 0.1, minval := -1.0, maxval := 1.0) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 32768, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(lower_quantile(x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 10000 * pow(z, 3) AS x FROM random_normal(10000, mean := 0.0, stddev := 0.1, minval := -1.0, maxval := 1.0) s(z)),
     perc AS (SELECT array_agg((i/100.0)::double precision) AS p FROM generate_series(1,99) s(i))
SELECT * FROM (
    SELECT
        p,
        a,
        LAG(a) OVER (ORDER BY p) AS b
    FROM (
        SELECT
            unnest((SELECT p FROM perc)) AS p,
            unnest(ddsketch_percentile(x, 0.001, 32768, (SELECT p FROM perc))) AS a
        FROM data
    ) foo ) bar WHERE a < b;

-- some basic tests to verify transforming from and to text work
-- 0.05 alpha
WITH data AS (SELECT i AS x FROM generate_series(1,10000) s(i)),
     intermediate AS (SELECT ddsketch(x, 0.05, 1024)::text AS intermediate_x FROM data),
     ddsketch_parsed AS (SELECT ddsketch_percentile(intermediate_x::ddsketch, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS a FROM intermediate),
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
    FROM ddsketch_parsed,
         pg_percentile
) foo;

WITH data AS (SELECT i - 5000 AS x FROM generate_series(1,10000) s(i)),
     intermediate AS (SELECT ddsketch(x, 0.05, 1024)::text AS intermediate_x FROM data),
     ddsketch_parsed AS (SELECT ddsketch_percentile(intermediate_x::ddsketch, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS a FROM intermediate),
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
    FROM ddsketch_parsed,
         pg_percentile
) foo;

WITH data AS (SELECT i - 10000 AS x FROM generate_series(1,10000) s(i)),
     intermediate AS (SELECT ddsketch(x, 0.05, 1024)::text AS intermediate_x FROM data),
     ddsketch_parsed AS (SELECT ddsketch_percentile(intermediate_x::ddsketch, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS a FROM intermediate),
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
    FROM ddsketch_parsed,
         pg_percentile
) foo;

-- test various cases of invalid ddsketch text representations

-- invalid flags
select 'flags 1 count 1000 alpha 0.050000 zero_count 0 maxbuckets 128 buckets 51 0 (-75, 1) (-72, 1) (-68, 1) (-58, 1) (-51, 1) (-48, 1) (-46, 1) (-44, 1) (-43, 1) (-42, 1) (-40, 2) (-39, 1) (-38, 3) (-37, 4) (-36, 1) (-35, 2) (-34, 4) (-33, 2) (-32, 6) (-31, 9) (-30, 5) (-29, 4) (-28, 6) (-27, 9) (-26, 12) (-25, 8) (-24, 4) (-23, 8) (-22, 14) (-21, 8) (-20, 10) (-19, 13) (-18, 23) (-17, 13) (-16, 25) (-15, 20) (-14, 22) (-13, 20) (-12, 28) (-11, 28) (-10, 40) (-9, 35) (-8, 45) (-7, 48) (-6, 54) (-5, 61) (-4, 68) (-3, 79) (-2, 75) (-1, 77) (0, 94)'::ddsketch;

-- invalid count and zero_count
select 'flags 0 count 0 alpha 0.050000 zero_count 0 maxbuckets 128 buckets 51 0 (-75, 1) (-72, 1) (-68, 1) (-58, 1) (-51, 1) (-48, 1) (-46, 1) (-44, 1) (-43, 1) (-42, 1) (-40, 2) (-39, 1) (-38, 3) (-37, 4) (-36, 1) (-35, 2) (-34, 4) (-33, 2) (-32, 6) (-31, 9) (-30, 5) (-29, 4) (-28, 6) (-27, 9) (-26, 12) (-25, 8) (-24, 4) (-23, 8) (-22, 14) (-21, 8) (-20, 10) (-19, 13) (-18, 23) (-17, 13) (-16, 25) (-15, 20) (-14, 22) (-13, 20) (-12, 28) (-11, 28) (-10, 40) (-9, 35) (-8, 45) (-7, 48) (-6, 54) (-5, 61) (-4, 68) (-3, 79) (-2, 75) (-1, 77) (0, 94)'::ddsketch;
select 'flags 0 count 1000 alpha 0.050000 zero_count -10 maxbuckets 128 buckets 51 0 (-75, 1) (-72, 1) (-68, 1) (-58, 1) (-51, 1) (-48, 1) (-46, 1) (-44, 1) (-43, 1) (-42, 1) (-40, 2) (-39, 1) (-38, 3) (-37, 4) (-36, 1) (-35, 2) (-34, 4) (-33, 2) (-32, 6) (-31, 9) (-30, 5) (-29, 4) (-28, 6) (-27, 9) (-26, 12) (-25, 8) (-24, 4) (-23, 8) (-22, 14) (-21, 8) (-20, 10) (-19, 13) (-18, 23) (-17, 13) (-16, 25) (-15, 20) (-14, 22) (-13, 20) (-12, 28) (-11, 28) (-10, 40) (-9, 35) (-8, 45) (-7, 48) (-6, 54) (-5, 61) (-4, 68) (-3, 79) (-2, 75) (-1, 77) (0, 94)'::ddsketch;
select 'flags 0 count 1000 alpha 0.050000 zero_count 1001 maxbuckets 128 buckets 51 0 (-75, 1) (-72, 1) (-68, 1) (-58, 1) (-51, 1) (-48, 1) (-46, 1) (-44, 1) (-43, 1) (-42, 1) (-40, 2) (-39, 1) (-38, 3) (-37, 4) (-36, 1) (-35, 2) (-34, 4) (-33, 2) (-32, 6) (-31, 9) (-30, 5) (-29, 4) (-28, 6) (-27, 9) (-26, 12) (-25, 8) (-24, 4) (-23, 8) (-22, 14) (-21, 8) (-20, 10) (-19, 13) (-18, 23) (-17, 13) (-16, 25) (-15, 20) (-14, 22) (-13, 20) (-12, 28) (-11, 28) (-10, 40) (-9, 35) (-8, 45) (-7, 48) (-6, 54) (-5, 61) (-4, 68) (-3, 79) (-2, 75) (-1, 77) (0, 94)'::ddsketch;

-- count mismatching buckets
select 'flags 0 count 1001 alpha 0.050000 zero_count 0 maxbuckets 128 buckets 51 0 (-75, 1) (-72, 1) (-68, 1) (-58, 1) (-51, 1) (-48, 1) (-46, 1) (-44, 1) (-43, 1) (-42, 1) (-40, 2) (-39, 1) (-38, 3) (-37, 4) (-36, 1) (-35, 2) (-34, 4) (-33, 2) (-32, 6) (-31, 9) (-30, 5) (-29, 4) (-28, 6) (-27, 9) (-26, 12) (-25, 8) (-24, 4) (-23, 8) (-22, 14) (-21, 8) (-20, 10) (-19, 13) (-18, 23) (-17, 13) (-16, 25) (-15, 20) (-14, 22) (-13, 20) (-12, 28) (-11, 28) (-10, 40) (-9, 35) (-8, 45) (-7, 48) (-6, 54) (-5, 61) (-4, 68) (-3, 79) (-2, 75) (-1, 77) (0, 94)'::ddsketch;

-- invalid alpha values
select 'flags 0 count 1000 alpha 0.00005 zero_count 0 maxbuckets 128 buckets 51 0 (-75, 1) (-72, 1) (-68, 1) (-58, 1) (-51, 1) (-48, 1) (-46, 1) (-44, 1) (-43, 1) (-42, 1) (-40, 2) (-39, 1) (-38, 3) (-37, 4) (-36, 1) (-35, 2) (-34, 4) (-33, 2) (-32, 6) (-31, 9) (-30, 5) (-29, 4) (-28, 6) (-27, 9) (-26, 12) (-25, 8) (-24, 4) (-23, 8) (-22, 14) (-21, 8) (-20, 10) (-19, 13) (-18, 23) (-17, 13) (-16, 25) (-15, 20) (-14, 22) (-13, 20) (-12, 28) (-11, 28) (-10, 40) (-9, 35) (-8, 45) (-7, 48) (-6, 54) (-5, 61) (-4, 68) (-3, 79) (-2, 75) (-1, 77) (0, 94)'::ddsketch;
select 'flags 0 count 1000 alpha 0.11 zero_count 0 maxbuckets 128 buckets 51 0 (-75, 1) (-72, 1) (-68, 1) (-58, 1) (-51, 1) (-48, 1) (-46, 1) (-44, 1) (-43, 1) (-42, 1) (-40, 2) (-39, 1) (-38, 3) (-37, 4) (-36, 1) (-35, 2) (-34, 4) (-33, 2) (-32, 6) (-31, 9) (-30, 5) (-29, 4) (-28, 6) (-27, 9) (-26, 12) (-25, 8) (-24, 4) (-23, 8) (-22, 14) (-21, 8) (-20, 10) (-19, 13) (-18, 23) (-17, 13) (-16, 25) (-15, 20) (-14, 22) (-13, 20) (-12, 28) (-11, 28) (-10, 40) (-9, 35) (-8, 45) (-7, 48) (-6, 54) (-5, 61) (-4, 68) (-3, 79) (-2, 75) (-1, 77) (0, 94)'::ddsketch;

-- invalid maxbuckets
select 'flags 0 count 1000 alpha 0.05 zero_count 0 maxbuckets 1 buckets 51 0 (-75, 1) (-72, 1) (-68, 1) (-58, 1) (-51, 1) (-48, 1) (-46, 1) (-44, 1) (-43, 1) (-42, 1) (-40, 2) (-39, 1) (-38, 3) (-37, 4) (-36, 1) (-35, 2) (-34, 4) (-33, 2) (-32, 6) (-31, 9) (-30, 5) (-29, 4) (-28, 6) (-27, 9) (-26, 12) (-25, 8) (-24, 4) (-23, 8) (-22, 14) (-21, 8) (-20, 10) (-19, 13) (-18, 23) (-17, 13) (-16, 25) (-15, 20) (-14, 22) (-13, 20) (-12, 28) (-11, 28) (-10, 40) (-9, 35) (-8, 45) (-7, 48) (-6, 54) (-5, 61) (-4, 68) (-3, 79) (-2, 75) (-1, 77) (0, 94)'::ddsketch;
select 'flags 0 count 1000 alpha 0.05 zero_count 0 maxbuckets 65536 buckets 51 0 (-75, 1) (-72, 1) (-68, 1) (-58, 1) (-51, 1) (-48, 1) (-46, 1) (-44, 1) (-43, 1) (-42, 1) (-40, 2) (-39, 1) (-38, 3) (-37, 4) (-36, 1) (-35, 2) (-34, 4) (-33, 2) (-32, 6) (-31, 9) (-30, 5) (-29, 4) (-28, 6) (-27, 9) (-26, 12) (-25, 8) (-24, 4) (-23, 8) (-22, 14) (-21, 8) (-20, 10) (-19, 13) (-18, 23) (-17, 13) (-16, 25) (-15, 20) (-14, 22) (-13, 20) (-12, 28) (-11, 28) (-10, 40) (-9, 35) (-8, 45) (-7, 48) (-6, 54) (-5, 61) (-4, 68) (-3, 79) (-2, 75) (-1, 77) (0, 94)'::ddsketch;

-- invalid number of bucket
select 'flags 0 count 1000 alpha 0.050000 zero_count 0 maxbuckets 128 buckets -1 0 (-75, 1) (-72, 1) (-68, 1) (-58, 1) (-51, 1) (-48, 1) (-46, 1) (-44, 1) (-43, 1) (-42, 1) (-40, 2) (-39, 1) (-38, 3) (-37, 4) (-36, 1) (-35, 2) (-34, 4) (-33, 2) (-32, 6) (-31, 9) (-30, 5) (-29, 4) (-28, 6) (-27, 9) (-26, 12) (-25, 8) (-24, 4) (-23, 8) (-22, 14) (-21, 8) (-20, 10) (-19, 13) (-18, 23) (-17, 13) (-16, 25) (-15, 20) (-14, 22) (-13, 20) (-12, 28) (-11, 28) (-10, 40) (-9, 35) (-8, 45) (-7, 48) (-6, 54) (-5, 61) (-4, 68) (-3, 79) (-2, 75) (-1, 77) (0, 94)'::ddsketch;
select 'flags 0 count 1000 alpha 0.050000 zero_count 0 maxbuckets 128 buckets 51 -1 (-75, 1) (-72, 1) (-68, 1) (-58, 1) (-51, 1) (-48, 1) (-46, 1) (-44, 1) (-43, 1) (-42, 1) (-40, 2) (-39, 1) (-38, 3) (-37, 4) (-36, 1) (-35, 2) (-34, 4) (-33, 2) (-32, 6) (-31, 9) (-30, 5) (-29, 4) (-28, 6) (-27, 9) (-26, 12) (-25, 8) (-24, 4) (-23, 8) (-22, 14) (-21, 8) (-20, 10) (-19, 13) (-18, 23) (-17, 13) (-16, 25) (-15, 20) (-14, 22) (-13, 20) (-12, 28) (-11, 28) (-10, 40) (-9, 35) (-8, 45) (-7, 48) (-6, 54) (-5, 61) (-4, 68) (-3, 79) (-2, 75) (-1, 77) (0, 94)'::ddsketch;
select 'flags 0 count 1000 alpha 0.050000 zero_count 0 maxbuckets 128 buckets 129 0 (-75, 1) (-72, 1) (-68, 1) (-58, 1) (-51, 1) (-48, 1) (-46, 1) (-44, 1) (-43, 1) (-42, 1) (-40, 2) (-39, 1) (-38, 3) (-37, 4) (-36, 1) (-35, 2) (-34, 4) (-33, 2) (-32, 6) (-31, 9) (-30, 5) (-29, 4) (-28, 6) (-27, 9) (-26, 12) (-25, 8) (-24, 4) (-23, 8) (-22, 14) (-21, 8) (-20, 10) (-19, 13) (-18, 23) (-17, 13) (-16, 25) (-15, 20) (-14, 22) (-13, 20) (-12, 28) (-11, 28) (-10, 40) (-9, 35) (-8, 45) (-7, 48) (-6, 54) (-5, 61) (-4, 68) (-3, 79) (-2, 75) (-1, 77) (0, 94)'::ddsketch;
select 'flags 0 count 1000 alpha 0.050000 zero_count 0 maxbuckets 128 buckets 51 52 (-75, 1) (-72, 1) (-68, 1) (-58, 1) (-51, 1) (-48, 1) (-46, 1) (-44, 1) (-43, 1) (-42, 1) (-40, 2) (-39, 1) (-38, 3) (-37, 4) (-36, 1) (-35, 2) (-34, 4) (-33, 2) (-32, 6) (-31, 9) (-30, 5) (-29, 4) (-28, 6) (-27, 9) (-26, 12) (-25, 8) (-24, 4) (-23, 8) (-22, 14) (-21, 8) (-20, 10) (-19, 13) (-18, 23) (-17, 13) (-16, 25) (-15, 20) (-14, 22) (-13, 20) (-12, 28) (-11, 28) (-10, 40) (-9, 35) (-8, 45) (-7, 48) (-6, 54) (-5, 61) (-4, 68) (-3, 79) (-2, 75) (-1, 77) (0, 94)'::ddsketch;

-- invalid bucket count
select 'flags 0 count 1000 alpha 0.050000 zero_count 0 maxbuckets 128 buckets 51 0 (-75, 1) (-72, 0) (-68, 1) (-58, 1) (-51, 1) (-48, 1) (-46, 1) (-44, 1) (-43, 1) (-42, 1) (-40, 2) (-39, 1) (-38, 3) (-37, 4) (-36, 1) (-35, 2) (-34, 4) (-33, 2) (-32, 6) (-31, 9) (-30, 5) (-29, 4) (-28, 6) (-27, 9) (-26, 12) (-25, 8) (-24, 4) (-23, 8) (-22, 14) (-21, 8) (-20, 10) (-19, 13) (-18, 23) (-17, 13) (-16, 25) (-15, 20) (-14, 22) (-13, 20) (-12, 28) (-11, 28) (-10, 40) (-9, 35) (-8, 45) (-7, 48) (-6, 54) (-5, 61) (-4, 68) (-3, 79) (-2, 75) (-1, 77) (0, 94)'::ddsketch;
select 'flags 0 count 1000 alpha 0.050000 zero_count 0 maxbuckets 128 buckets 51 0 (-75, 1) (-72, 1) (-68, 1) (-58, 1) (-51, 1) (-48, 1) (-46, 1) (-44, 1) (-43, 1) (-42, 1) (-40, 2) (-39, 1) (-38, 3) (-37, 4) (-36, 1) (-35, 2) (-34, 4) (-33, 2) (-32, 6) (-31, 9) (-30, 5) (-29, 4) (-28, 6) (-27, 9) (-26, 12) (-25, 8) (-24, 4) (-23, 8) (-22, 14) (-21, 8) (-20, 10) (-19, 13) (-18, 23) (-17, 13) (-16, 25) (-15, 20) (-14, 22) (-13, 20) (-12, 28) (-11, 28) (-10, 40) (-9, 35) (-8, 45) (-7, 48) (-6, 54) (-5, 61) (-4, 68) (-3, 79) (-2, 75) (-1, 77) (0, 94) (1, 10)'::ddsketch;
select 'flags 0 count 1000 alpha 0.050000 zero_count 0 maxbuckets 128 buckets 51 0 (-75, 1) (-72, 1) (-68, 1) (-58, 1) (-51, 1) (-48, 1) (-46, 1) (-44, 1) (-43, 1) (-42, 1) (-40, 2) (-39, 1) (-38, 3) (-37, 4) (-36, 1) (-35, 2) (-34, 4) (-33, 2) (-32, 6) (-31, 9) (-30, 5) (-29, 4) (-28, 6) (-27, 9) (-26, 12) (-25, 8) (-24, 4) (-23, 8) (-22, 14) (-21, 8) (-20, 10) (-19, 13) (-18, 23) (-17, 13) (-16, 25) (-15, 20) (-14, 22) (-13, 20) (-12, 28) (-11, 28) (-10, 40) (-9, 35) (-8, 45) (-7, 48) (-6, 54) (-5, 61) (-4, 68) (-3, 79) (-2, 75) (-1, 171)'::ddsketch;

-- invalid indexes in negative/positive count
select 'flags 0 count 1000 alpha 0.050000 zero_count 0 maxbuckets 128 buckets 50 50 (0, 95) (-1, 104) (-2, 93) (-3, 78) (-4, 60) (-5, 54) (-6, 48) (-7, 46) (-8, 40) (-9, 34) (-10, 28) (-11, 40) (-12, 26) (-13, 21) (-14, 27) (-15, 23) (-16, 14) (-17, 12) (-18, 16) (-19, 12) (-20, 10) (-21, 9) (-22, 10) (-23, 9) (-24, 6) (-25, 10) (-26, 5) (-27, 5) (-28, 5) (-29, 6) (-30, 10) (-31, 7) (-32, 3) (-33, 2) (-34, 3) (-35, 3) (-36, 6) (-37, 3) (-38, 2) (-39, 4) (-41, 1) (-44, 1) (-46, 1) (-47, 1) (-49, 1) (-70, 1) (-52, 1) (-57, 1) (-62, 2) (-75, 1)'::ddsketch;
select 'flags 0 count 1000 alpha 0.050000 zero_count 0 maxbuckets 128 buckets 51 0 (-75, 1) (-72, 1) (-50, 1) (-58, 1) (-51, 1) (-48, 1) (-46, 1) (-44, 1) (-43, 1) (-42, 1) (-40, 2) (-39, 1) (-38, 3) (-37, 4) (-36, 1) (-35, 2) (-34, 4) (-33, 2) (-32, 6) (-31, 9) (-30, 5) (-29, 4) (-28, 6) (-27, 9) (-26, 12) (-25, 8) (-24, 4) (-23, 8) (-22, 14) (-21, 8) (-20, 10) (-19, 13) (-18, 23) (-17, 13) (-16, 25) (-15, 20) (-14, 22) (-13, 20) (-12, 28) (-11, 28) (-10, 40) (-9, 35) (-8, 45) (-7, 48) (-6, 54) (-5, 61) (-4, 68) (-3, 79) (-2, 75) (-1, 77) (0, 94)'::ddsketch;

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
     intermediate AS (SELECT ddsketch_percentile(summary, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS a FROM intermediate_ddsketch),
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
        (SELECT ddsketch_percentile(summary, foo.p) FROM intermediate_ddsketch) AS a,
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
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.99])) AS a,
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
            unnest(ddsketch_percentile(input_data.val, 0.05, 1024, perc.percentiles)) v
        FROM perc, input_data
        GROUP BY perc.percentiles
    ) foo
) bar where v2 > v1;

-- <value,count> API

select trunc_value(ddsketch_percentile(value, count, 0.05, 1024, ARRAY[0.9, 0.95, 0.99])) as percentiles
from (values
  (47325940488,1),
  (15457695432,2),
  (6889790700,3),
  (4188763788,4),
  (2882932224,5),
  (2114815860,6),
  (1615194324,7),
  (2342114568,9),
  (1626471924,11),
  (1660755408,14),
  (1143728292,17),
  (1082582424,21),
  (911488284,26),
  (728863908,32),
  (654898692,40),
  (530198076,50),
  (417883440,62),
  (341452344,77),
  (274579584,95),
  (231921120,118),
  (184091820,146),
  (152469828,181),
  (125634972,224),
  (107059704,278),
  (88746120,345),
  (73135668,428),
  (61035756,531),
  (50683320,658),
  (42331824,816),
  (35234400,1012),
  (29341356,1255),
  (24290928,1556),
  (20284668,1929),
  (17215908,2391),
  (14737488,2964),
  (12692772,3674),
  (11220732,4555),
  (9787584,5647),
  (8148420,7000),
  (6918612,8678),
  (6015000,10758),
  (5480316,13336),
  (5443356,16532),
  (4535616,20494),
  (3962316,25406),
  (3914484,31495),
  (3828108,39043),
  (3583536,48400),
  (4104120,60000),
  (166024740,2147483647)) foo (count, value);

----------------------------------------------
-- nice data set with random data (uniform) --
----------------------------------------------

-- 0.05 alpha
WITH
 data AS (SELECT prng(10000) x, prng(10000, 29823218) cnt),
 data_expanded AS (SELECT x FROM (SELECT x, generate_series(1, (10 + 100 * cnt)::int) FROM data) foo ORDER BY random())
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(a) AS a,
        unnest(b) AS b
    FROM
       (SELECT lower_quantile(1.0 + x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) a FROM data_expanded) foo,
       (SELECT ddsketch_percentile(1.0 + x, (10 + 100 * cnt)::int, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) b FROM data) bar
) baz;

-- 0.01 alpha
WITH
 data AS (SELECT prng(10000) x, prng(10000, 29823218) cnt),
 data_expanded AS (SELECT x FROM (SELECT x, generate_series(1, (10 + 100 * cnt)::int) FROM data) foo ORDER BY random())
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(a) AS a,
        unnest(b) AS b
    FROM
       (SELECT lower_quantile(1.0 + x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) a FROM data_expanded) foo,
       (SELECT ddsketch_percentile(1.0 + x, (10 + 100 * cnt)::int, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) b FROM data) bar
) baz;

-- 0.001 alpha
WITH
 data AS (SELECT prng(10000) x, prng(10000, 29823218) cnt),
 data_expanded AS (SELECT x FROM (SELECT x, generate_series(1, (10 + 100 * cnt)::int) FROM data) foo ORDER BY random())
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(a) AS a,
        unnest(b) AS b
    FROM
       (SELECT lower_quantile(1.0 + x, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) a FROM data_expanded) foo,
       (SELECT ddsketch_percentile(1.0 + x, (10 + 100 * cnt)::int, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) b FROM data) bar
) baz;

-- test incremental API (adding values one by one)
CREATE TABLE t (d ddsketch);
INSERT INTO t VALUES (NULL);

-- check this produces the same result building the ddsketch at once
DO LANGUAGE plpgsql $$
DECLARE
  r RECORD;
BEGIN
    FOR r IN (SELECT i FROM generate_series(1,1000) s(i) ORDER BY md5(i::text)) LOOP
        UPDATE t SET d = ddsketch_add(d, r.i, 0.05, 1024);
    END LOOP;
END$$;

-- compare the results
WITH x AS (SELECT i FROM generate_series(1,1000) s(i) ORDER BY md5(i::text))
SELECT (SELECT ddsketch(d)::text FROM t) = (SELECT ddsketch(x.i, 0.05, 1024)::text FROM x) AS match;


-- now do the same thing, but add values with a count
TRUNCATE t;
INSERT INTO t VALUES (NULL);

-- check this produces the same result building the ddsketch at once
DO LANGUAGE plpgsql $$
DECLARE
  r RECORD;
BEGIN
    FOR r IN (SELECT i AS v, (1 + pow(mod(i,13), 2)::int) AS c FROM generate_series(1,1000) s(i) ORDER BY md5(i::text)) LOOP
        UPDATE t SET d = ddsketch_add(d, r.v, r.c, 0.05, 1024);
    END LOOP;
END$$;

-- compare the results
WITH x AS (SELECT i::double precision AS v, (1 + pow(mod(i,13), 2)::int) AS c FROM generate_series(1,1000) s(i) ORDER BY md5(i::text))
SELECT (SELECT ddsketch(d)::text FROM t) = (SELECT ddsketch(x.v, x.c, 0.05, 1024)::text FROM x) AS match;


-- now try the same thing with bulk incremental update (using arrays)
TRUNCATE t;
INSERT INTO t VALUES (NULL);

DO LANGUAGE plpgsql $$
DECLARE
  r RECORD;
BEGIN
    FOR r IN (SELECT a, array_agg(i::double precision) AS v FROM (SELECT mod(i,5) AS a, i FROM generate_series(1,1000) s(i) ORDER BY mod(i,5), md5(i::text)) foo GROUP BY a ORDER BY a) LOOP
        UPDATE t SET d = ddsketch_add(d, r.v, 0.05, 1024);
    END LOOP;
END$$;

-- compare the results
WITH x AS (SELECT mod(i,5) AS a, i::double precision AS d FROM generate_series(1,1000) s(i) ORDER BY mod(i,5), i)
SELECT (SELECT ddsketch(d)::text FROM t) = (SELECT ddsketch(x.d, 0.05, 1024)::text FROM x);


-- now try the same thing with bulk incremental update (using ddsketches)
TRUNCATE t;
INSERT INTO t VALUES (NULL);

DO LANGUAGE plpgsql $$
DECLARE
  r RECORD;
BEGIN
    FOR r IN (SELECT a, ddsketch(i, 0.05, 1024) AS d FROM (SELECT mod(i,5) AS a, i FROM generate_series(1,1000) s(i) ORDER BY mod(i,5), md5(i::text)) foo GROUP BY a ORDER BY a) LOOP
        UPDATE t SET d = ddsketch_union(d, r.d);
    END LOOP;
END$$;

-- compare the results
WITH x AS (SELECT a, ddsketch(i, 0.05, 1024) AS d FROM (SELECT mod(i,5) AS a, i FROM generate_series(1,1000) s(i) ORDER BY mod(i,5), md5(i::text)) foo GROUP BY a ORDER BY a)
SELECT (SELECT ddsketch(d)::text FROM t) = (SELECT ddsketch(x.d)::text FROM x);

-- percentile_of with an array of values
WITH data AS (SELECT i / 50.0 AS v FROM generate_series(-5000, 5000) s(i))
SELECT
    f,
    abs(a - b) < 0.05
FROM (
    SELECT
      unnest(ARRAY[-100.0, -75.0, -50.0, -25.0, 0.0, 25.0, 50.0, 75.0, 100.0]) f,
      unnest((SELECT ddsketch_percentile_of(data.v, 0.05, 1024, ARRAY[-100.0, -75.0, -50.0, -25.0, 0.0, 25.0, 50.0, 75.0, 100.0]) FROM data)) a,
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
        (SELECT ddsketch_percentile_of(data.v, 0.05, 1024, f) FROM data) a,
        (SELECT percent_rank(f) WITHIN GROUP (ORDER BY v) FROM data) AS b
    FROM unnest(ARRAY[-100.0, -75.0, -50.0, -25.0, 0, 25.0, 50.0, 75.0, 100.0]) AS f
) foo;

-- <value,count> API with percentile_of and individual values
WITH
    data AS (SELECT i / 50.0 AS v, 1 + abs(mod(i,13)) AS c FROM generate_series(-5000, 5000) s(i)),
    data_expanded AS (SELECT foo.v FROM (SELECT data.c, data.v FROM data) foo, LATERAL generate_series(1, c))
SELECT
    f,
    abs(a - b) < 0.05
FROM (
    SELECT
      f,
      (SELECT ddsketch_percentile_of(data.v, data.c, 0.05, 1024, f) FROM data) a,
      (SELECT percent_rank(f) WITHIN GROUP (ORDER BY v) FROM data_expanded) AS b
    FROM unnest(ARRAY[-100.0, -75.0, -50.0, -25.0, 0, 25.0, 50.0, 75.0, 100.0]) AS f
) foo;

-- <value,count> API with percentile_of and an array
WITH
    data AS (SELECT i / 50.0 AS v, 1 + abs(mod(i,13)) AS c FROM generate_series(-5000, 5000) s(i)),
    data_expanded AS (SELECT foo.v FROM (SELECT data.c, data.v FROM data) foo, LATERAL generate_series(1, c))
SELECT
    f,
    abs(a - b) < 0.05
FROM (
    SELECT
      unnest(ARRAY[-100.0, -75.0, -50.0, -25.0, 0.0, 25.0, 50.0, 75.0, 100.0]) f,
      unnest((SELECT ddsketch_percentile_of(data.v, data.c, 0.05, 1024, ARRAY[-100.0, -75.0, -50.0, -25.0, 0.0, 25.0, 50.0, 75.0, 100.0]) FROM data)) a,
      unnest((SELECT array_agg((SELECT percent_rank(f) WITHIN GROUP (ORDER BY v) FROM data_expanded)) foo FROM unnest(ARRAY[-100.0, -75.0, -50.0, -25.0, 0.0, 25.0, 50.0, 75.0, 100.0]) f)) AS b
) foo;

-- hypothetical-set aggregates
--
-- there's no guarantee for relative-errors of hypothetical-aggregates,
-- but for uniform distribution it's fairly close to the relative error
-- the ddsketch was defined with
WITH
  data AS (SELECT i / 10.0 AS v FROM generate_series(1,10000) s(i)), 
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
    (SELECT ddsketch_percentile_of(sketches.s, foo.v) FROM sketches) a,
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
  sketch_values AS (SELECT ddsketch_percentile_of(sketches.s, vals.v) AS v FROM sketches, vals),
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
  sketch_values AS (SELECT ddsketch_percentile_of(data.v, 0.05, 1024, vals.v) AS v FROM data, vals),
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
  sketch_values AS (SELECT ddsketch_percentile_of(data.v, 0.05, 1024, vals.v) AS v FROM data, vals),
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

CREATE TABLE src_data (v double precision);
INSERT INTO src_data SELECT z FROM random_normal(1000000, mean := 0.0, stddev := 0.1, minval := -1.0, maxval := 1.0) s(z);
ANALYZE src_data;

-- with parallelism
EXPLAIN (COSTS OFF) SELECT ddsketch(v, 0.05, 1024) FROM src_data;
SELECT ddsketch(v, 0.05, 1024) FROM src_data;

EXPLAIN (COSTS OFF) SELECT ddsketch(v, 0.05, 1024) FROM src_data;
SELECT trunc_value(ddsketch_percentile(v, 0.05, 1024, 0.9)) FROM src_data;

EXPLAIN (COSTS OFF) SELECT ddsketch(v, 0.05, 1024) FROM src_data;
SELECT ddsketch_percentile_of(v, 0.05, 1024, 0.9) FROM src_data;

-- without  parallelism
SET max_parallel_workers_per_gather = 0;
EXPLAIN (COSTS OFF) SELECT ddsketch(v, 0.05, 1024) FROM src_data;
SELECT ddsketch(v, 0.05, 1024) FROM src_data;

-- NULL handling

-- individual values, individual percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text))
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
    p,
    (SELECT ddsketch_percentile(data.v, 0.01, 1024, p) FROM data) AS a,
    (SELECT lower_quantile(v, p) FROM data) AS b
  FROM
    unnest((SELECT p FROM perc)) p
) foo;

-- <value,count> API, individual percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text)),
  data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c))
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
    p,
    (SELECT ddsketch_percentile(data.v, data.c, 0.01, 1024, p) FROM data) AS a,
    (SELECT lower_quantile(v, p) FROM data_exp) AS b
  FROM
    unnest((SELECT p FROM perc)) p
) foo;

-- <value,count> API, but count is NULL (should be treated as 1), individual percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text))
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
    p,
    (SELECT ddsketch_percentile(data.v, NULL, 0.01, 1024, p) FROM data) AS a,
    (SELECT lower_quantile(v, p) FROM data) AS b
  FROM
    unnest((SELECT p FROM perc)) p
) foo;

-- individual values, array of percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text))
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
  SELECT
    unnest((SELECT p FROM perc)) p,
    unnest((SELECT ddsketch_percentile(data.v, 0.01, 1024, perc.p) FROM data, perc)) AS a,
    unnest((SELECT lower_quantile(v, p) FROM data, perc)) AS b
) foo;

-- <value,count> API, array of percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text)),
  data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c))
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
  SELECT
    unnest((SELECT p FROM perc)) p,
    unnest((SELECT ddsketch_percentile(data.v, data.c, 0.01, 1024, perc.p) FROM data, perc)) AS a,
    unnest((SELECT lower_quantile(v, p) FROM data_exp, perc)) AS b
) foo;

-- <value,count> API, but count is NULL (should be treated as 1), array of percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text))
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
  SELECT
    unnest((SELECT p FROM perc)) p,
    unnest((SELECT ddsketch_percentile(data.v, NULL, 0.01, 1024, perc.p) FROM data, perc)) AS a,
    unnest((SELECT lower_quantile(v, p) FROM data, perc)) AS b
) foo;

-- NULL handling, but this time make sure the first value is NULL

-- individual values, individual percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT NULL AS v, 1 AS c UNION ALL SELECT * FROM (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text)) foo)
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
    p,
    (SELECT ddsketch_percentile(data.v, 0.01, 1024, p) FROM data) AS a,
    (SELECT lower_quantile(v, p) FROM data) AS b
  FROM
    unnest((SELECT p FROM perc)) p
) foo;

-- <value,count> API, individual percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT NULL AS v, 1 AS c UNION ALL SELECT * FROM (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text)) foo),
  data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c))
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
    p,
    (SELECT ddsketch_percentile(data.v, data.c, 0.01, 1024, p) FROM data) AS a,
    (SELECT lower_quantile(v, p) FROM data_exp) AS b
  FROM
    unnest((SELECT p FROM perc)) p
) foo;

-- <value,count> API, but count is NULL (should be treated as 1), individual percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT NULL AS v, 1 AS c UNION ALL SELECT * FROM (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text)) foo)
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
    p,
    (SELECT ddsketch_percentile(data.v, NULL, 0.01, 1024, p) FROM data) AS a,
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
    unnest((SELECT ddsketch_percentile(data.v, 0.01, 1024, perc.p) FROM data, perc)) AS a,
    unnest((SELECT lower_quantile(v, p) FROM data, perc)) AS b
) foo;

-- <value,count> API, array of percentiles
WITH
  perc AS (SELECT array_agg(i/10.0) AS p FROM generate_series(0,10) AS s(i)),
  data AS (SELECT NULL AS v, 1 AS c UNION ALL SELECT * FROM (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text)) foo),
  data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c))
SELECT
    round(p,2) AS p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
  SELECT
    unnest((SELECT p FROM perc)) p,
    unnest((SELECT ddsketch_percentile(data.v, data.c, 0.01, 1024, perc.p) FROM data, perc)) AS a,
    unnest((SELECT lower_quantile(v, p) FROM data_exp, perc)) AS b
) foo;

-- <value,count> API, but count is NULL (should be treated as 1), array of percentiles
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
    unnest((SELECT ddsketch_percentile(data.v, NULL, 0.01, 1024, perc.p) FROM data, perc)) AS a,
    unnest((SELECT lower_quantile(v, p) FROM data, perc)) AS b
) foo;


-- percentile-of values

-- individual values
WITH
    vals AS (SELECT (i*10.0) AS f FROM generate_series(-10,10) s(i)),
    data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text))
SELECT
    f,
    abs(a - b) < 0.05
FROM (
    SELECT
       f,
       (SELECT ddsketch_percentile_of(data.v, 0.01, 1024, f) FROM data) AS a,
       (SELECT percent_rank(f) WITHIN GROUP (ORDER BY v) FROM data WHERE data.v IS NOT NULL) AS b
    FROM vals
) foo;

WITH
    vals AS (SELECT (i*10.0) AS f FROM generate_series(-10,10) s(i)),
    data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text)),
    data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c))
SELECT
    f,
    abs(a - b) < 0.05
FROM (
    SELECT
       f,
       (SELECT ddsketch_percentile_of(data.v, data.c, 0.01, 1024, f) FROM data) AS a,
       (SELECT percent_rank(f) WITHIN GROUP (ORDER BY v) FROM data_exp WHERE data_exp.v IS NOT NULL) AS b
    FROM vals
) foo;

WITH
    vals AS (SELECT (i*10.0) AS f FROM generate_series(-10,10) s(i)),
    data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text))
SELECT
    f,
    abs(a - b) < 0.05
FROM (
    SELECT
       f,
       (SELECT ddsketch_percentile_of(data.v, NULL, 0.01, 1024, f) FROM data) AS a,
       (SELECT percent_rank(f) WITHIN GROUP (ORDER BY v) FROM data WHERE data.v IS NOT NULL) AS b
    FROM vals
) foo;

-- array of values

WITH data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text)),
     vals AS (SELECT array_agg(i) AS v FROM generate_series(-100,100,10) AS s(i)),
     sketch AS (SELECT ddsketch_percentile_of(data.v, data.c, 0.01, 1024, vals.v) s FROM data, vals),
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

WITH data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text)),
     data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c)),
     vals AS (SELECT array_agg(i) AS v FROM generate_series(-100,100,10) AS s(i)),
     sketch AS (SELECT ddsketch_percentile_of(data.v, data.c, 0.01, 1024, vals.v) s FROM data, vals),
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

WITH data AS (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text)),
     vals AS (SELECT array_agg(i) AS v FROM generate_series(-100,100,10) AS s(i)),
     sketch AS (SELECT ddsketch_percentile_of(data.v, NULL, 0.01, 1024, vals.v) s FROM data, vals),
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

WITH data AS (SELECT NULL AS v, 1 AS c UNION ALL SELECT * FROM (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text)) foo),
     vals AS (SELECT array_agg(i) AS v FROM generate_series(-100,100,10) AS s(i)),
     sketch AS (SELECT ddsketch_percentile_of(data.v, data.c, 0.01, 1024, vals.v) s FROM data, vals),
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

WITH data AS (SELECT NULL AS v, 1 AS c UNION ALL SELECT * FROM (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text)) foo),
     data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c)),
     vals AS (SELECT array_agg(i) AS v FROM generate_series(-100,100,10) AS s(i)),
     sketch AS (SELECT ddsketch_percentile_of(data.v, data.c, 0.01, 1024, vals.v) s FROM data, vals),
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

WITH data AS (SELECT NULL AS v, 1 AS c UNION ALL SELECT * FROM (SELECT (CASE WHEN mod(i,2) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v, 1 + mod(i,5) c FROM generate_series(1,10000) s(i) ORDER BY md5(i::text)) foo),
     vals AS (SELECT array_agg(i) AS v FROM generate_series(-100,100,10) AS s(i)),
     sketch AS (SELECT ddsketch_percentile_of(data.v, NULL, 0.01, 1024, vals.v) s FROM data, vals),
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
WITH data AS (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,10000) s(i)),
     sketches AS (SELECT ddsketch(data.v, 0.05, 1024) AS d FROM data GROUP BY c),
     perc AS (SELECT array_agg(i / 100.0) AS p FROM generate_series(0,100,10) s(i))
SELECT
  round(p, 2),
  check_relative_error(a, b, 0.05) AS check_error,
  print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile(sketches.d, p) FROM sketches) AS a,
     (SELECT lower_quantile(v, p) FROM data) AS b
  FROM unnest((SELECT p FROM perc)) p) foo;

WITH data AS (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,10000) s(i)),
     data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c)),
     sketches AS (SELECT ddsketch(data.v, data.c, 0.05, 1024) AS d FROM data GROUP BY c),
     perc AS (SELECT array_agg(i / 100.0) AS p FROM generate_series(0,100,10) s(i))
SELECT
  round(p, 2),
  check_relative_error(a, b, 0.05) AS check_error,
  print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile(sketches.d, p) FROM sketches) AS a,
     (SELECT lower_quantile(v, p) FROM data_exp) AS b
  FROM unnest((SELECT p FROM perc)) p) foo;

-- NULL is equal to count=1
WITH data AS (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,10000) s(i)),
     sketches AS (SELECT ddsketch(data.v, NULL, 0.05, 1024) AS d FROM data GROUP BY c),
     perc AS (SELECT array_agg(i / 100.0) AS p FROM generate_series(0,100,10) s(i))
SELECT
  round(p, 2),
  check_relative_error(a, b, 0.05) AS check_error,
  print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile(sketches.d, p) FROM sketches) AS a,
     (SELECT lower_quantile(v, p) FROM data) AS b
  FROM unnest((SELECT p FROM perc)) p) foo;

-- percentile_of
WITH data AS (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,10000) s(i)),
     sketches AS (SELECT ddsketch(data.v, 0.05, 1024) AS d FROM data GROUP BY c),
     vals AS (SELECT array_agg(i) AS p FROM generate_series(-100,100,10) s(i))
SELECT
  p,
  abs(a - b) < 0.05
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile_of(sketches.d, p) FROM sketches) AS a,
     (SELECT percent_rank(p) WITHIN GROUP (ORDER BY v) FROM data WHERE data.v IS NOT NULL) AS b
  FROM unnest((SELECT p FROM vals)) p) foo;

WITH data AS (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,10000) s(i)),
     data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c)),
     sketches AS (SELECT ddsketch(data.v, data.c, 0.05, 1024) AS d FROM data GROUP BY c),
     vals AS (SELECT array_agg(i) AS p FROM generate_series(-100,100,10) s(i))
SELECT
  p,
  abs(a - b) < 0.05
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile_of(sketches.d, p) FROM sketches) AS a,
     (SELECT percent_rank(p) WITHIN GROUP (ORDER BY v) FROM data_exp WHERE data_exp.v IS NOT NULL) AS b
  FROM unnest((SELECT p FROM vals)) p) foo;

WITH data AS (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,10000) s(i)),
     sketches AS (SELECT ddsketch(data.v, 0.05, 1024) AS d FROM data GROUP BY c),
     vals AS (SELECT array_agg(i) AS p FROM generate_series(-100,100,10) s(i))
SELECT
  p,
  abs(a - b) < 0.05
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile_of(sketches.d, p) FROM sketches) AS a,
     (SELECT percent_rank(p) WITHIN GROUP (ORDER BY v) FROM data WHERE data.v IS NOT NULL) AS b
  FROM unnest((SELECT p FROM vals)) p) foo;

-- sketches (first value is NULL)
WITH data AS (SELECT 1 AS c, NULL AS v UNION ALL SELECT * FROM (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,10000) s(i)) foo),
     sketches AS (SELECT ddsketch(data.v, 0.05, 1024) AS d FROM data GROUP BY c),
     perc AS (SELECT array_agg(i / 100.0) AS p FROM generate_series(0,100,10) s(i))
SELECT
  round(p, 2),
  check_relative_error(a, b, 0.05) AS check_error,
  print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile(sketches.d, p) FROM sketches) AS a,
     (SELECT lower_quantile(v, p) FROM data) AS b
  FROM unnest((SELECT p FROM perc)) p) foo;

WITH data AS (SELECT 1 AS c, NULL AS v UNION ALL SELECT * FROM (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,10000) s(i)) foo),
     data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c)),
     sketches AS (SELECT ddsketch(data.v, data.c, 0.05, 1024) AS d FROM data GROUP BY c),
     perc AS (SELECT array_agg(i / 100.0) AS p FROM generate_series(0,100,10) s(i))
SELECT
  round(p, 2),
  check_relative_error(a, b, 0.05) AS check_error,
  print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile(sketches.d, p) FROM sketches) AS a,
     (SELECT lower_quantile(v, p) FROM data_exp) AS b
  FROM unnest((SELECT p FROM perc)) p) foo;

WITH data AS (SELECT 1 AS c, NULL AS v UNION ALL SELECT * FROM (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,10000) s(i)) foo),
     sketches AS (SELECT ddsketch(data.v, NULL, 0.05, 1024) AS d FROM data GROUP BY c),
     perc AS (SELECT array_agg(i / 100.0) AS p FROM generate_series(0,100,10) s(i))
SELECT
  round(p, 2),
  check_relative_error(a, b, 0.05) AS check_error,
  print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile(sketches.d, p) FROM sketches) AS a,
     (SELECT lower_quantile(v, p) FROM data) AS b
  FROM unnest((SELECT p FROM perc)) p) foo;

-- NULL is equal to count=1
WITH data AS (SELECT 1 AS c, NULL AS v UNION ALL SELECT * FROM (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,10000) s(i)) foo),
     sketches AS (SELECT ddsketch(data.v, NULL, 0.05, 1024) AS d FROM data GROUP BY c),
     perc AS (SELECT array_agg(i / 100.0) AS p FROM generate_series(0,100,10) s(i))
SELECT
  round(p, 2),
  check_relative_error(a, b, 0.05) AS check_error,
  print_relative_error(a, b, 0.05) AS error_info
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile(sketches.d, p) FROM sketches) AS a,
     (SELECT lower_quantile(v, p) FROM data) AS b
  FROM unnest((SELECT p FROM perc)) p) foo;

-- percentile_of
WITH data AS (SELECT 1 AS c, NULL AS v UNION ALL SELECT * FROM (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,10000) s(i)) foo),
     sketches AS (SELECT ddsketch(data.v, 0.05, 1024) AS d FROM data GROUP BY c),
     vals AS (SELECT array_agg(i) AS p FROM generate_series(-100,100,10) s(i))
SELECT
  p,
  abs(a - b) < 0.05
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile_of(sketches.d, p) FROM sketches) AS a,
     (SELECT percent_rank(p) WITHIN GROUP (ORDER BY v) FROM data WHERE data.v IS NOT NULL) AS b
  FROM unnest((SELECT p FROM vals)) p) foo;

WITH data AS (SELECT 1 AS c, NULL AS v UNION ALL SELECT * FROM (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,10000) s(i)) foo),
     data_exp AS (SELECT data.* FROM data, lateral generate_series(1, data.c)),
     sketches AS (SELECT ddsketch(data.v, data.c, 0.05, 1024) AS d FROM data GROUP BY c),
     vals AS (SELECT array_agg(i) AS p FROM generate_series(-100,100,10) s(i))
SELECT
  p,
  abs(a - b) < 0.05
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile_of(sketches.d, p) FROM sketches) AS a,
     (SELECT percent_rank(p) WITHIN GROUP (ORDER BY v) FROM data_exp WHERE data_exp.v IS NOT NULL) AS b
  FROM unnest((SELECT p FROM vals)) p) foo;

WITH data AS (SELECT 1 AS c, NULL AS v UNION ALL SELECT * FROM (SELECT mod(i,10) AS c, (CASE WHEN mod(i,5) = 0 THEN NULL ELSE (i / 50.0 - 100.0) END) AS v FROM generate_series(1,10000) s(i)) foo),
     sketches AS (SELECT ddsketch(data.v, 0.05, 1024) AS d FROM data GROUP BY c),
     vals AS (SELECT array_agg(i) AS p FROM generate_series(-100,100,10) s(i))
SELECT
  p,
  abs(a - b) < 0.05
FROM (
  SELECT
     p,
     (SELECT ddsketch_percentile_of(sketches.d, p) FROM sketches) AS a,
     (SELECT percent_rank(p) WITHIN GROUP (ORDER BY v) FROM data WHERE data.v IS NOT NULL) AS b
  FROM unnest((SELECT p FROM vals)) p) foo;



SELECT ddsketch_percentile(NULL::ddsketch, 0.9);
SELECT ddsketch_percentile(NULL::ddsketch, ARRAY[0.5, 0.9]);
SELECT ddsketch_percentile_of(NULL::ddsketch, 0.9);
SELECT ddsketch_percentile_of(NULL::ddsketch, ARRAY[0.1, 0.9]);

SELECT ddsketch(NULL::ddsketch) FROM generate_series(1,10);

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
SELECT ddsketch_percentile(sketches.s, 0.9) FROM sketches;

WITH sketches AS (
    SELECT ddsketch(i/100.0, 0.01, 1024) AS s FROM generate_series(1,10000) s(i)
    UNION ALL
    SELECT ddsketch(i/100.0, 0.05, 1024) AS s FROM generate_series(1,10000) s(i)
)
SELECT ddsketch_percentile(sketches.s, ARRAY[0.95, 0.99]) FROM sketches;

WITH sketches AS (
    SELECT ddsketch(i/100.0, 0.01, 1024) AS s FROM generate_series(1,10000) s(i)
    UNION ALL
    SELECT ddsketch(i/100.0, 0.05, 1024) AS s FROM generate_series(1,10000) s(i)
)
SELECT ddsketch_percentile_of(sketches.s, 95) FROM sketches;

WITH sketches AS (
    SELECT ddsketch(i/100.0, 0.01, 1024) AS s FROM generate_series(1,10000) s(i)
    UNION ALL
    SELECT ddsketch(i/100.0, 0.05, 1024) AS s FROM generate_series(1,10000) s(i)
)
SELECT ddsketch_percentile_of(sketches.s, ARRAY[95, 99]) FROM sketches;


WITH sketches AS (
    SELECT ddsketch(i/100.0, 0.01, 1024) AS s FROM generate_series(1,10000) s(i)
)
SELECT ddsketch_union(sketches.s, NULL::ddsketch) FROM sketches;

WITH sketches AS (
    SELECT ddsketch(i/100.0, 0.05, 1024) AS s FROM generate_series(1,10000) s(i)
)
SELECT ddsketch_union(NULL::ddsketch, sketches.s) FROM sketches;

SELECT ddsketch_union(NULL::ddsketch, NULL::ddsketch);


-- test info functions
WITH data AS (SELECT 0.5 + i/10000.0 AS v FROM generate_series(1,10000) s(i)),
     sketch AS (SELECT ddsketch(v, 0.05, 128) s FROM data)
SELECT
  bytes,
  flags,
  alpha,
  count,
  zero_count,
  max_buckets,
  negative_buckets,
  positive_buckets,
  trunc_value(min_indexable) as min_indexable,
  trunc_value(max_indexable) as max_indexable
FROM ddsketch_info((SELECT s FROM sketch));

WITH data AS (SELECT 0.5 + i/10000.0 AS v FROM generate_series(1,10000) s(i)),
     sketch AS (SELECT ddsketch(v, 0.05, 128) s FROM data)
SELECT
  index,
  bucket_index,
  trunc_value(bucket_lower) as bucket_lower,
  trunc_value(bucket_upper) as bucket_upper,
  trunc_value(bucket_length) as bucket_length,
  bucket_count
FROM ddsketch_buckets((SELECT s FROM sketch));

SELECT
  trunc_value(min_indexable) as min_indexable,
  trunc_value(max_indexable) as max_indexable
FROM ddsketch_info(0.05);

SELECT
  index,
  bucket_index,
  trunc_value(bucket_min) as bucket_min,
  trunc_value(bucket_min) as bucket_max
FROM ddsketch_buckets(0.05, 0.5, 5);

SELECT
  index,
  bucket_index,
  trunc_value(bucket_min) as bucket_min,
  trunc_value(bucket_min) as bucket_max
FROM ddsketch_buckets(0.05, -5, -0.5);

CREATE TABLE random_data (v double precision, i int);
INSERT INTO random_data SELECT prng(10000, 45547, 34471541, 3, 1000000), generate_series(1,10000);

-- trimmed aggregates
SELECT ddsketch_sum(1000 * v, 0.05, 1024, 0.0, 1.0)   BETWEEN 4750000 AND 5250000 FROM random_data;
SELECT ddsketch_sum(1000 * v, 0.05, 1024, 0.0, 0.5)   BETWEEN 1200000 AND 1300000 FROM random_data;
SELECT ddsketch_sum(1000 * v, 0.05, 1024, 0.5, 1.0)   BETWEEN 3700000 AND 3800000 FROM random_data;
SELECT ddsketch_sum(1000 * v, 0.05, 1024, 0.25, 0.75) BETWEEN 2450000 AND 2550000 FROM random_data;

SELECT ddsketch_avg(1000 * v, 0.05, 1024, 0.0, 1.0)   BETWEEN 490 AND 510 FROM random_data;
SELECT ddsketch_avg(1000 * v, 0.05, 1024, 0.0, 0.5)   BETWEEN 240 AND 260 FROM random_data;
SELECT ddsketch_avg(1000 * v, 0.05, 1024, 0.5, 1.0)   BETWEEN 740 AND 760 FROM random_data;
SELECT ddsketch_avg(1000 * v, 0.05, 1024, 0.25, 0.75) BETWEEN 490 AND 510 FROM random_data;

-- trimmed aggregates, <value, count> API
SELECT ddsketch_sum(1000 * v, 1 + mod(i,5), 0.05, 1024, 0.0, 1.0)   BETWEEN 3 * 4750000 AND 3 * 5250000 FROM random_data;
SELECT ddsketch_sum(1000 * v, 1 + mod(i,5), 0.05, 1024, 0.0, 0.5)   BETWEEN 3 * 1200000 AND 3 * 1300000 FROM random_data;
SELECT ddsketch_sum(1000 * v, 1 + mod(i,5), 0.05, 1024, 0.5, 1.0)   BETWEEN 3 * 3700000 AND 3 * 3800000 FROM random_data;
SELECT ddsketch_sum(1000 * v, 1 + mod(i,5), 0.05, 1024, 0.25, 0.75) BETWEEN 3 * 2450000 AND 3 * 2550000 FROM random_data;

SELECT ddsketch_avg(1000 * v, 1 + mod(i,5), 0.05, 1024, 0.0, 1.0)   BETWEEN 490 AND 510 FROM random_data;
SELECT ddsketch_avg(1000 * v, 1 + mod(i,5), 0.05, 1024, 0.0, 0.5)   BETWEEN 240 AND 260 FROM random_data;
SELECT ddsketch_avg(1000 * v, 1 + mod(i,5), 0.05, 1024, 0.5, 1.0)   BETWEEN 740 AND 760 FROM random_data;
SELECT ddsketch_avg(1000 * v, 1 + mod(i,5), 0.05, 1024, 0.25, 0.75) BETWEEN 490 AND 510 FROM random_data;

-- trimmed aggregates when aggregating sketches
WITH sketches AS (SELECT mod(i,5) AS c, ddsketch(1000 * v, 0.05, 1024) AS s FROM random_data GROUP BY 1)
SELECT ddsketch_sum(s, 0.0, 1.0)   BETWEEN 4750000 AND 5250000 FROM sketches;

WITH sketches AS (SELECT mod(i,5) AS c, ddsketch(1000 * v, 0.05, 1024) AS s FROM random_data GROUP BY 1)
SELECT ddsketch_sum(s, 0.0, 0.5)   BETWEEN 1200000 AND 1300000 FROM sketches;

WITH sketches AS (SELECT mod(i,5) AS c, ddsketch(1000 * v, 0.05, 1024) AS s FROM random_data GROUP BY 1)
SELECT ddsketch_sum(s, 0.5, 1.0)   BETWEEN 3700000 AND 3800000 FROM sketches;

WITH sketches AS (SELECT mod(i,5) AS c, ddsketch(1000 * v, 0.05, 1024) AS s FROM random_data GROUP BY 1)
SELECT ddsketch_sum(s, 0.25, 0.75) BETWEEN 2450000 AND 2550000 FROM sketches;

WITH sketches AS (SELECT mod(i,5) AS c, ddsketch(1000 * v, 0.05, 1024) AS s FROM random_data GROUP BY 1)
SELECT ddsketch_avg(s, 0.0, 1.0)   BETWEEN 490 AND 510 FROM sketches;

WITH sketches AS (SELECT mod(i,5) AS c, ddsketch(1000 * v, 0.05, 1024) AS s FROM random_data GROUP BY 1)
SELECT ddsketch_avg(s, 0.0, 0.5)   BETWEEN 240 AND 260 FROM sketches;

WITH sketches AS (SELECT mod(i,5) AS c, ddsketch(1000 * v, 0.05, 1024) AS s FROM random_data GROUP BY 1)
SELECT ddsketch_avg(s, 0.5, 1.0)   BETWEEN 740 AND 760 FROM sketches;

WITH sketches AS (SELECT mod(i,5) AS c, ddsketch(1000 * v, 0.05, 1024) AS s FROM random_data GROUP BY 1)
SELECT ddsketch_avg(s, 0.25, 0.75) BETWEEN 490 AND 510 FROM sketches;

-- trimmed aggregates for a single sketch
SELECT ddsketch_sketch_sum(ddsketch(1000 * v, 0.05, 1024), 0.0, 1.0)   BETWEEN 4750000 AND 5250000 FROM random_data;
SELECT ddsketch_sketch_sum(ddsketch(1000 * v, 0.05, 1024), 0.0, 0.5)   BETWEEN 1200000 AND 1300000 FROM random_data;
SELECT ddsketch_sketch_sum(ddsketch(1000 * v, 0.05, 1024), 0.5, 1.0)   BETWEEN 3700000 AND 3800000 FROM random_data;
SELECT ddsketch_sketch_sum(ddsketch(1000 * v, 0.05, 1024), 0.25, 0.75) BETWEEN 2450000 AND 2550000 FROM random_data;

SELECT ddsketch_sketch_avg(ddsketch(1000 * v, 0.05, 1024), 0.0, 1.0)   BETWEEN 490 AND 510 FROM random_data;
SELECT ddsketch_sketch_avg(ddsketch(1000 * v, 0.05, 1024), 0.0, 0.5)   BETWEEN 240 AND 260 FROM random_data;
SELECT ddsketch_sketch_avg(ddsketch(1000 * v, 0.05, 1024), 0.5, 1.0)   BETWEEN 740 AND 760 FROM random_data;
SELECT ddsketch_sketch_avg(ddsketch(1000 * v, 0.05, 1024), 0.25, 0.75) BETWEEN 490 AND 510 FROM random_data;

-- check trim parameters
SELECT ddsketch_sum(1000 * v, 0.05, 1024, -0.1, 1.0)   BETWEEN 4750000 AND 5250000 FROM random_data;
SELECT ddsketch_sum(1000 * v, 0.05, 1024, 0.1, 1.1)   BETWEEN 4750000 AND 5250000 FROM random_data;
SELECT ddsketch_sum(1000 * v, 0.05, 1024, 0.9, 0.1)   BETWEEN 4750000 AND 5250000 FROM random_data;
SELECT ddsketch_sum(1000 * v, 0.05, 1024, 0.5, 0.5)   BETWEEN 4750000 AND 5250000 FROM random_data;
