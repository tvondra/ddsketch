\set ECHO none

-- disable the notices for the create script (shell types etc.)
SET client_min_messages = 'WARNING';
\i ddsketch--1.0.0.sql
SET client_min_messages = 'NOTICE';

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

    err := relative_error(estimated_value, actual_value);

    RETURN (err < allowed_error);

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

    IF err < allowed_error THEN
        RETURN NULL;
    END IF;

    RETURN format('estimate = %s, actual = %s, error = %s', estimated_value, actual_value, err);

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
WITH data AS (SELECT i AS x FROM generate_series(1,1000000) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM generate_series(1,1000000) s(i)),
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
WITH data AS (SELECT i AS x FROM generate_series(1,1000000) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM generate_series(1,1000000) s(i)),
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
WITH data AS (SELECT i AS x FROM generate_series(1,1000000) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM generate_series(1,1000000) s(i)),
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

------------------------------------------------------------
-- nice data set with ordered (desc) / evenly-spaced data --
------------------------------------------------------------

-- 0.05 alpha
WITH data AS (SELECT i AS x FROM generate_series(1000000,1,-1) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM generate_series(1000000,1,-1) s(i)),
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
WITH data AS (SELECT i AS x FROM generate_series(1000000,1,-1) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM generate_series(1000000,1,-1) s(i)),
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
WITH data AS (SELECT i AS x FROM generate_series(1000000,1,-1) s(i))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM generate_series(1000000,1,-1) s(i)),
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

----------------------------------------------------
-- nice data set with random / evenly-spaced data --
----------------------------------------------------

-- 10 centroids (tiny)
WITH data AS (SELECT i AS x FROM (SELECT generate_series(1,1000000) AS i, prng(1000000, 49979693) AS x ORDER BY x) foo)
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM (SELECT generate_series(1,1000000) AS i, prng(1000000, 49979693) AS x ORDER BY x) foo),
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

-- 100 centroids (okay-ish)
WITH data AS (SELECT i AS x FROM (SELECT generate_series(1,1000000) AS i, prng(1000000, 49979693) AS x ORDER BY x) foo)
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM (SELECT generate_series(1,1000000) AS i, prng(1000000, 49979693) AS x ORDER BY x) foo),
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

-- 1000 centroids (very accurate)
WITH data AS (SELECT i AS x FROM (SELECT generate_series(1,1000000) AS i, prng(1000000, 49979693) AS x ORDER BY x) foo)
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT i AS x FROM (SELECT generate_series(1,1000000) AS i, prng(1000000, 49979693) AS x ORDER BY x) foo),
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

-- 10 centroids (tiny)
WITH data AS (SELECT 1.0 + x AS x FROM prng(1000000) s(x))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + x AS x FROM prng(1000000) s(x)),
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

-- 100 centroids (okay-ish)
WITH data AS (SELECT 1.0 + x AS x FROM prng(1000000) s(x))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + x AS x FROM prng(1000000) s(x)),
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

-- 1000 centroids (very accurate)
WITH data AS (SELECT 1.0 + x AS x FROM prng(1000000) s(x))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + x AS x FROM prng(1000000) s(x)),
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

-- 10 centroids (tiny)
WITH data AS (SELECT 1.0 + sqrt(z) AS x FROM prng(1000000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + sqrt(z) AS x FROM prng(1000000) s(z)),
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

-- 100 centroids (okay-ish)
WITH data AS (SELECT 1.0 + sqrt(z) AS x FROM prng(1000000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + sqrt(z) AS x FROM prng(1000000) s(z)),
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

-- 1000 centroids (very accurate)
WITH data AS (SELECT 1.0 + sqrt(z) AS x FROM prng(1000000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + sqrt(z) AS x FROM prng(1000000) s(z)),
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

-- 10 centroids (tiny)
WITH data AS (SELECT 1.0 + sqrt(sqrt(z)) AS x FROM prng(1000000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + sqrt(sqrt(z)) AS x FROM prng(1000000) s(z)),
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

-- 100 centroids (okay-ish)
WITH data AS (SELECT 1.0 + sqrt(sqrt(z)) AS x FROM prng(1000000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + sqrt(sqrt(z)) AS x FROM prng(1000000) s(z)),
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

-- 1000 centroids (very accurate)
WITH data AS (SELECT 1.0 + sqrt(sqrt(z)) AS x FROM prng(1000000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + sqrt(sqrt(z)) AS x FROM prng(1000000) s(z)),
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

-- 10 centroids (tiny)
WITH data AS (SELECT 1.0 + pow(z, 2) AS x FROM prng(1000000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + pow(z, 2) AS x FROM prng(1000000) s(z)),
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

-- 100 centroids (okay-ish)
WITH data AS (SELECT 1.0 + pow(z, 2) AS x FROM prng(1000000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + pow(z, 2) AS x FROM prng(1000000) s(z)),
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

-- 1000 centroids (very accurate)
WITH data AS (SELECT 1.0 + pow(z, 2) AS x FROM prng(1000000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + pow(z, 2) AS x FROM prng(1000000) s(z)),
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

-- 10 centroids (tiny)
WITH data AS (SELECT 1.0 + pow(z, 4) AS x FROM prng(1000000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + pow(z, 4) AS x FROM prng(1000000) s(z)),
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

-- 100 centroids (okay-ish)
WITH data AS (SELECT 1.0 + pow(z, 4) AS x FROM prng(1000000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + pow(z, 4) AS x FROM prng(1000000) s(z)),
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

-- 1000 centroids (very accurate)
WITH data AS (SELECT 1.0 + pow(z, 4) AS x FROM prng(1000000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + pow(z, 4) AS x FROM prng(1000000) s(z)),
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

-- 10 centroids (tiny)
WITH data AS (SELECT 1.0 + pow(z, 4) AS x FROM random_normal(1000000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.05) AS check_error,
    print_relative_error(a, b, 0.05) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + pow(z, 4) AS x FROM random_normal(1000000) s(z)),
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

-- 100 centroids (okay-ish)
WITH data AS (SELECT 1.0 + pow(z, 4) AS x FROM random_normal(1000000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.01) AS check_error,
    print_relative_error(a, b, 0.01) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + pow(z, 4) AS x FROM random_normal(1000000) s(z)),
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

-- 1000 centroids (very accurate)
WITH data AS (SELECT 1.0 + pow(z, 4) AS x FROM random_normal(1000000) s(z))
SELECT
    p,
    check_relative_error(a, b, 0.001) AS check_error,
    print_relative_error(a, b, 0.001) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.001, 8192, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- make sure the resulting percentiles are in the right order
WITH data AS (SELECT 1.0 + pow(z, 4) AS x FROM random_normal(1000000) s(z)),
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

-- some basic tests to verify transforming from and to text work
-- 10 centroids (tiny)
WITH data AS (SELECT i AS x FROM generate_series(1,1000000) s(i)),
     intermediate AS (SELECT ddsketch(x, 0.05, 1024)::text AS intermediate_x FROM data),
     ddsketch_parsed AS (SELECT ddsketch_percentile(intermediate_x::ddsketch, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS a FROM intermediate),
     pg_percentile AS (SELECT percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x) AS b FROM data)
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

-- verify we can store ddsketch in a summary table
CREATE TABLE intermediate_ddsketch (grouping int, summary ddsketch);

WITH data AS (SELECT row_number() OVER () AS i, 1.0 + pow(z, 4) AS x FROM random_normal(1000000) s(z))
INSERT INTO intermediate_ddsketch
SELECT
    i % 10 AS grouping,
    ddsketch(x, 0.05, 1024) AS summary
FROM data
GROUP BY i % 10;

WITH data AS (SELECT 1.0 + pow(z, 4) AS x FROM random_normal(1000000) s(z)),
     intermediate AS (SELECT ddsketch_percentile(summary, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) AS a FROM intermediate_ddsketch),
     pg_percentile AS (SELECT percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY x) AS b FROM data)
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

-- verify 'extreme' percentiles for the dataset would not read out of bounds on the centroids
WITH data AS (SELECT x FROM generate_series(1,10) AS x)
SELECT
    p,
    check_relative_error(a, b, 0.055) AS check_error,
    print_relative_error(a, b, 0.055) AS error_info
FROM (
    SELECT
        unnest(ARRAY[0.01, 0.99]) AS p,
        unnest(ddsketch_percentile(x, 0.05, 1024, ARRAY[0.01, 0.99])) AS a,
        unnest(percentile_disc(ARRAY[0.01, 0.99]) WITHIN GROUP (ORDER BY x)) AS b
    FROM data
) foo;

-- check that the computed percentiles are perfectly correlated (don't decrease for higher p values)
-- first test on a tiny ddsketch with all centroids having count = 1
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

select ddsketch_percentile(value, count, 0.05, 1024, ARRAY[0.9, 0.95, 0.99])
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

-- 10 centroids (tiny)
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
       (SELECT percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY 1.0 + x) a FROM data_expanded) foo,
       (SELECT ddsketch_percentile(1.0 + x, (10 + 100 * cnt)::int, 0.05, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) b FROM data) bar
) baz;

-- 100 centroids (okay-ish)
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
       (SELECT percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY 1.0 + x) a FROM data_expanded) foo,
       (SELECT ddsketch_percentile(1.0 + x, (10 + 100 * cnt)::int, 0.01, 1024, ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) b FROM data) bar
) baz;

-- 1000 centroids (very accurate)
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
       (SELECT percentile_disc(ARRAY[0.01, 0.05, 0.1, 0.9, 0.95, 0.99]) WITHIN GROUP (ORDER BY 1.0 + x) a FROM data_expanded) foo,
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

-- compare the results, but do force a compaction of the incremental result
WITH x AS (SELECT i FROM generate_series(1,1000) s(i) ORDER BY md5(i::text))
SELECT (SELECT ddsketch(d)::text FROM t) = (SELECT ddsketch(x.i, 0.05, 1024)::text FROM x) AS match;

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

-- compare the results, but do force a compaction of the incremental result
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

-- compare the results, but do force a compaction of the incremental result
WITH x AS (SELECT a, ddsketch(i, 0.05, 1024) AS d FROM (SELECT mod(i,5) AS a, i FROM generate_series(1,1000) s(i) ORDER BY mod(i,5), md5(i::text)) foo GROUP BY a ORDER BY a)
SELECT (SELECT ddsketch(d)::text FROM t) = (SELECT ddsketch(x.d)::text FROM x);
