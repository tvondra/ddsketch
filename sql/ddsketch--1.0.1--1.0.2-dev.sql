DROP AGGREGATE ddsketch_percentile(ddsketch, double precision);
DROP AGGREGATE ddsketch_percentile(ddsketch, double precision[]);
DROP AGGREGATE ddsketch_percentile(double precision, bigint, double precision, integer, double precision);
DROP AGGREGATE ddsketch_percentile(double precision, bigint, double precision, integer, double precision[]);
DROP AGGREGATE ddsketch_percentile(double precision, double precision, integer, double precision);
DROP AGGREGATE ddsketch_percentile(double precision, double precision, integer, double precision[]);
DROP AGGREGATE ddsketch_percentile_of(ddsketch, double precision);
DROP AGGREGATE ddsketch_percentile_of(ddsketch, double precision[]);
DROP AGGREGATE ddsketch_percentile_of(double precision, bigint, double precision, integer, double precision);
DROP AGGREGATE ddsketch_percentile_of(double precision, bigint, double precision, integer, double precision[]);
DROP AGGREGATE ddsketch_percentile_of(double precision, double precision, integer, double precision);
DROP AGGREGATE ddsketch_percentile_of(double precision, double precision, integer, double precision[]);
DROP AGGREGATE ddsketch_sum(ddsketch, double precision, double precision);
DROP AGGREGATE ddsketch_sum(double precision, bigint, double precision, integer, double precision, double precision);
DROP AGGREGATE ddsketch_sum(double precision, double precision, integer, double precision, double precision);
DROP AGGREGATE ddsketch_avg(ddsketch, double precision, double precision);
DROP AGGREGATE ddsketch_avg(double precision, bigint, double precision, integer, double precision, double precision);
DROP AGGREGATE ddsketch_avg(double precision, double precision, integer, double precision, double precision);

DROP FUNCTION ddsketch_add_double                    (internal, double precision, double precision, integer, double precision);
DROP FUNCTION ddsketch_add_double_array              (internal, double precision, double precision, integer, double precision[]);
DROP FUNCTION ddsketch_add_double_array_count        (internal, double precision, bigint, double precision, integer, double precision[]);
DROP FUNCTION ddsketch_add_double_array_values       (internal, double precision, double precision, integer, double precision[]);
DROP FUNCTION ddsketch_add_double_array_values_count (internal, double precision, bigint, double precision, integer, double precision[]);
DROP FUNCTION ddsketch_add_double_count              (internal, double precision, bigint, double precision, integer, double precision);
DROP FUNCTION ddsketch_add_double_count_trimmed      (internal, double precision, bigint, double precision, integer, double precision, double precision);
DROP FUNCTION ddsketch_add_double_trimmed            (internal, double precision, double precision, integer, double precision, double precision);
DROP FUNCTION ddsketch_add_double_values             (internal, double precision, double precision, integer, double precision);
DROP FUNCTION ddsketch_add_double_values_count       (internal, double precision, bigint, double precision, integer, double precision);
DROP FUNCTION ddsketch_add_sketch                    (internal, ddsketch, double precision);
DROP FUNCTION ddsketch_add_sketch_array              (internal, ddsketch, double precision[]);
DROP FUNCTION ddsketch_add_sketch_array_values       (internal, ddsketch, double precision[]);
DROP FUNCTION ddsketch_add_sketch_trimmed            (internal, ddsketch, double precision, double precision);
DROP FUNCTION ddsketch_add_sketch_values             (internal, ddsketch, double precision);
DROP FUNCTION ddsketch_array_percentiles             (internal);
DROP FUNCTION ddsketch_array_percentiles_of          (internal);
DROP FUNCTION ddsketch_percentiles                   (internal);
DROP FUNCTION ddsketch_percentiles_of                (internal);
DROP FUNCTION ddsketch_trimmed_avg                   (internal);
DROP FUNCTION ddsketch_trimmed_sum                   (internal);

CREATE FUNCTION ddsketch_percentile(ddsketch, double precision)
    RETURNS double precision
    AS 'ddsketch', 'ddsketch_percentiles'
    LANGUAGE C IMMUTABLE PARALLEL SAFE;

CREATE FUNCTION ddsketch_percentile(ddsketch, double precision[])
    RETURNS double precision[]
    AS 'ddsketch', 'ddsketch_array_percentiles'
    LANGUAGE C IMMUTABLE PARALLEL SAFE;

CREATE FUNCTION ddsketch_percentile_of(ddsketch, double precision)
    RETURNS double precision
    AS 'ddsketch', 'ddsketch_percentiles_of'
    LANGUAGE C IMMUTABLE PARALLEL SAFE;

CREATE FUNCTION ddsketch_percentile_of(ddsketch, double precision[])
    RETURNS double precision[]
    AS 'ddsketch', 'ddsketch_array_percentiles_of'
    LANGUAGE C IMMUTABLE PARALLEL SAFE;

DROP FUNCTION ddsketch_sketch_avg                    (ddsketch, double precision, double precision);
DROP FUNCTION ddsketch_sketch_sum                    (ddsketch, double precision, double precision);

CREATE OR REPLACE FUNCTION ddsketch_sum(p_sketch ddsketch, p_low double precision = 0.0, p_high double precision = 1.0)
    RETURNS double precision
    AS 'ddsketch', 'ddsketch_sketch_sum'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION ddsketch_avg(p_sketch ddsketch, p_low double precision = 0.0, p_high double precision = 1.0)
    RETURNS double precision
    AS 'ddsketch', 'ddsketch_sketch_avg'
    LANGUAGE C IMMUTABLE STRICT;
