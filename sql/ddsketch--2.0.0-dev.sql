CREATE TYPE ddsketch;

CREATE OR REPLACE FUNCTION ddsketch_in(cstring)
    RETURNS ddsketch
    AS 'ddsketch', 'ddsketch_in'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION ddsketch_out(ddsketch)
    RETURNS cstring
    AS 'ddsketch', 'ddsketch_out'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION ddsketch_send(ddsketch)
    RETURNS bytea
    AS 'ddsketch', 'ddsketch_send'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION ddsketch_recv(internal)
    RETURNS ddsketch
    AS 'ddsketch', 'ddsketch_recv'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE TYPE ddsketch (
    INPUT = ddsketch_in,
    OUTPUT = ddsketch_out,
    RECEIVE = ddsketch_recv,
    SEND = ddsketch_send,
    INTERNALLENGTH = variable,
    STORAGE = extended,
    ALIGNMENT = double
);

CREATE OR REPLACE FUNCTION ddsketch_sketch(internal)
    RETURNS ddsketch
    AS 'ddsketch', 'ddsketch_sketch'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION ddsketch_add_double(internal, double precision, double precision, int)
    RETURNS internal
    AS 'ddsketch', 'ddsketch_add_double'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION ddsketch_combine(internal, internal)
    RETURNS internal
    AS 'ddsketch', 'ddsketch_combine'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION ddsketch_serial(internal)
    RETURNS bytea
    AS 'ddsketch', 'ddsketch_serial'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION ddsketch_deserial(bytea, internal)
    RETURNS internal
    AS 'ddsketch', 'ddsketch_deserial'
    LANGUAGE C IMMUTABLE STRICT;

CREATE AGGREGATE ddsketch(value double precision, alpha double precision, max_buckets int) (
    SFUNC = ddsketch_add_double,
    STYPE = internal,
    FINALFUNC = ddsketch_sketch,
    SERIALFUNC = ddsketch_serial,
    DESERIALFUNC = ddsketch_deserial,
    COMBINEFUNC = ddsketch_combine,
    PARALLEL = SAFE
);

CREATE OR REPLACE FUNCTION ddsketch_add_sketch(internal, ddsketch)
    RETURNS internal
    AS 'ddsketch', 'ddsketch_add_sketch'
    LANGUAGE C IMMUTABLE;

CREATE AGGREGATE ddsketch(sketch ddsketch) (
    SFUNC = ddsketch_add_sketch,
    STYPE = internal,
    FINALFUNC = ddsketch_sketch,
    SERIALFUNC = ddsketch_serial,
    DESERIALFUNC = ddsketch_deserial,
    COMBINEFUNC = ddsketch_combine,
    PARALLEL = SAFE
);

CREATE OR REPLACE FUNCTION ddsketch_count(sketch ddsketch)
    RETURNS bigint
    AS 'ddsketch', 'ddsketch_count'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;


CREATE OR REPLACE FUNCTION ddsketch_add_double_count(internal, double precision, bigint, double precision, int)
    RETURNS internal
    AS 'ddsketch', 'ddsketch_add_double_count'
    LANGUAGE C IMMUTABLE;

CREATE AGGREGATE ddsketch(value double precision, count bigint, alpha double precision, max_buckets int) (
    SFUNC = ddsketch_add_double_count,
    STYPE = internal,
    FINALFUNC = ddsketch_sketch,
    SERIALFUNC = ddsketch_serial,
    DESERIALFUNC = ddsketch_deserial,
    COMBINEFUNC = ddsketch_combine,
    PARALLEL = SAFE
);

CREATE OR REPLACE FUNCTION ddsketch_add(sketch ddsketch, "value" double precision, alpha double precision = NULL, max_buckets int = NULL)
    RETURNS ddsketch
    AS 'ddsketch', 'ddsketch_add_double_increment'
    LANGUAGE C IMMUTABLE PARALLEL SAFE;

CREATE OR REPLACE FUNCTION ddsketch_add(sketch ddsketch, "value" double precision, count bigint, alpha double precision = NULL, max_buckets int = NULL)
    RETURNS ddsketch
    AS 'ddsketch', 'ddsketch_add_double_count_increment'
    LANGUAGE C IMMUTABLE PARALLEL SAFE;

CREATE OR REPLACE FUNCTION ddsketch_add(sketch ddsketch, "values" double precision[], alpha double precision = NULL, max_buckets int = NULL)
    RETURNS ddsketch
    AS 'ddsketch', 'ddsketch_add_double_array_increment'
    LANGUAGE C IMMUTABLE PARALLEL SAFE;

CREATE OR REPLACE FUNCTION ddsketch_union(sketch1 ddsketch, sketch2 ddsketch)
    RETURNS ddsketch
    AS 'ddsketch', 'ddsketch_union_double_increment'
    LANGUAGE C IMMUTABLE PARALLEL SAFE;


CREATE OR REPLACE FUNCTION ddsketch_info(sketch ddsketch, out bytes bigint, out flags int, out alpha double precision, out count bigint, out zero_count bigint, out max_buckets int, out negative_buckets int, out positive_buckets int, out min_indexable double precision, out max_indexable double precision)
    RETURNS record
    AS 'ddsketch', 'ddsketch_sketch_info'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION ddsketch_buckets(sketch ddsketch, out index int, out bucket_index int, out bucket_lower double precision, out bucket_upper double precision, out bucket_length double precision, out bucket_count bigint)
    RETURNS SETOF record
    AS 'ddsketch', 'ddsketch_sketch_buckets'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION ddsketch_info(alpha double precision, out min_indexable double precision, out max_indexable double precision)
    RETURNS record
    AS 'ddsketch', 'ddsketch_param_info'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION ddsketch_buckets(alpha double precision, min_value double precision, max_value double precision, out index int, out bucket_index int, out bucket_min double precision, out bucket_max double precision)
    RETURNS SETOF record
    AS 'ddsketch', 'ddsketch_param_buckets'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION ddsketch_percentile(sketch ddsketch, percentile double precision)
    RETURNS double precision
    AS 'ddsketch', 'ddsketch_percentiles'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION ddsketch_percentile(sketch ddsketch, percentiles double precision[])
    RETURNS double precision[]
    AS 'ddsketch', 'ddsketch_array_percentiles'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION ddsketch_percentile_of(sketch ddsketch, "value" double precision)
    RETURNS double precision
    AS 'ddsketch', 'ddsketch_percentiles_of'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION ddsketch_percentile_of(sketch ddsketch, "values" double precision[])
    RETURNS double precision[]
    AS 'ddsketch', 'ddsketch_array_percentiles_of'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION ddsketch_sum(sketch ddsketch, low double precision = 0.0, high double precision = 1.0)
    RETURNS double precision
    AS 'ddsketch', 'ddsketch_sketch_sum'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION ddsketch_avg(sketch ddsketch, low double precision = 0.0, high double precision = 1.0)
    RETURNS double precision
    AS 'ddsketch', 'ddsketch_sketch_avg'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;
