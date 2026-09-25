CREATE TYPE ddsketch;

CREATE OR REPLACE FUNCTION ddsketch_in(cstring)
    RETURNS ddsketch
    AS 'ddsketch', 'ddsketch_in'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION ddsketch_out(ddsketch)
    RETURNS cstring
    AS 'ddsketch', 'ddsketch_out'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION ddsketch_send(ddsketch)
    RETURNS bytea
    AS 'ddsketch', 'ddsketch_send'
    LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION ddsketch_recv(internal)
    RETURNS ddsketch
    AS 'ddsketch', 'ddsketch_recv'
    LANGUAGE C IMMUTABLE STRICT;

CREATE TYPE ddsketch (
    INPUT = ddsketch_in,
    OUTPUT = ddsketch_out,
    RECEIVE = ddsketch_recv,
    SEND = ddsketch_send,
    INTERNALLENGTH = variable,
    STORAGE = extended,
    ALIGNMENT = 8
);

CREATE OR REPLACE FUNCTION ddsketch_sketch(p_pointer internal)
    RETURNS ddsketch
    AS 'ddsketch', 'ddsketch_sketch'
    LANGUAGE C IMMUTABLE;

CREATE AGGREGATE ddsketch(double precision, double precision, int) (
    SFUNC = ddsketch_add_double,
    STYPE = internal,
    FINALFUNC = ddsketch_sketch,
    SERIALFUNC = ddsketch_serial,
    DESERIALFUNC = ddsketch_deserial,
    COMBINEFUNC = ddsketch_combine,
    PARALLEL = SAFE
);

CREATE OR REPLACE FUNCTION ddsketch_add_sketch(p_pointer internal, p_element ddsketch)
    RETURNS internal
    AS 'ddsketch', 'ddsketch_add_sketch'
    LANGUAGE C IMMUTABLE;

CREATE AGGREGATE ddsketch(ddsketch) (
    SFUNC = ddsketch_add_sketch,
    STYPE = internal,
    FINALFUNC = ddsketch_sketch,
    SERIALFUNC = ddsketch_serial,
    DESERIALFUNC = ddsketch_deserial,
    COMBINEFUNC = ddsketch_combine,
    PARALLEL = SAFE
);

CREATE OR REPLACE FUNCTION ddsketch_count(ddsketch)
    RETURNS bigint
    AS 'ddsketch', 'ddsketch_count'
    LANGUAGE C IMMUTABLE STRICT;


CREATE OR REPLACE FUNCTION ddsketch_add_double_count(p_pointer internal, p_element double precision, p_count bigint, p_alpha double precision, p_buckets int)
    RETURNS internal
    AS 'ddsketch', 'ddsketch_add_double_count'
    LANGUAGE C IMMUTABLE;

CREATE AGGREGATE ddsketch(double precision, bigint, double precision, int) (
    SFUNC = ddsketch_add_double_count,
    STYPE = internal,
    FINALFUNC = ddsketch_sketch,
    SERIALFUNC = ddsketch_serial,
    DESERIALFUNC = ddsketch_deserial,
    COMBINEFUNC = ddsketch_combine,
    PARALLEL = SAFE
);

CREATE OR REPLACE FUNCTION ddsketch_add(p_sketch ddsketch, p_element double precision, p_alpha double precision, p_buckets int = NULL)
    RETURNS ddsketch
    AS 'ddsketch', 'ddsketch_add_double_increment'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION ddsketch_add(p_sketch ddsketch, p_element double precision, p_count bigint, p_alpha double precision, p_buckets int = NULL)
    RETURNS ddsketch
    AS 'ddsketch', 'ddsketch_add_double_count_increment'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION ddsketch_add(p_sketch ddsketch, p_elements double precision[], p_alpha double precision, p_buckets int = NULL)
    RETURNS ddsketch
    AS 'ddsketch', 'ddsketch_add_double_array_increment'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION ddsketch_union(p_sketch1 ddsketch, p_sketch2 ddsketch)
    RETURNS ddsketch
    AS 'ddsketch', 'ddsketch_union_double_increment'
    LANGUAGE C IMMUTABLE;


CREATE OR REPLACE FUNCTION ddsketch_info(p_sketch ddsketch, out bytes bigint, out flags bigint, out alpha double precision, out count bigint, out zero_count bigint, out max_buckets int, out negative_buckets int, out positive_buckets int, out min_indexable double precision, out max_indexable double precision)
    RETURNS record
    AS 'ddsketch', 'ddsketch_sketch_info'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION ddsketch_buckets(p_sketch ddsketch, out index int, out bucket_index int, out bucket_lower double precision, out bucket_upper double precision, out bucket_length double precision, out bucket_count bigint)
    RETURNS SETOF record
    AS 'ddsketch', 'ddsketch_sketch_buckets'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION ddsketch_info(p_alpha double precision, out min_indexable double precision, out max_indexable double precision)
    RETURNS record
    AS 'ddsketch', 'ddsketch_param_info'
    LANGUAGE C IMMUTABLE;

CREATE OR REPLACE FUNCTION ddsketch_buckets(p_alpha double precision, p_min_value double precision, p_max_value double precision, out index int, out bucket_index int, out bucket_min double precision, out bucket_max double precision)
    RETURNS SETOF record
    AS 'ddsketch', 'ddsketch_param_buckets'
    LANGUAGE C IMMUTABLE;

CREATE FUNCTION ddsketch_percentile(ddsketch, double precision)
    RETURNS double precision
    AS 'ddsketch', 'ddsketch_percentiles'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION ddsketch_percentile(ddsketch, double precision[])
    RETURNS double precision[]
    AS 'ddsketch', 'ddsketch_array_percentiles'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION ddsketch_percentile_of(ddsketch, double precision)
    RETURNS double precision
    AS 'ddsketch', 'ddsketch_percentiles_of'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION ddsketch_percentile_of(ddsketch, double precision[])
    RETURNS double precision[]
    AS 'ddsketch', 'ddsketch_array_percentiles_of'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION ddsketch_sum(p_sketch ddsketch, p_low double precision = 0.0, p_high double precision = 1.0)
    RETURNS double precision
    AS 'ddsketch', 'ddsketch_sketch_sum'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION ddsketch_avg(p_sketch ddsketch, p_low double precision = 0.0, p_high double precision = 1.0)
    RETURNS double precision
    AS 'ddsketch', 'ddsketch_sketch_avg'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;
