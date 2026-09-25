-----------------------------------------------------------
-- parameter validation
-----------------------------------------------------------

-- invalid percentile value
SELECT ddsketch_percentile(ddsketch(i / 1.0, 0.01, 1024), ARRAY[0.1, -0.1]) FROM generate_series(1,10) s(i);

-- alpha too low
SELECT ddsketch_percentile(ddsketch(i / 1.0, 0.00009, 1024), 0.5) FROM generate_series(1,10) s(i);

-- alpha too high
SELECT ddsketch_percentile(ddsketch(i / 1.0, 0.11, 1024), 0.5) FROM generate_series(1,10) s(i);

-- fewer than minimum number of buckets
SELECT ddsketch_percentile(ddsketch(i / 1.0, 0.01, 15), 0.5) FROM generate_series(1,10) s(i);

-- more than maximum number of buckets
SELECT ddsketch_percentile(ddsketch(i / 1.0, 0.01, 32769), 0.5) FROM generate_series(1,10) s(i);

-- too many buckets needed
SELECT ddsketch_percentile(ddsketch(i / 1.0, 0.01, 32), 0.5) FROM generate_series(1,10000) s(i);

-- invalid parameters
SELECT * FROM ddsketch_buckets(0.0::float8, 1.0, 2.0);
SELECT * FROM ddsketch_buckets('NaN'::float8, 1.0, 2.0);
SELECT * FROM ddsketch_buckets(1e-300::float8, 1.0, 2.0);

SELECT * FROM ddsketch_info(0.0::float8);
SELECT * FROM ddsketch_info(2.0::float8);
SELECT * FROM ddsketch_info(-1.0::float8);

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

-- test info functions
SELECT * FROM ddsketch_info(NULL::double precision);
SELECT * FROM ddsketch_info(NULL::ddsketch);
SELECT * FROM ddsketch_buckets(NULL::double precision, NULL::double precision, NULL::double precision);
SELECT * FROM ddsketch_buckets(NULL::ddsketch);

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


-- test sketches with negative buckets and zero counts

-- 5e-324 is below min_indexable_value, so it goes to the zero bucket.
-- The sketch below has count = 2, zero_count = 2, nbuckets = 0.
SELECT (ddsketch(v, 0.05, 1024))::text::ddsketch AS before
  FROM (VALUES (5e-324::float8), (-5e-324::float8)) t(v);

-- Adding one ordinary value must yield count = 3, zero_count = 2.
-- expected: "count 3 ... zero_count 2 ... buckets 1 0 (0, 1)"
SELECT ddsketch_add((SELECT ddsketch(v, 0.05, 1024)
                       FROM (VALUES (5e-324::float8), (-5e-324::float8)) t(v)),
                    1.0, 0.05, 1024)::text::ddsketch AS after;

-- a sketch with two negative buckets and no positive ones:
--   "... buckets 2 2 (7, 1) (0, 1)"   (nbuckets = 2, nbuckets_negative = 2)
SELECT (ddsketch(v, 0.05, 1024))::text::ddsketch AS before
  FROM (VALUES (-2.0::float8), (-1.0::float8)) t(v);

-- expected: "buckets 3 2 (7, 1) (0, 1) (17, 1)"
SELECT ddsketch_add((SELECT ddsketch(v, 0.05, 1024)
                       FROM (VALUES (-2.0::float8), (-1.0::float8)) t(v)),
                    5.0, 0.05, 1024)::text::ddsketch AS after;

-- same through ddsketch_union()'
SELECT ddsketch_union((SELECT ddsketch(v, 0.05, 1024)
                         FROM (VALUES (-2.0::float8), (-1.0::float8)) t(v)),
                      (SELECT ddsketch(v, 0.05, 1024)
                         FROM (VALUES (5.0::float8)) t(v)))::text::ddsketch;
