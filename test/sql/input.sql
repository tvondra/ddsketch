-- Reject overflowing fields instead of saturating or narrowing them.
SELECT 'flags 0 count 9223372036854775808 alpha 0.05 zero_count 0 maxbuckets 16 buckets 1 0 (0, 9223372036854775807)'::ddsketch;
SELECT 'flags 4294967296 count 1 alpha 0.05 zero_count 0 maxbuckets 16 buckets 1 0 (0, 1)'::ddsketch;
SELECT 'flags 0 count 9223372036854775807 alpha 0.05 zero_count 0 maxbuckets 16 buckets 1 0 (0, 9223372036854775808)'::ddsketch;
SELECT 'flags 0 count 1 alpha 0.05 zero_count 0 maxbuckets 16 buckets 1 0 (4294967296, 1)'::ddsketch;
SELECT 'flags 0 count 1 alpha 1e309 zero_count 0 maxbuckets 16 buckets 1 0 (0, 1)'::ddsketch;
SELECT 'flags 0 count 1 alpha 1e-400 zero_count 0 maxbuckets 16 buckets 1 0 (0, 1)'::ddsketch;
SELECT 'flags 0 count 1 alpha NaN zero_count 0 maxbuckets 16 buckets 1 0 (0, 1)'::ddsketch;

-- The sum must not wrap around to a seemingly valid header count.
SELECT 'flags 0 count 1 alpha 0.05 zero_count 0 maxbuckets 16 buckets 3 0 (0, 9223372036854775807) (1, 9223372036854775807) (2, 3)'::ddsketch;

SELECT ddsketch_count('flags 0 count 9223372036854775807 alpha 0.05 zero_count 0 maxbuckets 16 buckets 1 0 (0, 9223372036854775807)'::ddsketch);

-- Valid int32 indexes can still be outside the alpha-dependent mapping range.
SELECT 'flags 0 count 1 alpha 0.1 zero_count 0 maxbuckets 16 buckets 1 0 (2147483647, 1)'::ddsketch;
SELECT 'flags 0 count 1 alpha 0.1 zero_count 0 maxbuckets 16 buckets 1 1 (-2147483648, 1)'::ddsketch;
SELECT 'flags 0 count 1 alpha 0.1 zero_count 0 maxbuckets 16 buckets 1 0 (4000, 1)'::ddsketch;
SELECT 'flags 0 count 1 alpha 0.1 zero_count 0 maxbuckets 16 buckets 1 1 (-4000, 1)'::ddsketch;
SELECT ddsketch_count('flags 0 count 1 alpha 0.05 zero_count 0 maxbuckets 16 buckets 1 0 (4000, 1)'::ddsketch);

-- make sure the output is round-trip safe, regardless of extra_float_digits
WITH sketches AS (
  SELECT ddsketch(v, 0.01234567890123456::float8, 16) AS s FROM (VALUES (-1.0::float8), (0.0), (1.0)) vals(v)
  UNION ALL
  SELECT ddsketch(v, 0.09999999999999999::float8, 16) AS s FROM (VALUES (-1.0::float8), (0.0), (1.0)) vals(v)
  UNION ALL
  SELECT ddsketch(v, 0.05::float8, 16) AS s FROM (VALUES (-1.0::float8), (0.0), (1.0)) vals(v)
)
SELECT ddsketch_union(s, s::text::ddsketch) FROM sketches;

SET extra_float_digits = -3;
WITH sketches AS (
  SELECT ddsketch(v, 0.01234567890123456::float8, 16) AS s FROM (VALUES (-1.0::float8), (0.0), (1.0)) vals(v)
  UNION ALL
  SELECT ddsketch(v, 0.09999999999999999::float8, 16) AS s FROM (VALUES (-1.0::float8), (0.0), (1.0)) vals(v)
  UNION ALL
  SELECT ddsketch(v, 0.05::float8, 16) AS s FROM (VALUES (-1.0::float8), (0.0), (1.0)) vals(v)
)
SELECT ddsketch_union(s, s::text::ddsketch) FROM sketches;

SET extra_float_digits = 0;
WITH sketches AS (
  SELECT ddsketch(v, 0.01234567890123456::float8, 16) AS s FROM (VALUES (-1.0::float8), (0.0), (1.0)) vals(v)
  UNION ALL
  SELECT ddsketch(v, 0.09999999999999999::float8, 16) AS s FROM (VALUES (-1.0::float8), (0.0), (1.0)) vals(v)
  UNION ALL
  SELECT ddsketch(v, 0.05::float8, 16) AS s FROM (VALUES (-1.0::float8), (0.0), (1.0)) vals(v)
)
SELECT ddsketch_union(s, s::text::ddsketch) FROM sketches;

RESET extra_float_digits;
