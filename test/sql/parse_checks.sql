-- test parsing the text representation of sketches (ddsketch_in)
--
-- The input function has to reject malformed values - missing or unexpected
-- parts, fields that fail to parse or are out of range for their data type,
-- a list of buckets that does not match the header, and bucket counts that
-- would overflow. And it has to accept valid values, including values at the
-- boundaries of the allowed ranges.
--

-- don't print the (long) input values for each error
\set VERBOSITY terse

-- a valid sketch (with negative, zero and positive values)
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;

-- missing or misspelled keywords
SELECT ''::ddsketch;
SELECT 'flag 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 cnt 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 max_buckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 bucket 2 1 (0, 1) (7, 1)'::ddsketch;

-- missing and unexpected spaces
SELECT ' flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0  count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1(0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1)(7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1)  (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0 , 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1 ) (7, 1)'::ddsketch;

-- malformed buckets
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 [0, 1] (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0; 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1] (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1'::ddsketch;

-- truncated values
SELECT 'flags 0 count 3'::ddsketch;
SELECT 'flags 0 count 3 alpha'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7,'::ddsketch;

-- trailing characters
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1) '::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)x'::ddsketch;

-- fields that fail to parse
SELECT 'flags x count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count x alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha x zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count x maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets x buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets x 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 x (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (x, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, x) (7, 1)'::ddsketch;

-- fields out of range for their data type
SELECT 'flags 2147483648 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 9223372036854775808 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count -9223372036854775809 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 1e999 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 1e-999 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 9223372036854775808 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 2147483648 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2147483648 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 -2147483649 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (2147483648, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 9223372036854775808) (7, 1)'::ddsketch;

-- negative flags and total count
SELECT 'flags -1 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count -1 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;

-- bucket indexes just outside the indexable range for the alpha
SELECT 'flags 0 count 2 alpha 0.1 zero_count 0 maxbuckets 1024 buckets 2 0 (-3530, 1) (3537, 1)'::ddsketch;
SELECT 'flags 0 count 2 alpha 0.1 zero_count 0 maxbuckets 1024 buckets 2 0 (-3529, 1) (3538, 1)'::ddsketch;

-- a subnormal alpha (positive or negative) is parsed fine, but it's out of
-- the allowed range
SELECT 'flags 0 count 3 alpha 1e-310 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha -1e-310 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;

-- fewer or more buckets than declared in the header
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 3 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 1 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 2 maxbuckets 1024 buckets 0 0 (7, 1)'::ddsketch;

-- the number of negative buckets determines which buckets are negative, and
-- the indexes have to be descending in the negative part and ascending in
-- the positive part
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 0 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 2 (7, 1) (0, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 0 (7, 1) (0, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 2 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 0 (0, 1) (0, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 2 (0, 1) (0, 1)'::ddsketch;

-- total count has to match the zero bucket and the buckets
SELECT 'flags 0 count 2 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 4 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 2 1 (0, -1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 2 maxbuckets 1024 buckets 0 0'::ddsketch;

-- the sum of bucket counts (and the zero bucket) must not overflow, even if
-- it'd wrap around to match the total count
SELECT 'flags 0 count 9223372036854775807 alpha 0.050000 zero_count 1 maxbuckets 1024 buckets 1 0 (0, 9223372036854775807)'::ddsketch;
SELECT 'flags 0 count 9223372036854775807 alpha 0.050000 zero_count 0 maxbuckets 1024 buckets 2 1 (0, 9223372036854775807) (7, 1)'::ddsketch;
SELECT 'flags 0 count 1 alpha 0.050000 zero_count 0 maxbuckets 1024 buckets 3 0 (0, 9223372036854775807) (7, 9223372036854775807) (14, 3)'::ddsketch;

-- valid values at the boundaries of the allowed ranges
SELECT 'flags 0 count 3 alpha 0.0001 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.1 zero_count 1 maxbuckets 1024 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 1 maxbuckets 32768 buckets 2 1 (0, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 16 alpha 0.050000 zero_count 0 maxbuckets 16 buckets 16 8 (7, 1) (6, 1) (5, 1) (4, 1) (3, 1) (2, 1) (1, 1) (0, 1) (0, 1) (1, 1) (2, 1) (3, 1) (4, 1) (5, 1) (6, 1) (7, 1)'::ddsketch;
SELECT 'flags 0 count 3 alpha 0.050000 zero_count 3 maxbuckets 1024 buckets 0 0'::ddsketch;
SELECT 'flags 0 count 2 alpha 0.1 zero_count 0 maxbuckets 1024 buckets 2 0 (-3529, 1) (3537, 1)'::ddsketch;
SELECT 'flags 0 count 2 alpha 0.1 zero_count 0 maxbuckets 1024 buckets 2 2 (3537, 1) (-3529, 1)'::ddsketch;
SELECT 'flags 0 count 9223372036854775807 alpha 0.050000 zero_count 0 maxbuckets 1024 buckets 1 0 (0, 9223372036854775807)'::ddsketch;
SELECT 'flags 0 count 9223372036854775807 alpha 0.050000 zero_count 9223372036854775806 maxbuckets 1024 buckets 1 1 (0, 1)'::ddsketch;
SELECT 'flags 0 count 9223372036854775807 alpha 0.050000 zero_count 0 maxbuckets 1024 buckets 2 1 (0, 9223372036854775806) (7, 1)'::ddsketch;
