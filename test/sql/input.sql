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
