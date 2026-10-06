-- Tests ddsketch_add_double_count is correctly rejecting invalid counts
-- transition function for the ddsketch(value, count, alpha, buckets) aggregate.
--
-- The count has to be positive. Make sure zero and negative counts are rejected,
-- both in the first row (when the aggregate state gets created) and in subsequent
-- rows, for values mapped to the positive, negative and zero buckets.

-- zero count
SELECT ddsketch(v, c, 0.05, 1024) FROM (VALUES (1.0::double precision, 0::bigint)) AS t(v, c);
SELECT ddsketch(v, c, 0.05, 1024) FROM (VALUES (-1.0::double precision, 0::bigint)) AS t(v, c);
SELECT ddsketch(v, c, 0.05, 1024) FROM (VALUES (0.0::double precision, 0::bigint)) AS t(v, c);

-- negative count
SELECT ddsketch(v, c, 0.05, 1024) FROM (VALUES (1.0::double precision, -1::bigint)) AS t(v, c);
SELECT ddsketch(v, c, 0.05, 1024) FROM (VALUES (-1.0::double precision, -1::bigint)) AS t(v, c);
SELECT ddsketch(v, c, 0.05, 1024) FROM (VALUES (0.0::double precision, -1::bigint)) AS t(v, c);

-- the most negative bigint value
SELECT ddsketch(v, c, 0.05, 1024) FROM (VALUES (1.0::double precision, '-9223372036854775808'::bigint)) AS t(v, c);

-- invalid count after some valid rows
SELECT ddsketch(v, c, 0.05, 1024) FROM (VALUES (1.0::double precision, 10::bigint), (-1.0, 10), (0.0, 10), (1.0, 0)) AS t(v, c);
SELECT ddsketch(v, c, 0.05, 1024) FROM (VALUES (1.0::double precision, 10::bigint), (-1.0, 10), (0.0, 10), (-1.0, -1)) AS t(v, c);
SELECT ddsketch(v, c, 0.05, 1024) FROM (VALUES (1.0::double precision, 10::bigint), (-1.0, 10), (0.0, 10), (0.0, '-9223372036854775808')) AS t(v, c);

-- negative count must be rejected even if it'd just "undo" a preceding row
SELECT ddsketch(v, c, 0.05, 1024) FROM (VALUES (1.0::double precision, 10::bigint), (1.0, -10)) AS t(v, c);

-- positive counts are accepted, including the largest bigint value
SELECT ddsketch(v, c, 0.05, 1024) FROM (VALUES (1.0::double precision, 1::bigint), (-1.0, 2), (0.0, 3)) AS t(v, c);
SELECT ddsketch(v, c, 0.05, 1024) FROM (VALUES (1.0::double precision, 9223372036854775807::bigint)) AS t(v, c);
