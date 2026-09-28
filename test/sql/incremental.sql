-- creating a sketch still requires both initialization parameters
SELECT ddsketch_add(NULL::ddsketch, 1.0);
SELECT ddsketch_add(NULL::ddsketch, ARRAY[1.0]::float8[]);
SELECT ddsketch_add(NULL::ddsketch, 1.0, 3::bigint);
SELECT ddsketch_add(NULL::ddsketch, 1.0, p_alpha := 0.05);

-- existing sketches supply their own alpha and capacity
WITH sketch AS (SELECT ddsketch(1.0, 0.05, 16) AS s)
SELECT ddsketch_count(ddsketch_add(s, 2.0)) AS scalar_count,
       ddsketch_count(ddsketch_add(s, ARRAY[2.0, 3.0]::float8[])) AS array_count,
       ddsketch_count(ddsketch_add(s, 2.0, 3::bigint)) AS weighted_count
FROM sketch;

WITH sketch AS (SELECT ddsketch(1.0, 0.05, 16) AS s)
SELECT ddsketch_count(ddsketch_add(p_sketch := s, p_element := 2.0)) AS named_scalar,
       ddsketch_count(ddsketch_add(p_sketch := s, p_elements := ARRAY[2.0]::float8[])) AS named_array,
       ddsketch_count(ddsketch_add(p_sketch := s, p_element := 2.0, p_count := 3::bigint)) AS named_weighted
FROM sketch;

-- pre-existing three-argument alpha form must keep its meaning
SELECT ddsketch_count(ddsketch_add(ddsketch(1.0, 0.05, 16), 2.0, 0.05)) AS explicit_alpha;

-- NULL observations remain no-ops without initialization parameters
SELECT ddsketch_add(NULL::ddsketch, NULL::float8) IS NULL AS null_scalar,
       ddsketch_add(NULL::ddsketch, NULL::float8[]) IS NULL AS null_array,
       ddsketch_add(NULL::ddsketch, NULL::float8, 3::bigint) IS NULL AS null_weighted;
