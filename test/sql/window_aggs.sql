-- test ddsketch aggregates used as window functions
--
-- With a frame starting at the beginning of the partition, the final
-- function is called for each row, and then more values are added to the
-- same aggregate state, so the final function must not modify the state.
-- With a moving frame start, the aggregate has to be restarted for each
-- row. In both cases the result has to be the same as when building the
-- sketch from the rows in the frame directly.

CREATE TABLE window_aggs (id int, v double precision, c bigint, s ddsketch);

INSERT INTO window_aggs (id, v, c) VALUES
  (1, 1.0, 1), (2, -2.0, 2), (3, 0.0, 1), (4, 3.0, 3), (5, 1.0, 1),
  (6, -5.0, 2), (7, 100.0, 1), (8, 0.5, 2), (9, -0.5, 1), (10, 2.0, 3);

UPDATE window_aggs SET s = ddsketch_add(NULL::ddsketch, v, c, 0.05, 1024);

-- frames starting at the beginning of the partition
SELECT id,
       ddsketch_count(ddsketch(v, 0.05, 1024) OVER w) AS count,
       (ddsketch(v, 0.05, 1024) OVER w)::text =
         (SELECT ddsketch(x.v, 0.05, 1024)::text FROM window_aggs x WHERE x.id <= a.id) AS vals,
       (ddsketch(v, c, 0.05, 1024) OVER w)::text =
         (SELECT ddsketch(x.v, x.c, 0.05, 1024)::text FROM window_aggs x WHERE x.id <= a.id) AS counts,
       (ddsketch(s) OVER w)::text =
         (SELECT ddsketch(x.s)::text FROM window_aggs x WHERE x.id <= a.id) AS sketches
  FROM window_aggs a
WINDOW w AS (ORDER BY id)
 ORDER BY id;

SELECT id,
       ddsketch_count(ddsketch(v, 0.05, 1024) OVER w) AS count,
       (ddsketch(v, 0.05, 1024) OVER w)::text =
         (SELECT ddsketch(x.v, 0.05, 1024)::text FROM window_aggs x WHERE x.id % 2 = a.id % 2 AND x.id <= a.id) AS vals,
       (ddsketch(v, c, 0.05, 1024) OVER w)::text =
         (SELECT ddsketch(x.v, x.c, 0.05, 1024)::text FROM window_aggs x WHERE x.id % 2 = a.id % 2 AND x.id <= a.id) AS counts,
       (ddsketch(s) OVER w)::text =
         (SELECT ddsketch(x.s)::text FROM window_aggs x WHERE x.id % 2 = a.id % 2 AND x.id <= a.id) AS sketches
  FROM window_aggs a
WINDOW w AS (PARTITION BY id % 2 ORDER BY id)
 ORDER BY id;

-- moving frames
SELECT id,
       ddsketch_count(ddsketch(v, 0.05, 1024) OVER w) AS count,
       (ddsketch(v, 0.05, 1024) OVER w)::text =
         (SELECT ddsketch(x.v, 0.05, 1024)::text FROM window_aggs x WHERE x.id BETWEEN a.id - 2 AND a.id + 1) AS vals,
       (ddsketch(v, c, 0.05, 1024) OVER w)::text =
         (SELECT ddsketch(x.v, x.c, 0.05, 1024)::text FROM window_aggs x WHERE x.id BETWEEN a.id - 2 AND a.id + 1) AS counts,
       (ddsketch(s) OVER w)::text =
         (SELECT ddsketch(x.s)::text FROM window_aggs x WHERE x.id BETWEEN a.id - 2 AND a.id + 1) AS sketches
  FROM window_aggs a
WINDOW w AS (ORDER BY id ROWS BETWEEN 2 PRECEDING AND 1 FOLLOWING)
 ORDER BY id;

-- percentiles over a moving frame
SELECT id,
       round(ddsketch_percentile(ddsketch(v, 0.01, 1024) OVER w, 0.5)::numeric, 4) AS estimate,
       (SELECT lower_quantile(x.v, 0.5) FROM window_aggs x WHERE x.id BETWEEN a.id - 2 AND a.id + 1) AS exact
  FROM window_aggs a
WINDOW w AS (ORDER BY id ROWS BETWEEN 2 PRECEDING AND 1 FOLLOWING)
 ORDER BY id;

DROP TABLE window_aggs;
