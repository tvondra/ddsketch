-- table storing sketches with vastly different alpha values and bucket counts
CREATE TABLE digest_union_test(d1 ddsketch, d2 ddsketch);

-- two incompatible sketches
INSERT INTO digest_union_test
  SELECT ddsketch(v * random(), 0.01, 1000),
         ddsketch(v * random(), 0.02, 1000)
  FROM generate_series(1, 1000) v;

-- ddsketch_union(d1, d2)
EXPLAIN (COSTS OFF) SELECT ddsketch_union(d1, d2) FROM digest_union_test;
SELECT ddsketch_union(d1, d2) FROM digest_union_test;

DROP TABLE digest_union_test;
