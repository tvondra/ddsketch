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

-- do some randomized testing (but be careful not to generate sketches
-- with too few buckets)
DO $$
DECLARE
    v_scale_1    INT;
    v_scale_2    INT;
    v_alpha      DOUBLE PRECISION;
    v_buckets_1  INT;
    v_buckets_2  INT;
    v_rows_1     INT;
    v_rows_2     INT;
    v_rec        RECORD;
BEGIN

    FOR v_scale_1 IN 2..5 LOOP

        FOR v_scale_2 IN 2..5 LOOP

            FOR r IN 1..10 LOOP

                TRUNCATE digest_union_test;

                -- the same alpha value for both sketches
                v_alpha := random() * 0.1; -- max alpha = 0.1
                v_alpha := GREATEST(0.01, LEAST(v_alpha, 0.1));

                -- first random ddsketch
                v_buckets_1 := pow(10, v_scale_1) + (random() * pow(10, (v_scale_1 + 1)))::int;
                v_buckets_1 := GREATEST(500, LEAST(v_buckets_1, 32768));

                v_rows_1 := (random() * v_buckets_1 * 10);
                v_rows_1 := GREATEST(100, LEAST(v_rows_1, 100000));

                -- RAISE NOTICE '% % %', v_alpha, v_buckets_1, v_rows_1;

                -- second random ddsketch
                v_buckets_2 := pow(10, v_scale_1) + (random() * pow(10, (v_scale_1 + 1)))::int;
                v_buckets_2 := GREATEST(500, LEAST(v_buckets_2, 32768));

                v_rows_2 := (random() * v_buckets_2 * 10);
                v_rows_2 := GREATEST(100, LEAST(v_rows_2, 100000));

                -- RAISE NOTICE '% % %', v_alpha, v_buckets, v_rows;

                INSERT INTO digest_union_test VALUES
                ((SELECT ddsketch(v * random(), v_alpha, v_buckets_1) FROM generate_series(1, v_rows_1) v),
                (SELECT ddsketch(v * random(), v_alpha, v_buckets_2) FROM generate_series(1, v_rows_2) v));

                ANALYZE digest_union_test;

                -- FOR v_rec IN EXPLAIN (COSTS OFF) SELECT ddsketch_union(d1, d2) FROM digest_union_test LOOP
                --    RAISE NOTICE '%', v_rec."QUERY PLAN";
                -- END LOOP;

                PERFORM ddsketch_union(d1, d2) FROM digest_union_test;

            END LOOP;

        END LOOP;

    END LOOP;

END;
$$ LANGUAGE plpgsql;

DROP TABLE digest_union_test;
