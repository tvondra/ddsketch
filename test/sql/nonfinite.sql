-- Non-finite observations are invalid through every insertion interface.
SELECT ddsketch('NaN'::float8, 0.05, 16);
SELECT ddsketch('Infinity'::float8, 2::bigint, 0.05, 16);
SELECT ddsketch_add(NULL::ddsketch, '-Infinity'::float8, 0.05, 16);
SELECT ddsketch_add(NULL::ddsketch, 'NaN'::float8, 2::bigint, 0.05, 16);
SELECT ddsketch_add(NULL::ddsketch, ARRAY[1.0, 'NaN'::float8], 0.05, 16);

-- NaN must not bypass comparisons with the valid parameter ranges.
SELECT ddsketch(1.0, 'NaN'::float8, 16);
SELECT ddsketch_percentile(ddsketch(1.0, 0.05, 16), ARRAY[0.5, 'NaN'::float8]);
SELECT ddsketch_sum(ddsketch(1.0, 0.05, 16), 'NaN'::float8, 1.0);
SELECT ddsketch_avg(ddsketch(1.0, 0.05, 16), 0.0, 'NaN'::float8);
SELECT ddsketch_sum(ddsketch(1.0, 0.05, 16), '-Infinity'::float8, 1.0);
SELECT ddsketch_avg(ddsketch(1.0, 0.05, 16), 0.0, 'Infinity'::float8);
