-- Empty arrays and NULL elements must identify the offending argument.
SELECT ddsketch_add(NULL::ddsketch, ARRAY[]::float8[], 0.05, 16);
SELECT ddsketch_percentile(ddsketch(1.0, 0.05, 16), ARRAY[]::float8[]);
SELECT ddsketch_percentile_of(ddsketch(1.0, 0.05, 16), ARRAY[]::float8[]);
SELECT ddsketch_add(NULL::ddsketch, ARRAY[1.0, NULL]::float8[], 0.05, 16);
SELECT ddsketch_percentile(ddsketch(1.0, 0.05, 16), ARRAY[0.5, NULL]::float8[]);
SELECT ddsketch_percentile_of(ddsketch(1.0, 0.05, 16), ARRAY[1.0, NULL]::float8[]);
