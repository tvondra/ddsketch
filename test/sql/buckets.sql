-- A target-list SRF and LIMIT also keep the unfixed NaN cases bounded.
SELECT ddsketch_buckets(0.05, 'NaN'::float8, 1.0) LIMIT 1;
SELECT ddsketch_buckets(0.05, -1.0, 'NaN'::float8) LIMIT 1;
SELECT ddsketch_buckets(0.05, '-Infinity'::float8, 1.0) LIMIT 1;
SELECT ddsketch_buckets(0.05, -1.0, 'Infinity'::float8) LIMIT 1;

-- Name the failing endpoint and report its absolute magnitude.
SELECT ddsketch_buckets(0.1, -1.7e308::float8, 1.0) LIMIT 1;
SELECT ddsketch_buckets(0.1, -1.0, 1.7e308::float8) LIMIT 1;
