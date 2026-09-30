#!/usr/bin/env bash

mkdir -p corpus/in corpus/recv

for id in $(seq 1 100); do

	rows=$((1 + $RANDOM))

	alpha=$(psql -t -A test -c "select 0.0001 + random() * (0.1 - 0.0001)")
	buckets=$(psql -t -A test -c "select 16 + mod((random() * 1000000)::int, 32000)")

	psql -qAt -z -0 -c "select ddsketch(random(), $alpha, $buckets) from generate_series(1, $rows)" test > tmp
	mv tmp corpus/in/$(sha1sum tmp | awk '{print $1}')

	alpha=$(psql -t -A test -c "select 0.0001 + random() * (0.1 - 0.0001)")
	buckets=$(psql -t -A test -c "select 16 + mod((random() * 1000000)::int, 32000)")

	psql -qAt -c "select encode(ddsketch_send((select ddsketch(random(), $alpha, $buckets) from generate_series(1, $rows))),'base64')" test | base64 -d > tmp
	mv tmp corpus/recv/$(sha1sum tmp | awk '{print $1}')

done
