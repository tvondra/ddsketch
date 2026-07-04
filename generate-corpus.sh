#!/usr/bin/env bash

mkdir -p corpus/in corpus/recv

for id in $(seq 1 100); do

	buckets=$((16 + RANDOM % (32768 - 16)))

	rows=$((1 + $buckets * 10))

	psql -qAt -z -0 -c "select ddsketch(random(), greatest(0.0001, 0.1 * random()), $buckets) from generate_series(1, $rows)" test > corpus/in/$id

	psql -qAt -c "select encode(ddsketch_send((select ddsketch(random(), greatest(0.0001, 0.1 * random()), $buckets) from generate_series(1, $rows))),'base64')" test | base64 -d > corpus/recv/$id

done
