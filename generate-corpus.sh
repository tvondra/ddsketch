#!/usr/bin/env bash

mkdir corpus

for id in $(seq 1 100); do

	buckets=$((16 + RANDOM % (32768 - 16)))

	rows=$((1 + $buckets * 10))

	echo $buckets $rows

	psql -qAt -c "select encode(ddsketch_send((select ddsketch(random(), greatest(0.0001, 0.1 * random()), $buckets) from generate_series(1, $rows))),'base64')" test | base64 -d > corpus/$id

done
