#!/bin/sh

input_query="$1"
out_file="$2"


set -x
curl -H 'Accept: text/tab-separated-values' --data-urlencode 'query@'"$input_query" \
    https://sparql.rhea-db.org/sparql > "$out_file"

