#!/bin/sh

out_file="$1"


set -x
curl -H 'Accept: text/tab-separated-values' --data-urlencode 'query@rhea_sparql_query' https://sparql.rhea-db.org/sparql > "$out_file"

