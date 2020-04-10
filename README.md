# RHEA reactions extraction

Dependecies:
* R (> 3.5)
* tidyverse (> 1.2.1)
* curl (> 7.xx)


Clone the git repository
```bash
mkdir /path/to/repo
cd /path/to/repo
git clone https://github.com/metaganal/rhea_rmd.git .
```

Data are already included in the data directory, but you can update to the latest RHEA DB through the following procedure.


Obtain RHEA DB reactions from RHEA SPARQL end-point
```bash
curl -H 'Accept: text/tab-separated-values' --data-urlencode 'query@rhea_sparql_query' https://sparql.rhea-db.org/sparql > data/rhea_db.tsv
```

Wrangle RHEA DB with R to get a table that's easier to process.
Input: RHEA DB
Output: RHEA reactions, RHEA compound usage
```bash
Rscript R/rhea_wrangling.R data/rhea_db.tsv
```

Auto curation of RHEA based on frequency of compound usage.
Compounds with high frequency usage are removed as generic cofactors.
Some compounds are used to define the directions of enzymatic reaction.
```bash
Rscript R/rhea_table_curation.R data/rhea_db_reactions.tsv data/rhea_db_parsed.tsv
```