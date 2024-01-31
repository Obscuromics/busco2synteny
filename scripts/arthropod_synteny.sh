#!/usr/bin/env bash
# mamba install -c conda-forge ncbi-datasets-cli

GEN_ONE=GCA_947179485.1
GEN_TWO=GCA_030463065.1

#Get BUSCO results from https://a3cat.unil.ch/downloads.html
BUSCO_DB=arthropoda

BUSCO_ONE=${GEN_ONE}_${BUSCO_DB}.tar.gz
BUSCO_TWO=${GEN_TWO}_${BUSCO_DB}.tar.gz

wget https://a3cat.unil.ch/data/busco/${GEN_ONE}/$BUSCO_ONE
wget https://a3cat.unil.ch/data/busco/${GEN_ONE}/$BUSCO_TWO
gunzip $BUSCO_ONE; gunzip $BUSCO_TWO


datasets download genome accession GCA_947179485.1 --assembly-level chromosome --include genome | dataformat tsv genome --fields accession