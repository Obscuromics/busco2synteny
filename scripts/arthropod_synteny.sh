#!/usr/bin/env bash
# mamba install -c conda-forge ncbi-datasets-cli

GEN_ONE=GCA_947179485.1
GEN_TWO=GCA_030463065.1

BUSCO_DB=arthropoda

#Get BUSCO results from https://a3cat.unil.ch/downloads.html
#BUSCO_ONE=${GEN_ONE}_${BUSCO_DB}.tar.gz
#BUSCO_TWO=${GEN_TWO}_${BUSCO_DB}.tar.gz

#wget https://a3cat.unil.ch/data/busco/${GEN_ONE}/$BUSCO_ONE
#wget https://a3cat.unil.ch/data/busco/${GEN_ONE}/$BUSCO_TWO
#gunzip $BUSCO_ONE; gunzip $BUSCO_TWO

#Get genomefiles

N_CHR=${datasets summary genome accession ${GEN_ONE} --as-json-lines | dataformat tsv genome --fields assmstats-total-number-of-chromosomes}
datasets download genome accession ${GEN_ONE} --include seq-report
unzip ncbi_dataset.zip
cat ncbi_dataset/data/${GEN_ONE}/sequence_report.jsonl | dataformat tsv genome-seq | head -n "${N_CHR}" | tail -n "${N_CHR}"

