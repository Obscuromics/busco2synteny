#!/usr/bin/env bash
# mamba install -c conda-forge -c bioconda ncbi-datasets-cli docopt scipy numpy pandas matplotlib pyarrow

#Example command
#bash ../scripts/arthropod_synteny.sh lepidoptera GCA_905220365.1 GCA_905147765.2

get_busco_results () {
  #Get BUSCO results from https://a3cat.unil.ch/downloads.html
  GENB_ACC="$1"
  BUSCO_DB="$2"
  BUSCO_F=${GENB_ACC}_${BUSCO_DB}.tar.gz
  wget -q https://a3cat.unil.ch/data/busco/${GENB_ACC}/$BUSCO_F -P syn_downloads
  echo "Fetching ${GENB_ACC} ${BUSCO_DB} results..."
  tar -xzf syn_downloads/$BUSCO_F -C syn_downloads
  mv syn_downloads/${GENB_ACC}/${BUSCO_DB}/run_${BUSCO_DB}_*/full_table.tsv plot_input_files/${GENB_ACC}.${BUSCO_DB}.busco.tsv
}

get_genome_files () {
  GENB_ACC="$1"
  datasets download genome accession ${GENB_ACC} --include seq-report --filename syn_downloads/${GENB_ACC}/ncbi_dataset.zip
  unzip -qq -o syn_downloads/${GENB_ACC}/ncbi_dataset.zip -d syn_downloads/${GENB_ACC}
  cat "syn_downloads/${GENB_ACC}/ncbi_dataset/data/${GENB_ACC}/sequence_report.jsonl" | \
  dataformat tsv genome-seq | \
  tail -n +2 | \
  awk -F'\t' -v OFS='\t' '{ if ($11 == "assembled-molecule") { print $7,$12,"+",$7} else { exit }}' > plot_input_files/${GENB_ACC}.genomefile.tsv
}

mkdir plot_input_files

SCRIPT_DIR=$(dirname "$0")
BUSCO_DB="$1"

for ACC in ${@: 2}
do
	mkdir -p syn_downloads/${ACC}
	get_busco_results $ACC $BUSCO_DB
	get_genome_files $ACC
	rm -rf syn_downloads
done

ls plot_input_files/*.genomefile.tsv > plot_input_files/genomefile_paths.txt
ls plot_input_files/*.${BUSCO_DB}.busco.tsv > plot_input_files/busco_paths.${BUSCO_DB}.txt

python ${SCRIPT_DIR}/busco3synteny.py -a plot_input_files/genomefile_paths.txt -x plot_input_files/busco_paths.${BUSCO_DB}.txt

