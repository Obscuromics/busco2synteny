Synteny plotting with BUSCO genes
=========

## Dependencies (via [mamba](https://mamba.readthedocs.io/en/latest/))
-------
```
mamba create -n busco2synteny && \
mamba activate busco2synteny && \
mamba install -c conda-forge -c bioconda ncbi-datasets-cli docopt scipy numpy pandas matplotlib pyarrow
```

# Scripts
-------
Plotting two genomes - '[busco2synteny.py](https://github.com/Obscuromics/busco2synteny/blob/main/scripts/busco2synteny.py)'

Input: 
- Genome file for genome A
- Genome file for genome B
- BUSCO 'full_table.tsv' for genome A
- BUSCO 'full_table.tsv' for genome B

Example command:
'''
python scripts/busco2synteny.py -a genome_A.genomefile.tsv -b genome_B.genomefile.tsv -x genome_A.busco_results.tsv -y genome_B.busco_results.tsv
'''

-------
Plotting any number of genomes - '[busco3synteny.py](https://github.com/Obscuromics/busco2synteny/blob/main/scripts/busco3synteny.py)'

Input: 
- File of paths to genome files
- File of paths to busco files OR liftover file

Example command:
'''
python scripts/busco3synteny.py -a genomefile_paths.txt -x arthropod.buscofile_paths.txt
'''

-------
Plotting any number of genomes from '[A3Cat](https://a3cat.unil.ch/downloads.html)' example script - '[arthropod_synteny.sh]

Input: 
- BUSCO database (must be available for selected genomes)
- Genbank accession codes

Example command:
'''
bash scripts/arthropod_synteny.sh lepidoptera GCA_905220365.1 GCA_905147765.2
'''

-------

BUSCO runs should be against the same reference database.
See examples folder for example input files.
