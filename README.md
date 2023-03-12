# TranscriptM

A  metatranscriptome bioinformatics pipeline including metagenome contamination correction.



## Installation (conda)
```
git clone git@github.com:sternp/transcriptm.git

cd transcriptm

conda env create -n transcriptm -f transcriptm.yaml

conda activate transcriptm

pip install -e .
```


## Usage (QC, mapping, gDNA decontamination, and read counting pipeline)
```
transcriptm count
    -1 reads_R1.fastq.gz \
    -2 reads_R2.fastq.gz \
    -n 24 \
    --ref combined_reference.fna \
    --gff combined_reference.gff \
    -m 128 \
    -db /dir/to/bowtie2_1 /dir/to/bowtie2_2 \
    -o output_directory

Specifying -g will concatenate a directory of .fna genomes into a single ref sequence and annotate with prokka (time intensive)
Alternatively, you can use pre-contructed files via --ref <combined_reference.fna> and --gff <combined_reference.gff>

Please note, currently the required format for the contig's FASTA headers are <genomeID>_<number>.
For example: >Ardenticatenaceae-ID1234_00001, >Ardenticatenaceae-ID1234_00002...etc

positional arguments:
  {count}

optional arguments:
  -h, --help            show this help message and exit
  --version             Show version information.
  --verbosity VERBOSITY
                        1 = critical, 2 = error, 3 = warning, 4 = info, 5 = debug. Default = 4 (logging) (default: 4)
  --log LOG             Output logging information to file
  -1 <file1>            Forward FASTQ reads (default: none)
  -2 <file2>            Reverse FASTQ reads (default: none)
  --conda-prefix <path>
                        Path to the location of installed conda environments, or where to install new environments (default: /work/microbiome/users/sternesp/conda/envs/transcriptm-dev/envs/)
  -n <num>, --n-cores <num>
                        Maximum number of cores available for use. (default: 8)
  -m <num>, --max-memory <num>
                        Maximum memory for available usage in gigabytes, GB (default: 64)
  -o <dir>, --output <dir>
                        Output directory (default: ./)
  --dry-run [DRYRUN]    Perform snakemake dry run, tests workflow order and conda environments
  --conda-frontend <type>
                        Which conda frontend to use (default: mamba)
  -db <dir> [<dir> ...]
                        Location of one or more Bowtie2-formatted databases for contamination filtering during read QC (i.e. human, rRNA..etc) (default: none)
  --trimmomatic <options>
                        Apply custom trimmomatic values (default: SLIDINGWINDOW:4:20\ MINLEN:50)
  --sequencer-source <type>
                        NexteraPE, TruSeq2, TruSeq3, none (default: TruSeq3)
  --skip-qc             Skip the read QC step
  -g <dir>, --genome-dir <dir>
                        Directory containing FASTA files of each genome (default: none)
  --ref <file>          A single refernce .fna file of contigs/genomes (default: none)
  --gff <file>          GFF annotation of the reference sequence specified by --ref (default: none)
  -x <extension>, --genome-fasta-extension <extension>
                        File extension of genomes in the genome directory (default: fna)
  --kingdom <type>      For use in Prokka when constructing & annotating a new reference from .fna files (default: Bacteria)
  --min-read-aligned-percent <num>
                        Minimum read alignment percent for CoverM filtering (scale from 0-1) (default: 0.75)
  --min-read-percent-identity <num>
                        Minimum read percent identity for CoverM filtering (scale from 0-1) (default: 0.95)
  --gDNA <num>          Median x-fold gDNA coverage to enable gDNA contamination correction. (default: 1)
```

## Notes
You may need to download databases to filter contaminating reads from rRNA genes, the human genome...etc. Otherwise you can make your own Bowtie2-formatted database. You can specify multiple Bowtie2 databases in the transcriptm command, i.e. ```-db /dir/to/db1 /dir/to/db2```

Pre-made databases can be downloaded like so:
```
conda activate transcriptm

kneaddata_database --download human_genome bowtie2 /work/microbiome/db/human_grch37_bowtie2
```
Available pre-made databases:
```
KneadData Databases ( database : build = location )
human_genome : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1.tar.gz
human_genome : bmtagger = http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_BMTagger_v0.1.tar.gz
human_transcriptome : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg38_transcriptome_Bowtie2_v0.1.tar.gz
ribosomal_RNA : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/SILVA_128_LSUParc_SSUParc_ribosomal_RNA_v0.2.tar.gz
mouse_C57BL : bowtie2 = http://huttenhower.sph.harvard.edu/kneadData_databases/mouse_C57BL_6NJ_Bowtie2_v0.1.tar.gz

```
