UDiTaS v1.1 Generation Bio
==========================

Overview
--------

This repo is an updated version of the editasmedicine/uditas repo. It had been updated to work with python 3.9+ and docker. 

UDiTaS(TM) stands for UniDirectional Targeted Sequencing, a novel sequencing method useful for measuring small indels as well as
structural rearrangements, like translocations, in a single reaction.

See details of the method in Giannoukos et al. BMC Genomics (2018) 19:212, https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4561-9


Systems Requirements
--------------------

UDiTaS has been tested in python 3.9 and requires the python packages and tools listed in the files uditas_env_step1.yml and uditas_env_step2.yml. These are separated because of the challenges conda has with resolving installation requirements for environments with a large number of packages. 

Build the docker image with

`docker build -t uditas .`

This builds the image and tags it as `uditas:latest`. The code requires env vars, `BOWTIE2_INDEXES` and `GENOMES_2BIT`, that contain the location of genome indexes and 2bit files, as well as mounted access points for those files, and the input data. Assuming the indexes and 2bit files are in the same directory, these can be specified in the `docker run` call as follows

`docker run --rm -it -v /path/to/data:/data -v /path/to/genome/indexes:/genome -e BOWTIE2_INDEXES=/genome -e GENOMES_2BIT=/genome uditas:latest "/bin/bash"`

The conda environment should be activated automatically, but if not activate the conda environment with

`conda activate uditas_env`

Copy the code and data to an accessible directory and test the installation with

`cp -r /inst .`
`cd inst`
`pytest`

(You must have hg38 files present in the directories specified in the env vars above for this to work)

This will process the test data in

```
data/fig2a
data/fig2b
data/fig2c
```

These data are a subsample of the data displayed in Fig 2 of the paper.

Usage
-----
`uditas` is the main command to launch the UDiTaS analysis. The required argument is: `dir_sample`, the path of the directory with the data to be analyzed.

`dir_sample` should contain the fastq.gz files for R1, R2, I1 and I2. Used in the demultiplexing step.

`dir_sample` should also contain the analysis sheet `sample_info.csv` containing the description of all the samples, with their barcodes and guides used. See examples in the folder `data`

Once the setup has been run, the code can be run as

`uditas ./data/fig2a`

The full list of options are:

```
usage: uditas [-h] [-folder_genome_2bit FOLDER_GENOME_2BIT]
              [-skip_demultiplexing SKIP_DEMULTIPLEXING]
              [-skip_trimming SKIP_TRIMMING]
              [-skip_genome_local_alignment SKIP_GENOME_LOCAL_ALIGNMENT]
              [-skip_genome_global_alignment SKIP_GENOME_GLOBAL_ALIGNMENT]
              [-process_amplicon PROCESS_AMPLICON]
              [-skip_amplicon_global_alignment SKIP_AMPLICON_GLOBAL_ALIGNMENT]
              [-check_plasmid_insertions CHECK_PLASMID_INSERTIONS]
              [-skip_plasmid_alignment SKIP_PLASMID_ALIGNMENT] [-ncpu NCPU]
              [-window_size WINDOW_SIZE]
              [-default_amplicon_window_around_cut DEFAULT_AMPLICON_WINDOW_AROUND_CUT]
              [-min_MAPQ MIN_MAPQ] [-min_AS MIN_AS]
              [-process_AMP_seq_run PROCESS_AMP_SEQ_RUN]
              dir_sample

Process UDiTaS data

positional arguments:
  dir_sample            Directory with the sample to be processed

optional arguments:
  -h, --help            show this help message and exit
  -folder_genome_2bit FOLDER_GENOME_2BIT
                        Folder containing the 2bit file(s) with the reference
                        genome being used (default: GENOMES_2BIT)
  -skip_demultiplexing SKIP_DEMULTIPLEXING
                        Skip demultiplexing? Options: 0, 1 (skip) (default: 0)
  -skip_trimming SKIP_TRIMMING
                        Skip adapter trimming? Options: 0, 1 (skip) (default:
                        0)
  -skip_genome_local_alignment SKIP_GENOME_LOCAL_ALIGNMENT
                        Skip genome-wide local alignment? Options: 0 , 1
                        (skip) (default: 1)
  -skip_genome_global_alignment SKIP_GENOME_GLOBAL_ALIGNMENT
                        Skip genome-wide global alignment? Options: 0 , 1
                        (skip) (default: 0)
  -process_amplicon PROCESS_AMPLICON
                        Select row number (0-based) of amplicon to process,
                        set to all to process all amplicons (default: all)
  -skip_amplicon_global_alignment SKIP_AMPLICON_GLOBAL_ALIGNMENT
                        Skip amplicon global alignment? Options: 0, 1 (skip)
                        (default: 0)
  -check_plasmid_insertions CHECK_PLASMID_INSERTIONS
                        Check for plasmid insertions. Options: 0 (skip), 1
                        plamid_name and plasmid_sequence required in
                        sample_info.csv (default: 1)
  -skip_plasmid_alignment SKIP_PLASMID_ALIGNMENT
                        Skip plasmid alignment? Note, just alignment. Counts
                        still evaluated. Options: 0, 1 (skip) (default: 0)
  -ncpu NCPU            Number of CPUs to use (default: 4)
  -window_size WINDOW_SIZE
                        Window size around cut sites used to grab UDiTaS reads
                        (default: 15)
  -default_amplicon_window_around_cut DEFAULT_AMPLICON_WINDOW_AROUND_CUT
                        Window size around cut sites used to create amplicons
                        (default: 1000)
  -min_MAPQ MIN_MAPQ    Minimum mapping quality to include a read (default: 5)
  -min_AS MIN_AS        Minimum alignment score to include a read (default:
                        -180)
  -process_AMP_seq_run PROCESS_AMP_SEQ_RUN
                        Set to 1 to process an AMP-seq run using GUIDE-seq
                        adapters (default: 0)
```

Preprocessing
-------------

If the input data is already demultiplexed and missing UMIs, there is a utilty script that will preprocess the data so that it can be run by `uditas`.

`python /inst/uditas_utils.py preprocess /path/to/input/fastq /path/to/uditas/output`

The sample_info.csv file must be present in the output directory. The script looks for sample fastq files by concatenating the `Sample` column in the sample_info file with `_R[12]_001.fastq.gz`. It will then create the directory structure and copy fq files from the input directory to the output. It will also create umi files with a unique umi per read pair.
