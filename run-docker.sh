docker run --rm -it -v /uditas/data/30-702792638:/data -v /uditas/genomes/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index:/genome -e BOWTIE2_INDEXES=/genome -e GENOMES_2BIT=/genome uditas:latest
