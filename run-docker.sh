# docker run -it -v /Users/whoisthatdog/Documents/Consulting/design_pipeline/designLibraries/Homo_sapiens/Refseq:/genome -e BOWTIE2_INDEXES=/genome -e GENOMES_2BIT=/genome uditas:latest "/bin/bash"
#docker run -it -v /Users/whoisthatdog/Documents/Consulting/GenerationBio/Analysis/chr12:/genome -e BOWTIE2_INDEXES=/genome -e GENOMES_2BIT=/genome uditas:latest "/bin/bash"
docker run -it -v /uditas/data/30-702792638/uditas:/data -v /uditas/genomes/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index:/genome -e BOWTIE2_INDEXES=/genome -e GENOMES_2BIT=/genome uditas:latest "/bin/bash" 
#conda run -n uditas_env -skip_demultiplexing -skip_trimming -ncpu 8 /data"
