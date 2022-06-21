# docker run -it -v /Users/whoisthatdog/Documents/Consulting/design_pipeline/designLibraries/Homo_sapiens/Refseq:/genome -e BOWTIE2_INDEXES=/genome -e GENOMES_2BIT=/genome uditas:latest "/bin/bash"
docker run -it -v /Users/whoisthatdog/Documents/Consulting/GenerationBio/Analysis/chr12:/genome -e BOWTIE2_INDEXES=/genome -e GENOMES_2BIT=/genome uditas:latest "/bin/bash"
