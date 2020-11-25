# Download GRCh37.75 transcriptome (from https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/04_quasi_alignment_salmon.html)

dir.create('./transcriptome')
system(paste('cd',paste0('./transcriptome'),'&& { wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz ; }'))

# Unzip
system(paste('gunzip ./transcriptome/Homo_sapiens.GRCh37.75.cdna.all.fa.gz'))
