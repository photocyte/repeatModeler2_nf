##Sourcing a miniconda environment that has nextflow
#source /lab/solexa_weng/testtube/miniconda3/bin/activate
#mkdir results
#bsub -o results/output.log -e results/error.log "nextflow run repeatModeler.nf -resume --genome GCF_000005845.2_ASM584v2_genomic.fna"
nextflow run main.nf -resume --genome examples/U00096.3.fasta
