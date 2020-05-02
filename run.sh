##Actually, the bioconda install of RepeatModeler is F'ed.  But, tak has it installed by defailt. So maybe it would work?
#export PATH=${PATH}:/lab/solexa_weng/testtube/RECON-1.08/src/
source /lab/solexa_weng/testtube/miniconda3/bin/activate
mkdir results
bsub -o results/output.log -e results/error.log "nextflow run repeatModeler.nf -resume --genome GCF_000005845.2_ASM584v2_genomic.fna"
