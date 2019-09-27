
Channel.fromPath(params.genome).into{ genome_fasta_ch1 ; genome_fasta_ch2 ; genome_fasta_ch3 ; genome_fasta_ch4 }

process RepeatModeler_BuildDatabase {
  publishDir "results/db_dir"
  input:
     file fasta from genome_fasta_ch1
  output:
     file "*.translation" into db_translate_ch
     file "*.n*" into db_blastdb_ch
  tag "$fasta"
  script:
  """
  ##From database  
  THENAME=\$(basename ${fasta})
  THENAME=\${THENAME%.fasta}
  THENAME=\${THENAME%.fa}
  THENAME=\${THENAME%.fna}
  BuildDatabase -name \$THENAME -engine ncbi $fasta
  """
}

process RepeatModeler_execute {
  executor 'lsf'
  publishDir "results/rm_out"
//  stageInMode 'copy'
//  memory '180 GB'
//  scratch 'ram-disk'
  cpus 18
  queue '18'
  input:
     file db_translate from db_translate_ch
     file db_blastdb from db_blastdb_ch
  tag "$db_translate"
  output:
     file "RM_*"
     file "*-families.fa" into repeat_library_ch
     file "*-families.stk"     
  script:
  """
  ##From database  
  THENAME=\$(basename ${db_translate})
  THENAME=\${THENAME%.translation}
  RepeatModeler -engine ncbi -pa ${task.cpus} -database \$THENAME
  """
}

genome_fasta_ch2.splitFasta( by: 10 ).set{ genome_fasta_split }
genome_fasta_split.combine(repeat_library_ch).set{ repeat_masker_tuples }

process RepeatMasker_parallel_execute {
executor 'lsf'
cpus 4
input:
 set file(genome_chunk), file(repeat_library) from repeat_masker_tuples
output:
 file "*.out" into rm_out_chunk
script:
"""
 RepeatMasker -pa ${task.cpus} -gff -qq -lib ${repeat_library} ${genome_chunk}
"""
}

process convert_out_to_gff3 {
publishDir "repeatModeler_nextflow_results", mode:"copy"
conda "/lab/solexa_weng/testtube/miniconda3/envs/repeatmasker/"
input:
 file rm_out from rm_out_chunk.collectFile(keepHeader:true,skip:3,name:"combined_repeat_masker.outs")
 file genome from genome_fasta_ch3
output:
 file "*.gff3" into repeats_gff_ch
script:
"""
/lab/solexa_weng/testtube/miniconda3/envs/repeatmasker/share/RepeatMasker/util/rmOutToGFF3.pl ${rm_out} > tmp.gff
cat tmp.gff | gt gff3 -tidy -sort -retainids > ${genome}.repeats.gff3
"""
}

process soft_mask {
publishDir "repeatModeler_nextflow_results", mode:"copy"
input:
 file repeats_gff from repeats_gff_ch
 file genome from genome_fasta_ch4
output:
  file "softmasked.${genome}"
script:
"""
bedtools maskfasta -soft -fi ${genome} -bed ${repeats_gff} -fo softmasked.${genome}
"""

}
