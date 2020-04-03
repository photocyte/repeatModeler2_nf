Channel.fromPath(params.genome).set{ genome_fasta_ch1 }

process RepeatModeler_BuildDatabase {
  cache 'deep'
  publishDir "results/db_dir"
  input:
     file fasta from genome_fasta_ch1
  output:
     file "${fasta}" into cached_genome
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
  ##Put the path and/or version into the stdout
  which BuildDatabase
  BuildDatabase -h
  ##
  BuildDatabase -name \$THENAME -engine ncbi $fasta
  sleep 10 ##Filesystem latency issues
  """
}

cached_genome.into{ genome_fasta_ch2 ; genome_fasta_ch3 ; genome_fasta_ch4 }

process RepeatModeler_execute {
  publishDir "results/rm_out",mode:"copy"
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
     //file "RM_*" //I don't think we need these temporary files.
     file "unaligned.fa" optional true
     file "*-families.fa" into repeat_library_ch
     file "*-families.stk" into repeat_msa_ch    
  script:
  """
  ##From database  
  THENAME=\$(basename ${db_translate})
  THENAME=\${THENAME%.translation}
  ##Put the path and/or version into the stdout
  which RepeatModeler
  RepeatModeler -v
  ##
  RepeatModeler -engine ncbi -pa ${task.cpus} -database \$THENAME
  sleep 10 ##Filesystem latency issues
  """
}

process splitLibraryFasta {
conda "ucsc-fasplit"
input:
 file(inputFasta) from repeat_library_ch
output:
 file("split/*.fa") into library_fasta_split
script:
"""
mkdir split
faSplit about ${inputFasta} 100000 split/
sleep 10 ##Filesystem latency issues
"""
}

genome_fasta_ch2.combine(library_fasta_split.flatten()).set { repeat_masker_tuples }

process RepeatMasker_parallel_execute {
cpus 4
input:
 set file(genome_chunk), file(repeat_library) from repeat_masker_tuples
output:
 file "*.out" into rm_out_chunk
script:
"""
  ##Put the path and/or version into the stdout
  which RepeatMasker
  RepeatMasker -v
  ##
 RepeatMasker -pa ${task.cpus} -gff -qq -lib ${repeat_library} ${genome_chunk}
 sleep 10 ##Filesystem latency issues
"""
}

process convert_out_to_gff3 {
publishDir "results", mode:"copy"
conda "/lab/solexa_weng/testtube/miniconda3/envs/repeatmasker/"
input:
 file rm_out from rm_out_chunk.collectFile(keepHeader:true,skip:3,name:"combined_repeat_masker.outs")
 file genome from genome_fasta_ch3
output:
 file "*.gff3" into repeats_gff_ch
script:
"""
/lab/solexa_weng/testtube/miniconda3/envs/repeatmasker/share/RepeatMasker/util/rmOutToGFF3.pl ${rm_out} > tmp.gff
cat tmp.gff | grep -vP "^#" | gt gff3 -tidy -sort -retainids | uniq > ${genome}.repeats.gff3
sleep 10 ##Filesystem latency issues
"""
}

process soft_mask {
publishDir "results", mode:"copy"
input:
 file repeats_gff from repeats_gff_ch
 file genome from genome_fasta_ch4
output:
  file "softmasked.${genome}"
script:
"""
bedtools maskfasta -soft -fi ${genome} -bed ${repeats_gff} -fo softmasked.${genome}
sleep 10 ##Filesystem latency issues
"""
}

process rename_stockholm_record_ids {
input:
 file msaFile from repeat_msa_ch 
output:
 file "renamed.${msaFile}" into renamed_stockholm
script:
"""
#!/usr/bin/env python
import re
import os
wf = open('renamed.${msaFile}','w')
with open('${msaFile}','r') as rf:
    i=0
    for l in rf.readlines():
        m = re.search('.+:[0-9]+-[0-9]+.+',l)
        if m == None:
             wf.write(l)
        else:
             prefix = str(i)+"_"
             wf.write(prefix+m.group(0)+os.linesep)
        i+=1
wf.close()
"""
}

process convert_stockholm_to_fasta {
publishDir "results", mode:"copy"
input:
 file msaFile from renamed_stockholm
output:
 file "${msaFile}.msa.fa"
conda "hmmer"
script:
"""    
esl-reformat --informat stockholm -o ${msaFile}.msa.fa fasta ${msaFile}
sleep 10 ##Filesystem latency issues
"""
}
