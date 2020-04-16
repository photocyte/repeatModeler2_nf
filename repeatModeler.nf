Channel.fromPath(params.genome).set{ genome_fasta_ch1 }

process checksum_input {
executor 'local'
conda 'seqkit openssl coreutils'
publishDir "results",pattern:"input*.checksum.txt",mode:"copy",overwrite:"true"
input:
 file genome from genome_fasta_ch1
output:
 tuple env(FASCHK),file("*-*-*-*__${genome}") into checksummed_genome_ch
 file "input.*.checksum.txt"
tag "${genome}"
shell:
'''
f=!{genome}
FCHK=$(cat $f | openssl md5 | cut -f 2 -d " " | cut -c1-4) ##First 4 characters of a file contents md5 checksum
IDCHK=$(seqkit seq -n -i $f | sort | openssl md5 | cut -f 2 -d " " | cut -c1-6) ##First 6 characters of a md5 checksum of the sorted, concatenated FASTA IDs
SEQCHK=$(seqkit seq -u $f | seqkit sort -s | seqkit seq -s | openssl md5 | cut -f 2 -d " " | cut -c1-6) ##First 6 characters of a md5 checksum of the sorted, concatenated, uppercase, FASTA sequence
ESEQCHK=$(seqkit sort -s $f | seqkit seq -s | openssl md5 | cut -f 2 -d " " | cut -c1-4) ##First 4 characters of a md5 checksum of the sorted, concatenated FASTA sequence
FASCHK="${FCHK}-${IDCHK}-${SEQCHK}-${ESEQCHK}"
ln -s $f ${FASCHK}__${f}
echo "faschk:${FASCHK}__${f}"
echo "faschk:${FASCHK}__${f}" > input.${f}.checksum.txt
'''
}

process RepeatModeler_BuildDatabase {
  cache 'deep'
  publishDir "results/db_dir"
//  conda "repeatmodeler"
  input:
     tuple val(faschk),file(fasta) from checksummed_genome_ch
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
  ##Print the path and/or version into the stdout
  ##conda list > conda-env.txt
  which BuildDatabase
  ##
  BuildDatabase -name \$THENAME -engine ncbi $fasta
  sleep 10 ##Helps with rare filesystem latency issues
  """
}

cached_genome.into{ genome_fasta_ch2 ; genome_fasta_ch3 ; genome_fasta_ch4 ; genome_fasta_ch5 }

process RepeatModeler_execute {
  storeDir "results/RepeatModeler_out"
//  conda "repeatmodeler"
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
  ##Print the path and/or version into the stdout
  ##conda list > conda-env.txt
  which RepeatModeler
  ##
  RepeatModeler -engine ncbi -pa ${task.cpus} -LTRStruct -database \$THENAME 
  sleep 10 ##Helps with rare filesystem latency issues
  """
}

//Prefer to split this in a process, rather than the built-in "splitFasta" operator as
//having it in a process is more explicit for what is happening and where the files
//are stored
process splitLibraryFasta {
conda "ucsc-fasplit"
input:
 file(inputFasta) from repeat_library_ch
output:
 file("split/*.fa") into library_fasta_split
script:
"""
mkdir split
faSplit about ${inputFasta} 20000 split/
sleep 10 ##Helps with rare filesystem latency issues
"""
}

genome_fasta_ch2.combine(library_fasta_split.flatten()).set { repeat_masker_tuples }

process RepeatMasker_parallel_exec {
cpus 4
// conda "repeatmodeler"
input:
 set file(genome), file(repeat_lib_chunk) from repeat_masker_tuples
output:
 file "*.out" into rm_chunk_out
tag "${repeat_lib_chunk}, ${genome.baseName}"
script:
"""
  ##Print the path and/or version into the stdout
  ##conda list > conda-env.txt
  which RepeatMasker
  RepeatMasker -v
  ##
  RepeatMasker -nolow -no_is -norna -pa ${task.cpus} -gff -q -lib ${repeat_lib_chunk} ${genome}
"""
}

process RepeatMasker_simple_exec {
cpus 8
input:
 file genome from genome_fasta_ch5
output: 
 file "*.out" into rm_simple_out
tag "$genome"
script:
"""
RepeatMasker -noint -pa ${task.cpus} -gff -q ${genome}
"""
}

rm_simple_out.mix(rm_chunk_out).set{rm_outs}

process convert_out_to_gff {
// conda "repeatmodeler"
input:
 file rm_out from rm_outs.collectFile(keepHeader:true,skip:3,name:"combined_repeat_masker.outs")
output:
 file "tmp.gff" into repeats_gff_tmp_ch
script:
"""
#conda list > conda-env.txt
rmOutToGFF3.pl ${rm_out} > tmp.gff
"""
}

process tidy_to_gff3 {
storeDir "results"
conda "genometools-genometools"
input:
 file "tmp.gff" from repeats_gff_tmp_ch
 file genome from genome_fasta_ch3
output:
 file "${genome}.repeats.gff3.gz" into repeats_gff_ch
script:
"""
conda list > conda-env.txt
cat tmp.gff | grep -vP "^#" | gt gff3 -tidy -sort -retainids | uniq | gzip > ${genome}.repeats.gff3.gz
"""
}

process soft_mask {
storeDir "results"
conda "bedtools seqkit"
tag "${genome}"
input:
 file repeats_gff from repeats_gff_ch
 file genome from genome_fasta_ch4
output:
  file "softmasked.${genome}"
script:
"""
bedtools maskfasta -soft -fi <(seqkit seq -u ${genome}) -bed ${repeats_gff} -fo softmasked.${genome}
sleep 10 ##Helps with rare filesystem latency issues
"""
}

process rename_stockholm_record_ids {
input:
 file msaFile from repeat_msa_ch 
output:
 file "renamed.${msaFile}" into renamed_stockholm
shell:
'''
#!/usr/bin/env python
import re
import os
wf = open('renamed.!{msaFile}','w')
with open('!{msaFile}','r') as rf:
    i=0
    for l in rf.readlines():
        #g = re.search('^#=GF ID\\s+([^\\s]+)',l) ##For now only the "m" result is used
	#d = re.search('#=GF[\\s]+DE[\\s]+(.+)',l) ##For now, only the "m" result is used
        m = re.search('.+:[0-9]+-[0-9]+.+',l)
        if m == None:
             wf.write(l)
        else:
             prefix = str(i)+"_"
             wf.write(prefix+m.group(0)+os.linesep)
        i+=1
wf.close()
'''
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
sleep 10 ##Helps with rare filesystem latency issues
"""
}
