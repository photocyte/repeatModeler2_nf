nextflow.preview.dsl=2

process checksum_input {
executor 'local'
conda 'seqkit openssl coreutils'
publishDir "results",pattern:"input*.checksum.txt",mode:"copy",overwrite:"true"
input:
 path genome
output:
 tuple env(FASCHK),path("*-*-*-*__${genome}"), emit:checksum_and_genome
 path "input.*.checksum.txt"
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
  input:
     path fasta
  output:
     path "*.translation"
     path "*.n*"
  tag "$fasta"
  shell:
  '''
  ##From database  
  THENAME=$(basename !{fasta})
  THENAME=${THENAME%.fasta}
  THENAME=${THENAME%.fa}
  THENAME=${THENAME%.fna}
  ##Print the path and/or version into the stdout
  ##conda list > conda-env.txt
  which BuildDatabase
  ##
  BuildDatabase -name $THENAME -engine ncbi !{fasta}
  sleep 10 ##Helps with rare filesystem latency issues
  '''
}

process RepeatModeler_modelRepeatLibrary {
  //storeDir "results/RepeatModeler_out"
//  stageInMode 'copy'
//  memory '180 GB'
//  scratch 'ram-disk'
  cpus params.cpuNum
  //queue '18'
  input:
     path db_translate
     path db_blastdb
  tag "$db_translate"
  output:
     //path "RM_*" //I don't think we need these temporary files.
     path "*-families.fa", emit: repeat_library_ch
     path "*-families.stk", emit: repeat_msa_ch    
     path "unaligned.fa" optional true
  shell:
  '''
  ##From database  
  THENAME=$(basename !{db_translate})
  THENAME=${THENAME%.translation}
  ##Print the path and/or version into the stdout
  ##conda list > conda-env.txt
  which RepeatModeler
  ##
  RepeatModeler -engine ncbi -pa !{task.cpus} -LTRStruct -database \$THENAME 
  sleep 10 ##Helps with rare filesystem latency issues
  '''
}

//Prefer to split this in a process, rather than the built-in "splitFasta" operator as
//having it in a process is more explicit for what is happening and where the files
//are stored
process splitLibraryFasta {
conda "ucsc-fasplit"
input:
 path(inputFasta)
output:
 path "split/*.fa", emit: library_fasta_split
shell:
'''
mkdir split
faSplit about !{inputFasta} 20000 split/
sleep 10 ##Helps with rare filesystem latency issues
'''
}


process RepeatMasker_parallel_exec {
cpus params.cpuNum
// conda "repeatmodeler"
input:
 tuple path(genome), path(repeat_lib_chunk)
output:
 path "*.out"
tag "${repeat_lib_chunk}, ${genome.baseName}"
shell:
'''
  ##Print the path and/or version into the stdout
  ##conda list > conda-env.txt
  which RepeatMasker
  RepeatMasker -v
  ##
  RepeatMasker -nolow -no_is -norna -pa !{task.cpus} -gff -q -lib !{repeat_lib_chunk} !{genome}
'''
}

process RepeatMasker_simple_exec {
cpus params.cpuNum
input:
 path genome
output: 
 path "*.out", emit: rm_simple_out
tag "$genome"
script:
"""
RepeatMasker -noint -pa ${task.cpus} -gff -q ${genome}
"""
}

process convert_out_to_gff {
executor 'local'
// conda "repeatmodeler"
input:
 path rm_out
output:
 path "tmp.gff", emit: repeats_gff_tmp_ch
shell:
'''
#conda list > conda-env.txt
rmOutToGFF3.pl !{rm_out} > tmp.gff
'''
}

process tidy_to_gff3 {
executor 'local'
//storeDir "results"
conda "genometools-genometools"
input:
 path genome
 path "tmp.gff"
output:
 path "${genome}.repeats.gff3.gz", emit: repeats_gff_ch
shell:
'''
set -o pipefail
conda list > conda-env.txt
cat tmp.gff | grep -vP "^#" | gt gff3 -tidy -sort -retainids | uniq | gzip > !{genome}.repeats.gff3.gz
'''
}

process soft_mask {
executor 'local'
//storeDir "results"
conda "bedtools seqkit"
tag "${genome}"
input:
 path genome
 path repeats_gff
output:
  path "softmasked.${genome}"
shell:
"""
set -o pipefail
bedtools maskfasta -soft -fi <(seqkit seq -u !{genome}) -bed !{repeats_gff} -fo softmasked.!{genome}
sleep 10 ##Helps with rare filesystem latency issues
"""
}

process rename_stockholm_record_ids {
executor 'local'
input:
 path msaFile
output:
 path "renamed.${msaFile}", emit: renamed_stockholm
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
executor 'local'
publishDir "results", mode:"copy"
input:
 path msaFile
output:
 path "${msaFile}.msa.fa"
conda "hmmer"
shell:
'''   
esl-reformat --informat stockholm -o !{msaFile}.msa.fa fasta !{msaFile}
sleep 10 ##Helps with rare filesystem latency issues
'''
}

workflow modelAndMaskGenome_wf {
 take: genome
 main:
  checksum_input(genome)
  checksum_input.out.checksum_and_genome.map{ it[1] }.set{cached_genome}
  RepeatModeler_BuildDatabase(cached_genome)
  
  RepeatModeler_modelRepeatLibrary(RepeatModeler_BuildDatabase.out)
  rename_stockholm_record_ids(RepeatModeler_modelRepeatLibrary.out.repeat_msa_ch)
  convert_stockholm_to_fasta(rename_stockholm_record_ids.out.renamed_stockholm)

  splitLibraryFasta = RepeatModeler_modelRepeatLibrary.out.repeat_library_ch.splitFasta(size:20000)

  repeat_masker_tuples = cached_genome.combine(splitLibraryFasta)
  RepeatMasker_parallel_exec(repeat_masker_tuples)

  RepeatMasker_simple_exec(cached_genome)
  mergedOuts = RepeatMasker_simple_exec.out.mix(RepeatMasker_parallel_exec.out).collectFile(keepHeader:true,skip:3,name:"combined_repeat_masker.outs")
  convert_out_to_gff(mergedOuts)
  tidy_to_gff3(cached_genome,convert_out_to_gff.out)
  soft_mask(cached_genome,tidy_to_gff3.out)
}

workflow {
 log.info """\
 R E P E A T M O D E L E R 2 - N E X T F L O W
 =============================================
 genome        : ${params.genome}

 If you have set the --genome parameter, and
 have access to Miniconda / Singularity,
 and given the path to your tandem repeat 
 finder (TRF) executable under the 'runOptions' 
 in nextflow.config, this should just work.

 Intermediate results are stored in the
 'work' directory, in the current working dir

 The pipeline will publish results to the
 'results' dir, in the current working dir
 """
 modelAndMaskGenome_wf(Channel.fromPath(params.genome))
}
