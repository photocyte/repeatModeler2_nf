# repeatModeler2_nf
A Nextflow wrapper for RepeatModeler2.

Dependencies:
* Singularity , (typically has to be installed by an admin / system administrator)
* Miniconda , (can be installed by a user, see https://docs.conda.io/en/latest/miniconda.html)
* Nextflow , (once you have conda installed: `conda install -y nextflow`)

If you have set the `--genome` parameter, and have access to [Miniconda](https://docs.conda.io/en/latest/miniconda.html) & [Singularity](https://sylabs.io/singularity/), this pipeline should 'just work' and automatically pull the necessary RepeatModeler developer provided [Docker containers](https://hub.docker.com/r/dfam/tetools), and install other software via conda.

See here for more background on the dfam/tetools container: [https://github.com/Dfam-consortium/TETools](https://github.com/Dfam-consortium/TETools)

### Running the pipeline
```
git clone https://github.com/photocyte/repeatModeler2_nf.git
#After manually git cloning' the nextflow main.nf & nextflow.config into your working directory
nextflow run main.nf -resume -profile singularity --cpuNum 10 --genome examples/U00096.3.fasta
##-profile docker is also possible
```
or

```
#Letting nextflow manage the git cloning
nextflow pull https://github.com/photocyte/repeatModeler2_nf
nextflow run repeatModeler2_nf -latest -resume -profile singularity --cpuNum 10 --genome examples/U00096.3.fasta 
##-profile docker is also possible
```

### Local vs cluster/HPC execution
By default, the workflow runs locally. If you'd instead prefer Nextflow to submit the processes/jobs to a high-performance-computing (HPC) cluster, edit the first line of `nextflow.config`, to specify your desired cluster type. See https://www.nextflow.io/docs/latest/executor.html for more details. I.e., for running on a PBS cluster change:
`params.executor = "local"` to `params.executor = "pbs"`

### Results
 Look in `./results` once the pipeline is complete

### RepeatModeler2 citation
Flynn JM, Hubley R, Goubert C, Rosen J, Clark AG, Feschotte C, Smit AF. 2020. RepeatModeler2 for automated genomic discovery of transposable element families. PNAS. doi:10.1073/pnas.1921046117

### Directed acyclic graph of pipeline execution
(Note, DAG rendering is a little broken currently)

![Directed acyclic graph for program execution](./results/dag.svg)

### See also
https://github.com/darcyabjones/pante

