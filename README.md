# repeatModeler2_nf
A Nextflow wrapper for RepeatModeler2.

Dependencies:
* Singularity
* Nextflow
* Miniconda

If you have set the '--genome' parameter, and have access to [Miniconda](https://docs.conda.io/en/latest/miniconda.html) & [Singularity](https://sylabs.io/singularity/), this pipeline should 'just work' and automatically pull the necessary [Docker containers](https://hub.docker.com/r/dfam/tetools) and install software via conda.

See here for more background on the dfam/tetools container: [https://github.com/Dfam-consortium/TETools](https://github.com/Dfam-consortium/TETools)

### Running the pipeline
`nextflow run main.nf -resume --genome example.fa`

or

```
nextflow pull https://github.com/photocyte/repeatModeler2_nf
nextflow run repeatModeler2_nf -resume --genome example.fa
```

### Results
 Look in `./results` once the pipeline is complete

### RepeatModeler2 citation
Flynn JM, Hubley R, Goubert C, Rosen J, Clark AG, Feschotte C, Smit AF. 2020. RepeatModeler2 for automated genomic discovery of transposable element families. PNAS. doi:10.1073/pnas.1921046117

### Directed acyclic graph of pipeline execution
(Note, DAG rendering is a little broken currently)

![Directed acyclic graph for program execution](./results/dag.svg)

### See also
https://github.com/darcyabjones/pante

