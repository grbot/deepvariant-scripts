# deepvariant-scripts

## Set environment variables in ~/.bashrc
export PATH=$PATH:/cbio/soft/jdk-11.0.2/bin:/cbio/soft/nextflow/
export JAVA_CMD=/cbio/soft/jdk-11.0.2/bin/java/java
export JAVA_HOME=/cbio/soft/jdk-11.0.2
export NXF_OPTS="-Xms50m -Xmx24g"

## Workflow

1) Set sample sheet. See `NA12878.samplesheet.small.tsv` (contains Fastq reads)
2) Set parameters in `nextflow.config`
3) Run: `nextflow run main.nf -w /cbio/projects/020/gerrit/work -profile ilifu -resume`

