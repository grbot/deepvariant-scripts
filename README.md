# deepvariant-scripts

## Get the code
```
git clone https://github.com/grbot/deepvariant-scripts.git
```

## Set environment variables in `~/.bashrc`
```
export PATH=$PATH:/cbio/soft/jdk-11.0.2/bin:/cbio/soft/nextflow/
export JAVA_CMD=/cbio/soft/jdk-11.0.2/bin/java/java
export JAVA_HOME=/cbio/soft/jdk-11.0.2
export NXF_OPTS="-Xms50m -Xmx24g"
```

## Workflow

1) Set sample sheet. See `NA12878.samplesheet.small.tsv` (contains Fastq reads)
2) Set parameters in `nextflow.config`. Important set samplesheet and output directory.
3) Run: `nextflow run main.nf -w /cbio/projects/020/gerrit/work -profile ilifu -resume`. Change the working directory to where you have access.

# Note
- A samplesheet with 5 1000 genomes have also been prepared (`1kg.samplesheet.tsv`)
- Example output is in `/cbio/projects/031/gerrit/examples`
