# pfizer-opensource/nf-crispr-primer-design

## Introduction

**pfizer-opensource/nf-crispr-primer-design** is a bioinformatics pipeline that designs primers
for sequencing CRISPR-edited target sites. The pipeline takes CRISPR guide sequences as
input, finds their location in a reference genome, and then designs primers to amplify
the edited region. In addition to outputing target-specific primers, the pipeline can
also generate primers with NGS adapter sequences added on. This pipeline is a companion
to the **pfizer-opensource/nf-gene-editing-ngs** pipeline, and produces an `amplicons.yaml` file
that is compatible with that gene editing NGS analysis pipeline.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

The only required input file is a tab-delimited text file with guide identifiers and
guide sequences. The input data should look like this:

`guides.txt`:

```tsv
guide_id	guide_seq
HPRT-ctrl	AATTATGGGGATTACTAGGA
```

Each row represents a CRISPR guide, and a primer pair will be designed for that guide.

You will also need to tell the pipeline which reference genome to use for the primer
design and a prefix string to use for naming the output files.

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run pfizer-opensource/nf-crispr-primer-design \
   -profile <docker/singularity/.../institute> \
   --input guides.txt \
   --genome GRCh38 \
   --prefix my_guides \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

pfizer-opensource/nf-crispr-primer-design was originally written by Jason Arroyo.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
