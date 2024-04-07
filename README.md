# DENV-genomics
In this repository one can find curated dengue virus (DENV, serotypes 1-4) complete (or near complete) genomes and metadata pulled from GenBank (updated periodically) to construct global DENV phylogentic analysis. This repository also stores results of analyses performed using [nextstrain](https://nextstrain.org/) (see `auspice` directory) which can be visualized below:

* [DENV-1 Global](https://nextstrain.org/community/grubaughlab/DENV-genomics/d1)
* [DENV-2 Global](https://nextstrain.org/community/grubaughlab/DENV-genomics/d2)
* [DENV-3 Global](https://nextstrain.org/community/grubaughlab/DENV-genomics/d3)
* [DENV-4 Global](https://nextstrain.org/community/grubaughlab/DENV-genomics/d4)

## Sequences

Inside the directory `genomes` one can find alignments (`.fasta` files) of complete or near-complete DENV genomes. Viruses belonging to each serotype are split into individual files. Sequences are referenced with NCBI GenBank accession numbers (if previously published).

## Metadata

Also in the folder `genomes`, `.csv` files can be found, including metadata retrieved from public databases and from our own projects. In these files one can find the following data about the genomes:

* strain (GenBank accession number)
* Sample_ID (GenBank isolate name)
* Sample_Source (GenBank or other)
* date (collection date)
* Location (country or territory)
* Region (global region)
* Serotype
* Sero_Geno (genotype)
* Lineage (major and minor lineage)

## Other resouces

This repository will also contain new resources as they are developed. Currently it contains a workflow to help submit data to GenBank - `genomes`. 

---

**Grubaugh Lab** | Yale School of Public Health (YSPH) | [https://grubaughlab.com/](https://grubaughlab.com/)
