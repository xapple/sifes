# `sifes` version 2.0.1

The acronym `sifes` stands for **S**pecies **I**​dentification **F**​rom **E**​nvironmental **S**​equencing.

This project is a python pipeline handling the analysis -- from start to finish -- of microbial amplicon sequencing data (e.g. 16S rRNA or other regions).

## Introduction

At first, the main focus when developing `sifes` was to test the functioning of the new protocol we developed in our lab when switching from 454 to Illumina sequencers and to check the coherence and validity of the results obtained. Thus, the pipeline built fits our current needs and is designed to be easily used by the bioinformaticians in our company to quickly analyze the amplicon data that lots of researchers are generating.

The previous version of this pipeline was published under the name `illumitag` here:

[Microbial Community Composition and Diversity via 16S rRNA Gene Amplicons: Evaluating the Illumina Platform](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0116955)

Hence, the `sifes` project is *not* a biologist-oriented tool that supports all the possible use cases one could have with amplicon sequence reads out of the box. For instance, it does not have a graphical interface to operate, nor any bash/sh/csh commands. Indeed, as each sequencing experiment will have different goals and different scientific questions associated to it, there cannot be a standard set of procedures to apply to every dataset. To illustrate this, one could asks ourselves what should the following command do ?

    $ sifes --forward reads_fwd.fasta --reverse reads_rev.fasta

Hard to say. To solve all the underlying questions, the scientist would have to specify an endless list of options and the design of a tool supporting so many different cases would be greatly complicated.

    $ sifes --forward reads_fwd.fasta --reverse reads_rev.fasta --barcode_single TRUE --barcode_only_in_reverse_reads TRUE --use_illumina_i5 FALSE --discard_missmatch_barcode 2 --remove_sequences_from "Plastid, Mitochondrion, Thaumarchaeota" --seperate_phyla_in_graph_when_larger_than 3000 --version_of_silva_to_use SSURef111 etc...

Instead, the `sifes` project *is* a flexible and modular collections of packages written in proper, clean and commented object-oriented python which enables the user to survey, modify and extend the code-base easily -- provided he has a sufficient knowledge in programming. It is a basis upon which the scientist can set up the processing and analysis that he sees fit for his own data sparing him from having to develop lots of the infrastructure needed himself.

Many objects common to any analysis are provided such as a "FASTQ file pair", a "Sample", a "Collection of Samples", a "Cluster of sequences", a "Collection of OTUs", and so on. In addition you will find routines for sending these objects through well-known algorithms such as UCLUST, UPARSE, PandaSEQ, CREST classifier, Vegan NMDS, and so on. Lots of extra functionality is also present such as a multitude of visualizations in `matplotlib` and other things such as the ability to automatically distribute the computation on a network of computers (via SLURM). But here again, every cluster varies between each university and it would make no sense to provide all possible options in the list of command line arguments. Once again, this is why `sifes` is not a command-line tool.

## Installing

No automated installation has been developed for the `sifes` package yet. In the meantime, following this document and typing these commands on your bash prompt should get you started. It is designed so you don't need super user privileges at any step. If you cannot get a functional installation set up, contact the authors.

#### Step 1: Cloning the repository
Here you will download a copy of the code from github and place it somewhere in your home directory.

    $ cd ~
    $ mkdir repos
    $ cd repos
    $ git clone https://github.com/xapple/sifes.git

NB: The access to this repository is not public.

#### Step 2: Modify your search paths
Here you will edit either your ``~/.bashrc`` or ``~/.bash_profile`` to add a reference to the code you just downloaded.

    $ vim ~/.bash_profile
    export PYTHONPATH="$HOME/repos/sifes/":$PYTHONPATH

#### Step 3: Check installation

Check that you have a working installation of Python 3.

Check that you have a working installation of R. Either the system version, or a user version.

#### Step 4: Check you have all the required R dependencies
`sifes` will use some R packages that need to be installed. If you do not have them already, please install them:

    $ R install 'vegan'

#### Step 5: Check you have all the required executables
`sifes` will search for several different binaries as it processes your data. Please check all of these are available in your `$PATH`:

    $ which pandaseq27
    $ which usearch7
    $ which usearch6
    $ which fastqc
    $ which blastn
    $ which classify

## Flowchart
Below is drawn the flowchart describing the data processing along all the steps of `sifes`:

![Flowchart](/../master/documentation/flowchart.png?raw=true "Flowchart")
