# `sifes` version 2.0.0

The acronym `sifes` stands for **S**pieces **I**​dentification **F**​rom **E**​nvironmental **S**​equencing.

This project is a python pipeline handling the analysis -- from start to finish -- of 16S rRNA microbial amplicon sequencing data.

This source code is propriety of Lucas Sinclair <lucas@envonautics.com>, co-founder of Envonautics Ltd. (https://www.envonautics.com)

## Introduction

At first, the main focus when developing `sifes` was to test the functioning of the new protocol we developed in our lab when switching from 454 to Illumina sequencers and to check the coherence and validity of the results obtained. Thus, the pipeline built fits our current needs and is designed to be easily used by the bioinformaticians in our company to quickly analyze the 16S experiments that lots of researchers are generating.

The previous version of this pipeline was published under the name `illumitag` here:

[Microbial Community Composition and Diversity via 16S rRNA Gene Amplicons: Evaluating the Illumina Platform](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0116955)

Hence, the `sifes` project is *not* a biologist-oriented tool that supports all the possible use cases one could have with 16S rRNA sequence reads out of the box. For instance, it does not have a graphical interface to operate, nor any bash/sh/csh commands. Indeed, as each sequencing experiment will have different goals and different scientific questions associated to it, there cannot be a standard set of procedures to apply to every dataset. To illustrate this, one could asks ourselves what should the following command do ?

    $ sifes --forward reads_fwd.fasta --reverse reads_rev.fasta

Hard to say. To solve all the underlying questions, the scientist would have to specify an endless list of options and the design of a tool supporting so many different cases would be greatly complicated.

    $ sifes --forward reads_fwd.fasta --reverse reads_rev.fasta --barcode_single TRUE --barcode_only_in_reverse_reads TRUE --discard_missmatch_barcode 2 --remove_sequences_from "Plastid, Mitochondrion, Thaumarchaeota" --seperate_phyla_in_graph_when_larger_than 3000 --version_of_silva_to_use SSURef111 etc...

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

#### Step 3: Install your own version of python
Your system probably comes with a version of python installed. But the variations from system to system are too great to rely on any available python. We strongly suggest to just install our own version in your home directory.

For this we will be using this excellent project: https://github.com/yyuu/pyenv

To install it you may use this sister project: https://github.com/yyuu/pyenv-installer

Basically you just need to type this command:

    $ curl -L https://raw.githubusercontent.com/yyuu/pyenv-installer/master/bin/pyenv-installer | bash

These lines go into your `.bash_profile`:

    $ vim ~/.bash_profile
    export PYENV_ROOT="$HOME/.pyenv"
    export PATH="$PYENV_ROOT/bin:$PATH"
    eval "$(pyenv init -)"

Relaunch your shell and type these commands to get the right version of python now:

    pyenv install 2.7.11
    pyenv rehash
    pyenv global 2.7.11

#### Step 4: Install all required python packages
`sifes` uses many third party python libraries. You can get them by running these commands:

    $ pip install sh
    $ pip install decorator
    $ pip install biopython
    $ pip install threadpool
    $ pip install patsy
    $ pip install scipy
    $ pip install matplotlib
    $ pip install pandas
    $ pip install statsmodels
    $ pip install ipython
    $ pip install scikit-learn
    $ pip install rpy2
    $ pip install brewer2mpl
    $ pip install regex
    $ pip install ftputil
    $ pip install names
    $ pip install shell_command
    $ pip install pystache
    $ pip install tabulate
    $ pip install tqdm
    $ pip install humanfriendly
    $ pip install biom-format
    $ pip install future
    $ pip install scikit-bio

Don't forget to rehash the binary links at the end:

    $ pyenv rehash

#### Step 5: Check you have all the required executables
`sifes` will search for several different binaries as it processes your data. Please check all of these are available in your `$PATH`:

    $ which pandaseq27
    $ which usearch7
    $ which usearch6
    $ which fastqc
    $ which blastn
    $ which classify

#### Step 6: Check you have all the required R dependencies
`sifes` will use some R packages that need to be installed. If you do not have them already, please install them:

    $ R install 'vegan'

#### Step 7: Make a working directory with the raw data linked
By default, `sifes` will search for the sequence data in a directory called `SIFES` placed in your home directory. This can be modified of course for your own setup. Each specific collection of sequence data should have an associated `json` file placed in the `metadata/json` directory of the repository telling `sifes` exactly what the name of the raw input files are.

#### Step 8: Start typing python commands to analyze your data

    $ ipython -i -c "import sifes"

## Flowchart
Below is drawn the flowchart describing the data processing along all the steps of `sifes`:

![Flowchart](/../master/documentation/flowcharts/pipeline_overview.png?raw=true "Flowchart")
