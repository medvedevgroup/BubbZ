BubbZ 1.0.0
===============

Release date: 21th January 2019
=================

Authors
=======
* Ilia Minkin (Pennsylvania State University)
* Paul Medvedev (Pennsylvania State University)

Introduction
============
BubbZ is a whole-genome homology mapping pipeline. BubbZ works best in the case
when the user needs to find all (possibly overalpping) pairwise mappings in a 
collection of genomes. The mappings are output in GFF format.

Currently BubbZ does not support chromosomes in the input longer than
4294967296 bp, this will be fixed in the future releases.

Compilation and installation
============================
To compile the code, you need recent installations of the following software
(Linux only):

* Git
* CMake 
* A GCC compiler supporting C++11
* Intel TBB library properly installed on your system. In other words, G++
  should be able to find TBB libs (future releases will not depend on TBB)

Many systems will have this software already installed. If not, you will 
need sudo priviliges to install them. The easiest way to install the 
dependencies is to use a package management system; for APT on Debian systems
they can be installed by the following:

	sudo apt-get install git cmake g++ libtbb-dev

Once you installed the things above, do the following:

Clone the repository https://github.com/medvedevgroup/BubbZ by
running:

	git clone https://github.com/medvedevgroup/BubbZ 

Go to the root directory of the project and create the "build" folder by
executing:

	cd BubbZ
	mkdir build

Initialize dependencies by executing:

	git submodule update --init --recursive

Go to the "build" directory and compile and install the project by running:
	
	cd build
	cmake .. -DCMAKE_INSTALL_PREFIX=<path to install the binaries>
	make install

The make run will produce and install the executables of twopaco, BubbZ-lcb,
spoa and a wrapper script BubbZ which implements the pipeline.

BubbZ usage
==============
BubbZ takes a collection of FASTA file as an input. The simplest way to run
BubbZ is to run the following command:

	bubbz <input FASTA files>

For example:

	bubbz genome1.fa genome2.fa

By default, the output will be written in the directory "bubbZ_out" in
the current working directory. It will contain a GFF file "blocks_coords.gff"
containing coordinates of the mappings. BubbZ has several parameters that
affect the accuracy and performance, they are described below.

Output description
==================
The output directory will contain a GFF file called "blocks_coords.gff" with coordinates of the
locally-collinear blocks. Lines that have identical id fields correspond
to different copies of the same block.

Parameters affecting accuracy
=============================

The value of k
--------------
This parameter defines the order of the de Bruijn graph being used and controls
the tradeoff between the sensitivity on one hand, and speed and memory usage
on the other. The parameter is set by the key

	-k <an odd integer>

In general the lower the k, the slower and more sensitive the alignment is. For
small datasets, like bacteria, we recommend k=15, and for mammalian-sized
genomes k=21. The default is 21.

Vertices frequency threshold
----------------------------
Mammalian genomes contain many repeated elements that make the graph large and
convoluted. To deal with this issue, BubbZ removes all k-mers with frequency
more than a threshold, which is controlled by the option:

	-a <integer>

We recommend to set it to the twice of the maximum number of copies a homologous
block in the input genomes has. For example, if the largest gene family of the input
genomes has N members, set -a to at least N * 2. However, increasing this value may
significantly slow down the computation. The default value is 150.

Gap size threshold
---------------------
BubbZ analyzes the graph by looking for long chains of common vertices in it.
The gap size in a chain is limited by by parameter, which can be set using:

	-b <integer>

The default value of -b is 200. Increasing value may improve recall of divergent
sequences, but if -b is too high, it will decrease accuracy as well.

Mapping block size
----------------------------
BubbZ only output blocks longer than a specified threshold, which is set by

	-m <integer>

The default value is 200. 

Technical parameters
====================

Threads number
--------------
The maximum number of thread for BubbZ to use. This parameter is set by 

	-t <integer>

By default BubbZ tries to use as much threads as possible. You can limit
this number by using the above switch. Note that different stages of the
pipeline have different scalabilities. TwoPaCo will not use more than
16 threads, while graph analyzer BubbZ-lcb will use as much as possible.

Memory allocation
-----------------
The graph constructor TwoPaCo preallocates memory for Bloom filter. By default,
the Bloom filter size is thrice of the size of the input files. The Bloom
filter size can be set manually with the option:

	-f <memory amount in GB>

Output directory
----------------
The directory for the output files can be set by the argument

	-o <directory>

The default is "BubbZ_out" in the current working directory.

A note about the repeat masking
==============================
BubbZ and TwoPaCo currently do not recognize soft-masked characters (i.e. using
lowercase characters), so please convert soft-masked repeats to hard-maksed ones
(with Ns) if you would like to mask the repeats explicitly. However, it is not
necessary as BubbZ uses the abundance parameter -a to filter out high-copy
repeats.

License
=======
See LICENSE.txt

Contacts
========
E-mail your feedback at ivminkin@gmail.com.

Datasets used of analyses in the paper
======================================
See: https://github.com/medvedevgroup/BubbZ/blob/master/DATA.txt
