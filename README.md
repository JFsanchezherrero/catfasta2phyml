catfasta2phyml
==============

NAME

    catfasta2phyml_raxmlPartitions.pl -- Concatenate FASTA alignments PHYLIP, generates partition for 
    RaXML and estimates best model test fitting each partition.

SYNOPSIS

    catfasta2phyml_raxmlPartitions.pl [options] [files]

OPTIONS

    -h, -?, --help
            Print a brief help message and exits.

    -m, --man
            Prints the manual page and exits.

    -c, --concatenate
            Concatenate files even when number of taxa differ among
            alignments. Missing data will be filled with all gap (-)
            sequences.

    -f, --fasta
            Print output in FASTA format. Default is PHYML format.

    -p, --phylip
            Print output in a strict PHYLIP format. See
            http://evolution.genetics.washington.edu/phylip/doc/sequence.htm
            l.
            
            Note: the current output format is not entirely strict. Left
            to do is to efficiently print sequences in blocks of 10 characters.
            Check output if using this option.

    -s, --sequential
            Print output in sequential format. Default is interleaved.

	-output_file
			File to print alignment concatenated.

    -v, --verbose
            Be verbose by showing some useful output. See the combination
            with -n.

    -n, --noprint
            Do not print the concatenation, just check if all files have the
            same sequence lables and lengths. Program returns 1 on exit. See
            also the combination with -v.

	-protest_jar
			Java jar file for ProtTest3 (https://github.com/ddarriba/prottest3)
			in order to generate best model fitting alingment for each partition.
	
	-noProtTest
			Do not estimate model for each partition. Adds -model_partition [string] 
			to every partition generated.
	
	-model_partition            
			String provided will be 
            

DESCRIPTION

	This is a modification of the original script written by Johan A. A. Nylander, 
	catfasta2phym.pl.
	
	This version is intended to generate a phylip file with as many partitions as
	alignments provided each one containing the best protein model that fits according
	to ProtTest3.

	Please refer to https://github.com/nylander/catfasta2phyml for further details
	of the original file.
	
	The original catfasta2phyml.pl will concatenate FASTA alignments to one file
    (interleaved PHYML or FASTA format) after checking that all sequences
    are aligned (of same length). If there are sequence labels that are not
    present in all files, a warning will be issued. Sequenced can, however,
    still be concatenated (and missing sequences be filled with missing data
    (gaps)) if the argument --concatenate is used.

    Prints to STDOUT.

USAGE

	This is a modification of the original script written by Johan A. A. Nylander, 
	catfasta2phym.pl. This version is intended to generate a phylip file with as many 
	partitions as alignments provided each one containing the best protein model that 
	fits according to ProtTest3.
	
	Please refer to https://github.com/nylander/catfasta2phyml for the original file
	and original options.

    To concatenate fasta files to a phyml readable format, generate a partition for 
    each fasta file provided and estimate the best fitting model test.

        catfasta2phyml.pl -p file*.fas -prottest_jar /path/to/prottest_bin/prottest_java_file.jar >> out.phy
        
	To concatenate fasta files to a phyml readable format but do not generate a partition 
	for each fasta.

        catfasta2phyml.pl -p file*.fas -noProtTest -model_partition LG

AUTHOR

    Written by Johan A. A. Nylander
    Modified by Jose F. Sanchez Herrero (2018)

DEPENDENCIES

    Uses Perl modules Getopt::Long and Pod::Usage

LICENSE AND COPYRIGHT

    Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Johan Nylander. All
    rights reserved.

    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details. http://www.gnu.org/copyleft/gpl.html

DOWNLOAD

    https://github.com/nylander/catfasta2phyml
    https://github.com/JFsanchezherrero/catfasta2phyml

