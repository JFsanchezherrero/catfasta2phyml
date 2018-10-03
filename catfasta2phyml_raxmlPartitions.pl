#!/usr/bin/perl 


## 09/16/2015 04:57:30 PM
## TODO: print strict PHYLIP output if not -s
##   ./catfasta2phyml-branch.pl -v -c -p testing2/* > outfile

# 	Update: 09/05/2018 
#  	DESCRIPTION:  add RAXML partition file
#	AUTHOR:	 Jose F. Sanchez herrero

#	Update: 10/03/2018 
#	DESCRIPTION:  Add ProtTest for best model selection for each aln added to the partition
#	AUTHOR:	 Jose F. Sanchez herrero

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;
Getopt::Long::Configure("bundling_override");

#---------------------------------------------------------------------------
#  Global variables
#---------------------------------------------------------------------------
# my %HoH = (
#     'file' => {
#                 'aligned' => 1,
#                 'nseqs'   => 10,
#                 'nchars'  => 4,
#                 'seqs'    => {
#                               'apa' => 'ACGT',
#                               'bpa' => 'ACGT',
#                           }
# 
#               },
# 
# );
## ${HoH}{'file'}{'seqs'}{'apa'}
my %HoH              = ();
my %seqids           = ();   # 
my %nseq_hash        = ();   # key:infile, val:nseq
my %nchar_hash       = ();   # key:infile, val:nchar (for aligned data)
my %seqid_count_hash = ();   # key:seqid, val:count
my $term             = $/;   # input record separator
my @hash_ref_array   = ();   # array with hash references
my @all_labels_found = ();   # All unique labels in all files
my @infiles          = ();   # Infiles on ARGV
my $nfiles           = 0;    # count number of files
my $space            = "\t"; # spacer for aligned print
my $nchar            = 0;    # nchar for phyml header.
my $nseq             = 0;    # nseq for phyml header.
my $fasta            = 0;    # Print phyml format by default
my $concatenate      = 0;    # Concatenate after adding missing taxa
my $man              = 0;    # Manual
my $help             = 0;    # Help
my $verbose          = 0;    # Verbose
my $noprint          = 0;    # Do not print the concatenation
my $sequential       = 0;    # Print sequential with line breaks in in sequence (default is interleaved) 
my $strict_phylip    = 0;    # Print strict phylip format (http://evolution.genetics.washington.edu/phylip/doc/sequence.html)
my $lwidth           = 60;   # default line width for fasta
my $prottest_jar_file = 0;
my $Noprottest;
my $Model_partition;

#---------------------------------------------------------------------------
#  Handle arguments
#---------------------------------------------------------------------------
if (@ARGV < 1) {
    die "No arguments. Try:\n\n $0 -man\n\n";
}
else {
    GetOptions('h|help|?'         => sub { pod2usage(1) },
               'm|man'            => sub { pod2usage(-exitstatus => 0, -verbose => 2) },
               'f|fasta'          => \$fasta,
               'c|concatenate'    => \$concatenate,
               'v|verbose'        => \$verbose,
               'n|noprint'        => \$noprint,
               's|sequential'     => \$sequential,
               'p|phylip'         => \$strict_phylip,
               'protest_jar'	  => \$prottest_jar_file,
               'noProtTest'		  => \$Noprottest,
               'model_partition'  => \$Model_partition,
              );
}


#---------------------------------------------------------------------------
#  Read all infiles, count labels and get sequence lengths
#---------------------------------------------------------------------------
print STDERR "\nChecking sequences in infiles...\n\n" if ($verbose);

foreach my $infile (@ARGV) {

    my $seq_hash_ref = parse_fasta($infile); # key:seqid, value:sequence
    print STDERR "  File $infile: " if ($verbose);
    push (@infiles, $infile);
    $nfiles++;

    ## Are sequences aligned?
    my (@ret) = aligned($seq_hash_ref);
    if (scalar(@ret) > 1) {
        print STDERR "Error! Expecting aligned input sequences.\n";
        print STDERR "Sequences in $infile are not all of the same length ($ret[1] is $ret[0], $ret[2] is $ret[3])\n";
        exit(0);
    }
    elsif (scalar(@ret) == 1) {
        $HoH{$infile}{'aligned'} = 1;
        $HoH{$infile}{'nchars'} = $ret[0];
    }

    ## Save sequences from file(s) in one hash. Mem limit?
    foreach my $key (keys %$seq_hash_ref) {
    	my $string2;
    	if (length $key > 9) {
  	        $string2 = substr($key, 0, 9);
    	} else {
	    	$string2 = $key;
    	}
        $HoH{$infile}{'seqs'}{$key} = ${$seq_hash_ref}{$key};
    }
    ## Get nseqs for file
    $HoH{$infile}{'nseqs'} = scalar(keys %{${HoH}{$infile}{'seqs'}});
    $nseq_hash{$infile} = $HoH{$infile}{'nseqs'};

    ## Get length of sequence for all tax labels.
    foreach my $tax_key (keys %$seq_hash_ref) {
        $seqid_count_hash{$tax_key}++;
    }

    print STDERR " ntax=$HoH{$infile}{'nseqs'} nchar=$HoH{$infile}{'nchars'}\n" if ($verbose);

} # Done with file


#---------------------------------------------------------------------------
# Collect all taxon labels and find the length of longest name for printing
#---------------------------------------------------------------------------
(@all_labels_found) = sort { length($b) <=> length($a) } keys %seqid_count_hash;
$space = length($all_labels_found[0]) + 2;
#print Dumper @all_labels_found."\n";
#exit();


#---------------------------------------------------------------------------
# Get nseq. First check if nseqs are equal among input files
#---------------------------------------------------------------------------
my %string = map { $_, 1 } values %nseq_hash;
if (keys %string == 1) { # all values equal
    my (@nseqs) = values %nseq_hash;
    $nseq = shift(@nseqs);
}
elsif ($concatenate) {
    foreach my $file (keys %HoH) {
        my %second = map {$_=>1} (keys %{${HoH}{$file}{'seqs'}}); # Get all seq labels in file,
        my @only_in_all_labels_found = grep { !$second{$_} } @all_labels_found; # and compare with all labels found in all files
        if (@only_in_all_labels_found) {
            print STDERR "\n  Need to fill missing sequences for file $file\n" if ($verbose);
            my $allgapseq = '-' x $HoH{$file}{'nchars'};
            foreach my $seqid (@only_in_all_labels_found) {
                print STDERR "    Adding all gaps for seqid $seqid\n" if ($verbose);
                $HoH{$file}{'seqs'}{$seqid} = $allgapseq;
            }
        }
        $nseq = scalar(keys %{$HoH{$file}{'seqs'}});
        $HoH{$file}{'nseqs'} = $nseq;
    }
}
else {
    print STDERR "\n";
    foreach my $key (sort { $seqid_count_hash{$a} <=> $seqid_count_hash{$b} } (keys %seqid_count_hash)) {
        printf STDERR "%-${space}s --> %d\n", $key, $seqid_count_hash{$key};
    }
    print STDERR "\nError! Some sequence labels does not occur in all files.\n";
    print STDERR "That is, sequence labels needs to be identical for safe concatenation.\n";
    print STDERR "Use the --concatenate (or -c) to concatenate anyway.\nEmpty (all gap) sequences will be added where needed.\n\n";
    exit(1);
}


#---------------------------------------------------------------------------
# Check if names can be abbreviated
#---------------------------------------------------------------------------
if ($strict_phylip) {
    my %test = ();
    foreach my $label (@all_labels_found) {
        my $phylabel = phylip_label($label);
        if ($test{$phylabel}++) {
            print STDERR "\nWarning! Sctrict Phylip format results in duplicate labels for these data (e.g., $phylabel)!\n";
            exit(0);
        }
    }
}

#---------------------------------------------------------------------------
#  Get nchar
#---------------------------------------------------------------------------
my $start_position=1;
my $partition_file = "partitions.txt";
open (PART, ">$partition_file");

my $tmp_part = "tmp_folder";
mkdir $tmp_part, 0755;

## Proteins models available in IQTREE
# LG, DCMut, mtMAM, JTT, mtREV, cpREV, Dayhoff, WAG, mtART, rtREV, VT, HIVb, FLU, Blosum62, HIVw,
# PMB, JTTDCMut, GTR20, mtMet, mtVer, mtInv, Poisson, mtZOA,

# ProtTest: Include matrix (Amino-acid) = JTT LG DCMut MtREV MtMam MtArt Dayhoff WAG RtREV CpREV Blosum62 VT HIVb HIVw FLU 

## Share ProtTest and IQTREE: ATTENTION, different names!
# LG, DCMut, mtMAM, JTT, mtREV, mtART, cpREV, Dayhoff, WAG, rtREV,  VT, HIVb, FLU, Blosum62, HIVw,

## parse $output_tmp
my %conversion = (
## prottest => IQTREE
	'LG' => "LG",
	'JTT' => "JTT",
	'DCMut' => "DCMut",
	'MtREV' => "mtREV",
	'MtMam' => "mtMAM",
	'MtArt' => "mtART",			
	'Dayhoff' => "Dayhoff",
	'WAG' => "WAG",
	'RtREV' => "rtREV",			
	'VT' => "VT",			
	'HIVb' => "HIVw",			
	'Blosum62' => "Blosum62",			
	'CpREV' => "cpREV",			
	'RtREV' => "rtREV",			
);

foreach my $file (keys %HoH) {
    my $length = ($HoH{$file}{'nchars'});
    $nchar = $nchar + $length;

	#===  FUNCTION  ================================================================
	#      VERSION:  09/05/2018 
	#  DESCRIPTION:  add RAXML partition file
	#		AUTHOR:	 Jose F. Sanchez herrero
	#===============================================================================
    ## print partitions
    my $end_position=$start_position+$length-1;
    my @name = split(".fa", $file);
	
	if ($Noprottest) {
		
		print PART "$Model_partition, ".$name[0]." = "."$start_position"."-"."$end_position\n";
	} else {
	
		## add system call to java -jar prottest-3.4.2.jar -i $file
		my $output_tmp = $tmp_part."/".$name[0]."_tmp.txt";
		my $system_call = system("java -jar ".$prottest_jar_file." -i $file -o ".$output_tmp);
		if ($system_call == 1) {print "Error when calling ProtTest java jar file...\n"; exit(); }
		my $best_model;
		open (OUT, "<$output_tmp");
		while (<>) {
			chomp;
			if ($_ =~ /Best model according to BIC: (.*)/) {
				$best_model = $1; last;
		}}
		close (OUT);
		print PART "$conversion{$best_model}, ".$name[0]." = "."$start_position"."-"."$end_position\n";
	}
	$start_position += $length;
}
close (PART);

#---------------------------------------------------------------------------
#  Print everything to STDOUT
#---------------------------------------------------------------------------
if ($verbose) {
    print STDERR "\nChecked $nfiles files -- sequence labels and lengths seems OK.\n";
    print STDERR "Concatenated $nseq sequences, total length $nchar.\n";
    print STDERR "Printing concatenation to STDOUT.\n\n" unless $noprint;
}

if ($noprint) {
    print STDERR "\n\nEnd of script.\n\n" if ($verbose);
    exit(1);
}

if ($strict_phylip) {
    print STDOUT "   $nseq    $nchar\n";
}
else {
    print STDOUT "$nseq $nchar\n" unless $fasta;
}


#---------------------------------------------------------------------------
# Print the hash via intermediate hash (mem limit?)
# TODO: Try to circumvent the intermediate hash
#---------------------------------------------------------------------------
if ($fasta or $sequential) {
    ## First, concatenate all sequences from hashes
    my %print_hash = (); # key:label, value:sequence
    foreach my $file (@infiles) { # Keep input order
        die "Error: $file not in HoH\n" unless exists(${HoH}{$file});
        foreach my $seqid (keys %{${HoH}{$file}{'seqs'}}) {
            $print_hash{$seqid} .= $HoH{$file}{'seqs'}{$seqid}; # Concatenate seqs from files
        }
    }

    ## Then print, and add line breaks in sequences
    foreach my $label (sort keys  %print_hash) {
        if ($fasta) {
            print STDOUT ">$label\n";
            $print_hash{$label} =~ s/\S{$lwidth}/$&\n/gs; # replace word of size $lwidth with itself and "\n"
            print STDOUT $print_hash{$label}, "\n";
        }
        elsif ($strict_phylip) {
            my $phylip_label = phylip_label($label);
            printf STDOUT "%-10s ", $phylip_label;
            my $s = phylip_blocks($print_hash{$label});
            print $s, "\n";
        }
        else {
            printf STDOUT "%-${space}s ", $label;
            print STDOUT $print_hash{$label}, "\n";
        }
    }
}
else { # default: phyml interleaved
    my $did_first = 0;
    foreach my $file (@infiles) { # Keep input order
        die "Error: $file not in HoH\n" unless exists(${HoH}{$file});
        foreach my $seqid (sort keys %{$HoH{$file}{'seqs'}}) {
            ## create an array with the sequence in pieces of 10. If strict phylip, shift 5 pieces when printing the first time
            ## if printing the other times, shift 6 pieces
            ##    my @seq_array = unpack("(A10)*", $HoH{$file}{'seqs'}{$seqid});
            ##    while (my @five_pieces) = splice @seq_array, 0, 5) {
            ##        print Dumper(@seq_array);warn "\n HERE (hit return to continue)\n" and getc();
            ##    }
            if ($strict_phylip) {
                my $phylip_seqid = phylip_label($seqid);
                print STDOUT $phylip_seqid unless $did_first;
                ##my $s = phylip_blocks($HoH{$file}{'seqs'}{$seqid}); # 09/08/2015 01:58:15 PM: does not work with the sequence
                ##print STDOUT $s, "\n";
            }
            else {
                printf STDOUT "%-${space}s ", $seqid unless $did_first;
                ##print STDOUT "$HoH{$file}{'seqs'}{$seqid}\n";
            }
            ## Print sequence
            ## TODO: phylip strict printing of sequence in blocks of 10
            ## TODO: print length of 60
            print STDOUT "$HoH{$file}{'seqs'}{$seqid}\n"; # <<<<<<<<<<<<
        }
        print "\n";
        $did_first = 1;
    } 
}

print STDERR "\nEnd of script.\n\n" if ($verbose);


#===  FUNCTION  ================================================================
#         NAME:  aligned
#      VERSION:  08/31/2015 07:05:54 PM
#  DESCRIPTION:  ???
#   PARAMETERS:  ref to hash with seqs
#      RETURNS:  0 if aligned, array with names and lengths of the first encountered
#                seqs of unequal length otherwise.
#                "Sequences in $infile are not all of the same length ($lname is $length, $name is $l)"
#         TODO:  ???
#===============================================================================
sub aligned {

    my ($h_ref) = shift(@_);

    my $length;
    my $lname;
    my @aligned = ();

    foreach my $name (keys %$h_ref) {
        my $l = length($h_ref->{$name});
        if (defined $length) {
            if ($length != $l) {
                @aligned = ($length,$lname,$name,$l);
                last;
            }
        }
        else {
            $length = length($h_ref->{$name});
            $lname  = $name;
            @aligned = ($length);
        }
    }

    return @aligned;

} # end of aligned


#===  FUNCTION  ================================================================
#         NAME:  parse_fasta
#      VERSION:  11/07/2011 05:25:35 PM
#  DESCRIPTION:  ???
#   PARAMETERS:  filename
#      RETURNS:  hash ref
#         TODO:  ???
#===============================================================================
sub parse_fasta {

    my ($infile) = @_;

    my $term     = $/; # input record separator;
    my %seq_hash = (); # key:seqid, val:seq

    open my $INFILE, "<", $infile or die "could not open infile '$infile' : $! \n"; 
    $/ = ">";
    while(<$INFILE>) {
        chomp;
        next if($_ eq '');
        my ($id, @sequencelines) = split /\n/;
        foreach my $line (@sequencelines) {
            $seq_hash{$id} .= $line;
        }
    }
    $/ = $term;

    return(\%seq_hash);

} # end of parse_fasta


#===  FUNCTION  ================================================================
#         NAME:  phylip_label
#      VERSION:  02/18/2013 04:43:00 PM
#  DESCRIPTION:  manipulates input string to be of length 10, possibly padded with
#                white space.
#   PARAMETERS:  string
#      RETURNS:  string with the length of 10
#         TODO:  ???
#===============================================================================
sub phylip_label {

    my ($string) = @_;

    if (length($string) > 10) {
        $string = substr($string, 0, 10); # truncate label!
    }
    else {
        my $string_length = length($string);
        my $pad = ' ' x ((10 - $string_length)); # pad end with white space
        $string = $string . $pad;
    }

    return($string);

} # end of phylip_label


#===  FUNCTION  ================================================================
#         NAME:  phylip_blocks
#      VERSION:  09/03/2015 10:35:12 PM
#  DESCRIPTION:  return string in blocks of ten characters separated by spaces.
#                No more than 6 blocks wide, and 5 for the first row (giving space
#                to sequence label).
#   PARAMETERS:  string
#      RETURNS:  string
#         TODO:  ???
#===============================================================================
sub phylip_blocks {

    my ($string) = @_;

    my $ret_seq = '';
    my $first   = 1;
    my $i       = 0;

    my @foo = unpack("(A10)*", $string);

    foreach my $p (@foo) {
        $i++;
        $ret_seq .= $p;
        if ($i == 5) {
            if ($first) {
                $ret_seq .= "\n";
                $first = 0;
                $i = 0;
            }
            else {
               $ret_seq .= ' ';
            }
        }
        elsif ($i == 6) {
            $ret_seq .= "\n" unless $first;
            $i = 0;
        }
        else {
            $ret_seq .= ' ';
        }
    }

    return $ret_seq;

} # end of phylip_blocks


#===  POD DOCUMENTATION  =======================================================
#      VERSION:  10/12/2015 10:09:27 AM
#  DESCRIPTION:  Documentation
#         TODO:  ?
#===============================================================================
=pod

=head1 NAME

catfasta2phyml_raxmlPartitions.pl -- Concatenate FASTA alignments PHYLIP, generates partition for RaXML and estimates best model test fitting each partition.

=head1 SYNOPSIS

catfasta2phyml_raxmlPartitions.pl [options] [files]


=head1 OPTIONS

=over 8

=item B<-h, -?, --help>

Print a brief help message and exits.


=item B<-m, --man>

Prints the manual page and exits.


=item B<-c, --concatenate>

Concatenate files even when number of taxa differ among alignments.
Missing data will be filled with all gap (-) sequences.


=item B<-f, --fasta>

Print output in FASTA format. Default is PHYML format.


=item B<-p, --phylip>

Print output in a strict PHYLIP format.
See L<http://evolution.genetics.washington.edu/phylip/doc/sequence.html>.

B<Note:> the current output format is not entirely strict. Left to do is
to efficiently print sequences in blocks of 10 characters.
Check output if using this option.


=item B<-s, --sequential>

Print output in sequential format. Default is interleaved.


=item B<-v, --verbose>

Be verbose by showing some useful output. See the combination with B<-n>.


=item B<-n, --noprint>

Do not print the concatenation, just check if all files have the same sequence lables and lengths.
Program returns 1 on exit. See also the combination with B<-v>.

=item B<-protest_jar>

Java jar file for ProtTest3 (https://github.com/ddarriba/prottest3) in order to generate best model fitting alingment for each partition.

=item B<-noProtTest>

Do not estimate model for each partition. Adds -model_partition [string] to every partition generated in RaXML partition format style.

=item B<-model_partition>

Option provided will be added to every partition generated in RaXML partition format style.


=back

=head1 DESCRIPTION

This is a modification of the original script written by Johan A. A. Nylander, catfasta2phym.pl.

This version is intended to generate a phylip file with as many partitions as alignments provided each one containing the best protein model that fits according to ProtTest3.

Please refer to https://github.com/nylander/catfasta2phyml for further details of the original file.
	
The original catfasta2phyml.pl will concatenate FASTA alignments to one file (interleaved PHYML or FASTA format) after checking that all sequences are aligned (of same length). If there are sequence labels that are not present in all files, a warning will be issued. Sequenced can, however, still be concatenated (and missing sequences be filled with missing data (gaps)) if the argument --concatenate is used.     Prints to STDOUT.

=head1 USAGE

This is a modification of the original script written by Johan A. A. Nylander, catfasta2phym.pl. This version is intended to generate a phylip file with as many partitions as alignments provided each one containing the best protein model that fits according to ProtTest3.

Please refere to https://github.com/nylander/catfasta2phyml for the original file and original options.

To concatenate fasta files to a phyml readable format, generate a partition for each fasta file provided and estimate the best fitting model test.

        catfasta2phyml.pl -p file*.fas -prottest_jar /path/to/prottest_bin/prottest_java_file.jar >> out.phy
        
To concatenate fasta files to a phyml readable format but do not generate a partition for each fasta.

        catfasta2phyml.pl -p file*.fas -noProtTest -model_partition LG

=head1 AUTHOR

Written by Johan A. A. Nylander
Modified by Jose F. Sanchez Herrero (2018)


=head1 DEPENDENCIES

Uses Perl modules Getopt::Long and Pod::Usage


=head1 LICENSE AND COPYRIGHT

Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Johan Nylander. All rights reserved.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details. 
http://www.gnu.org/copyleft/gpl.html 


=head1 DOWNLOAD

https://github.com/nylander/catfasta2phyml
https://github.com/JFsanchezherrero/catfasta2phyml

=cut


__END__
