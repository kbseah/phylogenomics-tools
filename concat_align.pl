#!/usr/bin/perl

=head1 NAME

concat_align.pl - Concatenate alignment of marker genes

=head1 SYNOPSIS

    perl concat_align.pl --file <filename table> \
                         --wd <working directory> \
                         --markers <marker_table> \
                         --out <output prefix> \
                         --mask \
                         --model_select

    perl concat_align.pl --help

=head1 DESCRIPTION

Concatenate alignments of each marker gene listed in marker table to a single
alignment file, with output prefix for file name.

=head1 ARGUMENTS

=over 8

=item --file <file>

Shortnames and filenames table.

=item --wd <string>

Path to working directory

=item --markers <file>

Table of marker gene names for further processing.

=item --out <string>

Prefix for output file

=item --mask

Logical: Use Zorro to mask alignment columns by alignment quality? (Default: No)

=item --model_select

Logical: Use automated model selection? (Slow) (Default: No)

=item --help|-h

This help message

=back

=head1 OUTPUT

Alignment in Fasta format with filename given by --out, in working directory.

=head1 COPYRIGHT AND LICENSE

phylogenomics-tools. Copyright (C) 2013 Brandon Seah (kbseah@mpi-bremen.de)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

=cut

use strict;
use warnings;
use Cwd qw();
use Getopt::Long;
use Bio::Align::Utilities qw(:all);
use Bio::SeqIO;
use Bio::Seq;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Pod::Usage;

### Initialize variables

my @shortnames;
my @markers;

my $species_file;    # Species names file
my $marker_file;    # Marker genes list
my $path_to_wd;        # Working directory
my $output_aln_prefix;    # Prefix of concatenated alignment and partition file
my $use_mask = 0;    # Option whether to use the Zorro masking files
my $model_choose = 0;    # Option whether to perform automated model selection, using code from Alexis Stamatakis
my $raxmlExecutable = "raxmlHPC-PTHREADS -T 4";    # Which RAxML version to use?

my $path = Cwd::cwd();    # Print current working directory to $path

if ( ! @ARGV ) { pod2usage (-message=>"Insufficient input options",-exitstatus=>2); }

### Get options

GetOptions ("file=s" => \$species_file,
            "markers=s" => \$marker_file,
            "wd=s" => \$path_to_wd,
            "out=s" => \$output_aln_prefix,
            "mask" => \$use_mask,
            "model_select" => \$model_choose,
            "help|h" => sub {pod2usage (-exitstatus=>2,-verbose=>2); },
            ) or pod2usage (-message=>"Please check input options",-exitstatus=>2,-verbose=>2);

## MAIN ########################################################################

read_shortnames();
read_marker_names();
my $concat_len = scalar @markers;            # Get total number of genes in alignment

my $concat_aln;        # Bio::SimpleAlign object that will hold the concatenated alignment
my %concat_hash;    # temporary hash to store the alignments with marker names as keys
my %gene_positions;    # Hash to store the positions of each marker in the concatenated alignment, later to output partition file
my %best_model_choice;    # Hash to store the best models chosen by the automated model selection procedure, for each marker

my $out_fasta = Bio::AlignIO->new(-file=>">$path_to_wd/alignments/$output_aln_prefix\.aln",-format=>'fasta');    # define the Multifasta output file
my $out_phylip = Bio::AlignIO->new(-file=>">$path_to_wd/alignments/$output_aln_prefix\.phy",-format=>'phylip');    # define the Phylip output file

read_alignments();

if ($model_choose == 1) {
    print STDERR "Performing model test (could take a while...) \n";
    model_choice_wrapper();
}

print STDERR "Writing concatenated alignment to file \n";
write_concat_alignments();
$out_fasta -> write_aln($concat_aln);    # Write the alignments to file defined above
$out_phylip -> write_aln($concat_aln);    # Write phylip formatted alignment

print STDERR "*** Job complete ***\n";
print STDERR "Concatenated Fasta alignment written to ./$path_to_wd/alignments/$output_aln_prefix\.aln \n";
print STDERR "Concatenated Phylip alignment written to ./$path_to_wd/alignments/$output_aln_prefix\.phy \n";
print STDERR "Gene partitions file written \./$path_to_wd/alignments/$output_aln_prefix\.partitions \n";
print STDERR "************ \n";


### Subroutines

sub usage {                     # print usage statement
    print STDERR "******************************************************* \n";
    print STDERR "Concatenate alignments of marker genes \n";
    print STDERR "KBS 2013-11-11 \n";
    print STDERR "Usage: perl concat_align.pl \\ \n";
    print STDERR "\t --file GENOME_RENAME_TABLE \\ \n";
    print STDERR "\t --markers MARKER_TABLE \\ \n";
    print STDERR "\t --wd WORKING_FOLDER \\ \n";
    print STDERR "\t --mask [Use mask files produced by Zorro] \\ \n";
    print STDERR "\t --model_select [Use automated model selection (slow)] \\ \n";
    print STDERR "\t --out CONCATENATED_ALN_PREFIX \n";
    print STDERR "******************************************************* \n";
    exit;
}


sub read_shortnames {            # Read list of species shortnames
    open (INFO, "< $species_file") || die ("Cannot open list of species names: $!");
    while (<INFO>) {
        chomp;
        (my $throwaway, my $current_shortname) = split(/\t/,$_);
        push (@shortnames, $current_shortname);
    }
    close (INFO);
}

sub read_marker_names {            # Read list of marker genes for extraction
    open (INFO, "< $marker_file") || die ("Cannot open list of marker genes: $!");
    while (<INFO>) {
        chomp;
        push (@markers, $_);
    }
    close (INFO);
}



sub read_alignments {
    # Read alignments into hash %concat_hash with marker names as keys
    foreach my $marker (@markers) { # for each marker gene in the list
        # read the multi-Fasta alignment
        my $in = Bio::AlignIO->new(-file=>"$path_to_wd/alignments/$marker\.cat\.aln",-format=>'fasta');
        my $thein = $in->next_aln();
        # define the Phylip output file
        my $theout_phylip = Bio::AlignIO->new(-file=>">$path_to_wd/alignments/$marker\.cat\.phy",-format=>'phylip');
        $theout_phylip -> write_aln($thein);
        # add it to the temporary hash, indexed by marker name
        $concat_hash{$marker}=$thein;
    }
}


# This is a very ugly hack
sub write_concat_alignments {
    # Concatenate the alignments into new alignment file
    # Can't recursively cat on empty string, so must start from first element
    $concat_aln = cat($concat_hash{$markers[0]});
    $gene_positions{$markers[0]}= $concat_aln->length();

    # Initialize a file to contain the list of marker genes and their positions in the alignment
    open (my $PARTITION, "> $path_to_wd/alignments/$output_aln_prefix\.partitions") or die ("Cannot write partitions file: $!");
    if ( $model_choose == 0 ) {
    print $PARTITION join("" , "WAG" , ", " , $markers[0] , " = ", join("-","1",$gene_positions{$markers[0]})), "\n";
    } elsif ( $model_choose == 1 ) {
        print $PARTITION join("" , $best_model_choice{$markers[0]} , ", ", $markers[0], " = ", join("-","1",$gene_positions{$markers[0]})), "\n";
    }
    close ($PARTITION);

    if ( $use_mask == 1 ) {
        system ("cat $path_to_wd/alignments/$markers[0]\.cat\.mask >> $path_to_wd/alignments/$output_aln_prefix\.mask");
    } else {}

    for (my $x = 1; $x <= ($concat_len - 1); $x++){
        # Iterate through the hash and recursively concatenate the alignments
        $concat_aln = cat($concat_aln, $concat_hash{$markers[$x]});
        # store the end position of this marker gene in the alignment
        $gene_positions{$markers[$x]}= $concat_aln->length();
        # Add the marker and its positions to the partitions file for later use with phylogenetic software
        open (my $PARTITION, ">> $path_to_wd/alignments/$output_aln_prefix\.partitions") or die ("Cannot write partitions file: $!");
        if ( $model_choose == 0 ) {
            print $PARTITION join("", "WAG", ", " , $markers[$x], " = ", join("-",($gene_positions{$markers[$x-1]} + 1),$gene_positions{$markers[$x]}) ), "\n";
        } elsif ( $model_choose == 1 ) {
            print $PARTITION join("", $best_model_choice{$markers[$x]}, ", " , $markers[$x], " = ", join("-",($gene_positions{$markers[$x-1]} + 1),$gene_positions{$markers[$x]}) ), "\n";
        }
        close ($PARTITION);
        if ( $use_mask == 1 ) {
            system ("cat $path_to_wd/alignments/$markers[$x]\.cat\.mask >> $path_to_wd/alignments/$output_aln_prefix\.mask");
        }
    }
}


sub getLH        # Subroutine required for model test script by Alexis Stamatakis
  {
    my $fileID = $_[0];
    open(CPF, $fileID);
    my @lines = <CPF>;
    close(CPF);
    my $numIT = @lines;
    my $lastLH = pop(@lines);
    my $k = index($lastLH, '-');
    my $LH = substr($lastLH, $k);
    return $LH;
  }

#sub getTIME
#  {
#    my $fileID = $_[0];
#    open(CPF, $fileID);
#    my @lines = <CPF>;
#    close(CPF);
#    my $numIT = @lines;
#    my $lastLH = pop(@lines);
#    my $k = index($lastLH, '-');
#    my $TIME = substr($lastLH, 0, $k-1);
#    return $TIME;
#  }

sub choose_prot_model {        # Adapted from script by Alexis Stamatakis
    my $alignmentName = $_[0]; # Read the MARKER name
    my $UNLIKELY = -1.0E300;
    my @lh;

    my @AA_Models = ("DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG", "RTREV",
                     "CPREV", "VT", "BLOSUM62", "MTMAM", "LG", "MTART", "MTZOA",
                     "PMB", "HIVB", "HIVW", "JTTDCMUT", "FLU", "DAYHOFFF",
                     "DCMUTF", "JTTF", "MTREVF", "WAGF", "RTREVF", "CPREVF",
                     "VTF", "BLOSUM62F", "MTMAMF", "LGF", "MTARTF", "MTZOAF",
                     "PMBF", "HIVBF", "HIVWF", "JTTDCMUTF", "FLUF");


    #print "Determining AA model data\n";
    #print "Computing randomized stepwise addition starting tree number :".$i."\n";
    my $cmd = $raxmlExecutable." -y -p 12345 -m PROTCATJTT -s $path_to_wd/alignments/$alignmentName\.cat\.phy -n ST_".$alignmentName." \> ST_".$alignmentName."_out";
    system($cmd);

    my $numberOfModels = @AA_Models;

    for (my $i = 0; $i < $numberOfModels; $i++) {
        my $aa = "PROTGAMMA".$AA_Models[$i];
        my $cmd = $raxmlExecutable." -f e -m ".$aa." -s $path_to_wd/alignments/$alignmentName\.cat\.phy -t RAxML_parsimonyTree.ST_".$alignmentName." -n ".$AA_Models[$i]."_".$alignmentName."_EVAL \> ".$AA_Models[$i]."_".$alignmentName."_EVAL.out\n";
        #print($cmd);
        system($cmd);
    }


    for (my $i = 0; $i < $numberOfModels; $i++) {
        my $logFileName = "RAxML_log.".$AA_Models[$i]."_".$alignmentName."_EVAL";
        #print $logFileName."\n";
        $lh[$i] = getLH($logFileName);
    }

    my $bestLH = $UNLIKELY;
    my $bestI = -1;

    for (my $i = 0; $i < $numberOfModels; $i++) {
        #print "Model: ".$AA_Models[$i]." LH: ". $lh[$i]."\n";
        if($lh[$i] > $bestLH) {
            $bestLH = $lh[$i];
            $bestI = $i;
        }
    }

#    print "Best Model : ".$AA_Models[$bestI]."\n\n";
#    my $bestModel = $AA_Models[$bestI];
    $best_model_choice{$alignmentName} = $AA_Models[$bestI];    # Write best model chosen to the hash of best model choices
}

sub model_choice_wrapper {        # Wrapper routine to perform model test for each marker gene in the list
    foreach my $marker (@markers) {
        choose_prot_model($marker);
    }
}
