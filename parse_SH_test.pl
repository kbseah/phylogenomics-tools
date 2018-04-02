#!/usr/bin/perl

# Calls RAxML to infer multi-gene phylogeny from the concatenated alignment
# Infer individual gene trees from individual alignments
# For each gene alignment compare the best gene tree with concatenated tree with SH-test

## TODO: Integrate wtih CONSEL to see how they perform on other tests!
## Make it possible to specify number of starting trees
## Make it possible to choose between different tree analyses

### Packages ###
use strict;
use warnings;
use Cwd qw();
use File::Path qw(make_path);
use Getopt::Long;
use Bio::Align::Utilities qw(:all);
use Bio::SeqIO;
use Bio::Seq;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;

### Global variables ###

my $path = Cwd::cwd();
my $marker_file;
my $path_to_wd;
my $concat_aln_file;
my @markers;
my @shortnames;
my $NUMTHREADS = 4;

### Options ###

if ( !@ARGV ) { usage(); }

GetOptions ("markers=s" => \$marker_file,
            "wd=s" => \$path_to_wd,
            "concat_aln=s" => \$concat_aln_file,
            "numthreads=s" => \$NUMTHREADS,
) or usage();

## MAIN ########################################################################

read_marker_names();
#convert_phylip_aln();
#mkdir "$path_to_wd/trees";
#make_concat_tree();
#make_partitioned_tree();
#calc_indiv_gene_trees();
#perform_SH_test();
parse_SH_output();
#draw_gene_tree_bipartitions();
#draw_consensus_tree();

print STDERR "*** Job complete *** \n";
print STDERR "*** Happy Happy Joy Joy *** \n";

## SUBS ########################################################################

sub usage {
    # print usage statement
    print STDERR "******************************************************* \n";
    print STDERR "Infer trees and perform SH likelihood tests on pairs \n";
    print STDERR "Cite: RAxML, SH-test papers \n";
    print STDERR "KBS 2013-11-13 \n";
    print STDERR "Usage: perl tree_calculations.pl \\ \n";
    print STDERR "\t --markers MARKER_TABLE \\ \n";
    print STDERR "\t --wd WORKING_FOLDER \n";
    print STDERR "\t --concat_aln NAME_OF_CONCATENATED_ALIGNMENT_FILE \\ \n";
    print STDERR "\t --numthreads NUMBER_OF_RAXML_THREADS (at least 2!)  \n";
    print STDERR "******************************************************* \n";
    exit;
}

sub read_marker_names {
    # Read list of marker genes for extraction
    open (INFO, "< $marker_file") or die ("Cannot open list of marker genes: $!");
        while (<INFO>) {
            chomp;
            push (@markers, $_);
        }
    close (INFO);
}

sub parse_SH_output {
    # Parse SH test output from RAxML output files to a summary file
    ## TODO: Parse output into a table
    system ("cat /dev/null > $path_to_wd/trees/SH_TEST_OUTPUT");
    print STDERR "Copying SH test results to file $path_to_wd/trees/SH_TEST_OUTPUT \n";
    foreach my $marker (@markers) {
        open (my $INPUT, "< $path_to_wd/trees/RAxML_info\.$marker\_SH_test") or die ("Cannot open RAxML SH test result file: $!");
        my @data = <$INPUT>;
        close ($INPUT);
        my $totallines = scalar @data;
        my $totallines2 = $totallines - 1;
        my $totallines3 = $totallines - 2;
        open (my $OUTPUT, ">> $path_to_wd/trees/SH_TEST_OUTPUT");
        print $OUTPUT "gene: $marker \n";
        print $OUTPUT "$data[$totallines3]";
        print $OUTPUT "$data[$totallines2]\n";
        close ($OUTPUT);
    }
}

sub draw_gene_tree_bipartitions {
    # Writes the % of bipartitions found on set of all indiv gene trees onto the concat tree, much like a bootstrap %
    chdir "$path_to_wd/trees";
    # Concatenate all indiv gene trees into single tree file
    system ("cat *bestTree.*besttree > gene_trees_allset");
    print STDERR "Writing bipartitions found in individual gene trees to concatenated alignment tree... \n";
    system ("raxmlHPC-PTHREADS -T $NUMTHREADS -m PROTGAMMAWAG -s ../alignments/$concat_aln_file -n bipartitions_write -t RAxML_bestTree.concat_unpart -z gene_trees_allset -f b &> /dev/null");
    print STDERR "Bipartition writing done \n";
    print STDERR "************************************************** \n";
    chdir "$path";
}

sub draw_consensus_tree {
    # Draws an MRE consensus tree based on the individual gene trees
    chdir "$path_to_wd/trees";
    print STDERR "Calculating MRE consensus from gene trees...";
    system ("cat *bestTree.*besttree > gene_trees_allset");
    system ("raxmlHPC-PTHREADS -T $NUMTHREADS -m PROTGAMMAWAG -n consMRE -s ../alignments/concatenated.phy -J MRE -z gene_trees_allset &> /dev/null");
    print STDERR "MRE consensus calculated from individual gene trees \n";
    print STDERR "*************************************************** \n";
    chdir "$path";

}
