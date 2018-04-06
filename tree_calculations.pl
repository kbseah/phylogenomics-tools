#!/usr/bin/perl

=head1 NAME

tree_calculations.pl

=head1 DESCRIPTION

Calls RAxML to infer multi-gene phylogeny from the concatenated alignment
Infer individual gene trees from individual alignments
For each gene alignment compare the best gene tree with concatenated tree with SH-test


=head1 SYNOPSIS

perl tree_calculations.pl
        --markers MARKER_TABLE
	--wd WORKING_FOLDER 
	--concat_aln concatenated_alignment_file_prefix 
	--numthreads 4 

perl tree_calculations.pl --help

perl tree_calculations.pl --man

=cut 


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
use Pod::Usage;

### Global variables ###

my $path = Cwd::cwd();
my $marker_file;
my $path_to_wd;
my $concat_aln_prefix;
my @markers;
my @shortnames;
my $NUMTHREADS = 4;
my $use_mask = 0;
my $do_convert_phylip_aln = 0;
my $do_convert_nexus_aln = 0;
my $do_calc_indiv_gene_trees = 0;
my $do_perform_SH_test = 0;

### Options ###

if ( !@ARGV ) { pod2usage(-exitstatus=>2, -verbose=>0); }

GetOptions (
	"markers=s" => \$marker_file,
	"wd=s" => \$path_to_wd,
	"concat_aln=s" => \$concat_aln_prefix,
	"numthreads=s" => \$NUMTHREADS,
	"mask" => \$use_mask,
	"convert_phylip|p" => \$do_convert_phylip_aln,
	"convert_nexus|n" => \$do_convert_nexus_aln,
	"indiv|i" => \$do_calc_indiv_gene_trees,
	"sh_test|t" => \$do_perform_SH_test,
	"help|h" => sub {pod2usage(-exitstatus=>2, -verbose=>1)},
	"man|m" => sub {pod2usage(-exitstatus=>0, -verbose=>2)},
) or pod2usage(-exitstatus=>2, -verbose=>0);

=head1 ARGUMENTS

=over 8

=item --markers <file>

Tab-separated file containing list of all the gene markers contained in the
sequence alignment.

=item --wd <path>

Path to working directory (Default: Current folder)

=item --concat_aln <string>

Filename prefix for the concatenated sequence alignment

=item --numthreads <integer>

Number of threads for RAxML to use

=item --mask

Use alignment masking? (Default: No)

=item --convert_phylip|-p

Convert individual Fasta-formatted gene alignments to Phylip format 

=item --convert_nexus|-n

Convert individual Fasta-formatted gene alignments to Nexus format 

=item --indiv|-i

Compute individual gene trees, in addition to tree from concatenated alignment 

=item --sh_test|-t

=item --help|-h

Help message 

=item --man|-m

Manual page.

=back

=cut

### Main code block ###

read_marker_names();
if ($do_convert_phylip_aln == 1) { convert_phylip_aln(); }
if ($do_convert_nexus_aln == 1) { convert_nexus_aln(); }
mkdir "$path_to_wd/trees";
make_concat_tree();
make_partitioned_tree();
if ($do_calc_indiv_gene_trees == 1) {
	calc_indiv_gene_trees(); 
	draw_gene_tree_bipartitions(); 
	draw_consensus_tree();
}
if ($do_perform_SH_test == 1) { 
	perform_SH_test();
	parse_SH_output();
}

print STDERR "************************************ \n";
print STDERR "*********** JOB COMPLETE *********** \n";
print STDERR "******* Happy Happy Joy Joy! ******* \n";
print STDERR "************************************ \n";

### Subroutines ###

sub usage {                     # print usage statement
        print STDERR "******************************************************* \n";
        print STDERR "Infer trees and perform SH likelihood tests on pairs \n";
	print STDERR "Cite: RAxML, SH-test papers \n";
        print STDERR "KBS 2013-11-13 \n";
        print STDERR "Usage: perl tree_calculations.pl \\ \n";
        print STDERR "\t --markers MARKER_TABLE \\ \n";
        print STDERR "\t --wd WORKING_FOLDER \n";
	print STDERR "\t --concat_aln CONCATENATED_ALIGNMENT_PREFIX \\ \n";
	print STDERR "\t --numthreads NUMBER_OF_RAXML_THREADS (at least 2!) \\ \n";
	print STDERR "\t --mask [Use alignment mask?] \\ \n";
	print STDERR "\t -p [Convert indiv gene Fasta alignments to Phylip] \\ \n";
	print STDERR "\t -n [Convert indiv gene Fasta alignments to Nexus] \\ \n";
	print STDERR "\t -i [Compute individual gene trees] \\ \n";
	print STDERR "\t -t [Perform SH test] \\ \n";
        print STDERR "******************************************************* \n";
        exit;
}

sub read_marker_names {			# Read list of marker genes for extraction 
	open (INFO, "< $marker_file") || die ("Cannot open list of marker genes: $!");
		while (<INFO>) {
			chomp;
			push (@markers, $_);
		}
	close (INFO);
}


sub convert_phylip_aln {		# Convert Fasta alignments to Phylip format
	print STDERR "************************************ \n";
	print STDERR "Converting Fasta alignments to Phylip format for each gene: \n";
	print STDERR "************************************ \n";
        foreach my $marker (@markers) {                                                                 # for each marker gene in the list
		print STDERR "$marker \t";
                my $in = Bio::AlignIO->new(-file=>"$path_to_wd/alignments/$marker\.cat\.aln",-format=>'fasta'); # read the multi-Fasta alignment
                my $thein = $in->next_aln();
		my $out = Bio::AlignIO->new(-file=>">$path_to_wd/alignments/$marker\.cat\.phy",-format=>'phylip');
		$out -> write_aln($thein);
	}
	print STDERR "************************************ \n";
	print STDERR "Alignment format conversion complete \n";
	print STDERR "************************************ \n";
}

sub convert_nexus_aln {
	print STDERR "************************************ \n";
	print STDERR "Converting Fasta alignments to Nexus format for each gene: \n";
	print STDERR "************************************ \n";
	foreach my $marker (@markers) {
		print STDERR "$marker \t"; 
		my $in = Bio::AlignIO->new(-file=>"$path_to_wd/alignments/$marker\.cat\.aln",-format=>'fasta'); # read the multi-Fasta alignment
		my $thein = $in->next_aln();
		my $out = Bio::AlignIO->new(-file=>">$path_to_wd/alignments/$marker\.cat\.nxs",-format=>'nexus');
		$out -> write_aln($thein);
	}
}

sub make_concat_tree {			# Infer best tree for entire concatenated alignment
	chdir "$path_to_wd/trees";
	print STDERR "************************************ \n";
	print STDERR "Treeing concatenated alignment \n";
	print STDERR "This will take some time \n";
	print STDERR "************************************ \n";
	if ( $use_mask == 0 ) {
		my $make_concat_tree_exitcode = system ("raxmlHPC-PTHREADS -T $NUMTHREADS -m PROTCATWAG -s ../alignments/$concat_aln_prefix\.phy -n concat_unpart -N 10 -p 12345");}
	elsif ( $use_mask == 1 ) {
		my $make_concat_tree_exitcode = system ("raxmlHPC-PTHREADS -T $NUMTHREADS -m PROTCATWAG -s ../alignments/$concat_aln_prefix\.phy -n concat_unpart -N 10 -p 12345 -a ../alignments/$concat_aln_prefix\.mask");}
#	print STDERR "$make_concat_tree_exitcode \n";
	print STDERR "\n";
	print STDERR "************************************ \n";
	print STDERR "Concatenated alignment tree complete \n";
	print STDERR "************************************ \n";
	# Calculate support values for this tree
	print STDERR "Calculating SH-like support values for unpartitioned concatenated tree \n";
	if ( $use_mask == 0 ) {
		system ("raxmlHPC-PTHREADS -T $NUMTHREADS -m PROTCATWAG -s ../alignments/$concat_aln_prefix\.phy -n concat_unpart_SH -p 12345 -f J -t RAxML_bestTree.concat_unpart");}
	elsif ( $use_mask == 1 ) {
		system ("raxmlHPC-PTHREADS -T $NUMTHREADS -m PROTCATWAG -s ../alignments/$concat_aln_prefix\.phy -n concat_unpart_SH -p 12345 -f J -a ../alignments/$concat_aln_prefix\.mask -t RAxML_bestTree.concat_unpart"); }
	print STDERR "************************************ \n";
	print STDERR "SH-like support values calculated for unpartitioned tree \n";
	print STDERR "************************************ \n";
	print STDERR "\n";
	chdir "$path";
}

sub make_partitioned_tree { 		# Infer best tree under partitioned models for each gene in the concatenated alignment
	chdir "$path_to_wd/trees";
	print STDERR "************************************ \n";
	print STDERR "Treeing concatenated alignment, model partitioned by marker genes \n";
	print STDERR "Tree toplogy and branch lengths estimated jointly \n";
	print STDERR "This will take quite some time... \n";
	print STDERR "************************************ \n";
	if ( $use_mask == 0 ) {
		my $make_partitioned_tree_exitcode = system ("raxmlHPC-PTHREADS -T $NUMTHREADS -m PROTCATWAG -s ../alignments/$concat_aln_prefix\.phy -n concat_part -N 10 -p 12345 -q ../alignments/$concat_aln_prefix\.partitions"); }
	elsif ( $use_mask == 1 ) {
		my $make_partitioned_tree_exitcode = system ("raxmlHPC-PTHREADS -T $NUMTHREADS -m PROTCATWAG -s ../alignments/$concat_aln_prefix\.phy -n concat_part -N 10 -p 12345 -q ../alignments/$concat_aln_prefix\.partitions -a ../alignments/$concat_aln_prefix\.mask"); }
#	print STDERR "$make_partitioned_tree_exitcode \n";
	print STDERR "************************************ \n";
	print STDERR "Partitioned tree complete \n";
	print STDERR "************************************ \n";
	# Calculate support values for this tree
	print STDERR "************************************ \n";
	print STDERR "Calculating SH-like support values for partitioned concatenated tree \n";
	print STDERR "************************************ \n";
	if ( $use_mask == 0 ) {
		system ("raxmlHPC-PTHREADS -T $NUMTHREADS -m PROTCATWAG -s ../alignments/$concat_aln_prefix\.phy -n concat_part_SH -p 12345 -f J -t RAxML_bestTree.concat_part -q ../alignments/$concat_aln_prefix\.partitions");}
	elsif ( $use_mask == 1 ) {
		system ("raxmlHPC-PTHREADS -T $NUMTHREADS -m PROTCATWAG -s ../alignments/$concat_aln_prefix\.phy -n concat_part_SH -p 12345 -f J -t RAxML_bestTree.concat_part -q ../alignments/$concat_aln_prefix\.partitions -a ../alignments/$concat_aln_prefix\.mask");}
	print STDERR "************************************ \n";
	print STDERR "SH-like support values calculated for partitioned tree \n";
	print STDERR "************************************ \n";
	print STDERR "\n";
	chdir "$path";
}

sub calc_indiv_gene_trees {		# Infer best trees for each separate gene alignment
	chdir "$path_to_wd/trees";
	print STDERR "************************************ \n";
	print STDERR "Treeing individual genes \n";
	print STDERR "This will also take some time \n";
	print STDERR "************************************ \n";
	foreach my $marker (@markers) {
		print STDERR "************************************ \n";
		print STDERR "Calculating best tree for gene $marker \n";
		if ( $use_mask == 0 ) {
			system ("raxmlHPC-PTHREADS -T $NUMTHREADS -m PROTCATWAG -s ../alignments/$marker\.cat\.phy -n $marker\_besttree -N 10 -p 12345"); }
		elsif ( $use_mask == 1 ) {
			system ("raxmlHPC-PTHREADS -T $NUMTHREADS -m PROTCATWAG -s ../alignments/$marker\.cat\.phy -n $marker\_besttree -N 10 -p 12345 -a ../alignments/$marker\.cat\.mask"); }
	}
	print STDERR "************************************ \n";
	print STDERR "Individual gene trees complete \n";
	print STDERR "********************************* \n";
	print STDERR "\n";
	chdir "$path";
}

sub perform_SH_test {			# Perform SH test for each gene alignment to compare its own best tree and the best tree for the concatenated alignment
	chdir "$path_to_wd/trees";
	print STDERR "************************************ \n";
	print STDERR "Performing SH-test between best gene tree and concatenated tree for each gene \n";
	print STDERR "************************************ \n";
	foreach my $marker (@markers) {
	print STDERR "************************************ \n";
		print STDERR "Gene: $marker \n";
		if ( $use_mask == 0 ) {
			system ("raxmlHPC-PTHREADS -T $NUMTHREADS -m PROTGAMMAWAG -s ../alignments/$marker\.cat\.phy -n $marker\_SH_test -t RAxML_bestTree.$marker\_besttree -z RAxML_bestTree.concat_unpart -f H"); }
		elsif ( $use_mask == 1 ) { 
			system ("raxmlHPC-PTHREADS -T $NUMTHREADS -m PROTGAMMAWAG -s ../alignments/$marker\.cat\.phy -n $marker\_SH_test -t RAxML_bestTree.$marker\_besttree -z RAxML_bestTree.concat_unpart -f H -a ../alignments/$marker\.cat\.mask"); }
	}
	chdir "$path";
	print STDERR "************************************ \n";
	print STDERR "SH-tests complete\n";
	print STDERR "********************************* \n";
	print STDERR "\n";
}

sub parse_SH_output {			# Parse SH test output from RAxML output files to a summary file 
## TODO: Parse output into a table
	system ("cat /dev/null > $path_to_wd/trees/SH_TEST_OUTPUT");
	print STDERR "************************************ \n";
	print STDERR "Copying SH test results to file $path_to_wd/trees/SH_TEST_OUTPUT \n";
	print STDERR "************************************ \n";
	print STDERR "\n";
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

sub draw_gene_tree_bipartitions {	# Writes the % of bipartitions found on set of all indiv gene trees onto the concat tree, much like a bootstrap %
	chdir "$path_to_wd/trees";
	system ("cat *bestTree.*besttree > gene_trees_allset");		# Concatenate all indiv gene trees into single tree file
	print STDERR "************************************ \n";
	print STDERR "Writing bipartitions found in individual gene trees to concatenated alignment tree... \n";
	print STDERR "************************************ \n";
	system ("raxmlHPC-PTHREADS -T $NUMTHREADS -m PROTGAMMAWAG -s ../alignments/$concat_aln_prefix\.phy -n bipartitions_write -t RAxML_bestTree.concat_unpart -z gene_trees_allset -f b");
	print STDERR "************************************ \n";
	print STDERR "Bipartition writing done \n";
	print STDERR "************************************ \n";
	print STDERR "\n";
	chdir "$path";
}

sub draw_consensus_tree {		# Draws an MRE consensus tree based on the individual gene trees
	chdir "$path_to_wd/trees";
	print STDERR "************************************ \n";
	print STDERR "Calculating MRE consensus from gene trees...";
	print STDERR "************************************ \n";
	system ("cat *bestTree.*besttree > gene_trees_allset");
	system ("raxmlHPC-PTHREADS -T $NUMTHREADS -m PROTGAMMAWAG -n consMRE -s ../alignments/$concat_aln_prefix\.phy -J MRE -z gene_trees_allset");
	print STDERR "************************************ \n";
	print STDERR "MRE consensus calculated from individual gene trees \n";
	print STDERR "************************************ \n";
	print STDERR "\n";
	chdir "$path";

}
