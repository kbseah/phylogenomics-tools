#!/usr/bin/perl

=head1 NAME

count_output.pl - Tabulate how many markers detected per genome

=head1 SYNOPSIS

    perl count_output.pl --file <rename table> \
			 --wd <working directory> \
                         --markers <marker_table> \
    
    perl count_output.pl --help

=head1 DESCRIPTION

Tabulate how many marker genes detected per genome.

=head1 ARGUMENTS

=over 8

=item --file <file>

Shortnames and filenames table.

=item --wd <string>

Path to working directory

=item --markers <file>

Table of marker gene names for further processing.

=item --help|-h

This help message

=back

=head1 OUTPUT

Output files will be written to working directory supplied.

=over 8

=item counts_table_protein

Tab-separated table of genomes (rows) and marker genes (columns), giving
number of counts for each marker gene per genome.

=item counts_table_rna

Only reported if --16S switch is used.

=back

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
use Pod::Usage;

# Script to count number of marker genes per species
# After AMPHORA2 or Phyla-AMPHORA processing
# Given list of species names (same as used in previous scripts) 
# And list of marker genes (corresponding to file names of the marker .pep files)

# KBS 2013-11-11

my $path = Cwd::cwd();

#print "current path ",$path,"\n";

my %counts_by_sp;
my @spnames;
my @markers;
my $rename_table_file;
my $marker_table_file;
my $path_to_new_files;
my $include_16S = 0;

if ( ! @ARGV ) { pod2usage (-message=>"Insufficient input options",-exitstatus=>2); }

GetOptions (
        "file=s" => \$rename_table_file,
        "markers=s" => \$marker_table_file,
	"16S" => \$include_16S,
        "wd=s" => \$path_to_new_files,
        "help|h" => sub {pod2usage (-exitstatus=>2,-verbose=>2); }
    ) or pod2usage (-message=>"Please check input options",-exitstatus=>2,-verbose=>2);
# specify arguments, otherwise print usage statement

read_genome_shortnames();
read_marker_names();
if ( $include_16S == 1 ) { add_16S(); }
generate_output();
print STDERR "*** Job complete *** \n";

sub usage {                     # print usage statement
        print STDERR "******************************************************* \n";
        print STDERR "Count number of marker genes per species \n";
	print STDERR "Default outputs in working folder: \n";
	print STDERR "\t counts_table_protein and counts_table_rna \n";
        print STDERR "KBS 2013-11-11 \n";
        print STDERR "Usage: perl count_output.pl \\ \n";
        print STDERR "\t --file GENOME_RENAME_TABLE \\ \n";
        print STDERR "\t --markers MARKER_TABLE \\ \n";
	print STDERR "\t --16S \[include 16S in table?\] \\ \n";
        print STDERR "\t --wd WORKING_FOLDER \n";
        print STDERR "******************************************************* \n";
        exit;
}


sub read_genome_shortnames {
	open (INFO, "< $rename_table_file") || die("Cannot open file with renaming table: $!\n");
		while (<INFO>) {
			chomp;
			(my $throwaway, my $spname) = split(/\t/,$_);
			$spname =~ s/(\w*)\W*/$1/ge;
			push (@spnames, $spname);
		}
	close (INFO);
}

sub read_marker_names {
	open (TONFO, "< $marker_table_file") || die("Cannot open file with list of markers: $!\n");
		while (<TONFO>) {
			chomp;
			push (@markers, $_);
		}
	close(TONFO);
}

sub add_16S {
	my $theline = "16S";
	push (@markers, $theline);
}

sub generate_output {		# generates output to STDOUT
	my @sorted_markers = sort @markers;
	open (PEPFILE, "> $path_to_new_files/counts_table_protein") || die ("Cannot open file to write counts of protein marker genes: $! \n");
	if ($include_16S == 1) {open (RNAFILE, "> $path_to_new_files/counts_table_rna") || die ("Cannot open file to write counts of 16S rRNA marker genes: $! \n"); }
	print PEPFILE join ("\t","shortname", @sorted_markers), "\n";	# print header for the table to protein marker counts table, tab-separated
	if ($include_16S == 1) {print RNAFILE "shortname","\t","16S","\n";}
	foreach my $spname (@spnames){
		print PEPFILE $spname ,"\t";	# print species name
		my %counts_by_marker;
		chdir "$path_to_new_files/$spname";
			foreach my $marker(@markers) {
				$counts_by_marker{"$marker"} = 0;	# set associative array counter to zero
				open (MANFO, "< $marker\.pep") ;
					my @contents = <MANFO>;
					my @filtered = grep (/\>/,@contents);
					my $thecount = scalar @filtered;
					$counts_by_marker{"$marker"} = $thecount;
	#				while (<MANFO>) {
	#					if (/>/) {$counts_by_marker{"$marker"}++;}	# counts Fasta headers in marker .pep file to count instances of the gene
	#				}
				close(MANFO);
			}
		my @counts_values;
		foreach my $marker (sort keys %counts_by_marker) {
			push (@counts_values, $counts_by_marker{$marker});
		}
		print PEPFILE join ("\t", @counts_values), "\n";	# print counts for each marker gene, separated by tab
		if ( $include_16S == 1 ) {	# Do the same for 16S markers extracted by Metaxa if the -16S switch is on
			print RNAFILE $spname ,"\t";
			my $thecount = 0;
			open (CANFO, "< 16S\.fasta");
				my @contents = <CANFO>;
				my @filtered = grep (/\>/, @contents);
				$thecount = scalar @filtered;
			close (CANFO);
			print RNAFILE $thecount, "\n";
		}
		chdir $path;
	}
	close (PEPFILE);
	if ( $include_16S == 1) { close (RNAFILE); }
}
