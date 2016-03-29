#!/usr/bin/perl

=head1 NAME

repackage_fasta.pl - Reorganize Fasta files by marker gene

=head1 SYNOPSIS

    perl repackage_fasta.pl --file <filename table> \
                            --wd <working directory> \
                            --markers <marker_table> \
                            --mask
    
    perl repackage_fasta.pl --help

=head1 DESCRIPTION

Creates a new folder "alignments" in working folder, creates a separate Fasta
file for each marker gene in marker table and aligns them with Muscle.


=head1 ARGUMENTS

=over 8

=item --file <file>

Shortnames and filenames table.

=item --wd <string>

Path to working directory

=item --markers <file>

Table of marker gene names for further processing.

=item --16S

Logical: Include 16S gene? (Default: No)

=item --mask

Logical: Use Zorro to mask alignment columns by alignment quality? (Default: No)

=item --help|-h

This help message

=back

=head1 OUTPUT

New folder "alignments" in working directory. Fasta file for each marker gene
with the sequences from each genome. Alignments for each marker gene created
by Muscle.

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

# Script to read list of marker genes to align for phylogenetic analysis
# Repackage them into Multifasta files with species shortnames as headers
# Call MUSCLE for alignment

use strict;
use warnings;
use Cwd qw();
use File::Path qw(make_path);
use Getopt::Long;
use Pod::Usage;

my $path = Cwd::cwd();
my $species_file;
my $marker_file;
my $path_to_wd;
my @markers;
my @shortnames;
my $use_mask = 0;
my $include_16S =0;

if ( ! @ARGV ) { pod2usage (-message=>"Insufficient input options",-exitstatus=>2); }

GetOptions (
	"file=s" => \$species_file,
	"markers=s" => \$marker_file,
	"wd=s" => \$path_to_wd,
	"mask" => \$use_mask,
	"16S" => \$include_16S,
        "help|h" => sub {pod2usage (-exitstatus=>2,-verbose=>2); }
        ) or pod2usage (-message=>"Please check input options",-exitstatus=>2,-verbose=>2);

read_shortnames();
read_marker_names();
make_path("$path_to_wd/alignments");
call_Muscle();
if ($include_16S == 1) { align_16S(); }
print STDERR "***************************\n";
print STDERR "****** Job complete *******\n";
print STDERR "Alignments in folder $path_to_wd/alignments \n";
print STDERR "***************************\n";


sub usage {                     # print usage statement
        print STDERR "******************************************************* \n";
        print STDERR "Generate alignment for each marker gene using MUSCLE \n";
	print STDERR "Optionally: \n";
	print STDERR "\t Align 16S with MAFFT LINSI \n";
	print STDERR "\t Create alignment mask with Zorro \n";
	print STDERR "Cite: MUSCLE \n";
        print STDERR "KBS 2013-11-11 \n";
        print STDERR "Usage: perl repackage_fasta.pl \\ \n";
        print STDERR "\t --file GENOME_RENAME_TABLE \\ \n";
        print STDERR "\t --markers MARKER_TABLE \\ \n";
	print STDERR "\t --16S [Make 16S alignment?] \\ \n";
	print STDERR "\t --mask [Make alignment mask?] \\ \n";
        print STDERR "\t --wd WORKING_FOLDER \n";
        print STDERR "******************************************************* \n";
        exit;
}

sub read_shortnames {			# Read list of species shortnames
	open (INFO, "< $species_file") || die ("Cannot open list of species names: $!");
		while (<INFO>) {
			chomp;
			(my $throwaway, my $current_shortname) = split(/\t/,$_);
			push (@shortnames, $current_shortname);
		}
	close (INFO);
}

sub read_marker_names {			# Read list of marker genes for extraction 
	open (INFO, "< $marker_file") || die ("Cannot open list of marker genes: $!");
		while (<INFO>) {
			chomp;
			push (@markers, $_);
		}
	close (INFO);
}

sub call_Muscle {					# Perform alignment using MUSCLE and mask with Zorro
	foreach my $marker (@markers) {			# For each marker gene, concatenate all sequences to one multi-fasta file
		system("cat /dev/null > $path_to_wd/alignments/$marker\.cat.pep");
		foreach my $species (@shortnames) {
			open (my $origfile, "< $path_to_wd/$species/$marker\.pep") || die ("Cannot open .pep file to read: $!");
			open (my $catfile, ">> $path_to_wd/alignments/$marker\.cat.pep") || die ("Cannot open file to write multifasta: $!");
				while (<$origfile>) {
					my $newheader;
					if ($_ =~ /\>.*/) {$newheader = $species;}
					s/(\>).*/$1.$newheader/e;
					print $catfile $_;
				}
			close ($catfile);
			close ($origfile);
		}
		print STDERR "***************************\n";
		print STDERR "Created multi-Fasta file for marker ", $marker, "\n";
		print STDERR "Calling MUSCLE to perform alignment for ", $marker, "\n";
							# Call MUSCLE to perform alignment, output to multi-Fasta aln
		system("muscle -in $path_to_wd/alignments/$marker\.cat.pep -out $path_to_wd/alignments/$marker\.cat.aln");
		print STDERR "Alignment complete for marker ", $marker, "\n";
		if ( $use_mask == 1 ) {
								# Call Zorro to mask alignment 
			print STDERR "***************************\n";
			print STDERR "Calling Zorro to mask alignment by confidence for $marker \n";
			system("zorro $path_to_wd/alignments/$marker\.cat\.aln > $path_to_wd/alignments/$marker\.cat\.mask_raw");
			print STDERR "Alignment masking complete for marker $marker \n";
								# Reformat mask to suitable input format for RAxML
			open (MASKRAW, "< $path_to_wd/alignments/$marker\.cat\.mask_raw") || die ("Cannot open mask_raw file for reformatting: $!");
			open (MASKNEW, "> $path_to_wd/alignments/$marker\.cat\.mask") || die ("Cannot open file to write new alignment mask: $!");
				while (<MASKRAW>) {
				        my $rounded = int($_ + 0.5);
					print MASKNEW "$rounded ";
				}
			close (MASKNEW);
			close (MASKRAW);
		} else {}
	}
}

sub align_16S {
	system ("cat /dev/null > $path_to_wd/alignments/16S\.cat\.fasta");
	foreach my $species (@shortnames) {	# create concatenated 16S fasta file
		my $overcount = 0;
		open (my $origfile, "< $path_to_wd/$species/16S\.fasta") || die ("Cannot open .fasta file to read: $!");
		open (my $catfile, ">> $path_to_wd/alignments/16S\.cat\.fasta") || die ("Cannot open file to write multifasta: $!");
			while (<$origfile>) {
				my $newheader;
				if ($_ =~ /\>.*/) {
					$newheader = $species;
					$overcount = $overcount + 1;
				}
				if ($overcount > 1) { 
					print STDERR "There is more than one 16S sequence for species $species \n";
					print STDERR "Exiting... \n";
					exit;
				}
				s/(\>).*/$1.$newheader/e;
				print $catfile $_;
			}
		close ($catfile);
		close ($origfile);
	}
	print STDERR "***************************\n";
	print STDERR "Created multi-Fasta file for 16S sequences \n";
	print STDERR "Calling MAFFT LINSI to align 16S \n";
	system ("linsi $path_to_wd/alignments/16S\.cat\.fasta > $path_to_wd/alignments/16S\.cat\.aln");
	print STDERR "Alignment complete for 16S \n";
	print STDERR "***************************\n";
}
