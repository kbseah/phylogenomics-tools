#!/usr/bin/perl

# Calls AMPHORA2 scripts to call marker genes
# Cite AMPHORA2, HMMER, Metaxa

use strict;
use warnings;
use Cwd qw();
use Getopt::Long;

### Global variables

my $path = Cwd::cwd();
#print "$path\n";

my $rename_table = "rename_table";
my $path_to_working_directory = "new_fasta";
my $path_to_amphora;
my $use_metaxa = 0;
my $use_phyla = 0;
my $PHYLUM = 3;
my @dirnames;	# array to hold genome shortnames
my $NUMTHREADS = 1;

### Get options

if (! @ARGV) { usage(); } 	# print usage statement if no options specified
GetOptions (
	"file=s" => \$rename_table,
	"wd=s" => \$path_to_working_directory,
	"amphora_path=s" =>\$path_to_amphora, 
	"phyla_amphora" => \$use_phyla,
	"cpu:s" => \$NUMTHREADS,
	"phylum:s" => \$PHYLUM,
	"16S" => \$use_metaxa )
or usage();

### Main code block 

print STDERR "*** Job started *** \n";

read_rename_table();
if ( $use_phyla == 0 ) {
	run_amphora2();
}
elsif ( $use_phyla == 1 ) {
	run_phyla_amphora();
}
if ( $use_metaxa == 1 ) {
	run_metaxa();
	cleanup_metaxa();
}
print STDERR "*** Job complete *** \n";


### Subroutines

sub usage {                     # print usage statement
        print STDERR "******************************************************************************************** \n";
        print STDERR "Call AMPHORA2 or Phyla-AMPHORA scripts to identify markers in genomes \n";
	print STDERR "and Metaxa to extract 16S sequences \n";
	print STDERR "Cite: AMPHORA2, HMMer, and Metaxa papers \n";
        print STDERR "KBS 2014-01-21 \n";
        print STDERR "Usage: perl extract_markers.pl \\ \n";
        print STDERR "\t --file RENAME_TABLE \\ \n";
        print STDERR "\t --wd PATH_TO_WORKING_DIRECTORY \\ \n";
	print STDERR "\t --phyla_amphora [Use Phyla-AMPHORA instead of AMPHORA2] \\ \n";
	print STDERR "\t --phylum [Which phylum, for Phyla-AMPHORA? (Default: 3)] \\ \n";
	print STDERR "\t --amphora_path PATH_TO_AMPHORA_2 [or Phyla-AMPHORA] \\ \n";
	print STDERR "\t --cpu [Use parallelized Phyla-Amphora script with this no. of threads, slightly faster] \\ \n";
	print STDERR "\t --16S [Use Metaxa to extract 16S? ] \n";
        print STDERR "******************************************************************************************** \n";
        exit;
}

sub read_rename_table {
	open (INFO, "< $rename_table") || die("Cannot open renaming table: $!\n");
		while (my $text = <INFO>) {
			chomp $text;
			(my $origfile, my $newfile) = split(/\t/, $text);
			$newfile =~ s/(\w*)\W*/$1/ge;
			push(@dirnames, $newfile);
		}
	close(INFO);
}

sub run_amphora2 {
	print STDERR "Extracting AMPHORA2 markers... \n";
	foreach (@dirnames) {
		print STDERR "... Working on directory ", $_, " ... \n";
		chdir "$path_to_working_directory/$_";
		system("perl $path_to_amphora/Scripts/MarkerScanner.pl -DNA -Bacteria $_\.fasta");
	#	system("perl $path_to_amphora/Scripts/MarkerAlignTrim.pl -OutputFormat fasta");
		chdir $path;
	}
}

sub run_phyla_amphora {
	print STDERR "Extracting Phyla-AMPHORA markers... \n";
	foreach (@dirnames) {
		print STDERR "... Working on directory ", $_, " ... \n";
		chdir "$path_to_working_directory/$_";
		if ( $NUMTHREADS > 1 ) {
		system("perl $path_to_amphora/Scripts/MarkerScanner_souped.pl -DNA -Phylum $PHYLUM $_\.fasta --cpu $NUMTHREADS"); }
		elsif ( $NUMTHREADS == 1 ) {
		system("perl $path_to_amphora/Scripts/MarkerScanner.pl -DNA -Phylum $PHYLUM $_\.fasta"); }
#		system("perl $path_to_amphora/Scripts/MarkerAlignTrim.pl -OutputFormat fasta");
		chdir $path;
	}
}

sub run_metaxa {
	print STDERR "Extracting 16S with Metaxa ... \n";
	foreach (@dirnames) {
		print STDERR "... Working on directory ", $_, " ... \n";
		chdir "$path_to_working_directory/$_";
		system ("metaxa -i $_\.fasta -o $_\_metaxa -t b");
		chdir $path;
	}
}

sub cleanup_metaxa {
	print STDERR "Cleaning up Metaxa 16S hits ... \n";
	foreach (@dirnames) {
		print STDERR "... Working on directory ", $_, " ... \n";
		chdir "$path_to_working_directory/$_";
		system "mv $_\_metaxa.bacteria.fasta 16S.fasta";
		system "mv $_\_metaxa.graph 16S.graph";
		system "mv $_\_metaxa.extraction.results 16S.extraction.results";
		system "mv $_\_metaxa.summary.txt 16S.summary";
		chdir $path;
	}
}
