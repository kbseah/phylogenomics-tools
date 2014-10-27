#!/usr/bin/perl

# This script renames fasta files for further processing
# Input:
## Tab-delimited table with column 1 - original fasta filenames, column 2 - replacement filenames
## Fasta files with the above filenames
# Output:
## New fasta files with replacement filenames
## Fasta headers replaced by FILENAME_COUNTER where COUNTER enumerates header numbers in a multi-Fasta file

# KBS 2013-11-08

use strict;
use warnings;
use Getopt::Long;
use File::Path qw(make_path);

my %rename_hash;
my $rename_table = "rename_table";
my $path_to_original_files = "original_fasta";
my $path_to_new_files = "new_fasta";

if ( ! @ARGV ) {usage();}

GetOptions (
	"file=s" => \$rename_table,
	"orig=s" => \$path_to_original_files,
	"wd=s" => \$path_to_new_files,) or usage();	# specify arguments, otherwise print usage statement
read_rename_table();
check_for_duplicate_shortnames();
rename_fasta_files();
print "*** Job done *** \n";

sub usage {			# print usage statement
	print "******************************************************* \n";
	print "Rename fasta files by shortnames for further processing \n";
	print "KBS 2013-11-08 \n";
	print "Usage: perl rename_fasta_from_table.pl \\ \n";
	print "\t --file RENAME_TABLE \\ \n";
	print "\t --orig PATH_TO_ORIGINAL_FILES \\ \n";
	print "\t --wd NEW_WORKING_FOLDER \n";
	print "******************************************************* \n";
	exit;
}

sub read_rename_table {		# read rename_table file into rename_hash
	open(INFO, "< $rename_table") || die("Renaming table file not found: $!");
		while (my $text = <INFO>) {
			chomp($text);
			(my $origfile, my $newfile) = split(/\t/, $text);
			$newfile =~ s/(\w*)\W*/$1/ge; 
			$rename_hash{"$origfile"} = "$newfile";
		}
	close(INFO);
}

sub check_for_duplicate_shortnames {
	my %shortname_counts;
	my @shortnames = values %rename_hash;
	for (@shortnames) {
		$shortname_counts{$_}++;
	}
	my @duplicates = grep { $shortname_counts{$_} > 1} keys %shortname_counts;
	if (scalar @duplicates > 0) {
		print "*** ERROR *** \n";
		print "Duplicate shortnames found: ", join(", ", @duplicates), "\n";
		print "Exiting... \n";
		print "************* \n";
		exit;
	}
}

sub rename_fasta_files {	# for each fasta file, create new directory, fasta file, rename fasta headers by shortname
	foreach my $key (keys %rename_hash) {
		make_path("$path_to_new_files/$rename_hash{$key}");
		print "Original file: $key, New file: $rename_hash{$key}\n"; 
		open(my $editfile, "< $path_to_original_files/$key") || die("Cannot open input file: $!");
		open(my $outputfile, "> $path_to_new_files/$rename_hash{$key}/$rename_hash{$key}\.fasta") || die ("Cannot write new file: $!");
		my $counter;
		$counter=0;
		while (<$editfile>) {
			my $newname;
			if ($_ =~ /\>.*/) {
				$counter++; 
				$newname = join("_",$rename_hash{$key},$counter);
			}
			s/(\>).*/$1.$newname/e;
			print $outputfile $_;
		}
		close($editfile);
		close($outputfile);
	}
}


