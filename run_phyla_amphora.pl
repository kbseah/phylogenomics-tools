#!/usr/bin/perl

# Calls AMPHORA2 scripts to call marker genes
# Cite AMPHORA2, HMMER

use strict;
use warnings;
use Cwd qw();
use Getopt::Long;

my $path = Cwd::cwd();
#print "$path\n";

my $rename_table = "rename_table";
my $path_to_working_directory = "new_fasta";
my $path_to_phyla_amphora = "/home/kbseah/tools/Phyla_AMPHORA/Scripts";
my $PHYLUM = 3;
my @dirnames;  # array to hold genome shortnames
my $NUMTHREADS = 4;

# print usage statement if no options specified
if (! @ARGV) { usage(); }

GetOptions ("file=s" => \$rename_table,
            "wd=s" => \$path_to_working_directory,
            "phyla_path=s" => \$path_to_phyla_amphora,
            "phylum=s" => \$PHYLUM,
            "cpu=s" => \$NUMTHREADS,
            ) or usage();

## MAIN ########################################################################
read_rename_table();
run_phyla_amphora();
print "*** Job complete *** \n";

## SUBS ########################################################################

sub usage {
    # print usage statement
    print "******************************************************* \n";
    print "Call Phyla-AMPHORA to identify markers in genomes \n";
    print "Cite: Phyla-AMPHORA and HMMer papers \n";
    print "NB: Using souped-up version of MarkerScanner.pl called MarkerScanner_souped.pl \n";
    print "KBS 2013-11-14 \n";
    print "Usage: perl run_phyla_amphora.pl \\ \n";
    print "\t --file RENAME_TABLE \\ \n";
    print "\t --wd PATH_TO_WORKING_DIRECTORY \\ \n";
    print "\t --phyla_path PATH_TO_PHYLA_AMPHORA_SCRIPTS \\ \n";
    print "\t --phylum PHYLUM_NUMBER \\ \n";
    print "\t --cpu NUMBER_OF_THREADS \n";
    print "******************************************************* \n";
    exit;
}

sub read_rename_table {
    open (INFO, "< $rename_table") or die("Cannot open renaming table: $!\n");
    while (my $text = <INFO>) {
        chomp $text;
        (my $origfile, my $newfile) = split(/\t/, $text);
        $newfile =~ s/(\w*)\W*/$1/ge;
        push(@dirnames, $newfile);
    }
    close(INFO);
}

sub run_phyla_amphora {
    foreach (@dirnames) {
        print "... Working on directory ", $_, " ... \n";
        chdir "$path_to_working_directory/$_";
        system("perl $path_to_phyla_amphora/MarkerScanner_souped.pl -DNA -Phylum $PHYLUM $_\.fasta --cpu $NUMTHREADS");
    #    system("perl $path_to_phyla_amphora/MarkerAlignTrim.pl -OutputFormat fasta");
        chdir $path;
    }
}
