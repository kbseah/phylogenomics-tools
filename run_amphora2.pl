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
my $path_to_amphora2;
my @dirnames;    # array to hold genome shortnames

if (! @ARGV) { usage(); }     # print usage statement if no options specified
GetOptions ("file=s" => \$rename_table,
            "wd=s" => \$path_to_working_directory,
            "amphora_path=s" =>\$path_to_amphora2,
            ) or usage();

read_rename_table();
#run_amphora2();
#run_metaxa();
cleanup_metaxa();
print "*** Job complete *** \n";


sub usage {                     # print usage statement
    print "******************************************************* \n";
    print "Call AMPHORA2 scripts to identify markers in genomes \n";
    print "and Metaxa to extract 16S sequences \n";
    print "Cite: AMPHORA2, HMMer, and Metaxa papers \n";
    print "KBS 2014-01-21 \n";
    print "Usage: perl run_amphora2.pl \\ \n";
    print "\t --file RENAME_TABLE \\ \n";
    print "\t --wd PATH_TO_WORKING_DIRECTORY \\ \n";
    print "\t --amphora_path PATH_TO_AMPHORA_2 \n";
    print "******************************************************* \n";
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
    print "Extracting AMPHORA2 markers... \n";
    foreach (@dirnames) {
        print "... Working on directory ", $_, " ... \n";
        chdir "$path_to_working_directory/$_";
        system("perl $path_to_amphora2/Scripts/MarkerScanner.pl -DNA -Bacteria $_\.fasta");
    #    system("perl /home/kbseah/tools/AMPHORA2/Scripts/MarkerAlignTrim.pl -OutputFormat fasta");
        chdir $path;
    }
}

sub run_metaxa {
    print "Extracting 16S with Metaxa ... \n";
    foreach (@dirnames) {
        print "... Working on directory ", $_, " ... \n";
        chdir "$path_to_working_directory/$_";
        system ("metaxa -i $_\.fasta -o $_\_metaxa -t b");
        chdir $path;
    }
}

sub cleanup_metaxa {
    print "Cleaning up Metaxa 16S hits ... \n";
    foreach (@dirnames) {
        print "... Working on directory ", $_, " ... \n";
        chdir "$path_to_working_directory/$_";
        system "mv $_\_metaxa.bacteria.fasta 16S.fasta";
        system "mv $_\_metaxa.graph 16S.graph";
        system "mv $_\_metaxa.extraction.results 16S.extraction.results";
        system "mv $_\_metaxa.summary.txt 16S.summary";
        chdir $path;
    }
}
