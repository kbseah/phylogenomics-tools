#!/usr/bin/perl

# Check marker genes counts table and outputs list of marker genes which occur once and only once per genome
# I.e. where columns are all "1"s 

use strict;
use warnings;
use Getopt::Long;

my $inputfile;
my @headers;
my %tracking_hash;

if ( ! @ARGV ) { usage(); }

GetOptions (
    "table=s" => \$inputfile ) 
or die usage ();

open (my $INPUT, "< $inputfile") or die ("Cannot open counts table: $! \n");
    my $firstline = readline($INPUT);
    chomp $firstline;
    @headers = split (/\t/, $firstline);
    foreach my $theheader (@headers) {      # initialize the hash values
        $tracking_hash{$theheader} = 0;
    }
    while ( ! eof($INPUT) ) {
        my $currentline = readline($INPUT);
        my @splitcurrentline = split(/\t/, $currentline);
        my $counter = 0;
        my $length_of_line = scalar @splitcurrentline;
        while ($counter < $length_of_line) {
            if ($splitcurrentline[$counter] ne "1") {$tracking_hash{$headers[$counter]}++;}
            $counter++;
        }
    }
close ($INPUT);

foreach my $key (keys %tracking_hash) {
    if ($tracking_hash{$key} == 0) { print $key, "\n"; }
}


sub usage {
    print STDERR "********************************************* \n";
    print STDERR "Extract list of single copy genes found in \n";
    print STDERR "all genomes for analysis. \n";
    print STDERR "KBS 2013-11-15";
    print STDERR "Usage: perl counts_table_checker.pl \\ \n";
    print STDERR "\t --table COUNTS_TABLE_FILE \\ \n";
    print STDERR "\t > OUTPUT_LIST \n";
    print STDERR "********************************************** \n";
    exit;
}
