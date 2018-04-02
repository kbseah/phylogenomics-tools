#!/usr/bin/perl

# Check marker genes counts table and outputs list of marker genes which occur once and only once per genome
# I.e. where columns are all "1"s

=head1 NAME

counts_table_checker.pl - List single-copy marker genes extracted

=head1 SYNOPSIS

perl counts_table_checker.pl --table <counts_table_file> --out <shortlist>

perl counts_table_checker.pl --help

=head1 DESCRIPTION

Output shortlist of single-copy marker genes identified by the phylogenetics-tools
pipeline, from the table of marker gene counts produced by I<count_output.pl>.
That table lists how many of each marker gene was detected per genome. Those
occuring in multiple copies per genome are likely paralogs and have to be
manually examined to determine which copy is suitable for further analysis.

For phylogenetics, it is preferable to use single-copy genes only. If some genomes
are missing some genes that are otherwise single-copy in other genomes (e.g. if
some genomes in the analysis are incomplete), the --zeros option allows those
to be included in the shortlist. Otherwise, the shortlist contains only genes
that occur once and only once in all genomes.

=cut

use strict;
use warnings;
#use diagnostics;
use Getopt::Long;
use Pod::Usage;

my $inputfile;
my $outfile = "marker_list_singlecopy";
my $zeros;
my @headers;
my %tracking_hash;

if ( ! @ARGV ) { usage(); }

GetOptions ("table=s" => \$inputfile,
            "out=s"=> \$outfile,
            "zeros" => \$zeros,
            ) or die usage ();

=head1 ARGUMENTS

=over 8

=item --table <file>

Table of counts per marker gene per genome, produced by I<count_output.pl>

=item --zeros

Logical: Include marker genes that are missing in some genomes? (Default: No)

If this option is chosen, the missing genes will be represented by gaps in the
concatenated alignment. Use I<concat_align_with_gaps.pl> to concatenate
alignments where some species have missing genes.

=item --out <file>

Name for output file. (Default: "marker_list_singlecopy")

=back

=cut

open (my $INPUT, "< $inputfile") or die ("Cannot open counts table: $! \n");
    my $firstline = <$INPUT>;
    chomp $firstline;
    @headers = split (/\t/, $firstline);
    #shift @headers; # Remove the first column header which is not marker name
    foreach my $theheader (@headers) {
        # initialize the hash values
        $tracking_hash{$theheader} = 0;
    }
    while (my $currentline = <$INPUT> ) {
        chomp $currentline;
        my @splitcurrentline = split(/\t/, $currentline);
        #my $counter = 0;
        #my $length_of_line = scalar @splitcurrentline;
        for (my $counter=1; $counter <= $#headers; $counter++) {
            #print $splitcurrentline[0]."\t".$headers[$counter]."\t".$splitcurrentline[$counter]."\n";
            if (defined $splitcurrentline[$counter]) {
                if ($splitcurrentline[$counter] > 1) {
                    $tracking_hash{$headers[$counter]}++;
                } elsif ($splitcurrentline[$counter] == 0) {
                    $tracking_hash{$headers[$counter]}++ unless defined $zeros;
                }
            }

        }
    }
close ($INPUT);

open (my $fhout, ">", $outfile) or die ("Cannot write output file: $!");
foreach my $key (keys %tracking_hash) {
    if ($tracking_hash{$key} == 0) {
        print $fhout $key, "\n" unless $key eq 'shortname';
    }
}
close ($fhout);


sub usage {
    print STDERR "********************************************* \n";
    print STDERR "Extract list of single copy genes found in \n";
    print STDERR "all genomes for analysis. \n";
    print STDERR "KBS 2013-11-15";
    print STDERR "Usage: perl counts_table_checker.pl \\ \n";
    print STDERR "\t --table COUNTS_TABLE_FILE \\ \n";
    print STDERR "\t --zeros \\ # Allow missing genes \n";
    print STDERR "\t > OUTPUT_LIST \n";
    print STDERR "********************************************** \n";
    exit;
}

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
