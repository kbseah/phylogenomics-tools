#!/usr/bin/perl

=head1 NAME

perl extract_markers.pl - Extract protein-coding marker genes and 16S rRNA genes

=head1 SYNOPSIS

    perl extract_markers.pl --file <rename_table> \
                            --wd <working directory> \
                            --amphora_path <path to Amphora2>

    perl extract_markers.pl --help

=head1 DESCRIPTION

Run Amphora2 or Phyla-Amphora scripts to extract protein-coding marker genes
for each genome. Optionally extract 16S rRNA gene with Metaxa.

=head1 ARGUMENTS

=over 8

=item --file <file>

Shortnames and filenames table.

=item --wd <path>

Path to working directory

=item --seqtype <string>

Type of input sequence; either "DNA" or "protein" (Default: DNA)

=item --phyla_amphora

Use Phyla-Amphora instead of Amphora2 (Default: No)

=item --phylum <integer>

Which phylum-specific set of marker genes to use, for Phyla-Amphora (Default: 3)

=item --amphora_path <path>

Path to Amphora2 or Phyla-Amphora scripts and databases.

=item --archaea

Use Archaeal markers for Amphora2 (Default: no).
Does not work with Phyla-Amphora (which has only Bacterial markers)

=item --cpu <integer>

(For testing purposes only - doesn't work with standard Phyla-Amphora scripts)
Number of CPUs for modified Phyla-Amphora scripts

=item --16S

Logical: Extract 16S rRNA sequences too, with Metaxa?

=item --help|-h

This help message

=back

=head1 OUTPUT

For each genome specified in the list given in --file, a folder will be created
with the corresponding shortname, and the Amphora2 or Phyla-Amphora markers
extracted as Fasta files there.

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

# Calls AMPHORA2 scripts to call marker genes
# Cite AMPHORA2, HMMER, Metaxa

use strict;
use warnings;
use Cwd qw();
use Getopt::Long;
use Pod::Usage;

### Global variables

my $path = Cwd::cwd();
#print "$path\n";

my $rename_table = "rename_table";
my $path_to_working_directory = "new_fasta";
my $path_to_amphora;
my $use_metaxa = 0;
my $use_phyla = 0;
my $seqtype = "DNA";
my $use_Archaea = 0;
my $PHYLUM = 3;
my @dirnames;    # array to hold genome shortnames
my $NUMTHREADS = 1;

### Get options

if ( ! @ARGV ) { pod2usage (-message=>"Insufficient input options",-exitstatus=>2); }
GetOptions (
    "file=s" => \$rename_table,
    "wd=s" => \$path_to_working_directory,
    "seqtype=s" => \$seqtype,
    "amphora_path=s" =>\$path_to_amphora,
    "phyla_amphora" => \$use_phyla,
    "cpu:s" => \$NUMTHREADS,
    "phylum:s" => \$PHYLUM,
    "archaea" => \$use_Archaea,
    "16S" => \$use_metaxa,
    "help|h" => sub {pod2usage (-exitstatus=>2,-verbose=>2); }
    ) or pod2usage (-message=>"Please check input options",-exitstatus=>2,-verbose=>2);
# specify arguments, otherwise print usage statement


### Main code block

print STDERR "*** Job started *** \n";

read_rename_table();

if ( $use_phyla == 0 ) {
    run_amphora2();
}
elsif ( $use_phyla == 1 ) {
    if ($use_Archaea == 1) { # Catch invalid combination
        print STDERR "Error: Cannot use Archaeal markers with Phyla-Amphora \n";
        print STDERR "Exiting ...\n";
        exit;
    } else {
        if ($PHYLUM < 0 || $PHYLUM > 20) {
            print STDERR "Warning: Phylum number must be between 0 and 20 (see Phyla-Amphora manual)\n";
            print STDERR "Resetting Phylum parameter to 0 (use all markers)\n";
            $PHYLUM = 0;
        }
        run_phyla_amphora();
    }
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
            push(@dirnames, $newfile) unless ($text =~ m/^#/); # Ignore lines that are commented-out
        }
    close(INFO);
}

sub run_amphora2 {
    print STDERR "Extracting AMPHORA2 markers... \n";
    my @options; # Amphora2 options
    if ($seqtype eq "DNA") {
        push @options, "-DNA";
    }
    if ($use_Archaea == 0 ) {
        push @options, "-Bacteria";
    } elsif ($use_Archaea == 1) {
        push @options, "-Archaea";
    }
    my $options_string = join " ", @options;
    foreach (@dirnames) {
        my $curr_dir = $_;
        print STDERR "... Working on directory ", $curr_dir, " ... \n";
        chdir "$path_to_working_directory/$curr_dir";
        my $command_string = "perl $path_to_amphora/Scripts/MarkerScanner.pl ".$options_string." $curr_dir\.fasta";
        system($command_string);
        #system("perl $path_to_amphora/Scripts/MarkerScanner.pl -DNA -Bacteria $_\.fasta");
        #system("perl $path_to_amphora/Scripts/MarkerAlignTrim.pl -OutputFormat fasta");
        chdir $path;
    }
}

sub run_phyla_amphora {
    print STDERR "Extracting Phyla-AMPHORA markers... \n";
    my @options;
    if ($seqtype eq "DNA") {
        push @options, "-DNA";
    }
    push @options, "-Phylum $PHYLUM";
    my $options_string = join " ", @options;
    foreach (@dirnames) {
        my $curr_dir = $_;
        print STDERR "... Working on directory ", $curr_dir, " ... \n";
        chdir "$path_to_working_directory/$curr_dir";
        if ( $NUMTHREADS > 1 ) {
            my $command_string = "perl $path_to_amphora/Scripts/MarkerScanner_souped.pl ".$options_string." $curr_dir\.fasta --cpu $NUMTHREADS";
            system($command_string);
        } elsif ( $NUMTHREADS == 1 ) {
            my $command_string = "perl $path_to_amphora/Scripts/MarkerScanner.pl ".$options_string." $curr_dir\.fasta";
            system($command_string);
            #system("perl $path_to_amphora/Scripts/MarkerAlignTrim.pl -OutputFormat fasta");
        }
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
