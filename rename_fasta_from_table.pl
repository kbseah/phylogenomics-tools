#!/usr/bin/perl

=head1 NAME

rename_fasta_from_table.pl - Rename Fasta files by shortnames for phylogenomics-tools

=head1 SYNOPSIS

    perl rename_fasta_from_table.pl --file <filename table> \
                                    --orig <path to originals> \
                                    --wd <new folder name>

    perl rename_fasta_from_table.pl --help

=head1 DESCRIPTION

Rename Fasta files with shortnames (<10 characters) for phylogenomics pipeline.
See https://github.com/kbseah/phylogenomics-tools for more information.

=head1 ARGUMENTS

=over 8

=item --file <file>

Tab-delimited table with two columns: 1. original filenames of Fasta files
for microbial genomes, 2. shortnames (<10 characters alphanumeric and
underscore).

=item --orig <path>

Path to original files

=item --wd <string>

Name for new working directory (will overwrite if already exists)

=item --autoremove

Remove non-word characters automatically from shortnames

=item --strict

Check that all shortnames are 10 characters or fewer in length. This is necessary
if you want the final alignment in Phylip format, which limits the length of
sequence header names. (Default: No)

=item --help|-h

This help message

=back

=head1 OUTPUT

New fasta files with corresponding short names will be created in the new
folder specified by --wd.

Fasta headers will be replaced by FILENAME_COUNTER where COUNTER enumerates
header numbers in a multi-Fasta file

=cut

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
use Pod::Usage;

my %rename_hash;
my $rename_table = "rename_table";
my $path_to_original_files = "original_fasta";
my $path_to_new_files = "new_fasta";
my $strict;
my $autoremove;

if ( ! @ARGV ) { pod2usage (-message=>"Insufficient input options",-exitstatus=>2); }

# specify arguments, otherwise print usage statement
GetOptions ("file=s" => \$rename_table,
            "orig=s" => \$path_to_original_files,
            "strict" => \$strict,
            "autoremove" => \$autoremove,
            "wd=s" => \$path_to_new_files,
            "help|h" => sub {pod2usage (-exitstatus=>2,-verbose=>1); },
            "man" => sub { pod2usage(-exitstatus=>1,-verbose=>2); },
            ) or pod2usage (-message=>"Please check input options",-exitstatus=>2,-verbose=>2);

## MAIN ########################################################################
read_rename_table();
check_for_nonword_shortnames() unless $autoremove;
check_for_duplicate_shortnames();
check_for_shortname_length() if $strict;
rename_fasta_files();
print "*** Job done *** \n";

## SUBS ########################################################################

sub read_rename_table {
    # read rename_table file into rename_hash
    open(INFO, "< $rename_table") || die("Renaming table file not found: $!");
        while (my $text = <INFO>) {
            chomp($text);
            (my $origfile, my $newfile) = split(/\t/, $text);
            $newfile =~ s/(\w*)\W*/$1/ge if $autoremove; # Remove non-word characters
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

sub check_for_nonword_shortnames {
    my @invalid_shortnames;
    foreach my $shortname (values %rename_hash) {
        push @invalid_shortnames, $shortname if $shortname =~ m/\W/;
    }
    if (@invalid_shortnames) {
        print STDERR "The following shortnames in the table $rename_table have non-word characters\n";
        print STDERR "Please use only alphanumeric or underscore in shortnames\n";
        print STDERR join "\n", @invalid_shortnames;
        die ("\n\nExiting...\n");
    }
}

sub check_for_shortname_length {
    my @invalid_shortnames;
    foreach my $shortname (values %rename_hash) {
        push @invalid_shortnames, $shortname if length($shortname) > 10;
    }
    if (@invalid_shortnames) {
        print STDERR "The following shortnames in the table $rename_table are longer than 10 chars\n";
        print STDERR "You have chosen the --strict option. Please edit those names to be shorter than 10\n";
        print STDERR join "\n", @invalid_shortnames;
        die ("\n\nExiting ...\n");
    }
}

sub rename_fasta_files {
    # for each fasta file, create new directory, fasta file, rename fasta headers by shortname
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