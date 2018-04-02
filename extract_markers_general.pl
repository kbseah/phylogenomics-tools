#!/usr/bin/perl

# Calls AMPHORA2 scripts to call marker genes
# Cite AMPHORA2, HMMER, Metaxa

use strict;
use warnings;
use Cwd qw();
use Getopt::Long;
use Bio::SeqIO;

### Global variables

my $path = Cwd::cwd();
#print "$path\n";

my $rename_table = "rename_table";
my $path_to_working_directory = "new_fasta";
my $custom_hmm;
my @dirnames;    # array to hold genome shortnames
my $NUMTHREADS = 1;
my @marker_list;

### Get options

if (! @ARGV) { usage(); }     # print usage statement if no options specified
GetOptions ("file=s" => \$rename_table,
            "wd=s" => \$path_to_working_directory,
            "hmm:s" => \$custom_hmm,
            ) or usage();

### Main code block

print STDERR "*** Job started *** \n";

read_rename_table();
parse_marker_names();
extract_hmm_markers();


print STDERR "*** Job complete *** \n";


### Subroutines

sub usage {                     # print usage statement
    print STDERR "******************************************************************************************** \n";
    print STDERR "Extract protein sequences from genomes using HMM models \n";
    print STDERR "Cite: HMMer, EMBOSS \n";
    print STDERR "KBS 2014-02-18 \n";
    print STDERR "Usage: perl extract_markers_general.pl \\ \n";
    print STDERR "\t --file RENAME_TABLE \\ \n";
    print STDERR "\t --wd PATH_TO_WORKING_DIRECTORY \\ \n";
    print STDERR "\t --hmm HMM_FILE [Use custom HMM file with marker models] \\ \n";
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

sub parse_marker_names {    # Extract HMM model names from the HMM model file
    open (HMMFILE, "< $custom_hmm") || die ("Cannot open HMM model file: $!\n");
    while (<HMMFILE>) {
        chomp;
        if (/^NAME/) {
        my @split_nameline = split (' ', $_);
        push @marker_list, $split_nameline[1];
        }
        else {next;}
    }
    close (HMMFILE);
    open (MARKEROUT, "> marker_list_out") || die ("Cannot write marker list: $!\n");
        foreach my $markername (@marker_list) {
            print MARKEROUT $markername, "\n";
        }
    close (MARKEROUT);
}

sub extract_hmm_markers {
    print STDERR "Extracting markers... \n";
    foreach my $thegenome (@dirnames) {
        print STDERR "... Working on directory $thegenome ...\n";
        chdir "$path_to_working_directory/$thegenome";
        system ("getorf -sequence $thegenome.fasta -outseq $thegenome.orf -table 11");    # Extract ORFs usign EMBOSS getorf
        system ("hmmsearch --domtblout $thegenome.domtblout -o /dev/null $path\/$custom_hmm $thegenome.orf");    # Runn HMM search on each genome
        parse_hmm_domtblout($thegenome);
        chdir $path;
    }
}

sub parse_hmm_domtblout {
    my $thegenome = $_[0];
    my %hits;
    my %seq;
    open (INPUT, "< $thegenome\.domtblout") || die ("Cannot open hmmsearch output file: $!\n");
        while (<INPUT>) {
            chomp;
            next if /^#/;
            my @theline = split /\s+/;
            my $protein_cds = $theline[0];
            my $the_model = $theline[3];
            $hits{$the_model}{$protein_cds} = 1;    # There should be only one model hit per CDS
        }
    close (INPUT);

        my $seqin = new Bio::SeqIO('-file'=>$thegenome.".orf");
        while (my $seq = $seqin->next_seq) {
                $seq{$seq->id} = $seq;
        }

        for my $marker (keys %hits) {
            my $seqout = new Bio::SeqIO('-file'=>">>$marker.pep",'-format'=>'fasta');
            for my $seqid (keys %{$hits{$marker}}) {
                $seqout->write_seq($seq{$seqid});
            }
        }


}
