# SCRIPTS FOR AUTOMATED MULTI-GENE PHYLOGENIES

Have a bunch of microbial genomes that you want to make a multi-gene tree with? 

Daunted by the prospect of having to extract, align, and tree conserved marker genes?

Give these scripts a try!

Starting with genomic sequences in Fasta files, these scripts automate the process of:
 * Extracting marker genes using AMPHORA2 or Phyla-AMPHORA scripts
 * Extracting 16S genes using Metaxa
 * Aligning each gene (using MUSCLE for proteins, MAFFT LINSI for 16S) 
 * Generating alignment mask (weighted columns) using Zorro
 * Concatenating alignments into a single partitioned alignment
 * Reformatting input files for phylogeny reconstruction (with RAxML)
 * Choosing the best substitution model for protein coding genes (with RAxML)
 * Running different kinds of tree analysis (with RAxML)

## Software required

 * AMPHORA2 (https://github.com/martinwu/AMPHORA2) or Phyla-AMPHORA (http://wolbachia.biology.virginia.edu/WuLab/Software.html) to extract marker genes
 * Metaxa (http://microbiology.se/software/metaxa/) to extract 16S sequences (optional)
 * Perl and BioPerl
 * Muscle (http://www.drive5.com/muscle/) for sequence alignment
 * MAFFT (http://mafft.cbrc.jp/alignment/software/) for rRNA sequence alignment (optional)
 * Zorro (http://sourceforge.net/projects/probmask/) for alignment masking (optional)
 * RAxML (http://sco.h-its.org/exelixis/web/software/raxml/index.html) for tree inference (or you can use your favorite phylogenetic program)

## 0. Before you start

Retrieve genomes as Fasta files. Put all scaffolds into single multi-fasta file.

For each genome, specify a short name (10 characters or less, to be compatible with RAxML).

Create a plain text file with the original Fasta filename and new shortname for each genome separated by a tab, one line per genome.

Save this table (e.g. as `RENAME_TABLE`)

Make a list of all the marker genes used by AMPHORA2, save file as `MARKER_LIST`

Put all the original fasta files into a folder called `./original_fasta`

## 1. Rename Fasta files with short-names

```
 $ perl rename_fasta_from_table.pl \ 
         --file RENAME_TABLE \ 
         --orig PATH_TO_ORIGINAL_FILES \ 
         --wd NEW_WORKING_FOLDER 
```

Reads `RENAME_TABLE`.

Reads original Fasta files from `PATH_TO_ORIGINAL_FASTA`

Creates directory `NEW_WORKING_FOLDER` (silent overwrite)

Then generates new directory per genome, renames fasta file and directory according to shortname, as well as Fasta headers

## 2. Extract protein-coding marker genes and 16S sequences

```
 $ perl extract_markers.pl \
         --file RENAME_TABLE \
         --wd PATH_TO_WORKING_DIRECTORY \
         --phyla_amphora [Use Phyla-AMPHORA instead of AMPHORA2] \
         --phylum [Which phylum, for Phyla-AMPHORA? (Default: 3)] \
         --amphora_path PATH_TO_AMPHORA_2 [or Phyla-AMPHORA] \
         --16S [Use Metaxa to extract 16S? ]
```

Uses AMPHORA2 by default to extract conserved markers from each genome. Use --phyla_amphora switch to use Phyla-AMPHORA instead.

Regardless of which version, specify path to the AMPHORA2 or Phyla-AMPHORA directory (e.g. `/home/foobar/tools/AMPHORA2/`)

Calls the script MarkerScanner.pl for each genome given in RENAME_TABLE

If --16S switch is specified, also uses Metaxa to extract 16S sequences from each genome (preset to Bacterial mode - edit script if genomes are Archaea!)

## 3. Make table of how many markers detected per genome

```
 $ perl count_output.pl \ 
         --file GENOME_RENAME_TABLE \ 
         --markers MARKER_TABLE \ 
	 --16S [Include 16S in output?] \
         --wd WORKING_FOLDER \
```

The results are in two tab-separated tables: counts_table_protein and counts_table_rna for protein marker genes and 16S markers respectively.

Inspect the tables to choose which markers to use for phylogenetic analysis.

IMPORTANT: For 16S, if there is more than 1 16S sequence in a given genome, you will have to manually edit the 16S.fasta file to have only 1 paralog!

Draft genomes may contain incomplete 16S genes - check the 16S.graph file for each species to see if all the conserved regions (V1 to V9) are present; choose the sequence that has the most regions covered.

For the moment, select only columns which are all '1'. Can use the script counts_table_checker.pl to make a list of single-copy genes

In future, intend to implement new scripts which can deal with missing genes. 

Make a plain text list of marker names to be used for subsequent phylogenetic analysis: ANALYSIS_MARKERS

## 4. Reorganize the Fasta files by marker gene

```
 $ perl repackage_fasta.pl \ 
         --file GENOME_RENAME_TABLE \          
         --markers MARKER_TABLE \ 
	 --16S [Also make 16S alignment?] \
	 --mask [Mask alignment with Zorro?] \
         --wd WORKING_FOLDER 
```

Creates a new folder called "alignments" in working folder.

For each marker gene listed in `ANALYSIS_MARKERS`, creates a multi-Fasta file containing that marker sequence from each species. 

Optional: Also create 16S alignment

Optional: Use Zorro to mask parts of alignment that are more uncertain.

Each sequence header is renamed using the shortname for that species. 

If there is more than one 16S sequence for a given species, an error message will be shown.

For protein genes: Align with MUSCLE each marker file, outputs the alignment file.

For 16S genes: Align with MAFFT LINSI each marker file, outputs the alignment file. Known issue: Metaxa sometimes extracts 16S sequences that are too long (including parts of 23S). You'll be able to see this in the alignment; edit to remove non-16S regions.

At this point, you'll want to have a look at the alignments for each marker, remove poorly-aligned parts, obvious misalignments, etc.

## 5. Create concatenated alignment

```
 $ perl concat_align.pl \
         --file GENOME_RENAME_TABLE \
         --markers MARKER_TABLE \
         --wd WORKING_FOLDER \
         --mask [Use mask files produced by Zorro] \
         --model_select [Use automated model selection (slow)] \
         --out CONCATENATED_ALN_PREFIX
```

The `--model` switch implements perl script for finding best protein substitution model, written by Alexis Stamatakis. 

If you use this option, there will be a large number of working files generated by RAxML. Discard them with the commands:

```
 $ rm *EVAL*
 $ rm RAxML*
 $ rm ST_*
```

Otherwise, the partition file will just write "WAG" model for every partition.

Use the --mask switch if you want to include the alignment masking data produced by Zorro in the earlier step.

`CONCATENATED_ALN_PREFIX` is the file name prefix for the concatenated alignment and associated files.

## 6. Calculate trees

A. RAxML tree of the entire alignment as a single partition, using PROTCATWAG model for tree search, PROTGAMMAWAG model for optimization, find best tree from 10 randomized starting trees

```
raxmlHPC-PTHREADS -T 4 -m PROTCATWAG -s ./alignments/concatenated.aln.phy -n concat_unpart -N 10 -p 12345
```

Strip the best tree result of branch lengths to get only topology -> concat_constraint

B. RAxML best trees for each gene, single partition using PROTCATWAG model for tree search, PROTGAMMAWAG model for optimization

```
raxmlHPC-PTHREADS -T 4 -m PROTCATWAG -s ./alignments/$marker.cat.phy -n $marker_besttree -N 10 -p 12345 
```

C. Perform SH-test for each gene's best tree and constraint tree:

```
raxmlHPC-PTHREADS -T 4 -m PROTGAMMAWAG -s ./alignments/$marker.cat.phy -n $marker_SH_test -t RAxML_bestTree.$marker_besttree -z concat_constraint -f H
```

These (and more) automated by the script:

```
 $ perl tree_calculations.pl \
         --markers MARKER_TABLE \
         --wd WORKING_FOLDER
         --concat_aln CONCATENATED_ALIGNMENT_PREFIX \
         --numthreads NUMBER_OF_RAXML_THREADS (at least 2!) \
         --mask [Use alignment mask?]
```

All the output tree files are in the folder /trees in the working directory.

## Citing these tools

Please cite as: Brandon Seah (2014) Phylogenomics-tools. Online: https://github.com/kbseah/phylogenomics-tools
