# SCRIPTS FOR AUTOMATED MULTI-GENE PHYLOGENIES

[![DOI](https://zenodo.org/badge/10602/kbseah/phylogenomics-tools.svg)](https://zenodo.org/badge/latestdoi/10602/kbseah/phylogenomics-tools)

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

Retrieve genomes as Fasta files. Put all scaffolds into single multi-fasta file. Default is for genomic nucleotide sequences, but amino acid sequences can also be supplied (see part 1 below).

For each genome, specify a short name (10 characters or less, to be compatible with RAxML).

Create a plain text file with the original Fasta filename and new shortname for each genome separated by a tab, one line per genome.

Save this table (e.g. as `RENAME_TABLE`)

Make a list of all the marker genes used by AMPHORA2, save file as `MARKER_LIST`

Put all the original fasta files into a folder called `./original_fasta`

## 1. Rename Fasta files with short-names

```bash
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

```bash
 $ perl extract_markers.pl \
         --file RENAME_TABLE \
         --wd PATH_TO_WORKING_DIRECTORY \
         --phyla_amphora [Use Phyla-AMPHORA instead of AMPHORA2] \
         --phylum [Which phylum, for Phyla-AMPHORA? (Default: 3)] \
         --amphora_path PATH_TO_AMPHORA_2 [or Phyla-AMPHORA] \
         --16S [Use Metaxa to extract 16S? ]
```

Uses AMPHORA2 by default to extract conserved markers from each genome. Use --phyla_amphora switch to use Phyla-AMPHORA instead.

Regardless of which version, specify path to the AMPHORA2 or Phyla-AMPHORA directory (e.g. `/home/foobar/tools/AMPHORA2/`).

If using Phyla-AMPHORA, you can specify which phylum's marker set to use with the `--phylum` option (ignored if using AMPHORA2). If using AMPHORA2, you can use the Archaea marker set with the `--archaea` switch.

Calls the script MarkerScanner.pl for each genome given in RENAME_TABLE. By default it assumes the fasta files supplied are genomic nucleotide sequences.

Use the option `--seqtype protein` if you have amino acid sequences instead of genomic nucleotide sequences.

If `--16S` switch is specified, also uses Metaxa to extract 16S sequences from each genome (preset to Bacterial mode - edit script if genomes are Archaea!)

## 3. Make table of how many markers detected per genome

```bash
 $ perl count_output.pl \ 
         --file GENOME_RENAME_TABLE \ 
         --markers MARKER_TABLE \ 
	 --16S [Include 16S in output?] \
         --wd WORKING_FOLDER \
```

The results are in two tab-separated tables: counts_table_protein and counts_table_rna for protein marker genes and 16S markers respectively.

Inspect the tables to choose which markers to use for phylogenetic analysis, and also to remove genomes that are of poor quality. Many missing genes may indicate that genomes are incomplete, or evolutionarily reduced. Many duplicated genes may indicate mixed-strain genome assemblies, or errors introduced during the DNA sequencing process. Incomplete genomes with uneven duplication are common in single-cell MDA products, for example.

IMPORTANT: If there is more than one sequence for a given marker in a genome, you should either
 * Not include that marker for your analysis
 * Not include that genome in your analysis
 * Manually inspect the Fasta file to see which sequence is the "correct" one. For example, short hypothetical proteins or fragmentary ORFs may sometimes be mis-predicted as a marker, or the ORF for a given marker gene may be split into two fragments. In these cases you should manually edit the Fasta file to have only one sequence
 
Draft genomes may contain incomplete 16S genes - check the 16S.graph file for each species to see if all the conserved regions (V1 to V9) are present; choose the sequence that has the most regions covered.

Use the script `counts_table_checker.pl` to make a list of genes that are present in all genomes and only in single copy. Bear in mind that the more genomes you add to your analysis, the fewer will meet this criteria because of missing data, misprediction, etc.

In future, intend to implement new scripts which can deal with missing genes. 

Make a plain text list of marker names to be used for subsequent phylogenetic analysis: ANALYSIS_MARKERS

## 4. Reorganize the Fasta files by marker gene

```bash
 $ perl repackage_fasta.pl \ 
         --file GENOME_RENAME_TABLE \          
         --markers ANALYSIS_MARKERS \ 
	 --16S [Also make 16S alignment?] \
	 --mask [Mask alignment with Zorro?] \
         --wd WORKING_FOLDER 
```

Creates a new folder called "alignments" in working folder.

For each marker gene listed in `ANALYSIS_MARKERS`, creates a multi-Fasta file containing that marker sequence from each species. If a given species is missing a marker, they will be skipped. When concatenating the alignment, they will be replaced by gap characters.

Optional: Also create 16S alignment

Optional: Use Zorro to mask parts of alignment that are more uncertain.

Each sequence header is renamed using the shortname for that species. 

If there is more than one 16S sequence for a given species, an error message will be shown.

For protein genes: Align with MUSCLE each marker file, outputs the alignment file.

For 16S genes: Align with MAFFT LINSI each marker file, outputs the alignment file. Known issue: Metaxa sometimes extracts 16S sequences that are too long (including parts of 23S). You'll be able to see this in the alignment; edit to remove non-16S regions.

At this point, you'll want to have a look at the alignments for each marker, remove poorly-aligned parts, obvious misalignments, etc.

## 5. Create concatenated alignment

```bash
 $ perl concat_align_with_gaps.pl \
         --file GENOME_RENAME_TABLE \
         --markers MARKER_TABLE \
         --wd WORKING_FOLDER \
         --mask [Use mask files produced by Zorro] \
         --model_select [Use automated model selection (slow)] \
         --out PREFIX
```

This produces the following files in the `WORKING_FOLDER/alignments` folder:
 * `PREFIX.concat.fasta` Concatenated alignment in Fasta format
 * `PREFIX.concat.phy` Concatenated alignment in Phylip format
 * `PREFIX.partitions` File listing alignment positions corresponding to each marker (for partitioned phylogenetic analysis, for example)
 * Individual alignment files for each marker in Phylip format, with `PREFIX` in the filename

The `--model` switch implements perl script for finding best protein substitution model, written by Alexis Stamatakis. Otherwise, the partition file will just write "WAG" model for every partition.

Use the --mask switch if you want to include the alignment masking data produced by Zorro in the earlier step.


## 6. Calculate trees

A quick tree calculation can be performed with FastTree, which is useful e.g. to spot long-branching taxa, possible complications, before performing a full analysis with more sophisticated evolutionary models. In cases where alignments are very large, FastTree may be the most feasible option computationally.

```bash
fasttree < PREFIX.concat.phy > PREFIX.concat.fasttree
```


Alternatively, you could use RAxML for phylogenetic tree inference, because it offers more types of models and analyses: 

A. RAxML tree of the entire alignment as a single partition, using PROTCATWAG model for tree search, PROTGAMMAWAG model for optimization, find best tree from 10 randomized starting trees

```bash
raxmlHPC-PTHREADS -T 4 -m PROTCATWAG -s ./alignments/concatenated.aln.phy -n concat_unpart -N 10 -p 12345
```

Strip the best tree result of branch lengths to get only topology -> concat_constraint

B. RAxML best trees for each gene, single partition using PROTCATWAG model for tree search, PROTGAMMAWAG model for optimization

```bash
raxmlHPC-PTHREADS -T 4 -m PROTCATWAG -s ./alignments/$marker.cat.phy -n $marker_besttree -N 10 -p 12345 
```

C. Perform SH-test for each gene's best tree and constraint tree:

```bash
raxmlHPC-PTHREADS -T 4 -m PROTGAMMAWAG -s ./alignments/$marker.cat.phy -n $marker_SH_test -t RAxML_bestTree.$marker_besttree -z concat_constraint -f H
```

These (and more) automated by the script:

```bash
 $ perl tree_calculations.pl \
         --markers MARKER_TABLE \
         --wd WORKING_FOLDER
         --concat_aln CONCATENATED_ALIGNMENT_PREFIX \
         --numthreads NUMBER_OF_RAXML_THREADS (at least 2!) \
         --mask [Use alignment mask?]
```

All the output tree files are in the folder /trees in the working directory.

## Paper(s) using these tools

Petersen et al. 2016. [Nature Microbiology 2: 16195.](http://www.nature.com/articles/nmicrobiol2016195)

Rubin-Blum et al. 2017. [Nature Microbiology 2: 17093.](http://http://www.nature.com/articles/nmicrobiol201793)

Drop me a message if you have used these scripts and would like to add your publication to this list!

## Citing these tools

Please cite as: Brandon Seah (2014) Phylogenomics-tools. Online: https://github.com/kbseah/phylogenomics-tools

DOI available via Zenodo: [![DOI](https://zenodo.org/badge/10602/kbseah/phylogenomics-tools.svg)](https://zenodo.org/badge/latestdoi/10602/kbseah/phylogenomics-tools)
