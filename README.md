# msstitch - MS proteomics post-processing utilities

Shotgun proteomics has a number of bioinformatic tools available for identification 
and quantification of peptides, and the subsequent protein inference. These are a
collection of scripts to integrate a number of these tools, generating
ready to use result files.

If you need support for a specific tool, there is limited time but infinite gratitude :)

## Usage

An example command flow would first store mzML spectra data in an SQLite file:

```
msstitch storespectra --spectra file1.mzML file2.mzML \
  --setnames sampleset1 sampleset2 -o db.sqlite
```

Or, in case the lookup with some other spectra already exists:

```
msstitch storespectra --dbfile lookup.sqlite --spectra file1.mzML file2.mzML \
  --setnames sampleset1 sampleset2
```

Then store quantification data from kronik (MS1 precursor quant) and isobaric 
quantification from OpenMS together with the spectra:

```
msstitch storequant --dbfile db.sqlite --spectra file1.mzML file2.mzML \
  --kronik file1.kronik file2.kronik \
  --isobaric file1.consensusXML file2.consensusXML
```

Create a decoy database where peptides are reversed between tryptic residues:

```
msstitch makedecoy uniprot.fasta -o decoy.fasta --scramble tryp_rev --maxshuffle 10
```

Or without removing peptide sequences that match to the target DB:

```
msstitch makedecoy uniprot.fasta -o decoy.fasta --scramble tryp_rev --ignore-target-hits
```


After running two samples of MSGF and percolator, we can start making 
a more proper set of PSM tables:

```
# Add percolator data, filter 0.01 FDR
msstitch perco2psm -i psms1.txt \
  --perco percolator1.xml --mzid psms1.mzIdentML \
  --filtpsm 0.01 --filtpep 0.01
msstitch perco2psm -i psms2.txt \
  --perco percolator2.xml --mzid psms2.mzIdentML \
  --filtpsm 0.01 --filtpep 0.01
# Combine the two sets and split to a target and decoy file
msstitch concat -i psms1.txt psms2.txt -o allpsms.txt
msstitch split -i allpsms.txt --splitcol TD
```

Now refine the PSM tables, using the earlier created SQLite DB, 
adding more information (sample name, MS1 precursor quant,
isobaric quant, proteingroups, genes):

```
cp db.sqlite decoy_db.sqlite
msstitch psmtable -i target.tsv -o target_psmtable.txt --fasta uniprot.fasta \
  --dbfile db.sqlite --addmiscleav --addbioset --ms1quant --isobaric \
  --proteingroup --genes
msstitch psmtable -i decoy.tsv -o decoy_psmtable.txt --fasta decoy.fasta \
  --dbfile decoy_db.sqlite --proteingroup --genes --addbioset
```

If necessary (e.g. multiple TMT sample sets), split the table before making
protein/peptide tables:

```
msstitch split -i target_psmtable.txt --splitcol bioset
```

Create a peptide table, with summarized median isobaric quant ratios,
highest MS1 intensity PSM as the peptide MS1 quant intensity, and an additional
linear-modeled q-value column:

```
msstitch peptides -i set1_target_psms.txt -o set1_target_peptides.txt \
  --scorecolpattern svm --modelqvals --ms1quant \
  --isobquantcolpattern tmt10plex --denompatterns _126 _127C 
```

Create a protein table, with isobaric quantification as for peptides, the
average of the top-3 highest intensity peptides for MS1 quantification:

```
msstitch proteins -i set1_target_peptides.txt --decoyfn set1_decoy_peptides \
  -o set1_proteins.txt \
  --scorecolpattern '^q-value' --logscore \
  --ms1quant \
  --isobquantcolpattern tmt10plex --denompatterns _126 _127C 
```

Or the analogous process for genes

```
msstitch genes -i set1_target_peptides.txt --decoyfn set1_decoy_peptides \
  -o set1_genes.txt \
  --scorecolpattern '^q-value' --logscore \
  --ms1quant \
  --isobquantcolpattern tmt10plex --denompatterns _126 _127C 
```

Or when there are ENSEMBL entries in the fasta search database, even for ENSG:

```
msstitch ensg -i set1_target_peptides.txt --decoyfn set1_decoy_peptides \
  -o set1_ensg.txt \
  --scorecolpattern '^q-value' --logscore \
  --ms1quant \
  --isobquantcolpattern tmt10plex --denompatterns _126 _127C 
```

Finally, merge multiple sets of proteins (or genes/ENSG) into a single output.
Here we set an cutoff so that features with FDR > 0.01 are set to NA for the 
respective sample set.

```
msstitch merge -i set1_proteins.txt set2_proteins.txt \
  --setnames sampleset1 sampleset2 \
  --dbfile db.sqlite \
  --fdrcolpattern 'q-value' --mergecutoff 0.01 \
   --ms1quantcolpattern area --isobquantcolpattern plex --psmnrcolpattern quanted
```


## Some other useful commands are:
Trypsinize a fasta file (minimum retained peptide length, do cut K/RP, allow 1 missed cleavage)

```
msstitch trypsinize -i uniprot.fasta -o tryp_up.fasta --minlen 7 --cutproline --miscleav 1
```

Create an SQLite file with tryptic sequences for filtering out e.g. known-sequence data.
Options as for trypsinization, --insourcefrag builds lookup with support for 
in-source fragmented peptides that have lost some N-terminal residues:

```
msstitch storeseq -i canonical.fa --cutproline --minlen 7 --miscleav 1 --insourcefrag
```

Filter a percolator output file using the created SQLite, removing sequences
that match those stored in the SQLite. The below also removes sequences in the 
sample which are deamidated (i.e. D -> N), and sequences that have lost at most
2 N-terminal amino acids due to in-source fragmentation (DB must have been 
built with support for that).

```
msstitch filterperco -i perco.xml --dbfile tryptic.sqlite --insourcefrag 2 --deamidate -o filtered.xml
```

Create an SQLite file with full-protein sequences for filtering any peptide of 
a minimum length specified that matches to those. Slower than filtering tryptic 
sequences but more comprehensive:

```
msstitch storeseq -i canonical.fa --fullprotein --minlen 7
```

Filter a percolator output file on protein sequences using the SQLite, removing 
sequences in sample which match to anywhere in the protein. Sequences may be 
deamidated, and minimum length parameter must match the one the database is 
built with.

```
msstitch filterperco -i perco.xml --dbfile proteins.sqlite --fullprotein --deamidate --minlen 7 -o filtered.xml
```


Create an isobaric ratio table median-summarizing the PSMs by any column number 
you want in a PSM table. E.g. you have added a column with exons. The following 
uses average of two channels as denominator, outputs a new table with first column
the features found in column nr.20 of the PSM table:

```
msstitch isosummarize -i psm_table.txt --featcol 20 --isobquantcolpattern tmt10plex --denompatterns 126 127C
```
