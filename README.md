# msstitch - MS proteomics post-processing utilities

Shotgun proteomics has a number of bioinformatic tools available for identification 
and quantification of peptides, and the subsequent protein inference. `msstitch` is a 
tool to integrate a number of these tools, generating ready to use result files. Supported
formats are currently:

- mzML files
- Fasta peptides
- MSGF+ PSM tables (and mzIdentML)
- Sage PSM tables (not for quant yet)
- Percolator output
- OpenMS consensusXML isobaric quant
- Hardklor/Kronik MS1 quant output
- Dinosaur MS1 quant output

If you need support for another format, there is limited time but infinite gratitude :)

## Usage

### Storing data
An example command flow would first store mzML spectra data in an SQLite file:

```
msstitch storespectra --spectra file1.mzML file2.mzML \
  --setnames sampleset1 sampleset2 -o db.sqlite
```

Or, to add spectra to an existing SQLite lookup:

```
msstitch storespectra --dbfile lookup.sqlite --spectra file3.mzML file4.mzML \
  --setnames sampleset2 sampleset3
```

Then store quantification data from dinosaur (MS1 precursor quant) and isobaric 
quantification (including precursor purities, but use centroided MS1 for this)
from OpenMS together with the spectra. MS1 precursor quant features will be aligned
with MS2 spectra by selecting the best m/z match (an m/z window of 2x m/z 
tolerance), from all features within the retention time window of 2x RT tolerance.
--rttol is specified in seconds.

```
msstitch storequant --dbfile db.sqlite --spectra file1.mzML file2.mzML \
  --mztol 10 --mztoltype ppm --rttol 10 \
  --dinosaur file1.dinosaur file2.dinosaur \
  --isobaric file1.consensusXML file2.consensusXML
```

When using Hardklor/Kronik instead of Dinosaur, you can instead use:

```
msstitch storequant --dbfile db.sqlite --spectra file1.mzML file2.mzML \
  --kronik file1.kronik file2.kronik \
  --isobaric file1.consensusXML file2.consensusXML
```

For both Dinosaur and Kronik, the MS1 peak sum is used which theoretically would be more correct
when having differently shaped envelopes. If you'd rather use the envelope apex, pass `--apex`
in the above command.


### Handling MS search engines
Create a decoy database where peptides are reversed between tryptic residues. Decoy
peptides will be shuffled if they match a target sequence, but if they after shuffling
a number of times (max as in `--maxshuffle`) still match a target sequence, they will 
be removed. To avoid removal (e.g. for keeping database sizes identical), pass `--keep-target'.

```
msstitch makedecoy uniprot.fasta -o decoy.fasta --scramble tryp_rev --maxshuffle 10
```

Or without even trying to shuffle peptide sequences that match to the target DB:

```
msstitch makedecoy uniprot.fasta -o decoy.fasta --scramble tryp_rev --ignore-target-hits
```


After running two samples of MSGF and percolator, we can start making 
a more proper set of PSM tables by adding percolator data and filtering
on FDR. You need percolator data, a PSM table and a matching mzIdentML file.
The PSM table / mzIdentML should not have been filtered differently, as each PSM
is expected to be a single item in the mzIdentML file and vice versa. PSMs not found
in the percolator data are not output. The following adds percolator svm-score,
q-value (FDR), and posterior error as columns to the PSM table:

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

This also adds ExplainedIonCurrentRatio, which is parsed from the mzIdentML file as delivered by MSGF+.

If you're running Sage instead of MSGFplus, you will not need the mzIdentML, so the `perco2psm`
commands above can be written e.g.:
```
msstitch perco2psm -i psms1.sage.txt -o psms.sage.percolated.txt \
  --perco percolator1.xml \
  --filtpsm 0.01 --filtpep 0.01
... etc
```

Now refine the PSM tables, using the earlier created SQLite DB, 
adding more information (sample name, MS1 precursor quant,
isobaric quant, proteingroups, genes). In this example we set isobaric
quantitation intensities to NA if the precursor purity measured is <0.3.

```
cp db.sqlite decoy_db.sqlite
msstitch psmtable -i target.tsv -o target_psmtable.txt --fasta uniprot.fasta \
  --dbfile db.sqlite --addmiscleav --addbioset --ms1quant --isobaric \
  --min-precursor-purity 0.3 --proteingroup --genes
msstitch psmtable -i decoy.tsv -o decoy_psmtable.txt --fasta decoy.fasta \
  --dbfile decoy_db.sqlite --proteingroup --genes --addbioset
```

When using Sage results, you will have missed cleavages out of the box (and retention time, ion mobility)! So leave
out `--addmiscleav`. Otherwise msstitch will accept the Sage file format as well.


If necessary (e.g. multiple TMT sample sets), split the table before making
protein/peptide tables:

```
msstitch split -i target_psmtable.txt --splitcol bioset
```

### Summarizing PSMs
Create a peptide table, with summarized median isobaric quant ratios,
highest MS1 intensity PSM as the peptide MS1 quant intensity, and an additional
linear-modeled q-value column:

```
msstitch peptides -i set1_target_psms.txt -o set1_target_peptides.txt \
  --scorecolpattern svm --modelqvals --ms1quant \
  --isobquantcolpattern tmt10plex --denompatterns _126 _127C 
```

The same peptide table can also be made using median sweeping, which takes the median
intensity channel for each PSM as a denominator. Here is also exemplified how
to do channel-median centering of the ratios to normalize and use log2 intensity
values before calculating ratios:

```
msstitch peptides -i set1_target_psms.txt -o set1_target_peptides.txt \
  --scorecolpattern svm --modelqvals --ms1quant \
  --isobquantcolpattern tmt10plex --mediansweep --logisoquant --median-normalize
```

Or, if you only want the median PSM intensity per peptide summarized, use `--medianintensity`
Here is also illustrated that you can use the --keep-psms-na-quant flag to NOT
throw out the PSMs which have isobaric intensity below the mininum intensity 
(default 0, here 100) IN ANY channel:

```
msstitch peptides -i set1_target_psms.txt -o set1_target_peptides.txt \
  --scorecolpattern svm --modelqvals --ms1quant \
  --isobquantcolpattern tmt10plex --medianintensity \
  --minint 100 --keep-psms-na-quant
```

In case of analyzing peptides with PTMs, you may want to process a subset of
PSMs (those with the PTMs) to create a separate peptide table from. In
that case, there is an option to divide (or subtract for log2 data) isobaric 
quant values to a protein (or gene) table from a non-PTM search, often done on
another non-enriched sample. This allows discerning PTM-peptide 
differential expression from its respective protein differential expression in the sample.
The protein/gene table should obviously contain the same samples/channel,
and for example be from an `msstitch proteins` or `msstitch isosummarize` command,
using `--median-normalize` to get median centered ratios for the proteins or genes.
After that, create a peptide table from PTM-PSMs as follows:

```
msstitch peptides -i set1_ptm_psms.txt -o set1_ptm_peptides.txt \
  --scorecolpattern svm --isobquantcolpattern tmt10plex --denompatterns _126 _127C \
  --logisoquant --totalproteome set1_proteins.txt
```

For proper normalizing of this table (as it would otherwise be impacted by
sample differences per channel), you may want to median-center. In the case of 
small and possibly noisy PTM tables, it can be advisable to use another table
from a global search (or e.g. the full peptide or protein table from the PTM search)
for determining channel medians. This is possible by specifying the above command
plus:

```
  --median-normalize --normalization-factors-table /path/to/set1_global_proteins.txt
```

To create a protein table, with isobaric quantification as for peptides, the
average of the top-3 highest intensity peptides for MS1 quantification.
For all of these, summarizing isobaric PSM data to peptide, protein, gene features 
is done using medians of log2 PSM quantification values per feature (e.g. a protein). If you'd
rather use averages, use `--summarize-average` as below, where we also show log2
transformation of intensities before summarizing and subsequent median-centering.
FDR (q-values) for the protein table is here calculated 'classically', by ranking
target and decoy proteins, taking their best scoring peptide's q-value as a score.
`--logscore` is used since q-value is used (higher is better when comparing peptides).

```
msstitch proteins -i set1_target_peptides.txt --decoyfn set1_decoy_peptides \
  --psmtable set1_target_psms.txt \
  -o set1_proteins.txt \
  --scorecolpattern '^q-value' --logscore \
  --ms1quant \
  --isobquantcolpattern tmt10plex --denompatterns _126 _127C \
  --summarize-average --logisoquant --median-normalize
```

Or the analogous process for genes, using median sweeping to get intensity ratios instead of denominators:
As for peptides above, one can use the --keep-psms-na-quant flag to NOT
throw out the PSMs which have isobaric intensity below the mininum intensity
(default 0 used here) in any channel. Here we use picked FDR (Savitski et al. 2015 MCP)
to define q-values for the genes, for which you need a target and decoy fasta to form pairs.
The fasta files need to be analogous, i.e. the same order for matching T/D pair genes.

```
msstitch genes -i set1_target_peptides.txt --decoyfn set1_decoy_peptides \
  --psmtable set1_target_psms.txt \
  -o set1_genes.txt \
  --scorecolpattern '^q-value' --logscore \
  --fdrtype picked --targetfasta tdb.fa --decoyfasta decoy.fa \
  --ms1quant \
  --isobquantcolpattern tmt10plex --mediansweep \
  --keep-psms-na-quant
```

Or when there are ENSEMBL entries in the fasta search database, even for ENSG, here with summarized median PSM intensity per ENSG:

```
msstitch ensg -i set1_target_peptides.txt --decoyfn set1_decoy_peptides \
  --psmtable set1_target_psms.txt \
  -o set1_ensg.txt \
  --scorecolpattern '^q-value' --logscore \
  --ms1quant \
  --isobquantcolpattern tmt10plex --medianintensity \
  --median-normalize
```

Finally, merge multiple sets of proteins (or genes/ENSG) into a single output.
Here we set an cutoff so that features with FDR > 0.01 are set to NA for the 
respective sample set.

```
msstitch merge -i set1_proteins.txt set2_proteins.txt -o protein_table.txt \
  --setnames sampleset1 sampleset2 \
  --dbfile db.sqlite \
  --fdrcolpattern 'q-value' --mergecutoff 0.01 \
   --ms1quantcolpattern area --isobquantcolpattern plex
```


### Some other useful commands
Trypsinize a fasta file (minimum retained peptide length, do cut K/RP, allow 1 missed cleavage, 
also generate peptides from protein N-term with methionine loss).
This also treats stop codons as separators and will not read through them if they are in a protein.
Pass `--ignore-stop-codons` if you would like to treat an `*` as any other amino acid.

```
msstitch trypsinize -i uniprot.fasta -o tryp_up.fasta --minlen 7 \
  --cutproline --miscleav 1 --nterm-meth-loss
```

Create an SQLite file with tryptic sequences for filtering out e.g. known-sequence data.
Options and behaviour as for `trypsinize`, plus `--insourcefrag` which builds a lookup with support for 
in-source fragmented peptides that have lost some N-terminal residues. If one instead of 
in-source fragmentation only wishes to include protein N-term methionine loss, use `--nterm-meth-loss`

```
msstitch storeseq -i canonical.fa -o tryptic.sqlite --cutproline --minlen 7 \
  --miscleav 1 --insourcefrag
```

Filter a percolator output file, or a PSM file using the created SQLite,
removing sequences
that match those stored in the SQLite. The below also removes sequences in the 
sample which are deamidated (i.e. D -> N), and sequences that have lost at most
2 N-terminal amino acids due to in-source fragmentation (DB must have been 
built with support for that).

```
# Percolator:
msstitch filterperco -i perco.xml --dbfile tryptic.sqlite \
  --insourcefrag 2 --deamidate -o filtered.xml

# PSM file:
msstitch seqfilt -i psms.txt --dbfile tryptic.sqlite \
  --insourcefrag 2 --deamidate -o filtered.psms.txt
```

Now that you have filtered percolator, the svm-score populations will have changed,
and you may want to rerun its `qvality` component for the PSM and peptide score populations
to obtain fresh FDRs. To put those into the PSM table, you can use:

```
msstitch perco2psm -i psms.txt \
  --perco filtered.xml --mzid psms1.mzIdentML \
  --qvalitypsms qvality_psms.txt --qvalitypeps qvality_peptides.txt \
  --filtpsm 0.01 --filtpep 0.01
```
Note that the qvality PSMs and peptides files containing svm-scores and FDR, have to
be made of the filtered.xml percolator output, as their content will be zipped, thus
each line in the qvality PSMs file represents a PSM present in percolator XML, 
and vice versa. The same goes for peptides.

Create an SQLite file with full-protein sequences for filtering any peptide of 
a minimum length specified that matches to those. Stop codons are also here treated
as separating peptides. Slower than filtering tryptic sequences but more comprehensive:

```
msstitch storeseq -i canonical.fa -o proteins.sqlite --fullprotein --minlen 7
```

Filter a percolator output or PSM file on protein sequences using the SQLite, removing 
sequences in sample which match to anywhere in the protein. Sequences may be 
deamidated, and minimum length parameter must match the one the database is 
built with.

```
# Percolator
msstitch filterperco -i perco.xml --dbfile proteins.sqlite \
  --fullprotein --deamidate --minlen 7 -o filtered.xml

# PSM file:
msstitch seqfilt -i psms.txt --dbfile proteins.sqlite \
  --fullprotein --deamidate --minlen 7 -o filtered.psms.txt
```

Using a similar DB from `storeseq` with `--map-accessions`, we can add a column in a PSM file which contains
sequence-matches from a fasta file, so any peptide matching these sequences can
get annotated with the fasta ID for that sequences as provided in e.g. `external.fa`:

```
msstitch storeseq -i external.fa -o tryptic.sqlite --cutproline --minlen 7 \
  --miscleav 1 --insourcefrag --map-accessions

# Or use a full protein DB to match non-tryptic peptides:
msstitch storeseq -i external.fa -o fullprotein.sqlite --fullprotein \
  --minlen 7 --map-accessions

# Now add the annotations
msstitch seqmatch -i psms.txt --dbfile tryptic.sqlite --matchcolname COLUMN_NAME \
  --insourcefrag 3 --deamidate -o annotated_psms.txt

# For the fullprotein DB
msstitch seqmatch -i psms.txt --dbfile fullprotein.sqlite --matchcolname ANOTHER_COLUMN_NAME \
  --minlen 7 --fullprotein --enforce-tryptic -o annotated_psms.txt
```

Split a percolator file with PSMs and peptides into files with specific protein headers.
The below will split perco.xml and output two files: perco.xml_h0.xml, containing all
PSMs/peptides that have at least one mapping to a `novel_p` protein, and perco.xml_h1.xml,
which will contain all PSMs/peptides mapping to either `lncRNA` or `intergenic` proteins,
but which will __NOT__ contain PSMs/peptides mapping to also either/or `ENSP` and `sp`. 
More files, _h2.xml etc can be created by adding more headers.

```
msstitch splitperco -i perco.xml --protheaders "novel_p" "known:ENSP;sp|novel:lncRNA;intergenic"`
```

Remove duplicate peptides from percolator output (if youve introduced these somehow)
by only retaining the best peptides:

```
msstitch dedupperco -i perco.xml -o dedup_peptides.xml
```

Remove duplicate PSMs from the output if you've introduced those, for really edge cases 
like removing identical bare peptides where the only difference is 
a modification position.
```
msstitch deduppsms -i psms.txt --peptidecolpattern 'Peptide' -o dedup_psms.txt
```

Create an isobaric ratio table median-summarizing the PSMs by any column number 
you want in a PSM table. E.g. you have added a column with exons. The following 
uses average of two channels as denominator, outputs a new table with first column
the features found in column nr.20 of the PSM table:

```
msstitch isosummarize -i psm_table.txt --featcol 20 \
  --isobquantcolpattern tmt10plex --denompatterns 126 127C
```

We can also use this command to create a table multi-mapping PSMs are counted towards
identification and quant of all its mappings, e.g. proteins separated by `;`. This is
not recommended for conventional protein table construction, use at your own risk
in cases that benefit from this behaviour.

```
msstitch isosummarize -i psm_table.txt --featcol 20 \
  --isobquantcolpattern tmt10plex --denompatterns 126 127C \
  --split-multi-entries
```





Re-use an earlier PSM table and add PSMs from searched spectra files of a new or
re-searched sample. Saves time so you won't have to re-search all the spectra 
in case of a big analysis. In the example below, new PSMs are the result of a 
sample set that has been re-searched, (e.g. when MS reruns are done in case of 
bad spectra), so we delete the existing sample set before continuing. 
Protein grouping is done after regenerating the PSM table, to illustrate you 
can do protein grouping on the entire table
instead of only on the sample set. Since the new table is the one which supplies
the header, the columns not supplied in the command (here protein groups) will
be removed from the final result. This function assumes all PSMs presented are
in the same order in the table, so they should not have been inserted in parallel,
safest is to not generate the lookup table by hand.

```
msstitch deletesets -i old_psmtable.txt -o cleaned_psmtable.txt \
    --dbfile db.sqlite --setnames bad_set
msstitch psmtable -i rerun_target.tsv --oldpsms cleaned_psmtable.txt \
   -o new_almost_done_psmtable.txt --fasta uniprot.fasta \
  --dbfile db.sqlite --addmiscleav --addbioset --ms1quant --isobaric
msstitch psmtable -i new_almost_done_psmtable.txt -o new_target_psms.txt \
  --proteingroup
```

It is also possible to only pass a PSM table to `deletesets`:

```
msstitch deletesets -i old_psmtable.txt -o cleaned_psmtable.txt \
  --setnames bad_set
```
