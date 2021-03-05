# Changelog

## [3.6] - 2021-03-05
## Added
- Total protein normalization of PTM peptides
- PTM FLR stored and output in peptide table merges
- Precursor ion fraction (purity) from OpenMS to PSM table
- Precursor purity filter sets PSM quant to NA, keeps PSMs

## Fixed
- Do not create duplicate isobaric channels in DB and with that in PSM output
- Do not add standard fields to PSM refined header that are already there
- Do not crash on peptide table merge when a peptide has e.g. 600 protein matches
- Merging no-group peptide tables no longer output multiple identical protein accessions for a peptide
- Peptide table genecentric merge gets all genes per peptide, not only one


## [3.5] - 2020-09-23
## Changed
- Isobaric quant summarizing to features now throws out all PSMs with NA in a channel by default, keep them with --keep-psms-na-quant


## [3.4] - 2020-09-17
### Added
- Isobaric quant summarizing using denominators OR sweep OR median intensity
- Median normalization possible for isobaric quant
- Column in results "number of fully quanted PSMs", i.e. without missing values
- `msstitch deleteset` can delete a sample set from PSM tables and sqlite lookups
- `msstitch psmtable` can take `--oldpsms` to append new PSMs to an existing table (with identical search etc)

### Fixed
- Bug (and free speed increase) where merged tables the following fields would be all wrong and all the same for each set: amount PSMs, amount peptides, amount unique peptides
- MS1 quant bug where sometimes wrong features where fetched to align to a scan


## [3.3] - 2020-07-24
### Added
- Added `--dinosaur` to store results from MS1 quant using dinosaur, include FWHM data
- Added median-intensity (i.e. not ratio) summarization to `isosummarize` and `proteins` etc

### Changed
- Choose between MS1 apex and MS1 sum for quantification
- Choose between median and average when summarizing isobaric PSM data

### Fixed
- Bug in MS1 feat to PSM scan aligment, `ORDER BY` had been removed in `9717325837853de159ede5e77191b0c7de032820`


## [3.2] - 2020-07-09
### Changed
- `msstitch merge --no-group-annotation` does output a protein(s) column for whatever accession is in the DB

### Fixed
- Corrected typo in column to index on for ensg_proteins table (gene_id), which made merging ENSG very slow


## [3.1] - 2020-07-03
### Fixed
- Quant lookup on consensusXML stored channels in wrong order when >10 channels eg in TMT16
- --spectracol argument in msstitch psmtable did not work


## [3.0] - 2020-06-04
### Breaking
- Completely new interface, one command (msstitch), some subcommands, merged several
functionalities (e.g. SQLite creation/output) into fewer steps.

### Added
- Possible to add ion mobility value of scans to PSM table

### Changed
- Updated README
- Output column name changes (esp Gene ID/ Gene name from symbols)
- peptide sequence in first column of peptide table output
- trypsinize adds peptide index to fasta output as to not create duplicate entries
- PSM table creator does not longer take a mapfn from biomart, all info comes from fasta header
- Consequently, picked FDR can now work on ENSG ID and gene names
- SQLite creation for spectra has possiblity to specify output file 
- Percolator output splitprotein different behaviour: can specify "known" prot headers which excludes all scans with those annotations,
  when not specified it keeps any scan for a header (also if other headers are specified)

### Fixed
- msslookup quant no longer stores duplicates if >500.000 scans, those first 500.000 were stored twice

### Removed
- percolator split/merge/filter/qvality scoring, only kept protein identity splitting
- Protein error probability output


## [2.19] - 2019-12-16
### Fixed
- Protein FDR calculation for ENSG entries was correct but only if the fasta file contained a decoy_ prefix on the ENSG entries, which they did not
  in the case of makedecoy generated.


## [2.18] - 2019-12-04
### Changed
- makedecoy shuffles x times and then does no longer throw out the peptide but keeps the original decoy. Otherwise the decoy database gets smaller

### Fixed
- makedecoy with target checking works also on long proteins, which exceed the max parameter setting for SQLite 


## [2.17] - 2019-11-28
### Added
- FDR calculator from percolator output for PSM tables from MSGF, `msspsmtable percolator`
- Missed cleavage count column in PSM table as an option to `msspsmtable specdata --miscleav`

### Changed
- Linear modeled q-values from scores and bestpeptides for a feature are NA if there are fewer than 10 q-values above the threshold.
- Add biological set to PSM table is now optional with `msspsmtable specdata --addbioset` flag.
- Protein FDR calculation does not crash on ZeroDivisionError when T=0: D/T is then set to 1. 

### Fixed
- Bugs in `msslookup makedecoy` fixed, actually filter out all non-matching peptides and too short peptides
- Protein FDR calculation was for some proteins calculated as D/T+D where that should be D/T
- Peptide and PSM q-value calculation maxes at 1 (in case decoys > targets)


## [2.16] - 2019-10-21
### Added
- Trypsinization tool `msslookup trypsinize` to create trypsinized FASTA from proteins
- Biomart table header addition of "Gene stable ID" to comply with latest
- Params `--maxshuffle`, `--minlen`, and `--notrypsin` added to `msslookup seqspace` for more flexible/quicker operation if needed

### Changed
- Faster SQL for known_searchspace by UNIQUE constraint, accept text files also instead of .fa, 

## [2.15] - 2019-09-30
### Added
- Decoy creation tool to do protein-reverse and tryptic peptide reverse while shuffling decoy peptides that are identical to target peptides

### Fixed
- For some reason had missed that --noncentric also does not like having a coverage table (dependent on proteins, so irrelevant).

## [2.14] - 2019-08-19
### Changed
- Quicker XML parsing with big documents since there was non-cleared elements earlier

## [2.13] - 2019-06-04
### Changed
- Not querying SQLite for each feature/MS2 quant in msslookup quant, instead keep all data for given spectra file in memory, hopefully this speeds things up
### Fixed
- --noncentric flag for building combined output tables did not work if there was not a protein table in the lookup file. Now possible.

## [2.12] - 2018-12-07
### Changed
- Can now build lookups of multiple protein/gene/peptide tables without the nr-of-psms used for isobaric quant, to allow for other quant methods
- Output merged protein/gene tables with all mappings. Whereas protein tables earlier also got gene IDs and gene names, the reverse now is also true
### Fixed
- Header names in merged tables are now not only "Protein accession" anymore, but in case of e.g. ENSEMBL gene IDs the table is called gene ID
- Bug fixed, in merged tables we could get eg last line containing a protein with only NAs because it was not checked if the last merged feature was "empty" 

## [2.11] - 2018-12-07
### Fixed
- Bug fixed, now possible to load protein/gene/peptide tables in sqlite without storing "number of quanted PSMs"

## [2.10] - 2018-12-07
### Fixed
- Bug fixed where isobaric databases would also result in skipping the last line of the PSM table when running `msspsmtable quant`

## [2.9] - 2018-11-22
### Fixed
- Bug fixed where precursor-quant-only databases would result in skipping the last line of the PSM table when running `msspsmtable quant`

## [2.8] - 2018-04-19
### Fixed
- Bug which was not completely fixed in 2.7, now add NA correctly to all fields

## [2.7] - 2018-04-19
### Change
- Field outputs for peptable and prottable merge - used to keep set-fields together and interleave isobaric channels with their nr of psms.
  Now values are kept together S1-qvalue, S2-qvalue, S1-intensity, S2-intensity. Isobaric channels are no longer interleaved, kept in groups of the multiplex amount,
  and the nr of PSM groups follow at the very end.

### Fixed
- Bug where it tried to remove a dummy field from a list of NA
- Bug in removing `(pre=,post=)`, new regex

## [2.6] - 2018-03-15

### Added
- This changelog :)
- Ion injection time is now parsed from the mzML files and stored in the PSM table

### Changed
- Better, more thorough stripping of `(pre=,post=)` string from MSGF entries in peptide table. Useful eg when such are already in shoddy FASTA headers

### Fixed
- Strip also negative weight modifications from peptides when needing the unmodified sequence
