# Changelog

## [2.17] - 2019-11-
### Added
- FDR calculator from percolator output for PSM tables from MSGF, `msspsmtable percolator`
- Missed cleavage count column in PSM table as an option to `msspsmtable specdata --miscleav`

### Changed
- Linear modeled q-values from scores and bestpeptides for a feature are NA if there are fewer than 1 q-values above the threshold.
- Add biological set to PSM table is now optional with `msspsmtable specdata --addbioset` flag.
- Protein FDR calculation does not crash on ZeroDivisionError when T=0: D/T is then set to 1. 

### Fixed
- Bugs in `msslookup makedecoy` fixed, actually filter out all non-matching peptides and too short peptides
- Protein FDR calculation was for some proteins calculated as D/T+D where that should be D/T


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
