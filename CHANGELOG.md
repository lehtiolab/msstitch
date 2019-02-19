# Changelog

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
