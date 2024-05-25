# Fixture creation script
# Runs inside tests/
set -euo pipefail

# check if we are in tests/ and have a virtualenv with ansible
if [[ ! $(git rev-parse --show-prefix) = 'tests/' ]]
then
    echo You are not in the git repo tests folder, exiting
    exit 1
fi
echo Checking virtualenv installation
python3 -m venv .venv-fixtures
source .venv-fixtures/bin/activate
pip install -e ../
echo OK

# DB td generate
msstitch makedecoy -i base_fixtures/ens99_small.fasta -o fixtures/protrev_ens99_small.fasta --scramble prot_rev
rm decoychecker.sqlite
cat base_fixtures/ens99_small.fasta fixtures/protrev_ens99_small.fasta > /tmp/tmp_td.fa

# Run MSGF, perco, fake qvality
docker run -e UID=$(id -u) -e GID=$(id -g) -tiv $(pwd):/test -v /tmp:/tmp quay.io/biocontainers/msgf_plus:2020.03.14--0 bash -c \
	'mkdir -p /tmp/msgf && rm -rf /tmp/msgf/* && \
	msgf_plus -d /tmp/tmp_td.fa -s /test/base_fixtures/few_spectra.mzML -o /test/fixtures/few_spectra.mzid -mod /test/base_fixtures/mods.txt -tda 0 \
	    -t 10ppm -ti -1,2 -m 0 -inst 3 -e 1 -protocol 0 -ntt 2 -minLength 7 -maxLength 50 -minCharge 2 -maxCharge 6 -n 1 -addFeatures 1 && \
	    msgf_plus edu.ucsd.msjava.ui.MzIDToTsv -i /test/fixtures/few_spectra.mzid -o /tmp/msgf/msgf.tsv && \
	chown $UID:$GID /tmp/msgf/msgf.tsv'
docker run -e UID=$(id -u) -e GID=$(id -g) -tiv $(pwd):/test -v /tmp:/tmp quay.io/biocontainers/percolator:3.6.5--h6351f2a_0 bash -c \
	'mkdir -p /tmp/perco && rm -rf /tmp/perco/* && \
	msgf2pin -o /tmp/perco/msgf.pin -e trypsin -P "decoy_" /test/fixtures/few_spectra.mzid && \
	percolator -j /tmp/perco/msgf.pin -X /test/fixtures/perco.xml -N 500000 --decoy-xml-output && chown $UID:$GID /test/fixtures/perco.xml && \
	chown $UID:$GID /test/fixtures/perco.xml'

cat <(echo $'Score\tPEP\tq-value') \
       	<(grep -A3 '<psm p' fixtures/perco.xml| grep -v psm| grep -v '\-\-' | \
       	sed -E 's/\s*<[\/]*[a-z_]+>//g' | paste - - - | \
       	gawk -v FS='\t' -v OFS='\t' '{print $1, $3, $2+0.1}') > fixtures/qvality_psms.txt 
cat <(echo $'Score\tPEP\tq-value') \
       	<(grep -A3 '<peptide p' fixtures/perco.xml| grep -v peptide| grep -v '\-\-' | \
       	sed -E 's/\s*<[\/]*[a-z_]+>//g' | paste - - - | \
       	gawk -v FS='\t' -v OFS='\t' '{print $1, $3, $2}') > fixtures/qvality_peps.txt 

# TIMS MSGF-perco-quant
msstitch storespectra --spectra base_fixtures/few_spec_timstof.mzML --setnames Set1 -o fixtures/spectra_lookup_timstof.sqlite
docker run -e UID=$(id -u) -e GID=$(id -g) -tiv $(pwd):/test -v /tmp:/tmp quay.io/biocontainers/msgf_plus:2020.03.14--0 bash -c \
	'mkdir -p /tmp/tims && rm -rf /tmp/tims/* && \
	head -n2 /test/base_fixtures/mods.txt > /tmp/tims/timsmods && \
       	msgf_plus -d /tmp/tmp_td.fa -s /test/base_fixtures/few_spec_timstof.mzML \
        	-o /test/fixtures/few_spec_timstof.mzid \
		-mod /tmp/tims/timsmods \
		-tda 0 -t "10.0ppm" -ti "-1,2" -m 0 -inst 2 -e 1 -protocol 0 -ntt 2 \
		-minCharge 2 -maxCharge 6 -n 1 -addFeatures 1 && \
	msgf_plus edu.ucsd.msjava.ui.MzIDToTsv -i /test/fixtures/few_spec_timstof.mzid \
	    -o /test/fixtures/few_spec_timstof.tsv && \
	chown $GID:$UID /test/fixtures/few_spec_timstof.tsv &&
	chown $GID:$UID /test/fixtures/few_spec_timstof.mzid'

docker run -e UID=$(id -u) -e GID=$(id -g) -tiv $(pwd):/test -v /tmp:/tmp quay.io/biocontainers/percolator:3.6.5--h6351f2a_0 bash -c \
       'msgf2pin -o /tmp/tims.pin -e trypsin -P "decoy_" /test/fixtures/few_spec_timstof.mzid && \
       percolator -j /tmp/tims.pin -X /test/fixtures/perco_timstof.xml -N 500000 --decoy-xml-output && \
	chown $GID:$UID /test/fixtures/perco_timstof.xml'
docker run -e UID=$(id -u) -e GID=$(id -g) -tiv $(pwd):/test -v /tmp:/tmp quay.io/biocontainers/openms:3.1.0--h191ead1_4 bash -c \
 	'IsobaricAnalyzer -type tmt16plex -in /test/base_fixtures/few_spec_timstof.mzML -out /test/fixtures/few_spec_timstof.consXML \
 	    -extraction:select_activation "any" \
 	    -extraction:reporter_mass_shift 0.0013 \
 	    -extraction:min_precursor_intensity 1.0 \
 	    -extraction:keep_unannotated_precursor true \
 	    -quantification:isotope_correction true && \
 	chown $UID:$GID /test/fixtures/few_spec_timstof.consXML'

# Sage TIMS
docker run -e UID=$(id -u) -e GID=$(id -g) -tiv $(pwd):/test -v /tmp:/tmp ghcr.io/lazear/sage:v0.14.7 bash -c \
	'mkdir -p /tmp/timssage && rm -rf /tmp/timssage/* && \
	sage -f /tmp/tmp_td.fa -o /tmp/timssage /test/base_fixtures/sage_conf.tims.json /test/base_fixtures/few_spec_timstof.mzML && \
	chown $UID:$GID /tmp/timssage/results.sage.tsv'
cp /tmp/timssage/results.sage.tsv fixtures/few_spec_timstof.sage.tsv

# Run sage, percolator, fake qvality
docker run -e UID=$(id -u) -e GID=$(id -g) -tiv $(pwd):/test -v /tmp:/tmp ghcr.io/lazear/sage:v0.14.7 bash -c \
	'mkdir -p /tmp/sage && rm -rf /tmp/sage/* && \
	sage -f /tmp/tmp_td.fa -o /tmp/sage --write-pin /test/base_fixtures/sage_conf.json /test/base_fixtures/few_spectra.mzML && \
	chown $UID:$GID /tmp/sage/results.sage.tsv'
docker run -e UID=$(id -u) -e GID=$(id -g) -tiv $(pwd):/test -v /tmp:/tmp quay.io/biocontainers/percolator:3.6.5--h6351f2a_0 bash -c \
	'percolator -j /tmp/sage/results.sage.pin -X /test/fixtures/perco.sage.xml -N 500000 --decoy-xml-output && \
	chown $UID:$GID /test/fixtures/perco.sage.xml'

cat <(echo $'Score\tPEP\tq-value') \
       	<(grep -A3 '<psm p' fixtures/perco.sage.xml| grep -v psm| grep -v '\-\-' | \
       	sed -E 's/\s*<[\/]*[a-z_]+>//g' | paste - - - | \
       	gawk -v FS='\t' -v OFS='\t' '{print $1, $3, $2+0.1}') > fixtures/qvality_psms.sage.txt 
cat <(echo $'Score\tPEP\tq-value') \
       	<(grep -A3 '<peptide p' fixtures/perco.sage.xml| grep -v peptide| grep -v '\-\-' | \
       	sed -E 's/\s*<[\/]*[a-z_]+>//g' | paste - - - | \
       	gawk -v FS='\t' -v OFS='\t' '{print $1, $3, $2}') > fixtures/qvality_peps.sage.txt 


# post process add set
head -n1 /tmp/msgf/msgf.tsv | awk -F $'\t' '{OFS=FS ; print $0, "Biological set"}' > fixtures/few_spectra.tsv
tail -n+2 /tmp/msgf/msgf.tsv | awk -F $'\t' '{OFS=FS ; print $0, "Set1"}' >> fixtures/few_spectra.tsv
head -n1 /tmp/sage/results.sage.tsv | awk -F $'\t' '{OFS=FS ; print $0, "Biological set"}' > fixtures/few_spectra.sage.tsv
tail -n+2 /tmp/sage/results.sage.tsv | awk -F $'\t' '{OFS=FS ; print $0, "Set1"}' >> fixtures/few_spectra.sage.tsv


# Few spectra.tsv will have a bio set, add perco data
msstitch perco2psm --perco fixtures/perco.xml -i fixtures/few_spectra.tsv --mzid fixtures/few_spectra.mzid -o fixtures/few_spectra.tsv_fdr.tsv
msstitch perco2psm --perco fixtures/perco.sage.xml -i fixtures/few_spectra.sage.tsv -o fixtures/few_spectra.sage.fdr.tsv
msstitch perco2psm --perco fixtures/perco_timstof.xml -i fixtures/few_spec_timstof.tsv --mzid fixtures/few_spec_timstof.mzid -o fixtures/few_spec_timstof.tsv_fdr.tsv


# Add a set
cat fixtures/few_spectra.tsv_fdr.tsv <(sed 's/few_spectra/set2/;s/Set1/Set2/' < fixtures/few_spectra.tsv_fdr.tsv | tail -n+2) > /tmp/few_spectra.td.tsv
cat fixtures/few_spectra.sage.fdr.tsv <(sed 's/few_spectra/set2/;s/Set1/Set2/' < fixtures/few_spectra.sage.fdr.tsv | tail -n+2) > /tmp/few_spectra.sage.td.tsv
msstitch split -i /tmp/few_spectra.sage.td.tsv --splitcol TD
mv target.tsv fixtures/target.sage.tsv
mv decoy.tsv /tmp/decoy.sage.tsv
msstitch split -i /tmp/few_spectra.td.tsv --splitcol TD
mv decoy.tsv /tmp/decoy.tsv
mv target.tsv fixtures/target.tsv
msstitch splitperco -i fixtures/perco.xml -d fixtures --protheaders "^ENSP000002" "^ENSP000002|^ENSP000003"

# Quant etc:
rm -f /tmp/set2.mzML && ln -s $(pwd)/base_fixtures/few_spectra.mzML /tmp/set2.mzML
rm -f fixtures/spectra_lookup.sqlite && msstitch storespectra --spectra base_fixtures/few_spectra.mzML /tmp/set2.mzML base_fixtures/few_spec_timstof.mzML --setnames Set1 Set2 Set3 -o fixtures/spectra_lookup.sqlite
cp fixtures/spectra_lookup.sqlite fixtures/quant_lookup.sqlite
cp fixtures/spectra_lookup.sqlite fixtures/dinoquant_lookup.sqlite

docker run -e UID=$(id -u) -e GID=$(id -g) -tiv $(pwd):/test -v /tmp:/tmp quay.io/biocontainers/hardklor:2.3.2--hdbdd923_4 bash -c \
	'rm -rf /tmp/hardklor.out && \
	hardklor /test/base_fixtures/hardklor.conf'
docker run -e UID=$(id -u) -e GID=$(id -g) -tiv $(pwd):/test -v /tmp:/tmp quay.io/biocontainers/kronik:2.20--h4ac6f70_6 bash -c \
	'kronik -c 5 -d 3 -g 1 -m 8000 -n 600 -p 10 /tmp/hardklor.out /test/fixtures/few_spectra.kr && \
	chown $UID:$GID /test/fixtures/few_spectra.kr'
rm -f /tmp/set2.kr && ln -s $(pwd)/fixtures/few_spectra.kr /tmp/set2.kr
docker run -e UID=$(id -u) -e GID=$(id -g) -tiv $(pwd):/test -v /tmp:/tmp quay.io/biocontainers/proteowizard:3.0.9992 bash -c   \
	'rm -f /tmp/centroidms1.mzML && cd /tmp && \
	msconvert /test/base_fixtures/few_spectra.mzML --filter "peakPicking true 1" --outfile centroidms1.mzML'
docker run -e UID=$(id -u) -e GID=$(id -g) -tiv $(pwd):/test -v /tmp:/tmp quay.io/biocontainers/openms:2.5.0--h4afb90d_6 bash -c \
	'IsobaricAnalyzer -type tmt10plex -in /tmp/centroidms1.mzML -out /test/fixtures/few_spectra.consXML \
	    -extraction:select_activation "High-energy collision-induced dissociation" \
	    -extraction:reporter_mass_shift 0.0013 \
	    -extraction:min_precursor_intensity 1.0 \
	    -extraction:keep_unannotated_precursor true \
	    -quantification:isotope_correction true && \
	chown $UID:$GID /test/fixtures/few_spectra.consXML'

# For some reason I cant get the biocontainer thing to work, maybe it is too minimal and missing a font or something
docker run -e UID=$(id -u) -e GID=$(id -g) -tiv $(pwd):/test -v /tmp:/tmp msstitch_testdata:3.8 bash -c \
 	'rm -rf /tmp/dino && mkdir -p /tmp/dino && cd /tmp/dino && \
 	dinosaur --verbose --writeMsInspect --minCharge=2 --outDir=/tmp/dino /test/base_fixtures/few_spectra.mzML && \
 	mv /tmp/dino/few_spectra.features.tsv /test/fixtures/few_spectra.dino && \
 	chown $UID:$GID /test/fixtures/few_spectra.dino'

msstitch storequant --dbfile fixtures/quant_lookup.sqlite --kronik fixtures/few_spectra.kr /tmp/set2.kr --spectra base_fixtures/few_spectra.mzML /tmp/set2.mzML --mztol 20.0 --mztoltype ppm --rttol 5.0 --isobaric fixtures/few_spectra.consXML fixtures/few_spectra.consXML
msstitch storequant --dbfile fixtures/dinoquant_lookup.sqlite --dinosaur fixtures/few_spectra.dino --spectra base_fixtures/few_spectra.mzML --mztol 20.0 --mztoltype ppm --rttol 5.0 --isobaric fixtures/few_spectra.consXML


# PSM table MSGF
cat fixtures/quant_lookup.sqlite |tee fixtures/target_psms.sqlite fixtures/target_psms.sage.sqlite /tmp/decoylup /tmp/decoylup.sage > /dev/null
msstitch psmtable -i "fixtures/target.tsv" --dbfile "fixtures/target_psms.sqlite" --fasta base_fixtures/ens99_small.fasta --isobaric --ms1quant --proteingroup --genes --fastadelim pipe --genefield 2 -o fixtures/target_pg.tsv
msstitch psmtable -i "/tmp/decoy.tsv" --dbfile /tmp/decoylup --fasta fixtures/protrev_ens99_small.fasta --fastadelim pipe --genefield 2 --proteingroup --genes -o /tmp/decoy_pg.tsv
sed -i 's/0.0'\$'\\ttarget/0.001\\ttarget/' fixtures/target_pg.tsv
grep -v set2 fixtures/target_pg.tsv > /tmp/set1_target_pg.tsv
grep -v set2 /tmp/decoy_pg.tsv > /tmp/set1_decoy

# PSM table sage
msstitch psmtable -i "fixtures/target.sage.tsv" --dbfile "fixtures/target_psms.sage.sqlite" --fasta base_fixtures/ens99_small.fasta --isobaric --ms1quant --proteingroup --genes --fastadelim pipe --genefield 2 -o fixtures/target_pg.sage.tsv --spectracol 5
msstitch psmtable -i "/tmp/decoy.sage.tsv" --dbfile /tmp/decoylup.sage --fasta fixtures/protrev_ens99_small.fasta --fastadelim pipe --genefield 2 --proteingroup --genes -o /tmp/decoy_pg.sage.tsv --spectracol 5
grep -v set2 fixtures/target_pg.sage.tsv > /tmp/set1_sage
grep -v set2 /tmp/decoy_pg.sage.tsv > /tmp/set1_sage.decoy
# Interlace columns with msgf to have same quant, q values etc
# First cut the MSGF (longer) table:
head -n $(wc -l /tmp/set1_sage | sed 's/ .*//') /tmp/set1_target_pg.tsv > fixtures/set1_target_pg.tsv
head -n $(wc -l /tmp/set1_sage.decoy | sed 's/ .*//') /tmp/set1_decoy > /tmp/set1_decoy_cut
# Target peptides get quant and the actual peptide sequence, which we modify to look like sage output (almost, except the four digits become 3)
paste <(cut -f1,3,6 /tmp/set1_sage) <(cut -f12,14-19 fixtures/set1_target_pg.tsv | sed -E 's/(\+[0-9]+\.[0-9]+)/[\1]/g;s/^(\[[+0-9.]+\])/\1-/;s/Peptide/peptide/' ) <(cut -f12-42,49 /tmp/set1_sage) <(cut -f24- fixtures/set1_target_pg.tsv) > fixtures/set1_target_pg.sage.tsv
# Decoy peptides do not get quant of course
paste <(cut -f1,3,6 /tmp/set1_sage.decoy) <(cut -f12,14-19 /tmp/set1_decoy_cut | sed -E 's/(\+[0-9]+\.[0-9]+)/[\1]/g;s/^(\[[+0-9.]+\])/\1-/;s/Peptide/peptide/' ) <(cut -f12-42 /tmp/set1_sage.decoy) <(cut -f24- /tmp/set1_decoy_cut) > /tmp/set1_decoy_pg.sage.tsv

# make proteins for isosum feats
msstitch isosummarize -i fixtures/set1_target_pg.tsv -o fixtures/proteins_quantonly.txt --featcol 14 --isobquantcolpattern tmt10plex --denompatterns _126
msstitch isosummarize -i fixtures/set1_target_pg.tsv --featcol 11 -o fixtures/isosum_charge_column.txt --isobquantcolpattern plex --denompatterns 126 --keep-psms-na-quant
msstitch isosummarize -i fixtures/set1_target_pg.tsv -o fixtures/proteins_quantonly_splitmulti.txt --featcol 13 --isobquantcolpattern tmt10plex --denompatterns _126 --split-multi-entries

# peptide tables:
msstitch peptides -i fixtures/set1_target_pg.tsv -o fixtures/target_peptides.tsv --scorecolpattern svm --isobquantcolpattern plex --ms1quantcolpattern area --denompatterns 126
msstitch peptides -i fixtures/set1_target_pg.tsv -o fixtures/target_pep_quant_norm.tsv --scorecolpattern svm --isobquantcolpattern plex --ms1quantcolpattern area --denompatterns 126 --median-normalize
msstitch peptides -i fixtures/set1_target_pg.tsv -o fixtures/target_peptides_avg.tsv --scorecolpattern svm --summarize-average --isobquantcolpattern plex --ms1quantcolpattern area --denompatterns 126
msstitch peptides -i fixtures/set1_target_pg.tsv -o fixtures/target_peptides_sweep.tsv --scorecolpattern svm --isobquantcolpattern plex --ms1quantcolpattern area --mediansweep --keep-psms-na-quant

# sage peptides (for protein building etc)
msstitch peptides -i fixtures/set1_target_pg.sage.tsv -o fixtures/target_peptides.sage.tsv --scorecolpattern svm --isobquantcolpattern plex --ms1quantcolpattern area --denompatterns 126

# Decoy peptides
msstitch peptides -i /tmp/set1_decoy_cut -o fixtures/decoy_peptides.tsv --scorecolpattern svm
msstitch peptides -i /tmp/set1_decoy_pg.sage.tsv -o fixtures/decoy_peptides.sage.tsv --scorecolpattern svm

# protein tables
msstitch proteins -i fixtures/target_peptides.tsv -o fixtures/proteins.txt --psmtable fixtures/set1_target_pg.tsv --isobquantcolpattern plex --denompatterns 126 --ms1quant --scorecolpattern '^q-value' --logscore --decoyfn fixtures/decoy_peptides.tsv
msstitch proteins -i fixtures/target_peptides.tsv -o fixtures/proteins_isonorm_log.txt --psmtable fixtures/set1_target_pg.tsv --isobquantcolpattern plex --denompatterns 126 --logisoquant --median-normalize --ms1quant --scorecolpattern '^q-value' --logscore --decoyfn fixtures/decoy_peptides.tsv
msstitch proteins -i fixtures/target_peptides.tsv -o fixtures/proteins_isonorm_nolog.txt --psmtable fixtures/set1_target_pg.tsv --isobquantcolpattern plex --mediansweep --median-normalize --ms1quant --scorecolpattern '^q-value' --logscore --decoyfn fixtures/decoy_peptides.tsv
msstitch proteins -i fixtures/target_peptides.tsv -o fixtures/proteins_nonorm_log.txt --psmtable fixtures/set1_target_pg.tsv --isobquantcolpattern plex --denompatterns 126 --logisoquant --scorecolpattern '^q-value' --logscore --decoyfn fixtures/decoy_peptides.tsv
msstitch proteins -i fixtures/target_peptides.tsv -o fixtures/proteins_avg.txt --psmtable fixtures/set1_target_pg.tsv --isobquantcolpattern plex --denompatterns 126 --ms1quant --scorecolpattern '^q-value' --logscore --decoyfn fixtures/decoy_peptides.tsv --summarize-average
msstitch proteins -i fixtures/target_peptides.tsv -o fixtures/proteins_sweep.txt --psmtable fixtures/set1_target_pg.tsv --isobquantcolpattern plex --mediansweep --ms1quant --scorecolpattern '^q-value' --logscore --decoyfn fixtures/decoy_peptides.tsv --keep-psms-na-quant
msstitch proteins -i fixtures/target_peptides.tsv -o fixtures/proteins_intensities.txt --psmtable fixtures/set1_target_pg.tsv --isobquantcolpattern plex --medianintensity --ms1quant --scorecolpattern '^q-value' --logscore --decoyfn fixtures/decoy_peptides.tsv
msstitch proteins -i fixtures/target_peptides.tsv -o fixtures/proteins_denoms_pickfdr.txt --psmtable fixtures/set1_target_pg.tsv --isobquantcolpattern plex --denompatterns 126 --ms1quant --scorecolpattern '^q-value' --logscore --decoyfn fixtures/decoy_peptides.tsv --fdrtype picked --targetfasta base_fixtures/ens99_small.fasta --decoyfasta fixtures/protrev_ens99_small.fasta

# Peptides which need protein tables
msstitch peptides -i fixtures/set1_target_pg.tsv -o fixtures/target_peptides_totalprotnorm.txt --scorecolpattern svm --isobquantcolpattern plex --ms1quantcolpattern area --denompatterns 126 --logisoquant --totalproteome fixtures/proteins_isonorm_log.txt
msstitch peptides -i fixtures/set1_target_pg.tsv -o fixtures/target_peptides_totalp_sepnorm.txt --scorecolpattern svm --isobquantcolpattern plex --ms1quantcolpattern area --denompatterns 126 --logisoquant --totalproteome fixtures/proteins_nonorm_log.txt --normalization-factors-table fixtures/proteins_nonorm_log.txt --median-normalize
msstitch peptides -i fixtures/set1_target_pg.tsv -o fixtures/target_peptides_totalprotnorm_nolog.txt --scorecolpattern svm --isobquantcolpattern plex --denompatterns 126 --totalproteome fixtures/proteins.txt


# Genes
msstitch genes -i fixtures/target_peptides.tsv -o fixtures/genenames.txt --psmtable fixtures/set1_target_pg.tsv --isobquantcolpattern plex --denompatterns 126 --ms1quant \
  --scorecolpattern '^q-value' --logscore --fdrtype picked --decoyfn fixtures/decoy_peptides.tsv --targetfasta base_fixtures/ens99_small.fasta --decoyfasta fixtures/protrev_ens99_small.fasta

msstitch ensg -i fixtures/target_peptides.tsv -o fixtures/ensg.txt --psmtable fixtures/set1_target_pg.tsv --isobquantcolpattern plex --denompatterns 126 --ms1quant \
  --scorecolpattern 'svm' --logscore --fdrtype picked --decoyfn fixtures/decoy_peptides.tsv --targetfasta base_fixtures/ens99_small.fasta --decoyfasta fixtures/protrev_ens99_small.fasta

msstitch ensg -i fixtures/target_peptides.tsv -o /tmp/tmptxt --psmtable fixtures/set1_target_pg.tsv --isobquantcolpattern plex --denompatterns 126 --ms1quant \
  --scorecolpattern 'svm' --logscore --decoyfn fixtures/decoy_peptides.tsv --targetfasta base_fixtures/ens99_small.fasta --decoyfasta fixtures/protrev_ens99_small.fasta
cols=$(head -n1 /tmp/tmptxt | tr '\t' '\n' | grep -nv quanted | cut -f 1 -d ':' | tr '\n' ',' | sed 's/\,$//')
cut -f $cols /tmp/tmptxt > fixtures/ensg_nopsms.txt

# Fasta stuff
msstitch makedecoy -i base_fixtures/twoproteins.fasta -o fixtures/decoy_twoproteins.fasta --scramble prot_rev 
msstitch trypsinize -i base_fixtures/twoproteins.fasta -o fixtures/trypsinized_twoproteins.fasta
msstitch makedecoy -i fixtures/trypsinized_twoproteins.fasta -o fixtures/decoy_tryprev_pretryp_twoproteins.fasta --scramble tryp_rev --notrypsin 
mv decoychecker.sqlite fixtures/decoycheck.sqlite
msstitch  makedecoy -i base_fixtures/twoproteins.fasta -o fixtures/decoy_tryprev_twoproteins.fasta --scramble tryp_rev --ignore-target-hits
msstitch makedecoy -i base_fixtures/twoproteins.fasta -o fixtures/decoy_tryprev_targetcheck_twoproteins.fasta --scramble tryp_rev 
msstitch  makedecoy -i base_fixtures/twoproteins.fasta -o fixtures/decoy_tryprev_minlen_twoproteins.fasta --scramble tryp_rev --minlen 5
msstitch makedecoy -i base_fixtures/twoproteins.fasta -o fixtures/decoy_tryprev_tcheck_keeptarget_twoproteins.fasta --scramble tryp_rev --keep-target
msstitch makedecoy -i fixtures/trypsinized_twoproteins.fasta -o fixtures/decoy_tryprev_pretryp_keeptarget_twoproteins.fasta --scramble tryp_rev --notrypsin --keep-target

# duplicate thing
cat fixtures/few_spectra.tsv <(tail -n1 fixtures/few_spectra.tsv) <(head -n2 fixtures/few_spectra.tsv | tail -n1) > fixtures/few_spectra_duplicate.tsv

rm decoychecker.sqlite
rm normalization_factors_set1_target_pg.tsv
echo Done! Created all fixtures!
