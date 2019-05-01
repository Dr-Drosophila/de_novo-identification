# srun --partition=genetics_1 --nodes=1 --ntasks=1 --cpus-per-task=28 --mem=30G --time=02:00:00 --export=ALL --pty bash -i

species_directory=dana_RepeatModeler
species=Drosophila_ananassae
species_short=dana

# Drosophila_ananassae
# Drosophila_biarmipes
# Drosophila_bipectinata
# Drosophila_erecta
# Drosophila_eugracilis
# Drosophila_triauraria
# Drosophila_mauritiana
# Drosophila_mojavensis
# Drosophila_melanogaster
# Drosophila_persimilis
# Drosophila_pseudoobscura
# Drosophila_sechellia
# Drosophila_simulans
# Drosophila_virilis
# Drosophila_willistoni
# Drosophila_yakuba

# LOG FILE
log_file=dana_de_novo_log.txt
date_time=date

# BuildDatabase
# THIS FILE IS A VARIABLE INPUT FOR BuildDatabase
fasta_file=dana_nanopore_corrected
raw_file=nanopore/Dana.pass.minimap2.racon.x3.pilon.x3.fasta

# THIS IS THE NAME OF THE DATABASE THAT BuildDatabase makes.
db_name=dana_db

# RepeatModeler
# LOG FILE
log_file=dana_run.log

# UCLUST
# OUTPUT FILE NAME
uclust_output=uclusted_dana.fasta

# BLAST
# PEPTIDE SEQUENCE
peptide_seqs=melanogaster-all-translation-r6.24.fasta
#BLASTX OUTPUT
blastx_output=dana_blastx.out

# BEDTOOLS
# THE FILE THAT NEEDS TO BE PASSED INTO THE BEDTOOLS_MERGE; IT IS NOT SORTED.
bedtools_input=bedtools_merge_input.bed
#ALIGNMENTS THAT HAVE BEEN MERGED BY BEDTOOLS MERGE
merged_alignments=merged_alignments.bed
#FILE THAT HAS CORRECT LENGTH OF ALIGNMENT
shortened_merged=shortened_merged.bed

# PYTHON PARSING
# NOT GENES
uclust_no_genes=uclust_no_genes.fasta
# LIKELY GENES
uclust_genes=uclust_genes.fasta
# CUTOFF
cutoff=0.5

# REPEATMASKER DIRECTORY
dir_RepeatMasker=dana_RepeatMasker
repeatmasker_gff=dana_RepeatMasker.gff
repeatmasker_shortened_gff=dana_RepeatMasker_shortened.gff
repeatmasker_shortened_merged_bed=dana_RepeatMasker_shortened_merged.bed
uclust_postmatch=uclust_postmatch.fasta
uclust_summary=uclust_summary.txt
threshold=0.75

# SUMMARY FILE
family_table=family_table.txt
summary_file=dana_summary_file.txt
transposon_freq=transposon_freq.gff
dir_transposon_freq=dana_transposon_freq

echo "All variables written to memory."


# STEP 13
# DOWNSTREAM ANALYSIS
# 	ACCOUNTING FOR UNKNOWN SEQUENCES
# condense trf and rmask data to bed format
python3 ../trf_repmask_condensor.py -t $fasta_file.2.7.7.80.10.50.500.dat -r ./$dir_transposon_freq/$fasta_file.out -o pinko
# intersect two datasets and write
bedtools intersect -a rmask_out_sorted.bed -b trf_dat_sorted.bed -f 0.8 >> pinko
# merge intersection
cat pinko | sort -k1,1 -k2,2n | bedtools merge -c 4 -o distinct > merged_satellites.tab
# merge
python3 ../merged_Sat_to_out.py -p merged_satellites.tab -s shortened_rmask_file.out -o final_output
# condense classes
python3 ../class_condensor.py -i ./final_output -o ./../../identified_TEs/dana_condensed_classes.txt -s dana








cat ../../identified_TEs/dana_condensed_classes.txt | grep "DNA" | awk '{s += ($4)}END{print s}'
cat ../../identified_TEs/dana_condensed_classes.txt | grep "LINE" | awk '{s += ($4)}END{print s}'
cat ../../identified_TEs/dana_condensed_classes.txt | grep "LTR" | awk '{s += ($4)}END{print s}'
cat ../../identified_TEs/dana_condensed_classes.txt | grep "RC" | awk '{s += ($4)}END{print s}'
cat ../../identified_TEs/dana_condensed_classes.txt | grep "Simple_Repeat" | awk '{s += ($4)}END{print s}'
cat ../../identified_TEs/dana_condensed_classes.txt | grep "Satellite" | awk '{s += ($4)}END{print s}'
cat ../../identified_TEs/dana_condensed_classes.txt | grep "Unknown" | awk '{s += ($4)}END{print s}'

grep -v "species" dana_short | cut -f2,4 | sort -k1 | grep "DNA" | awk '{s += ($2)}END{print s}'
grep -v "species" dana_short | cut -f2,4 | sort -k1 | grep "LINE" | awk '{s += ($2)}END{print s}'
grep -v "species" dana_short | cut -f2,4 | sort -k1 | grep "LTR" | awk '{s += ($2)}END{print s}'
grep -v "species" dana_short | cut -f2,4 | sort -k1 | grep "RC" | awk '{s += ($2)}END{print s}'
grep -v "species" dana_short | cut -f2,4 | sort -k1 | grep "Satellite/Simple_Repeat" | awk '{s += ($2)}END{print s}'
grep -v "species" dana_short | cut -f2,4 | sort -k1 | grep "Unknown" | awk '{s += ($2)}END{print s}'

dana    DNA     6379774
dana    LINE    16595768
dana    LTR     29446771
dana    RC      8259686
dana    Simple_repeat   233808
dana    Unknown 9726699

-outfmt "7 qseqid qstart qend"

grep -v "# " blast_unknown | tr ":" "\t" | awk '{print $1"\t"$3"\t"$4}' | sort -k1,1 -k2,2n  | bedtools merge | wc -l


scp ../../trf_repmask_condensor.py cpr74@amarel.rutgers.edu:/home/cpr74
scp ../../class_condensor.py cpr74@amarel.rutgers.edu:/home/cpr74
scp <> cpr74@amarel.rutgers.edu:/home/cpr74
scp <> cpr74@amarel.rutgers.edu:/home/cpr74


grep ">" ./dana_nanopore_2019_03_05/outfiles/consensi.fa.classified | grep -v "?" | wc -l
dna="$(grep ">" ./dana_nanopore_2019_03_05/outfiles/consensi.fa.classified | tr " " "\t" | cut -f1 | grep -c "DNA")"
ltr="$(grep ">" ./dana_nanopore_2019_03_05/outfiles/consensi.fa.classified | tr " " "\t" | cut -f1 | grep -c "LTR")"
line="$(grep ">" ./dana_nanopore_2019_03_05/outfiles/consensi.fa.classified | tr " " "\t" | cut -f1 | grep -c "LINE")"
rc="$(grep ">" ./dana_nanopore_2019_03_05/outfiles/consensi.fa.classified | tr " " "\t" | cut -f1 | grep -c "RC")"
unk="$(grep ">" ./dana_nanopore_2019_03_05/outfiles/consensi.fa.classified | tr " " "\t" | cut -f1 | grep -c "Unknown")"
sr="$(grep ">" ./dana_nanopore_2019_03_05/outfiles/consensi.fa.classified | tr " " "\t" | cut -f1 | grep -c "Simple_repeat")"
sat="$(grep ">" ./dana_nanopore_2019_03_05/outfiles/consensi.fa.classified | tr " " "\t" | cut -f1 | grep -c "Satellite")"
expr $dna + $ltr + $line + $rc
expr $unk
expr $sr + $sat

grep ">" ./dana_nanopore_2019_03_05/outfiles/consensi.fa.classified | tr " " "\t" | cut -f1 | grep -v "Simple_repeat" | grep -v "DNA" | grep -v "LTR"| grep -v "LINE" | grep -v "RC" | grep -v "Unknown"


cat dbia_nanopore_2019_03_05/outfiles/dbia_RepeatMasker.gff | cut -f9 | tr " " "\t" | cut -f2 | grep "Unknown"


# ====================================================================================================================================================================================================================================================================================================================================================================================================

cat dana_nanopore_2019_03_05/outfiles/dana_RepeatMasker_shortened.gff | grep "Unknown"
cat dbia_nanopore_2019_03_05/outfiles/dbia_RepeatMasker_shortened.gff | grep "Unknown"
cat dbip_nanopore_2019_03_05/outfiles/dbip_RepeatMasker_shortened.gff | grep "Unknown"
cat dere_nanopore_2019_03_05/outfiles/dere_RepeatMasker_shortened.gff | grep "Unknown"
cat deug_nanopore_2019_03_05/outfiles/deug_RepeatMasker_shortened.gff | grep "Unknown"
cat dmau_nanopore_2019_03_05/outfiles/dmau_RepeatMasker_shortened.gff | grep "Unknown"
cat dmel_nanopore_2019_03_05/outfiles/dmel_RepeatMasker_shortened.gff | grep "Unknown"
cat dmoj_nanopore_2019_03_05/outfiles/dmoj_RepeatMasker_shortened.gff | grep "Unknown"
cat dper_nanopore_2019_03_05/outfiles/dper_RepeatMasker_shortened.gff | grep "Unknown"
cat dpse_nanopore_2019_03_05/outfiles/dpse_RepeatMasker_shortened.gff | grep "Unknown"
cat dsec_nanopore_2019_03_05/outfiles/dsec_RepeatMasker_shortened.gff | grep "Unknown"
cat dsim_nanopore_2019_03_05/outfiles/dsim_RepeatMasker_shortened.gff | grep "Unknown"
cat dtri_nanopore_2019_03_05/outfiles/dtri_RepeatMasker_shortened.gff | grep "Unknown"
cat dvir_nanopore_2019_03_05/outfiles/dvir_RepeatMasker_shortened.gff | grep "Unknown"
cat dwil_nanopore_2019_03_05/outfiles/dwil_RepeatMasker_shortened.gff | grep "Unknown"
cat dyak_nanopore_2019_03_05/outfiles/dyak_RepeatMasker_shortened.gff | grep "Unknown"

# ====================================================================================================================================================================================================================================================================================================================================================================================================

cat dana_nanopore_2019_03_05/outfiles/uclusted_dana.fasta | grep ">" | tr " " "\t" | cut -f1 | tr "#" "\t" | cut -f2 | tr "/" "\t" | cut -f1  | sort | uniq -c | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v "Other"
cat dbia_nanopore_2019_03_05/outfiles/uclusted_dbia.fasta | grep ">" | tr " " "\t" | cut -f1 | tr "#" "\t" | cut -f2 | tr "/" "\t" | cut -f1  | sort | uniq -c | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v "Other"
cat dbip_nanopore_2019_03_05/outfiles/uclusted_dbip.fasta | grep ">" | tr " " "\t" | cut -f1 | tr "#" "\t" | cut -f2 | tr "/" "\t" | cut -f1  | sort | uniq -c | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v "Other"
cat dere_nanopore_2019_03_05/outfiles/uclusted_dere.fasta | grep ">" | tr " " "\t" | cut -f1 | tr "#" "\t" | cut -f2 | tr "/" "\t" | cut -f1  | sort | uniq -c | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v "Other"
cat deug_nanopore_2019_03_05/outfiles/uclusted_deug.fasta | grep ">" | tr " " "\t" | cut -f1 | tr "#" "\t" | cut -f2 | tr "/" "\t" | cut -f1  | sort | uniq -c | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v "Other"
cat dmau_nanopore_2019_03_05/outfiles/uclusted_dmau.fasta | grep ">" | tr " " "\t" | cut -f1 | tr "#" "\t" | cut -f2 | tr "/" "\t" | cut -f1  | sort | uniq -c | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v "Other"
cat dmel_nanopore_2019_03_05/outfiles/uclusted_dmel.fasta | grep ">" | tr " " "\t" | cut -f1 | tr "#" "\t" | cut -f2 | tr "/" "\t" | cut -f1  | sort | uniq -c | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v "Other"
cat dmoj_nanopore_2019_03_05/outfiles/uclusted_dmoj.fasta | grep ">" | tr " " "\t" | cut -f1 | tr "#" "\t" | cut -f2 | tr "/" "\t" | cut -f1  | sort | uniq -c | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v "Other"
cat dper_nanopore_2019_03_05/outfiles/uclusted_dper.fasta | grep ">" | tr " " "\t" | cut -f1 | tr "#" "\t" | cut -f2 | tr "/" "\t" | cut -f1  | sort | uniq -c | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v "Other"
cat dpse_nanopore_2019_03_05/outfiles/uclusted_dpse.fasta | grep ">" | tr " " "\t" | cut -f1 | tr "#" "\t" | cut -f2 | tr "/" "\t" | cut -f1  | sort | uniq -c | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v "Other"
cat dsec_nanopore_2019_03_05/outfiles/uclusted_dsec.fasta | grep ">" | tr " " "\t" | cut -f1 | tr "#" "\t" | cut -f2 | tr "/" "\t" | cut -f1  | sort | uniq -c | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v "Other"
cat dsim_nanopore_2019_03_05/outfiles/uclusted_dsim.fasta | grep ">" | tr " " "\t" | cut -f1 | tr "#" "\t" | cut -f2 | tr "/" "\t" | cut -f1  | sort | uniq -c | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v "Other"
cat dtri_nanopore_2019_03_05/outfiles/uclusted_dtri.fasta | grep ">" | tr " " "\t" | cut -f1 | tr "#" "\t" | cut -f2 | tr "/" "\t" | cut -f1  | sort | uniq -c | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v "Other"
cat dvir_nanopore_2019_03_05/outfiles/uclusted_dvir.fasta | grep ">" | tr " " "\t" | cut -f1 | tr "#" "\t" | cut -f2 | tr "/" "\t" | cut -f1  | sort | uniq -c | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v "Other"
cat dwil_nanopore_2019_03_05/outfiles/uclusted_dwil.fasta | grep ">" | tr " " "\t" | cut -f1 | tr "#" "\t" | cut -f2 | tr "/" "\t" | cut -f1  | sort | uniq -c | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v "Other"
cat dyak_nanopore_2019_03_05/outfiles/uclusted_dyak.fasta | grep ">" | tr " " "\t" | cut -f1 | tr "#" "\t" | cut -f2 | tr "/" "\t" | cut -f1  | sort | uniq -c | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v "Other"

# ====================================================================================================================================================================================================================================================================================================================================================================================================

echo "dana" > sup.txt
cat dana_nanopore_2019_03_05/outfiles/dana_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -cv "Unknown" >> sup.txt
cat dana_nanopore_2019_03_05/outfiles/dana_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -c "Unknown" >> sup.txt
echo >> sup.txt
echo "dbia" >> sup.txt
cat dbia_nanopore_2019_03_05/outfiles/dbia_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -cv "Unknown" >> sup.txt
cat dbia_nanopore_2019_03_05/outfiles/dbia_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -c "Unknown" >> sup.txt
echo >> sup.txt
echo "dbip" >> sup.txt
cat dbip_nanopore_2019_03_05/outfiles/dbip_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -cv "Unknown" >> sup.txt
cat dbip_nanopore_2019_03_05/outfiles/dbip_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -c "Unknown" >> sup.txt
echo >> sup.txt
echo "dere" >> sup.txt
cat dere_nanopore_2019_03_05/outfiles/dere_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -cv "Unknown" >> sup.txt
cat dere_nanopore_2019_03_05/outfiles/dere_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -c "Unknown" >> sup.txt
echo >> sup.txt
echo "deug" >> sup.txt
cat deug_nanopore_2019_03_05/outfiles/deug_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -cv "Unknown" >> sup.txt
cat deug_nanopore_2019_03_05/outfiles/deug_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -c "Unknown" >> sup.txt
echo >> sup.txt
echo "dmau" >> sup.txt
cat dmau_nanopore_2019_03_05/outfiles/dmau_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -cv "Unknown" >> sup.txt
cat dmau_nanopore_2019_03_05/outfiles/dmau_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -c "Unknown" >> sup.txt
echo >> sup.txt
echo "dmel" >> sup.txt
cat dmel_nanopore_2019_03_05/outfiles/dmel_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -cv "Unknown" >> sup.txt
cat dmel_nanopore_2019_03_05/outfiles/dmel_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -c "Unknown" >> sup.txt
echo >> sup.txt
echo "dmoj" >> sup.txt
cat dmoj_nanopore_2019_03_05/outfiles/dmoj_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -cv "Unknown" >> sup.txt
cat dmoj_nanopore_2019_03_05/outfiles/dmoj_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -c "Unknown" >> sup.txt
echo >> sup.txt
echo "dper" >> sup.txt
cat dper_nanopore_2019_03_05/outfiles/dper_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -cv "Unknown" >> sup.txt
cat dper_nanopore_2019_03_05/outfiles/dper_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -c "Unknown" >> sup.txt
echo >> sup.txt
echo "dpse" >> sup.txt
cat dpse_nanopore_2019_03_05/outfiles/dpse_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -cv "Unknown" >> sup.txt
cat dpse_nanopore_2019_03_05/outfiles/dpse_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -c "Unknown" >> sup.txt
echo >> sup.txt
echo "dsec" >> sup.txt
cat dsec_nanopore_2019_03_05/outfiles/dsec_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -cv "Unknown" >> sup.txt
cat dsec_nanopore_2019_03_05/outfiles/dsec_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -c "Unknown" >> sup.txt
echo >> sup.txt
echo "dsim" >> sup.txt
cat dsim_nanopore_2019_03_05/outfiles/dsim_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -cv "Unknown" >> sup.txt
cat dsim_nanopore_2019_03_05/outfiles/dsim_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -c "Unknown" >> sup.txt
echo >> sup.txt
echo "dvir" >> sup.txt
cat dvir_nanopore_2019_03_05/outfiles/dvir_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -cv "Unknown" >> sup.txt
cat dvir_nanopore_2019_03_05/outfiles/dvir_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -c "Unknown" >> sup.txt
echo >> sup.txt
echo "dwil" >> sup.txt
cat dwil_nanopore_2019_03_05/outfiles/dwil_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -cv "Unknown" >> sup.txt
cat dwil_nanopore_2019_03_05/outfiles/dwil_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -c "Unknown" >> sup.txt
echo >> sup.txt
echo "dyak" >> sup.txt
cat dyak_nanopore_2019_03_05/outfiles/dyak_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -cv "Unknown" >> sup.txt
cat dyak_nanopore_2019_03_05/outfiles/dyak_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -c "Unknown" >> sup.txt
echo >> sup.txt
echo "dtri" >> sup.txt
cat dtri_nanopore_2019_03_05/outfiles/dtri_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -cv "Unknown" >> sup.txt
cat dtri_nanopore_2019_03_05/outfiles/dtri_RepeatMasker.gff | cut -f1,9 | tr " " "\t" | tr "\"" "\t" | tr ":" "\t" | tr "#" "\t" | cut -f2,6 | grep -v "buffer" | grep -v "?" | grep -v "rRNA" | grep -v ")n" | grep -c "Unknown" >> sup.txt
echo >> sup.txt

# ====================================================================================================================================================================================================================================================================================================================================================================================================

cp /scratch/cpr74/2019Spring/drosophila_classifications/dana_nanopore_2019_03_05/outfiles/uclust_postmatch.fasta /projects/genetics/ellison_lab/rele_c/de_novo/postmatch/dana_postmatch.fasta
cp /scratch/cpr74/2019Spring/drosophila_classifications/dbia_nanopore_2019_03_05/outfiles/uclust_postmatch.fasta /projects/genetics/ellison_lab/rele_c/de_novo/postmatch/dbia_postmatch.fasta
cp /scratch/cpr74/2019Spring/drosophila_classifications/dbip_nanopore_2019_03_05/outfiles/uclust_postmatch.fasta /projects/genetics/ellison_lab/rele_c/de_novo/postmatch/dbip_postmatch.fasta
cp /scratch/cpr74/2019Spring/drosophila_classifications/dere_nanopore_2019_03_05/outfiles/uclust_postmatch.fasta /projects/genetics/ellison_lab/rele_c/de_novo/postmatch/dere_postmatch.fasta
cp /scratch/cpr74/2019Spring/drosophila_classifications/deug_nanopore_2019_03_05/outfiles/uclust_postmatch.fasta /projects/genetics/ellison_lab/rele_c/de_novo/postmatch/deug_postmatch.fasta
cp /scratch/cpr74/2019Spring/drosophila_classifications/dmau_nanopore_2019_03_05/outfiles/uclust_postmatch.fasta /projects/genetics/ellison_lab/rele_c/de_novo/postmatch/dmau_postmatch.fasta
cp /scratch/cpr74/2019Spring/drosophila_classifications/dmel_nanopore_2019_03_05/outfiles/uclust_postmatch.fasta /projects/genetics/ellison_lab/rele_c/de_novo/postmatch/dmel_postmatch.fasta
cp /scratch/cpr74/2019Spring/drosophila_classifications/dmoj_nanopore_2019_03_05/outfiles/uclust_postmatch.fasta /projects/genetics/ellison_lab/rele_c/de_novo/postmatch/dmoj_postmatch.fasta
cp /scratch/cpr74/2019Spring/drosophila_classifications/dper_nanopore_2019_03_05/outfiles/uclust_postmatch.fasta /projects/genetics/ellison_lab/rele_c/de_novo/postmatch/dper_postmatch.fasta
cp /scratch/cpr74/2019Spring/drosophila_classifications/dpse_nanopore_2019_03_05/outfiles/uclust_postmatch.fasta /projects/genetics/ellison_lab/rele_c/de_novo/postmatch/dpse_postmatch.fasta
cp /scratch/cpr74/2019Spring/drosophila_classifications/dsec_nanopore_2019_03_05/outfiles/uclust_postmatch.fasta /projects/genetics/ellison_lab/rele_c/de_novo/postmatch/dsec_postmatch.fasta
cp /scratch/cpr74/2019Spring/drosophila_classifications/dsim_nanopore_2019_03_05/outfiles/uclust_postmatch.fasta /projects/genetics/ellison_lab/rele_c/de_novo/postmatch/dsim_postmatch.fasta
cp /scratch/cpr74/2019Spring/drosophila_classifications/dtri_nanopore_2019_03_05/outfiles/uclust_postmatch.fasta /projects/genetics/ellison_lab/rele_c/de_novo/postmatch/dtri_postmatch.fasta
cp /scratch/cpr74/2019Spring/drosophila_classifications/dvir_nanopore_2019_03_05/outfiles/uclust_postmatch.fasta /projects/genetics/ellison_lab/rele_c/de_novo/postmatch/dvir_postmatch.fasta
cp /scratch/cpr74/2019Spring/drosophila_classifications/dwil_nanopore_2019_03_05/outfiles/uclust_postmatch.fasta /projects/genetics/ellison_lab/rele_c/de_novo/postmatch/dwil_postmatch.fasta
cp /scratch/cpr74/2019Spring/drosophila_classifications/dyak_nanopore_2019_03_05/outfiles/uclust_postmatch.fasta /projects/genetics/ellison_lab/rele_c/de_novo/postmatch/dyak_postmatch.fasta
