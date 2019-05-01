#!/bin/bash

#SBATCH --partition=genetics_1          # Partition (job queue)
#SBATCH --requeue                    	# Return job to the queue if preempted
#SBATCH --job-name=dana(de_novo) 		# Assign an short name to your job
#SBATCH --nodes=1                    	# Number of nodes you require
#SBATCH --ntasks=1                   	# Total # of tasks across all nodes
#SBATCH --cpus-per-task=28            	# Cores per task (>1 if multithread tasks)
#SBATCH --mem=100G                   	# Real memory (RAM) required (MB)
#SBATCH --time=48:00:00              	# Total run time limit (HH:MM:SS)
#SBATCH --output=job_%j.out     		# STDOUT output file
#SBATCH --error=job_%j.err      		# STDERR output file (optional)
#SBATCH --export=ALL                	# Export you current env to the job env

# ...................................................................................
# ...................................................................................
# 	SHELL FOR FINDING DE NOVO REPEAT ELEMENTS IN DROSOPHILA SPECIES.
# 			SINGLE SPECIES AT A TIME
# ...................................................................................
# ...................................................................................
# CHECK STATUS OF JOB
# while true; do date; squeue -u cpr74; sleep 2;  var1="$(squeue -u cpr74 | wc -l)"; echo -n "number_of_jobs = "; expr $var1 - 1; printf "\n"; sleep 3; done

# OPEN INTERACTIVE NODE
# srun --partition=genetics_1 --nodes=1 --ntasks=1 --cpus-per-task=28 --mem=30G --time=01:00:00 --export=ALL --pty bash -i

#

# ...................................................................................
# 	CHECKLIST BEFORE RUNNING
# ...................................................................................


# STEP 0:
# MAKING AND CHANGING TO THE DIRECTORY OF INTEREST TO BE DONE BY HAND

# mkdir ./$species_directory/
# cd ./$species_directory/

# [ ] ARE ALL THE VARIABLES PROPERLY NAMED FOR THE SCRIPT TO WORK CORRECTLY?
# [ ] IS THE blastx_parser.py IN THE CORRECT SPECIES DIRECTORY?
# [ ] DOES THAT DIRECTORY HAVE THE CORRECT SPECIES FATA FILE IN IT?

# ...................................................................................
# 	VARIABLE NAMES
# ...................................................................................

#
#                                ,,            ,,          ,,
#                                db           *MM        `7MM
#                                              MM          MM
# `7M'   `MF' ,6"Yb.  `7Mb,od8 `7MM   ,6"Yb.   MM,dMMb.    MM   .gP"Ya  ,pP"Ybd
#   VA   ,V  8)   MM    MM' "'   MM  8)   MM   MM    `Mb   MM  ,M'   Yb 8I   `"
#    VA ,V    ,pm9MM    MM       MM   ,pm9MM   MM     M8   MM  8M"""""" `YMMMa.
#     VVV    8M   MM    MM       MM  8M   MM   MM.   ,M9   MM  YM.    , L.   I8
#      W     `Moo9^Yo..JMML.   .JMML.`Moo9^Yo. P^YbmdP'  .JMML. `Mbmmd' M9mmmP'
#
#

# Directory
# THE DIRECTORY THAT IS MADE TO HOSUE THE FUNCTIONS OF THE PARTICULAR SPECIES.
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
# Drosophila_persimilis
# Drosophila_pseudoobscura
# Drosophila_sechellia
# Drosophila_simulans
# Drosophila_virilis
# Drosophila_willistoni
# Drosophila_yakuba
# Drosophila_melanogaster

# LOG FILE
log_file=dana_de_novo_log.txt
date_time=date

# BuildDatabase
# THIS FILE IS A VARIABLE INPUT FOR BuildDatabase
fasta_file=dana_nanopore
raw_file=nanopore/dana.pass.minimap2.racon.x3.pilon.x3.fasta

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

#
#               ,,                        ,,    ,,
#               db                      `7MM    db
#                                         MM
#  `7MMpdMAo. `7MM  `7MMpdMAo.  .gP"Ya    MM  `7MM  `7MMpMMMb.   .gP"Ya
#    MM   `Wb   MM    MM   `Wb ,M'   Yb   MM    MM    MM    MM  ,M'   Yb
#    MM    M8   MM    MM    M8 8M""""""   MM    MM    MM    MM  8M""""""
#    MM   ,AP   MM    MM   ,AP YM.    ,   MM    MM    MM    MM  YM.    ,
#    MMbmmd'  .JMML.  MMbmmd'   `Mbmmd' .JMML..JMML..JMML  JMML. `Mbmmd'
#    MM               MM
#  .JMML.           .JMML.

sleep 3
cd /scratch/cpr74/2019Spring/drosophila_classifications/dana_nanopore_2019_03_05

# Deletion of temporary files required to get the directories and files set up
# rm ./../Amarel.zip ; sleep 5
# rm -r ./../__MACOSX/ ; sleep 5

ln -s /scratch/cpr74/2019Spring/drosophila_classifications/genomes/$raw_file ./$fasta_file

# IF NANOPORE FILE, THEN CONDENSE THE SEQID INTO 7-DIGIT CODE ONLY
#	IF THE NAME OF FASTA FILE HAS NANOPORE, THEN IT IS A NANOPORE SEQUENCE AND CAN RUN THE CONDENSOR ON IT
sleep 3
if [[ "$fasta_file" == *"nanopore"* ]]; then
	addition="_corrected"; # sleep 2
  	python3 nanopore_id_condensor.py -i $fasta_file -o $fasta_file$addition; # sleep 2
   	echo -e "Nanopore Sequence detected and corrected...\n"
   	fasta_file="$fasta_file$addition"; # sleep 5
else
   	echo -e "Illumina sequence detected; left unchanged...\n"
fi

echo -e "Fasta_File is now $fasta_file . \n\n\n"

# # STEP 1:
# # BUILDING DATABASE USING RepeatModeler
sleep 3
BuildDatabase -name $db_name -engine ncbi $fasta_file

# # STEP 2:
# RUNNING RepeatModeler
sleep 3
time ~/pkg/RepeatModeler-open-1.0.11/RepeatModeler -database $db_name -engine ncbi -pa 28 > $log_file; sleep 30

# STEP 3:
# RUNNING UCLUST
sleep 5
temp="$(find . -name "consensi.fa.classified")"
dir="$(cut -d "/" -f2 <<< "$temp")"
sleep 5
mv $dir outfiles
cp $peptide_seqs outfiles
cp $fasta_file outfiles
cd outfiles
sleep 5
usearch -cluster_fast consensi.fa.classified -id 0.9 -centroids $uclust_output
# FINDING THE NUMBER OF SEQUENCES THAT WERE CLUSTERED
echo -n 'The number of sequences from the {consensi.fa.classified} file = '
less consensi.fa.classified | grep -c '>'
echo -n 'The number of sequences from the {uclust_output}  file = '
less $uclust_output | grep -c '>'
var1="$(less consensi.fa.classified | grep -c '>')"
var2="$(less $uclust_output | grep -c '>')"
echo -n 'The number of sequeces that were clustered is: '
expr $var1 - $var2


# STEP4:
# MAKING A PEPTIDE DATABASE
mv ../$peptide_seqs ./
makeblastdb -in $peptide_seqs -dbtype prot
# RUNNING BLASTX (UCLUSTED dNTP VS PROTEINS) TO REMOVE ALL REPEATS THAT MIGHT BE GENES
#OUTPUTS ONLY THE QUERY ID, SEQUENCE ID QUERY LENGTH AND ALIGHNMENT LENGTH IN THAT ORDER
blastx -query $uclust_output -db $peptide_seqs -evalue 1e-10 -out $blastx_output -outfmt "7 qseqid sseqid qstart qend sstart send qlen"
#SORTING AND GROUPING THAT OUTPUT
awk '!/# /' $blastx_output | sort > sorted_blastx.out
# cat sorted_blastx.out | bedtools groupby -g 1,2 -c 3,4 -ops distinct,sum > grouped_blastx.out

# STEP 5:
# THIS BLASTX_PARSER.PY COMBINES THE NAME OF THE REPEAT AND THE NAME OF THE SEQUENCE IT ALIGNS TO.
# IT THEN PRINTS OUT THE NAME; ALIGHNMENT LENGTH; QUERY START; QUERY END; QUERY LENGTH
python3 ../blastx_parser.py -s sorted_blastx.out -u $uclust_output -b $bedtools_input
#SORTING AND THEN MERGING REPEATS THAT ALIGN TO THE SAME GENE MULTIPLE TIMES
cat $bedtools_input | sort -k1,1 -k2,2n | bedtools merge -c 4 -o distinct > $merged_alignments
#SUBTRACTING THE QSTART FROM QEND TO FIND THE LENGTH OF THE ALIGNMENT
awk '{ print $1, $3-$2, $4 }' $merged_alignments > $shortened_merged
# COMPARING THE QLENGTH AND THE CALCULATED Q LENGTH TO SEE WHICH OF THE GENES ARE REPEATS
awk '{ print $1, $2/$3 }' $shortened_merged > repeat_gene_alignments.bed

# STEP 6:
# PUT THIS repeat_gene_alignments.bed INTO ANOTHER PARSER THAT DOES THE FOLLOWING:
# 1. GOES THROUGH EACH OF THE ITEMS
# 2. SEPARATES THE NAME OF THE REPEAT AND THE NAME OF THE GENE
# 3. IF THERE ARE MULTIPLE PAIRS, THEN ADDS THE ALIGHNMENT SCORES AND RETURNS ONLY 1 PAIRS
# 4. RETURNS SEQUENCE OF REPEATS THAT ARE GENES AND NOT GENES (LIKELY)
# 5. CONTAINS THE CUTOFF VALUE AS A VARIABLE. CAN BE CHANGED AT A LATER DATE IF NEED BE.
python3 ../blastx_combined_name_separator.py -c $cutoff -s repeat_gene_alignments.bed -u $uclust_output -p $uclust_no_genes -g $uclust_genes

# STEP 9:
# REPEATMASKER
# TO FIGURE OUT IF RMASK ILBRARY FOR DROSOPHILA (REPBASE) IF ANY TE DESCRIBED WITHIN DROSOPHILA ARE IN THERE.
mkdir $dir_RepeatMasker
~/pkg/RepeatMasker/RepeatMasker -e ncbi -pa 28 -nolow -norna -species drosophila -no_is -dir $dir_RepeatMasker -gff $uclust_no_genes
# COPY THE REPEATMASKER DATA BACK TO THE MAIN DIRECTORY
cp ./$dir_RepeatMasker/$uclust_no_genes.out ./$repeatmasker_gff

# STEP 10:
# SHORTEN REPEATMASKER DATA AND RUN BEDTOOLS MERGE AGAIN TO MERGE THE SEQUENCES TOGETHER
# 	ONLY INCLUDE THE {NAME OF REPEAT_|_|_NAME OF ANNOTATION, START OF REPEAT ALIGNMENT, END OF REPEAT ALIGNMENT}
touch $repeatmasker_shortened_gff # TOUCHING SO CREATES A NEW EMPTY FILE AND DOESNT OUTPUT ERROR ON NEXT LINE.
python3 ../repeat_masker_shortener.py -i $repeatmasker_gff -o $repeatmasker_shortened_gff; sleep 2
#  	RUN BEDTOOLS TO MERGE OVERLAPPING ALIGNMENTS
cat $repeatmasker_shortened_gff | sort -k1,1 -k2,2n | sort -k1,1 -k2,2n | bedtools merge > $repeatmasker_shortened_merged_bed
#  	FINDING ALIGNMENTS THAT WERE MORE THAN _0.75_ TO BE ABLE TO ANNOTATE THEM TO KNOWN TES
python3 ../repeat_masker_parser.py -r $repeatmasker_shortened_merged_bed -t $threshold -u $uclust_no_genes -p $uclust_postmatch -s $uclust_summary; sleep 2

# STEP 10:
# MAKING A TABLE THAT HAS THE FREQENCIES OF EACH COUNTS OF FAMILIES FROM THE REPEATMODELER DATASET AS WELL AS THAT WHICH HAS BEEN 			CLASSIFIED BY REPEATMASKER.
# CONDENSING THE CLASSES OF EACH TE WITHIN THE SPECIES/RUN
python3 ../table_maker.py -c consensi.fa.classified -p $uclust_postmatch -f $family_table; sleep 5

# STEP 11:
# FINDING THE FREQUENCIES OF EACH OF THE TRANSPOSONS IDENTIFIED BY THE PIPELINE
mkdir $dir_transposon_freq # MAKE THE DIRECTORY (PREVENT CREATION OF MISLABELED DIRECTORY BY RepeatMasker)
~/pkg/RepeatMasker/RepeatMasker -e ncbi -pa 28 -nolow -norna -no_is -dir $dir_transposon_freq -lib $uclust_postmatch -gff ./../$fasta_file; sleep 3
cp ./$dir_transposon_freq/$fasta_file.out.gff ./$transposon_freq; sleep 3

# STEP 12
# RUNNING TRF TO FIND THE REPEATING SEQUENCES.
# 	INDEPENDENT OF THE PIPELINE, SO CAN RUN WHENEVER.
# 	CHOSE TO DO WITHIN PIPELINE BECAUSE CAN THEN RUN WITH EACH SPECIES.
trf ../$fasta_file 2 7 7 80 10 50 500 -h

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


# STEP -2
# CLEANUP
rm -r round-* ; sleep 30
rm ./../trfResults-* ; sleep 15

#
#
#
#
#  ,pP"Ybd `7MM  `7MM  `7MMpMMMb.pMMMb.  `7MMpMMMb.pMMMb.   ,6"Yb.  `7Mb,od8 `7M'   `MF'
#  8I   `"   MM    MM    MM    MM    MM    MM    MM    MM  8)   MM    MM' "'   VA   ,V
#  `YMMMa.   MM    MM    MM    MM    MM    MM    MM    MM   ,pm9MM    MM        VA ,V
#  L.   I8   MM    MM    MM    MM    MM    MM    MM    MM  8M   MM    MM         VVV
#  M9mmmP'   `Mbod"YML..JMML  JMML  JMML..JMML  JMML  JMML.`Moo9^Yo..JMML.       ,V
#                                                                               ,V
#                                                                            OOb"

# ...................................................................................
# 	SUMMARY FILE
# ...................................................................................

# STEP n:
# SUMMARIZING THE DATA FOR GENES, NOT GENES AND CLUSTERED FROM UCLUST.
# 	HEADER
echo "===========================================================================================" > $summary_file
echo "===========================================================================================" >> $summary_file
echo "" >> $summary_file
echo -ne " \t\t\t SUMMARY FILE FOR RUN OF " >> $summary_file
echo -e "$species ." >> $summary_file
echo "" >> $summary_file
echo "" >> $summary_file
echo -e "Things to note about this file: " >> $summary_file
echo -e "\t This file is NOT PARSABLE." >> $summary_file
echo -e "\t Has data of the number of copies of a repeat identified (de novo) by RepeatModeler." >> $summary_file
echo -e "\t Based on custom pipeline created by @AUTHOR." >> $summary_file
echo "" >> $summary_file
echo -e " Includes the following:" >> $summary_file
echo -e "\t 1. Summary of unclustered and clustered repeats via UCLUST." >> $summary_file
echo -e "\t 2. Summary of identified genes based on alignment scores using BLASTX." >> $summary_file
echo -e "\t 3. Summary of identified novel repeats based on alignment scores using BLASTX." >> $summary_file
echo "" >> $summary_file
echo -e "\t\t\t\t\t\t @AUTHOR: \t\t Chinmay Rele" >> $summary_file
echo -e "\t\t\t\t\t\t @LAST MODIFIED: \t $date_time" >> $summary_file
echo -e "\t\t\t\t\t\t @Date of run: \t\t $date_time" >> $summary_file
echo "" >> $summary_file
echo "===========================================================================================" >> $summary_file
echo -e "===========================================================================================\n\n" >> $summary_file
# CALCULATION AND SUMMARY
# UCLUST
echo -e "\t UCLUST" >> $summary_file
echo "" >> $summary_file
echo -n 'The number of sequences from the {consensi.fa.classified} file = ' >> $summary_file
less consensi.fa.classified | grep -c '>' >> $summary_file
echo -n "The number of sequences from the {$uclust_output} file = " >> $summary_file
less $uclust_output | grep -c '>' >> $summary_file
var1="$(less consensi.fa.classified | grep -c '>')" >> $summary_file
var2="$(less $uclust_output | grep -c '>')" >> $summary_file
echo -n 'The number of sequeces that were clustered is: ' >> $summary_file
expr $var1 - $var2 >> $summary_file
echo "" >> $summary_file
echo "" >> $summary_file
transposon_bp="$(cat transposon_freq.gff | grep -v "##" | cut -f4,5  | awk '{ print $2-$1 }' | awk '{total += $0} END {print total}')"
echo "Total TE bps in $species genome = $transposon_bp ." >> $summary_file
merged_TE_BP="$(cat transposon_freq.gff | grep -v "##" | cut -f1,4,5 | sort -k1,1 -k2,2n | bedtools merge | awk '{ print $3-$2 }' | awk '{total += $0} END {print total}')"
echo "Merged TE bps in $species genome = $merged_TE_BP ." >> $summary_file
trf_bp="$(cat dana_nanopore_corrected.2.7.7.80.10.50.500.dat | grep -v "Sequence" | grep -v "Parameters" | grep -v "\n" | grep -v -e '^$' | cut -f1,2 | awk '{ print $2-$1 }' | awk '{total += $0} END {print total}')"
echo "Simple Repeat bps in $species genome = $trf_bp ." >> $summary_file
assembly_size="$(grep -v ">" ../dana_nanopore_corrected | wc | awk '{print $3-$1}')"
echo "Assembly size of $species in bp = $assembly_size ." >> $summary_file
echo "" >> $summary_file
echo "" >> $summary_file


echo "===========================================================================================" >> $summary_file
echo "" >> $summary_file
# 	BLASTX
# 		NO GENES
echo -ne "\t NO GENES" >> $summary_file
echo -e "\t cutoff = $cutoff" >> $summary_file
echo "" >> $summary_file
echo -n "Size of {repeat_gene_alignments.bed} (uncondensed) = " >> $summary_file
a="$(less repeat_gene_alignments.bed | grep -c "_|_|_|_|_")"
echo $a >> $summary_file
echo -e "\t This is the unsorted and not-uniqued." >> $summary_file
echo -e "\t Means that there are multiples of the same gene-repeat alignment scores." >> $summary_file
echo -n "Size of {repeat_gene_alignments.bed} (condensed) = " >> $summary_file
b="$(less repeat_gene_alignments.bed | sort -k1,1n | uniq | grep -c "_|_|_|_|_" )"
echo $b >> $summary_file
echo -e "\t Condensed number that show as a single repeat-gene comparison." >> $summary_file
echo -n "Number of repeated Repeat-Gene alignments = " >> $summary_file
expr $a - $b >> $summary_file
echo -n "The number of sequences in the {$uclust_no_genes} file = " >> $summary_file
less $uclust_no_genes | grep -c '>' >> $summary_file
echo "" >> $summary_file
echo "NO Gene: IDs" >> $summary_file
cat $uclust_no_genes | grep ">" | sed 's/(.*//' | column -x -c 120 >> $summary_file
echo "" >> $summary_file
echo "" >> $summary_file
echo "===========================================================================================" >> $summary_file
echo "" >> $summary_file
#		GENES
echo -ne "\t GENES" >> $summary_file
echo -e "\t\t cutoff = $cutoff" >> $summary_file
echo "" >> $summary_file
echo -n "Size of {repeat_gene_alignments.bed} (uncondensed) = " >> $summary_file
a="$(less repeat_gene_alignments.bed | grep -c "_|_|_|_|_")"
echo $a >> $summary_file
echo -e "\t This is the unsorted and not-uniqued." >> $summary_file
echo -e "\t Means that there are multiples of the same gene-repeat alignment scores." >> $summary_file
echo -n "Size of {repeat_gene_alignments.bed} (condensed) = " >> $summary_file
b="$(less repeat_gene_alignments.bed | sort -k1,1n | uniq | grep -c "_|_|_|_|_" )"
echo $b >> $summary_file
echo -e "\t Condensed number that show as a single repeat-gene comparison." >> $summary_file
echo -n "Number of repeated Repeat-Gene alignments = " >> $summary_file
expr $a - $b >> $summary_file
echo -n "The number of sequences in the {$uclust_genes} file = " >> $summary_file
less $uclust_genes | grep -c '>' >> $summary_file
echo "" >> $summary_file
echo "Gene: IDs" >> $summary_file
cat $uclust_genes | grep ">" | sed 's/(.*//' | column -x -c 120 >> $summary_file
echo "" >> $summary_file
echo "" >> $summary_file

echo "===========================================================================================" >> $summary_file
echo "===========================================================================================" >> $summary_file
echo "" >> $summary_file
echo -e "\t\t\tDIRECTORY AND FILES SUMMARIES\n" >> $summary_file
echo -e "Working directory:" >> $summary_file
echo -ne "\t" >> $summary_file
pwd >> $summary_file
echo "" >> $summary_file
echo "===========================================================================================" >> $summary_file
echo -e "\t\t\tPython Scripts\n" >> $summary_file
echo -ne "nanopore_id_condensor.py \n\tCondenses seq ids of nanopore data. \n\tConsensus_Consensus_Consensus_utg000004l_pilon_pilon_pilon\n\t\tto\n\tutg000004l" >> $summary_file
echo -ne "blastx_parser.py\n\tTakes in uncommented BLASTX data and compares to UCLUST data.\n\tMerges the names of the query and subject so can use [bedtools merge] to merge overlapping alignments.\n\n" >> $summary_file
echo -ne "blastx_combined_name_separator.py \n\tGoes through each of the items. \n\tSeparates the name of the repeat and the name of the gene. \n\tIf there are multiple pairs, then adds the alighnment scores and returns only 1 pairs. \n\tReturns sequence of repeats that are genes and not genes (likely). \n\tContains the cutoff value as a variable. Can be changed at a later date if need be.\n\n" >> $summary_file
echo -ne "repeat_masker_shortener.py\n\tCombines data from RepeatMasker into a single sequenceid so we can run mergeBed on it.\n\n" >> $summary_file
echo -ne "repeat_masker_parser.py\n\tReads in the repeatmasker output to be able to align them and see sequence alignment similarities to the repeatmasker library.\n\n" >> $summary_file
echo -ne "table_maker.py\n\tTakes in consensi.fa.cassified and $uclust_post_match and finds the number of counts of the repeats in each family.\n\tWill be done for RepeatModeler data as well as RepeatMasker data.\n\n" >> $summary_file
echo "" >> $summary_file
echo "===========================================================================================" >> $summary_file
echo "" >> $summary_file
echo -e "\t\t\tFile Summaries\n" >> $summary_file; sleep 1
echo -e "\n$fasta_file\n\tRaw genetic data of $species. \n\tThe raw file that is the variable in put for BuildDatabase.\n \tGot from FTP client at FlyBase.org. \n\t[ftp://ftp.flybase.net/releases/current/] \n\tThe above line is the directory of all Drosophila genomes.\n" >> $summary_file; sleep 1
echo -e "$db_name \n \tThe database that BuildDatabase makes.\n" >> $summary_file; sleep 1
echo -e "$log_file \n \tRepeatModeler log file.\n \tCheck here for errors if job failed.\n" >> $summary_file; sleep 1
echo -e "$uclust_output \n \tUCLUST output.\n \tHas clustered repeats (Those that would occur multiple times in the RepeatModeler output).\n" >> $summary_file; sleep 1
echo -e "$blastx_output \n \tUCLUSTed file ran against $species peptide database." >> $summary_file; sleep 1

# ---------------------------
# 	END OF FILE
# ---------------------------
