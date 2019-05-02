#!/bin/bash

#######################################################################################################
##SET UP OPTIONS FOR MAKEBLASTDB AND BLASTP

while getopts a:b:c:d:e:f:g:h:k:l:m:n:o:p:q:r:s:t:u:v:x:y:z option
do
        case "${option}"
        in

                a) database=${OPTARG};;
                b) experimental=${OPTARG};;
                c) transcript_peps=${OPTARG};;
                d) gaf_db=${OPTARG};;
                f) max_matches=${OPTARG};;
                e) E_value=${OPTARG};;
                g) percID=${OPTARG};;
               # h) help=${OPTARG};; ADD USAGE STATEMENT
		m) perc_pos=${OPTARG};;
		o) out=${OPTARG};;
		s) bitscore=${OPTARG};;
                k) gapopen=${OPTARG};;
                l) gaps=${OPTARG};;
                q) qcovs=${OPTARG};;
		t) num_threads=${OPTARG};;
		u) assignedby=${OPTARG};;
		x) gaf_taxid=${OPTARG};;
        esac
done


ARGS=''
database="${database}"
experimental="${experimental}"
transcript_peps="${transcript_peps}"
trans_peps=$(basename "${transcript_peps}")

#IF STATEMENTS EXIST FOR EACH OPTIONAL BLAST PARAMETER
if [ -n "${E_value}" ]; then ARGS="$ARGS -evalue $E_value"; fi
if [ -n "${max_matches}" ]; then ARGS="$ARGS -max_target_seqs $max_matches"; fi
if [ -n "${num_threads}" ]; then ARGS="$ARGS -num_threads $num_threads"; else ARGS="$ARGS -num_threads 4"; fi
######################################################################################################
##CHOOSE BLAST DATABASE BASED ON WHETHER WE WANT EXPERIMENTAL ONLY OR ALL EVIDENCE CODES
##NEED TO FILTER BLAST DB FOR EXPONLY--RUN BLAST AGAINST ALL THEN FILTER GOA INFO FOR EXPONLY


if [[ "$experimental" = "yes" ]]; then database="$database"'_exponly'; fi
if [[ -z "$experimental" ]]; then database="$database"'_exponly'; fi
name="$database"
database='agbase_database/'"$database"'.fa'
Dbase="$name"'.fa'


##MAKE BLAST INDEX
#MAKE SOME SORT OF IF EXISTS $NAME.PQR (WHATEVER THE EXTENSIONS ARE FOR BLAST DB FILES) THEN SKIP MAKEBLASTDB STEP
#MAYBE MAKE ANOTHER OPTION TO PROVIDE PRE-MADE DATABASE INSTEAD?

makeblastdb -in agbase_database/$Dbase -dbtype prot -parse_seqids -out $name

##RUN BLASTP
##DO WE NEED MORE OPTIONS?
blastp -query $transcript_peps -db $name -out $out.asn -outfmt 11 $ARGS

##Blast_formatter will format stand-alone searches performed with an earlier version of a database if both the search and formatting databases are prepared so that fetching by sequence ID is possible. To enable fetching by sequence ID use the –parse_seqids flag when running makeblastdb, 
blast_formatter -archive $out.asn -out $out.html -outfmt 1 -html
#ALSO SEE WHERE WE CAN PULL PROTEIN NAME FROM--GOA COLUMN 10
#NEED TO FIGURE OUT HOW TO ADD STAXIDS--BLAST CAN'T PULL IT DIRECTLY BECAUSE IT ISN'T IN THE DATABASES
blast_formatter -archive $out.asn -out $out.tsv -outfmt '6 qseqid qstart qend sseqid sstart send evalue pident qcovs ppos gapopen gaps bitscore score'

##WANT THESE
##1. html - shows the alignments
##2. sliminput - this is the summary file of accessions and GO:IDs
##3. tsv - the Blast table (post filtering)
##4. GAF-like output (pull slim input from here)

#################################################################################################################
##ADD FILTERING BASED ON QCOVS (& %% ID & PPOS & BITSCORE & GAPS & GAPOPEN)
if [ -z "${perc_ID}" ]; then perc_ID="0"; fi
if [ -z "${qcovs}" ]; then qcovs="0"; fi
if [ -z "${perc_pos}" ]; then perc_pos="0"; fi
if [ -z "${bitscore}" ]; then bitscore="0"; fi
if [ -z "${gaps}" ]; then gaps="1000"; fi
if [ -z "${gapopen}" ]; then gapopen="100"; fi
awk -v x=$percID -v y=$qcovs -v z=$perc_pos -v w=$bitscore -v v=$gaps -v u=$gapopen '{ if(($8 > x) && ($9 > y) && ($10 > z) && ($13 > w) && ($12 < v) && ($11 < u)) { print }}' $out.tsv > tmp.tsv

##CALCULATE EXTRA COLUMNS AND ADD THEM TO OUTPUT
awk 'BEGIN { OFS = "\t" } {print $1, $3-$2, $2, $3, $4, $6-$5, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}' tmp.tsv > tmp2.tsv

##APPEND HEADER LINE TO TSV
echo -e "Query_ID\tQuery_length\tQuery_start\tQuery_end\tSubject_ID\tSubject_length\tSubject_start\tSubject_end\tE_value\tPercent_ID\tQuery_coverage\tPercent_positive_ID\tGap_openings\tTotal_gaps\tBitscore\tRaw_score" | cat - tmp2.tsv > temp && mv temp $out.tsv

##################################################################################################################
##PULL BLAST IDS FROM BLAST OUTPUT TSV

##THIS LINE IS CLOSE BUT STILL DOESN'T HANDLE HEADER QUITE RIGHT--MAY WORK RIGHT NOW
tail --lines=+2 $out.tsv | awk -F "\t" '{print $1, $5}' > "blstmp.txt"

awk 'BEGIN {OFS = "\t"} {sub(/_.*/, "", $2); print $1, $2}'  blstmp.txt > blastids.txt

if [ ! -d ./splitgoa ]; then mkdir "splitgoa"; fi

##SPLIT GOA DATABASE INTO SEVERAL TEMP FILES BASED ON THE NUMBER OF ENTRIES
if [[ "$experimental" = "no" ]]; then splitB.pl  "go_info/gene_association.goa_uniprot" "splitgoa"; else splitB.pl  "go_info/gene_association_exponly.goa_uniprot" "splitgoa"; fi

##PULL SUBSET OF GOA LINES THAT MATCHED BLAST RESULTS INTO GOA_ENTRIES
cyverse_blast2GO.pl "blastids.txt" "splitgoa"


#OUTGAF VARIABLES COUNT FROM 1 TO CORRESPOND TO THE GAF FILE SPEC
#THESE WILL ALWAYS BE THE SAME AND CAN BE DECLARED OUTSIDE THE LOOP

outgaf1="user_input_db"
outgaf15="user"
outgaf13="taxon:0000"
outgaf14=$(date +"%Y%m%d")
outgaf6="GO_REF:0000024"
outgaf7="ECO:0000247"
outgaf12="protein"
outgaf4=""
outgaf11=""
outgaf17=""
prefix="UniprotKB:"

if [ -n "${gaf_db}" ]; then outgaf1="$gaf_db"; fi
if [ -n "${assignedby}" ]; then outgaf15="$assignedby"; fi
if [ -n "${gaf_taxid}" ]; then outgaf13="taxon:""$gaf_taxid"; fi

#PULLING COLUMNS FROM BLASTIDS AND GOA_ENTRIES PRINTING TO NEW COMBINED FILE GOCOMBO; PULL INFO FROM GOCOMBO AND DECLARED VARIABLES ABOVE TO MAKE GAF OUTPUT
awk 'BEGIN {FS = "\t"}{OFS = "\t"} FNR==NR{a[$2]=$1;next}{ print a[$2], $0}' blastids.txt goa_entries.txt > gocombo_tmp.txt
awk  -v a="$outgaf1" -v b="$outgaf15" -v c="$outgaf13" -v d="$outgaf14" -v e="$outgaf6" -v f="$outgaf7" -v g="$outgaf12" -v h="$outgaf4" -v i="$outgaf11" -v j="$outgaf17" -v k="$prefix" 'BEGIN {FS = "\t"}{OFS = "\t"}{print a,$1,$1,h,$6,e,f,(k$3),$10,$1,i,g,c,d,b,$18,j}' gocombo_tmp.txt > goanna_gaf.txt

##APPEND HEADER TO GAF OUTPUT
sed -i '1 i\!gaf-version: 2.0' goanna_gaf.txt

##PULL COLUMNS FOR GO SLIM FILE
awk 'BEGIN {FS ="\t"}{OFS = "\t"} {print $2,$5,$9}' goanna_gaf.txt > slim_input.txt



##REMOVE FILES THAT ARE NO LONGER NECESSARY
#rm goa_entries.txt
#rm -r splitgoa
#rm gocombo_tmp.txt
#rm blastids.txt
#rm blstmp.txt
#rm tmp.tsv
#rm tmp2.tsv
#rm *.phr
#rm *.pin
#rm *.pog
#rm *.psd
#rm *.psi
#rm *.psq


