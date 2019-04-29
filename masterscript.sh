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
               # d) out_fmt=${OPTARG};;
                f) max_matches=${OPTARG};;
                e) E_value=${OPTARG};;
                g) percID=${OPTARG};;
               # h) help=${OPTARG};; ADD USAGE STATEMENT
		m) perc_pos=${OPTARG};;
		o) out=${OPTARG};;
		d) bitscore=${OPTARG};;
                k) gap_Open=${OPTARG};;
                l) gap_Extend=${OPTARG};;
                q) qcovs=${OPTARG};;
		t) num_threads=${OPTARG};;
		u) rawscore=${OPTARG};;
        esac
done


ARGS=''
database="${database}"
experimental="${experimental}"
transcript_peps="${transcript_peps}"
trans_peps=$(basename "${transcript_peps}")

#IF STATEMENTS EXIST FOR EACH OPTIONAL PARAMETER
if [ -n "${perc_pos}" ]; then ARGS="$ARGS -ppos $perc_pos"; fi
if [ -n "${E_value}" ]; then ARGS="$ARGS -eval $E_value"; fi
if [ -n "${percID}" ]; then ARGS="$ARGS -pident $percID"; fi
#if [ -n "${out}" ]; then ARGS="$ARGS -out $out"; fi
if [ -n "${gap_Open}" ]; then ARGS="$ARGS -gapopen $gap_Open"; fi
if [ -n "${gap_Extend}" ]; then ARGS="$ARGS -gapextend $gap_Extend"; fi
if [ -n "${max_matches}" ]; then ARGS="$ARGS -max_target_seqs $max_matches"; fi
if [ -n "${bitscore}" ]; then ARGS="$ARGS -bitscore $bitscore"; fi
if [ -n "${rawscore}" ]; then ARGS="$ARGS -score $rawscore"; fi
if [ -n "${num_threads}" ]; then ARGS="$ARGS -num_threads $num_threads"; else ARGS="$ARGS -num_threads 4"; fi

######################################################################################################
##CHOOSE BLAST DATABASE BASED ON WHETHER WE WANT EXPERIMENTAL ONLY OR ALL EVIDENCE CODES
##NEED TO FILTER BLAST DB FOR EXPONLY--RUN BLAST AGAINST ALL THEN FILTER GOA INFO FOR EXPONLY

if [[ "$experimental" = "yes" ]]; then database="$database"'_exponly'; fi
name="$database"
database='agbase_database/'"$database"'.fa'
Dbase="$name"'.fa'


##MAKE BLAST INDEX
#MAKE SOME SORT OF IF EXISTS $NAME.PQR (WHATEVER THE EXTENSIONS ARE FOR BLAST DB FILES) THEN SKIP MAKEBLASTDB STEP
#MAYBE MAKE ANOTHER OPTION TO PROVIDE PRE-MADE DATABASE INSTEAD?

#makeblastdb -in agbase_database/$Dbase -dbtype prot -parse_seqids -out $name

##RUN BLASTP
##DO WE NEED MORE OPTIONS?
#blastp -query $trans_peps -db $name -out $out.asn -outfmt 11 $ARGS

##Blast_formatter will format stand-alone searches performed with an earlier version of a database if both the search and formatting databases are prepared so that fetching by sequence ID is possible. To enable fetching by sequence ID use the –parse_seqids flag when running makeblastdb, 
#blast_formatter -archive $out.asn -out $out.html -outfmt 1 -html
#ALSO SEE WHERE WE CAN PULL PROTEIN NAME FROM--GOA COLUMN 10
#NEED TO FIGURE OUT HOW TO ADD STAXIDS--BLAST CAN'T PULL IT DIRECTLY BECAUSE IT ISN'T IN THE DATABASES
#blast_formatter -archive $out.asn -out $out.tsv -outfmt '6 qseqid qstart qend sseqid sstart send evalue pident qcovs ppos gapopen gaps bitscore score'

##WANT THESE
##1. html - shows the alignments
##2. sliminput - this is the summary file of accessions and GO:IDs
##3. tsv - the Blast table (post filtering)
##4. GAF-like output (pull slim input from here)

#################################################################################################################
##ADD FILTERING BASED ON QCOVS (& % ID & PPOS & BITSCORE)
#awk -v x=$percID -v y=$qcovs -v z=$perc_pos -v w=$bitscore '{ if(($8 > x) && ($9 > y) && ($10 > z) && ($13 > w)) { print }}' $out.tsv > tmp.tsv

##CALCULATE EXTRA COLUMNS AND ADD THEM TO OUTPUT
#awk 'BEGIN { OFS = "\t" } {print $1, $3-$2, $2, $3, $4, $6-$5, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}' tmp.tsv > tmp2.tsv

##APPEND HEADER LINE TO TSV
#echo -e "Query_ID\tQuery_length\tQuery_start\tQuery_end\tSubject_ID\tSubject_length\tSubject_start\tSubject_end\tE_value\tPercent_ID\tQuery_coverage\tPercent_positive_ID\tGap_openings\tTotal_gaps\tBitscore\tRaw_score" | cat - tmp2.tsv > temp && mv temp $out.tsv

##################################################################################################################
#script to transfer GOA info from gene_association.goa_uniprot file for the agbase uniprot entries -- maintained in the Data Store for public access

##PULL BLAST IDS FROM BLAST OUTPUT TSV

##THIS LINE IS CLOSE BUT STILL DOESN'T HANDLE HEADER QUITE RIGHT--MAY WORK RIGHT NOW
tail --lines=+2 $out.tsv | awk -F "\t" '{print $1, $5}' > "blstmp.txt"

#sed -i 's/_.*//'  blastids.txt ##NEEDS TO APPLY TO COL 5 BUT NOT COL 1
awk 'BEGIN {OFS = "\t"} {sub(/_.*/, "", $2); print $1, $2}'  blstmp.txt > blastids.txt

if [ ! -d ./splitgoa ]; then mkdir "splitgoa"; fi

##SPLIT GOA DATABASE INTO SEVERAL TEMP FILES BASED ON THE NUMBER OF ENTRIES
#if [[ "$experimental" = "yes" ]]; then splitB.pl  "go_info/gene_association_exponly.goa_uniprot" "splitgoa"; else splitB.pl  "go_info/gene_association.goa_uniprot" "splitgoa"; fi

##MERGE GOA INFO INTO BLAST RESULTS
#cyverse_blast2GO.pl "blastids.txt" "splitgoa"


##PULL NECESSARY COLUMNS FOR OUTPUT FILE
#NEED THESE COLUMNS:
#DB (of object--default to user_input)
#DB OBJECT = QUERY ID
#DB OBJECT SYMBOL =REPEAT DB OBJECT IF CAN'T GET SYMBOL
#QUALIFIER = MUST BE EMPTY BUT WE NEED TO HAVE THE EMPTY COLUMN
#GO ID = PULL FROM GOA GAF COLUMN 5
#DB:ref = ALWAYS USE GO_REF:0000024
#EVIDENCE CODE = ALWAYS USE ECO:0000247
#WITH OR FROM = UNIPROT:SUBJECT ID
#ASPECT = PULL FROM GOA COLUMN 9
#DB OBJECT NAME = EVERYTHING FROM QUERY FASTA HEADER
#DB OBJECT SYNONYM = EMPTY BUT NEED TO HAVE THE EMPTY COLUMN
#DB OBJECT TYPE = PROTEIN
#TAXON = TAXON ID OF QUERY (TAXON:####)--WILL BE USER INPUT?
#DATE = YYYYMMDD
#ASSIGNED BY = USER INPUT (DEFAULT TO 'USER INPUT')--OR MAKE IT REQUIRED
#ANNOTATION EXTENSION = PULL FROM GOA FILE COLUMN 16
#GENE PRODUCT FORM ID = EMPTY BUT NEED TO HAVE THE EMPTY COLUMN

##PULL COLUMNS FOR GO SLIM FILE

##REMOVE FILES THAT ARE NO LONGER NECESSARY
#rm gene_association.goa_uniprot
#rm goa_entries.txt

#echo "I think I'm done!"


