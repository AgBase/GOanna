#!/bin/bash

#######################################################################################################
##SET UP OPTIONS FOR MAKEBLASTDB AND BLASTP

while getopts a:b:c:d:e:f:g:hk:l:m:n:o:pq:r:s:t:u:v:x:y:z option
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
		m) perc_pos=${OPTARG};;
		o) out=${OPTARG};;
		s) bitscore=${OPTARG};;
                k) gapopen=${OPTARG};;
                l) gaps=${OPTARG};;
                q) qcovs=${OPTARG};;
		t) num_threads=${OPTARG};;
		u) assignedby=${OPTARG};;
		x) gaf_taxid=${OPTARG};;
		h) help=true ;;
		p) pdef=true ;;
        esac
done
#####################################################################################################
if [[ "$help" = "true" ]] ; then
  echo "Options:
    -a Blast database basename ('arthropod', 'bacteria', 'bird', 'crustacean', 'fish', 'fungi', 'human', 'insecta',
       'invertebrates', 'mammals', 'nematode', 'plants', 'rodents' 'uniprot_sprot', 'uniprot_trembl' or 'vertebrates')
    -c peptide fasta filename
    -o Blast output file basename
    [-b transfer GO with experimental evidence only ('yes' or 'no')]
    [-d database of query ID. Default: 'user_input_db']
    [-f Number of aligned sequences to keep. Default: 500]
    [-g Blast percent identity above which match should be kept. Default: keep all matches.]
    [-h help]
    [-m Blast percent positive identity above which match should be kept. Default: keep all matches.]
    [-s bitscore above which match should be kept. Default: keep all matches.]
    [-k Maximum number of gap openings allowed for match to be kept.Default: 100]
    [-l Maximum number of total gaps allowed for match to be kept. Default: 1000]
    [-q Minimum query coverage per subject for match to be kept. Default: keep all matches]
    [-t Number of threads.  Default: 4]
    [-u 'Assigned by' field of your GAF output file. Default: 'user']
    [-x Taxon ID of the peptides you are blasting. Default: 'taxon:0000']
    [-p parse_deflines. Parse query and subject bar delimited sequence identifiers]" 
  exit 0
fi
#####################################################################################################

ARGS=''
database="${database}"
experimental="${experimental}"
transcript_peps="${transcript_peps}"
trans_peps=$(basename "${transcript_peps}")

#IF STATEMENTS EXIST FOR EACH OPTIONAL BLAST PARAMETER
if [ -n "${E_value}" ]; then ARGS="$ARGS -evalue $E_value"; fi
if [ -n "${max_matches}" ]; then ARGS="$ARGS -max_target_seqs $max_matches"; fi
if [ -n "${num_threads}" ]; then ARGS="$ARGS -num_threads $num_threads"; fi
if [[ "$pdef" = "true" ]]; then ARGS="$ARGS -parse_deflines"; fi
######################################################################################################

##CHOOSE BLAST DATABASE BASED ON WHETHER WE WANT EXPERIMENTAL ONLY OR ALL EVIDENCE CODES
##DEFAULTS TO EXPERIMENTAL
if [[ "$experimental" = "yes" ]]; then database="$database"'_exponly'; fi
if [[ -z "$experimental" ]]; then database="$database"'_exponly'; fi
name="$database"
database='agbase_database/'"$database"'.fa'
Dbase="$name"'.fa'


##MAKE BLAST INDEX
test -f "/agbase_database/$Dbase" && makeblastdb -in /agbase_database/$Dbase -dbtype prot -parse_seqids -out $name
test -f "agbase_database/$Dbase" && makeblastdb -in agbase_database/$Dbase -dbtype prot -parse_seqids -out $name
    
##RUN BLASTP
blastp  -query $transcript_peps -db $name -out $out.asn -outfmt 11 $ARGS


##MAKE BLAST OUTPUT FORMATS 1 AND 6
blast_formatter -archive $out.asn -out $out.html -outfmt 0 -html
blast_formatter -archive $out.asn -out $out.tsv -outfmt '6 qseqid qstart qend sseqid sstart send evalue pident qcovs ppos gapopen gaps bitscore score'
#################################################################################################################

##FILTER BALST OUTPUT 6 (OPTIONALLY) BY %ID, QUERY COVERAGE, % POSITIVE ID, BITSCORE, TOTAL GAPS, GAP OPENINGS
if [ -z "${perc_ID}" ]; then perc_ID="0"; fi
if [ -z "${qcovs}" ]; then qcovs="0"; fi
if [ -z "${perc_pos}" ]; then perc_pos="0"; fi
if [ -z "${bitscore}" ]; then bitscore="0"; fi
if [ -z "${gaps}" ]; then gaps="1000"; fi
if [ -z "${gapopen}" ]; then gapopen="100"; fi
awk -v x=$percID -v y=$qcovs -v z=$perc_pos -v w=$bitscore -v v=$gaps -v u=$gapopen '{ if(($8 > x) && ($9 > y) && ($10 > z) && ($13 > w) && ($12 < v) && ($11 < u)) { print }}' $out.tsv > tmp.tsv

##CALCULATE QUERY AND SUBJECT LENGTH COLUMNS AND ADD THEM TO OUTPUT 6
awk 'BEGIN { OFS = "\t" } {print $1, $3-$2, $2, $3, $4, $6-$5, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}' tmp.tsv > tmp2.tsv

##APPEND HEADER LINE TO OUTPUT 6
echo -e "Query_ID\tQuery_length\tQuery_start\tQuery_end\tSubject_ID\tSubject_length\tSubject_start\tSubject_end\tE_value\tPercent_ID\tQuery_coverage\tPercent_positive_ID\tGap_openings\tTotal_gaps\tBitscore\tRaw_score" | cat - tmp2.tsv > temp && mv temp $out.tsv

##################################################################################################################
##PULL COLUMNS 1 AND 5 (QUERY ID AND SUBJECT ID) FOR ALL LINES EXCEPT HEADER
tail --lines=+2 $out.tsv | awk -F "\t" '{print $1, $5}' > "blstmp.txt"

##REMOVE THE _ AND EVERYTHING AFTER FROM THE SUBJECT ID SO THAT IT WILL MATCH THE GOA FILE
awk 'BEGIN {OFS = "\t"} {sub(/_.*/, "", $2); print $1, $2}'  blstmp.txt > blastids.txt

##MAKE KOABAS ANNOATATE INPUT FILE
#awk 'BEGIN {OFS = "\t"} {print $2}' blastids.txt | uniq > $out'_KOBAS_annotate_input.txt'

##SPLIT GOA DATABASE INTO SEVERAL TEMP FILES BASED ON THE NUMBER OF ENTRIES
if [ ! -d ./splitgoa ]; then mkdir "splitgoa"; fi

if [[ "$experimental" = "no" ]]
then
    test -f /go_info/gene_association.goa_uniprot && splitB.pl  "/go_info/gene_association.goa_uniprot" "splitgoa"
    test -f ./go_info/gene_association.goa_uniprot && splitB.pl  "/go_info/gene_association.goa_uniprot" "splitgoa"
elif [[ "$experimental" != "no" ]]
then
    test -f /go_info/gene_association_exponly.goa_uniprot && splitB.pl  "/go_info/gene_association_exponly.goa_uniprot" "splitgoa"
    test -f ./go_info/gene_association_exponly.goa_uniprot && splitB.pl  "./go_info/gene_association_exponly.goa_uniprot" "splitgoa"
fi

##PULL SUBSET OF GOA LINES THAT MATCHED BLAST RESULTS INTO GOA_ENTRIES.TXT
cyverse_blast2GO.pl "blastids.txt" "splitgoa"

#OUTGAF VARIABLES COUNT FROM 1 TO CORRESPOND TO THE GAF FILE SPEC
#THESE WILL ALWAYS BE THE SAME AND CAN BE DECLARED OUTSIDE THE AWK STATEMENT

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

#THESE ARE OPTIONALLY USER-SPECIFIED DEFAULTS IN LIST ABOVE
if [ -n "${gaf_db}" ]; then outgaf1="$gaf_db"; fi
if [ -n "${assignedby}" ]; then outgaf15="$assignedby"; fi
if [ -n "${gaf_taxid}" ]; then outgaf13="taxon:""$gaf_taxid"; fi

#PULLING COLUMNS FROM BLASTIDS.TXT AND GOA_ENTRIES.TXT AND  PRINTING TO NEW COMBINED FILE GOCOMBO; PULL INFO FROM GOCOMBO_TMP.TXT  AND DECLARED VARIABLES ABOVE TO MAKE GAF OUTPUT
awk 'BEGIN {FS = "\t"}{OFS = "\t"} FNR==NR{a[$2]=$1;next}{ print a[$2], $0}' blastids.txt goa_entries.txt > gocombo_tmp.txt
awk  -v a="$outgaf1" -v b="$outgaf15" -v c="$outgaf13" -v d="$outgaf14" -v e="$outgaf6" -v f="$outgaf7" -v g="$outgaf12" -v h="$outgaf4" -v i="$outgaf11" -v j="$outgaf17" -v k="$prefix" 'BEGIN {FS = "\t"}{OFS = "\t"}{print a,$1,$1,h,$6,e,f,(k$3),$10,$1,i,g,c,d,b,$18,j}' gocombo_tmp.txt > $out'_goanna_gaf.tsv'

##APPEND HEADER TO GAF OUTPUT
sed -i '1 i\!gaf-version: 2.0' $out'_goanna_gaf.tsv'

##PULL COLUMNS FOR GO SLIM FILE
awk 'BEGIN {FS ="\t"}{OFS = "\t"} {print $2,$5,$9}' $out'_goanna_gaf.tsv' > $out'_slim_input.txt'



##REMOVE FILES THAT ARE NO LONGER NECESSARY
if [ -s $out'_slim_input.txt' ]
then
    rm goa_entries.txt
    rm -r splitgoa
    rm gocombo_tmp.txt
    rm blstmp.txt
    rm blastids.txt 
    rm tmp.tsv
    rm tmp2.tsv
    rm *.phr
    rm *.pin
    rm *.pog
    rm *.psd
    rm *.psi
    rm *.psq
fi


