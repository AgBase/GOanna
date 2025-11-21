#!/bin/bash

#######################################################################################################
##SET UP OPTIONS FOR MAKEBLASTDB AND BLASTP

while getopts a:A:b:c:d:e:f:g:G:hk:l:m:n:o:pq:r:s:t:u:v:x:y:z option
do
        case "${option}"
        in

                a) database=${OPTARG};;
		A) agbase_db=${OPTARG};;
                b) experimental=${OPTARG};;
                c) transcript_peps=${OPTARG};;
                d) gaf_db=${OPTARG};;
                f) max_matches=${OPTARG};;
                e) E_value=${OPTARG};;
                g) percID=${OPTARG};;
		G) go_info=${OPTARG};;
                m) perc_pos=${OPTARG};;
                o) out=${OPTARG};;
                s) bitscore=${OPTARG};;
                k) gapopen=${OPTARG};;
                l) gaps=${OPTARG};;
                q) qcovs=${OPTARG};;
                r) ratioqs=${OPTARG};;
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
       'invertebrates', 'mammals', 'nematode', 'plants', 'rodents' 'uniprot_sprot', 'uniprot_trembl', 'vertebrates'
	or 'viruses')
    -c peptide fasta filename
    -o Blast output file basename
    [-A name of the agbase_database file to use (e.g. invertebrates_exponly.fa)]
    [-b transfer GO with experimental evidence only ('yes' or 'no'). Default = 'yes'.]
    [-d database of query ID. If your entry contains spaces either substitute and underscore (_) or, 
        to preserve the space, use quotes around your entry. Default: 'user_input_db']
    [-e Expect value (E) for saving hits. Default is 10.]
    [-f Number of aligned sequences to keep. Default: 5]
    [-g Blast percent identity above which match should be kept. Default: keep all matches.]
    [-G specify 'go_info/gene_association.goa_uniprot.gz' or 'go_info/gene_association_exponly.goa_uniprot.gz'
    [-h help]
    [-m Blast percent positive identity above which match should be kept. Default: keep all matches.]
    [-s bitscore above which match should be kept. Default: keep all matches.]
    [-k Maximum number of gap openings allowed for match to be kept.Default: 100]
    [-l Maximum number of total gaps allowed for match to be kept. Default: 1000]
    [-q Minimum query coverage per subject for match to be kept. Default: keep all matches]
    [-r Ratio of query length to subject length. Lengths should be comparable for match to be kept. Default: less than 1.2 so difference of up to 20% can be tolerated]
    [-t Number of threads.  Default: 8]
    [-u 'Assigned by' field of your GAF output file. If your entry contains spaces (eg. firstname lastname) 
        either substitute and underscore (_) or, to preserve the space, use quotes around your entry (eg. "firstname lastname")
        Default: 'user']
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
agdb="${agbase_db}"
goinfo="${go_info}"
echo -e "agdb is:" $agdb
echo -e "goinfo is:" $goinfo

#IF STATEMENTS EXIST FOR EACH OPTIONAL BLAST PARAMETER
if [ -n "${E_value}" ]; then ARGS="$ARGS -evalue $E_value"; fi
if [ -n "${max_matches}" ]; then ARGS="$ARGS -max_target_seqs $max_matches"; else ARGS="$ARGS -max_target_seqs 5"; fi
if [ -n "${num_threads}" ]; then ARGS="$ARGS -num_threads $num_threads"; else ARGS="$ARGS -num_threads 8"; fi
if [[ "$pdef" = "true" ]]; then ARGS="$ARGS -parse_deflines"; fi
######################################################################################################

##CHOOSE BLAST DATABASE BASED ON WHETHER WE WANT EXPERIMENTAL ONLY OR ALL EVIDENCE CODES
##DEFAULTS TO EXPERIMENTAL
if [[ -z "${agbase_db}" ]]
then
	if [[ "$experimental" = "yes" ]]; then database="$database"'_exponly'; fi
	if [[ -z "$experimental" ]]; then database="$database"'_exponly'; fi
	name="$database"
	database='agbase_database/'"$database"'.fa'
	Dbase="$name"'.fa'
	if [ -f "/agbase_database/$Dbase.gz" ]
	then
		echo "/agbase_database/$Dbase.gz"
		gunzip /agbase_database/$Dbase.gz
		makeblastdb -in /agbase_database/$Dbase -dbtype prot -parse_seqids -out $name
		echo 'database made' $(date)
	elif [ -f "agbase_database/$Dbase.gz" ]
	then
		echo "agbase_database/$Dbase.gz"
		gunzip "agbase_database/$Dbase.gz"
		makeblastdb -in agbase_database/$Dbase -dbtype prot -parse_seqids -out $name
		echo 'database made' $(date)
	elif [ -f "/agbase_database/$Dbase" ]
	then
		echo "/agbase_database/$Dbase"
		ls -l /agbase_database
		makeblastdb -in /agbase_database/$Dbase -dbtype prot -parse_seqids -out $name
		echo 'database made' $(date)
	elif [ -f "agbase_database/$Dbase" ]
	then
		echo "agbase_database/$Dbase"
		makeblastdb -in agbase_database/$Dbase -dbtype prot -parse_seqids -out $name
		echo 'database made' $(date)
	else
		echo "There is no agbase_database file."
	fi
else
	echo "You have specified an agbase_database file."
	name=$(basename "${agbase_db}")
	Dbase="${agbase_db}"
	makeblastdb -in $Dbase -dbtype prot -parse_seqids -out $name
fi

#GIVES ERRORS IF USERS TRY TO SELECT DBS THAT DON'T EXIST (BUT LOGICALLY SHOULD)
if [ $Dbase = "viruses_exponly.fa" ]; then echo "There are too few experimentally annotated viruses to perform this search. Please try all annotations instead (-b no)."; exit; fi
if [ $Dbase = "uniprot_sprot.fa" ]; then echo "This will search all of uniprot_sprot. To obtain high quality annotations please try experimental annotations only (-b yes)."; exit; fi
if [ $Dbase = "uniprot_trembl.fa" ]; then echo "This will search all of uniprot_trembl. To obtain high quality annotations please try experimental annotations only ( -b yes)."; exit; fi

#RUN BLASTP
blastp  -query $transcript_peps -db $name -out $out.asn -outfmt 11  $ARGS

echo 'blastp run' $(date)

##MAKE BLAST OUTPUT FORMATS 1 AND 6
blast_formatter -archive $out.asn -out $out.html -outfmt 0 -html
echo 'html formatter run' $(date)

blast_formatter -archive $out.asn -out $out.tsv -outfmt '6 qseqid qstart qend sseqid sstart send evalue pident qcovs ppos gapopen gaps bitscore score qlen slen'
echo 'tsv formatter run' $(date)
#################################################################################################################

##FILTER BLAST OUTPUT 6 (OPTIONALLY) BY %ID, QUERY COVERAGE, % POSITIVE ID, BITSCORE, TOTAL GAPS, GAP OPENINGS, RATIO OF QUERY LENGTH TO SUBJECT LENGTH
if [ -z "${perc_ID}" ]; then perc_ID="0"; fi
if [ -z "${qcovs}" ]; then qcovs="0"; fi
if [ -z "${perc_pos}" ]; then perc_pos="0"; fi
if [ -z "${bitscore}" ]; then bitscore="0"; fi
if [ -z "${gaps}" ]; then gaps="1000"; fi
if [ -z "${gapopen}" ]; then gapopen="100"; fi
if [ -z "${ratioqs}" ]; then ratioqs="1.2"; fi

awk -v x=$percID -v y=$qcovs -v z=$perc_pos -v w=$bitscore -v v=$gaps -v u=$gapopen -v r=$ratioqs '{ if(($8 > x) && ($9 > y) && ($10 > z) && ($13 > w) && ($12 < v) && ($11 < u) && ($15/$16 <= r)) { print }}' $out.tsv > tmp.tsv

##CALCULATE QUERY AND SUBJECT LENGTH COLUMNS AND ADD THEM TO OUTPUT 6
awk 'BEGIN { OFS = "\t" } {print $1, $3-$2, $2, $3, $4, $6-$5, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}' tmp.tsv > tmp2.tsv

##APPEND HEADER LINE TO OUTPUT 6
echo -e "Query_ID\tQuery_length\tQuery_start\tQuery_end\tSubject_ID\tSubject_length\tSubject_start\tSubject_end\tE_value\tPercent_ID\tQuery_coverage\tPercent_positive_ID\tGap_openings\tTotal_gaps\tBitscore\tRaw_score" | cat - tmp2.tsv > temp && mv temp $out.tsv

##################################################################################################################
##PULL COLUMNS 1 AND 5 (QUERY ID AND SUBJECT ID) FOR ALL LINES EXCEPT HEADER
tail --lines=+2 $out.tsv | awk -F "\t" '{print $1, $5}' > "blstmp.txt"

##REMOVE THE _ AND EVERYTHING AFTER FROM THE SUBJECT ID SO THAT IT WILL MATCH THE GOA FILE
awk 'BEGIN {OFS = "\t"} {sub(/_.*/, "", $2); print $1, $2}'  blstmp.txt > blastids.txt

##SPLIT GOA DATABASE INTO SEVERAL TEMP FILES BASED ON THE NUMBER OF ENTRIES
if [ ! -d splitgoa ]; then mkdir "splitgoa"; fi

if [[ -z "${go_info}" ]]
then
	#UNZIP NECESSARY  GO INFO DATABASES AND RUN SPLITB.PL
	if [[ "$experimental" = "no" ]]
	then
		echo "exponly=no"
		if [ -f /go_info/gene_association.goa_uniprot.gz ]
		then
			echo '/go_info/gene_association.goa_uniprot.gz'
			gunzip /go_info/gene_association.goa_uniprot.gz
			splitB.pl  "/go_info/gene_association.goa_uniprot" "splitgoa"
		elif [ -f /go_info/gene_association.goa_uniprot ]
		then
			echo '/go_info/gene_association.goa_uniprot'
			splitB.pl  "/go_info/gene_association.goa_uniprot" "splitgoa"
		elif [ -f go_info/gene_association.goa_uniprot.gz ]
		then
			echo 'go_info/gene_association.goa_uniprot.gz'
			gunzip go_info/gene_association.goa_uniprot.gz
			splitB.pl  "go_info/gene_association.goa_uniprot" "splitgoa"
		elif [ -f go_info/gene_association.goa_uniprot ]
		then
			echo 'go_info/gene_association.goa_uniprot'
			splitB.pl  "go_info/gene_association.goa_uniprot" "splitgoa"
		else
			echo 'There is no go_info file.'
		fi
	elif [[ "$experimental" != "no" ]]
	then
		echo "exponly=yes"
		if [ -f /go_info/gene_association_exponly.goa_uniprot.gz ]
		then
			echo '/go_info/gene_association_exponly.goa_uniprot.gz'
			gunzip /go_info/gene_association_exponly.goa_uniprot.gz
			splitB.pl  "/go_info/gene_association_exponly.goa_uniprot" "splitgoa"
		elif [ -f /go_info/gene_association_exponly.goa_uniprot ]
		then
			echo '/go_info/gene_association_exponly.goa_uniprot'
			splitB.pl  "/go_info/gene_association_exponly.goa_uniprot" "splitgoa"
		elif [ -f ./go_info/gene_association_exponly.goa_uniprot.gz ]
		then
			echo './go_info/gene_association_exponly.goa_uniprot.gz'
			gunzip ./go_info/gene_association_exponly.goa_uniprot.gz
			splitB.pl  "./go_info/gene_association_exponly.goa_uniprot" "splitgoa"
		elif [ -f ./go_info/gene_association_exponly.goa_uniprot ]
		then
			echo './go_info/gene_association_exponly.goa_uniprot'
			splitB.pl  "./go_info/gene_association_exponly.goa_uniprot" "splitgoa"
		else
			echo 'There is no go_info file.'
		fi
	fi
else
	echo "You have specified a go_info file:"
	echo "${go_info}"
	splitB.pl  "${go_info}" "splitgoa"
fi

##PULL SUBSET OF GOA LINES THAT MATCHED BLAST RESULTS INTO GOA_ENTRIES.TXT
cyverse_blast2GO.pl "blastids.txt" "splitgoa"

sed -i 's/\t/!/g' goa_entries.txt

#OUTGAF VARIABLES COUNT FROM 1 TO CORRESPOND TO THE GAF FILE SPEC
#THESE (IMMEDIATELY BELOW) WILL ALWAYS BE THE SAME AND CAN BE DECLARED OUTSIDE THE AWK STATEMENT

outgaf1="user_input_db"
outgaf15="user"
outgaf13="taxon:0000"
outgaf14=$(date +"%Y%m%d")
outgaf6="GO_REF:0000024"
outgaf7="ISA"
outgaf12="protein"
outgaf11=""
outgaf17=""
prefix="UniprotKB:"
outgaf4=''

#THESE ARE OPTIONALLY USER-SPECIFIED--DEFAULTS IN LIST ABOVE
if [ -n "${gaf_db}" ]; then outgaf1="$gaf_db"; fi
if [ -n "${assignedby}" ]; then outgaf15="$assignedby"; fi
if [ -n "${gaf_taxid}" ]; then outgaf13="taxon:""$gaf_taxid"; fi

#PULLING COLUMNS FROM BLASTIDS.TXT AND GOA_ENTRIES.TXT AND PRINTING TO NEW COMBINED FILE GOCOMBO
#PULL INFO FROM GOCOMBO_TMP.TXT  AND DECLARED VARIABLES ABOVE TO MAKE GAF OUTPUT
touch $out'_goanna_gaf.tsv'
sed -i 's/\t/!/g' blastids.txt
awk 'BEGIN {FS = "!"}{OFS = "!"} FNR==NR{a[$2]=$1;next}{ print a[$2], $0}' blastids.txt goa_entries.txt > gocombo_tmp.txt
cat gocombo_tmp.txt | while IFS="!" read -r xp unip id symbol qual goacc pmid evid empty aspect name sym prot tax date assdb empty2 empty3
do
	if [[ $aspect == P ]];
		then
		outgaf4="involved_in"
		echo -e "$outgaf1\t$xp\t$xp\t$outgaf4\t$goacc\t$outgaf6\t$outgaf7\t$prefix$id\t$aspect\t$xp\t$outgaf11\t$outgaf12\t$outgaf13\t$outgaf14\t$outgaf15\t$empty3\t$outgaf17" >> $out'_goanna_gaf.tsv'
	elif [[ $aspect == F ]];
		then
		outgaf4="enables"
		echo -e "$outgaf1\t$xp\t$xp\t$outgaf4\t$goacc\t$outgaf6\t$outgaf7\t$prefix$id\t$aspect\t$xp\t$outgaf11\t$outgaf12\t$outgaf13\t$outgaf14\t$outgaf15\t$empty3\t$outgaf17" >> $out'_goanna_gaf.tsv'
	elif [[ $aspect == C ]] && grep -q $goacc '/usr/GO0032991_and_children.json';
		then
		outgaf4="part_of"
		echo -e "$outgaf1\t$xp\t$xp\t$outgaf4\t$goacc\t$outgaf6\t$outgaf7\t$prefix$id\t$aspect\t$xp\t$outgaf11\t$outgaf12\t$outgaf13\t$outgaf14\t$outgaf15\t$empty3\t$outgaf17" >> $out'_goanna_gaf.tsv'
	elif [[ $aspect == C ]] && grep -v -q $goacc '/usr/GO0032991_and_children.json';
		then
		outgaf4="located_in"
		echo -e "$outgaf1\t$xp\t$xp\t$outgaf4\t$goacc\t$outgaf6\t$outgaf7\t$prefix$id\t$aspect\t$xp\t$outgaf11\t$outgaf12\t$outgaf13\t$outgaf14\t$outgaf15\t$empty3\t$outgaf17" >> $out'_goanna_gaf.tsv'
	else
		echo $goacc $aspect "ERROR: qualifier not set"
	fi
done
cat $out'_goanna_gaf.tsv' | uniq -u > uniq.gaf
mv uniq.gaf $out'_goanna_gaf.tsv'

#ENFORCE TAXON CONSTRAINTS
wget https://current.geneontology.org/ontology/imports/go-computed-taxon-constraints.obo -O constraints.obo
wget https://current.geneontology.org/ontology/imports/go-taxon-groupings.obo -O unions.obo

lineage=( $(echo "${gaf_taxid}" | taxonkit lineage -t | cut -f 3 | sed  's/;/,/g') )

python /usr/bin/taxon_constrain.py $out'_goanna_gaf.tsv' constraints.obo unions.obo $lineage .

mv new_keep_gaf.tsv $out'_goanna_gaf.tsv'

##APPEND HEADER TO GAF OUTPUT
header="!gaf-version: 2.2\n!date-generated:$(date +'%Y-%m-%d')\n!generated-by: AgBase\n\nDatabase\tDB_Object_ID\tDB_Object_Symbol\tQualifier\tGO_ID\tDB_Reference\tEvidence_Code\tWith_From\tAspect\tDB_Object_Name\tDB_Object_Synonyms\tDB_Object_Type\tTaxon\tDate\tAssigned_By\tAnnotation_Extension\tGene_Product_Form_Id"

printf "$header\n" | cat - $out'_goanna_gaf.tsv' > headergaf.tmp
mv headergaf.tmp $out'_goanna_gaf.tsv'

##REMOVE FILES THAT ARE NO LONGER NECESSARY
if [ -s $out'_goanna_gaf.tsv' ]
then
    rm goa_entries.txt
    rm -r splitgoa
    rm gocombo_tmp.txt
    rm blstmp.txt
    rm blastids.txt
    rm tmp.tsv
    rm tmp2.tsv
    rm $name'.pdb'
    rm $name'.pos'
    rm $name'.pot'
    rm $name'.ptf'
    rm $name'.pto'
    rm $name'.phr'
    rm $name'.pog'
    rm $name'.psq'
    rm $name'.pin'
    rm constraints.obo
    rm unions.obo
fi
