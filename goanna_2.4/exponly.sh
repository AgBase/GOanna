#! /bin/bash

while getopts i:o: option
do
        case "${option}"
        in

                i) uniprot=${OPTARG};;
                o) exponly=${OPTARG};;
        esac
done


awk '$6 == "EXP"; $6=="IDA"; $6=="IPI"; $6=="IMP"; $6=="IGI"; $6=="IEP"; $6=="HTP"; $6=="HMP"; $6=="HGI"; $6=="HDA"; $6=="HEP" {print}' $uniprot > $exponly

timestamp=$(date +"%F %T")

echo -e \
"!gaf-version: 2.1
!
!This file contains all GO annotations and gene product information for proteins in the UniProt KnowledgeBase (UniProtKB),
!ComplexPortal protein complexes, and RNAcentral identifiers.
!
!Generated: $timestamp
!GO-version: http://purl.obolibrary.org/obo/go/releases/2019-06-28/extensions/go-plus.owl
!
!Experimental only evidence codes were subset from the gene_association.goa_uniprot file using the following command:
!awk '"$6" == "EXP"; "$6"=="IDA"; "$6"=="IPI"; "$6"=="IMP"; "$6"=="IGI"; "$6"=="IEP"; "$6"=="HTP"; "$6"=="HMP"; "$6"=="HGI"; "$6"=="HDA"; "$6"=="HEP"'" | cat - $exponly > temp && mv temp $exponly

