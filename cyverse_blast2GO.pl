#!/usr/bin/perl
use Data::Dumper;
# ###############################################################
# script: uniprot2goa_gaf.pl
# desc: get GOA GAF records for blast hits
#
# command line arguments: blast_hits filename (blast_hits.txt)
#                         GAF records directory name (splitgoa)
#
# ############################################################
  
my ($blast_hit_filename,$gaf_directory)=@ARGV; 

my $hit_colnum=1;
my $gaf_obj_colnum=1;
my %blast_hits=();
my %EXP_EVD_CODES = ('EXP' => 1, 'IDA' => 1, 'IPI' => 1, 'IMP' => 1, 'IGI' => 1, 
     'IEP' => 1,'HTP' => 1,'HMP' => 1,'HGI' => 1,'HDA' => 1,'HEP' =>1);
	 
# ############################
# read in the blast hits 
# ############################
open(my $blast_fh,"<",$blast_hit_filename) or die "could not open blast hit file\n";
 while (my $line = <$blast_fh>) {
     chomp($line);
      $line =~ s/\r?\n//; 
    my @flds =split(/\t/, $line,-1);
    $blast_hits{$flds[$hit_colnum]}=1;
} # while input blast hits
close ($blast_fh);
my $blast_hit_cnt = keys %blast_hits;

open(my $out_fh, ">","goa_entries.txt");

 # ####################################
 # read every GAF file in tmp dir
 # searching for blast hits
 # if found print out experiment GAF
 #
 # first, create array of 
 # files in directory argument that
 # begins with string 'temp'
 # ###################################
opendir(DIR, $gaf_directory)  or die "Could not open $gaf_directory\n";
my @files = grep(/^temp/,readdir(DIR));
closedir(DIR);
  # ########################################
  #  foreach filename in array 
  # make sure a regular file with records
  #
  # pull in uniprotkb db, null qualifier,
  # experimental only that maps to EXPER evd  
  # ########################################
  my $found_cnt=0;
foreach $filename (@files) {
if ((!-f $filename) && (-s $filename)) { next; }  
 open(my $gaf_fh, "<","$gaf_directory/$filename") || die "can't open file";
 while (my $line = <$gaf_fh>) {
     chomp($line);
      $line =~ s/\r?\n//; 
    my @flds =split(/\t/, $line,-1);
	  # must be uniprotkb not rnacentral or intact
	if ( lc $flds[0] ne 'uniprotkb') {  next; }
	
	          # qualifier must be blank
  if ($flds[3] ne '') {  next; }

      # evd must be experimental
#  if (not exists $EXP_EVD_CODES{$flds[6]}) { next; }  
 
	if (exists $blast_hits{$flds[$gaf_obj_colnum]}) {
	     print $out_fh $line . "\n";
		 $found_cnt++;
    } # blast hit match in GAF 
  } # while input GAF records
  close ($gaf_fh);
} # foreach file in GAF directory
closedir(DIR);
close $out_fh;
 
 print "\n\nend of script. $blast_hit_filename contained $blast_hit_cnt unique records, these hits matched to $found_cnt GAF records from $gaf_directory\n";
