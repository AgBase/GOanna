#!/usr/bin/perl

# splitA.pl Part of ExView1 Pipeline, Roger Barthelson

$pslfile = $ARGV[0];
#$outroot = $ARGV[1];
$outroot = splitgoa;

if ( !defined $pslfile || !defined $outroot ) {
    die "Usage:  $0  pslfile  outfile\n";
}
open( PFIL, "$pslfile" ) or die "Cannot open pslfile\n";

$ct = 0;
$fl = 1;

$out = join "", $outroot, "/temp.", $fl;

open( OFIL, ">$out" ) or die "Cannot open outfile\n";

while ( $line = <PFIL> ) {
    @cole = split( /\t/, $line );
    if ( $ct < 19999999 ) {
        print OFIL "$line";
        $ct++;
        next;
    }
    if ( $ct == 19999999 ) {
        $read = $cole[1];
        print OFIL "$line";
        $ct++;
        next;
    }
    elsif ( $ct > 19999999 and $cole[1] eq $read ) {
        print OFIL "$line";
        $ct++;
        next;
    }
    else {
        $ct = 0;
        close(OFIL);
        $fl++;
        $out = join "", $outroot, "/temp.", $fl;
        open( OFIL, ">>$out" ) or die "Cannot open outfile\n";

        print OFIL "$line";
        next;
    }
}

close(PFIL);
close(OFIL);
