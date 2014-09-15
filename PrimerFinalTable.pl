#!/usr/bin/perl -w

use strict;

my ( $pIN, $strictTable, $out, $fa ) = @ARGV;
my %is;
my %state;
my %info;
open OUT,">$out";
print OUT "#ID\tSTATE\tALIGN_INFO\tLEFT_Primer\tRIGHT_Primer\tLEFT_TM\tRIGHT_TM\tLEFT_GC\tRIGHT_GC\tLEFT_ANY_TH\tRIGHT_ANY_TH\tLEFT_END_TH\tRIGHT_END_TH\tLEFT_HAIRPIN_TH\tRIGHT_HAIRPIN_TH\tLEFT_MISPRIMING\tRIGHT_MISPRIMING\tLEFT_POSITION\tRIGHT_POSITION\tPRODUCT_LENGHT(Mut)\tPRODUCT_LENGHT(Ref)\tSEQ\n";

open ( FA, $fa );
my %seq;
my ( $tit, $s );
$s = "";
while ( <FA> ){
	chomp;
	if ( /^>/ ){	
		$seq{$tit} = $s if ( $s ne "" );
		$_ =~ />(.*)/;
		$tit = $1;
		$s = "";
	}else{
		$s = "$s"."$_";
	}
}
close FA;
$seq{$tit} = $s;

open ( IN, $pIN );
while ( <IN> ){
	chomp;
	next if ( /^#/ );
	my @tmp = split /\s+/;
	$tmp[0] =~ s/^read//;
	$is{$tmp[0]} = "$tmp[2]";
	$state{$tmp[0]} = "$tmp[3]";
	$info{$tmp[0]} = "$tmp[1]";
}
close IN;

open ( IN, $strictTable );
my ( $id, $rawid );
my $i = 0;
my ($left_pos,$right_pos);
my $seq;
while ( <IN> ){
    chomp;
	if(/SEQUENCE_ID=(\S+)/){
		$rawid = $1;
		$id = $rawid;
		$i = 0;
	}elsif ( /^SEQUENCE_TEMPLATE=(.*)/ ){
		$seq = $1;
	}elsif ( /PRIMER_LEFT_NUM_RETURNED=0/ ){
        print OUT "$rawid\tMissing\tnull\tnull\tnull\tnull\tnull\tnull\tnull\tnull\tnull\tnull\tnull\tnull\tnull\tnull\t$seq{$rawid}\n";
	}elsif(/PRIMER_LEFT.*SEQUENCE=(\w+)/){
		
		$i ++;
		$id = "$rawid\_$i";
		print OUT "$id\t$state{$id}\t$info{$id}\t$1\t";
	}elsif(/PRIMER_RIGHT.*SEQUENCE=(\w+)/){	
        print OUT "$1\t";
	}elsif(/PRIMER_LEFT.*\d+=(\d+),(\d+)/){
		$left_pos=$1;
	}elsif(/PRIMER_RIGHT.*\d+=(\d+),(\d+)/){
        $right_pos=$1;
	}elsif(/PRIMER_LEFT.*TM=(\S+)/){
        print OUT "$1\t";
	}elsif(/PRIMER_RIGHT.*TM=(\S+)/){
        print OUT "$1\t";
	}elsif(/PRIMER_LEFT.*GC_PERCENT=(\S+)/){
        print OUT "$1\t";
    }elsif(/PRIMER_RIGHT.*GC_PERCENT=(\S+)/){
        print OUT "$1\t";
	}elsif (/PRIMER_LEFT_.*_SELF_ANY_TH=(\S+)/){
		print OUT "$1\t";
	}elsif (/PRIMER_RIGHT_.*_SELF_ANY_TH=(\S+)/){
		print OUT "$1\t";
	}elsif (/PRIMER_LEFT_.*_SELF_END_TH=(\S+)/){
        print OUT "$1\t";
    }elsif (/PRIMER_RIGHT_.*_SELF_END_TH=(\S+)/){
        print OUT "$1\t";
	}elsif (/PRIMER_LEFT_.*_HAIRPIN_TH=(\S+)/){
        print OUT "$1\t";
    }elsif (/PRIMER_RIGHT_.*_HAIRPIN_TH=(\S+)/){
        print OUT "$1\t";
	}elsif (/PRIMER_LEFT_.*_LIBRARY_MISPRIMING=(\S+)\,/){
        print OUT "$1\t";
    }elsif (/PRIMER_RIGHT_.*_LIBRARY_MISPRIMING=(\S+)\,/){
        print OUT "$1\t";
	}elsif(/PRIMER_PAIR_.*_PRODUCT_SIZE=(\d+)/){
	    print OUT "$left_pos\t$right_pos\t$1\t$is{$id}\t$seq{$rawid}\n";
    }
	

}
close IN;





