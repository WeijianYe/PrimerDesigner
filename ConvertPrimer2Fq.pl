#!/usr/bin/perl -w

my ( $in, $outprefix, $mis ) = @ARGV;
open ( IN, $in );
open ( F1, ">$outprefix.1.fq" );
open ( F2, ">$outprefix.2.fq" );
open ( M, ">$mis" );
my %data;
my $id;
my $i;
while ( <IN> ){
	chomp;
	if ( /^SEQUENCE_ID=(\S+)/ ){
		$id = $1;	
		$i = 0;
	}elsif ( /PRIMER_LEFT_NUM_RETURNED=0/ ){
		print M "$id\n";
	}elsif ( /PRIMER_.*_SEQUENCE=(\w+)/ ){
		$i ++;
		my $fq1 = $1;
		my $fq1Q = $1;
		$fq1Q =~ s/./g/g;
		chomp( my $tmp=<IN> );
		$tmp =~ /PRIMER_.*_SEQUENCE=(\w+)/;
		my $fq2 = $1;
		my $fq2Q = $1;
        $fq2Q =~ s/./g/g;
		print F1 "\@read$id\_$i\/1\n$fq1\n\+\n$fq1Q\n";
		print F2 "\@read$id\_$i\/2\n$fq2\n\+\n$fq2Q\n";
	}
}
close IN;
