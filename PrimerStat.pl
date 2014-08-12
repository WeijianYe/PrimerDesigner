#!/usr/bin/perl -w

my ( $bam, $info, $err, $stat, $r  ) = @ARGV; # $r acceptable region, approximately equal to variant size
my $is1 = 3 * $r; #acceptable insertsize, 3 times of variant size
my $is2 = $r / 3;
#die "$is1\t$is2";
open ( IN, "samtools view -X $bam|" );
open ( O, ">$info" ); #pass alignment info
print O "#ID\talignINFO\tInsertSize\tState\n";
open ( ERR, ">$err" );
open ( STAT, ">$stat" ); # only check the major primer situation
print STAT "#Tol\(Major\)\tPass\tRate\n";
my %data;
my $Tol;
my $Pass;
while ( <IN> ){
	chomp;
	my @tmp = split /\s+/;
	$Tol ++ if ( $tmp[0] =~ /1$/ );
	$data{$tmp[0]} ++;
}
close IN;
$Tol = $Tol / 2; # double output
open ( IN, "samtools view -X $bam|" );
while ( <IN> ){
    chomp;
    my @tmp = split /\s+/;
	my $l1 = $_;
	if ( $data{$tmp[0]} == 2 ){
		chomp( my $l2=<IN> );
		my @tmp2 = split /\s+/, $l2;
		$tmp[0] =~ /read(.*)\_(\d+)\_\d+\_\d+/;
		#$tmp[0] =~ /read(.*)\:\d+\-\d+\_\d+/; test
		my $chr = $1;
		my $p = $2;
		if ( $tmp[2] ne $tmp2[2] or $tmp[1] =~ /U/ or $tmp[1] =~ /u/ or $chr ne $tmp[2] ){
			print ERR "$l1\n$l2\n";
			print O "$tmp[0]\t$tmp[2]\,$tmp[3]\,$tmp[7]\,$tmp[8]\t$tmp[8]\tReferenceAlignmentProblem\n";
			next;
		}
		if ( $tmp[8] > $is1 or $tmp[8] < $is2 ){
			print ERR "$l1\n$l2\n";
			print O "$tmp[0]\t$tmp[2]\,$tmp[3]\,$tmp[7]\,$tmp[8]\t$tmp[8]\tReferenceAlignmentProblem\n";
            next;
		}
		#if ( $tmp[3] > ($p + $r) or $tmp[3] < ($p - $r) ){
		#	print ERR "$l1\n$l2\n";
		#	print O "$tmp[0]\t$tmp[2]\,$tmp[3]\,$tmp[7]\,$tmp[8]\t$tmp[8]\tReferenceAlignmentProblem\n";
        #    next;
		#}
		$Pass ++ if ( $tmp[0] =~ /1$/ );
		print O "$tmp[0]\t$tmp[2]\,$tmp[3]\,$tmp[7]\,$tmp[8]\t$tmp[8]\tPASS\n";
	}else{
		print ERR "$_\n";
		print O "$tmp[0]\tNA\tNA\tMultiAlignment\n";
	}
}
close IN;
my $rate = $Pass / $Tol;
print STAT "$Tol\t$Pass\t$rate\n";
