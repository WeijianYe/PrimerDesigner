#!/usr/bin/perl
use strict;
#Date: 20140801
#Author: Weijian Ye (BGI) & Siyang Liu (BGI&KU) based on a draft script from Liya Lin (BGI)
=head1 Description

This tool aims to offer an easy way to generate primer sequences for PCR using Primer 3 and check the specificity of the candidate primers. 
Mode1: Only the chromosome, start and end position are needed. 
Mode2: Provide your own fa sequence file with [] braces the variant - check the online Primer3 for details of the conventions.
Besides, you should specify the directory of the third tools.

=head1 Usage

  perl PrimerDesigner.pl [opts] input_for_primer3[o] primer_result_prefix[o]

  --region	input region file (or fa);
  --fa		input fasta file (or region);
  --ref <s>	reference version (hg18/[hg19]);
  --lt <int>	leftmost position of target [200];
  --len <int>	minimal length of target, starting from the leftmost [200];
  #sy 20140914: use --lt and --len when input is bed instead of fasta. target refers to the piece of DNA sequence that we aim to amplify, which not only contains the variant but also may contain a short stretch of the flanking region. Provided a bed encoding the sequence, --lt and --len defines the coordinates of the targetted region within this sequence. For example, if bed encodes 100bp insertion with 100bp upstream and 100bp downstream sequence, usually we use --lt 50 --len 150.
  --range <s>	possible range of the product size [200-600];
  --opt <f>	optimal TM [57.0];
  --max <f>	maximal TM [61.0];
  --min <f>	minimal TM [53.0];
  --optS <int> optimal primer size [23];
  --maxS <int> maximal primer size [27];
  --minS <int> minimal primer size [18];
  --diff <f>	maximal TM difference [5.0];
  --mag <int>	maximal GC content [60];
  --mig <int>	minimal GC content [40];
  --rep   the fa file of mispriming library
  --bwa        the directory of bwa binary ( newest bwa version, have to contain "mem" alignment pattern )
  --primer      the directory of primer3 binary
  --refD        you have to specify the directory contained your reference files ( must be splited by chromosomes and have a wholegenome fa named hg19.fa or hg18.fa )
  --thirdD      the directory of third tools( fastaDeal.pl, ConvertPrimer2Fq.pl, PrimerStat.pl, PrimerFinalTable.pl )
  --out     outdir
  --h		print this usage

  region.file:Name, chrN (do not use Chr/chromosome, please.) start end, tab separated.
  fasta file contains head and sequence for each region.
  Options 'lt' and 'len' are used to confine the region available for primer design. e.g. lt=200, len=200 means primers are not allowed to locate on the target starting at 200 leftmost and 200bp in size. More:Google it.

=head1 examples
  perl PrimerDesigner.pl -fa region.fa -refD refDir --bwa bwa --primer primer3 --thirdD demoDir --out outdir primerDesign.I primerDesign.O
  perl PrimerDesigner.pl -region region.txt --refD refDir --bwa bwa --primer primer3 --thirdD demoDir --out outdir primerDesign.I primerDesign.O

=cut

use Getopt::Long;

##get options from command line into variables and set default values
my ($region,$fasta,$ref,$lt,$len,$range,$opt_tm,$max_tm,$min_tm,$opt_size,$max_size,$min_size,$max_gc,$min_gc,$diff_tm,$bwa,$refD,$thirdD,$primer,$out,$rep,$help);

GetOptions(
	"region:s"=>\$region,
	"fa:s"=>\$fasta,
	"ref:s"=>\$ref,
	"lt:i"=>\$lt,
	"len:i"=>\$len,
	"range:s"=>\$range,
	"opt:s"=>\$opt_tm,
	"max:s"=>\$max_tm,
	"min:s"=>\$min_tm,
	"optS:s"=>\$opt_size,
    "maxS:s"=>\$max_size,
    "minS:s"=>\$min_size,
	"diff:s"=>\$diff_tm,
	"mag:i"=>\$max_gc,
	"mig:i"=>\$min_gc,
	"bwa:s"=>\$bwa,
	"refD:s"=>\$refD,
	"thirdD:s"=>\$thirdD,
	"primer:s"=>\$primer,
	"out:s"=>\$out,
	"rep:s"=>\$rep,
	"h"=>\$help,
);
$opt_size ||= 23;
$min_size ||= 18;
$max_size ||= 27;
$lt ||= 200;
$len ||= 200;
$range ||= "200-600";
$range =~ /(\d+)\-\d+/;
my $varS = $1;
$opt_tm ||= "57.0";
$max_tm ||= "61.0";
$min_tm ||= "53.0";
$diff_tm ||= "5.0";
$max_gc ||= 60;
$min_gc ||= 40;
$ref ||="hg19";
$bwa ||="/home/siyang/Bin/software_pip/bwa-0.7.5a/bwa";
$refD ||="/home/siyang/USER/yeweijian/Ref/hg19";
$thirdD ||="/faststorage/home/siyang/USER/yeweijian/Project/DanishPanGenome/2014July_Task/20140723_PrimerDesign/bin";
$primer ||="/faststorage/home/siyang/BACKUP/Bin/software_pip/primer3-2.3.6/src/primer3_core";
my $conf = `dirname $primer`;
$conf =~ s/\n//;
$conf = "$conf\/primer3_config\/";
$rep ||="/faststorage/home/siyang/USER/yeweijian/Project/DanishPanGenome/2014July_Task/20140723_PrimerDesign/bin/HumanMisPrimerRep.fa";
#$ref ||="/share/backup/jiawl/useful/human/hg19/splitbychr_hg19";

die `pod2text $0` if (@ARGV < 2 || $help);

my $input_fa;
if($fasta){
    $input_fa = $fasta;
}elsif($region){
	open POS,$region || die $!;
	my $files = "";
	while(<POS>){
		chomp;
		my @seq_pos = split /\t/,$_;
		my $ref_chr;
		#$seq_pos[0]=~s/-/_/g;
		$seq_pos[1]=~s/\W//g;
		$seq_pos[2]=~s/\D//g;
		$seq_pos[3]=~s/\D//g;
		if($ref eq "hg19"){
			$ref_chr = "$refD/$seq_pos[1].fa";
		}
		elsif($ref eq "hg18"){
			$ref_chr = "$refD/$seq_pos[1].fa";
		}
		else{
			die "Please choose one correct HG version\n";
		}
		`perl $thirdD/fastaDeal.pl -sub $seq_pos[2]-$seq_pos[3] $ref_chr > $out/$seq_pos[0]-$seq_pos[1]_$seq_pos[2]_$seq_pos[3].tmp`;
		$files .= " $out/$seq_pos[0]-$seq_pos[1]_$seq_pos[2]_$seq_pos[3].tmp";
	}
	close POS;

	`cat $out/*tmp > $out/all.input.fa`; ##cat all individual fa into one sinlge file.
	$input_fa = "$out/all.input.fa";
}

my ($target,$opt);
$target = "SEQUENCE_TARGET=$lt,$len\nPRIMER_PRODUCT_SIZE_RANGE=$range";
$opt = "PRIMER_OPT_SIZE=$opt_size\nPRIMER_MIN_SIZE=$min_size\nPRIMER_MAX_SIZE=$max_size\nPRIMER_OPT_TM=$opt_tm\nPRIMER_MIN_TM=$min_tm\nPRIMER_MAX_TM=$max_tm\nPRIMER_PAIR_MAX_DIFF_TM=$diff_tm\nPRIMER_MAX_GC=$max_gc\nPRIMER_MIN_GC=$min_gc\nPRIMER_FIRST_BASE_INDEX=1\nPRIMER_THERMODYNAMIC_PARAMETERS_PATH=$conf\nPRIMER_MISPRIMING_LIBRARY=$rep";

open IN,$input_fa || die $!;
open OUT,">",$ARGV[-2] || die $!;
my ($tmp,$id,$seq);
chomp ($tmp =<IN>);
$id = $1 if($tmp =~ /^>(.*)/);
while(<IN>){
	chomp;
	if(/^>(.*)/){
		

		if ( $seq =~ /</ and $seq =~ /\[/ ){
			( $lt,$len ) = &FindInclude ( $seq );
			my $r1 = $len - 100; # extend 100bp
			my $r2 = $len + 100;
			#$range = "$r1-$r2";
			$target = "SEQUENCE_TARGET=$lt,$len\nPRIMER_PRODUCT_SIZE_RANGE=$range";
			my @reg = &FindExclude ( $seq );
            my $exclude = "";
            for ( my $i = 0; $i < @reg; $i += 2 ){
                my $s = $reg[$i];
                my $n = $i + 1;
                my $l = $reg[$n];
                if ( $exclude eq "" ){
                    $exclude = "SEQUENCE_EXCLUDED_REGION=$s,$l";
                }else{
                    $exclude = "$exclude\nSEQUENCE_EXCLUDED_REGION=$s,$l";
                }
            }
			$seq =~ s/<//g;
            $seq =~ s/>//g;
			$seq =~ s/\[//g;
            $seq =~ s/\]//g;
			print OUT "SEQUENCE_ID=$id\nSEQUENCE_TEMPLATE=$seq\n$target\n$exclude\n$opt\n=\n";
		}elsif ( $seq =~ /\[/ ){
			( $lt,$len ) = &FindInclude ( $seq );
			my $r1 = $len - 100; # extend 100bp
            my $r2 = $len + 100;
			#$range = "$r1-$r2";
            $target = "SEQUENCE_TARGET=$lt,$len\nPRIMER_PRODUCT_SIZE_RANGE=$range";
			$seq =~ s/\[//g;
            $seq =~ s/\]//g;
			print OUT "SEQUENCE_ID=$id\nSEQUENCE_TEMPLATE=$seq\n$target\n$opt\n=\n";
		}elsif ( $seq =~ /</ ){
			my @reg = &FindExclude ( $seq );
			my $exclude = "";
			for ( my $i = 0; $i < @reg; $i += 2 ){
				my $s = $reg[$i];
				my $n = $i + 1;
				my $l = $reg[$n];
				if ( $exclude eq "" ){
					$exclude = "SEQUENCE_EXCLUDED_REGION=$s,$l";
				}else{
					$exclude = "$exclude\nSEQUENCE_EXCLUDED_REGION=$s,$l";
				}
			}
			$seq =~ s/<//g;
			$seq =~ s/>//g;
			print OUT "SEQUENCE_ID=$id\nSEQUENCE_TEMPLATE=$seq\n$target\n$exclude\n$opt\n=\n";
		}else{
			print OUT "SEQUENCE_ID=$id\nSEQUENCE_TEMPLATE=$seq\n$target\n$opt\n=\n";
		}
		$_ =~ />(.*)/;
		$id = $1;
		$seq = "";
		last;
	}
	else{
		#my $str = $1 if(/(\w+)/);
		my $str = $_;
		$seq = "$seq"."$str";
	}
}
### print the options.
while(<IN>){
	chomp;
	if(/^>(.*)/){
		if ( $seq =~ /</ and $seq =~ /\[/ ){
            ( $lt,$len ) = &FindInclude ( $seq );
            my $r1 = $len - 100; # extend 100bp
            my $r2 = $len + 100;
            #$range = "$r1-$r2";
            $target = "SEQUENCE_TARGET=$lt,$len\nPRIMER_PRODUCT_SIZE_RANGE=$range";
            my @reg = &FindExclude ( $seq );
            my $exclude = "";
            for ( my $i = 0; $i < @reg; $i += 2 ){
                my $s = $reg[$i];
                my $n = $i + 1;
                my $l = $reg[$n];
                if ( $exclude eq "" ){
                    $exclude = "SEQUENCE_EXCLUDED_REGION=$s,$l";
                }else{
                    $exclude = "$exclude\nSEQUENCE_EXCLUDED_REGION=$s,$l";
                }
            }
            $seq =~ s/<//g;
            $seq =~ s/>//g;
            $seq =~ s/\[//g;
            $seq =~ s/\]//g;
            print OUT "SEQUENCE_ID=$id\nSEQUENCE_TEMPLATE=$seq\n$target\n$exclude\n$opt\n=\n";
        }elsif ( $seq =~ /\[/ ){
            ( $lt,$len ) = &FindInclude ( $seq );
            my $r1 = $len - 100; # extend 100bp
            my $r2 = $len + 100;
            #$range = "$r1-$r2";
            $target = "SEQUENCE_TARGET=$lt,$len\nPRIMER_PRODUCT_SIZE_RANGE=$range";
            $seq =~ s/\[//g;
            $seq =~ s/\]//g;
            print OUT "SEQUENCE_ID=$id\nSEQUENCE_TEMPLATE=$seq\n$target\n$opt\n=\n";
		}elsif ( $seq =~ /</ ){
            my @reg = &FindExclude ( $seq );
            my $exclude = "";
            for ( my $i = 0; $i < @reg; $i += 2 ){
                my $s = $reg[$i];
                my $n = $i + 1;
                my $l = $reg[$n];
                if ( $exclude eq "" ){
                    $exclude = "SEQUENCE_EXCLUDED_REGION=$s,$l";
                }else{
                    $exclude = "$exclude\nSEQUENCE_EXCLUDED_REGION=$s,$l";
                }
	        }
			$seq =~ s/<//g;
            $seq =~ s/>//g;
            print OUT "SEQUENCE_ID=$id\nSEQUENCE_TEMPLATE=$seq\n$target\n$exclude\n$opt\n=\n";
		}else{
			print OUT "SEQUENCE_ID=$id\nSEQUENCE_TEMPLATE=$seq\n$target\n=\n";
		}
		$_ =~ />(.*)/;
		$id = $1;
		$seq = "";
	}else {
		#my $str = $1 if(/(\w+)/);
		my $str = $_;
		$seq = "$seq"."$str";
	}
}
if ( $seq =~ /</ and $seq =~ /\[/ ){
	( $lt,$len ) = &FindInclude ( $seq );
    my $r1 = $len - 100; # extend 100bp
    my $r2 = $len + 100;
    #$range = "$r1-$r2";
    $target = "SEQUENCE_TARGET=$lt,$len\nPRIMER_PRODUCT_SIZE_RANGE=$range";
    my @reg = &FindExclude ( $seq );
    my $exclude = "";
    for ( my $i = 0; $i < @reg; $i += 2 ){
    	my $s = $reg[$i];
        my $n = $i + 1;
        my $l = $reg[$n];
        if ( $exclude eq "" ){
        	$exclude = "SEQUENCE_EXCLUDED_REGION=$s,$l";
        }else{
            $exclude = "$exclude\nSEQUENCE_EXCLUDED_REGION=$s,$l";
        }
    }
    $seq =~ s/<//g;
    $seq =~ s/>//g;
    $seq =~ s/\[//g;
    $seq =~ s/\]//g;
    print OUT "SEQUENCE_ID=$id\nSEQUENCE_TEMPLATE=$seq\n$target\n$exclude\n$opt\n=\n";
}elsif ( $seq =~ /\[/ ){
    ( $lt,$len ) = &FindInclude ( $seq );
	my $r1 = $len - 100; # extend 100bp
	my $r2 = $len + 100;
    #$range = "$r1-$r2";
    $target = "SEQUENCE_TARGET=$lt,$len\nPRIMER_PRODUCT_SIZE_RANGE=$range";
    $seq =~ s/\[//g;
    $seq =~ s/\]//g;
    print OUT "SEQUENCE_ID=$id\nSEQUENCE_TEMPLATE=$seq\n$target\n$opt\n=\n";
}elsif ( $seq =~ /</ ){
	my @reg = &FindExclude ( $seq );
    my $exclude = "";
    for ( my $i = 0; $i < @reg; $i += 2 ){
    	my $s = $reg[$i];
        my $n = $i + 1;
        my $l = $reg[$n];
        if ( $exclude eq "" ){
        	$exclude = "SEQUENCE_EXCLUDED_REGION=$s,$l";
        }else{
            $exclude = "$exclude\nSEQUENCE_EXCLUDED_REGION=$s,$l";
        }
	}
    $seq =~ s/<//g;
    $seq =~ s/>//g;
    print OUT "SEQUENCE_ID=$id\nSEQUENCE_TEMPLATE=$seq\n$target\n$exclude\n$opt\n=\n";
}else{
	print OUT "SEQUENCE_ID=$id\nSEQUENCE_TEMPLATE=$seq\n$target\n=\n";
}
close IN;
close OUT;
print "========== Running Primer3 ==========\n";
`$primer -format_output <$ARGV[-2] >$ARGV[-1].format.out.txt`;
`$primer -strict_tags <$ARGV[-2] >$ARGV[-1].strict.out.txt`;
#`/faststorage/home/siyang/USER/yeweijian/Project/DanishPanGenome/2014July_Task/20140723_PrimerDesign/bin/Primer3_get_result.pl $ARGV[-1].strict.out.txt $ARGV[-1].strict.out.table.txt`;

print "========== Fq preparation ==========\n";
`perl $thirdD/ConvertPrimer2Fq.pl $ARGV[-1].strict.out.txt $ARGV[-1] $ARGV[-1].Mis.txt`;

print "========== Running BWA MEM ==========\n";
`$bwa mem -R '\@RG\tID:Primer\tPL:Illumina\tPI:500\tSM:Pilot' -T 0 -C -M -a $refD/$ref.fa $ARGV[-1].1.fq $ARGV[-1].2.fq | samtools view -Sb - > $ARGV[-1].bam` ;

print "========== STAT ==========\n";
`perl $thirdD/PrimerStat.pl $ARGV[-1].bam $ARGV[-1].alignInfo.txt $ARGV[-1].Failed.alignInfo.txt $ARGV[-1].Major.stat $varS`;
`perl $thirdD/PrimerFinalTable.pl $ARGV[-1].alignInfo.txt $ARGV[-1].strict.out.txt $ARGV[-1].Final.table $input_fa`;

`rm -f $out/*tmp` if ( $region );
print "========== All done ==========\n";

sub FindExclude{
    my ( $seq ) = @_;
    my @tmp = split //, $seq;
    my $seqOrder = 0;
    my @outRegion;
    my $start;
	my $flag1 = 0; #mark the '<' start position
	my $flag2 = 0;
    for ( my $i = 0; $i < @tmp; $i ++ ){
        if ( $tmp[$i] eq '<' ){
			if ( $flag1 ){
				$flag2 ++;
				next;		
			}
            $start = $seqOrder + 1;
            push ( @outRegion, $start );
			$flag1 ++;
        }elsif ( $tmp[$i] eq '>' ){
			if ( $flag2 ){
				$flag2 --;
				next;
			}
            my $len = $seqOrder - $start + 1;
            push ( @outRegion, $len );
			$flag1 = 0;
			$flag2 = 0; #reset
        }elsif ( $tmp[$i] ne '['  and $tmp[$i] ne ']' ){
            $seqOrder ++;
        }
    }
    return @outRegion;
}

sub FindInclude{
    my ( $seq ) = @_;
    my @tmp = split //, $seq;
    my $seqOrder = 0;
    my ( $start, $len );
    for ( my $i = 0; $i < @tmp; $i ++ ){
        if ( $tmp[$i] eq '[' ){
            $start = $seqOrder + 1;
        }elsif ( $tmp[$i] eq ']' ){
            $len = $seqOrder - $start + 1;
			last;
        }elsif ( $tmp[$i] ne '<' and $tmp[$i] ne '>' ){
            $seqOrder ++;
        }
    }
    return ( $start, $len );
}
