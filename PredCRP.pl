#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Cwd 'abs_path';
use File::Basename;
my $program = abs_path($0);
my ($filename, $PredCRPdir) = fileparse($program);
my $svmscale = $PredCRPdir.'/svm-scale';
my $svmpredict = $PredCRPdir.'/svm-predict';
my $model = $PredCRPdir.'/PredCRP_model';
my $tmp_location = $PredCRPdir.'/tmp';
my $result_location = $PredCRPdir.'/predict_result';
my $scale = $tmp_location."/crp_strong_169_12features_scl";
my $input = '';
my $help;
sub Usage(){
	print
"Usage: perl $program [Option]
Option:
-i FILE: input CRP binding sites information
-svmscale pathname: set svm-scale executable path and name
-svmpredict pathname: set svm-predict executable path and name
-model pathname: set PredCRP_model path and name
-h, -help
";
}
GetOptions(
	'i=s'	=>\$input,
	'svmscale=s'	=>\$svmscale,
	'svmpredict=s'	=>\$svmpredict,
	'model=s'	=>\$model,
	'h|help'	=>\$help,
);

if($help){
	&Usage();
	exit;
}
if(!$input){
	print STDERR "No input file\n";
	&Usage();
	exit;
}
if($input !~ /\.csv/){
	print STDERR "The format of input file must be csv\n";
	exit;
}
`perl -p -i -e "s/\r/\n/g" $input`;
my $svm_12features = fileparse($input);
$svm_12features =~ s/\.csv//g;
#my $input_tmp = $tmp_location.'/'.$svm_12features."_tmp";
my $location_svm_12features = $tmp_location.'/'.$svm_12features."_svm";
my $location_svm_12features_scl = $tmp_location.'/'.$svm_12features."_svm_scl";
my $location_predict_result_tmp = $result_location.'/'.$svm_12features."_predict_tmp";
my $location_predict_result = $result_location.'/'.$svm_12features."_predict";
my $location_predict_final = $result_location.'/'.$svm_12features."_PredictResult.csv";
my ($sequence, $tss, $unit, $class_tmp, $lengthSeq, $promoter_flag, $bubble_flag, $ar_overlap, $ca_overlap);
my ($aacg, $catt, $gaac, $gagc, $tgcg, $ttac, $ttat, $tttt);
my $header = "CRPBS,Distance of Center Position of CRPBS to TSS,Transcription Unit,Regulatory Role";
open SVM,">",$location_svm_12features;
open FILE,"<",$input;
my $line=<FILE>;
chomp $line;
	if($line !~ /$header/){
		#print "$line\n$header\n";
		print STDERR "No header in input file\n";
		exit;
	}
my $num=0;
while($line=<FILE>){
	chomp $line;
	my @ele = split(/,/,$line);
	$ele[0] =~ tr/A-Z/a-z/;
	$sequence = $ele[0];
	$tss = $ele[1];
	$unit = $ele[2];
	$class_tmp = $ele[3];
	$lengthSeq = length($sequence);
	if($lengthSeq != 42){
		print STDERR "input sequence error, 10bp + CRPBS(22bp) + 10bp\n";
		exit;
	}
	feature_extraction($sequence, $tss, $unit, $class_tmp);	
}
`$svmscale -r $scale $location_svm_12features > $location_svm_12features_scl`;
`$svmpredict -b 1 $location_svm_12features_scl $model $location_predict_result_tmp`;
open Pred,">",$location_predict_result;
print Pred ",Pred role,Pred probability\n";
open Predtmp,"<",$location_predict_result_tmp;
$line=<Predtmp>;
while($line=<Predtmp>){
	chomp $line;
	my @ele = split(/ /,$line);
	if($ele[0] == 0){
		print Pred ",activation,$ele[1]\n";
	}else{
		print Pred ",repression,$ele[1]\n";
	}	
}
`paste $input $location_predict_result > $location_predict_final`;
`rm $location_predict_result_tmp`;
`rm $location_predict_result`;
close FILE;
close SVM;
close Predtmp;
close Pred;

sub feature_extraction{
	my ($sequence, $tss, $unit, $class_tmp) = @_;
	$aacg=0;
	$catt=0;
	$gaac=0;
	$gagc=0;
	$tgcg=0;
	$ttac=0;
	$ttat=0;
	$tttt=0;
	my @split_seq = split(//,$sequence);
	foreach my $i (@split_seq){
		if($i !~ m/[atcg]/){
			print STDERR "input sequence error, only ATCG\n";	
			exit;		
		}
	}
	for(my $i=0; $i<@split_seq-3;$i++){
		my $four_motif = substr($sequence,$i,4);
		if($four_motif eq 'aacg'){
			$aacg++;
		}elsif($four_motif eq 'catt'){
			$catt++;
		}elsif($four_motif eq 'gaac'){
			$gaac++;
		}elsif($four_motif eq 'gagc'){
			$gagc++;
		}elsif($four_motif eq 'tgcg'){
			$tgcg++;
		}elsif($four_motif eq 'ttac'){
			$ttac++;
		}elsif($four_motif eq 'ttat'){
			$ttat++;
		}elsif($four_motif eq 'tttt'){
			$tttt++;
		}
	}	
	my $class = 0;
	if($class_tmp eq '-'){
		$class = 1;
	}elsif($class_tmp eq '+'){
		$class = 0;
	}else{
		$class = 0;
	}
	#=======8 motifs =====================
	my @motif;
	$promoter_flag = 0;
	$bubble_flag = 0;
	$ar_overlap = 0;
	$ca_overlap = 0;
	#======= -35 ~ -10 ==============
	if(($tss >= -35) and ($tss <= -10)){
		$promoter_flag = 1;
	}elsif((($tss-11) > -35) and (($tss-11) < -10)){
		$promoter_flag = 1;
	}elsif((($tss+11) > -35) and (($tss+11) < -10)){
		$promoter_flag = 1;
	}
	#===== -10 ~ 2 ===================
	if(($tss >= -10) and ($tss <= 2)){
		$bubble_flag = 1;
	}elsif((($tss-11) > -10) and (($tss-11) < 2)){
		$bubble_flag = 1;
	}elsif((($tss+11) > -10) and (($tss+11) < 2)){
		$bubble_flag = 1;
	} 
	#===== -60 ~ 60 ==================
	if(($tss >= -60) and ($tss <= 60)){
		my $right = abs($tss+60);
		my $left = abs(60-$tss);
		if($right>11){
			$right = 11;
		}
		if($left>11){
			$left = 11;
		}
		$ar_overlap = abs($right+$left);	
	}elsif((($tss-11) > -60) and (($tss-11) < 60)){
		$ar_overlap = abs(60-($tss-11));	
	}elsif((($tss+11) > -60) and (($tss+11) < 60)){
		$ar_overlap = abs(($tss+11)+60);
	}
	#===== -95 ~ -35 ===================
	if(($tss >= -95) and ($tss <= -35)){
		my $right = abs($tss+95);
		my $left = abs(-35-$tss);
		if($right > 11){
			$right = 11;
		}
		if($left > 11){
			$left = 11;
		}
		$ca_overlap = abs($left + $right);		
	}elsif((($tss-11) > -95) and (($tss-11) < -35)){
		$ca_overlap = abs(-35-($tss-11));
	}elsif((($tss+11) > -95) and (($tss+11) < -35)){
		$ca_overlap = abs(($tss+11)+95);
	}
	print SVM "$class 1:$aacg 2:$catt 3:$gaac 4:$gagc 5:$tgcg 6:$ttac 7:$ttat 8:$tttt 9:$promoter_flag 10:$bubble_flag 11:$ar_overlap 12:$ca_overlap\n";
	#die;
}


