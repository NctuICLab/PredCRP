#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long qw(:config no_ignore_case);
my $program = abs_path($0);
my $location_independent = 'location_independent_feature_extract.pl';
my ($regulonDB,$chooseTF,$evidence,$length,$help);
sub Usage(){
	print STDERR "Usage $program [Option]
	Option:
		-input		[FILE] The BindingSitesSet.txt from RegulonDB
		-TF		[STR]	The interested TF (Ex: CRP)
		-evidence	[No]	Evidence level (0:Weak, 1:Strong, 2:Both)
		-length		[No]	10bp+BindingSites+10bp (Ex: The length of CRP-BS is 42)
		-h		Show the usage
		";
}
GetOptions(
	'input=s'	=>\$regulonDB,
	'TF=s'	=>\$chooseTF,
	'evidence=s'	=>\$evidence,
	'length=s'	=>\$length,
	'h'	=>\$help,
);
if($help || !$regulonDB || !$chooseTF || !$evidence || !$length){
	&Usage();die;
}
if(! -e $regulonDB){
	print STDERR "The regulonDB file is not exist\n";
	&Usage();
	die;
}
if(!$chooseTF){
	print STDERR "Please input the interested TF name\n";
	&Usage();
	die;
}
if(($evidence != 0)and($evidence != 1)and($evidence != 2)){
	print STDERR "Please inpput the correct evidence level\n";
	&Usage();
	die;
}
if(!$length){
	print STDERR "Please input the lengh of TF-BS\n";
	&Usage();
	die;
}
=pod
my $regulonDB = $ARGV[0];
my $chooseTF = $ARGV[1];
my $evidence = $ARGV[2];
my $length = $ARGV[3];
if(@ARGV<3){
	print STDERR "Usage: perl feature_extraction_BindingSites.pl RegulonDB TF Evidence(0:weak, 1:strong, 2:Both Length_BS(ex:CRP 42))\n";
	die;
}
=cut
my $output;
if($evidence == 0){
	$output = $regulonDB."_".$chooseTF."_weak.txt";
}elsif($evidence == 1){
	$output = $regulonDB."_".$chooseTF."_strong.txt";
}elsif($evidence == 2){
	$output = $regulonDB."_".$chooseTF."_all.txt";
}
open OUT,">",$output;
my $no = 0;
my $alltf_num = 0;
my ($bl_overlap,$br_overlap,@b_overlap,@b_count);
my (@up_down_poly,$l_overlap,$r_overlap,@overlap,@count,@dis_poly);
my ($c1l_overlap,$c1r_overlap,@c1_overlap,@c1_count);
my ($c2l_overlap,$c2r_overlap,@c2_overlap,@c2_count);
my ($arl_overlap,$arr_overlap,@ar_overlap,@rl_overlap,@cr_overlap);
my ($rrl_overlap,$rrr_overlap,@rr_overlap,@rr_count);
my ($ral_overlap,$rar_overlap,@ra_overlap,@ra_count);
my ($cal_overlap,$car_overlap,@ca_overlap,@ca_count);
my ($perl_overlap,$perr_overlap,@per_overlap,@per_count);
my ($bs_seq,$len_bs,$tf,$relative_dis,$role,$direction);
my (@ar_count,@rl_count,@cr_count,@a_r_number);
my ($role_class,$dir_class);
my ($is_overlap,$is_b_overlap,$is_c1_overlap,$is_c2_overlap,$is_ar_overlap,$is_rr_overlap);
my ($is_rl_overlap,$is_ra_overlap,$is_ca_overlap,$is_cr_overlap,$is_per_overlap);
#read all TF data
my (@start, @end, @class, @sequence, @tf_name);
my ($start, $end);
open REGDB,"<",$regulonDB;
while(my $line=<REGDB>){
	chomp $line;
	if($line =~ /^#/){
		next;
	}
	my @stemp=split(/\t/,$line);
	my $pass = 1;
	foreach my $i (@stemp){
		if(!$i){
			$pass = 0;
		}elsif(($i eq '+-')or($i eq '?')){
			$pass = 0;
		}elsif($i eq 'Strong'){
			if($evidence == 0){
				$pass = 0;
			}
		}elsif($i eq 'Weak'){
			if($evidence == 1){
				$pass = 0;
			}
		}
	}
	if(!$pass){
		next;
	}
	($tf,$start,$end,$role) = ($stemp[1],$stemp[3],$stemp[4],$stemp[8]);
	if($role eq '-'){
		$role_class = 1;
	}elsif($role eq '+'){
		$role_class = 0;
	}
	$start[$alltf_num] = $start;
	$end[$alltf_num] = $end;
	$tf_name[$alltf_num] = $tf;
	$class[$alltf_num] = $role_class;
	$alltf_num++;
}
close REGDB;





open REGDB,"<",$regulonDB;
while(my $line=<REGDB>){
	chomp $line;
	if($line =~ /^#/){
		next;
	}
	my @stemp=split(/\t/,$line);
	my $pass = 1;
	foreach my $i (@stemp){
		if(!$i){
			$pass = 0;
		}elsif(($i eq '+-')or($i eq '?')){
			$pass = 0;
		}elsif($i eq 'Strong'){
			if($evidence == 0){
				$pass = 0;
			}
		}elsif($i eq 'Weak'){
			if($evidence == 1){
				$pass = 0;
			}
		}
	}
	if(!$pass){
		next;
	}
	($start,$end) = ($stemp[3],$stemp[4]);
	($tf,$direction,$role,$relative_dis,$bs_seq) = ($stemp[1],$stemp[5],$stemp[8],$stemp[10],$stemp[11]);
	$len_bs = length($bs_seq);
	#print "length:$len_bs\n";
	if(length($bs_seq) != $length){
		next;
	}
	if($tf ne $chooseTF){
		next;
	}
	$a_r_number[0]=0;
	$a_r_number[1]=0;
	for(my $i=0; $i<$alltf_num; $i++){
		if($tf_name[$i] eq $chooseTF){
			if(($start>=$start[$i])and($end<=$end[$i])){
				$a_r_number[$class[$i]]++;
			}elsif(($start<=$start[$i])and($end>$start[$i])and($end!=$end[$i])){
				$a_r_number[$class[$i]]++;
			}elsif(($start<$end[$i])and($end>=$end[$i])and($start!=$start[$i])){
				$a_r_number[$class[$i]]++;
			}
		}else{
			if(($start>=$start[$i])and($end<=$end[$i])){
				$a_r_number[$class[$i]]++;
			}elsif(($start<=$start[$i])and($end>$start[$i])){
				$a_r_number[$class[$i]]++;
			}elsif(($start<$end[$i])and($end>=$end[$i])){
				$a_r_number[$class[$i]]++;
			}
		}
	}

	#print $line."\n";
	#print "dis:$relative_dis\n";
	if($role eq '+'){
		$role_class = 0;
	}elsif($role eq '-'){
		$role_class = 1;
	}
	if($direction eq 'forward'){
		$dir_class = 0;
	}elsif($direction eq 'reverse'){
		$dir_class = 1;
	}
	#======== -10~2 =========================
	if(($relative_dis >= -10) && ($relative_dis <= 2))
	{
		$bl_overlap=abs($relative_dis+10);
		$br_overlap=abs(2-$relative_dis);
		if($bl_overlap>11)
		{
			$bl_overlap=11;
		}
		if($br_overlap>11)
		{
			$br_overlap=11;
		}
		$b_overlap[$no]=abs($bl_overlap+$br_overlap);
		$b_count[$role_class]++;


	}elsif((($relative_dis-11)> -10)&&(($relative_dis-11)< 2))
	{
		$b_overlap[$no]=abs(2-($relative_dis-11));
		$b_count[$role_class]++;

	}elsif((($relative_dis+11)> -10)&&(($relative_dis+11)< 2))
	{
		$b_overlap[$no]=abs(($relative_dis+11)+10);
		$b_count[$role_class]++;
	}else
	{
		$b_overlap[$no]=0;
	}

	#also include Activation by DNA conformational change
	#and Repression by steric hindrance
	@up_down_poly=();
	if(($relative_dis >= -35) && ($relative_dis <= -10))
	{
		#$most[$temp[0]]++;
		$l_overlap=abs($relative_dis+35);
		#print $l_overlap."\n";
		$r_overlap=abs(-10-$relative_dis);
		#print $r_overlap."\n";
		if($l_overlap < $r_overlap)
		{
			$up_down_poly[$no]=0;
		}elsif($l_overlap > $r_overlap)
		{
			$up_down_poly[$no]=1;
		}else
		{
			$up_down_poly[$no]=2;
		}
		if($l_overlap>11)
		{
			$l_overlap=11;
		}
		if($r_overlap>11)
		{
			$r_overlap=11;
		}

		$overlap[$no]=abs($l_overlap+$r_overlap);
		#print $overlap[$no]."\n";
		#die;
		$count[$role_class]++;
		$dis_poly[$no]=0;

	}elsif((($relative_dis-11)> -35)&&(($relative_dis-11)< -10))
	{
		$overlap[$no]=abs(-10-($relative_dis-11));
		$count[$role_class]++;
		$dis_poly[$no]=0;
		$up_down_poly[$no]=1;

	}elsif((($relative_dis+11)> -35)&&(($relative_dis+11)< -10))
	{
		$overlap[$no]=abs(($relative_dis+11)+35);
		$count[$role_class]++;
		$dis_poly[$no]=0;
		$up_down_poly[$no]=0;

	}elsif(($relative_dis-11)>-10)
	{
		$overlap[$no]=0;
		$dis_poly[$no]=abs(-10-($relative_dis-11));
		$up_down_poly[$no]=1;

	}elsif(($relative_dis+11)<-35)
	{
		$overlap[$no]=0;
		$dis_poly[$no]=abs(($relative_dis+11)+35);
		$up_down_poly[$no]=0;

	}

	#calculate class I activation -95 ~ -60#
	if(($relative_dis >= -95) && ($relative_dis <= -60))
	{
		#$most[$temp[0]]++;
		$c1l_overlap=abs($relative_dis+95);
		#print $l_overlap."\n";
		$c1r_overlap=abs(-60-$relative_dis);
		#print $r_overlap."\n";

		if($c1l_overlap>11)
		{
			$c1l_overlap=11;
		}
		if($c1r_overlap>11)
		{
			$c1r_overlap=11;
		}

		$c1_overlap[$no]=abs($c1l_overlap+$c1r_overlap);
		#print $overlap[$no]."\n";
		#die;
		$c1_count[$role_class]++;


	}elsif((($relative_dis-11)> -95)&&(($relative_dis-11)< -60))
	{
		$c1_overlap[$no]=abs(-60-($relative_dis-11));
		$c1_count[$role_class]++;

	}elsif((($relative_dis+11)> -95)&&(($relative_dis+11)< -60))
	{
		$c1_overlap[$no]=abs(($relative_dis+11)+95);
		$c1_count[$role_class]++;
	}else
	{
		$c1_overlap[$no]=0;
	}

	#calculate class II activation -50~-35 #
	if(($relative_dis >= -50) && ($relative_dis <= -35))
	{
		#$most[$temp[0]]++;
		$c2l_overlap=abs($relative_dis+50);
		#print $l_overlap."\n";
		$c2r_overlap=abs(-35-$relative_dis);
		#print $r_overlap."\n";

		if($c2l_overlap>11)
		{
			$c2l_overlap=11;
		}
		if($c2r_overlap>11)
		{
			$c2r_overlap=11;
		}

		$c2_overlap[$no]=abs($c2l_overlap+$c2r_overlap);
		#print $overlap[$no]."\n";
		#die;
		$c2_count[$role_class]++;


	}elsif((($relative_dis-11)> -50)&&(($relative_dis-11)< -35))
	{
		$c2_overlap[$no]=abs(-35-($relative_dis-11));
		$c2_count[$role_class]++;

	}elsif((($relative_dis+11)> -50)&&(($relative_dis+11)< -35))
	{
		$c2_overlap[$no]=abs(($relative_dis+11)+50);
		$c2_count[$role_class]++;
	}else
	{
		$c2_overlap[$no]=0;
	}

	#calculate Activation by repressor modulation -60~60
	#Repression by DNA looping
	#Cooperative repression
	if(($relative_dis >= -60) && ($relative_dis <= 60))
	{
		#$most[$temp[0]]++;
		$arl_overlap=abs($relative_dis+60);
		#print $l_overlap."\n";
		$arr_overlap=abs(60-$relative_dis);
		#print $r_overlap."\n";

		if($arl_overlap>11)
		{
			$arl_overlap=11;
		}
		if($arr_overlap>11)
		{
			$arr_overlap=11;
		}

		$ar_overlap[$no]=abs($arl_overlap+$arr_overlap);
		$rl_overlap[$no]=abs($arl_overlap+$arr_overlap);
		$cr_overlap[$no]=abs($arl_overlap+$arr_overlap);
		#print "ar_overlap[$no]:$ar_overlap[$no]\n";die;
		#print $overlap[$no]."\n";
		#die;
		$ar_count[$role_class]++;
		$rl_count[$role_class]++;
		$cr_count[$role_class]++;


	}elsif((($relative_dis-11)> -60)&&(($relative_dis-11)< 60))
	{
		$ar_overlap[$no]=abs(60-($relative_dis-11));
		$ar_count[$role_class]++;
		$rl_overlap[$no]=abs(60-($relative_dis-11));
		$rl_count[$role_class]++;
		$cr_overlap[$no]=abs(60-($relative_dis-11));
		$cr_count[$role_class]++;

	}elsif((($relative_dis+11)> -60)&&(($relative_dis+11)< 60))
	{
		$ar_overlap[$no]=abs(($relative_dis+11)+60);
		$ar_count[$role_class]++;
		$rl_overlap[$no]=abs(($relative_dis+11)+60);
		$rl_count[$role_class]++;
		$cr_overlap[$no]=abs(($relative_dis+11)+60);
		$cr_count[$role_class]++;
	}else
	{
		$ar_overlap[$no]=0;
		$rl_overlap[$no]=0;
		$cr_overlap[$no]=0;
	}

	#calculate Repression by roadblock -10~60#
	if(($relative_dis >= -10) && ($relative_dis <= 60))
	{
		#$most[$temp[0]]++;
		$rrl_overlap=abs($relative_dis+10);
		#print $l_overlap."\n";
		$rrr_overlap=abs(60-$relative_dis);
		#print $r_overlap."\n";

		if($rrl_overlap>11)
		{
			$rrl_overlap=11;
		}
		if($rrr_overlap>11)
		{
			$rrr_overlap=11;
		}

		$rr_overlap[$no]=abs($rrl_overlap+$rrr_overlap);
		#print $overlap[$no]."\n";
		#die;
		$rr_count[$role_class]++;


	}elsif((($relative_dis-11)> -10)&&(($relative_dis-11)< 60))
	{
		$rr_overlap[$no]=abs(60-($relative_dis-11));
		$rr_count[$role_class]++;

	}elsif((($relative_dis+11)> -10)&&(($relative_dis+11)< 60))
	{
		$rr_overlap[$no]=abs(($relative_dis+11)+10);
		$rr_count[$role_class]++;
	}else
	{
		$rr_overlap[$no]=0;
	}

	#calculate Repression by activator modulation -95~-10#
	if(($relative_dis >= -95) && ($relative_dis <= -10))
	{
		#$most[$temp[0]]++;
		$ral_overlap=abs($relative_dis+95);
		#print $l_overlap."\n";
		$rar_overlap=abs(-10-$relative_dis);
		#print $r_overlap."\n";

		if($ral_overlap>11)
		{
			$ral_overlap=11;
		}
		if($rar_overlap>11)
		{
			$rar_overlap=11;
		}

		$ra_overlap[$no]=abs($ral_overlap+$rar_overlap);
		#print $overlap[$no]."\n";
		#die;
		$ra_count[$role_class]++;


	}elsif((($relative_dis-11)> -95)&&(($relative_dis-11)< -10))
	{
		$ra_overlap[$no]=abs(-10-($relative_dis-11));
		$ra_count[$role_class]++;

	}elsif((($relative_dis+11)> -95)&&(($relative_dis+11)< -10))
	{
		$ra_overlap[$no]=abs(($relative_dis+11)+95);
		$ra_count[$role_class]++;
	}else
	{
		$ra_overlap[$no]=0;
	}

	#calculate Cooperative activation -95~-35#
	if(($relative_dis >= -95) && ($relative_dis <= -35))
	{
		#$most[$temp[0]]++;
		$cal_overlap=abs($relative_dis+95);
		#print $l_overlap."\n";
		$car_overlap=abs(-35-$relative_dis);
		#print $r_overlap."\n";

		if($cal_overlap>11)
		{
			$cal_overlap=11;
		}
		if($car_overlap>11)
		{
			$car_overlap=11;
		}

		$ca_overlap[$no]=abs($cal_overlap+$car_overlap);
		#print $overlap[$no]."\n";
		#die;
		$ca_count[$role_class]++;


	}elsif((($relative_dis-11)> -95)&&(($relative_dis-11)< -35))
	{
		$ca_overlap[$no]=abs(-35-($relative_dis-11));
		$ca_count[$role_class]++;

	}elsif((($relative_dis+11)> -95)&&(($relative_dis+11)< -35))
	{
		$ca_overlap[$no]=abs(($relative_dis+11)+95);
		$ca_count[$role_class]++;
	}else
	{
		$ca_overlap[$no]=0;
	}

	#calculate Promoter escape regulation -10~10#
	if(($relative_dis >= -10) && ($relative_dis <= 10))
	{
		#$most[$temp[0]]++;
		$perl_overlap=abs($relative_dis+10);
		#print $l_overlap."\n";
		$perr_overlap=abs(-10-$relative_dis);
		#print $r_overlap."\n";

		if($perl_overlap>11)
		{
			$perl_overlap=11;
		}
		if($perr_overlap>11)
		{
			$perr_overlap=11;
		}

		$per_overlap[$no]=abs($perl_overlap+$perr_overlap);
		#print $overlap[$no]."\n";
		#die;
		$per_count[$role_class]++;


	}elsif((($relative_dis-11)> -10)&&(($relative_dis-11)< 10))
	{
		$per_overlap[$no]=abs(10-($relative_dis-11));
		$per_count[$role_class]++;

	}elsif((($relative_dis+11)> -10)&&(($relative_dis+11)< 10))
	{
		$per_overlap[$no]=abs(($relative_dis+11)+10);
		$per_count[$role_class]++;
	}else
	{
		$per_overlap[$no]=0;
	}

	if($overlap[$no]>0)
	{
		$is_overlap=1;
	}else
	{
		$is_overlap=0;
	}
	print OUT "$role_class\t$relative_dis\t$bs_seq\t$is_overlap\t$overlap[$no]\t";
	print OUT "$dis_poly[$no]\t$up_down_poly[$no]\t";
	if($b_overlap[$no]>0)
	{
		$is_b_overlap = 1;
	}else
	{
		$is_b_overlap = 0;
	}
    print OUT "$is_b_overlap\t$b_overlap[$no]\t";#5,6
    print OUT "$a_r_number[0]\t$a_r_number[1]\t";#7,8
		print OUT "$dir_class\t";#9
	if($c1_overlap[$no]>0)
	{
		$is_c1_overlap = 1;
	}else
	{
		$is_c1_overlap = 0;
	}
    print OUT "$is_c1_overlap\t$c1_overlap[$no]\t";#10,11
	if($c2_overlap[$no]>0)
	{
		$is_c2_overlap = 1;
	}else
	{
		$is_c2_overlap = 0;
	}
    print OUT "$is_c2_overlap\t$c2_overlap[$no]\t";#12,13
	if($ar_overlap[$no]>0)
	{
		$is_ar_overlap = 1;
	}else
	{
		$is_ar_overlap = 0;
	}
    print OUT "$is_ar_overlap\t$ar_overlap[$no]\t";#14,15
	if($rr_overlap[$no]>0)
	{
		$is_rr_overlap = 1;
	}else
	{
		$is_rr_overlap = 0;
	}
    print OUT "$is_rr_overlap\t$rr_overlap[$no]\t";#16,17
	if($rl_overlap[$no]>0)
	{
		$is_rl_overlap = 1;
	}else
	{
		$is_rl_overlap = 0;
	}
    print OUT "$is_rl_overlap\t$rl_overlap[$no]\t";#18,19

	if($ra_overlap[$no]>0)
	{
		$is_ra_overlap = 1;
	}else
	{
		$is_ra_overlap = 0;
	}
    print OUT "$is_ra_overlap\t$ra_overlap[$no]\t";#20,21
	if($ca_overlap[$no]>0)
	{
		$is_ca_overlap = 1;
	}else
	{
		$is_ca_overlap = 0;
	}
    print OUT "$is_ca_overlap\t$ca_overlap[$no]\t";#22,23
	if($cr_overlap[$no]>0)
	{
		$is_cr_overlap = 1;
	}else
	{
		$is_cr_overlap = 0;
	}
    print OUT "$is_cr_overlap\t$cr_overlap[$no]\t";#24,25
	if($per_overlap[$no]>0)
	{
		$is_per_overlap = 1;
	}else
	{
		$is_per_overlap = 0;
	}
    print OUT "$is_per_overlap\t$per_overlap[$no]\n";#26,27
		$no++;
}
close REGDB;
`perl $location_independent $output`;
`rm $output`;
