#! usr/bin/perl
my $input = $ARGV[0];
if(@ARGV<1){
	print "useage ATCG_feature_extract.pl feature_dependent.txt\n";die;
}
my $output = $input;
$output =~ s/.txt//g;
$output .= "_svmformat.txt";
open DATA,">",$output;
#useage ATCG_feature_extract.pl feature_dependent.txt
#print "111";
#   A      C      G      T
#0  259    271    253    267
#1  491.2  467.2  507.2  482.2
#2  15200  9300   13700  9600
$num=0;
open IFILE,"<",$input;
while(my $line =<IFILE>)
{
	chomp $line;
	$line =~ tr/A-Z/a-z/;
	#print $_."\n";
	#die;
	@temp=split(/\t/,$line);
	#print "\n";
	#print $temp[0]."\t".$temp[1]."\t".$temp[2]."\t".$temp[3]."\t".$temp[4]."\t".$temp[5]."\n";
	#print $temp[6]."\t".$temp[7]."\t".$temp[8]."\t".$temp[9]."\t".$temp[10]."\n";
	#die;

	push(@class,$temp[0]);
	push(@distance,$temp[1]);
	@allsequence[$num]=$temp[2];
	# print $temp[0]."\n";
	# print $temp[1]."\t";
	# print $temp[2];
	# die;
	@{$sequence[$num]}=split(//,$temp[2]);
	#pop(@{$sequence[$num]});
	#print @{$sequence[$num]};
	#print "\n";
	#die;
	push(@d_overlap,$temp[3]);
	push(@num_overlap,$temp[4]);
	push(@distance_poly,$temp[5]);
	push(@up_down,$temp[6]);
	push(@d_overlap_bubble,$temp[7]);
	push(@number_overlap_bubble,$temp[8]);
	push(@repressor_o_tf,$temp[9]);
	push(@activator_o_tf,$temp[10]);
	push(@forward_reverse,$temp[11]);

	push(@classI,$temp[12]);
	push(@classI_number,$temp[13]);
	push(@classII,$temp[14]);
	push(@classII_number,$temp[15]);
	push(@a_r,$temp[16]);
	push(@a_r_number,$temp[17]);
	push(@r_r,$temp[18]);
	push(@r_r_number,$temp[19]);
	push(@r_l,$temp[20]);
	push(@r_l_number,$temp[21]);
	push(@r_a,$temp[22]);
	push(@r_a_number,$temp[23]);
	push(@c_a,$temp[24]);
	push(@c_a_number,$temp[25]);
	push(@c_r,$temp[26]);
	push(@c_r_number,$temp[27]);
	push(@per,$temp[28]);
	push(@per_number,$temp[29]);

	#@atemp=split(//,$temp[30]);
	#print "$atemp[0]\n$atemp[1]";
	#pop(@atemp);
	#print $atemp[0];
	#push(@CRP_mbs_number,$atemp[0]);
	#die;
	$num++;
	#print @{$sequence[0]}."\n";
	#print $sequence[0][41]."\n";
	#die;
}


open(MOTIF,"tmp/4-mer.txt")||die("open 4-mer.txt error\n");
while(<MOTIF>)
{
	chomp;
	tr/A-Z/a-z/;
	@temp_m=split(/\t/);
	#print $temp_m[1]."\n";
	#print $_."\n";
	push(@motif_4,$temp_m[1]);

}
close(MOTIF);
open(MOTIF2,"tmp/3-mer.txt")||die("open 3-mer.txt error\n");
while(<MOTIF2>)
{
    chomp;
    tr/A-Z/a-z/;
    @temp2_m=split(/\t/);
    #print $temp2_m[1]."\n";
    #die;
    #print $_."\n";
    push(@motif_3,$temp2_m[1]);

}
close(MOTIF2);
#print $motif_4[0]."\n";

#print $allsequence[0]."\n";

#$ss=substr($allsequence[0],0,4);
#$ss1=substr($allsequence[0],1,4);
#print $ss."\t$ss1\n";
#die;


for($i=0;$i<$num;$i++)
{
	$a_c=0;
	$c_c=0;
	$g_c=0;
	$t_c=0;
	$tra_a=0;
	$tra_c=0;
	$tra_g=0;
	$tra_t=0;
	$trc_c=0;
	$trc_g=0;
	$trc_t=0;
	$trg_g=0;
	$trg_t=0;
	$trt_t=0;
	@distribution_a=();
	@distribution_c=();
	@distribution_g=();
	@distribution_t=();
	%motif_c=();
    %motif_c2=();
	@feature=();
	#print $allsequence[$i]."\n";
	for($k=0;$k<@{$sequence[$i]}-3;$k++)
	{
		#print "111\n";
		$string_temp=substr($allsequence[$i],$k,4);
		$motif_c{"$string_temp"}++;
		#print $string_temp."\n";
	}
    for($k=0;$k<@{$sequence[$i]}-2;$k++)
    {
        #print "111\n";
        $string_temp2=substr($allsequence[$i],$k,3);
        $motif_c2{"$string_temp2"}++;
        #print $string_temp2."\n";
        #die;
    }

	for($j=0;$j<=@{$sequence[$i]}+0;$j++)
	{
		if($sequence[$i][$j] eq 'a')
		{
			push(@distribution_a,$j);
			$a_c++;
			if(($j+1) == @{$sequence[$i]}+0)
			{
				last;
			}
			if($sequence[$i][$j+1] eq 'a')
			{
				$tra_a++;
			}elsif($sequence[$i][$j+1] eq 'c')
			{
				$tra_c++;
			}elsif($sequence[$i][$j+1] eq 'g')
			{
				$tra_g++;
			}elsif($sequence[$i][$j+1] eq 't')
			{
				$tra_t++;
			}

		}elsif($sequence[$i][$j] eq 'c')
		{
			push(@distribution_c,$j);
			$c_c++;
			if(($j+1) == @{$sequence[$i]}+0)
			{
				last;
			}
			if($sequence[$i][$j+1] eq 'a')
			{
				$tra_c++;
			}elsif($sequence[$i][$j+1] eq 'c')
			{
				$trc_c++;
			}elsif($sequence[$i][$j+1] eq 'g')
			{
				$trc_g++;
			}elsif($sequence[$i][$j+1] eq 't')
			{
				$trc_t++;
			}
		}elsif($sequence[$i][$j] eq 'g')
		{
			push(@distribution_g,$j);
			$g_c++;
			if(($j+1) == @{$sequence[$i]}+0)
			{
				last;
			}
			if($sequence[$i][$j+1] eq 'a')
			{
				$tra_g++;
			}elsif($sequence[$i][$j+1] eq 'c')
			{
				$trc_g++;
			}elsif($sequence[$i][$j+1] eq 'g')
			{
				$trg_g++;
			}elsif($sequence[$i][$j+1] eq 't')
			{
				$trg_t++;
			}
		}elsif($sequence[$i][$j] eq 't')
		{
			push(@distribution_t,$j);
			$t_c++;
			if(($j+1) == @{$sequence[$i]}+0)
			{
				last;
			}
			if($sequence[$i][$j+1] eq 'a')
			{
				$tra_t++;
			}elsif($sequence[$i][$j+1] eq 'c')
			{
				$trc_t++;
			}elsif($sequence[$i][$j+1] eq 'g')
			{
				$trg_t++;
			}elsif($sequence[$i][$j+1] eq 't')
			{
				$trt_t++;
			}
		}else
		{
			print $i."\t".$j."\t".$sequence[$i][$j]."\n";
			print "Wrong acgt sequence \n";
			print "1111";
			die;
		}
	}
	$seqlen=@{$sequence[$i]}+0;
	$qa=$a_c/$seqlen;
	$qc=$c_c/$seqlen;
	$qg=$g_c/$seqlen;
	$qt=$t_c/$seqlen;

	#print @motif_4."\n";
	#die;
	for($k=0;$k<@motif_4+0;$k++)
	{
		#print $motif_4[$k]."\t".$motif_c{"$motif_4[$k]"}."\t";
		$feature[$i][$k]=0;
		$feature[$i][$k]=0+$motif_c{"$motif_4[$k]"};
		#print $feature[$i][$k]."\n";

	}
	#print $motif_4[0]."\t".$feature[0][0]."\n$motif_4[255]\t$feature[$i][255] \n";
	#print "\n\n\n";

	#die;

	$feature[$i][256]=(($a_c*259)+($c_c*271)+($g_c*253)+($t_c*267))/$seqlen;
	$feature[$i][257]=(($a_c*491.2)+($c_c*467.2)+($g_c*507.2)+($t_c*482.2))/$seqlen;
	$feature[$i][258]=(($a_c*15200)+($c_c*9300)+($g_c*13700)+($t_c*9600))/$seqlen;

	$feature[$i][259]=$qa;
	$feature[$i][260]=$qc;
	$feature[$i][261]=$qg;
	$feature[$i][262]=$qt;

	$feature[$i][263]=($qa**2)+($qc**2)+($qg**2)+($qt**2);
	$feature[$i][264]=-(($qa * log($qa)/log(10))+($qc * log($qc)/log(10))+($qg * log($qg)/log(10))+($qt * log($qt)/log(10)));# There are some problem #
	#print $qa."\t".$qc."\t".$qg."\t".$qt."\n\n";
	#$feature[$i][4]=0.518826;
	$feature[$i][265]=(-1/$feature[$i][264])*($qa * log($qa)/log(10));
	$feature[$i][266]=(-1/$feature[$i][264])*($qc * log($qc)/log(10));
	$feature[$i][267]=(-1/$feature[$i][264])*($qg * log($qg)/log(10));
	$feature[$i][268]=(-1/$feature[$i][264])*($qt * log($qt)/log(10));

	$feature[$i][269]=$tra_a/($seqlen-1);
	$feature[$i][270]=$tra_c/($seqlen-1);
	$feature[$i][271]=$tra_g/($seqlen-1);
	$feature[$i][272]=$tra_t/($seqlen-1);
	$feature[$i][273]=$trc_c/($seqlen-1);
	$feature[$i][274]=$trc_g/($seqlen-1);
	$feature[$i][275]=$trc_t/($seqlen-1);
	$feature[$i][276]=$trg_g/($seqlen-1);
	$feature[$i][277]=$trg_t/($seqlen-1);
	$feature[$i][278]=$trt_t/($seqlen-1);

	if((@distribution_a+0) >= 5)
	{
		$mid=0;
		$mid25=0;
		$mid75=0;

		$mid=int((@distribution_a+0)/2);
		$mid25=int($mid/2);
		$mid75=$mid25+$mid;
		#print $mid_a."\t".$mid_a25."\t".$mid_a75."\n";

		$feature[$i][279]=(($distribution_a[0]+1)/$seqlen)*100;
		$feature[$i][280]=(($distribution_a[$mid25]+1)/$seqlen)*100;
		$feature[$i][281]=(($distribution_a[$mid]+1)/$seqlen)*100;
		$feature[$i][282]=(($distribution_a[$mid75]+1)/$seqlen)*100;
		$feature[$i][283]=(($distribution_a[(@distribution_a-1)]+1)/$seqlen)*100;

	#	for($k=0;$k<@distribution_a+0;$k++)
	#	{
	#		print $distribution_a[$k]."\n";
	#	}
	}else
	{
		if((@distribution_a+0) == 4)
		{
			$feature[$i][279]=(($distribution_a[0]+1)/$seqlen)*100;
			$feature[$i][280]=(($distribution_a[1]+1)/$seqlen)*100;
			$feature[$i][281]=(($distribution_a[2]+1)/$seqlen)*100;
			$feature[$i][282]=0;
			$feature[$i][283]=(($distribution_a[3]+1)/$seqlen)*100;
		}elsif((@distribution_a+0) == 3)
		{
			$feature[$i][279]=(($distribution_a[0]+1)/$seqlen)*100;
			$feature[$i][280]=0;
			$feature[$i][281]=(($distribution_a[1]+1)/$seqlen)*100;
			$feature[$i][282]=0;
			$feature[$i][283]=(($distribution_a[2]+1)/$seqlen)*100;
		}elsif((@distribution_a+0) == 2)
		{
			$feature[$i][279]=(($distribution_a[0]+1)/$seqlen)*100;
			$feature[$i][280]=0;
			$feature[$i][281]=0;
			$feature[$i][282]=0;
			$feature[$i][283]=(($distribution_a[1]+1)/$seqlen)*100;

		}else
		{
			$feature[$i][279]=(($distribution_a[0]+1)/$seqlen)*100;
			$feature[$i][280]=0;
			$feature[$i][281]=0;
			$feature[$i][282]=0;
			$feature[$i][283]=0;
		}
	}
	#print
	#print (@distribution_c+0)."\n";

	if((@distribution_c+0) >= 5)
	{

		$mid=0;
		$mid25=0;
		$mid75=0;

		$mid=int((@distribution_c+0)/2);
		$mid25=int($mid/2);
		$mid75=$mid25+$mid;
		#print $mid_a."\t".$mid_a25."\t".$mid_a75."\n";

		$feature[$i][284]=(($distribution_c[0]+1)/$seqlen)*100;
		$feature[$i][285]=(($distribution_c[$mid25]+1)/$seqlen)*100;
		$feature[$i][286]=(($distribution_c[$mid]+1)/$seqlen)*100;
		$feature[$i][287]=(($distribution_c[$mid75]+1)/$seqlen)*100;
		$feature[$i][288]=(($distribution_c[(@distribution_c-1)]+1)/$seqlen)*100;

		#for($k=0;$k<@distribution_a+0;$k++)
		#{
		#	print $distribution_a[$k]."\n";
		#}
	}else
	{
		if((@distribution_c+0) == 4)
		{
			$feature[$i][284]=(($distribution_c[0]+1)/$seqlen)*100;
			$feature[$i][285]=(($distribution_c[1]+1)/$seqlen)*100;
			$feature[$i][286]=(($distribution_c[2]+1)/$seqlen)*100;
			$feature[$i][287]=0;
			$feature[$i][288]=(($distribution_c[3]+1)/$seqlen)*100;
		}elsif((@distribution_c+0) == 3)
		{
			$feature[$i][284]=(($distribution_c[0]+1)/$seqlen)*100;
			$feature[$i][285]=0;
			$feature[$i][286]=(($distribution_c[1]+1)/$seqlen)*100;
			$feature[$i][287]=0;
			$feature[$i][288]=(($distribution_c[2]+1)/$seqlen)*100;
		}elsif((@distribution_c+0) == 2)
		{
			$feature[$i][284]=(($distribution_c[0]+1)/$seqlen)*100;
			$feature[$i][285]=0;
			$feature[$i][286]=0;
			$feature[$i][287]=0;
			$feature[$i][288]=(($distribution_c[1]+1)/$seqlen)*100;

		}else
		{
			$feature[$i][284]=(($distribution_c[0]+1)/$seqlen)*100;
			$feature[$i][285]=0;
			$feature[$i][286]=0;
			$feature[$i][287]=0;
			$feature[$i][288]=0;
		}
	}

	if((@distribution_g+0) >= 5)
	{
		$mid=0;
		$mid25=0;
		$mid75=0;

		$mid=int((@distribution_g+0)/2);
		$mid25=int($mid/2);
		$mid75=$mid25+$mid;
		#print $mid25."\t".$mid."\t".$mid75."\n";

		$feature[$i][289]=(($distribution_g[0]+1)/$seqlen)*100;
		$feature[$i][290]=(($distribution_g[$mid25]+1)/$seqlen)*100;
		$feature[$i][291]=(($distribution_g[$mid]+1)/$seqlen)*100;
		$feature[$i][292]=(($distribution_g[$mid75]+1)/$seqlen)*100;
		$feature[$i][293]=(($distribution_g[(@distribution_g-1)]+1)/$seqlen)*100;

		#$tg=($distribution_g[0]+1);
		#$tg1=($distribution_g[$mid25]+1);
		#$tg2=($distribution_g[$mid]+1);
		#$tg3=($distribution_g[$mid75]+1);
		#$tg4=($distribution_g[(@distribution_g-1)]+1);
		#print $tg."\t$tg1\t$tg2\t$tg3\t$tg4\n";
		#for($k=0;$k<@distribution_a+0;$k++)
		#{
		#	print $distribution_a[$k]."\n";
		#}
	}else
	{
		if((@distribution_g+0) == 4)
		{
			$feature[$i][289]=(($distribution_g[0]+1)/$seqlen)*100;
			$feature[$i][290]=(($distribution_g[1]+1)/$seqlen)*100;
			$feature[$i][291]=(($distribution_g[2]+1)/$seqlen)*100;
			$feature[$i][292]=0;
			$feature[$i][293]=(($distribution_g[3]+1)/$seqlen)*100;
		}elsif((@distribution_g+0) == 3)
		{
			$feature[$i][289]=(($distribution_g[0]+1)/$seqlen)*100;
			$feature[$i][290]=0;
			$feature[$i][291]=(($distribution_g[1]+1)/$seqlen)*100;
			$feature[$i][292]=0;
			$feature[$i][293]=(($distribution_g[2]+1)/$seqlen)*100;
		}elsif((@distribution_g+0) == 2)
		{
			$feature[$i][289]=(($distribution_g[0]+1)/$seqlen)*100;
			$feature[$i][290]=0;
			$feature[$i][291]=0;
			$feature[$i][292]=0;
			$feature[$i][293]=(($distribution_g[1]+1)/$seqlen)*100;

		}else
		{
			$feature[$i][289]=(($distribution_g[0]+1)/$seqlen)*100;
			$feature[$i][290]=0;
			$feature[$i][291]=0;
			$feature[$i][292]=0;
			$feature[$i][293]=0;
		}
	}

	if((@distribution_t+0) >= 5)
	{
		$mid=0;
		$mid25=0;
		$mid75=0;

		$mid=int((@distribution_t+0)/2);
		$mid25=int($mid/2);
		$mid75=$mid25+$mid;
		#print $mid_a."\t".$mid_a25."\t".$mid_a75."\n";

		$feature[$i][294]=(($distribution_t[0]+1)/$seqlen)*100;
		$feature[$i][295]=(($distribution_t[$mid25]+1)/$seqlen)*100;
		$feature[$i][296]=(($distribution_t[$mid]+1)/$seqlen)*100;
		$feature[$i][297]=(($distribution_t[$mid75]+1)/$seqlen)*100;
		$feature[$i][298]=(($distribution_t[(@distribution_t-1)]+1)/$seqlen)*100;

		#for($k=0;$k<@distribution_a+0;$k++)
		#{
		#	print $distribution_a[$k]."\n";
		#}
	}else
	{
		if((@distribution_t+0) == 4)
		{
			$feature[$i][294]=(($distribution_t[0]+1)/$seqlen)*100;
			$feature[$i][295]=(($distribution_t[1]+1)/$seqlen)*100;
			$feature[$i][296]=(($distribution_t[2]+1)/$seqlen)*100;
			$feature[$i][297]=0;
			$feature[$i][298]=(($distribution_t[3]+1)/$seqlen)*100;
		}elsif((@distribution_t+0) == 3)
		{
			$feature[$i][294]=(($distribution_t[0]+1)/$seqlen)*100;
			$feature[$i][295]=0;
			$feature[$i][296]=(($distribution_t[1]+1)/$seqlen)*100;
			$feature[$i][297]=0;
			$feature[$i][298]=(($distribution_t[2]+1)/$seqlen)*100;
		}elsif((@distribution_t+0) == 2)
		{
			$feature[$i][294]=(($distribution_t[0]+1)/$seqlen)*100;
			$feature[$i][295]=0;
			$feature[$i][296]=0;
			$feature[$i][297]=0;
			$feature[$i][298]=(($distribution_t[1]+1)/$seqlen)*100;

		}else
		{
			$feature[$i][294]=(($distribution_t[0]+1)/$seqlen)*100;
			$feature[$i][295]=0;
			$feature[$i][296]=0;
			$feature[$i][297]=0;
			$feature[$i][298]=0;
		}
	}



	if($distance[$i] > 0)
	{
		$feature[$i][299]=0;
		$feature[$i][300]=$distance[$i];
	}else
	{
		$feature[$i][299]=1;
		$feature[$i][300]=abs($distance[$i]);
	}



	$feature[$i][301]=$d_overlap[$i];#302
	$feature[$i][302]=$num_overlap[$i];#303
	$feature[$i][303]=$distance_poly[$i];#304
	$feature[$i][304]=$up_down[$i];#305
	$feature[$i][305]=$d_overlap_bubble[$i];#306
	$feature[$i][306]=$number_overlap_bubble[$i];#307
	$feature[$i][307]=$repressor_o_tf[$i];#308
	$feature[$i][308]=$activator_o_tf[$i];#309
	$feature[$i][309]=$forward_reverse[$i];#310
	$feature[$i][310]=@classI[$i];#311
	$feature[$i][311]=@classI_number[$i];#312
	$feature[$i][312]=@classII[$i];#313
	$feature[$i][313]=@classII_number[$i];#314
	$feature[$i][314]=@a_r[$i];#315
	$feature[$i][315]=@a_r_number[$i];#316
	$feature[$i][316]=@r_r[$i];#317
	$feature[$i][317]=@r_r_number[$i];#318
	$feature[$i][318]=@r_l[$i];#319
	$feature[$i][319]=@r_l_number[$i];#320
	$feature[$i][320]=@r_a[$i];#321
	$feature[$i][321]=@r_a_number[$i];#322
	$feature[$i][322]=@c_a[$i];#323
	$feature[$i][323]=@c_a_number[$i];#324
	$feature[$i][324]=@c_r[$i];#325
	$feature[$i][325]=@c_r_number[$i];#326
	$feature[$i][326]=@per[$i];#327
	$feature[$i][327]=@per_number[$i];#328
    for($k=0;$k<@motif_3+0;$k++)
    {
        $no = $k+328;
        $feature[$i][$no]=0;
        $feature[$i][$no]=0+$motif_c2{"$motif_3[$k]"};
        #print $feature[$i][$k]."\n";
    }

	print DATA $class[$i]." ";
	for($l=0;$l<@{$feature[$i]}+0;$l++)
	{
		$f_number=$l+1;
		print DATA $f_number.":".$feature[$i][$l]." ";
		#$f_number++;
	}
	print DATA "\n";

}
