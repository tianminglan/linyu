#!/usr/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;

#######################
###USAGE
########################

sub USAGE{
	my $usage=<<"USAGE";

Version:1.0
Wed Mar 27 16:27:23 HKT 2019   linyu\@genomics.cn


	Usage:perl $0 -i fq.gz -name samplename -formatdb formatdb -blastall blastall -Rpath R -database mtDNA.fa -coveragemin num -coveragemax num -similaritymin num -similaritymax num -delete yes/no -o outputpath

              eg.perl $0 -i /ldfssz1/ST_HEALTH/Population_Genomics/PMO/Ancient_DNA/main/test.fq.gz -name sample1 -formatdb /zfssz3/ST_HEALTH/Population_Genomics/USER/lantianming/bin/blast-2.2.26/bin/formatdb -blastall /zfssz3/ST_HEALTH/Population_Genomics/USER/lantianming/bin/blast-2.2.26/bin/blastall -Rpath /hwfssz1/ST_HEALTH/Population_Genomics/User/xuwenhao/software/link/app/R-3.2.1/bin/R -database /ldfssz1/ST_HEALTH/Population_Genomics/PMO/Ancient_DNA/main/animal_mtDNA.fa -coveragemin 0.1 -coveragemax 1.00 -similaritymin 10 -similaritymax 100  -mutationsite 10 -ct_mutationnum 1 -ga_mutationnum 1 -delete yes -o /ldfssz1/ST_HEALTH/Population_Genomics/PMO/Ancient_DNA/main
	
	where:

	Options
		-i <file>               : Fastq file after adapter filtration.
		-name <character>       : Sample name that used for identify different samples.
		-formatdb <file>        : Absolute path of Formatdb. If you have formatdbed database,then ignore this.
		-blastall <file>        : Blastall Absolute path.
		-Rpath    <file>	: R Absolute path.
		-database <file>        : MtDNA database.
		-coveragemin <num>      : The minimal query coverage.
		-coveragemax <num>      : The maximal query coverage.
		-similaritymin  <num>   : The minimal sequence similarity.
		-similaritymax  <num>   : The maximal sequence similarity.
		-mutationsite	<num>	: The first and last N bases that used for screening reads with the C->T and G->A changes.
		-ct_mutationnum <num>	: The number of C to T changes within the -mutationsite <num>. This should be used with the parameter -mutationsite.  
		-ga_mutationnum <num>   : The number of G to A changes within the -mutationsite <num>. This should be used with the parameter -mutationsite.
		-delete <character>     : Deleting the temp files, yes or no.
		-o <path>               : The output directory.
		-h                      : Show this help.
USAGE
print $usage;
exit;
}

######### Check Parameters And Write Into Memory ######

my ($fq,$sample,$formatdb,$blastall,$Rpath,$database,$coveragemin,$coveragemax,$similaritymin,$similaritymax,$CT,$c_t_mutation,$g_a_mutation,$delete,$output);
GetOptions("help|?" =>\&USAGE,
	   "i:s" => \$fq,
	   "name:s" => \$sample,
	   "formatdb:s" => \$formatdb,
	   "blastall:s" => \$blastall,
	   "Rpath:s" => \$Rpath,
	   "database:s" => \$database,
	   "coveragemin:s" => \$coveragemin,
	   "coveragemax:s" => \$coveragemax,
	   "similaritymin:s" => \$similaritymin,
	   "similaritymax:s" => \$similaritymax,
	   "mutationsite:s" => \$CT,
	   "ct_mutationnum:s" => \$c_t_mutation,
	   "ga_mutationnum:s" => \$g_a_mutation,
	   "delete:s" => \$delete,
	   "o:s" => \$output
	  ) or &USAGE;
&USAGE unless ($fq && $sample && $blastall && $database && $coveragemin && $coveragemax && $similaritymin && $delete && $output);

######### Check The Coverage and Similarity Format ######

if ($coveragemin>$coveragemax){
	print "ERROR:coveragemin cannot be greater than coveragemax\n";
	exit;
}
if ($similaritymin>$similaritymax){
	print "ERROR:similaritymin cannot be greater than similaritymax\n";
	exit;
}
if ($coveragemin>1){
	print "Message:You must ensure the coveragemin fomat it must less then 1\n        eg. 98% ,you should input 0.98,not 98\n";
	exit;
} 
if ($coveragemax>1){
	print "Message:You must ensure the coveragemax fomat it must less then 1\n        eg. 98% ,you should input 0.98,not 98\n";
	exit;
}
if ($similaritymin>100 || $similaritymax>100){
	print "ERROR:Similarity must less than 100\n";
	exit;
} 
if ($similaritymin<1){
	print "Message:You must ensure the similaritymin fomat,eg. 0.98 not mean 98%\n        if you want to input 98%,you must input 98\n        if you want to input 0.98%,you must input 0.98\n";
	print "        if you ensure you are right,please ignore this message\n";
}
if ($similaritymax<1){
	print "Message:You must ensure the similaritymax fomat,eg. 1 not mean 100%,if you want to input 100%,you must input 100\n        if you want to input 1%,you must input 1 \n";
	print "        if you ensure you are right,please ignore this message\n";
}

######### Creat The Database Index ######
if (defined $formatdb) {
	system("$formatdb \-i $database \-p F \-a F \-o T");
}
######### Filter Low Quality Reads ######

open FQ, "gzip -dc $fq|";
open FilterFQ,">$output/$sample.filterfq.fq";

while(my $line1=<FQ>, my $line2=<FQ>,my $line3=<FQ>,my $line4=<FQ>){
	chomp($line1);chomp($line2);chomp($line3);chomp($line4);
	my @base=split//,$line2;
	my @basequality=split//,$line4;
	my $len=length($line2);
	next if $len<30;
	my $sum_topbasequality=0;
	foreach my $position (0..5){
		$sum_topbasequality+=ord($basequality[$position])-33;
	}
	my $sum_lastbasequality=0;
	foreach my $position (($#basequality-5+1)..$#basequality){
		$sum_lastbasequality+=ord($basequality[$position])-33;
	}
	my $avetop=$sum_topbasequality/5;my $avelast=$sum_lastbasequality/5;
	next if ($avetop<15 || $avelast<15);
	my $numN=0;
	foreach my $N (@base) {
		if ($N eq "N"){$numN +=1;}
	}
	next if ($numN/$len>0.1);
	my ($i,$m,$n,$k);
	for ($i=0;$i<=$len-1;$i++){
		last if ($base[$i]=~/[^N]/);
	}
	for ($k=$len-1;$k>=0;$k--){
		last if ($base[$k]=~/[^N]/);
	}
	for ($m=0;$m<=$len-1;$m++){
		my $c=ord($basequality[$m])-33;
		last if ($c>15);
	}
	for ($n=$len-1;$n>=0;$n--){
		my $d=ord($basequality[$n])-33;
		last if ($d>15);
	}
	my ($min,$max);
	if ($i>$m) {
		$min = $i;
	} else { 
		$min = $m;
	}
	if ($k<$n) {
		$max = $k;
	} else {
		$max = $n;
	}
	my $base=join("", @base[$min..$max]);
	my $basequality=join("", @basequality[$min..$max]);
	my $len1=length($base);
	if ($len1>30){
		print  FilterFQ "$line1\n$base\n$line3\n$basequality\n";
	}
}

close FQ;
close FilterFQ;

######### Change Fastq To Fasta ######

open FilterFQ,"<$output/$sample.filterfq.fq";
open FASTA,">$output/$sample.filterfasta.fa";

while (my $line1=<FilterFQ>,my $line2=<FilterFQ>,my $line3=<FilterFQ>,my $line4=<FilterFQ>) {
        chomp($line1);chomp($line2);chomp($line3);chomp($line4);
	$line1 =~ s/^@/>/;
	print FASTA "$line1\n$line2\n";
}
						
close FASTA;
close FilterFQ;

######### Add Length Message ######

open FASTA,"<$output/$sample.filterfasta.fa";
open Addlength,">$output/$sample.adlengthfasta.fa";

while (my $line1=<FASTA>,my $line2=<FASTA>){
        chomp($line1);chomp($line2);
	my $len=length($line2);
	$line1 .=" length=$len";
	print Addlength "$line1\n$line2\n";
}

close FASTA;
close Addlength;

######### Blast Align ######

my $input="$sample.adlengthfasta.fa";
system("$blastall \-i $output/$input \-d $database \-o $output/$input.out \-p blastn \-F F \-e 1e\-10 \-m 3 \-a 8");

######### Add Blast Message In The End ######

open Blastout,"<$output/$sample.adlengthfasta.fa.out";
open AddBLASTN,">$output/$sample.adlengthfasta.fa.out.jiablast";

while(my $line=<Blastout>) {
        chomp($line);
	last if ($line =~ /  Database:/); 
	print  AddBLASTN "$line\n";
}

close Blastout;
close AddBLASTN;

######### Choose The aligned Reads ######

open AddBLASTN,"<$output/$sample.adlengthfasta.fa.out.jiablast";
open Choose,">$output/$sample.adlengthfasta.fa.out.jiablast.choose";

my $seq="";
while(my $line=<AddBLASTN>){
	if ($line !~ /^BLASTN/) {
		$seq .= $line;
		next;
	} else {
		$seq .=$line;
		if ($seq =~ /No hits found/) {
			$seq = "";
			next;
		} else {
			my @line = split /\n/, $seq;
			for (my $i=0; $i<@line; $i++) {
				if ($line[$i] =~ /^Query/) {
					print Choose "$line[$i] $line[($i+1)]\n";
				} elsif ($line[$i] =~ /Score    E/) {
					for (my $j=$i; $j<$#line; $j++) {
						print Choose "$line[$j]\n";
					}
					print Choose "#####\n";
					$seq="";
					last;
				}
			}
		}
	}
}

close AddBLASTN;
close Choose;

######### Stitching information ######

open Choose,"<$output/$sample.adlengthfasta.fa.out.jiablast.choose";
open Jiewei,">$output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei";

my $sequence="";
my $j=0;
my %h1;
my %h2;
while(my $line=<Choose>) {
	if ($line!~/#####/) {
		$sequence .=$line;
		next;			
	} else {
		my @a =split /\n/,$sequence;
		for (my $i=0;$i<@a;$i++){
			chomp($a[$i]);
			if ($j==0) {
				if (length($a[$i])==0) {
					$j++;
					print  Jiewei "$a[$i]\n";
					next;
				}
				print  Jiewei "$a[$i]\n";
			} elsif ($j==1) {
				if (length($a[$i])==0) {
                                        $j++;
					next;
                                }
				my @gbl=split/\|/,$a[$i];
				my $id = (split /\./, $gbl[1])[0];
				$h1{$id}=$gbl[2];
			} elsif ($j==2) {
				if (length($a[$i])==0) {
                                        $j++;
					next;
                                }
				my @blast1=split/\s+/,$a[$i];
				if ($blast1[1]==0) {        #0 delete
					next;
				}
				$h2{$blast1[0]." ".$blast1[3]}=$a[$i];
			} elsif ($j==3) {
				if (length($a[$i])==0) {
                                        $j++;
					next;
				}
				my @blast2=split/\s+/,$a[$i];
				if ($blast2[1]==0) {        #0 delete
                                 	 next;
                                  }
				my $s=$blast2[1]-1;
				my $e=$blast2[1]+1;
				if (($blast2[3]-$blast2[1])>=0) {
					my $second1 =join(" ",@blast2[2..3]);
					if (exists $h2{$blast2[0]." ".$s}) {
						my @cc =split /\s+/,$h2{$blast2[0]." ".$s};
						my @dd =split /\S+/,$h2{$blast2[0]." ".$s};
						$h2{$blast2[0]." "."new"." ".$blast2[3]}=$cc[0].$dd[1].$cc[1].$dd[2].$cc[2].$second1;
						delete $h2{$blast2[0]." ".$s};
					} else {
						next;
					}
				} elsif (($blast2[3]-$blast2[1])<=0) {
					my $second2 =join(" ",@blast2[2..3]);
					if (exists $h2{$blast2[0]." ".$e}) {
                                                my @ee =split /\s+/,$h2{$blast2[0]." ".$e};
                                                my @ff =split /\S+/,$h2{$blast2[0]." ".$e};
                                                $h2{$blast2[0]." "."new"." ".$blast2[3]}=$ee[0].$ff[1].$ee[1].$ff[2].$ee[2].$second2;
						 delete $h2{$blast2[0]." ".$e};
                                        } else {
                                                next;
                                        }
				}
			} elsif ($j==4) {
				if (length($a[$i])==0) {
                                        $j++;
                                        next;
                                }
				my @blast3=split/\s+/,$a[$i];
				if ($blast3[1]==0) {        #0 delete
                                           next;
                                }
				my $t=$blast3[1]-1;
				my $d=$blast3[1]+1;
				if (($blast3[3]-$blast3[1])>=0) {
					my $third1 =join(" ",@blast3[2..3]);
					if (exists $h2{$blast3[0]." "."new"." ".$t}) {
						my @gg =split /\s+/,$h2{$blast3[0]." "."new"." ".$t};
                                        	my @hh =split /\S+/,$h2{$blast3[0]." "."new"." ".$t};
                                        	$h2{$blast3[0]." "."new1"." ".$blast3[3]}=$gg[0].$hh[1].$gg[1].$hh[2].$gg[2].$third1;
						delete $h2{$blast3[0]." "."new"." ".$t};
					} else {
                                	        next;
                               		}
				} elsif (($blast3[3]-$blast3[1])<=0) {
					my $third2 =join(" ",@blast3[2..3]);
					if (exists $h2{$blast3[0]." "."new"." ".$d}) {
						my @ii =split /\s+/,$h2{$blast3[0]." "."new"." ".$d};
						my @jj =split /\S+/,$h2{$blast3[0]." "."new"." ".$d};
						$h2{$blast3[0]." "."new1"." ".$blast3[3]}=$ii[0].$jj[1].$ii[1].$jj[2].$ii[2].$third2;
						delete $h2{$blast3[0]." "."new"." ".$d};
					} else {
						next;
					}
				}
			} elsif ($j==5) {
				if (length($a[$i])==0) {
                                        $j++;
                                        next;
                                }
				my @blast4=split/\s+/,$a[$i];
				if ($blast4[1]==0) {        #0 delete
                                           next;
                                }
                                my $fo=$blast4[1]-1;
                                my $ur=$blast4[1]+1;
                                if (($blast4[3]-$blast4[1])>=0) {
                                        my $four1 =join(" ",@blast4[2..3]);
                                        if (exists $h2{$blast4[0]." "."new1"." ".$fo}) {
                                                my @kk =split /\s+/,$h2{$blast4[0]." "."new1"." ".$fo};
                                                my @ll =split /\S+/,$h2{$blast4[0]." "."new1"." ".$fo};
                                                $h2{$blast4[0]." "."new2"." ".$blast4[3]}=$kk[0].$ll[1].$kk[1].$ll[2].$kk[2].$four1;
						delete $h2{$blast4[0]." "."new1"." ".$fo};
					} else {
                                                next;
                                        }
                                } elsif (($blast4[3]-$blast4[1])<=0) {
                                        my $four2 =join(" ",@blast4[2..3]);
                                        if (exists $h2{$blast4[0]." "."new1"." ".$ur}) {
						my @mm =split /\s+/,$h2{$blast4[0]." "."new1"." ".$ur};
                                                my @nn =split /\S+/,$h2{$blast4[0]." "."new1"." ".$ur};
                                                $h2{$blast4[0]." "."new2"." ".$blast4[3]}=$mm[0].$nn[1].$mm[1].$nn[2].$mm[2].$four2;
						delete $h2{$blast4[0]." "."new1"." ".$ur};
                                        } else {
                                                next;
                                        }
                                }
			} elsif ($j==6) {
				my @blast5=split/\s+/,$a[$i];
				if ($blast5[1]==0) {        #0 delete
                                           next;
                                }
                                my $fi=$blast5[1]-1;
                                my $ve=$blast5[1]+1;
                                if (($blast5[3]-$blast5[1])>=0) {
                                        my $five1 =join(" ",@blast5[2..3]);
                                        if (exists $h2{$blast5[0]." "."new2"." ".$fi}) {
                                                my @oo =split /\s+/,$h2{$blast5[0]." "."new2"." ".$fi};
                                                my @pp =split /\S+/,$h2{$blast5[0]." "."new2"." ".$fi};
                                                $h2{$blast5[0]." "."new2"." ".$fi}=$oo[0].$pp[1].$oo[1].$pp[2].$oo[2].$five1;
                                        } else {
                                                next;
                                        }
                                } elsif (($blast5[3]-$blast5[1])<=0) {
                                        my $five2 =join(" ",@blast5[2..3]);
                                        if (exists $h2{$blast5[0]." "."new2"." ".$ve}) {
                                                my @qq =split /\s+/,$h2{$blast5[0]." "."new2"." ".$ve};
                                                my @rr =split /\S+/,$h2{$blast5[0]." "."new2"." ".$ve};
                                                $h2{$blast5[0]." "."new2"." ".$ve}=$qq[0].$rr[1].$qq[1].$rr[2].$qq[2].$five2;
                                        } else {
                                                next;
                                        }
                                }
			}			
		}
		my @k = sort keys %h2;
		foreach my $k (@k) {
			my $k2 = (split /\s+/, $k)[0];
			if (exists $h1{$k2}) {
       				print  Jiewei "$h2{$k}\t$h1{$k2}\n";
			} else {
				print Jiewei "$h2{$k}\n";
			}
		}
		print Jiewei  "#####\n";
		$sequence="";
		$j=0;
		%h1=();
		%h2=();
	}
}

close Choose;
close Jiewei;

######### Calculate Similarity ######

open Jiewei,"<$output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei";
open Inden,">$output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei.inden";

my $seqone="";
my $jk=0;
my $m=0;
while (my $line=<Jiewei>) {
	if ($line !~ /#####/) {
		$seqone .=$line;
		next;
	} else {
		my @a=split /\n/,$seqone;
		my $len=0;
		for (my $i=0;$i<@a;$i++) {
			chomp($a[$i]);
			if ($jk==0) {
				if (length($a[$i])==0) {
					$jk++;
					print Inden "$a[$i]\n";
					next;
				}	
				print Inden "$a[$i]\n";
			} elsif ($jk==1) {
				if ($m==0) {
					$m++;
					my @b=split/\s+/,$a[$i];
					$len=length($b[2]);
					print Inden "identy $a[$i]\n";
				} else {
					my @c=split/\s+/,$a[$i];
					my @b2=split //,$c[2];
					my $dian=0;
					foreach my $xulie (@b2) {
						if ($xulie eq ".") {
							$dian++;
						}
					}
					my $identity=$dian/$len;
					$identity=sprintf "%.4f",$identity;
					$a[$i]=$identity." ".$a[$i];
					print Inden "$a[$i]\n";
				}
			}			
		}
		print Inden "#####\n";
		$jk=0;
		$m=0;
		$seqone="";
	}
}
close Jiewei;
close Inden;	

######### Calculate Mutation message ######

open Inden,"<$output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei.inden";
open Mutation,">$output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei.inden.mutation";

my %hash;
my $length;
my $start;
my $end;
my %hashbar;
my $llll;

while (my $line=<Inden>){
	chomp($line);
	if($line=~/^Query|^Sequences|^ |^$/){
		print Mutation "$line\n";
		next;
	}elsif($line=~/^###/){
		%hash=();
		%hashbar=();
		$length=0;
		$llll=0;
		$start=0;
		$end=0;
		print Mutation "$line\n";
		next;
	}
	if($line=~/identy/){
		print Mutation "Compute $line\n";
		my @temp=split /\s+/,$line;
		$start=$temp[2];
		$end=$temp[4];
		$llll=length($temp[3]);
#########store - ###############################
		my @tmp=split //,$temp[3];
		foreach my $n (0..$#tmp){
			if($tmp[$n] =~ /-/){
				my $key=$start+$n;
				$hashbar{$key}=1;
			}
		}
########################################
		my @temp1=split /\d+/,$line;
		my $seqandspace=$temp1[3];
		my @temp2=split /a|t|c|g/,$seqandspace;
		my @temp3=split //,$temp2[0];
		$length=length($start)+$#temp3+1;#important 
		my @temp4=split //,$temp[3];
		foreach my $num (0..$#temp4){#Actural position on the read
			my $keynew=$start+$num;
			$hash{$keynew}=$temp4[$num];
		}
	}
	my $change;
	my $identity;
	if($line=~/^1|^0/){
		my @temp=split /\s+/,$line;
#######################################
		my $count=0;
		my @t=split //,$temp[3];
		foreach my $le (@t){
			if($le eq "."){
				$count++;
			}
		}
		my $id=$count/$llll;	
		$identity=sprintf "%.4f",$id;
#######################################
		
		my $numlen=length($temp[2]);#length of the third colume
		my @temp1=split /\d+/,$line;
		my @temp2=split //,$temp1[4];
		my $startnum=$length-$numlen;
		foreach my $num (0..$#temp2){
			if($temp2[$num]=~/a|t|c|g/){
################# counting - ############################################
				my $newkey=$start+$num-$startnum;
				my $numbar=0;
				foreach my $nu ($start..$newkey){
					if(exists $hashbar{$nu}){
						$numbar++;
					}
				}
				my $truekey;
				if ($temp[4]-$temp[2]>0) {
					$truekey=$newkey-$numbar;
				} elsif ($temp[4]-$temp[2]<0) {
					$truekey=$end-($newkey-$numbar)+1;
				}
###########################################################	
				if($hash{$newkey}!~/-/){
					$change.="$end:$truekey:$hash{$newkey}:$temp2[$num] ";
				}
			}
		}
	}
	if(defined $change){
		print Mutation "$identity $line\tRL:MP:RB:MB\t$change\n";#RL: read length; MP: mutation position; RB: Raw base on the read; MB: base after mutation
	}else{
		if($line!~/^ide/ and defined $identity){
			print Mutation "$identity $line\n";
		}
	}
}

close Mutation;
close Inden;

######### Change To M8 Format ######

open Mutation,"<$output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei.inden.mutation";
open M8,">$output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei.inden.mutation.m8";

my $name;
my $lengthone;
my $startone;
my $endone;

while(my $line=<Mutation>){
	chomp($line);
	if($line=~/^Query/){
		my @temp=split /=| /,$line;
		$name=$temp[2];
		$lengthone=(split/\s+/,(split/length=/,$line)[1])[0];
	}
	if($line=~/^Compute/){
		my @temp=split /\s+/,$line;
		$startone=$temp[3];
		$endone=$temp[5];
	}
	if($line=~/^0|^1/){
		my @temp=split /\s+/,$line;
		my $similarity=$temp[0]*100;
		if ($line=~/RL:MP:RB:MB/) {
			my @xinxi=split/\t/,$line;
			print M8 "$name\t$temp[2]\t$similarity\t$lengthone\t$startone\t$endone\t$temp[3]\t$temp[5]\t$temp[6] $temp[7]\t$xinxi[2]\t$xinxi[3]\n";
		} else {
			print M8 "$name\t$temp[2]\t$similarity\t$lengthone\t$startone\t$endone\t$temp[3]\t$temp[5]\t$temp[6] $temp[7]\n";
		}
	}
}
close Mutation;
close M8;

########## CT mutation #######

if(defined $CT ) {
	system("awk \-F \'\\t\' \'\$10\~\/RL\:MP\:RB\:MB\/\' $output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei.inden.mutation.m8 \>$output/$sample.onlymutation.result");
	open OnlyMutation,"<$output/$sample.onlymutation.result";
	open OnlyMutationonly,">$output/$sample.onlymutation.result1";
	my %hashmutationuniq;
	while(my $line=<OnlyMutation>){
		chomp($line);
		my @temp=split /\t/,$line;
		my $key="$temp[0]###$temp[1]";
		if(exists $hashmutationuniq{$key}){
			my $num=(split/\t/,$hashmutationuniq{$key})[0];
			if($temp[2]=~/^\d/ and $num<=$temp[2]){
				$hashmutationuniq{$key}=$temp[2]."\t".$temp[3]."\t".$temp[4]."\t".$temp[5]."\t".$temp[6]."\t".$temp[7]."\t".$temp[8]."\t".$temp[9]."\t".$temp[10];
			}
		}else{
			if($temp[2]=~/^\d/){
				$hashmutationuniq{$key}=$temp[2]."\t".$temp[3]."\t".$temp[4]."\t".$temp[5]."\t".$temp[6]."\t".$temp[7]."\t".$temp[8]."\t".$temp[9]."\t".$temp[10];
			}
		}
	}
	close OnlyMutation;
	foreach my $key (keys %hashmutationuniq){
		my @temp=split /###/,$key;
		print OnlyMutationonly "$temp[0]\t$temp[1]\t$hashmutationuniq{$key}\n";
	}
	close OnlyMutationonly;
	open OnlyMutationonly,"<$output/$sample.onlymutation.result1";
	open MutationResult2,">$output/$sample.onlymutation.result2";
	while (my $line=<OnlyMutationonly>){
		chomp($line);
		my ($seq,$species,$inden,$length,$start,$end,$bstart,$estart,$specie,$zhushi,$mutation)=(split/\t/,$line)[0,1,2,3,4,5,6,7,8,9,10];
		my @allmutation=split/\s/,$mutation;
		my $n_5_c_t=0;my $n_5_g_a=0;my $n_3_c_t=0;my $n_3_g_a=0;my $mutationsite="";
		foreach my $each (@allmutation) {
			my @juti=split/\:/,$each;
			if ($each=~/t:c/ && ($juti[1]<=$CT)) {
				$n_5_c_t+=1;
				$mutationsite.="\t".$each;
			}
			if ($each=~/t:c/ && ($juti[1]>=($juti[0]-$CT+1)) && ($juti[1]<=$juti[0])) {
				$n_3_c_t+=1;
				$mutationsite.="\t".$each;
			}
			if ($each=~/a:g/ && ($juti[1]<=$CT)) {
				$n_5_g_a+=1;
				$mutationsite.="\t".$each;
			}
			if ($each=~/a:g/ && ($juti[1]>=($juti[0]-$CT+1)) && ($juti[1]<=$juti[0])) {
				$n_3_g_a+=1;
				$mutationsite.="\t".$each;
			}
		}
		print MutationResult2 "$seq\t$species\t$inden\t$length\t$start\t$end\t$bstart\t$estart\t$specie\t$n_5_c_t\t$n_3_c_t\t$n_5_g_a\t$n_3_g_a\t$mutationsite\n";
	}
	close OnlyMutationonly;
	close MutationResult2;
		if (defined $c_t_mutation &&  !defined $g_a_mutation) {
			system("awk \-F \'\\t\' \'\(\(\$6\-\$5\+1\)\/\$4\)\>\=$coveragemin \&\& \(\(\$6\-\$5\+1\)\/\$4\)\<=$coveragemax \&\&  \$3\>\=$similaritymin \&\& \$3\<\=$similaritymax \&\& \(\$10\>\=$c_t_mutation \|\| \$11\>\=$c_t_mutation\) \' $output/$sample.onlymutation.result2 \> $output/$sample.onlymutation.result2_ct");
			my $in="$output/$sample.onlymutation.result2_ct";
			my $out="$output/$sample.onlymutation.result2_ct2";
			&result1($in,$out);
			my $input_result2="$output/$sample.onlymutation.result2_ct2";
			my $output_result2_one="$output/$sample.onlymutation.result2_ct3";
			my $output_result2_two="$output/$sample.onlymutation.result2_ct4";
			&result2($input_result2,$output_result2_one,$output_result2_two);
			my $input_result3="$output/$sample.onlymutation.result2_ct4";
			my $output_result3="$output/$sample.Result.txt";
			my $output_result4="$output/Result.$sample.txt";
			&result3($input_result3,$output_result3,$output_result4);
			my $barfile="$sample.Result.txt";
			if (defined $Rpath && $output_result3){&barplot($output,$sample,$barfile);}
			&delete;
		} elsif (defined $g_a_mutation &&  !defined $c_t_mutation) {
			system("awk \-F \'\\t\' \'\(\(\$6\-\$5\+1\)\/\$4\)\>\=$coveragemin \&\& \(\(\$6\-\$5\+1\)\/\$4\)\<=$coveragemax \&\&  \$3\>\=$similaritymin \&\& \$3\<\=$similaritymax \&\& \(\$12\>\=$g_a_mutation \|\| \$13\>\=$g_a_mutation\) \' $output/$sample.onlymutation.result2 \> $output/$sample.onlymutation.result2_ga");
			my $in="$output/$sample.onlymutation.result2_ga";
			my $out="$output/$sample.onlymutation.result2_ga2";
			&result1($in,$out);
			my $input_result2="$output/$sample.onlymutation.result2_ga2";
			my $output_result2_one="$output/$sample.onlymutation.result2_ga3";
			my $output_result2_two="$output/$sample.onlymutation.result2_ga4";
			&result2($input_result2,$output_result2_one,$output_result2_two);
			my $input_result3="$output/$sample.onlymutation.result2_ga4";
			my $output_result3="$output/$sample.Result.txt";
			my $output_result4="$output/Result.$sample.txt";
			&result3($input_result3,$output_result3,$output_result4);
			my $barfile="$sample.Result.txt";
			if (defined $Rpath && $output_result3){&barplot($output,$sample,$barfile);}
			&delete;
		} elsif (defined $g_a_mutation && defined $c_t_mutation) {
			system("awk \-F \'\\t\' \'\(\(\$6\-\$5\)\/\$4\)\>\=$coveragemin \&\& \(\(\$6\-\$5\)\/\$4\)\<=$coveragemax \&\&  \$3\>\=$similaritymin \&\& \$3\<\=$similaritymax  \&\& \(\$10\>\=$c_t_mutation \|\| \$11\>\=$c_t_mutation \|\| \$12\>\=$g_a_mutation \|\| \$13\>\=$g_a_mutation\)\' $output/$sample.onlymutation.result2 \> $output/$sample.onlymutation.result2_ct_ga");
			my $in="$output/$sample.onlymutation.result2_ct_ga";
			my $out="$output/$sample.onlymutation.result2_ct_ga2";
			&result1($in,$out);
			my $input_result2="$output/$sample.onlymutation.result2_ct_ga2";
			my $output_result2_one="$output/$sample.onlymutation.result2_ct_ga3";
			my $output_result2_two="$output/$sample.onlymutation.result2_ct_ga4";
			&result2($input_result2,$output_result2_one,$output_result2_two);
			my $input_result3="$output/$sample.onlymutation.result2_ct_ga4";
			my $output_result3="$output/$sample.Result.txt";
			my $output_result4="$output/Result.$sample.txt";
			&result3($input_result3,$output_result3,$output_result4);
			my $barfile="$sample.Result.txt";
			if (defined $Rpath && $output_result3){&barplot($output,$sample,$barfile);}
			&delete;
		} else {
			print "if you set mutationsite,you must set ct_mutationnum or ga_mutationnum or you will get the wrong result,please try again\n";
			&delete;
			exit;
		}
}else {
	######### Get Important Message ######

	system"awk \-F \'\\t\' \'\{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9\}\' $output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei.inden.mutation.m8 \>$output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei.inden.mutation.m8.part";

	######### Get The Uniq Sequence ######

	open Part,"<$output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei.inden.mutation.m8.part";
	open Only,">$output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei.inden.mutation.m8.part.only";

	my %hashuniq;
	while(my $line=<Part>){
        	chomp($line);
		my @temp=split /\t/,$line;
		my $key="$temp[0]###$temp[1]";
		if(exists $hashuniq{$key}){
			my $num=(split/\t/,$hashuniq{$key})[0];
			if($temp[2]=~/^\d/ and $num<=$temp[2]){
				$hashuniq{$key}=$temp[2]."\t".$temp[3]."\t".$temp[4]."\t".$temp[5]."\t".$temp[6]."\t".$temp[7]."\t".$temp[8];
			}
		}else{
			if($temp[2]=~/^\d/){
				$hashuniq{$key}=$temp[2]."\t".$temp[3]."\t".$temp[4]."\t".$temp[5]."\t".$temp[6]."\t".$temp[7]."\t".$temp[8];
			}
		}
	}
	foreach my $key (keys %hashuniq){
        	my @temp=split /###/,$key;
		print Only "$temp[0]\t$temp[1]\t$hashuniq{$key}\n";
	}

	close Part;
	close Only;

	######### Select Coverage and Similarity ######
	my $coveragemin2=$coveragemin*100;
	my $coveragemax2=$coveragemax*100;

	if ((defined $coveragemin) && (defined $coveragemax) && (defined $similaritymin) && (defined $similaritymax)) {
		system"awk \-F \'\\t\' \'\(\(\$6\-\$5\+1\)\/\$4\)\>\=$coveragemin \&\& \(\(\$6\-\$5\+1\)\/\$4\)\<=$coveragemax \&\&  \$3\>\=$similaritymin \&\& \$3\<\=$similaritymax \' $output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei.inden.mutation.m8.part.only \> $output/$sample.coverage$coveragemin.$coveragemax.similarity$similaritymin.$similaritymax";
		system"mv $output/$sample.coverage$coveragemin.$coveragemax.similarity$similaritymin.$similaritymax $output/$sample.coverage$coveragemin2.$coveragemax2.similarity$similaritymin.$similaritymax"
	} else {
		die (print "Lack of coveragemin,coveragemax or similaritymin similaritymax\n");
	}	

	open Resultzero,"<$output/$sample.coverage$coveragemin2.$coveragemax2.similarity$similaritymin.$similaritymax";
	open Resultone,">$output/$sample.coverage$coveragemin2.$coveragemax2.similarity$similaritymin.$similaritymax.result1";


	my %hashresult1;

	while(my $line=<Resultzero>){
        	chomp($line);
		my @temp=split /\t/,$line;
		my $key="$temp[8]";
		$hashresult1{$key}+=1;
		}
	close Resultzero;
	foreach my $key (keys %hashresult1){
		print Resultone "$key\t$hashresult1{$key}\n";
	}

	close Resultone;

	######### Get Result ######

	open Resultone,"<$output/$sample.coverage$coveragemin2.$coveragemax2.similarity$similaritymin.$similaritymax.result1";

	my %hashresult2;
	my $zongnum=0;

	while (my $line=<Resultone>) {
        	chomp($line);
		my $num=(split/\t/,$line)[1];
		$zongnum+=$num;
	}
	$hashresult2{1}=$zongnum;

	close Resultone;

	open Resultone,"<$output/$sample.coverage$coveragemin2.$coveragemax2.similarity$similaritymin.$similaritymax.result1";
	open Resulttwo,">$output/$sample.coverage$coveragemin2.$coveragemax2.similarity$similaritymin.$similaritymax.result2";

	while (my $line=<Resultone>) {
		chomp($line);
		my $xinxi=(split/\t/,$line)[0];
		my $num=(split/\t/,$line)[1];
		my $bili=$num/($hashresult2{1});
		print Resulttwo "$xinxi\t$num\t$bili\t$hashresult2{1}\n";
	}
	system"sort \-t \$\'\\t\' \-r \-n \-k 2 $output/$sample.coverage$coveragemin2.$coveragemax2.similarity$similaritymin.$similaritymax.result2 \> $output/$sample.coverage$coveragemin2.$coveragemax2.similarity$similaritymin.$similaritymax.result3";

	close Resultone;
	close Resulttwo;


	open Resultthree,"<$output/$sample.coverage$coveragemin2.$coveragemax2.similarity$similaritymin.$similaritymax.result3";
	open Resultfour,">$output/$sample.Result.txt";

	my $rank=1;
#	print Resultfour "Result: Coverage form $coveragemin2% to $coveragemax2% and Similarity from $similaritymin% to $similaritymax%\n";
	print Resultfour "Species Ranking\tSpecies\tVMH\tPoVMH\tTotal hits\n";
	while (my $line=<Resultthree>) {
		chomp($line);
		print Resultfour "$rank\t$line\n";
		$rank++;
	}
	close Resultthree;
	close Resultfour;

	open Resultfour,"<$output/$sample.Result.txt";
	open Resultfive,">$output/Result.$sample.txt";
	my %hash1result;
	my $R;
	while (my $line2=<Resultfour>) {
		chomp($line2);
		my $ranknum=(split/\t/,$line2)[0];
		my $PoVMH=(split/\t/,$line2)[3];
		$hash1result{$ranknum}=$PoVMH;
		print Resultfive"$line2\n";
	}
	if (exists $hash1result{2}) {
		$R=$hash1result{1}/$hash1result{2};
	}else {
		$R="NA";
	}
	print Resultfive"Rvalue=$R\n";
	close Resultfive;
	close Resultfour;
	if (defined $Rpath) {
		system("awk \'\$1\!\~\/Result\:\/\' $output/$sample.Result.txt \| head -n 10 \> $output/$sample.coverage$coveragemin2.$coveragemax2.similarity$similaritymin.$similaritymax.result4.txt");
		my $barfile="$sample.coverage$coveragemin2.$coveragemax2.similarity$similaritymin.$similaritymax.result4.txt";
		&barplot($output,$sample,$barfile);
	}
	&delete;
}
#################### Creat R shell and plot###############

sub barplot {
	my ($output,$sample,$file)=@_;
	my $Rline=<<"Rline";
	library(ggplot2)
	data<-read.table("$output/$file",header = T,sep ="\\t")
	Rank<-data\$Species.Ranking
	Mappedratio<-data\$PoVMH
	Mappedcounts<-data\$VMH
	Species<-data\$Species
	pdf(file="$output/Result.$sample.Mappedratio_Rank.pdf",w=8,h=6)
	ggplot(data=data,mapping=aes(x=factor(Rank),y=Mappedratio,fill=Species))+
	geom_bar(stat="identity",width=0.3)+
	theme(axis.title.x =element_text(size=14),axis.text.x =element_text(size=12), axis.title.y=element_text(size=14),axis.text.y=element_text(size=12),legend.text= element_text(size=12))+xlab("Species Ranking")+ylab("PoVMH")
	dev.off()
	pdf(file="$output/Result.$sample.Mappedcounts_Rank.pdf",w=8,h=6)
	ggplot(data=data,mapping=aes(x=factor(Rank),y=Mappedcounts,fill=Species))+
	geom_bar(stat="identity",width=0.3)+
	theme(axis.title.x =element_text(size=14),axis.text.x =element_text(size=12), axis.title.y=element_text(size=14),axis.text.y=element_text(size=12),legend.text= element_text(size=12))+xlab("Species Ranking")+ylab("VMH")
	dev.off()
Rline
	open (ROUT,">$output/$sample.Rplot.shell");
	print ROUT "$Rline";
	close (ROUT);
	system("$Rpath CMD BATCH $output/$sample.Rplot.shell");
}

######### Decide Whether To Delete The Intermediate File ######

sub delete {
	my $panduan1="yes";
	my $panduan2="no";
	if ($panduan1 =~/^$delete$/i) {
	system"rm $output/$sample\*";

#	system"rm $output/$sample.filterfq.fq $output/$sample.filterfasta.fa $output/$sample.adlengthfasta.fa $output/$sample.adlengthfasta.fa.out $output/$sample.adlengthfasta.fa.out.jiablast $output/$sample.adlengthfasta.fa.out.jiablast.choose $output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei $output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei.inden $output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei.inden.mutation $output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei.inden.mutation.m8 $output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei.inden.mutation.m8.part $output/$sample.adlengthfasta.fa.out.jiablast.choose.jiewei.inden.mutation.m8.part.only $output/$sample.coverage$coveragemin2.$coveragemax2.similarity$similaritymin.$similaritymax $output/$sample.coverage$coveragemin2.$coveragemax2.similarity$similaritymin.$similaritymax.result1  $output/$sample.coverage$coveragemin2.$coveragemax2.similarity$similaritymin.$similaritymax.result2 $output/$sample.coverage$coveragemin2.$coveragemax2.similarity$similaritymin.$similaritymax.result3";
	} elsif ($panduan2 =~/^$delete$/i) {
	} else {
	}
}
################## Result Format Sub ##############

sub result1 {
	my ($input,$output)=@_;
	my %hashresult1;
	open Resultzero,"<$input";
	while(my $line=<Resultzero>){
	        chomp($line);
		my @temp=split /\t/,$line;
		my $key="$temp[8]";
		$hashresult1{$key}+=1;
		}
	close Resultzero;
	open Resultone,">$output";
	foreach my $key (keys %hashresult1){
		print Resultone "$key\t$hashresult1{$key}\n";
	}
	close Resultone;
}

sub result2 {
	my ($input,$output,$output2)=@_;
	open Resultone,"<$input";
	my %hashresult2;
	my $zongnum=0;
while (my $line=<Resultone>) {
	chomp($line);
	my $num=(split/\t/,$line)[1];
	$zongnum+=$num;
}
$hashresult2{1}=$zongnum;
close Resultone;
open Resultone,"<$input";
open Resulttwo,">$output";
while (my $line=<Resultone>) {
	chomp($line);
	my $xinxi=(split/\t/,$line)[0];
	my $num=(split/\t/,$line)[1];
	my $bili=$num/($hashresult2{1});
	print Resulttwo "$xinxi\t$num\t$bili\t$hashresult2{1}\n";
}
close Resultone;
close Resulttwo;
system"sort \-t \$\'\\t\' \-r \-n \-k 2 $output \> $output2";
}

sub result3 {
	my ($input,$output,$output2)=@_;
	open Resultthree,"<$input";
	open Resultfour,">$output";
	my $rank=1;
	print Resultfour "Species Ranking\tSpecies\tVMH\tPoVMH\tTotal hits\n";
	while (my $line=<Resultthree>) {
		chomp($line);
		print Resultfour "$rank\t$line\n";
		$rank++;
	}
	close Resultthree;
	close Resultfour;
	open Resultfour,"<$output";
	open Resultfive,">$output2";
	my %hash1;
	my $R;
	while (my $line2=<Resultfour>) {
		chomp($line2);
		my $ranknum=(split/\t/,$line2)[0];
		my $PoVMH=(split/\t/,$line2)[3];
		$hash1{$ranknum}=$PoVMH;
		print Resultfive"$line2\n";
	}
	if (exists $hash1{2}) {
		$R=$hash1{1}/$hash1{2};
	}else {
		$R="NA";
	}
	print Resultfive"Rvalue=$R\n";

			
close Resultfive;
close Resultfour;
}











