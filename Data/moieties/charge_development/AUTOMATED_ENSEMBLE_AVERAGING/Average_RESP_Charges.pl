#!/usr/bin/perl
# Version 1.05 : 2009-04-08
# Usage:			ARGV[0]		ARGV[1]	ARGV[2]
#	Average_RESP_Charges.pl file_prefix output_file num_atoms
#
#use Data::Types qw(:all);
use strict;
use warnings;
use Cwd;

use File::Find;

my $file_pattern  =$ARGV[0];
my $output_name=$ARGV[1];
my $num_atoms = $ARGV[2];
my $file_extension = ".log.resp.pch";
my $search_pattern = "Point charges before & after optimization";
my $i=0;
my $j=0;
my $k=0;
my $LINE;
my $temp=0;
my (@RAWARRAY,@DATA,@inputATOMID);
my @ATOMLIST = ("e-","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Me","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn"); # Elements by element number

print "\nGenerating Ensemble Averaged Charges...\n";
find(\&d, cwd);

sub d {
  my $file = $File::Find::name;

#  $file =~ s,/,\\,g;

#  return unless -f $file;
#  print "$file\n";
  return  unless $file =~ m/$file_pattern/;
  return  unless $file =~ m/$file_extension/;
#  print $file." j number is ".$j."\n";
  open (INFILE, $file) or print "couldn't open $file\n" && return;
  my @LINES = <INFILE>;
  close(INFILE);
#  print $search_pattern." is the searchable patter\n";
  foreach $LINE (@LINES){
  	if($LINE =~ /$search_pattern/){
#		print $LINE." k is ".$k."\n";
		$k=1;
#		print "now k is ".$k."\n";
	}
	if($k>2){
		if($j==0){
	                $temp=substr($LINE,10,2);
			$temp=~ s/ //g;
	                $inputATOMID[$i]=$ATOMLIST[$temp];
#			print $inputATOMID[$i]."\n";
		}
		my $temp2=substr($LINE,29,10);
		$temp2=~ s/ //g;
		$temp = ("$temp2");
		$RAWARRAY[$i][$j]=$temp;
#		print $RAWARRAY[$i][$j]."\n";
		$i++;
	}
	if($i==($num_atoms)){
	$i=0;
	$j++;
	$k=0;
	return;
	}
	if($k>0){ $k++; }
  } 
}
print "Number of files ".$j."\nNumber of atoms ".($#RAWARRAY+1)."\n";

sub stdev {
        my ($a,$b,$c); ### (netcharge of values squared, population(n), average value)
        ($a,$b,$c) = ($_[0],$_[1],$_[2]);
	my $d = (($a/$b)**(1/2));
	if($a > 0 && $b > 0 && (sprintf("% .3f",$c) != sprintf("% .3f",$d))){ # Check for positive squared values, positive population, and non-zero standard deviation
	        return (sqrt(( $a - ( $b * (  $c ** 2 ))) / ($b - 1)));
	}
	else{ return(0); }
}

my $no_files = $j;
my $netcharge = 0;
my $poscharge = 0;
my $negcharge = 0;

open (OUTPUT,">$output_name");
for($i=0;$i<=$#RAWARRAY;$i++){ # Scrolls through each atom's charge set
	my $avg = 0;
	my $sq = 0;
	for($j=0;$j<$no_files;$j++){ # Scrolls through each charge value
		$temp = ( ($avg * $j) + $RAWARRAY[$i][$j] ) / ($j + 1);
		$avg = $temp;
		$temp = $sq + ($RAWARRAY[$i][$j] ** 2);
		$sq = $temp;
	}
	$DATA[$i][0] = sprintf("%.3f",$avg);
	$DATA[$i][1] = sprintf("%.3f",&stdev($sq,$no_files,$avg));
	print OUTPUT sprintf("%4s",($i+1))."   ".$inputATOMID[$i]."    Charge:  ".sprintf("% .3f",$DATA[$i][0])."    Std Dev:  ".sprintf("% .3f",$DATA[$i][1])."\n";
	$temp = $netcharge + $DATA[$i][0];
        $netcharge = $temp;
	if($DATA[$i][0]>0){
		$temp = $poscharge + $DATA[$i][0];
		$poscharge = $temp;
	}
	if($DATA[$i][0]<0){
                $temp = $negcharge + $DATA[$i][0];
                $negcharge = $temp;
	}
}
print OUTPUT "Net Charge: ".sprintf("% .3f",$netcharge)."\n";
print OUTPUT "\n# Charge Magnitude #\nPos Charge: ".sprintf("% .3f",$poscharge)."\n";
print OUTPUT "Neg Charge: ".sprintf("% .3f",$negcharge)."\n";

if (sprintf("%.3f",$netcharge) == sprintf("%.0f",$netcharge)){
	print "Net Charge: ".sprintf("% .3f",$netcharge)."\n";
}
else{
	print "\n!!!!!!!!\nWARNING: Net Charge is not an integer: ".sprintf("% .3f",$netcharge)."\nCorrect before adding to PREP/MOL2 file\n!!!!!!!!\n";
}
close (OUTPUT);
print "\n";
