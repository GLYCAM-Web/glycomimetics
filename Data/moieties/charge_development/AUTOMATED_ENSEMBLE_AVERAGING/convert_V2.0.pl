#!/usr/bin/perl
#
#                   @ARGV[0]  @ARGV[1] @ARGV[2] @ARGV[3]
# usage : convert.pl beginnum endnum interval prefix_name
#
#
#


for($x=@ARGV[0];$x<=@ARGV[1];$x=@ARGV[2]+$x)
{


   $datafile = @ARGV[3].".pdb.".$x ;


   open(DATA,"$datafile") || die "Can't open standard PDB file : $datafile\n";

   $outfile = @ARGV[3]."_".$x.".pdb";
   open(OUT,">>$outfile") || die "Can't open output file : $outfile\n";

        while($line=<DATA>){
        print OUT $line
}
      print OUT TER;
   close(DATA);
   close(OUT);


}

