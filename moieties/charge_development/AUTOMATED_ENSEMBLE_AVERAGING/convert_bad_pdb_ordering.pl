#!/usr/bin/perl
############# User adjusted input files goes here ##################
# The prefix for the copied and transfered files
my $prefix = "0IB-OME_ensembleMD";

# The Template file name
my $template = "template.pdb";   

# The number of lines in the pdb file
my $num_lines = 40;

################# Working Portion of the Script ###############

  for($x=1;$x<=100;$x=$x+1){ 
          
 	$datafile=$prefix.".pdb.".$x;
 	open(TEST,"$datafile") || die "couldn't open $datafile";
		
	open(MODCRDS,">".$prefix.$x.".pdb") || die "oops: $!";
		
        #get atom and residue info
        open(OUTPUT,"$template") || die "not open: $!";
	my @pdb_array;
        for($j=1;$j<=$num_lines;$j++){
      	 	chomp($partfile=<OUTPUT>);
        	chomp($dataline=<TEST>);
                                                                                                                             
        #$num gets the coordinates and $num2 gets res and atom info
		$test=substr($partfile,0,6);
		if($test eq "HETATM" or $test eq "ATOM  "){
	              $num3=substr($partfile,6,5);
		      $num2=substr($partfile,11,21);
        	      $num= substr($dataline,32,34);
	              $pdb_array[$num3]="$test"."$num3"."$num2"."$num";

		}
	      #concatenate $num2 and $num and print to new file      	
       }
        for($j=1;$j<@pdb_array;$j++){
		print MODCRDS $pdb_array[$j]."\n";
	}
       close(TEST);
       close(OUTPUT);
       close(MODCRDS);

}
