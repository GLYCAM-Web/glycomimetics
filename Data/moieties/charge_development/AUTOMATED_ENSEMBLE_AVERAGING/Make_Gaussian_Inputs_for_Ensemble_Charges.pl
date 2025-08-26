#!/usr/bin/perl
# Version 1.0.0.1
# 2011-7-19 M.B. Tessier
# Converts, builds and prepares Gaussian input files from Ensemble MD simulations to be used for RESP charges
#_____________________________________________________________________________________________________#


#___# Usage #____________________________________#
#                                   ARGV[0]
# ./Ensemble_Charges.pl Settings_Input_File
#________________________________________________#


#___# REQUIREMENTS #_________________________________________________________________________________________________#
#	1. The bonds and angles will come directly from the template file.  The torsions that are frozen and non-
#		frozen-copied will be the only values not copied from this file.

#	2. The template and input files must have the same file format and atom-order.  The dihedral values in the
#		z-matrix may not be setup through the same atoms, so an option exists to use PTRAJ or VMD torsion
#		outputs (the torsions needs to be in the same sequential order as the files, timestep values will be
#		ignored).

#	3. If the torsion reassignment option is chosen, files must be named with the torsion mask (ex. d11) followed
#		by the '.dat' extension.  So the file 'd11.dat' must be in the working directory.  If it is not present,
#		an output warning will be issued and the value of 'd11' will be extracted from the input files.

#	4. The dihedral (torsion) identifiers in the input files must be unique so a search can extract the values.
#		For example searching for 'd11' will reveal at least two entries in a g94 file but 'd11=' will be
#		specific to the torsion value.  The torsion mask below automatically adds the equals sign unless the
#		variable is altered within the code.

#	5. The non-frozen and frozen torsions can contain the same torsion identifiers but the default action will be
#		to freeze any torsion found in both masks.

#____________________________________________________________________________________________________________________#



#______# Input Control File Keywords  #_______________________________________#
#TEMPLATE_GAUSS_FILE		 template_OPT.com
#PARTIAL_TEMPLATE_END		 0
#INT_START                       1
#INT_FINISH                      100
#INT_FREQ			 1
#INPUT_FILE_PREFIX               0IA-OME_ensembleMD
#INPUT_FILE_EXTENSION            .pdb.com
#OUTPUT_FILE_PREFIX              0IA-OME_ensMDgauss
#FROZEN_COPIED_TORSIONS          d5,d8,d12,d15,d26,d27
#NONFROZEN_COPIED_TORSIONS       d16,d17,d18,d24,d25
#REASSIGN_TORSIONS               YES
#GAUSSIAN_MEM                    128MW
#GAUSSIAN_QM_OPT_COM             hf/6-31g(d)
#GAUSSIAN_QM_RESP		 hf/6-31g(d)
#GAUSSIAN_COMP_SYS               pcluster
#GAUSSIAN_NUM_PROCS              4
#_____________________________________________________________________________#


#____________________# Variable Descriptions #______________________________________________________________________#
#TEMPLATE_GAUSS_FILE             The template Gaussian input file that contains the optimized bond, angle and torsion
#                               values.  This file will serve as the template for all torsions that are not replaced
#				from the FROZEN_COPIED_TORSIONS and NONFROZEN_COPIED_TORSIONS masks.  All bonds and
#				and angles will be assigned from the geometries in this file.
#				
#				If the file DOES NOT exist, i.e. NULL, then it is assumed you don't want to used a
#				template and will instead use the original files as templates.  Essentially this just
#				allows you to freeze torsions.  Non-frozen torsions definitions should ONLY be used
#				from external files if the bond, angle and torsion assignments are all equivalent in
#				each file.
#______________________________________________________________________
#PARTIAL_TEMPLATE_END
#				 This allows selection of a contiguous (from 1 to ##)list of atoms which will be used
#				from the template.  For now it assumes that from position 1 to # will come from the
#				template.  A zero (0) value will mean that the entire template is used.
#______________________________________________________________________
#INT_START                       The integer value for the lowest filename value (INT)
#				 The filenames should be standard to match the 'prefixINT.extension' naming system
#______________________________________________________________________
#INT_FINISH                      The integer value for the highest filename value (INT)
#                                The filenames should be standard to match the 'prefixINT.extension' naming system
#______________________________________________________________________
#INT_FREQ                        The integer value for the frequency of values between the START and FINISH filenames
#______________________________________________________________________
#INPUT_FILE_PREFIX               The 'prefix' for the input filenames
#                                The filenames should be standard to match the 'prefixINT.extension' naming system
#______________________________________________________________________
#INPUT_FILE_EXTENSION            The 'extension' for the input filenames
#                                The filenames should be standard to match the 'prefixINT.extension' naming system
#______________________________________________________________________
#OUTPUT_FILE_PREFIX              The 'prefix' for the gaussian output ('gauss_ext' is 'com') filenames
#                                The filenames should be standard to match the 'prefixINT.extension' naming system
#______________________________________________________________________
#FROZEN_COPIED_TORSIONS          The frozen torsion(s) to be copied from the input file separated by commas.  The
#				program	assumes that the torsion values will be found in the INPUT files by searching
#				for string input found here followed by an equal sign.  This construct can be changed
#				by changing the $search_motif variable.
#______________________________________________________________________
#NONFROZEN_COPIED_TORSIONS       The non-frozen torsion(s) to be copied from the input file separated by commas.
#				These torsions are rotatable and vital to proper geometry but will not be frozen for
#				optimization purposes.  Like the FROZEN_COPIED_TORSIONS variable(s), these values are
#				searched for by the string provided followed by an equal sign.  This construct is tied
#				to the same $search_motif variable from FROZEN_COPIED_TORSIONS.
#______________________________________________________________________
#REASSIGN_TORSIONS              The YES/NO option to reassign the torsions copied from the input files from a PTRAJ or VMD
#                               output file (two column format;first is timestep, second is torsion value). The
#                               expected filenames will be the torsion value specified in the FROZEN and NONFROZEN
#                               sections followed by a '.dat' extension.  Example: d14.dat
#______________________________________________________________________
#GAUSSIAN_MEM                    The Gaussian %mem variable value.  This is a string that should contain the integer
#				memory value immediately followed by the 'MW' unit. Example: 128MW
#______________________________________________________________________
#GAUSSIAN_QM_OPT_COM             The Gaussian theory level used for geometry optimization. Example: hf/6-31g(d)
#______________________________________________________________________
#GAUSSIAN_QM_RESP                The Gaussian theory level used for RESP charge generation. Example: hf/6-31g(d)
#______________________________________________________________________
#GAUSSIAN_COMP_SYS               The computer system where Gaussian will be run.
#				Available systems:
#					pcluster	The UGA PCluster system (uses %NprocLinda=1
#							and %NprocShared=GAUSSIAN_NUM_PROCS)
#					generic		Any generic system - any value will be used that isn't
#							recognized as a specific system listed above (uses 
#							%nproc=GAUSSIAN_NUM_PROCS)
#______________________________________________________________________
#GAUSSIAN_NUM_PROCS              The number of processors the Gaussian job will run on.  Must be an integer of 1 or
#				more.  Default is 1 processor.

#_____________________________________________________________________________________________________________________#


use strict;
use Cwd;
use File::Copy;
use Tie::File;

###################### Assigned Variables #############################################################################
my $dir = getcwd;
my (@frozenTor,@nonFroTor,@array,@tempary,@FORMATEDTEMP);
my (%LOG_FILES);
my ($i,$j,$x,$y)=0;
# The extension for the torsion files $search_motif
my $search_motif="dat";
# The output file Gaussian98/03 extension definition
my $gauss_ext="com";
my ($backup_dir) = "backup";
##################### End of variables ################################################################################

print "\n";

##################### Parsing the input control file ##################################################################
open(SETTINGS_FILE,@ARGV[0]) || die "Cannot open the settings file\nSettings file requirements are listed at the top of this program.\nNext time read between the lines.\n      .--.\n   .'';  |.-.\n   |  |  |  |\n   |  |  |  |\n   |  |  |  |\n   |  |  |  |_ \n   |  | _|  / `,\n   }  /``) /  /\n   |`/   /:__/ \\\n   |/   /      |\n   (   '\\      |\n    \\    `.   /\n     |       |\n     |       |\nHave a nice day!\n";

# Build the hash table for the settings file
our %hashtable_input=();
our @keys_input =();
our @values_input=();
    while(my $read_input = <SETTINGS_FILE>){
       if($read_input=~ /^\s*(\S+)\s+(.*)$/)#Matching the config values
            {
               my $first_input = $1;#extracting the key
               my $second_input = $2;#extracting the value
                if( ($second_input =~ /^\S*\s*$/)|| ($second_input =~ /^\S.*\s+$/))
                 {
                   $second_input =~ s/\s+$//;#removing the space if present in the value part
                  }   ################# End of "If" statement
                       #Creating the hash table with the help of key & value
                        $hashtable_input{$first_input} = "$second_input";
               } ##################### End of "If" statment
         }########### Closing "While" for input file

# Assign the variables
our $templateFile  = $hashtable_input{'TEMPLATE_GAUSS_FILE'};
our $templateEndMask = $hashtable_input{'PARTIAL_TEMPLATE_END'};
our $nStart  = $hashtable_input{'INT_START'};
our $nStop  = $hashtable_input{'INT_FINISH'};
if ($nStart > $nStop) { die "ERROR: The INT_START integer, $nStart, must be smaller than the INT_FINISH integer, $nStop \nFix your input settings file before continuing\n"; }
if ($nStart =~ /\D/ || $nStop =~ /\D/) { die "ERROR: Expected an interger for both INT_START and INT_STOP\n";}
our $nFreq  = $hashtable_input{'INT_FREQ'};
if ($nFreq =~ /\D/) { die "ERROR: Expected an interger for INT_FREQ\n";}
$j = ($nStop - $nStart + 1) / $nFreq;
if ( $j =~ /\D/) {die "ERROR: There is not an integer number of points between $nStart and $nStop using $nFreq steps\n";}
our $inPrefix  = $hashtable_input{'INPUT_FILE_PREFIX'};
our $inExt  = $hashtable_input{'INPUT_FILE_EXTENSION'};
our $outPrefix  = $hashtable_input{'OUTPUT_FILE_PREFIX'};
our $frozenTors = $hashtable_input{'FROZEN_COPIED_TORSIONS'};
our $nonFroTors = $hashtable_input{'NONFROZEN_COPIED_TORSIONS'};
our $reassignTors = $hashtable_input{'REASSIGN_TORSIONS'};
$reassignTors =~ tr/a-z/A-Z/;
our $gauMem = $hashtable_input{'GAUSSIAN_MEM'};
our $gauOptTheory = $hashtable_input{'GAUSSIAN_QM_OPT_COM'};
our $gauRespTheory = $hashtable_input{'GAUSSIAN_QM_RESP'};
our $gauCompSys = $hashtable_input{'GAUSSIAN_COMP_SYS'};
$gauCompSys =~ tr/a-z/A-Z/;
our $gauProcs = $hashtable_input{'GAUSSIAN_NUM_PROCS'};
our @masterArray;

if($templateEndMask eq '') { $templateEndMask = 0;}
# Diagnostic output

## Create PCLUSTER Gaussian98/03 header ##########################################
sub gauHeaderPcluster{
        return "%NprocLinda=1\n%NprocShared=".$gauProcs."\n";
}#### End of create PCLUSTER Gaussian98/03 header ################################
## Create generic Gaussian98/03 header ###########################################
sub gauHeaderGeneric{
        return "%nproc=".$gauProcs."\n";
}#### End of create generic Gaussian98/03 header #################################

# Template File setup subroutine
sub formTemplate{
	my ($templateStart, $endName, $templateEnd) = @_;
	my (@LINES,@sLINES,@eLINES); # The final output array
        open(TEMPLATESTART,"<$templateStart") || die "ERROR: Could not open $templateEnd. $!";
        @sLINES = <TEMPLATESTART>; # The starting part of the Gaussian input file
        close(TEMPLATESTART);

	if($endName!=0 or $endName ne "0"){	
		open(TEMPLATEEND,"<$templateEnd") || die "ERROR: Could not open $templateEnd. $!";
		@eLINES = <TEMPLATEEND>; # The end of the Gaussian input file
		close(TEMPLATEEND);
	}

        $y=0;
        while($y==0){ # This loop runs through the template file array and removes the template header information
                if(substr($sLINES[0],0,1) eq "#") {$y = 1;}
                shift(@sLINES);
        }
        if ($#LINES == 0) { die "ERROR: The template is not a standard Gaussian input file with a control line beginning with the pound (#) sign\n"; }

	$y=0;
	my $z=0;
	if($endName !=0 or $endName ne "0"){
		while($y==0){ # This loop runs through the template file array and removes everything after finding the $endName
			$LINES[$z]=$sLINES[$z];
			if(substr($sLINES[$z],0,length($endName)) eq "$endName"){
				$y=1;
			}
			$z++;
			
		}
		my $yi=0;
		$y=0;
		for($yi=0;$yi<$#eLINES;$yi++){
			if ($y==1){
				$LINES[$z]=$eLINES[$yi];
				$z++;
			}
			if(substr($eLINES[$yi],0,length($endName)) eq "$endName"){
				$y=1;
			}
		}
		$LINES[$z++]="\n";
	}
	else{ @LINES=@sLINES;}	

        if($gauCompSys eq "PCLUSTER"){ $x=&gauHeaderPcluster(); } # Need to tell these two subroutines which file to process
        else { $x=&gauHeaderGeneric(); }
#       print $x;
        unshift(@LINES, "# opt=z-matrix ".$gauOptTheory." scf=tight\n");
        unshift(@LINES, "$x");
        unshift(@LINES, "%mem=$gauMem\n");
#       foreach (@LINES){print $_;}
        return @LINES; 
}#### End of template reformat ###################################################

if (-e "$templateFile"){
	print "##### Input Parameters #####\n";
        print "Template File:          ".$templateFile."\nInput files:            ".$inPrefix." ".$nStart." to ".$nStop." every ".$nFreq." ".$inExt."\n";
        print "Output files:           ".$outPrefix." ".$nStart." to ".$nStop." every ".$nFreq." .com\nFrozen torsions:        ".$frozenTors."\n";
# The Frozen torsion array diagnostic
#       print "Index    Torsion\n";
#       for($i=0;$i<=$#fTorLeng;$i++){
#               print "$i         $frozenTor[$i]\n";
#       }
        print "Non-frozen torsions:    ".$nonFroTors."\n";
# The non-Frozen torsion array diagnostic
#        print "Index    Torsion\n";
#        for($i=0;$i<=$#nFTorLeng;$i++){
#                print "$i         $nonFroTor[$i]\n";
#        }
        print "Torsion reassignment?:  ".$reassignTors."\n";
        print "Gaussian parameters:\n Memory: ".$gauMem."\n OPT Theory: ".$gauOptTheory."\n RESP Theory: ".$gauRespTheory."\n Processors: ".$gauProcs."\n System: ".$gauCompSys."\n";
	print "############################\n\n";
	if ($templateEndMask eq "0"){
		print "The entire Template File, $templateFile, will be used.\n";
		@masterArray = &formTemplate("$templateFile","$templateEndMask","0");
	}
        if (-e "$templateFile" && $templateEndMask ne "0"){
	        print "Part of the Template File, $templateFile, will be used until $templateEndMask is reached and then the target structure is used.\n";
	}
}
else { 
#	die "ERROR: The template file named ".$templateFile." does not exist\n";
        print "##### Input Parameters #####\n";
        print "IMPORTANT WARNING: INPUT FILES will be as the Template because the template file, $templateFile, doesn't exist\nInput files:            ".$inPrefix." ".$nStart." to ".$nStop." every ".$nFreq." ".$inExt."\n";
        print "Output files:           ".$outPrefix." ".$nStart." to ".$nStop." every ".$nFreq." .com\nFrozen torsions:        ".$frozenTors."\n";
# The Frozen torsion array diagnostic
#       print "Index    Torsion\n";
#       for($i=0;$i<=$#fTorLeng;$i++){
#               print "$i         $frozenTor[$i]\n";
#       }
        print "Non-frozen torsions:    ".$nonFroTors."\n";
# The non-Frozen torsion array diagnostic
#        print "Index    Torsion\n";
#        for($i=0;$i<=$#nFTorLeng;$i++){
#                print "$i         $nonFroTor[$i]\n";
#        }
        print "Torsion reassignment?:  ".$reassignTors."\n";
        print "Gaussian parameters:\n Memory: ".$gauMem."\n OPT Theory: ".$gauOptTheory."\n RESP Theory: ".$gauRespTheory."\n Processors: ".$gauProcs."\n System: ".$gauCompSys."\n";
        print "############################\n\n";

}


# Placing the Frozen Torsions mask into @frozenTor
my @tempary = split('\,',$frozenTors);
my $fTorLeng = $#tempary;
for($i=0;$i<=$fTorLeng;$i++){
	$frozenTor[$i][0] = $tempary[$i];
	$frozenTor[$i][1] = 0;
#        print "Testing f 0 ".$frozenTor[$i][0]."\n";
#        print "Testing f 1 ".$frozenTor[$i][1]."\n";
}
undef @tempary;
# Placing the non-Frozen Torsions mask into @nonFroTor
my @tempary = split('\,',$nonFroTors);
my $nFTorLeng = $#tempary;
for($i=0;$i<=$nFTorLeng;$i++){
        $nonFroTor[$i][0] = $tempary[$i];
        $nonFroTor[$i][1] = 0;
#        print "Testing nf 0 ".$nonFroTor[$i][0]."\n";
#        print "Testing nf 1 ".$nonFroTor[$i][1]."\n";
}
undef @tempary;

# Check for the value of the length of the Torsion arrays
#print $fTorLeng."\n";
#print $nFTorLeng."\n";

$j = 0;

# Check for the Frozen Torsion Reassignment files
if ($reassignTors eq 'YES'){
	print "##### Torsion reassignment #####\n";
	for($i=0;$i<=$fTorLeng;$i++){
		if (-e "$frozenTor[$i][0].$search_motif"){
			$frozenTor[$i][1] = 1;
			print "Frozen torsion ".$frozenTor[$i][0]." will be replaced using ".$frozenTor[$i][0].".dat\n";
#			print $frozenTor[$i][1]."\n";
			$j = 1;
		}
		else {print "Frozen torsion ".$frozenTor[$i][0]." will be replaced using the value in the INPUT files\n";}
	}
        print "\n";
}

# Check for the Non-Frozen Torsion Reassignment files
if ($reassignTors eq 'YES'){
        for($i=0;$i<=$nFTorLeng;$i++){
                if (-e "$nonFroTor[$i][0].$search_motif"){
                        $nonFroTor[$i][1] = 1;
                        print "Relaxed torsion ".$nonFroTor[$i][0]." will be replaced using ".$nonFroTor[$i][0].".dat\n";
#                       print $nonFroTor[$i][1]."\n";
			$j = 1;
                }
		else { print "Relaxed torsion ".$nonFroTor[$i][0]." will be replaced using the value in the INPUT files\n";}
        }
	if($j==0){ print "WARNING: Torsion Reassignment selected but reassignment files are not found for the format [torsion].dat\nProceeding with processing using the values in the input files...\n"; }
	print "################################\n";
}

############################ End of parsing #############################################################


##################### Subroutines #######################################################################
#
## Subroutine to extract a line from a file based on the line number ####################
sub extract {
   my ($filename, $line_no)=@_;
   my $line;
   open (FILE, $filename) || die "$filename can't be opened $! ";
   if ($line_no =~ /\D/) {
        while ($line=<FILE>) {
            if ($line =~ /$line_no/) {
                return $line;
            }
        }
   }
   else {
     foreach (1..$line_no) {
           $line = <FILE>;
        }
        return $line;
   }
}
#### End of the extract subroutine ##############################################

## Find and extract line by string name #########################################
sub findExtract{
	my ($filename,$search)=@_;
	open(TORFILE,"<$filename") || die "ERROR: Could not open: ".$filename."\n$!";
	my @LINES = <TORFILE>;
	close(TORFILE);
	my ($LINE);
	foreach $LINE (@LINES){
		if($LINE =~ /$search/){
#			print "Found $search in the following line:\n$LINE";
			return($LINE);
		}
	}
	die "ERROR: I have not found the pattern: ".$search."\n";
}#### End of find and extract line by string name ################################


## Strip the timestep from the special input file torsion value ##################
sub formatTor{ # Input format is "dihedral_name=","timestep dihedral_value"
	my ($dihName,$fixThis)=@_;
	my @tempary = split(/\s+/,$fixThis);
	my $fixedDih = $dihName." ".$tempary[$#tempary]."\n";
	return($fixedDih);
}### End of reformatting torsion input from special input file ###################


## Find torsion index in template array ##########################################
sub findIndex{
	my ($search) = @_;
	for($y=0;$y<=$#FORMATEDTEMP;$y++){
                if($FORMATEDTEMP[$y] =~ /$search/){
#                        print "Found $search in the following line:\n$FORMATEDTEMP[$y]";
                        return($y);
                }
        }
	die "ERROR: Could not find search pattern ".$search." in the template file\n";
}#### End of non-frozen torsion replacement ######################################


## Delete frozen torsion values from OUTPUT file and add frozen torsion to the end of the file ##
sub processFrozenTor{

}#### End of frozen torsion processing ####
#
############################ End of subroutines #########################################################


##################### File processing ###################################################################
# Build the correct header into the template
for($i=$nStart;$i<=$nStop;$i++){ # Run through the output files
	if (-e "$templateFile" && $templateEndMask ne "0"){
		@masterArray = &formTemplate("$templateFile","$templateEndMask","$inPrefix$i$inExt");
	}
	unless (-e "$templateFile"){
		@masterArray = &formTemplate("$inPrefix$i$inExt","0","0");
	}
	@FORMATEDTEMP = @masterArray;
#	print $#FORMATEDTEMP."\n";
	for($j=0;$j<=$nFTorLeng;$j++){
		if ($nonFroTor[$j][1]==1){
			$x=&extract("$nonFroTor[$j][0].$search_motif","$i");
			$nonFroTor[$j][2]=&formatTor($nonFroTor[$j][0]."=",$x);
		}
		else {$nonFroTor[$j][2]=&findExtract("$inPrefix$i$inExt",$nonFroTor[$j][0]."=");}
#		print "Torsion: ".$nonFroTor[$j][0].", File indicator: ".$nonFroTor[$j][1].", Line replacement value: ".$nonFroTor[$j][2];
		$x = &findIndex($nonFroTor[$j][0]."=");
		$FORMATEDTEMP[$x] = $nonFroTor[$j][2];
#		print $FORMATEDTEMP[$x];
	}
        for($j=0;$j<=$fTorLeng;$j++){
                if ($frozenTor[$j][1]==1){
                        $x=&extract("$frozenTor[$j][0].$search_motif","$i");
                        $frozenTor[$j][2]=&formatTor($frozenTor[$j][0]."=",$x);
                }
                else {$frozenTor[$j][2]=&findExtract("$inPrefix$i$inExt",$frozenTor[$j][0]."=");}
#                print "Torsion: ".$frozenTor[$j][0].", File indicator: ".$frozenTor[$j][1].", Line replacement value: ".$frozenTor[$j][2];
                $x = &findIndex($frozenTor[$j][0]."=");
#		print "Before: X is: ".$x." and array line value is: ".$FORMATEDTEMP[$x];
                splice(@FORMATEDTEMP,$x,1);
		push(@FORMATEDTEMP,"$frozenTor[$j][2]");
#                print "After: X is: ".$x." and array line value is: ".$FORMATEDTEMP[$x];
	}
	push(@FORMATEDTEMP,"\n");
	push(@FORMATEDTEMP,"--link1--\n");
	push(@FORMATEDTEMP,"%chk=$outPrefix$i.chk\n");
	push(@FORMATEDTEMP,"#".$gauRespTheory." geom=allcheck guess=read pop=chelpg iop(6/33=2) scf=tight\n");
	push(@FORMATEDTEMP,"\n");
	push(@FORMATEDTEMP,$inPrefix.$i." charge calculation\n");
	push(@FORMATEDTEMP,"\n");
	$x = &findIndex("#");
	push(@FORMATEDTEMP,$FORMATEDTEMP[$x+4]."\n");
        unshift(@FORMATEDTEMP, "%chk=$outPrefix$i.chk\n");
	open(OUTPUT,">$outPrefix$i.$gauss_ext");
	foreach(@FORMATEDTEMP){ print OUTPUT $_;}
	close(OUTPUT);
} # End of output file index read through
print "Done file processing\n";
############################ End of file processing #####################################################

