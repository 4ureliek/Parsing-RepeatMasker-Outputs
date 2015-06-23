#!/usr/bin/perl -w

#######################################################
# Author  :  Aurelie Kapusta
#			 - adapted from another script (Qi Wang, former student in C. Feschotte Lab)
# email   :  4urelie.k@gmail.com
# version :  1.4 (see updates below)
#######################################################
# Date    : v1.0, Dec 2012
# Updated : v1.1, Nov 2013
#		- made loop on columns dynamic to take in account different outputs (with correction of fragmentation or not for ex.)
#		- updated filtering out headers of RMoutput
# Updated : v1.2, Mar 2014
#		- changed structure of the code; subs, loops, $nesting hash => start reducing code length and variable amounts + make the code more clear
#		- added option TEclass to correct class/fam
#		- This was needed to implement counting of nesting cases based on classes
# Updated : v1.3, Apr 2014
#		- Bug fix for TEc (was not defined if there was not a / in the classfam column or RMoutput)
#		- Bug fix 
######################################################
use strict;
use warnings;
use Getopt::Long;

my $v = "1.3";
my $usage = "Usage, v$v:
     perl <scriptname.pl> -RMout <MaskedGenome.out> [-TEage <TEage>] [-TEs <TEclass>]

	(note that order of these options doesn't matter)
	
	MANDATORY ARGUMENTS:
	 -RMout <MaskedGenome.out>  => <MaskedGenome.out> = file to analyse
	
	[OPTIONAL ARGUMENTS]:
	 -TEage  = tabulated file with TEs information including age
	           with columns as follow: Rname \\t Rclass \\t Rfamily \\t Rclass/Rfam \\t avg_%div \\t lineage \\t etc
	           if a TE is not found, age category will be \"undet\", for undetermined
			   If this option is chosen, there will be an extra output based on age.
			   If this option is chosen, -TEs is not necessary and even if provided, won't be used (redundant info)
	 -TEs    = tabulated file with TEs information, can be only:
	           with columns as follow: Rname \\t Rclass \\t Rfamily \\t [etc; other columns won't be used even if there are some]
	           this is to allow cleaner output for the classes (useful if some repeat names are not formatted like: Rname#Rclass/Rfamily)
			   To build this file, all repeat names + current class / fam can be obtained easily using output of parseRM_vX.X_ak.pl script (and then checked and corrected)

	WHAT IT DOES: 
		This script reads a Repeat Masker output (.out) and find nested groups, see below for more details on method and outputs
		Note that the TE that suffered the insertion is the nesting TE. The TE that inserted is the nested TE.
		It provides 3 outputs:
			- only nested/nesting blocks
			- original lines of RM output, with annotations of nesting/nested
			- all lines, but with corrections of coordinates. TEs fragmented by nesting events are merged in one line (noted by additional column with number of frags in it)
			  note that coordinates of the nesting TE will be true, but WRONG TO CALCULATE LENGTH MASKED BY IT. For that, use additional column with real lenght
	 
	METHOD (from Qi Wang)
		For each TE, the previous and next TE in the file are compared. If the previous and next TEs have:
			- the same Repeat Name
			- the same strand
			- genomic End of the previous TE is within 50 bp of the genomic Start of the nested TE
			- genomic Start of the next TE is within 50 bp of the genomic End of the nested TE
			- repeat End (in consensus) of the previous TE is within +/- 20 bp of the repeat Start (in consensus) of the next TE
		Then the TE is determined to be nested within the previous/next TE.
	 
	SUMMARY OF THE NESTED/NESTING STRUCTURES FOUND BY THE SCRIPT (frg = fragment):
		Frg A in B:                           [BBBBBB][AAAAAA][BBBBBB]
		Frg A in B in C:              [CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC]
		Frg A in B in C in D: [DDDDDD][CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC][DDDDDD]
		
		Two indep. nested frg in C:                      [CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC]
		Two indep. nested frg in C, nested in D: [DDDDDD][CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC][DDDDDD]	
		
		Two frg in C:                      [CCCCCC][BBBBBB][AAAAAA][CCCCCC]
		Two frg in C, nested in D: [DDDDDD][CCCCCC][BBBBBB][AAAAAA][CCCCCC][DDDDDD]
		   
		Three frg in C:                      [CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC]
		Three frg in C, nested in D: [DDDDDD][CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC][DDDDDD]\n\n";

my $RMout;
my $TEclass;
my $TEage;
GetOptions ('RMout=s' => \$RMout, 'TEs=s' => \$TEclass, 'TEage=s' => \$TEage);

#check that RMout file is provided
if (! $RMout){
	die $usage;
}

#log file
my $log = "$RMout.Nest.GetNesting.log";
open LOG, ">$log" or die "could not open $log $!\n";
print LOG "\n --- Script started, v$v:\n";
print LOG "     RM output file to analyze = $RMout\n";
(($TEage)?(print LOG "     TEclass file defined = $TEclass, but since TE age ($TEage) is defined it will be used instead to correct class of the elements\n"):(print LOG "     TEclass file (to correct class of the elements) defined = $TEclass\n\n")) if ($TEclass);
print LOG "     TEage file (for additional output, split on age) defined = $TEage\n" if ($TEage);
print LOG "\n\n";


# Get TE infos if provided => will allow to count SINEs in LINEs etc
print LOG " --- Getting TE info from defined file...\n" if (($TEage) || ($TEclass));
my $TE;
$TE = get_TEs_infos($TEage) if ($TEage);
$TE = get_TEs_infos($TEclass) if ((! $TEage) && ($TEclass));
# 0       1       2         3        4       5			6
#Rname	Rclass	Rfam	Rclassfam	%div	~Age	Acient/LineageSpe


# Loop on RMout
# Structure =
# 0		1	2	3	4		5		6		7		8		9		10			11			12		13		14	[15			16			17]
# Score	Div	Del	Ins	Gname	Gstart	Gend	Gleft	strand	Rname	Rclass/fam	Repstart	Rend	Repleft	ID	[Ifoverlap	IfPrevCut	Nb-of-frg]
my @RMarray = ();
my $col;
print LOG " --- Storing RMout infos...\n";
open RMOUT, "<$RMout" or die print LOG "ERROR: could not open $RMout $!\n"; 
LINE: while(<RMOUT>) {
	chomp (my $line = $_);	
	next LINE if (($line =~ /position|[sS]core|Gname/) || ($line !~ /\w/)); #skip headers and white lines
	$line =~ s/^\s+//; #remove spaces in beginning of lines
	my @line = split('\s+',$line);
	$line[8] = "-" if ($line[8] eq "C");
	
	#In case there are some low complexity stuff, we don't want them (are not TEs)
	if (($TEclass) || ($TEage)) {
		my $TEname = lc ($line[9]);
		next LINE if ($TE->{$TEname} =~ "nonTE");
	}
	next LINE if ((lc ($line[10]) =~ /low[.]+complexity/) || ($line[10] =~ /simple[.]+repeat/) || ($line[10] =~ /satellite/) || ($line[10] =~ /sat/) || ($line[10] =~ /rna\/.*rna/));
	
	#now, for relevant lines:
	push @RMarray, \@line;
	$col = $#line-1 unless ($col);
}
close(RMOUT);


# Create output files
open(NEWOUT, ">$RMout".".Nest.Corr.out") or die "Can not create >$RMout".".Nest.Corr.out $!\n";
print NEWOUT "Gname\tGstart\tGend\tRealMaskedLen\tStrand\tRname\tRclass/fam\tRepstart\tRend\tRepleft\t%Div\t%Del\t%Ins\tIf-nest(ed/ing)\tID(blockNumber)\n\n";
open(NEST, ">$RMout".".Nest.Only.out") or die "Can not create >$RMout".".Nest.Only.out $!\n";
print NEST "Gname\tGstart\tGend\tRealMaskedLen\tStrand\tRname\tRclass/fam\tRepstart\tRend\tRepleft\t%Div\t%Del\t%Ins\tIf-nest(ed/ing)\tID(blockNumber)\n\n";
open(ORIG, ">$RMout".".Nest.OriginalLines.out") or die "Can not create >$RMout".".OriginalLines.out $!\n";
print ORIG "\tScore\t%Div\t%Del\t%Ins\tGname\tGstart\tGend\tGleft\tStrand\tRname\tRclass/fam\tRepstart\tRend\tRepleft\tID\tIfoverlap\tIfPrevCut\tNb-of-frg\tIf-nest(ing/ed)\tID(blockNumber)\n\n";


###############################################################################
# Loop through the array and look for nestings
###############################################################################
print LOG " --- parsing RM output to get nesting events...\n";
my $n = 1;
my %countnest = ();
my %nesting=(
      '101' => 0,                   #[BBBBBB][AAAAAA][BBBBBB]
	 '21012' => 0,          #[CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC]
	'3210123' => 0, #[DDDDDD][CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC][DDDDDD]

	 '2X2X2' => 0,          #[CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC]
	'32X2X23' => 0, #[DDDDDD][CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC][DDDDDD]
	
	 '2XXX2' => 0,          #[CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC]
	'32XXX23' => 0, #[DDDDDD][CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC][DDDDDD]
	
	 'A01A' => 0,          #[CCCCCC][BBBBBB][AAAAAA][CCCCCC]
	'CA01AC' => 0, #[DDDDDD][CCCCCC][BBBBBB][AAAAAA][CCCCCC][DDDDDD]
);

for (my $i = 1; $i < $#RMarray; $i++){
	#----------------------------------------------------------
	# store original lines of the window I am going to look at
	# +Correct TE class/fam
	my ($TEcf,$TEc,$line_ref);
	for (my $c = -3; $c <= 3; $c++) {
		($line_ref->{$c},$TEc,$TEcf) = get_line_and_TEc(\@RMarray,$c,$col,$i,$TEc,$TEcf,$TE,$TEclass,$TEage);
	}
	#----------------------------------------------------------

	#values to write lines
	for (my $c = 1; $c <= 3; $c++) {
		if (($i >= $c) && ($i < $#RMarray-$c) && ($RMarray[$i-$c]->[4] eq $RMarray[$i+$c]->[4])) {
			$nesting{'len'}{$c} = ($RMarray[$i-$c]->[6]-$RMarray[$i-$c]->[5]+1) + ($RMarray[$i+$c]->[6]-$RMarray[$i+$c]->[5]+1);
			$nesting{'div'}{$c} = (($RMarray[$i-$c]->[$c] * ($RMarray[$i-$c]->[6]-$RMarray[$i-$c]->[5]+1)) +($RMarray[$i+$c]->[$c] * ($RMarray[$i+$c]->[6]-$RMarray[$i+$c]->[5]+1)))/$nesting{'len'}{$c};
			$nesting{'del'}{$c} = (($RMarray[$i-$c]->[2] * ($RMarray[$i-$c]->[6]-$RMarray[$i-$c]->[5]+1)) +($RMarray[$i+$c]->[2] * ($RMarray[$i+$c]->[6]-$RMarray[$i+$c]->[5]+1)))/$nesting{'len'}{$c};
			$nesting{'ins'}{$c} = (($RMarray[$i-$c]->[3] * ($RMarray[$i-$c]->[6]-$RMarray[$i-$c]->[5]+1)) +($RMarray[$i+$c]->[3] * ($RMarray[$i+$c]->[6]-$RMarray[$i+$c]->[5]+1)))/$nesting{'len'}{$c};
		}
	}
	my $len = $RMarray[$i]->[6] - $RMarray[$i]->[5] + 1;
	
	##########################################################################################################################################
	# SUMMARY OF THE ACHITECTURE OF THE FOLLOWING PART OF THE SCRIPT:
	# 1) check +2 surrounding
	#	IF [CCCCCC][XXXXXX][XXXXXX][XXXXXX][CCCCCC]
	#			IF [CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC], ie check +1 surrounding
	#				check +3 surrounding
	#					IF [DDDDDD][CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC][DDDDDD] => print +3 corresponding lines
	#					else => print +2 corresponding lines
	#			IF [CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC], ie wouldn't be in previous case
	#				check +3 surrounding
	#					IF [DDDDDD][CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC][DDDDDD] => print +3 corresponding lines
	#					else => print +2 corresponding lines
	#			IF [CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC]
	#				check +3 surrounding
	#					IF [DDDDDD][CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC][DDDDDD] => print +3 corresponding lines
	#					else => print +2 corresponding lines
	
	# 2) ELSE check +1 surrounding 
	# 	IF same TE, ie [BBBBBB][AAAAAA][BBBBBB] 
	#								 -1		 0		+1											-1		 0	    +1
	#		BUT, need to check IF [BBBBBB][AAAAAA][BBBBBB][XXXXXX][BBBBBB] OR [BBBBBB][XXXXXX][BBBBBB][AAAAAA][BBBBBB] => NO PRINT, because already detected
	#		ELSE => print +1 corresponding lines
	
	# 3) check particular case of
	#	IF [CCCCCC][BBBBBB][AAAAAA][CCCCCC] 
	# 		check surrounding
	#			IF [DDDDDD][CCCCCC][BBBBBB][AAAAAA][CCCCCC][DDDDDD] 
	
	# ELSE => print current line.
	##########################################################################################################################################
	

	###########################################################
	# 1) check +2 surrounding 
	###########################################################
	#		-2		-1		0		+1		+2
	# IF [CCCCCC][XXXXXX][XXXXXX][XXXXXX][CCCCCC]
	if (($i > 1) && ($i < $#RMarray-1) # only if this is not first or last round of loop
	&& ($RMarray[$i-2]->[4] eq $RMarray[$i+2]->[4]) && ($RMarray[$i-2]->[9] eq $RMarray[$i+2]->[9]) && ($RMarray[$i-2]->[8] eq $RMarray[$i+2]->[8]) 
	&& (($RMarray[$i-1]->[5] - $RMarray[$i-2]->[6]) <= 50) && (($RMarray[$i+2]->[5] - $RMarray[$i+1]->[6]) <= 50)) { #[CCCCCC] could be nesting TE
		
		#check if same TE cut in 2 (ie stuff inside are other TEs and could be other nesting event
		if ((($RMarray[$i-2]->[8] eq '+') && (($RMarray[$i+2]->[11] - $RMarray[$i-2]->[12]) <= 20) && (($RMarray[$i+2]->[11] - $RMarray[$i-2]->[12]) >= -20)) 
		|| (($RMarray[$i-2]->[8] eq '-') && (($RMarray[$i-2]->[13] - $RMarray[$i+2]->[12]) <= 20) && (($RMarray[$i-2]->[13] - $RMarray[$i+2]->[12]) >= -20))) {

			#		-1		0		+1
			# IF [BBBBBB][AAAAAA][BBBBBB]
			if (($RMarray[$i-1]->[4] eq $RMarray[$i+1]->[4]) && ($RMarray[$i-1]->[9] eq $RMarray[$i+1]->[9]) && ($RMarray[$i-1]->[8] eq $RMarray[$i+1]->[8])
			&& ((($RMarray[$i]->[5] - $RMarray[$i-1]->[6]) <= 50) && (($RMarray[$i+1]->[5] - $RMarray[$i]->[6]) <= 50) 
			&& ((($RMarray[$i-1]->[8] eq '+') && (($RMarray[$i+1]->[11] - $RMarray[$i-1]->[12]) <= 20) && (($RMarray[$i+1]->[11] - $RMarray[$i-1]->[12]) >= -20)) 
			|| (($RMarray[$i-1]->[8] eq '-') && (($RMarray[$i-1]->[13] - $RMarray[$i+1]->[12]) <= 20) && (($RMarray[$i-1]->[13] - $RMarray[$i+1]->[12]) >= -20))))) {

				#				  -3	  -2	  -1	  0		  +1	  +2	 +3
				# CHECK +3; IF [DDDDDD][CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC][DDDDDD]
				if (($i > 2) && ($i < $#RMarray-2) 
				&& ($RMarray[$i-3]->[4] eq $RMarray[$i+3]->[4]) && ($RMarray[$i-3]->[9] eq $RMarray[$i+3]->[9]) && ($RMarray[$i-3]->[8] eq $RMarray[$i+3]->[8])
				&& ((($RMarray[$i-2]->[5] - $RMarray[$i-3]->[6]) <= 50) && (($RMarray[$i+3]->[5] - $RMarray[$i+2]->[6]) <= 50) 
				&& ((($RMarray[$i-3]->[8] eq '+') && (($RMarray[$i+3]->[11] - $RMarray[$i-3]->[12]) <= 20) && (($RMarray[$i+3]->[11] - $RMarray[$i-3]->[12]) >= -20)) 
				|| (($RMarray[$i-3]->[8] eq '-') && (($RMarray[$i-3]->[13] - $RMarray[$i+3]->[12]) <= 20) && (($RMarray[$i-3]->[13] - $RMarray[$i+3]->[12]) >= -20))))) {
					
					$nesting{'3210123'}++;
					$countnest{$TEc->{3}}{$TEc->{2}}++;
					
					# 1) print original lines
					print ORIG "$line_ref->{-3}\tnesting\t$n\n$line_ref->{-2}\tnesting/nested\t$n\n$line_ref->{-1}\tnesting/nested\t$n\n$line_ref->{0}\tnested\t$n\n$line_ref->{1}\tnesting/nested\t$n\n$line_ref->{2}\tnesting/nested\t$n\n$line_ref->{3}\tnesting\t$n\n";
					
					# 2) print new files
					# print nesting TEs
					my $Rend;
					#print [DDDDDD]
					if ($RMarray[$i-3]->[8] eq '+') { 
						$Rend = $RMarray[$i+3]->[12];
					} else {
						$Rend = $RMarray[$i-3]->[12];
					}				
					print NEWOUT "$RMarray[$i-3]->[4]\t$RMarray[$i-3]->[5]\t$RMarray[$i+3]->[6]\t$nesting{'len'}{3}\t$RMarray[$i-3]->[8]\t$RMarray[$i-3]->[9]\t$TEcf->{-3}\t$RMarray[$i-3]->[11]\t$Rend\t$RMarray[$i+3]->[13]\t$nesting{'div'}{3}\t$nesting{'del'}{3}\t$nesting{'ins'}{3}\tnesting\t$n\n";
					print NEST "$RMarray[$i-3]->[4]\t$RMarray[$i-3]->[5]\t$RMarray[$i+3]->[6]\t$nesting{'len'}{3}\t$RMarray[$i-3]->[8]\t$RMarray[$i-3]->[9]\t$TEcf->{-3}\t$RMarray[$i-3]->[11]\t$Rend\t$RMarray[$i+3]->[13]\t$nesting{'div'}{3}\t$nesting{'del'}{3}\t$nesting{'ins'}{3}\tnesting\t$n\n";
					#print [CCCCCC]
					if ($RMarray[$i-2]->[8] eq '+') { 
						$Rend = $RMarray[$i+2]->[12];
					} else {
						$Rend = $RMarray[$i-2]->[12];
					}	
					print NEWOUT "$RMarray[$i-2]->[4]\t$RMarray[$i-2]->[5]\t$RMarray[$i+2]->[6]\t$nesting{'len'}{2}\t$RMarray[$i-2]->[8]\t$RMarray[$i-2]->[9]\t$TEcf->{-2}\t$RMarray[$i-2]->[11]\t$Rend\t$RMarray[$i+2]->[13]\t$nesting{'div'}{2}\t$nesting{'del'}{2}\t$nesting{'ins'}{2}\tnesting/nested\t$n\n";
					print NEST "$RMarray[$i-2]->[4]\t$RMarray[$i-2]->[5]\t$RMarray[$i+2]->[6]\t$nesting{'len'}{2}\t$RMarray[$i-2]->[8]\t$RMarray[$i-2]->[9]\t$TEcf->{-2}\t$RMarray[$i-2]->[11]\t$Rend\t$RMarray[$i+2]->[13]\t$nesting{'div'}{2}\t$nesting{'del'}{2}\t$nesting{'ins'}{2}\tnesting/nested\t$n\n";
					#print [BBBBBB]
					if ($RMarray[$i-1]->[8] eq '+') { 
						$Rend = $RMarray[$i+1]->[12];
					} else {
						$Rend = $RMarray[$i-1]->[12];
					}	
					print NEWOUT "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i+1]->[6]\t$nesting{'len'}{1}\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$Rend\t$RMarray[$i+1]->[13]\t$nesting{'div'}{1}\t$nesting{'del'}{1}\t$nesting{'ins'}{1}\tnesting/nested\t$n\n";
					print NEST "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i+1]->[6]\t$nesting{'len'}{1}\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$Rend\t$RMarray[$i+1]->[13]\t$nesting{'div'}{1}\t$nesting{'del'}{1}\t$nesting{'ins'}{1}\tnesting/nested\t$n\n";
					# print nested TE 
					#print [AAAAAA]
					print NEWOUT "$RMarray[$i]->[4]\t$RMarray[$i]->[5]\t$RMarray[$i]->[6]\t$len\t$RMarray[$i]->[8]\t$RMarray[$i]->[9]\t$TEcf->{0}\t$RMarray[$i]->[11]\t$RMarray[$i]->[12]\t$RMarray[$i]->[13]\t$RMarray[$i]->[1]\t$RMarray[$i]->[2]\t$RMarray[$i]->[3]\tnested\t$n\n";
					print NEST "$RMarray[$i]->[4]\t$RMarray[$i]->[5]\t$RMarray[$i]->[6]\t$len\t$RMarray[$i]->[8]\t$RMarray[$i]->[9]\t$TEcf->{0}\t$RMarray[$i]->[11]\t$RMarray[$i]->[12]\t$RMarray[$i]->[13]\t$RMarray[$i]->[1]\t$RMarray[$i]->[2]\t$RMarray[$i]->[3]\tnested\t$n\n";	
					$n++;
				} 
				
				else {  #ie no +3 => [CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC] 
					$nesting{'21012'}++;
					$countnest{$TEc->{2}}{$TEc->{1}}++;
					
					# 1) print original lines
					print ORIG "$line_ref->{-2}\tnesting\t$n\n$line_ref->{-1}\tnesting/nested\t$n\n$line_ref->{0}\tnested\t$n\n$line_ref->{1}\tnesting/nested\t$n\n$line_ref->{2}\tnesting\t$n\n";
				
					# 2) print new files
					# print nesting TEs
					my $Rend;
					#print [CCCCCC]
					if ($RMarray[$i-2]->[8] eq '+') { 
						$Rend = $RMarray[$i+2]->[12];
					} else {
						$Rend = $RMarray[$i-2]->[12];
					}
					print NEWOUT "$RMarray[$i-2]->[4]\t$RMarray[$i-2]->[5]\t$RMarray[$i+2]->[6]\t$nesting{'len'}{2}\t$RMarray[$i-2]->[8]\t$RMarray[$i-2]->[9]\t$TEcf->{-2}\t$RMarray[$i-2]->[11]\t$Rend\t$RMarray[$i+2]->[13]\t$nesting{'div'}{2}\t$nesting{'del'}{2}\t$nesting{'ins'}{2}\tnesting/nested\t$n\n";
					print NEST "$RMarray[$i-2]->[4]\t$RMarray[$i-2]->[5]\t$RMarray[$i+2]->[6]\t$nesting{'len'}{2}\t$RMarray[$i-2]->[8]\t$RMarray[$i-2]->[9]\t$TEcf->{-2}\t$RMarray[$i-2]->[11]\t$Rend\t$RMarray[$i+2]->[13]\t$nesting{'div'}{2}\t$nesting{'del'}{2}\t$nesting{'ins'}{2}\tnesting/nested\t$n\n";
					#print [BBBBBB]
					if ($RMarray[$i-1]->[8] eq '+') {
						$Rend = $RMarray[$i+1]->[12];
					} else {
						$Rend = $RMarray[$i-1]->[12];
					}
					print NEWOUT "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i+1]->[6]\t$nesting{'len'}{1}\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$Rend\t$RMarray[$i+1]->[13]\t$nesting{'div'}{1}\t$nesting{'del'}{1}\t$nesting{'ins'}{1}\tnesting/nested\t$n\n";
					print NEST "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i+1]->[6]\t$nesting{'len'}{1}\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$Rend\t$RMarray[$i+1]->[13]\t$nesting{'div'}{1}\t$nesting{'del'}{1}\t$nesting{'ins'}{1}\tnesting/nested\t$n\n";
					#print [AAAAAA]
					# print nested TE
					print NEWOUT "$RMarray[$i]->[4]\t$RMarray[$i]->[5]\t$RMarray[$i]->[6]\t$len\t$RMarray[$i]->[8]\t$RMarray[$i]->[9]\t$TEcf->{0}\t$RMarray[$i]->[11]\t$RMarray[$i]->[12]\t$RMarray[$i]->[13]\t$RMarray[$i]->[1]\t$RMarray[$i]->[2]\t$RMarray[$i]->[3]\tnested\t$n\n";
					print NEST "$RMarray[$i]->[4]\t$RMarray[$i]->[5]\t$RMarray[$i]->[6]\t$len\t$RMarray[$i]->[8]\t$RMarray[$i]->[9]\t$TEcf->{0}\t$RMarray[$i]->[11]\t$RMarray[$i]->[12]\t$RMarray[$i]->[13]\t$RMarray[$i]->[1]\t$RMarray[$i]->[2]\t$RMarray[$i]->[3]\tnested\t$n\n";
					$n++;
				}
			}	
		} elsif 
		# Now, +2 is potentially nesting TE, but not +1; is it 2 indep nestings in +2? check if IF [CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC]
		(
		($RMarray[$i-2]->[4] eq $RMarray[$i]->[4]) && ($RMarray[$i-2]->[9] eq $RMarray[$i]->[9]) && ($RMarray[$i-2]->[8] eq $RMarray[$i]->[8])
		&& (($RMarray[$i]->[5] - $RMarray[$i-1]->[6]) <= 50) && (($RMarray[$i+1]->[5] - $RMarray[$i]->[6]) <= 50) 
		&& ((($RMarray[$i-2]->[8] eq '+') && (($RMarray[$i]->[11] - $RMarray[$i-2]->[12]) <= 20) && (($RMarray[$i]->[11] - $RMarray[$i-2]->[12]) >= -20) && (($RMarray[$i+2]->[11] - $RMarray[$i]->[12]) <= 20) && (($RMarray[$i+2]->[11] - $RMarray[$i]->[12]) <= 20))
		|| (($RMarray[$i-2]->[8] eq '-') && (($RMarray[$i-2]->[13] - $RMarray[$i]->[12]) <= 20) && (($RMarray[$i-2]->[13] - $RMarray[$i]->[12]) >= -20) && (($RMarray[$i]->[13] - $RMarray[$i+2]->[12]) <= 20) && (($RMarray[$i]->[13] - $RMarray[$i+2]->[12]) >= -20)))
		) {
			# [CCCCCC] will require only one line, and recalculation of the %div del ins and len masked 	
			my $nesting2X2_len = ($RMarray[$i-2]->[6]-$RMarray[$i-2]->[5]+1) + ($RMarray[$i+2]->[6]-$RMarray[$i+2]->[5]+1) + ($RMarray[$i]->[6]-$RMarray[$i]->[5]+1);
			my $nesting2X2_div = (($RMarray[$i-1]->[1] * ($RMarray[$i-1]->[6]-$RMarray[$i-1]->[5]+1)) + ($RMarray[$i+1]->[1] * ($RMarray[$i+1]->[6]-$RMarray[$i+1]->[5]+1)) + ($RMarray[$i]->[1] * ($RMarray[$i]->[6]-$RMarray[$i]->[5]+1)))/$nesting2X2_len;
			my $nesting2X2_del = (($RMarray[$i-1]->[2] * ($RMarray[$i-1]->[6]-$RMarray[$i-1]->[5]+1)) + ($RMarray[$i+1]->[2] * ($RMarray[$i+1]->[6]-$RMarray[$i+1]->[5]+1)) + ($RMarray[$i]->[2] * ($RMarray[$i]->[6]-$RMarray[$i]->[5]+1)))/$nesting2X2_len;
			my $nesting2X2_ins = (($RMarray[$i-1]->[3] * ($RMarray[$i-1]->[6]-$RMarray[$i-1]->[5]+1)) + ($RMarray[$i+1]->[3] * ($RMarray[$i+1]->[6]-$RMarray[$i+1]->[5]+1)) + ($RMarray[$i]->[3] * ($RMarray[$i]->[6]-$RMarray[$i]->[5]+1)))/$nesting2X2_len;
			
			# CHECK +3; IF [DDDDDD][CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC][DDDDDD]
			if (($i > 2) && ($i < $#RMarray-2) 
			&& ($RMarray[$i-3]->[4] eq $RMarray[$i+3]->[4]) && ($RMarray[$i-3]->[9] eq $RMarray[$i+3]->[9]) && ($RMarray[$i-3]->[8] eq $RMarray[$i+3]->[8])
			&& ((($RMarray[$i-2]->[5] - $RMarray[$i-3]->[6]) <= 50) && (($RMarray[$i+3]->[5] - $RMarray[$i+2]->[6]) <= 50) 
			&& ((($RMarray[$i-3]->[8] eq '+') && (($RMarray[$i+3]->[11] - $RMarray[$i-3]->[12]) <= 20) && (($RMarray[$i+3]->[11] - $RMarray[$i-3]->[12]) >= -20)) 
			|| (($RMarray[$i-3]->[8] eq '-') && (($RMarray[$i-3]->[13] - $RMarray[$i+3]->[12]) <= 20) && (($RMarray[$i-3]->[13] - $RMarray[$i+3]->[12]) >= -20))))) { 
				
				$nesting{'32X2X23'}++;
				$countnest{$TEc->{3}}{$TEc->{2}}++;
				
				# 1) print original lines
				print ORIG "$line_ref->{-3}\tnesting\t$n\n$line_ref->{-2}\tnesting/nested\t$n\n$line_ref->{-1}\tnested\t$n\n$line_ref->{0}\tnesting/nested\t$n\n$line_ref->{1}\tnested\t$n\n$line_ref->{2}\tnesting/nested\t$n\n$line_ref->{3}\tnesting\t$n\n";
				
				# 2) print new files
				# print nesting TEs
				#print [DDDDDD]
				my $Rend;
				if ($RMarray[$i-3]->[8] eq '+') {
					$Rend = $RMarray[$i+3]->[12];
				} else {
					$Rend = $RMarray[$i-3]->[12];
				}
				print NEWOUT "$RMarray[$i-3]->[4]\t$RMarray[$i-3]->[5]\t$RMarray[$i+3]->[6]\t$nesting{'len'}{3}\t$RMarray[$i-3]->[8]\t$RMarray[$i-3]->[9]\t$TEcf->{-3}\t$RMarray[$i-3]->[11]\t$Rend\t$RMarray[$i+3]->[13]\t$nesting{'div'}{3}\t$nesting{'del'}{3}\t$nesting{'ins'}{3}\tnesting\t$n\n";
				print NEST "$RMarray[$i-3]->[4]\t$RMarray[$i-3]->[5]\t$RMarray[$i+3]->[6]\t$nesting{'len'}{3}\t$RMarray[$i-3]->[8]\t$RMarray[$i-3]->[9]\t$TEcf->{-3}\t$RMarray[$i-3]->[11]\t$Rend\t$RMarray[$i+3]->[13]\t$nesting{'div'}{3}\t$nesting{'del'}{3}\t$nesting{'ins'}{3}\tnesting\t$n\n";
				if ($RMarray[$i-2]->[8] eq '+') {
					$Rend = $RMarray[$i+2]->[12];
				} else {
					$Rend = $RMarray[$i-2]->[12];
				}
				print NEWOUT "$RMarray[$i-2]->[4]\t$RMarray[$i-2]->[5]\t$RMarray[$i+2]->[6]\t$nesting2X2_len\t$RMarray[$i-2]->[8]\t$RMarray[$i-2]->[9]\t$TEcf->{-2}\t$RMarray[$i-2]->[11]\t$Rend\t$RMarray[$i+2]->[13]\t$nesting2X2_div\t$nesting2X2_del\t$nesting2X2_ins\tnesting/nested\t$n\n";
				print NEST "$RMarray[$i-2]->[4]\t$RMarray[$i-2]->[5]\t$RMarray[$i+2]->[6]\t$nesting2X2_len\t$RMarray[$i-2]->[8]\t$RMarray[$i-2]->[9]\t$TEcf->{-2}\t$RMarray[$i-2]->[11]\t$Rend\t$RMarray[$i+2]->[13]\t$nesting2X2_div\t$nesting2X2_del\t$nesting2X2_ins\tnesting/nested\t$n\n";
				# print nested TEs
				print NEWOUT "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i-1]->[6]\t$len\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$RMarray[$i-1]->[12]\t$RMarray[$i-1]->[13]\t$RMarray[$i-1]->[1]\t$RMarray[$i-1]->[2]\t$RMarray[$i-1]->[3]\tnested\t$n\n";
				print NEWOUT "$RMarray[$i+1]->[4]\t$RMarray[$i+1]->[5]\t$RMarray[$i+1]->[6]\t$len\t$RMarray[$i+1]->[8]\t$RMarray[$i+1]->[9]\t$TEcf->{1}\t$RMarray[$i+1]->[11]\t$RMarray[$i+1]->[12]\t$RMarray[$i+1]->[13]\t$RMarray[$i+1]->[1]\t$RMarray[$i+1]->[2]\t$RMarray[$i+1]->[3]\tnested\t$n\n";
				print NEST "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i-1]->[6]\t$len\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$RMarray[$i-1]->[12]\t$RMarray[$i-1]->[13]\t$RMarray[$i-1]->[1]\t$RMarray[$i-1]->[2]\t$RMarray[$i-1]->[3]\tnested\t$n\n";
				print NEST "$RMarray[$i+1]->[4]\t$RMarray[$i+1]->[5]\t$RMarray[$i+1]->[6]\t$len\t$RMarray[$i+1]->[8]\t$RMarray[$i+1]->[9]\t$TEcf->{1}\t$RMarray[$i+1]->[11]\t$RMarray[$i+1]->[12]\t$RMarray[$i+1]->[13]\t$RMarray[$i+1]->[1]\t$RMarray[$i+1]->[2]\t$RMarray[$i+1]->[3]\tnested\t$n\n";	
				$n++;
			} 
			
			else { # ie no +3 => [CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC]
				$nesting{'2X2X2'}++;
				$countnest{$TEc->{2}}{$TEc->{-1}}++;
				$countnest{$TEc->{2}}{$TEc->{1}}++; #Note that can be an over estimation, but cases like that are not the most frequent
				
				# 1) print original lines
				print ORIG "$line_ref->{-2}\tnesting\t$n\n$line_ref->{-1}\tnested\t$n\n$line_ref->{0}\tnesting\t$n\n$line_ref->{1}\tnested\t$n\n$line_ref->{2}\tnesting\t$n\n";
				
				# 2) print new files
				# print nesting TEs	
				my $Rend;
				if ($RMarray[$i-2]->[8] eq '+') {
					$Rend = $RMarray[$i+2]->[12];
				} else {
					$Rend = $RMarray[$i-2]->[12];
				}
				print NEWOUT "$RMarray[$i-2]->[4]\t$RMarray[$i-2]->[5]\t$RMarray[$i+2]->[6]\t$nesting2X2_len\t$RMarray[$i-2]->[8]\t$RMarray[$i-2]->[9]\t$TEcf->{-2}\t$RMarray[$i-2]->[11]\t$Rend\t$RMarray[$i+2]->[13]\t$nesting2X2_div\t$nesting2X2_del\t$nesting2X2_ins\tnesting/nested\t$n\n";
				print NEST "$RMarray[$i-2]->[4]\t$RMarray[$i-2]->[5]\t$RMarray[$i+2]->[6]\t$nesting2X2_len\t$RMarray[$i-2]->[8]\t$RMarray[$i-2]->[9]\t$TEcf->{-2}\t$RMarray[$i-2]->[11]\t$Rend\t$RMarray[$i+2]->[13]\t$nesting2X2_div\t$nesting2X2_del\t$nesting2X2_ins\tnesting/nested\t$n\n";
				# print nested TEs
				print NEWOUT "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i-1]->[6]\t$len\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$RMarray[$i-1]->[12]\t$RMarray[$i-1]->[13]\t$RMarray[$i-1]->[1]\t$RMarray[$i-1]->[2]\t$RMarray[$i-1]->[3]\tnested\t$n\n";
				print NEWOUT "$RMarray[$i+1]->[4]\t$RMarray[$i+1]->[5]\t$RMarray[$i+1]->[6]\t$len\t$RMarray[$i+1]->[8]\t$RMarray[$i+1]->[9]\t$TEcf->{1}\t$RMarray[$i+1]->[11]\t$RMarray[$i+1]->[12]\t$RMarray[$i+1]->[13]\t$RMarray[$i+1]->[1]\t$RMarray[$i+1]->[2]\t$RMarray[$i+1]->[3]\tnested\t$n\n";
				print NEST "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i-1]->[6]\t$len\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$RMarray[$i-1]->[12]\t$RMarray[$i-1]->[13]\t$RMarray[$i-1]->[1]\t$RMarray[$i-1]->[2]\t$RMarray[$i-1]->[3]\tnested\t$n\n";
				print NEST "$RMarray[$i+1]->[4]\t$RMarray[$i+1]->[5]\t$RMarray[$i+1]->[6]\t$len\t$RMarray[$i+1]->[8]\t$RMarray[$i+1]->[9]\t$TEcf->{1}\t$RMarray[$i+1]->[11]\t$RMarray[$i+1]->[12]\t$RMarray[$i+1]->[13]\t$RMarray[$i+1]->[1]\t$RMarray[$i+1]->[2]\t$RMarray[$i+1]->[3]\tnested\t$n\n";	
				$n++;
			}
		} elsif 
		#so it is not double nesting or 2 indep nesting - could it be 3 other pieces of 3 indep TEs? => [CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC]
		((($RMarray[$i-1]->[5] - $RMarray[$i-2]->[6]) <= 50) && (($RMarray[$i]->[5] - $RMarray[$i-1]->[6]) <= 50)  && (($RMarray[$i+1]->[5] - $RMarray[$i]->[6]) <= 50) && (($RMarray[$i+2]->[5] - $RMarray[$i+1]->[6]) <= 50) ) {

			# CHECK +3; IF [DDDDDD][CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC][DDDDDD]
			if (($i > 2) && ($i < $#RMarray-2) 
			&& ($RMarray[$i-3]->[4] eq $RMarray[$i+3]->[4]) && ($RMarray[$i-3]->[9] eq $RMarray[$i+3]->[9]) && ($RMarray[$i-3]->[8] eq $RMarray[$i+3]->[8])
			&& ((($RMarray[$i-2]->[5] - $RMarray[$i-3]->[6]) <= 50) && (($RMarray[$i+3]->[5] - $RMarray[$i+2]->[6]) <= 50) 
			&& ((($RMarray[$i-3]->[8] eq '+') && (($RMarray[$i+3]->[11] - $RMarray[$i-3]->[12]) <= 20) && (($RMarray[$i+3]->[11] - $RMarray[$i-3]->[12]) >= -20)) 
			|| (($RMarray[$i-3]->[8] eq '-') && (($RMarray[$i-3]->[13] - $RMarray[$i+3]->[12]) <= 20) && (($RMarray[$i-3]->[13] - $RMarray[$i+3]->[12]) >= -20))))) { 
				
				$nesting{'32XXX23'}++;
				$countnest{$TEc->{3}}{$TEc->{2}}++;
				
				# 1) print original lines
				print ORIG "$line_ref->{-3}\tnesting\t$n\n$line_ref->{-2}\tnesting/nested\t$n\n$line_ref->{-1}\tnested\t$n\n$line_ref->{0}\tnested\t$n\n$line_ref->{1}\tnested\t$n\n$line_ref->{2}\tnesting/nested\t$n\n$line_ref->{3}\tnesting\t$n\n";
				
				# 2) print new files
				# print nesting TEs
				my $Rend;
				# print [DDDDDD]
				if ($RMarray[$i-3]->[8] eq '+') {
					$Rend = $RMarray[$i+3]->[12];
				} else {
					$Rend = $RMarray[$i-3]->[12];
				}
				print NEWOUT "$RMarray[$i-3]->[4]\t$RMarray[$i-3]->[5]\t$RMarray[$i+3]->[6]\t$nesting{'len'}{3}\t$RMarray[$i-3]->[8]\t$RMarray[$i-3]->[9]\t$TEcf->{-3}\t$RMarray[$i-3]->[11]\t$Rend\t$RMarray[$i+3]->[13]\t$nesting{'div'}{3}\t$nesting{'del'}{3}\t$nesting{'ins'}{3}\tnesting\t$n\n";
				print NEST "$RMarray[$i-3]->[4]\t$RMarray[$i-3]->[5]\t$RMarray[$i+3]->[6]\t$nesting{'len'}{3}\t$RMarray[$i-3]->[8]\t$RMarray[$i-3]->[9]\t$TEcf->{-3}\t$RMarray[$i-3]->[11]\t$Rend\t$RMarray[$i+3]->[13]\t$nesting{'div'}{3}\t$nesting{'del'}{3}\t$nesting{'ins'}{3}\tnesting\t$n\n";
				# print [CCCCCC]
				if ($RMarray[$i-2]->[8] eq '+') {
					$Rend = $RMarray[$i+2]->[12];
				} else {
					$Rend = $RMarray[$i-2]->[12];
				}
				print NEWOUT "$RMarray[$i-2]->[4]\t$RMarray[$i-2]->[5]\t$RMarray[$i+2]->[6]\t$nesting{'len'}{2}\t$RMarray[$i-2]->[8]\t$RMarray[$i-2]->[9]\t$TEcf->{-2}\t$RMarray[$i-2]->[11]\t$Rend\t$RMarray[$i+2]->[13]\t$nesting{'div'}{2}\t$nesting{'del'}{2}\t$nesting{'ins'}{2}\tnesting/nested\t$n\n";
				print NEST "$RMarray[$i-2]->[4]\t$RMarray[$i-2]->[5]\t$RMarray[$i+2]->[6]\t$nesting{'len'}{2}\t$RMarray[$i-2]->[8]\t$RMarray[$i-2]->[9]\t$TEcf->{-2}\t$RMarray[$i-2]->[11]\t$Rend\t$RMarray[$i+2]->[13]\t$nesting{'div'}{2}\t$nesting{'del'}{2}\t$nesting{'ins'}{2}\tnesting/nested\t$n\n";
				# print nested TEs (all 3 insides)
				print NEWOUT "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i-1]->[6]\t$len\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$RMarray[$i-1]->[12]\t$RMarray[$i-1]->[13]\t$RMarray[$i-1]->[1]\t$RMarray[$i-1]->[2]\t$RMarray[$i-1]->[3]\tnested\t$n\n";
				print NEWOUT "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i-1]->[6]\t$len\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$RMarray[$i-1]->[12]\t$RMarray[$i-1]->[13]\t$RMarray[$i-1]->[1]\t$RMarray[$i-1]->[2]\t$RMarray[$i-1]->[3]\tnested\t$n\n";
				print NEWOUT "$RMarray[$i+1]->[4]\t$RMarray[$i+1]->[5]\t$RMarray[$i+1]->[6]\t$len\t$RMarray[$i+1]->[8]\t$RMarray[$i+1]->[9]\t$TEcf->{1}\t$RMarray[$i+1]->[11]\t$RMarray[$i+1]->[12]\t$RMarray[$i+1]->[13]\t$RMarray[$i+1]->[1]\t$RMarray[$i+1]->[2]\t$RMarray[$i+1]->[3]\tnested\t$n\n";
				print NEST "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i-1]->[6]\t$len\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$RMarray[$i-1]->[12]\t$RMarray[$i-1]->[13]\t$RMarray[$i-1]->[1]\t$RMarray[$i-1]->[2]\t$RMarray[$i-1]->[3]\tnested\t$n\n";
				print NEST "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i-1]->[6]\t$len\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$RMarray[$i-1]->[12]\t$RMarray[$i-1]->[13]\t$RMarray[$i-1]->[1]\t$RMarray[$i-1]->[2]\t$RMarray[$i-1]->[3]\tnested\t$n\n";	
				print NEST "$RMarray[$i+1]->[4]\t$RMarray[$i+1]->[5]\t$RMarray[$i+1]->[6]\t$len\t$RMarray[$i+1]->[8]\t$RMarray[$i+1]->[9]\t$TEcf->{1}\t$RMarray[$i+1]->[11]\t$RMarray[$i+1]->[12]\t$RMarray[$i+1]->[13]\t$RMarray[$i+1]->[1]\t$RMarray[$i+1]->[2]\t$RMarray[$i+1]->[3]\tnested\t$n\n";	
				$n++;
			} else {
				# Still +2 surrounding
				#else => [CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC]
				$nesting{'2XXX2'}++;
				$countnest{$TEc->{2}}{$TEc->{0}}++;
				$countnest{$TEc->{2}}{$TEc->{-1}}++ unless ($TEc->{1} eq $TEc->{0}); #trying to avoid counting multitude stuff that might be same element
				$countnest{$TEc->{2}}{$TEc->{1}}++ unless ($TEc->{1} eq $TEc->{0}); #trying to avoid counting multitude stuff that might be same element
				
				# 1) print original lines
				print ORIG "$line_ref->{-2}\tnesting\t$n\n$line_ref->{-1}\tnested\t$n\n$line_ref->{0}\tnested\t$n\n$line_ref->{1}\tnested\t$n\n$line_ref->{2}\tnesting\t$n\n";
				
				# 2) print new files
				# print nesting TE
				my $Rend;
				# print [CCCCCC]
				if ($RMarray[$i-2]->[8] eq '+') {
					$Rend = $RMarray[$i+2]->[12];
				} else {
					$Rend = $RMarray[$i-2]->[12];
				}
				print NEWOUT "$RMarray[$i-2]->[4]\t$RMarray[$i-2]->[5]\t$RMarray[$i+2]->[6]\t$nesting{'len'}{2}\t$RMarray[$i-2]->[8]\t$RMarray[$i-2]->[9]\t$TEcf->{-2}\t$RMarray[$i-2]->[11]\t$Rend\t$RMarray[$i+2]->[13]\t$nesting{'div'}{2}\t$nesting{'del'}{2}\t$nesting{'ins'}{2}\tnesting\t$n\n";
				print NEST "$RMarray[$i-2]->[4]\t$RMarray[$i-2]->[5]\t$RMarray[$i+2]->[6]\t$nesting{'len'}{2}\t$RMarray[$i-2]->[8]\t$RMarray[$i-2]->[9]\t$TEcf->{-2}\t$RMarray[$i-2]->[11]\t$Rend\t$RMarray[$i+2]->[13]\t$nesting{'div'}{2}\t$nesting{'del'}{2}\t$nesting{'ins'}{2}\tnesting\t$n\n";
				# print nested TEs (all 3 insides)
				print NEWOUT "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i-1]->[6]\t$len\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$RMarray[$i-1]->[12]\t$RMarray[$i-1]->[13]\t$RMarray[$i-1]->[1]\t$RMarray[$i-1]->[2]\t$RMarray[$i-1]->[3]\tnested\t$n\n";
				print NEWOUT "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i-1]->[6]\t$len\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$RMarray[$i-1]->[12]\t$RMarray[$i-1]->[13]\t$RMarray[$i-1]->[1]\t$RMarray[$i-1]->[2]\t$RMarray[$i-1]->[3]\tnested\t$n\n";
				print NEWOUT "$RMarray[$i+1]->[4]\t$RMarray[$i+1]->[5]\t$RMarray[$i+1]->[6]\t$len\t$RMarray[$i+1]->[8]\t$RMarray[$i+1]->[9]\t$TEcf->{1}\t$RMarray[$i+1]->[11]\t$RMarray[$i+1]->[12]\t$RMarray[$i+1]->[13]\t$RMarray[$i+1]->[1]\t$RMarray[$i+1]->[2]\t$RMarray[$i+1]->[3]\tnested\t$n\n";
				print NEST "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i-1]->[6]\t$len\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$RMarray[$i-1]->[12]\t$RMarray[$i-1]->[13]\t$RMarray[$i-1]->[1]\t$RMarray[$i-1]->[2]\t$RMarray[$i-1]->[3]\tnested\t$n\n";
				print NEST "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i-1]->[6]\t$len\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$RMarray[$i-1]->[12]\t$RMarray[$i-1]->[13]\t$RMarray[$i-1]->[1]\t$RMarray[$i-1]->[2]\t$RMarray[$i-1]->[3]\tnested\t$n\n";	
				print NEST "$RMarray[$i+1]->[4]\t$RMarray[$i+1]->[5]\t$RMarray[$i+1]->[6]\t$len\t$RMarray[$i+1]->[8]\t$RMarray[$i+1]->[9]\t$TEcf->{1}\t$RMarray[$i+1]->[11]\t$RMarray[$i+1]->[12]\t$RMarray[$i+1]->[13]\t$RMarray[$i+1]->[1]\t$RMarray[$i+1]->[2]\t$RMarray[$i+1]->[3]\tnested\t$n\n";	
				$n++;
			}	
		}				
	} elsif
	###########################################################
	# 2) no +2, check if only a +1 surrounding 
	###########################################################
	# IF [BBBBBB][AAAAAA][BBBBBB] 
	(($RMarray[$i-1]->[4] eq $RMarray[$i+1]->[4]) && ($RMarray[$i-1]->[9] eq $RMarray[$i+1]->[9]) && ($RMarray[$i-1]->[8] eq $RMarray[$i+1]->[8])
	&& ((($RMarray[$i]->[5] - $RMarray[$i-1]->[6]) <= 50) && (($RMarray[$i+1]->[5] - $RMarray[$i]->[6]) <= 50) 
	&& ((($RMarray[$i-1]->[8] eq '+') && (($RMarray[$i+1]->[11] - $RMarray[$i-1]->[12]) <= 20) && (($RMarray[$i+1]->[11] - $RMarray[$i-1]->[12]) >= -20)) 
	|| (($RMarray[$i-1]->[8] eq '-') && (($RMarray[$i-1]->[13] - $RMarray[$i+1]->[12]) <= 20) && (($RMarray[$i-1]->[13] - $RMarray[$i+1]->[12]) >= -20))))) {	
#								 -1		 0		+1											-1		 0	    +1
#		BUT, need to check IF [BBBBBB][AAAAAA][BBBBBB][XXXXXX][BBBBBB] OR [BBBBBB][XXXXXX][BBBBBB][AAAAAA][BBBBBB] => NO PRINT, because already detected when $i will be or was central [BBBBBB]
		unless (
		((($RMarray[$i+1]->[4] eq $RMarray[$i+3]->[4]) && ($RMarray[$i+1]->[9] eq $RMarray[$i+3]->[9]) && ($RMarray[$i+1]->[8] eq $RMarray[$i+3]->[8])) 
		&& ((($RMarray[$i+1]->[5] - $RMarray[$i+2]->[6]) <= 50) && (($RMarray[$i+3]->[5] - $RMarray[$i+2]->[6]) <= 50) 
		&& ((($RMarray[$i+1]->[8] eq '+') && (($RMarray[$i+3]->[11] - $RMarray[$i+1]->[12]) <= 20) && (($RMarray[$i+3]->[11] - $RMarray[$i+1]->[12]) >= -20)) 
		|| (($RMarray[$i+1]->[8] eq '-') && (($RMarray[$i+1]->[13] - $RMarray[$i+3]->[12]) <= 20) && (($RMarray[$i+1]->[13] - $RMarray[$i+3]->[12]) >= -20)))))
		|| 
		((($RMarray[$i-1]->[4] eq $RMarray[$i-3]->[4]) && ($RMarray[$i-1]->[9] eq $RMarray[$i-3]->[9]) && ($RMarray[$i-1]->[8] eq $RMarray[$i-3]->[8]))
		&& ((($RMarray[$i-2]->[5] - $RMarray[$i-3]->[6]) <= 50) && (($RMarray[$i-1]->[5] - $RMarray[$i-2]->[6]) <= 50) 
		&& ((($RMarray[$i-3]->[8] eq '+') && (($RMarray[$i+1]->[11] - $RMarray[$i-3]->[12]) <= 20) && (($RMarray[$i+1]->[11] - $RMarray[$i-3]->[12]) >= -20)) 
		|| (($RMarray[$i-3]->[8] eq '-') && (($RMarray[$i-3]->[13] - $RMarray[$i-1]->[12]) <= 20) && (($RMarray[$i-3]->[13] - $RMarray[$i-1]->[12]) >= -20)))))
		) {
			$nesting{'101'}++;
			$countnest{$TEc->{1}}{$TEc->{0}}++;
			
			# 1) print original lines
			print ORIG "$line_ref->{-1}\tnesting\t$n\n$line_ref->{0}\tnested\t$n\n$line_ref->{1}\tnesting\t$n\n";
			
			# 2) print new files
			my $Rend;
			# print [BBBBBB]
			if ($RMarray[$i-1]->[8] eq '+') {
				$Rend = $RMarray[$i+1]->[12];
			} else {
				$Rend = $RMarray[$i-1]->[12];
			}
			print NEWOUT "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i+1]->[6]\t$nesting{'len'}{1}\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$Rend\t$RMarray[$i+1]->[13]\t$nesting{'div'}{1}\t$nesting{'del'}{1}\t$nesting{'ins'}{1}\tnesting\t$n\n";
			print NEST "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i+1]->[6]\t$nesting{'len'}{1}\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$Rend\t$RMarray[$i+1]->[13]\t$nesting{'div'}{1}\t$nesting{'del'}{1}\t$nesting{'ins'}{1}\tnesting\t$n\n";
			# print nested TE
			# print [AAAAAA]
			print NEWOUT "$RMarray[$i]->[4]\t$RMarray[$i]->[5]\t$RMarray[$i]->[6]\t$len\t$RMarray[$i]->[8]\t$RMarray[$i]->[9]\t$TEcf->{0}\t$RMarray[$i]->[11]\t$RMarray[$i]->[12]\t$RMarray[$i]->[13]\t$RMarray[$i]->[1]\t$RMarray[$i]->[2]\t$RMarray[$i]->[3]\tnested\t$n\n";
			print NEST "$RMarray[$i]->[4]\t$RMarray[$i]->[5]\t$RMarray[$i]->[6]\t$len\t$RMarray[$i]->[8]\t$RMarray[$i]->[9]\t$TEcf->{0}\t$RMarray[$i]->[11]\t$RMarray[$i]->[12]\t$RMarray[$i]->[13]\t$RMarray[$i]->[1]\t$RMarray[$i]->[2]\t$RMarray[$i]->[3]\tnested\t$n\n";	
			$n++;
		}
	} elsif
	###########################################################
	#                                -1      0       +1      +2
	# 3) check particular case of [CCCCCC][BBBBBB][AAAAAA][CCCCCC]  
	###########################################################
	(($i < $#RMarray-2) && ($RMarray[$i-1]->[4] eq $RMarray[$i+2]->[4]) && ($RMarray[$i-1]->[9] eq $RMarray[$i+2]->[9]) && ($RMarray[$i-1]->[8] eq $RMarray[$i+2]->[8])
	&& ((($RMarray[$i]->[5] - $RMarray[$i-1]->[6]) <= 50) && (($RMarray[$i+1]->[5] - $RMarray[$i]->[6]) <= 50) && (($RMarray[$i+2]->[5] - $RMarray[$i+1]->[6]) <= 50) 
	&& ((($RMarray[$i-1]->[8] eq '+') && (($RMarray[$i+2]->[11] - $RMarray[$i-1]->[12]) <= 20) && (($RMarray[$i+2]->[11] - $RMarray[$i-1]->[12]) >= -20)) 
	|| (($RMarray[$i-1]->[8] eq '-') && (($RMarray[$i-1]->[13] - $RMarray[$i+2]->[12]) <= 20) && (($RMarray[$i-1]->[13] - $RMarray[$i+2]->[12]) >= -20))))) {
		
		my $nesting2x_len = ($RMarray[$i-1]->[6]-$RMarray[$i-1]->[5]+1) + ($RMarray[$i+2]->[6]-$RMarray[$i+2]->[5]+1);
		my $nesting2x_div = (($RMarray[$i-1]->[1] * ($RMarray[$i-1]->[6]-$RMarray[$i-1]->[5]+1)) + ($RMarray[$i+2]->[1] * ($RMarray[$i+2]->[6]-$RMarray[$i+2]->[5]+1)))/$nesting2x_len;
		my $nesting2x_del = (($RMarray[$i-1]->[2] * ($RMarray[$i-1]->[6]-$RMarray[$i-1]->[5]+1)) + ($RMarray[$i+2]->[2] * ($RMarray[$i+2]->[6]-$RMarray[$i+2]->[5]+1)))/$nesting2x_len;
		my $nesting2x_ins = (($RMarray[$i-1]->[3] * ($RMarray[$i-1]->[6]-$RMarray[$i-1]->[5]+1)) + ($RMarray[$i+2]->[3] * ($RMarray[$i+2]->[6]-$RMarray[$i+2]->[5]+1)))/$nesting2x_len;
		
		#			-2		-1      0       +1      +2	   +3
		# CHECK [DDDDDD][CCCCCC][AAAAAA][BBBBBB][CCCCCC][DDDDDD]
		if (($i > 2) && ($i < $#RMarray-2) 
		&& ($RMarray[$i-2]->[4] eq $RMarray[$i+3]->[4]) && ($RMarray[$i-2]->[9] eq $RMarray[$i+3]->[9]) && ($RMarray[$i-2]->[8] eq $RMarray[$i+3]->[8])
		&& ((($RMarray[$i-1]->[5] - $RMarray[$i-2]->[6]) <= 50) && (($RMarray[$i+3]->[5] - $RMarray[$i+2]->[6]) <= 50) 
		&& ((($RMarray[$i-2]->[8] eq '+') && (($RMarray[$i+3]->[11] - $RMarray[$i-2]->[12]) <= 20) && (($RMarray[$i+3]->[11] - $RMarray[$i-2]->[12]) >= -20)) 
		|| (($RMarray[$i-2]->[8] eq '-') && (($RMarray[$i-2]->[13] - $RMarray[$i+3]->[12]) <= 20) && (($RMarray[$i-2]->[13] - $RMarray[$i+3]->[12]) >= -20))))) { 
			
			my $nesting2x3_len = ($RMarray[$i-2]->[6]-$RMarray[$i-2]->[5]+1) + ($RMarray[$i+3]->[6]-$RMarray[$i+3]->[5]+1);
			my $nesting2x3_div = (($RMarray[$i-2]->[1] * ($RMarray[$i-2]->[6]-$RMarray[$i-2]->[5]+1)) + ($RMarray[$i+3]->[1] * ($RMarray[$i+3]->[6]-$RMarray[$i+3]->[5]+1)))/$nesting2x3_len;
			my $nesting2x3_del = (($RMarray[$i-2]->[2] * ($RMarray[$i-2]->[6]-$RMarray[$i-2]->[5]+1)) + ($RMarray[$i+3]->[2] * ($RMarray[$i+3]->[6]-$RMarray[$i+3]->[5]+1)))/$nesting2x3_len;
			my $nesting2x3_ins = (($RMarray[$i-2]->[3] * ($RMarray[$i-2]->[6]-$RMarray[$i-2]->[5]+1)) + ($RMarray[$i+3]->[3] * ($RMarray[$i+3]->[6]-$RMarray[$i+3]->[5]+1)))/$nesting2x3_len;
				
			$nesting{'CA01AC'}++;
			$countnest{$TEc->{-2}}{$TEc->{-1}}++;
			
			# 1) print original lines
			print ORIG "$line_ref->{-2}\tnesting\t$n\n$line_ref->{-1}\tnesting/nested\t$n\n$line_ref->{0}\tnested\t$n\n$line_ref->{1}\tnested\t$n\n$line_ref->{2}\tnesting/nested\t$n\n$line_ref->{3}\tnesting\t$n\n";
			
			# 2) print new files
			# print nesting TEs
			my $Rend;
			# print [DDDDDD]
			if ($RMarray[$i-2]->[8] eq '+') {
				$Rend = $RMarray[$i+3]->[12];
			} else {
				$Rend = $RMarray[$i-2]->[12];
			}
			print NEWOUT "$RMarray[$i-2]->[4]\t$RMarray[$i-2]->[5]\t$RMarray[$i+3]->[6]\t$nesting2x3_len\t$RMarray[$i-2]->[8]\t$RMarray[$i-2]->[9]\t$TEcf->{-2}\t$RMarray[$i-2]->[11]\t$Rend\t$RMarray[$i+3]->[13]\t$nesting2x3_div\t$nesting2x3_del\t$nesting2x3_ins\tnesting\t$n\n";
			print NEST "$RMarray[$i-2]->[4]\t$RMarray[$i-2]->[5]\t$RMarray[$i+3]->[6]\t$nesting2x3_len\t$RMarray[$i-2]->[8]\t$RMarray[$i-2]->[9]\t$TEcf->{-2}\t$RMarray[$i-2]->[11]\t$Rend\t$RMarray[$i+3]->[13]\t$nesting2x3_div\t$nesting2x3_del\t$nesting2x3_ins\tnesting\t$n\n";
			# print [CCCCCC]
			if ($RMarray[$i-1]->[8] eq '+') {
				$Rend = $RMarray[$i+2]->[12];
			} else {
				$Rend = $RMarray[$i-1]->[12];
			}
			print NEWOUT "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i+2]->[6]\t$nesting2x_len\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$Rend\t$RMarray[$i+2]->[13]\t$nesting2x_div\t$nesting2x_del\t$nesting2x_ins\tnesting/nested\t$n\n";
			print NEST "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i+2]->[6]\t$nesting2x_len\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$Rend\t$RMarray[$i+2]->[13]\t$nesting2x_div\t$nesting2x_del\t$nesting2x_ins\tnesting/nested\t$n\n";
			# print nested TEs
			print NEWOUT "$RMarray[$i]->[4]\t$RMarray[$i]->[5]\t$RMarray[$i]->[6]\t$len\t$RMarray[$i]->[8]\t$RMarray[$i]->[9]\t$TEcf->{0}\t$RMarray[$i]->[11]\t$RMarray[$i]->[12]\t$RMarray[$i]->[13]\t$RMarray[$i]->[1]\t$RMarray[$i]->[2]\t$RMarray[$i]->[3]\tnested\t$n\n";
			print NEWOUT "$RMarray[$i+1]->[4]\t$RMarray[$i+1]->[5]\t$RMarray[$i+1]->[6]\t$len\t$RMarray[$i+1]->[8]\t$RMarray[$i+1]->[9]\t$TEcf->{1}\t$RMarray[$i+1]->[11]\t$RMarray[$i+1]->[12]\t$RMarray[$i+1]->[13]\t$RMarray[$i+1]->[1]\t$RMarray[$i+1]->[2]\t$RMarray[$i+1]->[3]\tnested\t$n\n";
			print NEST "$RMarray[$i]->[4]\t$RMarray[$i]->[5]\t$RMarray[$i]->[6]\t$len\t$RMarray[$i]->[8]\t$RMarray[$i]->[9]\t$TEcf->{0}\t$RMarray[$i]->[11]\t$RMarray[$i]->[12]\t$RMarray[$i]->[13]\t$RMarray[$i]->[1]\t$RMarray[$i]->[2]\t$RMarray[$i]->[3]\tnested\t$n\n";
			print NEST "$RMarray[$i+1]->[4]\t$RMarray[$i+1]->[5]\t$RMarray[$i+1]->[6]\t$len\t$RMarray[$i+1]->[8]\t$RMarray[$i+1]->[9]\t$TEcf->{1}\t$RMarray[$i+1]->[11]\t$RMarray[$i+1]->[12]\t$RMarray[$i+1]->[13]\t$RMarray[$i+1]->[1]\t$RMarray[$i+1]->[2]\t$RMarray[$i+1]->[3]\tnested\t$n\n";
			$n++;
		} else {
			# no D surrounding
			$nesting{'A01A'}++;
			$countnest{$TEc->{-1}}{$TEc->{0}}++;
			$countnest{$TEc->{-1}}{$TEc->{1}}++ unless ($TEc->{1} eq $TEc->{0}); #trying to avoid counting multitude stuff that might be same element
			
			# 1) print original lines
			print ORIG "$line_ref->{-1}\tnesting\t$n\n$line_ref->{0}\tnested\t$n\n$line_ref->{1}\tnested\t$n\n$line_ref->{2}\tnesting\t$n\n";
				
			# 2) print new files
			my $Rend;
			if ($RMarray[$i-1]->[8] eq '+') {
				$Rend = $RMarray[$i+1]->[12];
			} else {
				$Rend = $RMarray[$i-1]->[12];
			}
			print NEWOUT "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i+2]->[6]\t$nesting2x_len\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$Rend\t$RMarray[$i+2]->[13]\t$nesting2x_div\t$nesting2x_del\t$nesting2x_ins\tnesting\t$n\n";
			print NEST "$RMarray[$i-1]->[4]\t$RMarray[$i-1]->[5]\t$RMarray[$i+2]->[6]\t$nesting2x_len\t$RMarray[$i-1]->[8]\t$RMarray[$i-1]->[9]\t$TEcf->{-1}\t$RMarray[$i-1]->[11]\t$Rend\t$RMarray[$i+2]->[13]\t$nesting2x_div\t$nesting2x_del\t$nesting2x_ins\tnesting\t$n\n";
			# print nested TEs
			print NEWOUT "$RMarray[$i]->[4]\t$RMarray[$i]->[5]\t$RMarray[$i]->[6]\t$len\t$RMarray[$i]->[8]\t$RMarray[$i]->[9]\t$TEcf->{0}\t$RMarray[$i]->[11]\t$RMarray[$i]->[12]\t$RMarray[$i]->[13]\t$RMarray[$i]->[1]\t$RMarray[$i]->[2]\t$RMarray[$i]->[3]\tnested\t$n\n";
			print NEWOUT "$RMarray[$i+1]->[4]\t$RMarray[$i+1]->[5]\t$RMarray[$i+1]->[6]\t$len\t$RMarray[$i+1]->[8]\t$RMarray[$i+1]->[9]\t$TEcf->{1}\t$RMarray[$i+1]->[11]\t$RMarray[$i+1]->[12]\t$RMarray[$i+1]->[13]\t$RMarray[$i+1]->[1]\t$RMarray[$i+1]->[2]\t$RMarray[$i+1]->[3]\tnested\t$n\n";
			print NEST "$RMarray[$i]->[4]\t$RMarray[$i]->[5]\t$RMarray[$i]->[6]\t$len\t$RMarray[$i]->[8]\t$RMarray[$i]->[9]\t$TEcf->{0}\t$RMarray[$i]->[11]\t$RMarray[$i]->[12]\t$RMarray[$i]->[13]\t$RMarray[$i]->[1]\t$RMarray[$i]->[2]\t$RMarray[$i]->[3]\tnested\t$n\n";
			print NEST "$RMarray[$i+1]->[4]\t$RMarray[$i+1]->[5]\t$RMarray[$i+1]->[6]\t$len\t$RMarray[$i+1]->[8]\t$RMarray[$i+1]->[9]\t$TEcf->{1}\t$RMarray[$i+1]->[11]\t$RMarray[$i+1]->[12]\t$RMarray[$i+1]->[13]\t$RMarray[$i+1]->[1]\t$RMarray[$i+1]->[2]\t$RMarray[$i+1]->[3]\tnested\t$n\n";	
			$n++;
		}	
	} else {
		#does not meet any condition of nesting -> just print line in the new output
		print NEWOUT "$RMarray[$i]->[4]\t$RMarray[$i]->[5]\t$RMarray[$i]->[6]\t$len\t$RMarray[$i]->[8]\t$RMarray[$i]->[9]\t$TEcf->{0}\t$RMarray[$i]->[11]\t$RMarray[$i]->[12]\t$RMarray[$i]->[13]\t$RMarray[$i]->[1]\t$RMarray[$i]->[2]\t$RMarray[$i]->[3]\tna\tna\n";
	}
}
close(NEWOUT);
close(ORIG);
close(NEST);

print LOG "
	$nesting{'101'} cases like [BBBBBB][AAAAAA][BBBBBB]
		(TE A nested in TE B)
	$nesting{'21012'} cases like [CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC]
		(TE A nested in TE B, nested in TE C)
	$nesting{'3210123'} cases like [DDDDDD][CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC][DDDDDD] 
		(TE A nested in TE B, nested in TE C, nested in TE D)
		
	$nesting{'2XXX2'} cases like [CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC]
		(TEs A, E and B are nested in TE C - A B E could be pieces of same TE badly masked, since they are close)
	$nesting{'32XXX23'} cases like [DDDDDD][CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC][DDDDDD]
		(TEs A, E and B are nested in TE C, nested in TE D - A B E could be pieces of same TE badly masked, since they are close)
		
	$nesting{'2X2X2'} cases like [CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC]
		(TEs X are 2 independent nesting events in TE C)
	$nesting{'32X2X23'} cases like [DDDDDD][CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC][DDDDDD]
		(TEs X are 2 independent nesting events in TE C, nesting in D)
		
	$nesting{'A01A'} cases like [CCCCCC][AAAAAA][BBBBBB][CCCCCC] 
		(TEs A and B are nested in TE C - A B could be pieces of same TE badly masked, since they are close)
	$nesting{'CA01AC'} cases like [DDDDDD][CCCCCC][AAAAAA][BBBBBB][CCCCCC][DDDDDD]
		(TEs A and B are nested in TE C, nested in TE D - A B could be pieces of same TE badly masked, since they are close)\n\n";

print LOG "COUNTING OF NESTING EVENTS:\nNestingTE_class\tNestingTE_classes and nb of events:\n";
foreach my $nestingTE (keys %countnest) {
	print LOG "\n$nestingTE\t";
	foreach my $nestedTE (keys %{$countnest{$nestingTE}}) {
		print LOG "\t$nestedTE\t$countnest{$nestingTE}{$nestedTE}\n\t";
	}
}
exit;

##################################################################################################################################################################
# SUBROUTINES
##################################################################################################################################################################
######################################################
# QUITE GENERAL PARSE RM STUFF
######################################################
#----------------------------------------------------------------------------
# get TE infos
#----------------------------------------------------------------------------
sub get_TEs_infos {
	my $input = shift; #file name
	my %TEs = ();
	open(my $input_fh, "<", $input) or die print "ERROR: could not open $input!\n";
	LINE: while(<$input_fh>) {
		chomp (my $line = $_);
		next LINE if ($line !~ /\w/);
		my @TEs = split('\t', $line); 
		my $lcRname = lc ($TEs[0]);
		$TEs{$lcRname} = \@TEs;
	}	
	return \%TEs;
}

#----------------------------------------------------------------------------
# get line_ref class and classfam array
# ($line_ref->{$c},$TEcf,$TEc) = get_line_and_TEc(\@RMarray,$c,$col,$i,$TEc,$TEcf,$TE,$TEclass,$TEage);
#----------------------------------------------------------------------------
sub get_line_and_TEc {
	my ($RMarray,$c,$col,$i,$TEc,$TEcf,$TE,$TEclass,$TEage) = @_;
	my $line_ref = "";
		
	#avoid undef stuff
	return ($line_ref,$TEc,$TEcf) if (($i < abs($c)) || ($c >= ($#RMarray - $i)));
	
	#if OK carry on
	#modify $i to get the relevant info (from -3 to +3)
	$i = $i + $c;
		
	#line ref
	for (my $j = 0; $j <= $col ; $j++){
		$line_ref = "$line_ref\t$RMarray->[$i][$j]";
	}
	#classfam
	my $Rname = lc($RMarray->[$i][9]);
	my $classfam = $RMarray->[$i][10];
	if ((($TEclass) || ($TEage)) && ($TE->{$Rname}->[1])) {
		$TEcf->{$c} = $TE->{$Rname}->[1]."/".$TE->{$Rname}->[2];
		$TEc->{$c} = $TE->{$Rname}->[1];
	} else {
		$TEcf->{$c} = $classfam;
		$TEc->{$c} = $classfam;
		$TEc->{$c} =~ s/(.*)\/.*/$1/;
	}
	return ($line_ref,$TEc,$TEcf)
}	


