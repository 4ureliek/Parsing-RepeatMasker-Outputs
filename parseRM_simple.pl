#!/usr/bin/perl -w
#######################################################
# Author  :  Aurelie K
# date    :  04 May 2011 = v1.0 [other script, updated until 1.9]
# version :  2.3
# email   :  4urelie.k@gmail.com
# Pupose  :  Parse RepeatMasker outputs to get summary info for each repeat as well as masked amounts by class and family.
#			 For the repeats, it will provide:
#				- fragment number (frg nb): all, complete frg nb (from start to end), frg nb corrected for interupted repeats
#				- %div, ins, del: pondered average, median
#				- length masked + %genome by this repeat
#				- amount of DNA that is masked several times (usually 2) by this element 
#			 Plus the script provides a set of files with overlap info to help giving real amounts/%
# Steps   :
# 		- calcul of genome length, with an optional step to remove the Ns in sequences (if the genome used was already without Ns, this does not matter, just a useless step)
# 		- calcul of length and % of genome masked by each TE
# 		- evaluate percentage of genome masked two times (overlaps) -> details in files
#######################################################
# Updated : v2.0, 03-04 Oct 2013
#		- possibility of providing a "TEclass" file with TE names, class and family as first 3 columns [additional columns can be there but won't be used]
#		- corrected error at a if statement that prevented writing overlap infos in all-repeats output file
#		- calculation of length was missing a +1 (2 spots)
#		- OVERLAP info are more complete (it was not giving all overlap by DNA class for ex, but only by DNA class when family was different)
#		- add class and family amount counts in a summary file
#		- changed getting options [! the command line is different]
#		- removed concatenation of inputs, it's useless, cat is easy to do
#		- added a log file (-> to keep with output files)
#		- added follow up of version used by writing it in the log file
#		- changed names of output files
#		- changed overlap_all and diff output columns => easier to filter or check in excel
# Updated : v2.1, 29 Oct 2013
#		- just changed $Rmout reading, from: 
#				unless ($line =~ /perc|[sS]core|Gname/) 
#				to
#				unless ($line =~ /position|[sS]core|Gname/)
#		   Because perc is also included in "supercontig"...
# Updated : v2.2, 05 Nov 2013
#		- optional addition of a column that will contain the length of the consensus repeat used to mask (if library is provided)
#		- added some values in the summary file, gain of time not having to do it in excel
# Updated : v2.3, 19 Nov 2013
#		- Fixed bug, Class/Fam column was from RM output and not from the TEclass file


### TO DO


# Updated : v2.4
#		- re write with subs
#		- create nonTE when possible => filtered out in GetLandscape
#		- improve class det (ex. if "hAT" in name, class = DNA)
#		- Option -TEageOUT to create output file readily usable as a TEage file for parseRM_GetLandscape
# 			 -TEageOUT
#         	 => create output file readily usable as a TEclass file for parseRM_GetLandscape (e.g. Rname \\t Rclass \\t Rfamily \\t Rclass/Rfam \\t avg_%div \\t avg_%div_bin)
######################################################

use strict;
use warnings;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Getopt::Long;
use Array::Utils qw(:all);

my $version = "v2.3";
my $usage = "\nUsage [$version]: 
    perl <scriptname.pl> -genfile <FastaFile> [-noNrem] -RMout <FastaFile.out> [-TEs <TEclass>] [-lib <repeat_library>] [-fast] [-norw]
    OR
    perl <scriptname.pl> -genlen LengthOfGenome(nt) -RMout <Genome.out> [-TEs <TEclass>] [-lib <repeat_library>] [-fast] [-norw]
	
	Order of options/arguments doesn't matter. 
	
    WITH MANDATORY ARGUMENTS:
    -genfile <FastaFile> OR -genlen LengthOfGenome(nt)
        => <GenomeFile> = fasta file THAT WAS MASKED, typically a genome. If Ns were NOT removed from it before masking, also add the option -noNrem
		=> LengthOfGenome = amount of nucleotides in the genome (e.g. genome size in nt; it's the only info that is required from the genome file)
    -RMout <FastaFile.out>
        => <FastaFile.out> = Repeat masker output file (typically, Genome.out).
           If several .out files to parse, concatenate them before using this script, using following command:
           cat file1.out file2.out file3.out > FastaFile.out
           OR, if all .out in current folder are to parse:
           cat *.out > FastaFile.out
	 
    OPTIONAL ARGUMENTS (flagged by brackets [...] around it)
    -TEs <TEclass>
        => optional input file, with TE information as follow, 3 first columns mandatory: Rname \\t Rclass \\t Rfam \\t Rclass/Rfam \\t etc (only first 3 columns will be used)
	       an easy way to get it is to run this script a first time, copy the first columns of the output all-repeats.tab, modify them accordingly and then copy this in a text file
    -lib <repeat_library> 
        => <repeat_library> = optional input file corresponding to the library used to mask. It will allow printing length of consensus sequences in the output.
    -norw
        => do not rewrite in a file the whole output, with additional length column. Will gain space.
    -fast
        => do not calculate averages and medians for %div, %ins, %del, Len_masked.
    -noNrem
        => do not check for N removal in genome (chose that if repeat masker was done on a genome with Ns, e.g. with assembly gaps)
	 
    Note that parsing won't be as relevant if names are not formatted consistently (Rname#Rclass/Rfam)
      the use of a TEclass file can help correcting that without changing the repeat library itself
      (files XXX.parseRM.all-repeats.tab and XXX.parseRM.overlaps-all.tab only will be the only really useful ones if class/fam are not right)\n\n";

################################################################################################################################################################
# Get arguments, options and paths
################################################################################################################################################################
my $GenFile;
my $GenLen;
my $RMout;
my $TEs;
my $fast;
my $norw;
my $noNrem;
my $repfile;
GetOptions ('fast' => \$fast, 'norw'=> \$norw, 'noNrem'=> \$noNrem, 'genfile=s' => \$GenFile, 'genlen=s' => \$GenLen, 'RMout=s' => \$RMout, 'TEs=s' => \$TEs, 'lib=s' => \$repfile);

#check step to see if mandatory arguments are provided and if OK
if (((! $GenFile) && (! $GenLen)) || (! $RMout) || (($GenLen) && ($GenLen =~ /\D/))){
	die $usage;
}

#get path of RM output file, to be used as output location after
my $path = $RMout;
if ($path =~ m/^.*\/.*/){
	$path =~ s/(.*)\/.*$/$1/; #get the path if there are directories in the argument
} else {
	$path = "."; #if no "/" then path = . (ie current folder)
}

#Create folder to parse results
my $filetmp = $RMout;
$filetmp =~ s/.*\/(.*)$/$1/; #extract name of file
my $j=1;
until (!-d "$path/$filetmp.parsed$j"){
	$j++;
}
mkdir "$path/$filetmp.parsed$j";
my $newpath = "$path/$filetmp.parsed$j";

#Log file
my $log = "$path/$filetmp.parseRM.$j.log";
open LOG, ">$log" or die "\t !! Failed to create file $log $!\n";

#print in log details of pipeline
print LOG "\n----------------------------\nPipeline started ($version):\n----------------------------\n\t-> analysis of file = $RMout\n";
print LOG "\t-> with genome file = $GenFile\n" if ($GenFile);
print LOG "\t-> with genome length = $GenLen\n" if ($GenLen);
print LOG "\t-> with TE names, class and fam info = $TEs\n" if ($TEs);
print LOG "\t-> with library fasta file = $repfile\n" if ($repfile);
print LOG "\t-> option -fast chosen  => no stats\n" if ($fast);
print LOG "\t-> option -norw chosen  => do not rewrite RM output with additional length column\n" if ($norw);
print LOG "\t-> option -nonrem chosen  => do not remove Ns / check for N removal from genome\n" if ($noNrem);

#get Rclass/fam + %div and age if relevant
my %TEs = ();
if ($TEs) {
	open(TE, "<$TEs") or die print " !! $TEs defined but can't be opened $!\n";
	while(<TE>) {
		chomp(my $te = $_);
		my @TEs = split('\t', $te); 
		# 0       1       2         3        4       5			6
		#Rname	Rclass	Rfam	Rclassfam	%div	~Age	Acient/LineageSpe
		my $lcRname = lc ($TEs[0]);
		$TEs{$lcRname} = \@TEs;
	}
	close (TE);
}

#get repeats length if relevant
my %rep_len = ();
if ($repfile) {
	my $lib = Bio::SeqIO->new(-file => $repfile, -format => "fasta") or die print " !! $repfile defined but can't be opened $!\n";
	while( my $seq = $lib->next_seq() ) {
		my $id = $seq->display_id;
		my $len = $seq->length;
		$rep_len{$id} = $len;
	}
}


################################################################################################################################################################
# Calculate length of the genome to have the total length (to get percentages)
# (includes check steps to avoid repeating this if length already calculated)
# ONLY IF OPTION -genfile
################################################################################################################################################################
if ($GenFile) {
	print LOG "\n--- Getting length of the genome from $GenFile...\n";
	my $genome;
	#removing the N in the genome if relevant
	unless ($noNrem) {
		my $genometmp = "$GenFile.Nrem";
		$genome = "$GenFile.Nrem.fa";
		if (-e $genome) {
			print LOG "    N already removed from $GenFile ($genome exists) - skipping N removal step...\n";
		} else {
			print LOG "    N may be not removed from $GenFile ($genome does not exists) - Removing N...\n";
			my $Gentmp = Bio::SeqIO->new(-file => $GenFile, -format => "fasta") or die print LOG "\t !! Failed to create Bio::SeqIO object from $GenFile $!\n";
			open GMINN, ">$genometmp" or die print LOG "\t !! Failed to create $genometmp file $!\n";
			while( my $seq = $Gentmp->next_seq() ) {
				my $sequence = $seq->seq;
				$sequence =~ s/[Nn]//g;
				my $id = $seq->display_id."\t".$seq->desc;
				print GMINN ">$id\n$sequence\n";
			}
			close GMINN;
			# Rewrite sequences in fasta format just to be sure
			my $Gen = Bio::SeqIO->new(-file => $genometmp, -format => "fasta") or die print LOG "\t !! Failed to create Bio::SeqIO object from $genometmp $!\n";
			my $GenMinusN = Bio::SeqIO->new(-file => ">$genome", -format => "fasta") or die print LOG "\t !! Failed to create Bio::SEQIO outputfile $genome $!\n";
			while( my $seq2 = $Gen->next_seq() ) {
				$GenMinusN->write_seq($seq2);		
			}
		}
	} else {
		$genome = $GenFile;
	}
	
	# index the genome and connect to the fasta file with N removed
	my $reindex;
	my $indexfile = "$genome.index";
	if (-e $indexfile) {
		$reindex = 0;
		print LOG "    Genome previously indexed ($indexfile exists) - Skipping indexing step...\n";
	} else {
		$reindex = 1;
		print LOG "    Genome not indexed ($indexfile does not exists) - Indexing...\n";
	}
	my $db = Bio::DB::Fasta->new($genome, -reindex=>$reindex) or die print LOG "\t !! Failed to create Bio::DB::Fasta object from $genome $!\n";

	#create list of the ID of the genome file
	my @dbIDs = $db->get_all_ids();

	#calculating total length of the genome (of the database)
	my $lengthfile = "$genome.length";
	
	unless (-e $lengthfile) {
		print LOG "    Total length of genome not known ($lengthfile does not exists) - Calculating total length...";
		open LENGTH, ">$lengthfile" or die print LOG "\t !! could not create file $lengthfile $!\n";
		$GenLen = 0;
		foreach my $ID (@dbIDs) {
			my $obj = $db->get_Seq_by_id($ID);
			my $len = $obj->length;
			$GenLen += $len;
		}
		print LENGTH $GenLen;
		close (LENGTH);
	} else {
		print LOG "    Total length of genome has been previously calculated ($lengthfile exists)";
	}
	
	open LENGTH, "<$lengthfile" or die print LOG "\t !! could not open $lengthfile $!\n";
	while (<LENGTH>) {
		$GenLen = $_;
	}	
	print LOG " = $GenLen nucleotides\n";
	close(LENGTH);
}



################################################################################################################################################################
# Parsing the RM output
################################################################################################################################################################
print LOG "\n--- Parsing Repeatmasker output (this can be long)...\n";

#open the RMout file and put it into a hash
open INFILE, "<$RMout" or die print LOG "\t !! could not open $RMout $!\n";
my @RMout = ();
while(<INFILE>) {
	chomp(my $line = $_);
	unless ($line =~ /position|[sS]core|Gname/){ #remove column names
		if ($line =~ /\w/){ #remove blank lines
			$line =~ s/^\s+//; #remove spaces in beginning of lines
			push(@RMout, $line);
		}
	}
}
close(INFILE);


# Declare variables needed to parse
my %rep_len_fullID = ();
my %class = ();
my %fam = ();
my %elem = ();
my %blocks = ();
my %double = ();
my %doubleE = ();
my %doubleF = ();
my %doubleFnotE = ();
my %doubleC = ();
my %doubleCnotF = ();
my %doubleD = ();
my %length = ();
my %ListLen = ();
my %nb = ();
my %Cnb = ();
my %div = ();
my %ListDiv = ();
my %del = ();
my %ListDel = ();
my %ins = ();
my %ListIns = ();
my %errorTE = ();
my $totalmasked;
my $doublemasked;
my $totaldouble;
my $MaskedDiff;
my $MaskedSameF;
my $MaskedSameFnotE;
my $MaskedSameE;
my $MaskedSameC;
my $MaskedSameCnotF;
my $output;

#prepare output with all lines rewritten, except if -fast option
unless ($norw){
	$output = "$newpath/$filetmp.length.tab";
	open FILE, ">$output" or die print LOG "\t !! could not open $output $!\n";
	print FILE "Score\t%div\t%del\t%ins\tGname\tGstart\tGend\tGleft\tGmasked\tStrand\tRname\tRclassfam\tRclass\tRfam\tRfullname\tRstart\tRend\tRleft\tblock\tIfOverlap\n\n";
}

for (my $i = 0; $i <= $#RMout; $i++){ #read the input file line by line
	my ($score,$div,$del,$ins,$Gname,$Gstart,$Gend,$Gleft,$Strand,$Rname,$Rclassfam,$Repstart,$Rend,$Repleft,$block,$IfOverlap) = split(/\s+/, $RMout[$i]);
	#now deal with repeat names, best is $TEs file but otherwise at least try to get Rclass/Rfam fom the RMout
	my ($Rclass,$Rfam);
	my $lcRname = lc ($Rname);
	
	if ($TEs{$lcRname}->[0]) { #ie if TE is in the list of TEclass
		# 0       1       2         3
		#Rname	Rclass	Rfam	Rclassfam
		print " => in -TEs";
		$Rclass = $TEs{$lcRname}->[1];
		$Rfam = $TEs{$lcRname}->[2];
		$Rclassfam = "$Rclass/$Rfam";
	} else {
		print LOG "\t! $Rname is not in $TEs\n" if (($TEs) && (! $errorTE{$Rname})); #eg $TEs{$lcname} not defined despite $TEs file in input
		$errorTE{$Rname} = 1 if ($TEs);
		if ($Rclassfam =~ /\//) {
			($Rclass,$Rfam) = split(/\//, $Rclassfam);
		} else {
			$Rfam = $Rclassfam;
			$Rfam=~ s/^(.*)\..*$/$1/;
			$Rclass = $Rclassfam;
			$Rclass =~ s/^.*\.(.*)$/$1/;
		}
	}
	#make sure key of the hash will be unique; $Gname\t$Gstart in the beginning allows to sort by query name + order (=input)
	my $uniqname = "$Gname\t$Gstart\t$Gend\t$Rname\t$Rclass\t$Rfam"; #[v2.0]
	
	#calculate length in genome masked by this fragment
	my $Gmasked = ($Gend - $Gstart + 1); #[v2.0 add the +1]
	$totalmasked+=$Gmasked;
	
	#check if there is the * for overlaps
	if ($IfOverlap) {
		$IfOverlap =~ s/\*/yes/;			
	} else {
		$IfOverlap = ".";
	}
	
	#Storing datas for each element (Name#Rclass_Rfam)
	my $Rfullname = "$Rname"."#"."$Rclass"."/"."$Rfam"; #[v2.0] => put back the /
	my $key = "$Rname\t$Rclass\t$Rfam\t$Rclassfam\t$Rfullname";
	
	#get length of sequence if relevant
	if ($repfile) {
		my $Rlen = $rep_len{$Rname."#".$Rclassfam};
		$rep_len_fullID{$key} = $Rlen;
	}	
	
	#counting length masked (add length masked)
	$length{$key} += $Gmasked;
	$class{$Rclass} += $Gmasked;
	$fam{$Rclassfam} += $Gmasked;
	#counting fragments number (add 1 each time)
	$nb{$key} ++;			
	
	#Getting hash to count nr fragments (ie if same TE fragmented) 
	if ($blocks{$key}) {
		$blocks{$key}="$blocks{$key}#$Gname-$block";
	} else {
		$blocks{$key}="$Gname-$block";
	}

	#Counting nb of 1 to end fragments
	my $Rstart;
	my $Rleft;
	if ($Strand eq '+'){
		$Rstart = $Repstart;
		$Rleft = $Repleft;
	} else {
		$Rstart = $Repleft;
		$Rleft = $Repstart;
	}
	$Rleft =~ s/\(//;
	$Rleft =~ s/\)//;	
	if (($Rstart == 1) && ($Rleft == 0)) {
		$Cnb{$key} ++;
	}
	
	#Do the maths unless option fast
	unless ($fast){
		$div{$key} += ($div*$Gmasked);		#adding values of div percentage, x length to ponderate 
		$del{$key} += ($del*$Gmasked);		#adding values of del percentage, x length to ponderate 
		$ins{$key} += ($ins*$Gmasked);		#adding values of ins percentage, x length to ponderate 

		#putting values in hashes to be able to calculate medians further by splitting again into lists
		if (exists $ListLen{$key}) {
			$ListLen{$key} = "$ListLen{$key} $Gmasked";
		} else {
			$ListLen{$key} = $Gmasked;
		}	

		if (exists $ListDiv{$key}) {
			$ListDiv{$key} = "$ListDiv{$key} $div";
		} else {
			$ListDiv{$key} = $div;
		}	

		if (exists $ListDel{$key}) {
			$ListDel{$key} = "$ListDel{$key} $ins";
		} else {
			$ListDel{$key} = $del;
		}

		if (exists $ListIns{$key}) {
			$ListIns{$key} = "$ListIns{$key} $del";
		} else {
			$ListIns{$key} = $ins;
		}
	}
	
	
	##########################################################################################################
	#Take care of overlap
	unless ($i == 0) {
		my ($scoreU,$divU,$delU,$insU,$GnameU,$GstartU,$GendU,$GleftU,$StrandU,$RnameU,$RclassfamU,$RepstartU,$RendU,$RepleftU,$blockU,$IfoverlapU,$IfPrevCutU,$NboffrgU) = split(/\s+/, $RMout[$i-1]);
		
		#now deal with repeat names, best is $TEs file but otherwise at least try to get Rclass/Rfam fom the RMout
		my ($RclassU,$RfamU);
		my $lcRnameU = lc ($RnameU);
		if ($TEs{$lcRnameU}->[0]) { #ie if TE is in the list of TEclass
			# 0       1       2         3
			#Rname	Rclass	Rfam	Rclassfam
			$RclassU = $TEs{$lcRnameU}->[1];
			$RfamU = $TEs{$lcRnameU}->[2];
		} else {
			print LOG "\t! $RnameU is not in $TEs\n" if (($TEs) && (! $errorTE{$RnameU})); #eg $TEs{$lcnameU} not defined despite $TEs file in input
			$errorTE{$RnameU} = 1 if ($TEs);
			if ($RclassfamU =~ /\//) {
				($RclassU,$RfamU) = split(/\//, $RclassfamU);
			} else {
				$RfamU = $RclassfamU;
				$RfamU =~ s/^(.*)\..*$/$1/;
				$RclassU = $RclassfamU;
				$RclassU =~ s/^.*\.(.*)$/$1/;
			}
		}	
		my $RfullnameU = $RnameU."#".$RclassU."/".$RfamU; #[v2.0] => put back the / + fixed name; it was $Rname instead of $RnameU
		if (($Gname eq $GnameU) && ($Gstart < $GendU) && ($Gstart > $GstartU)){	#ie same sequence of genome + there is overlapping but still progression of Gstart
			$doublemasked = ($GendU - $Gstart +1); # length which is masked twice for this element (Name#Class_Fam) #[v2.0] add the +1
			$totaldouble += $doublemasked; # total value
			$double{$uniqname} = "$GstartU\t$GendU\t$RnameU\t$RclassU\t$RfamU\t$doublemasked";	#hash to list all overlapping #[v2.0]
			
			if ($Rfullname eq $RfullnameU){ #same element => will be in stat output file
				$doubleE{$key} += $doublemasked;
				$MaskedSameE += $doublemasked;
			} 
			
			if ($Rclassfam eq $RclassfamU) { #same class/family [will include same elements]
				$doubleF{$Rclassfam} += $doublemasked;
				$MaskedSameF += $doublemasked;
			} 
			
			if ($Rclass eq $RclassU) {	#same class (LTR, DNA etc) [will include same elements, same families]
				$doubleC{$Rclass} += $doublemasked;
				$MaskedSameC += $doublemasked;
			} 
			
			if (($Rclass eq $RclassU) && ($Rclassfam ne $RclassfamU)) {	#same class (LTR, DNA etc) but NOT same family
				$doubleCnotF{$Rclass} += $doublemasked;
				$MaskedSameCnotF += $doublemasked;
			} 
			
			if (($Rclassfam eq $Rclassfam) && ($Rfullname ne $Rfullname)) {	#same class/family but NOT same element
				$doubleFnotE{$Rclassfam} += $doublemasked;
				$MaskedSameFnotE += $doublemasked;
			} 
			
			#get overlaps done by different repeats and different classes
			unless (($Rfullname eq $RfullnameU) || ($Rclassfam eq $RclassfamU) || ($Rclass eq $RclassU)) { #what's left basically, and class is different
				$doubleD{$uniqname} = "$GstartU\t$GendU\t$RnameU\t$RclassU\t$RfamU\t$doublemasked"; #[v2.0]
				$MaskedDiff += $doublemasked;	
			}
		}
	}
	
	
	##########################################################################################################
	#rewrite RM ouput file with the addition of length and a column with Fullname (easier to sort by query names in excel)
	#and use this step to separate columns with tabs, and introduce Fullname (=Query name with Name#Rclassfam)
	unless ($norw){
		print FILE "$score\t$div\t$del\t$ins\t$Gname\t$Gstart\t$Gend\t$Gleft\t$Gmasked\t$Strand\t$Rname\t$Rclassfam\t$Rclass\t$Rfam\t$Rfullname\t$Repstart\t$Rend\t$Repleft\t$block\t$IfOverlap\n";
		#NB if IfOverlapp = yes, ie there is a higher-scoring match whose domain partly (<80%) includes the domain of this match. 	
	}
}
close FILE;



##########################################################################################################
# Read hashes and write in files
##########################################################################################################
# ALL ELEMENTS
###################################################
my $stats = "$newpath/$filetmp.parseRM.all-repeats.tab";
my $totalmaskedper;

my %nrCounts = ();
foreach my $TE (sort keys %blocks) {
	my @blocks = split ("#",$blocks{$TE});
	my @unique_blocks = unique(@blocks);
	$nrCounts{$TE} = $#unique_blocks + 1;
}
	
unless ($fast){
	# File of data per element of query
	my @ListLen;
	my %hashMedLen = ();
	foreach my $klen (sort keys %ListLen) {
		@ListLen = split (/\s+/, $ListLen{$klen});
		$hashMedLen{$klen} = &median(\@ListLen);
	}

	my @ListDiv;
	my %hashMedDiv = ();
	foreach my $kdiv (sort keys %ListDiv) {
		@ListDiv = split (/\s+/, $ListDiv{$kdiv});
		$hashMedDiv{$kdiv} = &median(\@ListDiv);
	}

	my @ListDel;
	my %hashMedDel = ();
	foreach my $kdel (sort keys %ListDel) {
		@ListDel = split (/\s+/, $ListDel{$kdel});
		$hashMedDel{$kdel} = &median(\@ListDel);
	}

	my @ListIns;
	my %hashMedIns = ();
	foreach my $kins (sort keys %ListIns) {
		@ListIns = split (/\s+/, $ListIns{$kins});
		$hashMedIns{$kins} = &median(\@ListIns);
	}
	
	open STAT, ">$stats" or die print LOG "\t !! could not open $stats $!\n";
	if ($repfile) {
		print STAT "Rname\tRclass\tRfam\tRclassfam\tRfullname\tRlen\tFRG_NB\tFRG_NB_StartToEnd\tNR_FRG_NB\tAVG_%DIV\tMED_%DIV\tAVG_%DEL\tMED_%DEL\tAVG_%INS\tMED_%INS\tLEN_MASKED\tAVG_LEN_MASKED\tMED_LEN_MASKED\t%_GENOME\tLEN_OVERLAP\t%_OVERLAP_(GENOME)\t%_OVERLAP_(LEN_MASKED)\n\n";
	} else {
		print STAT "Rname\tRclass\tRfam\tRclassfam\tRfullname\tFRG_NB\tFRG_NB_StartToEnd\tNR_FRG_NB\tAVG_%DIV\tMED_%DIV\tAVG_%DEL\tMED_%DEL\tAVG_%INS\tMED_%INS\tLEN_MASKED\tAVG_LEN_MASKED\tMED_LEN_MASKED\t%_GENOME\tLEN_OVERLAP\t%_OVERLAP_(GENOME)\t%_OVERLAP_(LEN_MASKED)\n\n";
	}	
	foreach my $keyStats (sort keys %length) {
		my $percentageS = $length{$keyStats} / $GenLen * 100;
		$totalmaskedper+=$percentageS;
		my $meanL = $length{$keyStats} / $nb{$keyStats};
		my $meanDiv = $div{$keyStats} / $length{$keyStats};
		my $meanDel = $del{$keyStats} / $length{$keyStats};
		my $meanIns = $ins{$keyStats} / $length{$keyStats};
		$Cnb{$keyStats}=0 unless (defined $Cnb{$keyStats});	
		if ($repfile) {
			$rep_len_fullID{$keyStats} = "na" unless ($rep_len_fullID{$keyStats});
			print STAT "$keyStats\t$rep_len_fullID{$keyStats}\t$nb{$keyStats}\t$Cnb{$keyStats}\t$nrCounts{$keyStats}\t$meanDiv\t$hashMedDiv{$keyStats}\t$meanDel\t$hashMedDel{$keyStats}\t$meanIns\t$hashMedIns{$keyStats}\t$length{$keyStats}\t$meanL\t$hashMedLen{$keyStats}\t$percentageS\t";
		} else { 
			print STAT "$keyStats\t$nb{$keyStats}\t$Cnb{$keyStats}\t$nrCounts{$keyStats}\t$meanDiv\t$hashMedDiv{$keyStats}\t$meanDel\t$hashMedDel{$keyStats}\t$meanIns\t$hashMedIns{$keyStats}\t$length{$keyStats}\t$meanL\t$hashMedLen{$keyStats}\t$percentageS\t";
		}

		#Now check if for this query (element#Rclassfam) there is a value associated for a overlapping
		if (exists $doubleE{$keyStats}){ #[v2.0]
			my $EachQDblPerG = $doubleE{$keyStats} / $GenLen * 100;
			my $EachQDblPerE = $doubleE{$keyStats} / $length{$keyStats} * 100;
			print STAT "$doubleE{$keyStats}\t$EachQDblPerG\t$EachQDblPerE\n";
		} else {
			print STAT"\t\t\n";
		}
	}
close STAT;
} elsif ($fast){
	$stats = "$newpath/$filetmp.parseRM.all-repeats.light.tab";
	open STATF, ">$stats" or die print LOG "	 !!could not open $stats $!\n";
	if ($repfile) {
		print STATF "Rname\tRclass\tRfam\tRclassfam\tRfullname\tRlen\tFRG_NB\tFRG_NB_StartToEnd\tNR_FRG_NB\tAVG_LEN_MASKED\tMED_LEN_MASKED\t%_GENOME\tLEN_OVERLAP\t%_OVERLAP_(GENOME)\t%_OVERLAP_(LEN_MASKED)\n\n";
	} else {
		print STATF "Rname\tRclass\tRfam\tRclassfam\tRfullname\tFRG_NB\tFRG_NB_StartToEnd\tNR_FRG_NB\tAVG_LEN_MASKED\tMED_LEN_MASKED\t%_GENOME\tLEN_OVERLAP\t%_OVERLAP_(GENOME)\t%_OVERLAP_(LEN_MASKED)\n\n";
	}
	foreach my $keyStatsF (sort keys %length) {
		my $percentageF = $length{$keyStatsF} / $GenLen * 100;
		$totalmaskedper+=$percentageF;
		my $meanLF = $length{$keyStatsF} / $nb{$keyStatsF};
		$Cnb{$keyStatsF}=0 unless ($Cnb{$keyStatsF});
		if ($repfile) {
			$rep_len_fullID{$keyStatsF} = "na" unless ($rep_len_fullID{$keyStatsF});
			print STATF "$keyStatsF\t$rep_len_fullID{$keyStatsF}\t$nb{$keyStatsF}\t$Cnb{$keyStatsF}\t$nrCounts{$keyStatsF}\t$length{$keyStatsF}\t$meanLF\t$percentageF\t";
		} else {
			print STATF "$keyStatsF\t$nb{$keyStatsF}\t$Cnb{$keyStatsF}\t$nrCounts{$keyStatsF}\t$length{$keyStatsF}\t$meanLF\t$percentageF\t";
		}
		#Now check if for this query (element#Rclassfam) there is a value associated for a overlapping
		if (exists $doubleE{$keyStatsF}){ #[v2.0]
			my $EachQDblPerGF = $doubleE{$keyStatsF} / $GenLen * 100;
			my $EachQDblPerQF = $doubleE{$keyStatsF} / $length{$keyStatsF} * 100;
			print STATF "$doubleE{$keyStatsF}\t$EachQDblPerGF\t$EachQDblPerQF\n";
		} else {
			print STATF"\n";
		}
	}
close STATF;
}

###################################################
# SUMMARY FILE WITH OVERLAP VALUES
###################################################
# output files:
my $DoubleStatFile = "$newpath/$filetmp.parseRM.Summary.tab";
my $lengthDouble = "$newpath/$filetmp.parseRM.overlaps-all.tab";
my $lengthDoubleD = "$newpath/$filetmp.parseRM.overlaps-if-diff-TEs.tab";

if (defined $totaldouble){
	# STAT FILE
	open DOUBLESTATS, ">$DoubleStatFile" or die print LOG "\t !! could not open $DoubleStatFile $!\n";
	print DOUBLESTATS "#This file gives summary info for total, by class and by class/fam\n#See file $filetmp.parseRM.all-repeats.light.tab for more info about repeats one by one\n#Overlap or double corresponds to DNA fragments that are masked by several elements. These amounts need to be subtracted in order to get more accurate TE amount\n#In file $filetmp.parseRM.overlaps-all.tab, all overlaps are listed.\n#In file $filetmp.parseRM.overlaps-if-diff-TEs.tab, all overlaps by DIFFERENT CLASSES are listed.\n#Note that these files may not be parsed correctly if names are not formatted consistently (Rname#Rclass/Rfam)\n\n";
	
	#ALL
	my $totdblper = $totaldouble / $GenLen * 100;
	my $tot_nodouble = $totalmasked - $totaldouble;
	my $tot_nodouble_per = $totalmaskedper - $totdblper;
	print DOUBLESTATS "#TOTAL\nAMOUNT\t\tPERCENTAGE\nmasked\tin-double\tmasked_minus_in-double\t%masked\t%in-double\t%masked_minus_%in-double\n$totalmasked\t$totaldouble\t$tot_nodouble\t$totalmaskedper\t$totdblper\t$tot_nodouble_per\n";
	
	
	#SUMMARY OVERLAP
	print DOUBLESTATS "\n\n#SUMMARY OF OVERLAPS\nMasked_by:\tAmount\t%\n";
	if (defined $MaskedSameC){
		my $SameCPer = $MaskedSameC / $GenLen * 100;
		print DOUBLESTATS "Same: class (all)\t$MaskedSameC\t$SameCPer\n";
	}
	if (defined $MaskedSameCnotF){
		my $SameCnotFPer = $MaskedSameCnotF / $GenLen * 100;
		print DOUBLESTATS "Same: class (but NOT same family)\t$MaskedSameCnotF\t$SameCnotFPer\n";
	}
	if (defined $MaskedSameF){
		my $SameFper = $MaskedSameF / $GenLen * 100;
		print DOUBLESTATS "Same: class/fam\t$MaskedSameF\t$SameFper\n";
	}
	if (defined $MaskedSameFnotE){
		my $SameFnotEper = $MaskedSameFnotE / $GenLen * 100;
		print DOUBLESTATS "Same: class/fam (but NOT same element)\t$MaskedSameFnotE\t$SameFnotEper\n";
	}
	if (defined $MaskedSameE){
		my $SameEper = $MaskedSameE / $GenLen * 100;
		print DOUBLESTATS "Same: name#class/fam\t$MaskedSameE\t$SameEper\n";
	}
	if (defined $MaskedDiff){
		my $Diffper = $MaskedDiff / $GenLen * 100;
		print DOUBLESTATS "Different: class\t$MaskedDiff\t$Diffper\n";
	}
	
	
	#BY CLASS
	print DOUBLESTATS "\n\n#DETAILS BY CLASS\n\tAMOUNT\t\t\tPERCENTAGE_OF_GENOME\nclass_name\tmasked\tin_double[same_class]\tin_double[same_class/family]\tmasked\tin_double[same_class]\tin_double[same_class/family]\n";
	foreach my $keyC (sort keys %class) {
		$doubleC{$keyC} = 0 unless ($doubleC{$keyC});
		$doubleCnotF{$keyC} = 0 unless ($doubleCnotF{$keyC});
		my $samefam = $doubleC{$keyC} - $doubleCnotF{$keyC};
		my $perCA = $class{$keyC} / $GenLen * 100;
		my $perCD = $doubleC{$keyC} / $GenLen * 100;
		my $perCDD = $samefam / $GenLen * 100;
		print DOUBLESTATS "$keyC\t$class{$keyC}\t$doubleC{$keyC}\t$samefam\t$perCA\t$perCD\t$perCDD\n";
	}
	
	#BY CLASS/FAM
	print DOUBLESTATS "\n\n#DETAILS BY CLASS/FAMILY\n\tAMOUNT\t\t\tPERCENTAGE_OF_GENOME\nclass/family_name\tmasked\tin_double[same_class/family]\tin_double[same_repeat]\tmasked\tin_double[same_class/family]\tin_double[same_repeat]\n\n";
	foreach my $keyF (sort keys %fam) {
		$doubleF{$keyF} = 0 unless ($doubleF{$keyF});
		$doubleFnotE{$keyF} = 0 unless ($doubleFnotE{$keyF});
		my $samerep = $doubleF{$keyF} - $doubleFnotE{$keyF};
		my $perFA = $fam{$keyF} / $GenLen * 100;
		my $perFD = $doubleF{$keyF} / $GenLen * 100;
		my $perFDD = $samerep / $GenLen * 100;
		print DOUBLESTATS "$keyF\t$fam{$keyF}\t$doubleF{$keyF}\t$samerep\t$perFA\t$perFD\t$perFDD\n";
	}

	
	###################################################
	# OVERLAP FILES
	###################################################
	# all
	open DOUBLE, ">$lengthDouble" or die print LOG "\t !! could not open $lengthDouble $!\n";
	print DOUBLE "\n#See below details about ALL overlaps, e.g. everytime the same piece of DNA is masked twice\n
#Scaffold/Chr:\tThis element:\t\t\t\t\tHas overlap with this upstream element:\t\t\t\t\tOVERLAP:
#Gname\tGstart\tGend\tRname\tRclass\tRfam\tGstart\tGend\tRname\tRclass\tRfam\tLength\n\n";
	foreach my $keyD (sort keys %double) {
		print DOUBLE "$keyD\t$double{$keyD}\n";
	}
	close DOUBLE;
	
	#details for different class (oter file)
	open DOUBLED, ">$lengthDoubleD" or die print LOG "\t !! could not open $lengthDoubleD $!\n";
	print DOUBLED "\n#See below details about overlapping by different TE of different class\n
#Scaffold/Chr:\tThis element:\t\t\t\t\tHas overlap with this upstream element:\t\t\t\t\tOVERLAP:
#Gname\tGstart\tGend\tRname\tRclass\tRfam\tGstart\tGend\tRname\tRclass\tRfam\tLength\n\n";
	foreach my $keyDD (sort keys %doubleD) {
		print DOUBLED "$keyDD\t$doubleD{$keyDD}\n";
	}
	close DOUBLED;
}
###################################################



print LOG "\n--- Parsing finished, see:\n\t$stats\n";

unless (defined $norw){
	print LOG "\n    Plus RM out rewritten in:\n\t$output\n";
}

if (defined $totaldouble){
	print LOG "\n    Plus some amounts/% per class and family in files:\n\t$DoubleStatFile\n\t$lengthDouble\n\t$lengthDoubleD\n----------------------------\n";
}
close LOG;
exit;




sub median { 
	@_ == 1 or die ('Sub usage: $median = median(\@array);'); 
	my ($array_ref) = @_; 
	my $count = scalar @$array_ref;
	# Sort a COPY of the array, leaving the original untouched 
	my @array = sort { $a <=> $b } @$array_ref; 
	if ($count % 2) { 
		return $array[int($count/2)]; 
	} else { 
		return ($array[$count/2] + $array[$count/2 - 1]) / 2; 
	} 
} 
