#!/usr/bin/perl -w

#######################################################
# Author  :  Aurelie K
# date    :  Jan 2013 = version 1.0
# version :  1.6
# email   :  4urelie.k@gmail.com
# Pupose  :  Parse RepeatMasker output
# 	 			-> get big tables with % divergence breakdown in bins to have dynamic landscape of a given family
#				   as on Repeat Masker website download pages
######################################################
# Updated : v1.1, 25 Feb 2013
#		- allow to skip Nrem step (ie if masked genome from UCSC or RM wensite)
# Updated : v1.2, 28 Feb 2013
#		- changes for TEclass stuff
#		- output with percentages grouped by class -> faster to create graph
# Updated : v1.3, Nov 2013
#		- changed TEclass format (to match input used in other scripts)
#		- updated way of getting arguments
#		- follow up of version
# Updated : v1.4, Feb 2014
#		- -TEage option => age outputs
#		- %div correcteed using -300/4*ln(1-%div*4/300)
#		- changes => subroutines
# Updated : v1.5, Feb 2014
#		- -filter option
#		- when filter, no output for class, useless
#		- 100 bins, too much, even when sometimes there are artifiicially high %div (rm blast does that) => stop at 50.
# Updated : v1.6, Feb 2014
#		- Option -recent => will provide outputs only up to 5% divergence, but bins of 0.25 (e.g. [0-0.25[, [0.25;0.50[ etc



### TO DO

#		- allow filtering on %div, more general that "recent", and more dynamic
#		- check for nonTE stuff and filter them out as an option


######################################################
use strict;
use warnings;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Getopt::Long;
my $version = "1.6";
my $usage = "Usage, v.1.5:
     perl <scriptname.pl> -genfile <GenomeFile> -RMout <MaskedGenome.out> [-TEage <TEage>] [-TEs <TEclass>] [-filter [type,word>] [-noNrem] [-recent]
     OR
     perl <scriptname.pl> -genlen LengthOfGenone(nt) -RMout <MaskedGenome.out> [-TEage <TEage>] [-TEs <TEclass>] [-recent]
	
	(note that order of these options doesn't matter)
	
	MANDATORY ARGUMENTS:
	 -RMout <MaskedGenome.out>  => <MaskedGenome.out> = file to analyse
	 -genfile OR -genlen:
	     -genfile <GenomeFile>      => <GenomeFile> = genome file (path of the file)
	     -genlen LengthOfGenone(nt) => provide length of genome instead of genome file
	
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
	 -filter = apply a filter (on class or family through \"type\" argument)
	           typically, to get output only on DNA TEs: -filter class,DNA
	                      to get output only for Alus elements: -filter family,ALU
	           Note that: - the filter is case insentitive (DNA is same as dna, LINE1 same as line1)
	                      - this might be super wrong to use unless you provide at least a TEclass file to make sure all DNA TEs have a class set as DNA for ex
	                      - this is mostly relevant if a TEage file is provided
	                      - exact match needed (you don't want to get nonLTR when LTR is chosen for ex)
	 -noNrem = only relevant if -genfile is used
			   chosing this means \"do not check for removal of Ns in genome / do not remove Ns\". If not chosen, Ns will be removed from genome before getting length.
	           (using genome with or without Ns will change % of genome, but obviously not raw amounts)
	 -recent = script will provide outputs only up to 5% divergence, but bins of 0.25 (e.g. [0-0.25[, [0.25;0.50[ etc
			   These 2 parameters (max %div printed + size of bins) can be changed in the subroutine get_boundaries_bins (just find it, toward the end of the script)
			   
	WHAT IT DOES: 
	 For each fragment (line of repeat masker output), amount of DNA masked is put in a bin, based on the % of divergence to consenus of this fragment.
	 This is done per repeat and per class, as well as per lineage is -TEage is used (see above)
	 This script generates outputs that are tabulated text files. 
	
	SUGGESTION:
	 If you have -aln files, correct the %div with CG correction (Kimura)
	 How: check you RM version
	    If > 4.0.5, should be already included in the -aln files... Is it in .out - not released yet??
	    If < 4.0.5, look in the RepeatMasker directory, Utils directory, called calcDivergenceFromAlign.pl and run it
	    			then run this script on the modified .out\n\n";

my $GenFile;
my $GenLen;
my $RMout;
my $TEclass;
my $TEage;
my $noNrem;
my $filter;
my $recent;
GetOptions ('noNrem'=> \$noNrem, 'recent'=> \$recent, 'genfile=s' => \$GenFile, 'genlen=s' => \$GenLen, 'RMout=s' => \$RMout, 'TEs=s' => \$TEclass, 'TEage=s' => \$TEage, 'filter=s' => \$filter);

#check that all options are provided
if (((! $GenFile) && (! $GenLen)) || (! $RMout) || (($GenLen) && ($GenLen =~ /\D/))){
	die $usage;
}

#check if -recent option
my $ifrecent;
($recent)?($ifrecent = 1):($ifrecent = 0);

#get filter, in lower case + check if OK
my ($filter_t, $filter_w) = split (",",lc($filter)) if ($filter);
if (($filter) && ((! $filter_t) || (! $filter_w) || (($filter_t ne "class") && ($filter_t ne "family")))){
	print "\t\t$filter was defined but filter type is neither class or family => filter can't be applied\n";
	die $usage;
}

my $log;
($filter)?($log = "$RMout.GetLandscape.$filter_w.log"):($log = "$RMout.GetLandscape.log");
open LOG, ">$log" or die "could not open $log $!\n";
print LOG "\n --- Script started, v$version:\n";
print LOG "     RM output file to analyze = $RMout\n";
print LOG "     About the genome:\n";
print LOG "          Option -noNrem chosen...\n" if ($noNrem);
print LOG "          Genome file defined = $GenFile\n" if ($GenFile);
print LOG "          Length of genome given = $GenLen\n" if ($GenLen);
(($TEage)?(print LOG "     TEclass file defined = $TEclass, but since TE age ($TEage) is defined it will be used instead to correct class of the elements\n"):(print LOG "     TEclass file (to correct class of the elements) defined = $TEclass\n\n")) if ($TEclass);
print LOG "     TEage file (for additional output, split on age) defined = $TEage\n" if ($TEage);
print LOG "     Filter $filter will be applied => keep only $filter_t that match $filter_w (case insentitive)\n" if ($filter);
print LOG "\n";
close LOG;

# GET LENGTH OF GENOME FROM FILE IF RELEVANT
unless ($GenLen) {
	my $genome;
	# 1) check if Ns should be removed / check for removal
	($noNrem)?($genome = $GenFile):($genome = Nrem($GenFile,$log));
	# 2) calculate length of the genome, to get percentages
	$GenLen = get_len($genome,$log);
}

# Parsing the RM output
# Get TE infos if provided
open LOG, ">>$log" or die "could not open $log $!\n";
print LOG " --- Getting TE info from defined file...\n" if (($TEage) || ($TEclass));
my $TE;
$TE = get_TEs_infos($TEage) if ($TEage);
$TE = get_TEs_infos($TEclass) if ((! $TEage) && ($TEclass));
#TE:
# 0       1       2         3        4       5			6
#Rname	Rclass	Rfam	Rclassfam	%div	~Age	Acient/LineageSpe


#####################################################
#Now core of script
#####################################################
print LOG " --- Parsing Repeatmasker output to load all the data...\n";

#prepare outputs; use $t to make output with or without corrected class and fam
my ($out_Len,$out_Per,$out_LenClass,$out_PerClass);
($filter)?($out_Len = "$RMout.Landscape.$filter_w.Len.all.tab"):($out_Len = "$RMout.Landscape.Len.all.tab");
($filter)?($out_Per = "$RMout.Landscape.$filter_w.Perc.all.tab"):($out_Per = "$RMout.Landscape.Perc.all.tab");
($filter)?($out_LenClass = ""):($out_LenClass = "$RMout.Landscape.Len.class.tab");
($filter)?($out_PerClass = ""):($out_PerClass = "$RMout.Landscape.Perc.class.tab");
my $t = 0;
$t = 1 if (($TEclass) || ($TEage));
prep_output($out_Len,$t,$ifrecent);
prep_output($out_Per,$t,$ifrecent);
prep_output($out_LenClass,2,$ifrecent) unless ($out_LenClass eq "");
prep_output($out_PerClass,2,$ifrecent) unless ($out_PerClass eq "");
my ($out_LenAge,$out_PerAge);
if ($TEage) {
	($filter)?($out_LenAge = "$RMout.Landscape.$filter_w.Len.age.tab"):($out_LenAge = "$RMout.Landscape.Len.age.tab");
	($filter)?($out_PerAge = "$RMout.Landscape.$filter_w.Perc.age.tab"):($out_PerAge = "$RMout.Landscape.Perc.age.tab");
}	
prep_output($out_LenAge,3,$ifrecent) if ($TEage);
prep_output($out_PerAge,3,$ifrecent) if ($TEage);


#parse
my %hash_of_elem = ();
my %hash_of_counts = ();
my %hash_of_class = ();
my %hash_of_classcounts = ();
my %hash_of_age = ();
my %hash_of_agecounts = ();

open RMOUT, "<$RMout" or die print LOG "ERROR: could not open $RMout $!\n"; 
LINE: while(<RMOUT>) {
	chomp (my $line = $_);	
	next LINE if (($line =~ /position|[sS]core|Gname/) || ($line !~ /\w/)); #skip headers and white lines
	$line =~ s/^\s+//; #remove spaces in beginning of lines
	my ($div,$Gstart,$Gend,$Rname,$Rclassfam,$Rclass,$Rfam,$Gmasked) = get_info_from_RMout_file($line);
	
	# correct the % of divergence
	$div = -300/4*log(1-$div*4/300); #note that in perl, log in ln; log10 
	
	#Storing data for each element (Fullname)
	my $Rfullname = $Rname."#".$Rclass."/".$Rfam;
	my @infos = ",$div,$Gmasked";
	
	#check if there was a TEclass or TEage file provided, if not fake the hash
	unless (($TEclass) || ($TEage)){
		my @TE = ($Rname,$Rclass,$Rfam,"$Rclass/$Rfam","na","undet");
		$TE->{lc($Rname)} = \@TE;
		#Rname	Rclass	Rfam	Rclassfam	%div	~Age	Acient/LineageSpe
	}
	
	#get class, fam + also data if relevant [if corrected class => will be used]; also note if $TE{lc($Rname) should exist but doesn't
	unless($TE->{lc($Rname)}) {
		my @TE = ($Rname,$Rclass,$Rfam,"$Rclass/$Rfam","na","undet");
		$TE->{lc($Rname)} = \@TE;
		print LOG "\t\t$Rname not found in $TEclass => RM output class / fam used\n" if (($TEclass) && (! $TEage));
		print LOG "\t\t$Rname not found in $TEage => RM output class / fam used, age = undet\n" if ($TEage);
	}
	
	#If there was a filter defined, process only lines that go through it
	#Apply filter on corrected Rclass and Rfam
	if ($filter) {
		if ($filter_t eq "class") {
			next LINE unless ($filter_w eq lc($TE->{lc($Rname)}->[1]));
		} elsif ($filter_t eq "family") {
			next LINE unless ($filter_w eq lc($TE->{lc($Rname)}->[2]));
		}
	}
	
	#load hashes now
	$hash_of_elem{$Rfullname} .= "@infos";
	($hash_of_counts{$Rfullname})?($hash_of_counts{$Rfullname}++):($hash_of_counts{$Rfullname} = 1);
	
	my $Rclasscorr = $TE->{lc($Rname)}->[1];
	unless ($filter) {
		$hash_of_class{$Rclasscorr} .= "@infos";
		($hash_of_classcounts{$Rclasscorr})?($hash_of_classcounts{$Rclasscorr}++):($hash_of_classcounts{$Rclasscorr} = 1);
	}	
	
	#deal with age split if relevant
	if ($TEage){
		my $age = $TE->{lc($Rname)}->[5];
		$hash_of_age{$age} .= "@infos";
		($hash_of_agecounts{$age})?($hash_of_agecounts{$age}++):($hash_of_agecounts{$age} = 1);
	}
}
close RMOUT;

##########################################################################################################
#Step 3. Read hashes and write in files [could still be "subroutined" but already better]
##########################################################################################################
print LOG " --- Repeats are being processed => write outputs...\n";

# ALL ELEMENTS
open (OUTLEN,">>$out_Len") or die print LOG "ERROR: could not open $out_Len $!\n";
open (OUTPER,">>$out_Per") or die print LOG "ERROR: could not open $out_Per $!\n";
foreach my $key (sort keys %hash_of_elem) {
	my ($len_by_div,$per_by_div) = get_bins($key,$hash_of_counts{$key},\%hash_of_elem,$GenLen,$ifrecent); #it gets references of tables @len_by_div
	my ($Rname,$Rclassfam)= split (/#/,$key);
	my ($Rclass,$Rfam)= split (/\//,$Rclassfam);
	my $firstcols = "$Rname\t$Rclass\t$Rfam";
	$firstcols = "$Rname\t$Rclass\t$Rfam\t$TE->{lc($Rname)}->[1]\t$TE->{lc($Rname)}->[2]" if (($TEclass) || ($TEage));
	print_bins($firstcols,$out_Len,$len_by_div,$ifrecent);
	print_bins($firstcols,$out_Per,$per_by_div,$ifrecent);
}

# CLASS
unless ($filter) {
	open (OUTLEN,">>$out_LenClass") or die print LOG "ERROR: could not open $out_LenClass $!\n";
	open (OUTPER,">>$out_PerClass") or die print LOG "ERROR: could not open $out_PerClass $!\n";
	foreach my $key (sort keys %hash_of_class) {
		my ($len_by_div,$per_by_div) = get_bins($key,$hash_of_classcounts{$key},\%hash_of_class,$GenLen,$ifrecent);
		print_bins($key,$out_LenClass,$len_by_div,$ifrecent);
		print_bins($key,$out_PerClass,$per_by_div,$ifrecent);
	}
}

# AGE
if ($TEage) {
	open (OUTLEN,">>$out_LenAge") or die print LOG "ERROR: could not open $out_LenAge $!\n";
	open (OUTPER,">>$out_PerAge") or die print LOG "ERROR: could not open $out_PerAge $!\n";
	foreach my $key (sort keys %hash_of_age) {
		my ($len_by_div,$per_by_div) = get_bins($key,$hash_of_agecounts{$key},\%hash_of_age,$GenLen,$ifrecent);
		print_bins($key,$out_LenAge,$len_by_div,$ifrecent);
		print_bins($key,$out_PerAge,$per_by_div,$ifrecent);
	}
}

print LOG " --- Parsing finished. \n     Outputs for all elements";
(($TEclass) || ($TEage))?(print LOG " (with class corrected using provided file):\n"):(print LOG ":\n");
print LOG "     \t $out_Len\n     \t $out_Per\n";
unless ($filter) {
	print LOG "     Outputs for all classes:";
	(($TEclass) || ($TEage))?(print LOG " (with class corrected using provided file):\n"):(print LOG "\n");
	print LOG "     \t $out_LenClass\n     \t $out_PerClass\n";
}	
print LOG "     Outputs for lineages:\n     \t $out_LenAge\n     \t $out_PerAge" if ($TEage);

print "\n --- Script GetLandscape done\n     Check the log file = $log\n     for steps + if some repeats not found in provided TE files (if relevant)\n\n";

exit;




##################################################################################################################################################################
# SUBROUTINES
##################################################################################################################################################################
######################################################
# QUITE GENERAL FASTA FILES TREATMENT
######################################################
#----------------------------------------------------------------------------
# remove Ns in a genome file
#----------------------------------------------------------------------------
#removing the N in the genome
sub Nrem {
	$GenFile = shift;
	$log = shift;
	my $genometmp = "$GenFile.Nrem";
	my $genome = "$GenFile.Nrem.fa";
	
	#open log file
	open(my $log_fh, ">>", $log) or print "ERROR: could not open $log $!\n";	
	if (-e $genome) {
		print $log_fh "\n --- N already removed from $GenFile ($genome exists) - skipping N removal step...\n";
	} else {
		#remove Ns
		print $log_fh "\n --- N may be not removed from $GenFile ($genome does not exists) - Removing N...\n";
		my $Gentmp = Bio::SeqIO->new(-file => $GenFile, -format => "fasta") or die "Failed to create Bio::SeqIO object from $GenFile $!\n";
		open(my $GminN_fh, ">", $genometmp) or die print $log "ERROR: could not create $genometmp $!\n";	
		while( my $seq = $Gentmp->next_seq() ) {
			my $sequence = $seq->seq;
			$sequence =~ s/[Nn]//g;
			my $id = $seq->display_id."\t".$seq->desc;
			print $GminN_fh ">$id\n$sequence\n";
		}
		# Rewrite sequences in fasta format just to be sure
		my $Gen = Bio::SeqIO->new(-file => $genometmp, -format => "fasta") or die print $log "ERROR: could not create Bio::SeqIO object from $genometmp $!\n";
		my $GenMinusN = Bio::SeqIO->new(-file => ">$genome", -format => "fasta") or die print $log "ERROR: could not create Bio::SEQIO outputfile $genome $!\n";
		while( my $seq2 = $Gen->next_seq() ) {
			$GenMinusN->write_seq($seq2);		
		}
	}
	return $genome;
}
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# calculate length of the genome to have the total length (to get percentages) - includes check steps to avoid repeating this if length already calculated
# => my $GenLen = get_len($genome,$log);
# 	 note that log file need to be closed in the main before calling the subroutine, and reopened after
#----------------------------------------------------------------------------
sub get_len {
	my $genome = shift;
	my $log = shift;

	#open log file
	open(my $log_fh, ">>", $log) or print "ERROR: could not open $log $!\n";	

	#calculating total length of the genome (of the database)
	my $lengthfile = "$genome.length";
	
	unless (-e $lengthfile) {
		print $log_fh " --- Total length of genome not known ($lengthfile does not exists) - Calculating total length...\n";
		# index the genome and connect to the fasta file
		my $reindex;
		my $indexfile = "$genome.index";
		if (-e $indexfile) {
			$reindex = 0;
			print $log_fh "     Genome previously indexed ($indexfile exists) - Skipping indexing step...\n";
		} else {
			$reindex = 1;
			print $log_fh "     Genome not indexed ($indexfile does not exists) - Indexing...\n";
		}
		my $db = Bio::DB::Fasta->new($genome, -reindex=>$reindex) or die print $log "     ERROR: could not create Bio::DB::Fasta object from $genome $!\n";
		
		#create list of the ID of the genome file
		my @dbIDs = $db->get_all_ids();
		
		open (my $len_fh, ">", $lengthfile) or die print $log "     ERROR: could not create file $lengthfile $!\n\n";
		$GenLen = 0;
		foreach my $ID (@dbIDs) {
			my $obj = $db->get_Seq_by_id($ID);
			my $len = $obj->length;
			$GenLen += $len;
		}
		print $len_fh $GenLen;
	} else {
		print $log_fh " --- Total length of genome has been previously calculated ($lengthfile exists)\n";
	}
	
	open (my $len_fh2, "<", $lengthfile) or die print $log_fh "ERROR: could not open $lengthfile $!\n";
	while (<$len_fh2>) {
		$GenLen = $_;
		}	
	print $log_fh "     => total length = $GenLen nt\n";
	
	return ($GenLen);
}
#----------------------------------------------------------------------------



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


#----------------------------------------------------------------------------
# get values from RM output line, while (<RMOUT>)(=> split in subroutine)
#----------------------------------------------------------------------------
sub get_info_from_RMout_file {
	my $line = shift;
	my @line = split(/\s+/,$line);
	#get infos needed
	my $div = $line[1];
	my $Gstart = $line[5];
	my $Gend = $line[6];
	my $Rname = $line[9];
	my $Rclassfam = $line[10];
	#getting other values
	my ($Rclass,$Rfam);
	if ($Rclassfam =~ /\//) {
		($Rclass,$Rfam) = split(/\//, $Rclassfam);
	} else {
		$Rfam = $Rclassfam;
		$Rfam=~ s/^(.*)\..*$/$1/;
		$Rclass = $Rclassfam;
		$Rclass =~ s/^.*\.(.*)$/$1/;
	}
	my $Gmasked = ($Gend - $Gstart + 1);
	return ($div,$Gstart,$Gend,$Rname,$Rclassfam,$Rclass,$Rfam,$Gmasked);
}
#----------------------------------------------------------------------------

######################################################
# QUITE SPECIFIC TO THIS SCRIPT
######################################################
#----------------------------------------------------------------------------
# to get the boundaries of the bins
# my ($max, $bin) = get_boundaries_bins($ifrecent);
#----------------------------------------------------------------------------
sub get_boundaries_bins {
	my $ifrecent = shift;
	my ($max,$bin);
	if ($ifrecent == 0) { #this means the option "ifrecent" is NOT chosen
		$max = 50; #e.g. bins won't go higher than 50% divergence
		$bin = 1; #e.g. bins will be of size 1 => [0;1[, [1;2[ etc
	} elsif ($ifrecent == 1){ #this means the option "ifrecent" is chosen
		$max = 5; #e.g. bins won't go higher than 50% divergence
		$bin = 0.25; #e.g. bins will be of size 0.25 => [0;0.25[, [0.25;0.5[ etc
	}	
	return($max,$bin);
}

#----------------------------------------------------------------------------
# prep outputs
#----------------------------------------------------------------------------
sub prep_output {
	my $out = shift; #file name
	my $t = shift; #0, 1, 2 or 3 here; det behavior
	my $ifrecent = shift; #0 or 1
	open(my $out_fh, ">", $out) or die "ERROR: could not create $out $!\n";
	print $out_fh "\t% of divergence bins:\nLineage" if ($t == 3);
	print $out_fh "\t% of divergence bins:\nClass" if ($t == 2);
	print $out_fh "\t\t\t\t\t% of divergence bins\nRname\tRclass\tRfam\tRclassCorrect\tRfamCorrect" if ($t == 1);
	print $out_fh "\t\t\t% of divergence bins\nRname\tRclass\tRfam" if ($t == 0);
	my ($max, $bin) = get_boundaries_bins($ifrecent);
	for (my $i = 0; $i < $max; $i+=$bin) {
		my $i2 = $i+$bin;
		print $out_fh "\t\[$i;$i2\[";
	}
	print $out_fh "\n\n";
}
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# get values from RM output line (simple array > split in subroutine)
# foreach my $key (sort keys %hash_of_elem) {
	# my ($len_by_div,$per_by_div) = get_bins($key,$hash_of_counts{$key},\%hash_of_elem,$GenLen,$ifrecent);
#----------------------------------------------------------------------------
sub get_bins {
	my ($key,$n,$refhash,$GenLen,$ifrecent) = @_; 
		#$n = number of elements
		#$refhash = info from RM
	#deal with size of bins
	my ($max, $bin) = get_boundaries_bins($ifrecent);
	#prep values
	my %len_by_div = ();
	my %per_by_div = ();
	$refhash->{$key} =~ s/\s/,/g;
	$refhash->{$key} =~ s/,//;
	my @infos = split(/,/,$refhash->{$key}); #=> RMoutput info, 2 values per element, => ",$div,$Gmasked" added
	#loop in values stored
	for (my $i = 0; $i < $n*2; $i=$i+2) {
		my $div = $infos[$i]; #%div
		my $len = $infos[$i+1]; #length masked
		my $per = $len/$GenLen*100;
		unless ($div > $max) {
			FINDBIN: for (my $j = $bin; $j <= $max; $j+=$bin) {
				my $coord = $j-$bin; 
				if (($div >= $coord) && ($div < $j)) {
					$len_by_div{$coord}+=$len;
					$per_by_div{$coord}+=$per;
					last FINDBIN;
				}
			}
		}	
	}
	return(\%len_by_div,\%per_by_div);
}
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# to print the bins
#----------------------------------------------------------------------------
sub print_bins {
	my ($key,$out,$bins,$ifrecent) = @_;
	open(my $out_fh, ">>", $out) or die print "ERROR: could not open $out $!\n";
	print $out_fh $key;
	my ($max, $bin) = get_boundaries_bins($ifrecent);
	for (my $k=0; $k<=$max; $k+=$bin) {
		($bins->{$k})?(print $out_fh "\t$bins->{$k}"):(print $out_fh "\t0");
	}
	print $out_fh "\n";
}
#----------------------------------------------------------------------------


