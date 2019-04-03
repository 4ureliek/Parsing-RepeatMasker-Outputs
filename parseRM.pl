#!/usr/bin/perl -w
#----------------------------------------------------------------------------
# Author  :  Aurelie Kapusta
# email   :  4urelie.k@gmail.com
# PURPOSE :  parsing Repeat Masker outputs, .out and/or .align
#----------------------------------------------------------------------------
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Bio::SeqIO;
#use Data::Dumper;

my $VERSION = "5.8.2";
my $CHANGELOG;
set_chlog();
sub set_chlog {
	$CHANGELOG = "
#	- v1.0 = 03 Apr 2015
#            Script that was only for parsing age
#	- v1.1 = 08 Apr 2015
#            added -align and changed way of counting on this script
#	- v1.2 = 09 Apr 2015	
#            Filter nonTE stuff
#	- v3.0 = 04 May 2015 (MayThe4th!)
#            Make this a global parseRM script => up version to match the fact that previous parseRM was 2.4
#            and changed name from parseRM_age.pl -> parseRM.pl
#	- v3.1 = 05 May 2015
#            counter for number of Gnames processed, 10 by 10, to see progress
#	- v3.2 = 08 May 2015
#            Modify counter + this is taking way to long, so change of the data structure to store %div
#	- v3.3 = 22 May 2015
#            Log correction
#	- v4.0 = 20 Jun 2016
#            Fixed the -parse option; changes in storing of %div, %del and %ins 
#               to remove redundant info and make it run faster
#            Count fragments that are start to end of consensus, +/- 2
#            Count fragments that are > 50% and > 90 % of the consensus (different than above)
#            Get median of fragments length for .out but removed for .align since it is
#               a fake median (not on all lengths but on all values for all positions)
#	- v4.1 = 09 Sep 2016
#            Fixed bug that would crash when -p is set without -f
#	- v5.0 = 04 Jul 2017
#            Clean up the usage, the code and the subs
#            Bug fix - last sequence was not processed
#            Bug fix - counts hash check for age was wrong (fullid / Rid)
#            Bug fix - counts age hash when print was without the 'a'...?
#            Bug fix - one sub was calling itself... It was probably the issue of -p non progressing
#	- v5.1 = 09 Jul 2017
#            Bug fix - hash per genome sequence were not emptied
#            Speed up the parsing - thanks to GitHub user WuChangCheng (BioWu), who used NYTProf to identify 
#               that the '\@div = sort \@div if (\$DIV[0])' line was super long, I changed way of doing things           
#	- v5.2 = 19 Jul 2017
#            Bug fix in the 'double_masking' & thus total_masked amounts (nr were OK)
#            Bug fix in loading the kimura corrected %div from .align
#            Check step in case no repeat to parse in a file (for example if only nonTE stuff)
#            Usage update
#	- v5.3 = 19 Jul 2017
#            Bug fix when -a used alone
#	- v5.4 = Sep 13 2017
#            convention of variables uc/lc
#            changes in the subs to use more references
#            address the speed issue -> need to still improve but it is already 2x faster on a 4000 lines .out file
#	- v5.5 = Sep 26 2017
#            bug fix for age parsing, when -d set and -a INT (and not a file)
#            bug fix genome files -> ls does not have full path
#            bug fix load and access repeat lengths from library file
#	- v5.6 = Oct 3-5 2017
#            bug fix for -d: 
#                in log, \"XXX Gname done\" counted for all when -d used => counters passed as local in sub
#                split age was hash not reinitialized => numbers added up => add sub clean_up_hashes
#            bug fix finding fa file when -d not set
#            bug fix for N removal
#            bug fix in printing all splitage in one file when -d set
#            bug fix landscape
#            bug fix length of genomes
#	- v5.7 = Oct 12 2017
#            Won't crash parsing a .align if some repeats don't have the # (that defines Rname#Rclassfam)
#               it prints warnings if -v is set, so user can check if it was unintentional (typo, etc)
#	- v5.8 = Feb 28 2018
#            Change for the total length calculation of the fasta file
#	- v5.8.1 = Jul 16 2018
#            Bug fix for -s all    
#	- v5.8.2 = Apr 03 2019
#            Cosmetic & usage edits   

# TO DO:
# dig into using intervals with a start and end in an array instead of position by position...?
#    
\n";
	return 1;
}

my $USAGE;
set_usage();
my $FULLUSAGE;
set_help();

#-----------------------------------------------------------------------------
#-------------------------- LOAD AND CHECK OPTIONS ---------------------------
#-----------------------------------------------------------------------------
my $NONTE = "no";
my ($IN,$P,$AA,$LAND);
my ($GFILE,$GLEN,$NREM,$LIB);
my ($FILTER,$CONT,$TES);
my ($MY,$DIR,$K);
my ($HELP,$V,$CHLOG);
GetOptions ('in=s'     => \$IN, 
            'dir'      => \$DIR, 
            'parse'    => \$P, 
            'age=s'    => \$AA, 
            'land=s'   => \$LAND, 
            'my=s'     => \$MY,
            'fa'       => \$GFILE, 
            'nrem'     => \$NREM, 
            'glen=s'   => \$GLEN, 
            'rlib=s'   => \$LIB, 
            'kim'      => \$K, 
            'simple=s' => \$NONTE, 
            'edit=s'   => \$TES, 
            'what=s'   => \$FILTER, 
            'contain'  => \$CONT, 
            'updates'  => \$CHLOG, 
            'help'     => \$HELP, 
            'v'        => \$V);

#check steps on the options
my $TEAGE; #flag if --age points to a file
check_opt();
#print some log if $V
print_log("1") if ($V);

#make STDERR buffer flush immediately
select((select(STDERR), $|=1)[0]); 

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
#----- Get list of input files if -dir, or just load the one from -in
my @FILES = ();
($DIR)?(@FILES = `ls $IN`):(push(@FILES,$IN));
my $PATH = get_path($IN);

#----- Prep --parse stuff if relevant
my $GENLEN = ();
my $LIBLEN = ();
prep_parse() if ($P);

#----- Prep --age stuff if relevant
my $AGE = ();
prep_age() if ($AA);

#----- Prep --land stuff if relevant
my ($MAX,$BIN) = split(",",$LAND) if ($LAND);

#----- Prep steps not specific to p, a or l
#Load substitution rates from $MY if relevant
print STDERR " --- Loading substitution rates info ($MY)\n" if ($MY && $V);
my $SRATES = ();
load_val($MY,"rates") if ($MY);
#Deal with $TES if relevant
my $TE = load_TE_info($TES) if ($TES);

#----- Now parse all the RM files
print STDERR " --- Now parsing RM output(s)\n" if ($V);
my $NB;
my $TOT = ();
$TOT->{'nr'}= 0;
$TOT->{'double'}= 0;
#main parsers, per Gname:
my %LOWDIV = ();
my %ID = ();
my %INDEL = ();
my %DOUBLE = ();
#processed data:
my $PARSED = ();
my $MASKED = ();
my $COUNTS = ();
my $NR_CHECK = ();
my $LANDSCAPE = ();
#Loading data:
my $F;
my $BIG = ();
my ($DIV,$DEL,$INS,$GNAME,$GST,$GEN,$STRAND,$RFULLID,$RID); #'columns' in $BIG
my ($RNAME,$CLASSFAM);
my %IFCLASS = (); #check if no class or fam
my ($BLOCK,$RFULLNAME);
my ($RCLASS,$RFAM,$RCLASSFAM);
my ($SKIPPED,$PREVSKIP);
my ($FNAME,$FANAME,$ALLREP);
FILE: foreach my $file (@FILES) {
	chomp $file;
	next FILE unless ($file); #ls is weird sometimes, blank values??
	$F = $IN."/".filename($file) if ($DIR);
	$F = $PATH."/".filename($file) if (! $DIR);
	$FNAME = filename($F);
	print STDERR "     -> $F..\n" if ($V);
	
	#Check if f can/should be parsed
	if (! -e $F) {
		print STDERR "        ..skipped, does not exist?\n" if ($V);
		#note: this looks stupid but that way if a file is deleted while this is running, won't die
		next FILE;
	}
	if ($F !~ /\.align$/ && $F !~ /\.out$/) {
		print STDERR "        ..skipped, not .out or .align?\n" if ($V);
		next FILE;	
	}
	
	#Now parse:
	#filter & load files in an array of array
	print STDERR "        ..loading in array..\n" if ($V);
	load_RM_in_array();
	print STDERR "          ..done\n" if ($V);

	print STDERR "          WARN: no repeat to parse in $file? Skipping\n" if (! $BIG && $V);
	next FILE if (! $BIG);
			
	#loop through the array to store all %div per nt piece so that nt can be split 
	#into age category and bins [not the best memory usage wise]
	print STDERR "        ..Looping through array (parsing)..\n" if ($V);		
	$NB = scalar(@{$BIG});	
	parseRM_table();
		
	#Now, print all data contained in hash tables
	print STDERR "        ..printing parsed data..\n" if ($V);	
	$TOT->{'double_na'} = 0; #value needed for now in print_split
	if ($P) {
		$ALLREP = "$F.parseRM.all-repeats.tab";
		print_parsed_summary();	
		print_parsed_allrep();
	}	
	if ($AA) {
		my $fall = $IN.".splitage_all.tab" if ($DIR);
		open my $fhall, ">", $fall or confess "\nERROR (main): could not open to write $fall $!\n" if ($DIR);
		print_age($fhall);
		close $fhall if ($DIR); 
	}	
	print_landscape() if ($LAND);
	#clean up all the hashes now - empty for each file
	clean_up_hashes();
	print STDERR "          .. done\n" if ($V);
}
#----- Done - log & exit
print_log("2") if ($V);
exit;

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

#--------------------------------- MAIN SUBS ---------------------------------
#----------------------------------------------------------------------------
sub load_RM_in_array {
	my $i = 0;
	$PREVSKIP = "";
	$SKIPPED = 0; #for RMout (if no block)
	open (my $fh, "<", $F) or confess "\nERROR (sub get_RMalign_array): can't open to read $F $!\n";
	LINE: while(defined(my $l = <$fh>)) {	
		chomp $l;
		my $next = check_line(\$l);
		next LINE if ($next eq "y");
	
		#Load these values:
		if ($F=~ /\.align/) {
			get_RMalign_val(\$l,\$i);
		} elsif ($F =~ /\.out/) {
			get_RMout_val(\$l);	
			next LINE if (! $BLOCK);
		}
		
		#get, + correct if needed, the Rclass and Rfam
		($RCLASS,$RFAM,$RCLASSFAM) = get_Rclass_Rfam();
		
		#filter stuff (if relevant)
		my $skip = "no";
		$skip = get_RM_array_skip() unless ($NONTE eq "all" && ! $FILTER);
		$PREVSKIP = "no";
		$PREVSKIP = $skip if ($skip);
		next LINE if ($skip eq "yes");
		
		#Reconstruct large array
		if ($F=~ /\.align/) {
			$RFULLID = $RNAME."-#-".$GNAME.":".$GST."_".$GEN."_".$STRAND;
		} else {
			#For .out, use blocks:
			$RFULLID = $RNAME."-#-".$BLOCK.":".$GNAME;
		}	
		$RID = $RNAME."-#-".$GNAME.":".$GST."-".$GEN;
		
		#Now load
		if (substr($l,0,1) =~ /[0-9]/) {
			my @new = ($DIV,$DEL,$INS,$GNAME,$GST,$GEN,$STRAND,$RFULLID,$RID);
			push(@{$BIG},\@new);
			$i++;
		}	
	}		
	close $fh;
	print STDERR "          WARN: $SKIPPED lines skipped because they had no block info\n" if ($V && $SKIPPED > 0);
	return 1;
}		

#----------------------------------------------------------------------------
	sub get_RMalign_val {
	my $l = shift;
	my $i = shift;
	#Load values:
	if (substr($$l,0,1) =~ /[0-9]/) {
		#kindda the RM.out line, without block info
		#773 12.31 10.00 0.00 Scaffold10260 147 276 (1002) Cannrnd-3_family-5#LTR/ERVL 440 582 (2) 1 m_b1349s001i0
		#453 15.99 4.90 1.90 Scaffold10260 1174 1275 (3) C Cannrnd-5_family-2606#LTR/ERV1 (0) 580 476 2 m_b1349s001i1
		#=> apparently there is only a strand column when it is C. Interesting.
		my @l = split(/\s+/,$$l);
		my $sc;
		($sc,$DIV,$DEL,$INS,$GNAME,$GST,$GEN) = ($l[0],$l[1],$l[2],$l[3],$l[4],$l[5],$l[6]);
		if (($l[8] eq "C") || ($l[8] eq "+")) { #just in case it gets changed in a later version
			($STRAND,$RFULLNAME) = ($l[8],$l[9]);
		} else {
			($STRAND,$RFULLNAME) = ("+",$l[8]);
		}
		if ($RFULLNAME =~ /#/) {
			($RNAME,$CLASSFAM)=split("#",$RFULLNAME);
		} else {
			$RNAME = $RFULLNAME;
			$CLASSFAM = "Unknown";
			print STDERR "          WARN: $RNAME had no class or family defined -> set as Unknown\n" unless (defined $IFCLASS{$RNAME} && $V);
			$IFCLASS{$RNAME} = 1;
		}	
	#check the %div:	
	} elsif (substr($l,0,6) eq "Kimura" && $PREVSKIP eq "no") {	
		# (NB: missing for many - seems to be simple repeats...??)
		my ($whatever,$Kimdiv) = split("=",$$l);
		$Kimdiv =~ s/\s//g;
		$BIG->[$$i-1]->[0]=$Kimdiv;
	}
	return 1;
}

#----------------------------------------------------------------------------
	sub get_RMout_val {
	my $l = shift;
	$$l =~ s/^\s+//; #remove spaces in beginning of lines
	my @l = split(/\s+/,$$l);			
	my $sc;
	($sc,$DIV,$DEL,$INS) = ($l[0],$l[1],$l[2],$l[3]);
	($GNAME,$GST,$GEN) = ($l[4],$l[5],$l[6]);
	($STRAND,$RNAME,$CLASSFAM,$BLOCK) = ($l[8],$l[9],$l[10],$l[14]);
	$SKIPPED++ unless ($BLOCK);						
	# correct the % of divergence if $K
	$DIV = -300/4*log(1-$DIV*4/300) if ($K);
}

#-----------------------------------------------------------------------------
sub parseRM_table {
	my $len = 0;
	my ($nbscaff,$nbscafftmp) = (0,0);
	LINE: for (my $i = 0; $i < $NB; $i++){
		($DIV,$DEL,$INS) = ($BIG->[$i][0],$BIG->[$i][1],$BIG->[$i][2]);
		($GNAME,$GST,$GEN) = ($BIG->[$i][3],$BIG->[$i][4],$BIG->[$i][5]);
		($STRAND,$RFULLID,$RID) = ($BIG->[$i][6],$BIG->[$i][7],$BIG->[$i][8]);			
		
		#Process data
		#First, see if %div should be converted in My using substitution rate provided and just replace the value
		if ($MY) {
			my $subs = $SRATES->{filename($F)};
			my $myears = $DIV / 100 / ($subs * 2);
			$DIV = $myears;
		}
		#Now store lowest %div(or My) per position as well as the associated repeat info
		for (my $n = $GST; $n <=$GEN; $n++) {
			$TOT->{'tot'}++;
			$DOUBLE{$n}++ if ($LOWDIV{$n});
			$TOT->{'double'}++ if ($LOWDIV{$n});
			if (! $LOWDIV{$n} || ($LOWDIV{$n} && $LOWDIV{$n} > $DIV)) {
				$LOWDIV{$n} = $DIV;
				$ID{$n}{$DIV}{'f'}=$RFULLID;
				$ID{$n}{$DIV}{'r'}=$RID;
				if ($P) {
					$INDEL{$n}{$DIV}[0]=$DEL;
					$INDEL{$n}{$DIV}[1]=$INS;
				}	
			}		
		}	
					
		#Now parse all this if this is the last line of a Gname, or of file
		if ($i == $NB-1 || $GNAME ne $BIG->[$i+1][3]) {
			my $TOTdb = keys %DOUBLE;
			foreach my $n (sort keys %LOWDIV) {
				my $lowest = $LOWDIV{$n};
				$RFULLID = $ID{$n}{$lowest}{'f'}; #get corresponding repeat
				my $rest;
				($RNAME,$rest) = split("-#-",$RFULLID);
				$RID = $ID{$n}{$lowest}{'r'}; #get corresponding repeat
				my $lcname = lc($RNAME);
				($RCLASS,$RFAM) = ($TE->{$lcname}[1],$TE->{$lcname}[2]);
				parse_all_parse(\$n,\$lowest) if ($P);			
				my $type = "na";
				$type = parse_all_age(\$n,\$lowest) if ($AA);
				parse_all_land(\$type,\$lowest) if ($LAND);
			}			
			$nbscaff++;
			$nbscafftmp++;
			#empty these hash for each Gname
			%LOWDIV = ();
			%ID = (); 
			%INDEL = ();
			%DOUBLE = ();
			if ($nbscafftmp == 100) {
				print STDERR "          ..$nbscaff Gname done\n" if ($V);
				$nbscafftmp = 0;
			}
		}
	}
	print STDERR "          ..done\n" if ($V);
	return 1;		
}



#---------------------------------- GENERAL ---------------------------------
#----------------------------------------------------------------------------
sub load_val {
	my $data = shift;
	my $type = shift; #My or div
	if (-e $data) { #then it's a file
		open (my $fh, "<", $data) or confess "\nERROR (sub load_val): could not open to read $data $!\n"; 
 		while (defined(my $l = <$fh>)) {
 			chomp $l;
 			my ($name,$d) = split(/\s+/,$l);
			if ($type eq "div") {
				load_val_hash(\$d,\$name);
			} else {
				$SRATES->{$name}=$d if ($type eq "rates");
				$GENLEN->{$name}=$d if ($type eq "len");
			}
 		}
 		close ($fh);
	} else { #not a file
		my $name = filename($IN);
		if ($type eq "div") {
			load_val_hash(\$data,\$name);
		} else {
			$SRATES->{$name}=$data if ($type eq "rates");
			$GENLEN->{$name}=$data if ($type eq "len");
		}
	}
	return 1;
}

#----------------------------------------------------------------------------
sub load_val_hash {	
	my $d = shift;
	my $name = shift;
	my ($m,$n);
	($$d =~ /,/)?(($m,$n) = split(",",$$d)):(($m,$n)=($$d,$$d));
	my @v = (); #values
	($m < $n)?(@v=($m,$n)):(@v=($n,$m));
	$AGE->{$$name}=\@v;
	return 1;
}	

#----------------------------------------------------------------------------
sub get_Rclass_Rfam {
	my ($Rc,$Rf);
	if ($CLASSFAM =~ /\//) {
		($Rc,$Rf) = split(/\//, $CLASSFAM);
	} else {
		$Rf = $CLASSFAM;
		$Rf=~ s/^(.*)\..*$/$1/;
		$Rc = $CLASSFAM;
		$Rc =~ s/^.*\.(.*)$/$1/;
	}
	my $Rcf = "$Rc/$Rf";
	
	#Correct class anf fam if relevant
	my $lcRname = lc($RNAME);
	($Rc,$Rf) = ($AGE->{$lcRname}[1],$AGE->{$lcRname}[2]) if ($TEAGE && $AGE->{$F}{$lcRname}[1]);
	if ($TE->{$lcRname}[1]) {
		($Rc,$Rf) = ($TE->{$lcRname}[1],$TE->{$lcRname}[2]);
	} else {
		($TE->{$lcRname}[1],$TE->{$lcRname}[2]) = ($Rc,$Rf);
	}
	
	return ($Rc,$Rf,$Rcf);
}

#----------------------------------------------------------------------------
sub check_line {
	my $l = shift;
	if ($F=~ /\.align/) {
		return "y" if ($$l !~ /\w/ || substr($$l,0,1) =~ /\s/ || substr($$l,0,2) eq "C ");
		return "y" if (substr($$l,0,6) eq "Matrix" || substr($$l,0,11) eq "Transitions" || substr($$l,0,8) eq "Gap_init");	
	} else {	
		return "y" if (($$l =~ /position|[sS]core|Gname/) || ($$l !~ /\w/));
	}
	return "n";
}

#----------------------------------------------------------------------------
sub get_RM_array_skip {
	my $skip = "no";
	#filter out non TE stuff unless asked not
	unless (($NONTE eq "nonTE") || ($NONTE eq "all")) {
		$skip = "yes" if ($RCLASS eq "nonTE");
	}
	unless ($NONTE eq "all") {
		$skip = "yes" if (($RCLASS eq "Simple_repeat") 
		               || ($RCLASS eq "Low_complexity") 
					   || ($RCLASS eq "Satellite") 
					   || ($RCLASS =~ /RNA$/) 
					   || ($RCLASS =~ /omeric$/) 
					   || ($RCLASS eq "ARTEFACT"));
	}
	#filter out stuff if relevant		
	if ($FILTER) {
		my ($F_type,$F_name) = split(",",$FILTER);
		my ($lcRname,$lcRclass,$lcRfam,$lcf_name) = (lc($RNAME),lc($RCLASS),lc($RFAM),lc($F_name));
		if ($CONT) {
			#check if what's det as the filter is included in the names
			$skip = "yes" unless ((($F_type eq "name") && ($lcRname =~ /$lcf_name/))
						      || (($F_type eq "class") && ($lcRclass =~ /$lcf_name/)) 
						      || (($F_type eq "family") && ($lcRfam =~ /$lcf_name/)));
		} else {
			$skip = "yes" unless ((($F_type eq "name") && ($lcf_name eq $lcRname))
						      || (($F_type eq "class") && ($lcf_name eq $lcRclass)) 
						      || (($F_type eq "family") && ($lcf_name eq $lcRfam)));
		}	
	}
	return ($skip);
}

#----------------------------------------------------------------------------
sub get_path {
	my($file) = shift;
	($file =~ /\//)?($file =~ s/(.*)\/.*$/$1/):($file = ".");
	return $file;
}

#----------------------------------------------------------------------------
sub filename {
	my($file) = shift;
	$file =~ s/.*\/(.*)$/$1/;
	return $file;
}

#----------------------------------------------------------------------------
sub clean_up_hashes {
	$TOT = ();
	$TOT->{'nr'}= 0;
	$TOT->{'double'}= 0;
	#main parsers, per Gname:
	%LOWDIV = ();
	%ID = ();
	%INDEL = ();
	%DOUBLE = ();
	#processed data:
	$PARSED = ();
	$MASKED = ();
	$COUNTS = ();
	$NR_CHECK = ();
	$LANDSCAPE = ();
	#Loading data:
	$BIG = ();
	return 1;
}



#----------------------------- RELATED TO PARSE -----------------------------
#----------------------------------------------------------------------------
sub prep_parse {
	print STDERR " --- Prepping steps for --parse...\n" if ($V);
	#genome length
	if ($GLEN) {
		print STDERR "     - Loading sequence length(s) from $GLEN...\n" if ($V);
		load_val($GLEN,"len");
	} elsif ($GFILE) {
		print STDERR "     - Getting sequence length(s) from genome file(s)...\n" if ($V);
		get_tot_length(); 
	}
	#lib lengths
	if ($LIB) {
		print STDERR "     - Getting lengths of consensus sequences from $LIB...\n" if ($V);
		get_all_lengths();
	}	
	return 1;
}

#----------------------------------------------------------------------------
sub get_tot_length {
	my @list;	
	if (! $DIR) {
		my $fa = $IN;
		$fa = $1 if ($IN =~ /(.*)\.out/ || $IN =~ /(.*)\.align/);
		$fa = $fa.".fa" unless ($fa =~ /.fa(sta|$)/);
		@list = `find $PATH/$fa`;
		if (! $list[0]) {
			print STDERR "       WARN: $fa could not be found, % genome won't be determined\n";
			return 1;
		}
	} else {
		@list = `find $IN/*fa`;
	}
	#Now get length & print it	
	FA: foreach my $fa (@list) {
		chomp $fa;
		next FA if ($fa eq $IN);
		next FA if ($fa !~ /.fa(sta|$)/);
		next FA unless (substr(`head -n 1 $fa`,0,1) eq ">");
		print STDERR "        -> $fa\n" if ($V);	
		my $flen = "$fa.length";		
		print STDERR "           obtaining total length\n" if ($V);
		my $len = 0;
		if ($NREM) {
			$len = `cat $fa | grep -v ">" | awk '{ total += length(\$0) }; END { print total }'`;
		} else {
			$len = `cat $fa | grep -v ">" | sed 's/n\|N//g' | awk '{ total += length(\$0) }; END { print total }'`;
		}
		chomp $len;
		$GENLEN->{filename($fa)} = $len;
		print STDERR "           => ".$GENLEN->{filename($fa)}." nt\n" if ($V);

	}	
	return 1;
}

#----------------------------------------------------------------------------
sub get_all_lengths {
	#Get lengths of all sequences and store that by file and sequence ID. 
	#Note that if some are not unique, it just replaces by last length.
	#Specifically for TEs.
	print STDERR "        -> $LIB\n" if ($V);
	unless (substr(`head -n 1 $LIB`,0,1) eq ">") {
		print STDERR "       WARN: $LIB not a fasta file? => skipped\n" if ($V);
		return 1;
	}	
	if (! -e $LIB) {
		print STDERR "       WARN: $LIB could not be found, % genome won't be determined\n";
		return 1;
	}
	my $lfile = "$LIB.lengths";
	my $do = 0;
	if (-e $lfile) {
		print STDERR "           lengths have been previously calculated ($lfile exists) => extracting\n" if ($V);
		#extract lengths now
		open (my $lfh, "<", $lfile) 
		     or warn "           WARN: could not open to read $lfile, but lengths will be recalculated ($!)\n" && $do++;
		unless ($do == 1) {
		while (defined(my $l = <$lfh>)) {
				chomp $l;
				my ($ID,$ln) = split(/\s+/,$l);
				my ($RNAME,$RCLASSFAM) = split('#',$ID);
				$LIBLEN->{lc($RNAME)}=$ln;
			}	
			close ($lfh);
			return 1;
		}
	}
	if ((! -e $lfile) || ($do == 1)) {
		print STDERR "           obtaining lengths\n" if ($V);
		my $id = "";
		my $ln = 0;
		my $c = 0;
		open (my $fa_fh, "<", $LIB) 
		     or warn "           WARN: could not open to read $LIB, repeat lengths won't be in the output $!\n";
		open (my $len_fh, ">", $lfile) 
		     or warn "           WARN: could not open to write $lfile, but lengths will be calculated $!\n";
		while (defined(my $l = <$fa_fh>)) {
			chomp $l;
			if (substr($l,0,1) eq ">") {
				#first get and print unless first header
				unless ($c == 0) {
					print $len_fh "$id\t$ln\n";
					my ($Rn,$Rc) = split('#',$id);
					$LIBLEN->{lc($Rn)}=$ln;
				}
				$c=1;
				#store header and reinitialize length
				my @id = split (/\s+/,$l);
				$id = $id[0];
				$id =~ s/>//;
				$ln = 0;
			} else {
				#get length; could be more than one line so increment
				$ln+=length($l);
			}
		}
		#get and print len last sequence
		print $len_fh "$id\t$ln\n";
		my ($Rn,$Rcf) = split('#',$id);
		$LIBLEN->{lc($Rn)}=$ln;
		close ($fa_fh);
		close ($len_fh);
		print STDERR "           lengths are now extracted in $lfile\n" if ($V);
	}	
	return 1;
}

#----------------------------------------------------------------------------
sub Nrem {
	my $fa = shift;
	my $genome = $fa;
	$genome = $fa.".Nrem.fa" if ($fa !~ /\.Nrem\.fa/);
	if (-e $genome) {
		print STDERR "           Ns already removed from $fa ($genome exists) - skipping N removal step\n";
	} else {
		print STDERR "           Ns may be not removed from $fa ($genome does not exists) - Removing N...\n";
		my $faobj = Bio::SeqIO->new(-file => $fa, -format => "fasta") or die "Failed to create Bio::SeqIO object from $fa $!\n";
		open(my $fh, ">", $genome) or confess "\nERROR (sub Nrem): could not open to write $genome $!\n";	
		while( my $seq = $faobj->next_seq() ) {
			my $sequence = $seq->seq;
			$sequence =~ s/[Nn]//g;
			my $id = $seq->display_id."\t".$seq->desc;
			print $fh ">$id\n$sequence\n";
		}
		close $fh;
	}
	return $genome;
}

#----------------------------------------------------------------------------				
sub parse_all_parse {
	my $n = shift;
	my $ldiv = shift;	

	my $dbl = $DOUBLE{$$n};	
	
	#Store values, basically increment for each base => nr
	$MASKED->{'pn'}{$RNAME}{'nr'}++;
	$MASKED->{'pn'}{$RNAME}{'double'}++ if ($dbl);
	
	$MASKED->{'pc'}{$RCLASS}{'nr'}++;
	$MASKED->{'pc'}{$RCLASS}{'double'}++ if ($dbl);
	
	$MASKED->{'pf'}{$RFAM}{'nr'}++;
	$MASKED->{'pf'}{$RFAM}{'double'}++ if ($dbl);

	#count fragments using the full ID; with .align files it will be a bit off mostly for LINEs
	if (! $NR_CHECK->{'pnr'}{$RFULLID}) {
		$COUNTS->{'pn'}{$RNAME}{'nr'}++;
		$COUNTS->{'pc'}{$RCLASS}{'nr'}++;
		$COUNTS->{'pf'}{$RFAM}{'nr'}++;
		$NR_CHECK->{'pnr'}{$RFULLID}=1;
	}
	if (! $NR_CHECK->{'ptot'}{$RID}) {
		$COUNTS->{'pn'}{$RNAME}{'tot'}++;
		$COUNTS->{'pc'}{$RCLASS}{'tot'}++;
		$COUNTS->{'pf'}{$RFAM}{'tot'}++;
		$NR_CHECK->{'ptot'}{$RID}=1;
	}

	#make the list of %div/My, %del, %ins to be able to do median for each repeat, and average can be done at the same time.
	#Increment parsed hashes:
	$PARSED->{$RNAME}{'div'}.=$$ldiv.",";
	$PARSED->{$RNAME}{'del'}.=$INDEL{$$n}{$$ldiv}[0].",";
	$PARSED->{$RNAME}{'ins'}.=$INDEL{$$n}{$$ldiv}[1].",";
	return 1;	
}

#----------------------------------------------------------------------------
sub print_parsed_summary {		
	my $out = "$F.parseRM.summary.tab";
	open(my $fh, ">", $out) or confess "\nERROR (sub print_parsed_summary): could not open to write $out $!\n";
	print $fh  "#This file gives some summary info about the masking, total, by class and by family\n";
	print $fh  "#(all repeats are detailed in the file $ALLREP)\n";
	print $fh  "#Overlap or double corresponds to DNA fragments that are masked by several elements.\n";
	print $fh  "#These amounts need to be subtracted in order to get more accurate TE amount.\n";
	print $fh  "#If a .align file was parsed, these amounts will be much higher than for the associated .out\n";
	print $fh  "#Note that overlaps may not be parsed correctly if names are not formatted consistently (Rname#Rclass/Rfam)\n";	
	
	$FANAME = $1 if ($FNAME =~ /(.*)\.out/ || $FNAME =~ /(.*)\.align/);
	$FANAME = $FANAME.".fa" unless ($FANAME =~ /.fa(sta|$)/);
	
	$GENLEN->{$FANAME} = "nd" unless ($GENLEN->{$FANAME});
	my $pertotnr = "nd";
	$TOT->{'nr'} = $TOT->{'tot'} - $TOT->{'double'};
	$pertotnr = $TOT->{'nr'} / $GENLEN->{$FANAME} * 100 if ($GENLEN->{$FANAME} ne "nd");

	print $fh "\n#TOTAL:\n";
	print $fh "#nt_total_in_genome\tnt_masked-minus-double\t%_masked\tFYI:nt_masked_double\n";
	print $fh "$GENLEN->{$FANAME}\t$TOT->{'nr'}\t$pertotnr\t$TOT->{'double'}\n";	
		
	print $fh "\n#BY CLASS\n";
	print $fh "#class\tnt_masked-minus-double\t%_masked\tFYI:nt_masked_double\n";	
	print_parsed_summary_details('pc',$fh);
	
	print $fh "\n#BY FAMILY\n";
	print $fh "#family\tnt_masked-minus-double\t%_masked\tFYI:nt_masked_double\n";
	print_parsed_summary_details('pf',$fh);
	close $fh;
	print STDERR "          -> $out\n" if ($V);
	return 1;
}	

#----------------------------------------------------------------------------
sub print_parsed_summary_details {	
	my $type = shift;
	my $fh = shift;
	foreach my $key (keys %{$MASKED->{$type}}) {
		my $nt = $MASKED->{$type}{$key}{'nr'};
		my $per = "nd";
		$per = $nt / $GENLEN->{$FANAME} * 100 if ($GENLEN->{$FANAME} ne "nd");
		my $db = 0;
		$db = $MASKED->{$type}{$key}{'double'} if ($MASKED->{$type}{$key}{'double'});
		my $TOT = $nt + $db;
		print $fh "$key\t$nt\t$per\t$db\n";	
	}
}	

#----------------------------------------------------------------------------
sub print_parsed_allrep {			
	open(my $fh, ">", $ALLREP) or confess "\nERROR (sub print_parsed_allrep): could not open to write $ALLREP $!\n";
	print $fh "#Rname\tRclass\tRfam\tRlen\tFRG_NB_all";
#	print $fh "\tFRG_NB_StartToEnd";
	print $fh "\tFRG_NB_Reconstructed_repeats";
	print $fh "\tLEN_MASKED_NR\tAVG_%DIV\tMED_%DIV\tAVG_%DEL\tMED_%DEL\tAVG_%INS\tMED_%INS\tAVG_LEN_MASKED\t%_GENOME\n\n";
	foreach my $name (keys %{$MASKED->{'pn'}}) {
		my $lc = lc($name);
		my $rlen = "nd";
		$rlen = $LIBLEN->{$lc} if (defined $LIBLEN->{$lc});
		my ($dela,$delm,$diva,$divm,$insa,$insm) = get_avg_med_from_list(\$name);
		my $len = $MASKED->{'pn'}{$name}{'nr'};
		my $avglen = $len / $TOT->{'nr'};
		my $perlen = "nd";
		$perlen = $len / $GENLEN->{$FANAME} if ($GENLEN->{$FANAME} ne "nd");
		print $fh "$name\t$TE->{$lc}[1]\t$TE->{$lc}[2]\t$rlen\t$COUNTS->{'pn'}{$name}{'tot'}";
#		print $fh "\tnd";
		print $fh "\t$COUNTS->{'pn'}{$name}{'nr'}\t$len\t$$diva\t$$divm\t$$dela\t$$delm\t$$insa\t$$insm\t$avglen\t$perlen\n";
	}
	close $fh;
	print STDERR "          -> $ALLREP\n" if ($V);
	return 1;
}

#----------------------------------------------------------------------------
sub get_avg_med_from_list {
	my $name = shift;
	my @med = ();
	@med = ('nd','nd','nd') if ($F =~ /\.align/);
	my @avg = ();
	foreach my $t (sort keys %{$PARSED->{$$name}}) { #del div ins in that order since sorted
		my @list = split (',', $PARSED->{$$name}{$t});
		push(@med,median(\@list)) unless ($F =~ /\.align/);
		push(@avg,average(\@list));
	}	
	return(\$avg[0],\$med[0],\$avg[1],\$med[1],\$avg[2],\$med[2]);
}

#----------------------------------------------------------------------------
sub median {
	my ($array_ref) = @_; 
	my $count = scalar @$array_ref;
	my @array = sort { $a <=> $b } @$array_ref; 
	if ($count % 2) { 
		return $array[int($count/2)]; 
	} else { 
		return ($array[$count/2] + $array[$count/2 - 1]) / 2; 
	} 
} 

#----------------------------------------------------------------------------
sub average { 
	my ($array_ref) = @_; 
	my $count = scalar @$array_ref;
	my $total = 0;
	foreach my $val (@$array_ref) {
		$total+=$val;
	}
	return($total/$count);
} 



#----------------------------- RELATED TO LAND ------------------------------
#----------------------------------------------------------------------------
sub parse_all_land {
	my $type = shift;
	my $div = shift;
	unless ($$div > $MAX) {
		FINDBIN: for (my $j = $BIN; $j <= $MAX; $j+=$BIN) {
			my $coord = $j-$BIN; 
			if ($$div >= $coord && $$div < $j) {
				$LANDSCAPE->{"Rname"}{"$RNAME\t$RCLASS\t$RFAM"}{$coord}++; #since it's per position here, simple increment
				$LANDSCAPE->{"Rclass"}{$RCLASS}{$coord}++; 
				$LANDSCAPE->{"Rfam"}{"$RCLASS\t$RFAM"}{$coord}++; 
				$LANDSCAPE->{"Rage"}{$type}{$coord}++ if ($AA && $TEAGE);
				last FINDBIN;
			}
		}	
	}	
	return 1;
}

#----------------------------------------------------------------------------
sub print_landscape {
	my $n = "Div";
	$n = "My" if ($MY);
	foreach my $type (keys %{$LANDSCAPE}) {
		my $out = $F.".landscape.$n.$type.tab";
		open (my $fh, ">", $out) or confess "\nERROR (sub print_landscape): could not open to write $out $!\n";
		prep_landscape_out($n,$type,$fh);			
		foreach my $key (keys %{$LANDSCAPE->{$type}}) {
			print $fh "$key";
			for (my $K=0; $K<$MAX; $K+=$BIN) {
				if ($LANDSCAPE->{$type}{$key}{$K}) {
					print $fh "\t$LANDSCAPE->{$type}{$key}{$K}";
				} else {
					print $fh "\t0";
				}	
			}
			print $fh "\n";	
		}
		close $fh;
		print STDERR "          -> $out\n" if ($V);
	}	
	return 1;
}

#----------------------------------------------------------------------------
sub prep_landscape_out {
	my $n = shift;
	my $type = shift;
	my $fh = shift;
	my $header = "";

	#prep
	$n = "Million years" if ($n eq "My");
	$n = "% of divergence" if ($n eq "Div");
	$header = "\t\t\t$n bins:\nRname\tRclass\tRfam" if ($type eq "Rname");
	$header = "\t\t$n bins:\nRclass\tRfam" if ($type eq "Rfam");
	$header = "\t$n bins:\nRclass" if ($type eq "Rclass");
	$header = "\t$n bins:\nLineage" if ($type eq "Rage");
	
	#Now print
	print $fh "$header";
	for (my $i = 0; $i < $MAX; $i+=$BIN) {
		my $i2 = $i+$BIN;
		print $fh "\t\[$i;$i2\[";
	}
	print $fh "\n\n";
	return 1;
}



#------------------------------ RELATED TO AGE ------------------------------
#----------------------------------------------------------------------------
sub prep_age {
	if ($TEAGE) {
		print STDERR "\n --- Loading age info from $AA\n" if ($V);
		load_TE_info($AA);
	} else {
		print STDERR "\n --- Loading age value ($AA)\n" if ($V);
		load_val($AA,"div");
	}	
	return 1;
}

#----------------------------------------------------------------------------
sub load_TE_info {
	my $in = shift;
	open(my $fh, "<", $in) or confess "\nERROR (sub get_TEs_infos): could not open to read $in $!\n"; 
	LINE: while(defined(my $l = <$fh>)) {
		chomp $l;
		next LINE if ($l !~ /\w/);
		my @TEs = split(/\t/, $l); 
		my $lcRname = lc($TEs[0]);
		$TE->{$lcRname} = \@TEs;
	}	
	close ($fh);
	return 1;
}

#----------------------------------------------------------------------------				
sub det_age_type {
	my $lowdiv = shift;	
	my $type = "na";
	if ($TEAGE) { #then it was a file, so use $rname
		my $lcRname = lc($RNAME);
		$type = $AGE->{$lcRname}[3];
	} else {
		my $name = $FNAME;
		$name = filename($IN) if ($DIR); #several files, but one value
		my ($min,$max)=($AGE->{$name}[0],$AGE->{$name}[1]);
		$type = "LS" if ($$lowdiv <= $min);
		$type = "A" if ($$lowdiv >= $max);
		$type = "nd" if (($$lowdiv < $max) && ($$lowdiv > $min));
	}
	return \$type;
}

#----------------------------------------------------------------------------				
sub parse_all_age {
	my $n = shift;
	my $lowdiv = shift;
	my $dbl = $DOUBLE{$$n};
	
	my $type = det_age_type($lowdiv);
	#Now store values, basically increment for each base => nr
	$MASKED->{'a'}{$$type}{'nr'}++;
	$MASKED->{'a'}{$$type}{'double'}++ if ($dbl);	

	#count fragments using the full ID; not really possible to do nr when .align, but check counting full IDs
	if (! $NR_CHECK->{'a'}{$RFULLID}) {
		$COUNTS->{'a'}{$$type}{'nr'}++;
		$NR_CHECK->{'a'}{$RFULLID}=1;
	}
	#get total counts with repeat ID; note that with .align files it will be a bit off			
	if (! $NR_CHECK->{'atot'}{$RID}) {
		$COUNTS->{'a'}{$$type}{'tot'}++;
		$NR_CHECK->{'atot'}{$RID}=1;
	}
	return $type;
}

#----------------------------------------------------------------------------
sub print_age {
	my $fhall = shift;
	my $out = $F.".splitage.tab";
	print STDERR "          -> $out\n" if ($V);	

	#Prep outputs if directory	
	open my $fh, ">", $out or confess "\nERROR (sub print_age): could not open to write $out $!\n";
	prep_age_out_headers($fh);
	prep_age_out_headers($fhall) if ($DIR);

	foreach my $type (keys %{$MASKED->{'a'}}) {
		$TOT->{'nr'} = $TOT->{'tot'} - $TOT->{'double'};
		my $nr_per = $MASKED->{'a'}{$type}{'nr'}/$TOT->{'nr'}*100;
		my $nr_c = "na";
		$nr_c = $COUNTS->{'a'}{$type}{'nr'} if ($COUNTS->{'a'}{$type}{'nr'});
		my $TOT_c = "na";
		$TOT_c = $COUNTS->{'a'}{$type}{'tot'} if ($COUNTS->{'a'}{$type}{'tot'});		
		print $fh    "$F\t$TOT->{'nr'}\t$TOT->{'tot'}\t$TOT->{'double'}\t$type";
		print $fhall "$F\t$TOT->{'nr'}\t$TOT->{'tot'}\t$TOT->{'double'}\t$type" if ($DIR);
		print $fh    "\t$TOT_c\t$nr_c\t$MASKED->{'a'}{$type}{'nr'}\t$nr_per\n";		
		print $fhall "\t$TOT_c\t$nr_c\t$MASKED->{'a'}{$type}{'nr'}\t$nr_per\n" if ($DIR);
	}	
	close $fh;
	return 1;
}

#----------------------------------------------------------------------------
sub prep_age_out_headers {
	my $fh = shift;
	print $fh "#nr_masked = amount of masked nt to consider\n";
	print $fh "#tot_masked = total amount of masked nt, including redundant maskins / overlaps\n";
	print $fh "#double_masked = total amount of nt that were masked by at least 2 repeats = redundant masking / overlaps\n\n";	
	print $fh "#Input_file\tnt_masked-minus-double\ttot_masked\tnt_masked_double\tAgeCat";
	print $fh "\tCounts_this_age\tnr_Counts_this_age\tnr_masked_this_age\t%nr_masked_this_age\n\n";	
	return 1;
}
	


#---------------------------- LOG, USAGE & HELP -----------------------------
#----------------------------------------------------------------------------
sub check_opt {
	die "\n Script parseRM.pl version $VERSION\n\n" if (! $IN && ! $P && ! $AA && ! $LAND && ! $HELP && ! $CHLOG && $V);
	die $CHANGELOG if ($CHLOG);
	die $FULLUSAGE if ($HELP);
	die $USAGE if (! $IN && ! $P && ! $AA && ! $LAND);
	die "\n please chose one of --parse, --age or --land (use -h to see the full usage)\n\n" if (! $P && ! $AA && ! $LAND);
	die "\n $IN does not exist?\n\n" if ($IN && ! -e $IN);
	die "\n $LIB does not exist?\n\n" if ($LIB && ! -e $LIB);
	check_file(\$GLEN) if ($GLEN);
	check_file(\$AA) if ($AA);
	$TEAGE = 1 if ($AA && -e $AA); #if it's a file, flag
	check_file(\$MY) if ($MY);	
	#avoid / at the end of path in case it's a directory
	$IN = $1 if ($DIR && $IN =~ /^(.*)\/$/);	
	return 1;
}

#----------------------------------------------------------------------------
sub check_file {
	my $file = shift;
	die "\n $file does not exist as a file but also not numerical?\n\n" if (! -e $$file && $$file !~ /[0-9\.,]+/);
	return:
}

#----------------------------------------------------------------------------
sub print_log {	
	my $type = shift;
	if ($type == 1) {
		print STDERR "\n--------------------------------------------------\n";
		print STDERR " --- Script parseRM.pl started (v$VERSION), with:\n";
		print STDERR "      - Directory containing input files = $IN (-d chosen)\n" if ($DIR);
		print STDERR "      - Input file = $IN\n" if (! $DIR);	
		print STDERR "      - All non TE repeats will be filtered out\n" if ($NONTE eq "no");
		print STDERR "      - Elements with class = nonTE won't be filtered out (--simple nonTE chosen)\n" if ($NONTE eq "nonTE");
		print STDERR "      - Non TE repeats won't be filtered out (--simple all chosen)\n" if ($NONTE eq "all");
		print STDERR "      - Repeats will be filtered based on $FILTER (--what)\n" if ($FILTER);
		print STDERR "        regexp will be used and not exact match (--contain)\n" if ($FILTER && $CONT);
		print STDERR "      - %div will be corrected using -300/4*log(1-%div*4/300)\n" if ($K);
		print STDERR "      => global parsing will be done, with following options:\n" if ($P);
		print STDERR "         - genome file(s) will be looked for (--fa chosen)\n" if ($GFILE);
		print STDERR "         - Ns will be removed from the genome before getting its length (--nrem chosen)\n" if ($NREM);
		print STDERR "         - length(s) are given instead (--glen $GLEN)\n" if ($GLEN);
		print STDERR "         - library file = $LIB\n" if ($LIB);
		print STDERR "         - no options\n" if ($P && ! $GFILE && ! $NREM && ! $GLEN && ! $LIB);
		print STDERR "      => age split will be done, with following options:\n" if ($AA);
		print STDERR "         - Age info will be based on $AA\n" if ($AA);
		print STDERR "         - Age will be in My, using substitution rates ($MY)\n" if ($AA && $MY);
		print STDERR "      => landscape data will be generated, with following options:\n" if ($LAND);
		print STDERR "         - With max_bin,bin_len = $LAND\n" if ($LAND);
		print STDERR "         - Age will be in My, using substitution rates ($MY)\n" if ($LAND && $MY);
		print STDERR "--------------------------------------------------\n";
	} else {
		print STDERR " --- Script done\n";
		if ($DIR) {
			print STDERR "    -> age split files: *.agesplit.tab for all files\n" if ($AA);
			print STDERR "                        + all in $IN.splitage_all.tab\n" if ($AA);
			print STDERR "    -> parsed files: *.parseRM.all-repeats.tab\n" if ($P);
			print STDERR "                     *.parseRM.summary.tab\n" if ($P);
			print STDERR "    -> landscape files: *.landscape.*.Rname.tab\n" if ($LAND);
			print STDERR "                        *.landscape.*.Rfam.tab\n" if ($LAND);
			print STDERR "                        *.landscape.*.Rclass.tab\n" if ($LAND);
			print STDERR "                        *.landscape.*.Rage.tab\n" if ($LAND && $TEAGE && $AA);
		}
		print STDERR "--------------------------------------------------\n\n";
	}	
	return 1;
}

#----------------------------------------------------------------------------
sub set_usage {
$USAGE = "
   PLEASE CITE: Kapusta, Suh & Feschotte (2017) PNAS (doi: 10.1073/pnas.1616702114)  
   And link to parseRM.pl, vX.X, https://github.com/4ureliek/Parsing-RepeatMasker-Outputs 
 
   Usage [v$VERSION]:   
   perl parseRM.pl -i <genome.(out|align)>
            [-p] [-f] [-n] [-g <X>] [-r <repeat_library.fa>]
            [-a <file.txt> OR <X,X>] [-m <file.txt> OR <X>] [-l <max,bin>] 
            [-d] [-k] [-s <type>] [-e <file>] [-w <type,name>] [-c] 
            [-v] [-u] [-h]
	
   ==> To print full help and details of all options as well as more examples, 
   type: parseRM.pl -h
   	
   This script will process RepeatMasker outputs .out or .align file(s), 
   with 3 non exclusive behaviors that can all be set together.
   If all 3 options -a, -p and -l are set, there will be an additional output 
   with bins by %div or My, but by age categories (specified in -a) .
      
   Examples below are for one input file, but note that -m option can also load a file,
   in case of several input files (load files in a directory with -d). 
   For more examples, type: parseRM.pl -h

   FOR SUMMARY BY REPEAT:
   Set -p to get a summary of the masking, as well as amount or DNA, 
   counts of fragments, etc, for each repeat name (all-repeats.tab file), 
   family, class and total amount (summary file) - if sequence names were 
   formatted as per RepeatMasker nomenclature: Rname#Rclass/Rfamily, or
   at least Rname#Rclass (Rfamily will equal Rclass if no / ).
   Typically:
      perl parseRM.pl -i <MySpecies.align> -p -v
   Or, with all options:
      perl parseRM.pl -i <MySpecies.align> -p -f MySpecies.fa -n -r repeat_library.fa
   Providing the genome (file or total length) and the repeat library file 
   will add columns in the output, but they are not really necessary. 
   The -n option is to remove Ns before calculating % of the genome. 
	
   EVOLUTIONARY LANDSCAPE:
   Use -l to set behavior to split the amount of DNA by bins of %div or My, 
   allowing to generate landscape graphs for each repeat name, family or class.
   For examples - if input file is named MySpecies.align:
      To get the numbers in bins of 1% of divergence to consensus, up to 50%:
         perl parseRM.pl -i MySpecies.align -l 50,1 -v
      To get the numbers in bins of 0.25% of divergence to consensus, up to 5%:
         perl parseRM.pl -i MySpecies.align -l 5,0.25 -v
      To get the numbers in bins of 1M, up to 50, and with a substitution rate of 0.0021:
         perl parseRM.pl -i MySpecies.align -l 50,1 -m 0.0021 -v
         	
   AMOUNTS SPLIT BY AGE:
   Use -a to determine the amounts of DNA in a genome that is masked by repeats 
   of different lineages / %divergence categories.
   For examples - if the input file is named MySpecies.align:
      To split at 10% divergence to consensus:
         perl parseRM.pl -i MySpecies.align -a 10 -v
      To split at 25 My, if the the substitution rate is 0.0021:
         perl parseRM.pl -i MySpecies.align -a 25 -m 0.0021 -v

   ==> To print full help and details of all options as well as more examples, 
   type: parseRM.pl -h\n\n";      
	return 1;
}

#----------------------------------------------------------------------------
sub set_help {
    $FULLUSAGE = "
   PLEASE CITE: Kapusta, Suh & Feschotte (2017) PNAS (doi: 10.1073/pnas.1616702114)  
   And link to parseRM.pl, vX.X, https://github.com/4ureliek/Parsing-RepeatMasker-Outputs 
 
   Usage [v$VERSION]:   
   perl parseRM.pl -i <genome.(out|align)>
            [-p] [-f] [-n] [-g <X>] [-r <repeat_library.fa>]
            [-a <file.txt> OR <X,X>] [-m <file.txt> OR <X>] [-l <max,bin>] 
            [-d] [-k] [-s <type>] [-e <file>] [-w <type,name>] [-c] 
            [-v] [-u] [-h]
		
   This script will process RepeatMasker outputs .out or .align file(s), 
   with 3 non exclusive behaviors that can all be set together.
   If all 3 options -a, -p and -l are set, there will be an additional output 
   with bins by %div or My, but by age categories (specified in -a).
   
   See the see end of this help for examples of command lines.

   PARSING:
   Use -p to get a summary of the masking, as well as amount or DNA, 
   counts of fragments, etc, for each repeat name (all-repeats file), 
   family, class and total amount (summary file), if sequence names were 
   formatted as per RepeatMasker nomenclature: Rname#Rclass/Rfamily, or
   at least Rname#Rclass (Rfamily will equal Rclass if no / ).
	
   LANDSCAPE:
   Use -l <max,bin> to set behavior to split the amount of DNA by bins of %div or My, 
   allowing to generate landscape graphs for each repeat name, family or class.
         	
   AGE SPLIT:
   Set -a to set behavior to determine the amounts of DNA in a genome that is
   masked by repeats of different lineages / %divergence categories.
    
   MANDATORY ARGUMENTS:	
     -i,--in (STRING) 
         RepeatMasker output .out or RepeatMasker output .align 
         (use -d and provide a directory path here to parse several files)
         Use of .align is best for -l, because the %div corrected 
         for higher mutation rate at CpG sites and the Kimura 2-Parameter 
         divergence metric is in the .align files (thus, the graphs are smoother)
         However, the script will treat separately repeats that overlap 
         (merged in the .out, with one name kept): this means that for -p, 
         the amount of DNA masked by each of the repeats will be much higher 
         than if a .out file is parsed.

   OPTIONAL ARGUMENTS RELATED TO --parse
     -p,--parse (BOOL)
         To get a summary file with number of counts, average length, 
         number of reconstructed repeats, total length annotated as each repeat, etc
         The options below are only relevant if this is chosen.
         Note that for a .align, the median is not really median
         since it is done by position, which is why it is skipped
     -f,--fa (BOOL)
         If the file THAT WAS MASKED, typically a genome, is in 
         the same directory as the .out or .align file(s). 
         If not provided, the % of genome masked won't be calculated.
         Names should correspond before .align and .fa(sta), and/or .out and .fa(sta),
         for ex: '-i hg38.out -f' will look for hg38.fa
         and 'hg38.out_RM405_RB20140131 -f' will also look for hg38.fa
     -g,--glen (INT or STRING)
         Alternatively to genome file(s), you can provide the total length of the genome (in nt)
         If several genomes are being looked at (-d chosen for example) this can be a file,
         containing 2 columns: filename \\t X                           
            filename = name of the file to be parsed
            X = the value (in nt)
     -n,--nrem (BOOL)
         To remove Ns from the genome file before getting its length 
         (to calculate the percentages on the amount of DNA that is not Ns)
     -r,--rlib (STRING)
         To add the length of the consensus sequence included in the output,
         set here the library of consensus sequences used to mask the genome 
         If several maskings are parsed, you can concatenate the different libraries
         (same if RM was run with -lib, concatenate the RM library and the user one)   
         If a repeat could not be found, or if -r is not chosen, the value 
         in the output will be \"nd\". 

   OPTIONAL ARGUMENTS RELATED TO --land
     -l,--land (STRING)
         To generate additional outputs that can be used to make 
         landscape graphs by repeat, by family and by class (3 files).
         Two values should be given, separated by a comma: <max>,<bin>
            <max> is the value of the last bin, where everything 
                  higher than this %div or My will be placed
            <bin> is the size of the bins
         If -m is set, then numbers here HAVE TO BE in My instead of %divergence
         Typically: -l 50,1 (if %div are corrected for CpGs then values tend to be higher)
         See the end of this help for more examples.
     -m,--my (INT or STRING)
         To output bins in My and not %div: substitution rates need to be provided:
         if -d not chosen, use -m substitution_rate (ex. -m 0.00199)
         if -d is chosen, then different values can be provided for each input files,
            if you use -m subst-rates.tab, and subst-rates.tab is a file with 2 columns: 
            filename_as_in_dir \\t substitution_rate
         
   OPTIONAL ARGUMENTS RELATED TO --age 
     -a,--age (STRING)
         Option1: load a file
            File should contain TE age data, tabular delimited. Columns should be:
            Rname \\t Rclass \\t Rfam \\t Lineage
            The Rname needs to be an exact match to the repeat names in the RM output; 
            class and family can be changed
            If -d is chosen, the file must contain info for ALL repeats of all files
         Option2: 1 or 2 values (INT or FLOAT). Use a comma to sepearate 2 values;
            their order does not matter (-a 13.5,15 and -a 15,13.5 are equivalent)
            Any DNA masked with %div < the smallest number given (here 13.5) will be 
            placed as lineage specific. Any DNA masked with %div > the largest number 
            given (here 15) will be placed as ancient, and the rest will be placed as \"nd\". 
            If -m is set, then numbers here HAVE TO BE in My instead of %divergence
            If -d is chosen and different values should be used for each input files
               you can provide a file containing 2 columns: filename \\t X                           
               filename = name of the file to be parsed
               X = the value (in nt)
     -m,--my (INT or STRING)
            Same as above for --land
                              
   OTHER OPTIONAL ARGUMENTS
     -d,--dir (BOOL)
         Chose this if -i points to a directory containing several .out or .align files                       
         Note that the directory should contain no other file besides previous 
         outputs of this script (symbolic links created with ln -s are OK). 
            If genomes are provided, they need to end by .fa and names should be 
            identical before .align and .fa(sta), and/or .out and .fa(sta)
            No technical objection to have a mixture of .out and .align in the folder 
            (how to read them is based on having .out or .align in their names)
     -k,--k (BOOL)
         If no .align file available, use this to correct %div with 
         the Jukes Cantor equation: %div = -300/4*log(1-%div*4/300)
         Do not chose this if you are parsing any .align files, as the 
         %div are already corrected, with the kimura equation (CpG corrected)
     -s,--simple (STRING)
         If you wish to KEEP (some) nonTE sequences (as far as annotation goes of course)
         By default, all non TE sequences are filtered out.
         Chose between all or nonTE
            all   = KEEP ALL nonTE stuff (Low_complexity, Simple_repeat, Satellites, etc)
            nonTE = KEEP only non TE sequences when class = nonTE 
            (other stuff will be filtered out)                   
     -e,--edit (STRING)
         File with TE information as follow: Rname \\t Rclass \\t Rfam
         This can be used to correct some class and/or families of some repeats.
         It may complement a -a file. Not all repeats need to be in it, 
         just the ones to correct. Rname needs to match the repeat masker output; 
         an easy way to get a list that can be edited is to run the parseRM_simple.pl script
         and use the first columns of the output *.all-repeats.tab
         However, this can be proven difficult when parsing .align files, 
         because of the 3' and 5' separation of LINE1 elements (XXX_3end or XXX_5end)                   
     -w,--what (STRING)
         Typically: -w <type,name> to run the script on only a subset of repeats. 
         Not case sensitive.
         The type can be: name, class or family 
         and it will done on EXACT MATCH unless -c is chosen as well
         ex: -w name,nhAT1_ML => only fragments corresponding to the repeat 
                                 named exactly nhAT1_ML will be looked at
             -w class,DNA     => all repeats with class named exactly DNA, 
                                 as in ...#DNA/hAT or ...#DNA/Tc1
             -w family,hAT    => all repeats with family named exactly hAT,
                                 so NOT ...#DNA/hAT-Charlie for example)
     -c,--contain (BOOL) 
         To check if the \"name\" determined with -w is included in the 
         value of the RM output, instead of an exact match
         ex: name,HERVK => all fragments containing HERVK in their name
             family,hAT => all repeats with family containing hAT, such as 
                           ...#DNA/hAT, ...#DNA/hAT-Charlie, etc
     -v,--v (BOOL)
         Verbose mode, make the script talks to you
         Print the version if only option
     -u,--updates (BOOL)
         Print change log (updates)
     -h,--help (BOOL)
         Print this usage
         

   EXAMPLES:
     ==> EXAMPLES OF --land USAGE, if input file is named MySpecies.align:
         To get the numbers in bins of 1% of divergence to consensus, up to 50%:
            perl parseRM.pl -i MySpecies.align -l 50,1 -v
         To get the numbers in bins of 0.22% of divergence to consensus, up to 5%:
            perl parseRM.pl -i MySpecies.align -l 5,0.25 -v
         To get the numbers in bins of 1M, up to 50, and with a substitution rate of 0.0021:
            perl parseRM.pl -i MySpecies.align -l 1,50 -m 0.0021 -v
	     Alternatively, to get the numbers in bins of 1M, up to 50, and
	     MyStudy is a directory containg MySpecies1.align & MySpecies1.align:
			perl parseRM.pl -i MyStudy -d -l 50,1 -m subst-rates.tab -v
			where subst-rates.tab is a tabulated file that contains the substitution rates:
				MySpecies1.align\t0.0025
				MySpecies2.align\t0.0041      

     ==> EXAMPLES OF --age USAGE, if input file is named MySpecies.align:
         To split at 10% divergence to consensus:
            perl parseRM.pl -i MySpecies.align -a 10 -v
         To split at 25 My, if the the substitution rate is 0.0021:
            perl parseRM.pl -i MySpecies.align -a 25 -m 0.0021 -v
	     Alternatively, to split on 25My, and MyStudy is a directory 
	     containg MySpecies1.align & MySpecies1.align:
			perl parseRM.pl -i MyStudy -d -a 25 -m file.txt -v
			where file.txt contains the substitution rates:
				MySpecies1.align\t0.0025
				MySpecies2.align\t0.0041
			note that MyStudy may contain a mixture of .out and .align files,
			but that to split on age, .align files are advised.

     ==> EXAMPLES OF --parse USAGE, if input file is named MySpecies.align:
         Typically, simple:
            perl parseRM.pl -i MySpecies.align -p -v
         Or, with all options:
            perl parseRM.pl -i MySpecies.align -p -f MySpecies.fa -n -r repeat_library.fa
         If several files in a directory named MyStudy:
         	perl parseRM.pl -i MyStudy -d -p -v
         	where files in MyStudy may as follow (check usage of -d):
         	   MySpecies1.out.align
         	   MySpecies1.out.fasta
         	   MySpecies2.RM.out
         	   MySpecies2.RM.fa
\n";      
	return 1;
}
