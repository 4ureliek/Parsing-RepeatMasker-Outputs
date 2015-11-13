#!/usr/bin/perl -w
#######################################################
# Author  :  Aurelie Kapusta
# version :  see below / see changelog
# email   :  4urelie.k@gmail.com  
# github  :  https://github.com/4ureliek?tab=repositories
#######################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Bio::SeqIO; #required if -lib, comment if issues and don't use -lib
use Array::Transpose::Ragged qw/transpose_ragged/; #required if -big, comment if issues and don't use -big
use Statistics::R; #required if -R, comment if issues and don't use -R

my $version = "3.0";
my $scriptname = "TE-analysis_Coverage.pl";
my $changelog = "
#	- v1.0 = March 2012
#		Script written with Xiaoyu Zhuo
#	- v2.0 = May 2014
#		Re wrote the script with subroutines
#		-lib or most common length is kept for consensus length
#		Added several options (filter, min_len, min_frg, R etc)
#   - v2.1 = May 2014
#		R plots through Statistics::R module
#	- v3.0 = 09-13 Nov 2015
#		Input file can now also be the TEjoin file from the TE-analysis-pipeline script, v4+
#		Input file can now also be the TrInfo file from the TE-analysis-pipeline script, v4+
#       -update of the -filter option (add -contain)
#       -cat option
#       -cons, -skip options
#       Some other changes to updtae the code (more subroutines, mostly)
#       No intermediate files

# TO DO: differenciate sense and antisense when relevant
# TO DO: read .align files...?
\n";

my $usage = "\nUsage [$version]: 
	perl $scriptname -in <file.tab> -type <X>
                    [-RM <X>] [-lib <repeats.fa>] [-big] [-TEs <TEclass>] [-cons] [-force <X>] [-R] [-Rfile]
	                [-filter <type,filter>] [-contain] [-skip <X>] [-min_frg <X>] [-min_len <X>] [-cat <X>] 
	                [-v] [-h] [-chlog]
	
	This script will output the coverage of each repeat consensus in the input file
	Use -filter to to restrict this to some repeats
	This can be in the genome (RM.out) or after intersection with some features (TE-analysis_pipeline: https://github.com/4ureliek/TEanalysis)
	Output files can be directly used to plot the coverage in R (all repeats in one file with -big), use -Rfile to get command lines
	
	/!\\ Because of indels, convert coordinates from genome to consensus is not super precise. 
    This is why even features with 1 nt (TSS, etc) may span more than 1 nt.
    You can use -force Rstart or -force Rend to convert both coordinates from the start or the end in consensus respectively

    MANDATORY ARGUMENT:	
     -in      =>   (STRING) input file, see below for the 3 possibilities:
     -type    =>   (STRING) 3 options:
                               -type RMout   if the input file is a repeat masker .out output file
                               -type TEjoin  if the input file is a TEjoin.bed file (from TE-analysis_pipeline, v4+)
                               -type TrInfo if the input file is a exons.TEjoin.TrInfos.tab file (from TE-analysis_pipeline, v4+)
                                             /!\\ you also need to provide the original repeat masker .out file used for the analysis using -RM
     
    OPTIONAL ARGUMENTS 
     -RM      => (STRING)   If -type TrInfo, then use this to provide the repeat masker output file .out (needed for coordinates in consensus)
     -lib     => (STRING)   Provide the library file (fasta) to obtain consensus lengths (otherwise calculated from the .out info)
                            Requries Bio::SeqIO
     -big     =>   (BOOL)   add an output: a big tabulated file with all repeats in it, in addition to a file per repeat
                            Requires Array::Transpose::Ragged
     -TEs     => (STRING)   Optional input file to correct class and/or family of (some) TEs
                            with TE information as follow, 3 first columns mandatory: 
                            Rname \\t Rclass \\t Rfam \\t Rclass/Rfam \\t etc (only first 3 columns will be used)
     -cons    =>   (BOOL)   Print an output file that has the consensus lengths for all repeats
                            This can be important because of the merging of repeats between .align and .out of Repeat Masker
                            (coordinates in consensus could be wrong; use this option to check how many instances,
                            and if it could have a great impact skip the inconsistencies with -skip in filtering options)
     -force   =>   (STRING) Not relevant for -type RMout.
                            use -force Rstart to convert both coordinates from the start in consensus 
                            use -force Rend to convert both coordinates from the end in consensus
                            Note this can create issues (values outside of the consensus length)
     -R       =>   (BOOL)   To directly get the images of the plots (pdf)
                            Requires Statistics::R
     -Rfile   =>   (BOOL)   To print a file with R command lines to plot the coverage graphs
                            Behavior will be different if bigtable or not

    OPTIONAL FILTERING
     -filter  => (STRING)   To run the script on only a subset of repeats. Not case sensitive.
                            The type can be: name, class or family and it will be EXACT MATCH unless -contain is chosen as well
                              ex: name,nhAT1_ML => only fragments corresponding to the repeat named exactly nhAT1_ML will be looked at
                                   class,DNA => all repeats with class named exactly DNA (as in ...#DNA/hAT or ...#DNA/Tc1)
                                   family,hAT => all repeats with family named exactly hAT (so NOT ...#DNA/hAT-Charlie for example)
     -contain =>   (BOOL)   To check if the \"name\" determined with -filter is included in the value in Repeat Masker output, instead of exact match
                            Note that presence of [ or ] in names will create errors
                               ex: name,HERVK => all fragments containing HERVK in their name
                                   family,hAT => all repeats with family containing hAT (...#DNA/hAT, ...#DNA/hAT-Charlie, etc)
     -skip    =>    (INT)   To skip all instances of consensus length issues when > Xnt, that can arise from the merging of repeats 
                            between .align and .out of Repeat Masker (coordinates in consensus could be wrong)
                            Typically, -skip 1
     -min_frg =>    (INT)   Filter for a minimum number of fragments (e.g. coverage not so useful if less than 3 fragments for example)
                            Typically, -min_frg 3   
     -min_len =>    (INT)   Filter for a minimum length of the repeat consensus (in nt)
                            Can be useful if lots of unclasified short repeats
                            Typically, -min_len 80         
     -cat     => (STRING)   Relevant only if -TrInfo. Will filter on the category column to plot some features.
                            If there are several transcripts with the same feature it won't be counted several times.
                            5 possibilities: 
                               -cat exonized [default] to plot all exonized pieces of the TE. 
                                              -> This includes the \"exonized\" category as well as all partial overlaps from other categories
                               -cat TSS      to plot all TSS located in the TE (i.e. for the 3 categories: TSS, TSS_5SPL, TSS_polyA)
                               -cat polyA    to plot all polyA sites located in the TE (i.e. for the 3 categories: polyA, 3SPL_polyA, TSS_polyA)
                               -cat 5SPL     to plot all 5' SPL sites located in the TE (i.e. for the 3 categories: 5SPL, 3SPL_exon_5SPL, TSS_5SPL)
                               -cat 3SPL     to plot all 3' SPL sites located in the TE (i.e. for the 3 categories: 3SPL, 3SPL_exon_5SPL, 3SPL_polyA)

    OTHER OPTIONAL ARGUMENTS
     -v       =>   (BOOL)   verbose mode, makes the script talk to you (print in STDERR)
     -v       =>   (BOOL)   print version if only option
     -h,-help =>   (BOOL)   print this usage
     -chlog   =>   (BOOL)   print changelog
\n";


################################################################################
### Get arguments/options and check some of them
################################################################################
my $cat = "exonized";
my ($in,$type,$RM,$cons,$force,$TEclass,$filter,$f_regex,$min_frg,$min_len,$skip,$big,$lib,$R,$Rfile,$v,$help,$chlog);
GetOptions ('in=s' => \$in, 'type=s' => \$type, 'RM=s' => \$RM, 'filter=s' => \$filter, 'contain' => \$f_regex, 'cons' => \$cons, 'force=s' => \$force, 'skip=s' => \$skip, 'min_frg=s' => \$min_frg, 'min_len=s' => \$min_len, 'cat=s' => \$cat, 'big' => \$big, 'TEs=s' => \$TEclass, 'lib=s' => \$lib, 'R' => \$R, 'Rfile' => \$Rfile, 'chlog' => \$chlog, 'h' => \$help, 'help' => \$help, 'v' => \$v);

#check step to see if mandatory argument is provided + if help
die "\n $scriptname version: $version\n (or you forgot to specify -in, use -h to see the usage)\n\n" if ((! $in) && (! $type) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ($help);
die "\n Please provide input file (-in, use -h to see the usage)\n\n" if (! $in);
die "\n $in does not exist?\n\n" if (! -e $in);
die "\n Please set/check the type of input file (-type, use -h to see the usage)\n\n" if ((! $type) || (($type ne "RMout") && ($type ne "TEjoin") && ($type ne "TrInfo")));
die "\n -RM is mandatory when -type TrInfo is chosen (use -h to see the usage)\n\n" if (($type eq "TrInfo") && (! $RM));
die "\n $RM does not exist?\n\n" if (($RM) && (! -e $RM));
die "\n -lib set, but $lib does not exist?\n\n" if (($lib) && (! -e $lib));
die "\n -TEs set, but $TEclass does not exist?\n\n" if (($TEclass) && (! -e $TEclass));
die "\nPlease check the usage of -skip, not an integer\n\n" if (($skip) && ($skip !~ /\d/));
die "\nPlease check the usage of -min_frg, not an integer\n\n" if (($min_frg) && ($min_frg !~ /\d/));
die "\nPlease check the usage of -min_len, not an integer\n\n" if (($min_len) && ($min_len !~ /\d/));
my ($f_type,$f_name) = split(",",$filter) if ($filter);
die "\nPlease check the usage of -filter, unknown type chosen ($filter)\n\n" if (($filter) && ($f_type ne "name") && ($f_type ne "class") && ($f_type ne "family"));
die "\nPlease check the usage of -cat, unknown type chosen ($cat)\n\n" if (($cat ne "exonized") && ($cat ne "TSS") && ($cat ne "polyA") && ($cat ne "5SPL") && ($cat ne "3SPL") && ($cat ne "SPL"));
die "\nPlease check the usage of -force, unknown type chosen ($force)\n\n" if (($force) && ($force ne "Rstart") && ($force ne "Rend"));


################################################################################
### MAIN
################################################################################
########################################
# VERBOSE STUFF
########################################
print STDERR "\n--------------------------------------------------\n" if ($v);
print STDERR " --- $scriptname started (v$version), with:\n" if ($v);
print STDERR "       -in      => $in\n" if ($v);
print STDERR "                   is a repeat masker *.out file (-type RMout chosen)\n" if (($type eq "RMout") && ($v));
print STDERR "                   is a *.TEjoin.bed file (-type TEjoin chosen)\n" if (($type eq "TEjoin") && ($v));
print STDERR "                   is a *.exons.TEjoin.TrInfo.tab (-type TrInfo chosen)\n" if (($type eq "TrInfo") && ($v));
print STDERR "                    -> repeat masker output = $RM\n" if (($type eq "TrInfo") && ($v));
print STDERR "       -lib     => $lib\n" if (($lib) && ($v));
print STDERR "       -TEs     => $TEclass\n" if (($TEclass) && ($v));
print STDERR "       -cons    => additional output with all consensus lengths\n" if (($cons) && ($v));
print STDERR "       -force   => coordinates will be converted from $force only\n" if (($force) && ($v));
print STDERR "       -R       => R plots will be generated (this may not work, depends on your R install)\n" if (($R) && ($v));
print STDERR "       -Rfile   => R command lines will be printed out\n" if (($Rfile) && ($v));
print STDERR "       -big     => big output with all repeats in it to load in R\n" if (($big) && ($v));
print STDERR " --- Filtering options:\n" if ((($filter) || ($min_frg) || ($min_len)) && ($v));
print STDERR "       -filter  => $filter\n" if (($filter) && ($v));
print STDERR "                    -> using exact match\n" if (($filter) && (! $f_regex));
print STDERR "                    -> using regular expression\n" if (($filter) && ($f_regex));
print STDERR "       -skip    => instances when consensus length is more than +/- $skip nt different than the major length will be skipped\n" if (($skip) && (! $lib) && ($v));
print STDERR "       -skip    => instances when consensus length is more than +/- $skip nt different than the major length (or length calculated from $lib) will be skipped\n" if (($skip) && ($lib) && ($v));
print STDERR "       -min_frg => TEs with fewer than $min_frg fragments/intersections will be filtered out\n" if (($min_frg) && ($v));
print STDERR "       -min_len => repeats with consensus < $min_len nt will be filtered out\n" if (($min_len) && ($v));
print STDERR "       -cat     => !! -cat chosen but not -type TrInfo => ignored\n" if (($type ne "TrInfo") && ($cat ne "exonized") && ($v));
undef ($cat) if (($cat ne "exonized") && ($type ne "TrInfo"));
print STDERR "       -cat     => only $cat TEs will be mapped on the consensus\n" if (($cat) && ($v));
print STDERR "       -force   => !! -force chosen but -type RMout => ignored\n" if (($force) && ($type eq "RMout") && ($v));
undef ($force) if (($force) && ($type eq "RMout"));
print STDERR "       -force   => coordinates will be converted to start AND end in consensus based only on $force\n" if (($force) && ($v));


########################################
# ACTUAL PARSING
########################################
print STDERR "--------------------------------------------------\n" if ($v);
#Get name and path of the RMout file + create folder for outputs
my $path = path($in);
my $name = filename($in);
my $outcore = "coverage";
$outcore = $cat.".".$outcore if (($cat) && ($cat ne "exonized"));
$outcore = $f_name.".".$outcore if ($filter);
$outcore =~ s/[\[\]]//g if ($filter); #problem if brackets in folder name
print STDERR " --- Making output directory\n" if ($v);
$path = make_outdir($path,$name,$outcore);
print STDERR "      -> $path\n" if ($v);

#Get consensus lengths if -lib provided
my $replen = (); #ref of hash that will contain consensus lengths
$replen = get_fasta_lengths($lib,$v) if ($lib);

#Get TE infos if provided
my $TE = (); #will contain TE info
$TE = get_TEs_infos($TEclass,$v) if ($TEclass);

#define stuff for the subroutines
$TEclass = "n" unless ($TEclass);
$skip = "n" unless ($skip);
($cons)?($cons = "y"):($cons = "n");
($big)?($big = "y"):($big = "n");
($R)?($R = "y"):($R = "n");
($Rfile)?($Rfile = "y"):($Rfile = "n");
$filter = "n" unless ($filter);
($f_regex)?($f_regex = "y"):($f_regex = "n");
$min_frg = "n" unless ($min_frg);
$min_len = "n" unless ($min_len);
$cat = "n" unless ($cat);
$force = "n" unless ($force);

#Get the RMline when -type = TrInfo
my $RMlines = (); #will contain coords in consensus for each copy unless filtered
$RMlines = load_RM_coords($RM,$filter,$f_regex,$TEclass,$v) if ($type eq "TrInfo");

#Read input file
my ($frgs,$lib_incomplete,$which_len,$coords);
($frgs,$replen,$lib_incomplete,$which_len,$coords,$TE) = parse_in($in,$path,$type,$TEclass,$TE,$replen,$filter,$f_regex,$lib,$RMlines,$cat,$force,$v);

#Now make decision for $Rlen using the $which_len hash => fill up or update replen hash
#Also print it if -cons set
my $cout = $path."/".$name.".cons.len.tab"; 
$replen = check_Rlen($which_len,$replen,$lib,$lib_incomplete,$cons,$cout,$v);

#Now get coverage for each Rname
my ($covfiles,$bigt) = get_coverage($coords,$path,$replen,$skip,$frgs,$min_frg,$min_len,$big,$R,$v);

#now if bigtable => get all repeats in one file
my $bigout = "$path/_ALL.$outcore.tab";
get_bigtable($bigt,$bigout,$v) if ($big eq "y");

#print R command lines if relevant
get_Rcmdlines($bigt,$bigout,$replen,$covfiles,$TE,$v) if ($Rfile eq "y");

#get plots through R if relevant
get_Rplots($bigt,$replen,$covfiles,$TE,$path,$v) if ($R eq "y");

print STDERR " --- Script done -> see files in $path\n" if ($v); 
print STDERR "--------------------------------------------------\n\n" if ($v);
exit;


##########################################################################################################
# SUBROUTINES
##########################################################################################################
#----------------------------------------------------------------------------
# get a filename from a full path
# my $name = filename($filename);
#----------------------------------------------------------------------------
sub filename {
	my($file) = shift;
	$file =~ s/.*\/(.*)$/$1/;
	return $file;
}

#----------------------------------------------------------------------------
# from a filename keep only its path - independently of the current dir
# my $path = path($filename);
#----------------------------------------------------------------------------
sub path {
	my($file) = shift;
	($file =~ /\//)?($file =~ s/(.*)\/.*$/$1/):($file = ".");
	return $file;
}

#----------------------------------------------------------------------------
# make output dir
# $path = make_outdir($path,$name,$outcore);
#----------------------------------------------------------------------------
sub make_outdir {
	my ($path,$name,$outcore) = @_;
	my $c=0;
	my $full = "$path/$name.$outcore";
	until (!-d "$full.$c"){
		$c++;
	}
	mkdir ("$full.$c", 0755) or confess "ERROR (sub make_outdir): could not mkdir $path $!\n";
	return "$full.$c";
}

#----------------------------------------------------------------------------
# get_fasta_lengths (returns ref of hash, fasta IDs as keys)
# $replen = get_fasta_lengths($lib,$v);
#----------------------------------------------------------------------------
sub get_fasta_lengths {
	my ($repfile,$v) = @_;
	print STDERR " --- getting TE consensus lengths from $repfile\n" if ($v);
	my $lib = Bio::SeqIO->new(-file => $repfile, -format => "fasta") or confess "ERROR (sub get_fasta_lengths): can't open to read $repfile $!\n";
	my %replen = ();
	while( my $seq = $lib->next_seq() ) {
		my $Rname = $seq->display_id;
		my $len = $seq->length;
		$Rname =~ s/^(.*)#/$1/;
		$replen{$Rname} = $len;
	}
	undef $lib;
	return \%replen;
}

#----------------------------------------------------------------------------
# get TE infos
# $TE = get_TEs_infos($TEclass,$v) if ($TEclass);
#----------------------------------------------------------------------------
sub get_TEs_infos {
	my ($input,$v) = @_;
	print STDERR " --- getting TE infos from $input\n" if ($v);
	my %TEs = ();
	open(my $ifh, "<", $input) or confess print "ERROR (sub get_TEs_infos): could not open $input $!\n";
	LINE: while(<$ifh>) {
		chomp (my $line = $_);
		next LINE if ($line !~ /\w/);
		my @TEs = split('\t', $line); 
		$TEs{$TEs[0]} = \@TEs;
	}
	close ($ifh);
	return \%TEs;
}

#----------------------------------------------------------------------------
# load RM coords in cons if relevant
# $RMlines = load_RM_coords($RM,$filter,$f_regex,$TEclass,$v) if ($type eq "TrInfo");
#----------------------------------------------------------------------------
sub load_RM_coords {
	my ($RM,$filter,$f_regex,$TEclass,$v) = @_;
	print STDERR " --- loading coordinates in consensus for each fragment from $RM\n" if ($v);
	print STDERR "      (unless filtered out based on $filter)\n" if ($v);
	my %RMlines = ();
	open(my $ifh, "<", $RM) or confess print "ERROR (sub load_RM_coords): could not open to read $RM $!\n";
	LINE: while(<$ifh>) {
		chomp (my $line = $_);
		next LINE if (($line =~ /position|[sS]core|Gname/) || ($line !~ /\w/)); #skip headers and white lines
		$line =~ s/^\s+//; #remove spaces in beginning of lines
		my @line = split(/\s+/,$line);
		my ($Gname,$Gstart,$Gend,$Rname,$classfam,$Rend,$Rblock) = ($line[4],$line[5],$line[6],$line[9],$line[10],$line[12],$line[14]);
		
		my $TEs = (); #just so I can use same sub, but at that stage, no load in memory of the TE infos
		$TEs = get_Rclass_Rfam($Rname,$classfam,$TEclass,$TEs); 
		my ($Rclass,$Rfam) = ($TEs->{$Rname}[1],$TEs->{$Rname}[2]);
		
		my $skip = filter_TE($Rname,$Rclass,$Rfam,$filter,$f_regex); #this is to avoid loading in memory non relevant info
		next LINE if ($skip == 1);
		
		my ($Rstrand,$Rstart,$Rleft) = get_TE_cov(\@line);
		my $id = join "#",$Rname,$Gname,$Gstart,$Gend,$Rblock;
		$RMlines{$id} = \@line;
	}
	close ($ifh);
	return(\%RMlines);
}

#----------------------------------------------------------------------------
# filter out nonTE stuff and if filter
# my $skip = filter_TE($Rname,$Rclass,$Rfam,$filter,$f_regex);
#----------------------------------------------------------------------------
sub filter_TE {
	my ($Rname,$Rclass,$Rfam,$filter,$f_regex) = @_;
	my $skip = 0;
	#filter out nonTE stuff
	$skip = 1 if (($Rclass eq "nonTE") || ($Rclass eq "Simple_repeat") || ($Rclass eq "Low_complexity") 
			   || ($Rclass eq "Satellite") || ($Rclass =~ /RNA$/) || ($Rclass =~ /omeric$/));

	#filter out stuff from -filter if relevant
	my ($f_type,$f_name) = split(",",$filter) if ($filter ne "n");
	my ($lcRname,$lcRclass,$lcRfam,$lcf_name) = (lc($Rname),lc($Rclass),lc($Rfam),lc($f_name));
	if ($filter ne "n") {
		if ($f_regex eq "y") {
			$skip = 1 unless ((($f_type eq "name") && ($lcRname =~ /$lcf_name/))
						   || (($f_type eq "class") && ($lcRclass =~ /$lcf_name/)) 
						   || (($f_type eq "family") && ($lcRfam =~ /$lcf_name/)));
		} else {
			$skip = 1 unless ((($f_type eq "name") && ($lcf_name eq $lcRname))
						   || (($f_type eq "class") && ($lcf_name eq $lcRclass)) 
						   || (($f_type eq "family") && ($lcf_name eq $lcRfam)));
		}		
	}
	return ($skip)
}

#----------------------------------------------------------------------------
# Read input file => get coordinates in consensus
# ($frgs,$replen,$lib_incomplete,$which_len,$coords,$TE) = parse_in($in,$path,$type,$TEclass,$TE,$replen,$filter,$f_regex,$lib,$RMlines,$cat,$v);
#----------------------------------------------------------------------------
sub parse_in {
	my ($in,$path,$type,$TEclass,$TE,$replen,$filter,$f_regex,$lib,$RMlines,$cat,$force,$v) = @_;
	print STDERR " --- Reading $in to get coordinates in consensus\n" if ($v);
	my $lib_incomplete = 0;
	my %check = ();
	my %which_len = ();
	my %frgs = ();
	my %coords = ();
	open my $infh, "<$in" or confess "ERROR (sub parse_in): could not open to read $in $!\n";	
	LINE: while (<$infh>){
		chomp (my $line = $_);	
		next LINE if ($line !~ /\w/); #skip white lines
		my @line = ();
		my @RMline = ();
		if ($type eq "RMout") {
			next LINE if ($line =~ /position|[sS]core|Gname/); #skip headers
			$line =~ s/^\s+//; #remove spaces in beginning of lines
			@RMline = split(/\s+/,$line);
		} elsif ($type eq "TEjoin") {
			#TE join output will always look like this:
			#chr9	33165266	33165415	335;32.9;9.6;1.7;chr9;33165266;33165415;(108048016);-;MIRb;SINE/MIR;(31);237;71;2546253	.	-	chr	st	en	[uniqID	.	-]
			@line = split(/\s+/,$line);
			@RMline = split(";",$line[3]);
		} else {
			#Structure of the TrInfos file
			#0			1			2				3		4		5			6		7		8		9			10		11			12		13			14		15		16		17	
			#Gene_ID	Gene_name	Gene_biotype	Tr_ID	Tr_name	Tr_biotype	TrChr	TrStart	TrEnd	TrStrand	Trlen	TrMatureLen	ExonNb	ExonType	ExChr	ExStart	ExEnd	ExLen	
			#18					19				20					21						22						23			24		25		26		27		28		29		30			
			#len_thisTE-Exon	%thisTE-Exon	%thisTE-MatureTr	len_allTEs-MatureTr(nr)	%allTEs-MatureTr(nr)	Category	TEname	TEclass	TEfam	TEchr	TEstart	TEend	TEstrand		
			#31			32			33			34
			#RMblock	TEage(1)	TEage(2)	TE_avg_%div
			next LINE if (substr($line,0,1) eq "#"); #skip headers			
			@line = split(/\s+/,$line);
			next LINE unless (($cat eq "exonized") || ($line[23] =~ /$cat/)); #the $cat value will be contained in the value if should be kept. filtering on "exonized" means keeping all.
			#avoid counting same feature several times
			my ($st,$en,$strand) = get_cat_coord($cat,\@line);
			my $feature = join "#",$st,$en,$strand,$line[23],$line[24]; #feature coordinate and strand, TEname and category [just to be sure]
			next LINE if ($check{$feature});
			$check{$feature}=1;
			#Now carry on
			my $id = join "#",$line[24],$line[27],$line[28],$line[29],$line[31];
			#load from hash; if TE filtered out when loading infos then just fake it here (avoids putting 2 times the code to filter out)
			($RMlines->{$id})?(@RMline = @{$RMlines->{$id}}):(@RMline = ("na","na","na","na",$line[27],$line[28],$line[29],"na",$line[30],$line[24],$line[25]."/".$line[26],"na","na"));
		}
		
		#get Rclass and Rfam, to filter
		my ($Gstart,$Gend,$Rname,$classfam,$Rend) = ($RMline[5],$RMline[6],$RMline[9],$RMline[10],$RMline[12]);
		$TE = get_Rclass_Rfam($Rname,$classfam,$TEclass,$TE); #complete the hash if not all repeats were loaded from the TEclass file
		my ($Rclass,$Rfam) = ($TE->{$Rname}[1],$TE->{$Rname}[2]);

		#filter out stuff
		my $skip = filter_TE($Rname,$Rclass,$Rfam,$filter,$f_regex);
		next LINE if ($skip == 1);
				
		#Now get coordinates for the coverage, first get Rstart and Rleft based on Rstrand, valid for all
		my ($Rstrand,$Rstart,$Rleft) = get_TE_cov(\@RMline);
		my $strand = "+";
		#"correct" the $Rstart,$Rleft and $Rend for when it is -type TEjoin
		($strand,$Rstart,$Rend,$Rleft) = get_feat_cov(\@RMline,\@line,$strand,$Rstart,$Rend,$Rleft,$force) if ($type eq "TEjoin");
		#"correct" the $Rstart,$Rleft and $Rend for when it is -type TrInfo
		($strand,$Rstart,$Rend,$Rleft) = get_cat_cov(\@RMline,\@line,$Rstart,$Rend,$Rleft,$cat,$force) if ($type eq "TrInfo");
				
		#store length of consensus unless $repinfo has value for this TE (eg was in the library provided)
		my $Rlen = $Rend + $Rleft; #brackets got replaced
		unless ($replen->{$Rname}) {
			print STDERR "      !! $Rname not found in $lib => most common reconstructed length will be used\n" if (($lib) && ($v) && (! $check{$Rname}));
			$check{$Rname}=1;
			$lib_incomplete = 1;
			($which_len{$Rname}{$Rlen})?($which_len{$Rname}{$Rlen}++):($which_len{$Rname}{$Rlen}=1); #keep in mind number of "lengths" -> merging issue .align to .out
		}	
		
		#count number of fragments
		($frgs{$Rname})?($frgs{$Rname}++):($frgs{$Rname}=1);
		
		#store coordinates for each TE
		my @c = ();
		@c = @{$coords{$Rname}} unless (! $coords{$Rname}); #de-reference array that was value
		push(@c,($Rstart,$Rend,$Rleft));
		$coords{$Rname} = \@c;
	}
	close $infh;
	undef $RMlines;
	undef $TEclass;
	return ($frgs,$replen,$lib_incomplete,\%which_len,\%coords,$TE);
}

#----------------------------------------------------------------------------
# Get Rclassfam from RMout, correct if relevant
# my ($Rclass,$Rfam,$Rclassfam) = get_Rclass_Rfam($Rname,$classfam,$TEclass,$TE);
#----------------------------------------------------------------------------
sub get_Rclass_Rfam {
	my($Rname,$classfam,$TEclass,$TE) = @_;
# 	my %check = ();
	my ($Rclass,$Rfam);
	if ($classfam =~ /\//) {
		($Rclass,$Rfam) = split(/\//, $classfam);
	} else {
		$Rfam = $classfam;
		$Rfam=~ s/^(.*)\..*$/$1/;
		$Rclass = $classfam;
		$Rclass =~ s/^.*\.(.*)$/$1/;
	}
	my $Rclassfam = "$Rclass/$Rfam";
	
	#Deal with class and family, if replacement from TEclass file required
	if (($TEclass ne "n") && ($TE->{$Rname})) {
		$Rclass = $TE->{$Rname}->[1];
		$Rfam = $TE->{$Rname}->[2];
		$Rclassfam = $TE->{$Rname}->[3];	
 	} 
 	my @TE = ($Rname,$Rclass,$Rfam,$Rclassfam);
	$TE->{$Rname} = \@TE;
	return $TE;
}

#----------------------------------------------------------------------------
# Get the correct rep start and rep left / correct strand
# ($Rstrand,$Rstart,$Rleft) = get_TE_cov(\@RMline);
#----------------------------------------------------------------------------
sub get_TE_cov {
	my $l = shift;
	#($score,$div,$del,$ins,$Gname,$Gstart,$Gend,$Gleft,$strand,$Rname,$classfam,$repstart,$Rend,$repleft,$block)
	my ($Rstart,$Rleft,$Rstrand);
	$Rstrand = $l->[8];
	$Rstrand = "-" unless ($Rstrand eq "+"); #to avoid the "C"
	($Rstrand eq "+")?($Rstart = $l->[11]):($Rstart = $l->[13]);
	($Rstrand eq "+")?($Rleft = $l->[13]):($Rleft = $l->[11]);
	$Rleft =~ s/[\(\)]//g; #remove parentheses
	return ($Rstrand,$Rstart,$Rleft);
}

#----------------------------------------------------------------------------
# Get coordinates of feature in consensus if -type TEjoin
# ($strand,$Rstart,$Rend,$Rleft) = get_feat_cov(\@RMline,\@line,$strand,$Rstart,$Rend,$Rleft) if ($type eq "TEjoin");
#----------------------------------------------------------------------------
sub get_feat_cov {
	my ($RMline,$line,$strand,$Rst,$Ren,$Rlf) = @_;
	$strand = $line->[11] if ($line->[11]);
	my ($Gst,$Gen,$st,$en) = ($RMline->[5],$RMline->[6],$line->[7],$line->[8]);	
	$Rst = $Rst + ($st-$Gst) - 1 if ($Gst > $st); #if $st < $Gst then $Rst is good
	$Ren = $Ren - ($Gen-$en) if ($en < $Gen); #if $en > $Gen then $Ren is good
	$Rlf = $Rlf + ($Gen-$en) if ($en < $Gen); #if $en > $Gen then $Rleft is good
	return($strand,$Rst,$Ren,$Rlf);
}

#----------------------------------------------------------------------------
# Get coordinates of the cat stuff if -type TrInfo
# my ($st,$en,$strand) = get_cat_coord($fcat,$line);
#----------------------------------------------------------------------------
sub get_cat_coord {
	my ($fcat,$line) = @_;
	my $strand = $line->[9]; #strand of the transcript - will affect coordinates of the features
	my ($st,$en) = ($line->[15],$line->[16]); #coords of the exon
	#Categories that should not be plotted have been filtered out, so all of them need to be considered now
	#Correct for the "punctual" ones; possible values of the category are: exonized, TSS, TSS_5SPL, 5SPL, 3SPL_exon_5SPL, 3SPL, 3SPL_polyA, TSS_polyA	
	($st,$en) = ($st,$st) if ((($fcat eq "TSS") && ($strand eq "+")) || (($fcat eq "polyA") && ($strand eq "-")) || (($fcat eq "5SPL") && ($strand eq "-")) || (($fcat eq "3SPL") && ($strand eq "+")));
	($st,$en) = ($en,$en) if ((($fcat eq "TSS") && ($strand eq "-")) || (($fcat eq "polyA") && ($strand eq "+")) || (($fcat eq "5SPL") && ($strand eq "+")) || (($fcat eq "3SPL") && ($strand eq "-")));
	return($st,$en,$strand);	
}

#----------------------------------------------------------------------------
# Get coordinates of feature in consensus if -type TrInfo
# ($strand,$Rstart,$Rend,$Rleft) = get_cat_cov(\@RMline,\@line,$Rstart,$Rend,$Rleft,$cat,$force) if ($RMout eq "TrInfo");
#----------------------------------------------------------------------------
sub get_cat_cov {
	my ($RMline,$line,$Rst,$Ren,$Rlf,$fcat,$force) = @_;
	my ($st,$en,$strand) = get_cat_coord($fcat,$line);
	
	#now get the values in consensus
	my ($Gst,$Gen,$Rstrand) = ($RMline->[5],$RMline->[6],$RMline->[8]);	
	if ($Rstrand eq "+") {
		$Rst = $Rst + ($st-$Gst) if (($Gst < $st) && ($force ne "Rend")); #if $st <= $Gst then $Rst is good
		$Rst = $Ren - ($Gen-$st) if (($Gst < $st) && ($force eq "Rend")); #if $st <= $Gst then $Rst is good

		$Ren = $Ren - ($Gen-$en) if (($en < $Gen) && ($force ne "Rstart")); #if $en >= $Gen then $Ren is good				
		$Ren = $Rst + ($en-$Gst) if (($en < $Gen) && ($force eq "Rstart")); #if $en >= $Gen then $Ren is already good	

		$Rlf = $Rlf + ($Gen-$en) if (($en < $Gen) && ($force ne "Rstart"));
		$Rlf = $Rlf - ($en-$Gst) if (($en < $Gen) && ($force eq "Rstart"));		
	} else {
		$Rst = $Rst + ($Gen-$en) if (($en < $Gen) && ($force ne "Rend")); 
		$Rst = $Ren - ($Gen-$en) if (($en < $Gen) && ($force eq "Rend")); 

		$Ren = $Ren - ($st-$Gst) if (($Gst < $st) && ($force ne "Rstart"));
		$Ren = $Rst + ($Gen-$st) if (($Gst < $st) && ($force eq "Rstart"));		
		
		$Rlf = $Rlf + ($st-$Gst) if (($Gst < $st) && ($force ne "Rstart"));
		$Rlf = $Rlf - ($Gen-$st) if (($Gst < $st) && ($force eq "Rstart"));		
	}
	return($strand,$Rst,$Ren,$Rlf);
}

#----------------------------------------------------------------------------
# get Rlen from hash, returns ref of hash
# $replen = check_Rlen($which_len,$replen,$lib,$lib_incomplete,$cons,$cout,$v);
#----------------------------------------------------------------------------
sub check_Rlen {
	my ($which_len,$replen,$lib,$lib_incomplete,$cons,$cout,$v) = @_;	
	print STDERR " --- Chosing most common consensus lengths\n" if (($v) && (! $lib));
	print STDERR " --- Chosing most common consensus lengths [not all sequences were in $lib]\n" if (($v) && ($lib) && ($lib_incomplete == 1));
	print STDERR "     And printing all lengths from the repeat masking info in:\n     $cout\n" if (($v) && ($cons eq "y"));
	open my $cfh, ">", $cout or confess "ERROR (sub check_Rlen): Failed to open to write $cout $!\n" if ($cons eq "y");
	foreach my $rep (sort keys %{$which_len}){
		my @lengths = sort {$which_len->{$rep}{$b} <=> $which_len->{$rep}{$a}} keys %{$which_len->{$rep}};
		$replen->{$rep}=$lengths[0] unless ($replen->{$rep}); #most common length
		print $cfh "$rep\t@lengths\n" if ($cons eq "y");
	}
	close $cfh if ($cons eq "y");
	return ($replen);
}

#----------------------------------------------------------------------------
# Coverage for each Rname
# my ($covfiles,$big) = get_coverage($coords,$path,$replen,$skip,$frgs,$min_frg,$min_len,$big,$R,$v);
#----------------------------------------------------------------------------
sub get_coverage {
	my ($coords,$path,$replen,$skip,$frgs,$min_frg,$min_len,$big,$R,$v) = @_;
	print STDERR " --- Going through arrays of coordinates to get coverages\n" if ($v);
	my %big = ();# if ($big eq "y");
	my %covfiles = ();# if ($R eq "y");
	RNAME: foreach my $Rname (keys %{$coords}) {
		my $clen = $replen->{$Rname}; # get consensus length
		next RNAME if (($min_frg ne "n") && ($frgs->{$Rname} < $min_frg)); #skip this one if not enough fragments to see a coverage trend [defined by user]
		next RNAME if (($min_len ne "n") && ($clen < $min_len)); #skip this one if not long enough [defined by user]

		#initialize coverage array
		my @cov = ();
		for (my $c = 0; $c < $clen; $c++){
			push(@cov, 0);
		}
		
		#get the list and fill+print coverage
		my @Rcoords = @{ $coords->{$Rname} }; #de-reference array that was value
		COORDS: for (my $i=0, my $j=1; $i<$#Rcoords+1; $i+=3,$j+=3){ #$j<$#Rcoords+2
			my ($s,$e) = ($Rcoords[$i], $Rcoords[$j]);
			#check if should be skipped
			if ($skip ne "n") {
				my $l = $Rcoords[$j+1];
				my $Rlen = $e+$l;
				next COORDS if (($Rlen < $clen - $skip) || ($Rlen > $clen + $skip));
			}	
			foreach my $k ($s-1..$e-1){
				$cov[$k]++;
			}
		}	

		my $out = "$path/$Rname.coverage.tab";
		$covfiles{$Rname} = $out;
		open my $ofh, ">$out" or confess "ERROR (sub get_coverage): could not open to write $out $!\n";
		for (my $i=0; $i<=$#cov; $i++) { 
			print $ofh "$cov[$i]\n";		
		}
		close $ofh;
		$big{$Rname}=\@cov if ($big eq "y"); 	#if bigtable option then memorize all coverages
		
	}
	return(\%covfiles,\%big);
}

#----------------------------------------------------------------------------
# Get bigtable output if relevant
# get_bigtable($bigt,$bigout,$v) if ($big eq "y");
#----------------------------------------------------------------------------
sub get_bigtable {
	my ($big,$bigout,$v) = @_;
	print STDERR " --- Getting big file with coverage for all repeats\n" if ($v);
	
	#get big table with everything
	my @coverage = ();
	foreach my $key (sort keys %{$big}) {
		my @cov = @{ $big->{$key} }; #de-reference array that was value
		push(@coverage,[$key,@cov]); #now push that
	}
	
	#Transpose table to invert columns and rows - use ragged because it is an irregular matrix
	my @TR_coverage = transpose_ragged(\@coverage);
		
	#Get dimensions of the matrix
	my $catnb = @coverage;
	my $maxlength = @TR_coverage;

	my $col; #= number cat for this repeat
	my $line;#= number of max length
		
	#create a big outputfile	
	open my $ofh, ">$bigout" or confess "ERROR (sub get_bigtable): Failed to open to write $bigout $!\n";
	for ($line=0; $line<=$maxlength; $line+=1) { #for each line,
		for ($col=0; $col<=$catnb; $col+=1) { #print all columns
			if (defined ($TR_coverage[$line][$col])) {
				print $ofh "$TR_coverage[$line][$col]\t";			
			} else {
				print $ofh "0\t";
			}
		}
		print $ofh "\n";
	}
	close $ofh;
}

#----------------------------------------------------------------------------
# print R command lines to plot coverages
# get_Rcmdlines($bigt,$bigout,$replen,$covfiles,$TE,$v) if ($Rfile);
#----------------------------------------------------------------------------
sub get_Rcmdlines {
	my ($big,$bigout,$replen,$covfiles,$TE,$v) = @_;
	print STDERR " --- Printing R command lines\n" if ($v); 
	#output file	
	my $out = "$path/_ALL.R.cmdlines.txt";
	open my $ofh, ">$out" or confess "ERROR (sub get_Rcmdlines): could not open to write $out $!\n";
	print $ofh "setwd(\"\") #fill with location of files on computer\n\n";
	print $ofh "dat<-read.csv(\"_ALL.R.cmdlines.txt\", sep=\"\t\", header = TRUE)\n" unless ($bigout eq "n");
	#now loop
	my $nb = 1;
	my $i = 1;
	foreach my $Rname (keys %{ $big }) {
		my $filename = $covfiles->{$Rname};
		#deal with consensus length and tick marks
		my $len = $replen->{$Rname};
		my $ticks = int( $len/100 );
		$ticks = 10 if ($ticks < 10);
		my $fullname = $Rname."#".$TE->{$Rname}[1]."/".$TE->{$Rname}[2];
		#now print
		print $ofh "\npar(mfrow=c(4,2))\n" if ($nb == 1);		
		if ($bigout eq "n") {
			print $ofh "dat<-read.csv(\"$filename\", header = FALSE)\n";
			print $ofh "plot(dat\$V1, col=\"blue\",pch = 18, xlab = \"position in consensus\", ylab = \"coverage\", main=\"$fullname\", xaxp = c(0,$len,$ticks))\n";
		} else {
			print $ofh "plot(dat\$V$i, col=\"blue\",pch = 18, xlab = \"position in consensus\", ylab = \"coverage\", main=\"$fullname\", xaxp = c(0,$len,$ticks))\n";
		}
		$i++;	
		$nb = 0 if ($nb == 8);
		$nb++;
	}	
	close $ofh;
#Table[!is.na(Table$Col),])
}

#----------------------------------------------------------------------------
# get R plots
# get_Rplots($bigt,$replen,$covfiles,$TE,$path,$v) if ($R);
#----------------------------------------------------------------------------
sub get_Rplots {
	my ($big,$replen,$filelist,$TE,$path,$v) = @_;	
	print STDERR " --- Getting R plots\n" if ($v); 
	#Start R bridge
	my $R = Statistics::R->new() ;
	$R->startR ;

	foreach my $Rname (keys %{ $big }) {
		my $filename = $filelist->{$Rname};
		
		#get infos
		my $fullname = $Rname."#".$TE->{$Rname}[1]."/".$TE->{$Rname}[2];		
		my $len = $replen->{$Rname};
		
		#plot
		$R->send(qq`dat <- read.table("$filename", header = FALSE)`);
		my $out = "$filename.pdf"; 
		#$R->run(qq`pdf("$out")`);
		$R->run(qq`pdf("$out")`);	
		
# 		# PDF Format  
# 		# The 'pdf' device allows for PDF output.  
# 		#  
# 		# Height and Width: Dimensions of the image in inches  
# 		# Onefile: If true, it plots multiple figures on a page. If false,  
# 		#      each figure gets its own numbered page.  
# 		# Family: The font family used to render text. You can use most fonts  
# 		#     available to your computer.  
# 		# Paper: The type of paper determines the pdf size.  
# 		# Pointsize: Font size.  
# 		pdf(file='PDF Output.pdf', height=6, width=6, onefile=TRUE, family='Helvetica', paper='letter', pointsize=12) 
		
		
		my $V1 = "\$V1";
		$R->send(qq`plot(dat$V1, col="blue", pch = 18, xlab = "position in consensus (tot = $len)", ylab = "coverage", main="$fullname")`);
		$R->run(q`dev.off()`);
	}	
	#End R bridge
	$R->stopR() ;
}

