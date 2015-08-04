#!/usr/bin/perl -w
#######################################################
# Author  :  Aurelie K
# version :  2.0 (see below)
# email   :  4urelie.k@gmail.com  
# Purpose :  It reads a repeat masker output file and extracts sequences (all or subset based on filter and/or random)
#			 Check Usage for the many options
#######################################################
BEGIN{
   #what to do on: kill -s ALRM <pid> so I can check where it is if it stalls
   $SIG{ALRM}  = sub {print STDERR "SIGALRM received\n"; print STDERR Carp::longmess; print "\n";};
	unshift(@INC, "~/bin/BioPerl-1.6.901");
}

#always load forks before anything else
use forks;
use forks::shared;
use Thread::Semaphore;
#now load rest
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Bio::Perl;
use Bio::DB::Fasta;
use List::Util 'shuffle';

my $version = "2.13";

# UPDATES
my $changelog = "
#   - v1.0 = Jun 2011
#   - v1.1 = May 2014 (other script named parseRM-ExtractSeq2-ak.pl)
#		Bug fix - avoid extractions error because of chr17 included in chr17_gl000203_random for ex
#   - v2.0 = Jun 2014
#		Restructuration (subs) + have Rname in fasta_id, better
#		Options: -append, -filter, -min_frg, -min_len, -TEs, -v, multiple RMout files...
#		Parallelization
#   - v2.1 = Jul 2014
#		Return in thread loop that shouldn't be there -> changed as a next FILE
#		Put extracted sequences in a folder
#   - v2.2 = 10 Jul 2014
#		Corrected cleaning previous outputs
#		Little bug fixed at several places
#		chlog option
#   - v2.3 = 15 Jul 2014
#		-extract for a subset of sequences was printing duplicates -> fixed
#		-extract modified => possibility to also get complementary set than the one extracted
#       Avoid having ZZZ:XX-XX structure in SeqIDs -> problems with Bio::DB::Fasta later
#   - v2.4 = 18 Jul 2014
#		Bug fix in -flank option
#		   -> correct coordinates for ends of scaffold + keep track of actual fl length in description
#		   -> had to introduce changes in how sequence is store to append fragments to each others
#		Bug fix in extract_random sub
#   - v2.5 = 29 Jul 2014
#		Typos; \$se (semaphore that I don't use anymore) remained in getting arguments of sub extract_one_seq
#	- v2.6 = 31 Jul 2014
#		extract_one_seq was modified to avoid regexp against whole db header list, it was a long useless step
#		\"random\" in name of complementary set file name, removed
#	- v2.7 = 27 Aug 2014
#		option to keep stuff with the class named \"nonTE\"
#       option cat to cat the extracted fasta files
#		Separation of fa files by RMout file (were all put in the same directory...)
#   - v2.8 = 07 Jan 2015
#       option -rc
#   - v2.9 = 09 Jan 2015
#       options -f5 and -f3 [TO DO: change it to when extraction is done so that filtering is done after reconstructing repeats]
#   - v2.10 = 15 Jan 2015
#       Bio::DB isssues, index during threading -> has to be done before
#   - v2.11 = 23 Jan 2015
#       -min_div and -max_div
#   - v2.12 = 03 Feb 2015
#       Correct starting threads, was starting number asked +1
#       remove filtering out nested TEs (was a cc leftover...)
#   - v2.14 = 04 Mar 2015
#       Bug fix in filtering + added -contain
\n";

my $usage = "\nUsage [$version]: 
    perl <scriptname.pl> -dir <dir_with_all_files> [-min_frg <X>] [-min_len <X>] [-f5 <X>] [-f3 <X>] [-min_div <X>] [-max_div <X>] [-TEs <TEclass>] [-filter <type,name>] [-contain] [-nonTE] [-flank <X>] [-rc] [-presplit] [-append] [-extract <type,X>] [-cat] [-cpu <X>] [-v]
	
	This script will create a file per repeat (Rname) with extracted sequences + a concatenated fasta file
	Names will be >Rname#Rclass/Rfam_Chr-start-end genome_file_of_origin

    MANDATORY ARGUMENT:	
    -dir <dir_with_all_files>
        directory that should contain:
           - at least 1 repeat masker output [there can be several]
           - Corresponding genomes (exact same ones that are masked) 
               /!\\ Same name is required before extensions = .out and .fa; note that .out.XXXX is OK if you need more info in the file name
               NOTE: sym links will work (ln -s path/file new_name) so no need to actually copy the genome files
           - previous repeats.fa files with extracted sequences with this script if you want to append new sequences 
             /!\\ Use -append or they will be deleted
        For ex. minimal list of files in this directory = species1.out, species1.fa
        Outputs will be located in this folder

	  
    OPTIONAL ARGUMENTS (flagged by brackets [...] around each of them)
    FILTERING OPTIONS
    -min_frg <X>
        filter for a minimum number of fragments (useful if masking with unclassified things that have very low copy nb)
           ex: -min_frg 3   
    -min_len <X>
        X = minimum length of the fragment (in nt) to be retained (after reconstruction based on interrupted repeats / blockID of Repeat Masker)
           Note that this reconstruction in RM software requires repeats to be named like: Rname#classfam
           ex: -min_len 80
    -max_div <X>
        will extract only FRAGMENT that have %div < X [TO DO: do this on average of all fragments, after reconstruction]
    -min_div <X>
        will extract only FRAGMENT that have %div > X [TO DO: do this on average of all fragments, after reconstruction]
    -f5 <X> 
        to skip Repeat Masker output line when start in repeat is more than X
        For example, to extract fragment only if the 5' end of the consensus was found, type -f5 1
        !!! If you set both -f5 and -f3, interrupted repeats that indeed have both ends present will still be skipped
    -f3 <X>
        to skip Repeat Masker output line when what's left of the consensus of the repeat is more than X
        For example, to extract fragment only if the 3' end of the consensus was found, type -f3 0
        /!\\ If you set both -f5 and -f3, interrupted repeats that indeed have both ends present will still be skipped
    -TEs <TEclass>
        This optional input file allows to correct class and/or family for (some) repeats, before filtering
           Note that you would need to \"cat\" the various files, if you had more than one (e.g. for several species etc)
           Should contain TE information as follow, 3 first columns mandatory: 
           Rname \\t Rclass \\t Rfam \\t Rclass/Rfam \\t etc (only first 3 columns will be used)
           an easy way to get it is to run my other script parseRM.pl.
           -> copy the first columns of the output all-repeats.tab, modify them accordingly and then copy this in a text file
    -filter <type,name>
        run the script on only a subset of repeats. Not case sensitive.
        This is advised instead of grepping the Repeat Masker output, because elements will be skipped when nested (e.g. flankings are repeated as well).
        The type can be: name, class or family and it will be EXACT MATCH unless -contain is chosen as well
        ex: name,nhAT1_ML => only fragments corresponding to the repeat named exactly nhAT1_ML will be looked at
            class,DNA => all repeats with class named exactly DNA (as in ...#DNA/hAT or ...#DNA/Tc1)
            family,hAT => all repeats with family named exactly hAT (so NOT ...#DNA/hAT-Charlie for example)
    -contain
        to check if the \"name\" determined with -filter is included in the value in Repeat Masker output, instead of exact match
        ex: name,HERVK => all fragments containing HERVK in their name
            family,hAT => all repeats with family containing hAT (...#DNA/hAT, ...#DNA/hAT-Charlie, etc)
    -nonTE <type>
        chose this option if you wish to keep (some) nonTE sequences
        type --> all, nonTE
          all   = KEEP ALL nonTE stuff (Low_complexity, Simple_repeat, Satellites, etc)
          nonTE = KEEP only non TE sequences when class = nonTE (other stuff will be filtered out)

    OTHER OPTIONS
    -flank <X>
        add X nt in 5' and X nt in 3' of the sequences 
    -rc
        to extract reverse complement of the sequence when TE is not on strand +
    -presplit
    	if RMoutput files have been split already => do not do it again
    -append
    	existing fasta files with extracted sequences lcoated in the folder won't be deleted, and newly extracted sequences will be added
    -extract <type,X>
    	To extract a subset (if not chosen, all sequences are extracted)
        - type --> all, only, compl
          all = extract all AND extract X in additional files
          only = extract only these X sequences
          compl = extract all AND extract X in additional files AND extract complementary set
        - X = number of a subset of sequences per repeat and per RM output that will be extracted
          Note that ~5% highest %div as well as ~5% lowest %div will be extracted to capture diversity; rest will be random
          ex: -extract only,100  
    -cat 
        concatenate extracted fasta files
    -cpus <X>
        X = max number of parallel jobs running (CPUs used)
        default = 1 (e.g. no parallelization)
    -chlog
    	print updates
    -v
        verbose mode, make the script talks to you
        print version if only option
    -h|help 
        Print this help\n\n";
        
        
################################################################################
# Get arguments/options, check some of them, define shared stuff and start threads
################################################################################
my ($flank,$min_frg,$min_len,$filter,$filter5,$filter3,$extract,$nonTE,$maxdiv,$mindiv) = ("na","na","na","na","na","na","na,na","na","na","na");
my ($dir,$presplit,$append,$f_regexp,$TEclass,$help,$v,$chlog,$cat,$rc);
my $cpus = 1;
GetOptions ('dir=s' => \$dir, 'presplit' => \$presplit, 'append' => \$append, 'extract=s' => \$extract,'flank=s' => \$flank,'rc' => \$rc,'filter=s' => \$filter, 'contain' => \$f_regexp, 'f5=s' => \$filter5, 'f3=s' => \$filter3, 'nonTE=s' => \$nonTE, 'min_frg=s' => \$min_frg, 'min_len=s' => \$min_len, 'max_div=s' => \$maxdiv, 'min_div=s' => \$mindiv, 'TEs=s' => \$TEclass, 'cat' => \$cat, 'cpus=s' => \$cpus, 'chlog' => \$chlog, 'h' => \$help, 'help' => \$help, 'v' => \$v);

($f_regexp)?($f_regexp = "y"):($f_regexp="n");

#check step to see if mandatory argument is provided + if help/changelog
die "\n version $version\n\n" if ((! $dir) && (! $help)  && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ((! $dir) || ($help));

print STDERR "\n --- Script parseRM_ExtractSeqs.pl started (v$version)\n" if ($v);

#avoid / at the end of paths
$dir = $1 if ($dir =~ /^(.*)\/$/);
print STDERR "      - Directory containing input files = $dir\n" if ($v);

print STDERR "      - WARNING => -nonTE option chosen but -presplit chosen as well: previous files will be used.\n         Don't chose -presplit if you want to change the filtering options\n" if (($filter ne "na") && ($presplit));

#if relevant, check if filter has the ,
die "\nERROR (main): check -filter option usage (use -h if you need to see the usage)\n" if (($filter ne "na") && ($filter !~ /,/));
die "\nERROR (main): check -extract option usage (use -h if you need to see the usage)\n" if (($extract ne "na") && ($extract !~ /,/));

#clean previous files unless -append or -presplit + hive them values
unless ($append) {
	print STDERR "      - Clean previous output directories\n" if ($v);
	my @fa = `ls $dir | grep Extract`;
	my $checkfa = @fa;
	if ($checkfa > 0) {
		print STDERR "        rm -Rf $dir/*.Extract\n" if ($v);
		`rm -Rf $dir/*.Extract`; #or die "\n      ERROR (main, prep): can not rm -Rf $dir/*.Extract $!\n\n";
	}
	my @checksplit = `ls $dir | grep split`;
	my $checksplitnb = @checksplit;
	if ($checksplitnb > 0) {
		my @checktab = `ls $dir/*.split | grep tab`;
		my $checktabnb = @checktab;
		if ($presplit) {
			if ($checktabnb > 0) {
				print STDERR "        rm -f $dir/*.split/*tab\n" if ($v);
				`rm -f $dir/*.split/*tab`; #or die "\n      ERROR (main, prep): can not rm -f tab files in $dir/*.split directories $!\n\n";
			}	
		} else {
			print STDERR "        rm -Rf $dir/*.split\n" if ($v);
			`rm -Rf $dir/*.split`; #or die "\n      ERROR (main, prep): can not rm -Rf $dir/*.split $!\n\n";
		}	
	}
} else {
	print STDERR "      - Option -append chosen, no cleaning previous outputs\n" if ($v);
}

#make index files for the genomes
print STDERR "      - Indexing genome files\n" if ($v);
my @Genomes = `ls $dir/*.fa` or die "\n      ERROR (main, prep): can't list files .fa in $dir $!\n\n";
foreach my $g (@Genomes) {
	chomp ($g);
	# index the genomes if necessary
	my $index = "$g.index";
	if (-e $index) {
		print STDERR "         -> $index file exists, skipped\n" if ($v);
	} else {
		my $db = Bio::DB::Fasta->new($g, -reindex => 1) or confess "\n      ERROR (Main): Failed to open create Bio::DB object (to index) from $g $!\n\n";
		print STDERR "         -> $index file created\n" if ($v);
		undef ($db);
	}	
}		

#keep STDOUT and STDERR from buffering
select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately
select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately

#Initialize and open Thread stuff
my @RMout_list :shared; 
my @RMout_done :shared;
my @posi_list :shared; 
my @posi_done :shared; 
my %frgs_all :shared;
my %frgs_nr :shared;
my $finished :shared; 

################################################################################
# MAIN
################################################################################
print STDERR "      - Max number of CPUs used = $cpus\n" if ($v);

($rc)?($rc="y"):($rc="n");

#start threads
print STDERR "       => Starting $cpus threads\n" if ($v);
for(my $i = 1; $i < $cpus; $i++){
    threads->create({ 'context' => 'scalar' }, \&thread, \@RMout_list, \@RMout_done, \%frgs_all, \%frgs_nr, \@posi_list, \@posi_done, \$finished, \$flank, \$rc, \$min_frg, \$min_len, \$extract, \$dir);
}

#Get TE infos if provided
print STDERR "\n --- getting TE infos from $TEclass...\n" if (($TEclass) && ($v));
my $TE->{"na"} = "na"; #will be replaced by TE info if TEclass provided
$TE = get_TEs_infos($TEclass) if ($TEclass);
print STDERR "     ...getting TE infos done\n" if (($TEclass) && ($v));
$TEclass = "na";

#Split RMout files per repeat and get list of RMout files
if ($presplit) {
	print STDERR " --- no splitting RM output files (-presplit chosen) => gathering previous split RMout in $dir/*split\n" if ($v);
	my @RMout_files = `find $dir/*split -type f -name \"*.out*##*\"` or die "\n      ERROR (main): can't find files *.out*##* in $dir $!\n\n";
	my $list = check_for_genome(\@RMout_files,$dir,$v); #check that corresponding genome exists before doing anything
	@RMout_list = @{$list};
} else {
	print STDERR " --- Splitting RM output files and getting list of RM output files to process...\n" if ($v);
	my @RMout_files = `ls $dir/*.out*` or die "\n      ERROR (main): can't list files .out in $dir $!\n\n";
	my $split_RMouts = split_RM_files($dir,\@RMout_files,$filter,$f_regexp,$filter5,$filter3,$maxdiv,$mindiv,$TE,$TEclass,$nonTE,$v); #checking that corresponding genome exists is done in sub (gain of time, less files)
	@RMout_list = @{$split_RMouts};
	print STDERR "     ...Splitting done\n" if ($v);
}

#run threads
print STDERR " --- Processing RM output files [threading if -cpus chosen]...\n" if ($v);
print STDERR "     (Note that repeats with same repeat masker block ID will be reconstructed when extracted)\n" if ($v);
print STDERR "     With filter = $filter\n" if (($v) && ($filter ne "na"));
print STDERR "     With $flank nt in 5' and 3'\n" if (($v) && ($flank ne "na"));
print STDERR "     With -f5 set ($filter5)\n" if (($v) && ($filter5 ne "na"));
print STDERR "     With -f3 set ($filter3)\n" if (($v) && ($filter3 ne "na"));
print STDERR "     Only when repeat has more than $min_frg fragments\n" if (($v) && ($min_frg ne "na"));
print STDERR "     Only when reconstructed length is > $min_len nt\n" if (($v) && ($min_len ne "na"));
print STDERR "     Only when % divergence of the fragment is < $maxdiv\n" if (($v) && ($maxdiv ne "na"));
print STDERR "     Only when % divergence of the fragment is > $mindiv\n" if (($v) && ($mindiv ne "na"));
print STDERR "     Reverse complement will be extracted when TE is on minus strand (-rc chosen)\n" if (($v) && ($rc eq "y"));

#flag finished now so that "log" above is printed correctly
$finished = 1;
thread(\@RMout_list, \@RMout_done, \%frgs_all, \%frgs_nr, \@posi_list, \@posi_done, \$finished, \$flank, \$rc, \$min_frg, \$min_len, \$extract, \$dir, $v);

#clean threads
print STDERR "   - Cleaning the $cpus threads\n" if ($v);
my $totflag = 0;
foreach my $thr (threads->list){
    my $flag = $thr->join();
    $totflag+=$flag;
}
my $totRMdone = 0;
foreach my $done (@RMout_done) {
	$totRMdone++;
}
my $totP = 0;
foreach my $p (@posi_list) {
	$totP++;
}
my $totPdone = 0;
foreach my $p (@posi_done) {
	$totPdone++;
}
#Just check that files were processed and same number of position files generated
print STDERR "     => $totRMdone RMout files processed\n" if ($v);
print STDERR "        $totP position files generated\n" if ($v);
print STDERR "     => total of $totPdone position files processed to extract sequences\n" if ($v);
print STDERR "     !! Some files were not processed ($totflag)\n" unless ($totflag == 0);

my @Extract = `ls $dir | grep Extract`;
cat_fa(\@Extract,$dir,$v) if ($cat);

print STDERR "\n --- Script done\n" if ($v);
print STDERR "     --> see files .fa in $dir/Extract\n" if ($v);
print STDERR "         + see concatenated fasta file(s) in $dir\n" if (($cat) && ($v));
print STDERR "\n" if ($v);
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
# from a filename or a directory keep only its path - independently of the current dir
# my $path = path($filename);
#----------------------------------------------------------------------------
sub get_path {
	my($file) = shift;
	($file =~ /\//)?($file =~ s/(.*)\/.*$/$1/):($file = ".");
	return $file;
}

#----------------------------------------------------------------------------
# get TE infos
# $TE = get_TEs_infos($TEclass)
#----------------------------------------------------------------------------
sub get_TEs_infos {
	my $input = shift; #file name
	my %TEs = ();
	open(my $input_fh, "<", $input) or confess print STDERR "ERROR (sub get_TEs_infos): could not open $input $!\n";
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
# get TE infos
# my $list = check_for_genome(\@RMout_files,$dir,$v)
#----------------------------------------------------------------------------
sub check_for_genome {
	my ($files,$dir,$v) = @_;
	my @list = ();
	foreach my $file (@$files) {
		chomp $file;
		my $fa = filename($file);
		$fa = $1.".fa" if $fa =~ /^(.*)\.out/;
		unless (-e "$dir/$fa") {
			print STDERR "     ERROR (sub check_for_genome): $fa doesn't exist - check genome file name (see usage)\n" if ($v);
			print STDERR "       -> $file skipped for analysis\n" if ($v);
		} else {
			push(@list,$file);
		}
	}
	return (\@list);
}

#----------------------------------------------------------------------------
# split RMoutput files and filter them - so that it's faster to load in arrays
# my $split_RMouts = split_RM_files($dir,\@RMout_files,$filter,$f_regexp,$filter5,$filter3,$maxdiv,$mindiv,$TE,$TEclass,$nonTE,$v);
#----------------------------------------------------------------------------
sub split_RM_files {
	my ($dir,$RMouts,$filter,$f_regexp,$filter5,$filter3,$maxdiv,$mindiv,$TE,$TEclass,$nonTE,$v) = @_;	
	my @list = ();

	#get filter if relevant
	my ($f_type,$f_name) = split(",",$filter) unless ($filter eq "na");
	my %check = ();	
	#check that corresponding genome exists, if not exclude from list
	my $RMs = check_for_genome($RMouts,$dir,$v);
	
	print STDERR "     [Elements with class = nonTE won't be filtered out (-nonTE nonTE chosen)]\n" if (($v) && ($nonTE eq "nonTE"));
	print STDERR "     [Non TE repeats won't be filtered out (-nonTE all chosen)]\n" if (($v) && ($nonTE eq "all"));
	
	FILE: foreach my $file (@{$RMs}) {
		chomp $file;
		next FILE unless (-e $file); #for some reason there was an extra empty value in the list. Better check on file existence.
		print STDERR "     -> $file...\n" if ($v);
		#get genome name
		my $fa = $file;
		$fa = $1.".fa" if $fa =~ /^(.*)\.out/;
		
		#make outputdir
		my $filename = filename($file);
		my $splitpath = $dir."/".$filename.".split";
		system "mkdir $splitpath" unless (-e $splitpath);
		
		#put RMout in array to be able to loop on it (because need to access to previous and next $block)
		my @RMout = get_RMout_array($file);		
		#now loop on array
		my %check = ();
		LINE: for (my $i = 0; $i <= $#RMout; $i++){
			my ($Rname,$classfam,$block) 
			   = ($RMout[$i]->[9],$RMout[$i]->[10],$RMout[$i]->[14]);	
			next LINE unless ($block); #apparently, that happens	
			my ($Rclass,$Rfam,$Rclassfam) = get_Rclass_Rfam($classfam);
 			
			#Deal with class and family, if replacement required
			if ($TEclass ne "na") {
				if ($TE->{lc($Rname)}) {
					$Rclass = $TE->{lc($Rname)}->[1];
					$Rfam = $TE->{lc($Rname)}->[2];
					$Rclassfam = $TE->{lc($Rname)}->[3];	
				} else {
					print STDERR "        !! $Rname not found in $TEclass => RM output class / fam used instead\n" if (($TEclass) && ($v) && (! $check{lc($Rname)}{"class"}));
					$check{lc($Rname)}{"class"}=1;
				}	
			} 
			#filter out non TE stuff unless asked not
			unless (($nonTE eq "nonTE") ||  ($nonTE eq "all")) {
				next LINE if ($Rclass eq "nonTE");
			}
			unless ($nonTE eq "all") {
				next LINE if (($Rclass eq "Simple_repeat") || ($Rclass eq "Low_complexity") 
							   || ($Rclass eq "Satellite") || ($Rclass =~ /RNA$/) || ($Rclass =~ /omeric$/) || ($Rclass eq "ARTEFACT"));
			}

			#filter out stuff if relevant		
			unless ($filter eq "na") {
				my ($lcRname,$lcRclass,$lcRfam,$lcf_name) = (lc($Rname),lc($Rclass),lc($Rfam),lc($f_name));
				if ($f_regexp eq "y") {
					#check if what's det as the filter is included in the names
					next LINE unless ((($f_type eq "name") && ($lcRname =~ /$lcf_name/))
					               || (($f_type eq "class") && ($lcRclass =~ /$lcf_name/)) 
					               || (($f_type eq "family") && ($lcRfam =~ /$lcf_name/)));
				} else {
					next LINE unless ((($f_type eq "name") && ($lcf_name eq $lcRname))
					               || (($f_type eq "class") && ($lcf_name eq $lcRclass)) 
					               || (($f_type eq "family") && ($lcf_name eq $lcRfam)));
				}	
			}
			
			#filter max and/or min div
			my $div = $RMout[$i]->[1];
			next LINE if ((($mindiv ne "na") && ($div < $mindiv)) || (($maxdiv ne "na") && ($div > $maxdiv)));
			
			#filter for presence of ends (relatively to consensus)
			my ($strand,$repstart,$repleft) = ($RMout[$i]->[8],$RMout[$i]->[11],$RMout[$i]->[13]);	
			my $Rstart = $repstart;
			my $Rleft = $repleft;
			unless ($strand eq "+") {
				$Rstart = $repleft;
				$Rleft = $repstart;
			}
			$Rleft = $1 if $Rleft =~ /\(([0-9]+)\)/;
			next LINE if ((($filter5 ne "na") && ($Rstart > $filter5)) || (($filter3 ne "na") && ($Rleft > $filter3)));
			
			#now print in split files what went through
			my $split = $splitpath."/".$filename."##$Rname";
			push(@list,$split) unless (-e $split);
			my $line = "$RMout[$i]->[0]\t$RMout[$i]->[1]\t$RMout[$i]->[2]\t$RMout[$i]->[3]\t$RMout[$i]->[4]\t$RMout[$i]->[5]\t$RMout[$i]->[6]\t$RMout[$i]->[7]\t$RMout[$i]->[8]\t$RMout[$i]->[9]\t$RMout[$i]->[10]\t$RMout[$i]->[11]\t$RMout[$i]->[12]\t$RMout[$i]->[13]\t$RMout[$i]->[14]\n";
			open (my $fh_split,">>",$split) or confess "\n      ERROR (Sub split_RM_files): Failed to open $split $!\n\n";
			print $fh_split "$line";
			close($fh_split);
		}
		print STDERR "     ...done\n" if ($v);
	}
	return \@list;
}

#----------------------------------------------------------------------------
# Threaded actions, loop on .out files split by repeat name
# thread(\@RMout_list, \@RMout_done, \%frgs_all, \%frgs_nr, \@posi_list, \@posi_done, \$finished, \$flank, \$min_frg, \$min_len, \$extract, \$dir, \$v);
#----------------------------------------------------------------------------
sub thread {
    my ($RMout_list,$RMout_done,$frgs_all,$frgs_nr,$posi_list,$posi_done,$finished,$flank,$rc,$min_frg,$min_len,$extract,$dir,$v) = @_; 
	
	my $allflag = 0;
    FILE: while((my $RMout = shift @$RMout_list) || (! $$finished)){
		if(!$RMout){
			sleep 1;
			next;
		}
		chomp ($RMout);
		print STDERR "     -> $RMout in progress...\n" if ($v);
		
		#get positions
		my $posi;
		($posi_list,$posi) = get_posi_from_RMout($RMout,$frgs_all,$frgs_nr,$posi_list,$v);
		push(@{$RMout_done},$RMout);
		
		#check posifiles for frg numbers
		my $posiname = filename($posi);
		unless (($$min_frg eq "na") || ($frgs_nr->{$posiname} < $$min_frg)) {
			$allflag += 1;
			next FILE;
		}
		
		#OK, extract sequences now
		my $name = $1 if ($RMout =~ /^(.*?)\.split/);
		my $faloc = $name.".Extract";
		`mkdir $faloc` unless (-e $faloc);
		extract_sequences($dir,$faloc,$posi,$flank,$rc,$min_len,$extract,$frgs_all,$v);		
		push(@{$posi_done},$posi);
		
		#done for this RMoutput file
		print STDERR "     ...$RMout done\n" if ($v);
		
	}
	return ($allflag);
}

#----------------------------------------------------------------------------
# Get position files to extract repeats
# ($posi_list,$posi) = get_posi_from_RMout($RMout,$frgs_all,$frgs_nr,$posi_list,$v);
#----------------------------------------------------------------------------
sub get_posi_from_RMout {
	my ($RMout,$frgs_all,$frgs_nr,$posi_list,$v) = @_; 
	print STDERR "        ..in progress: get positions from $RMout\n" if ($v);
		
	#put RMout in array to be able to loop on it (because need to access to previous and next $block)
	my @RMout = get_RMout_array($RMout);
	
	#now loop on array
	my %blockcheck = ();
	my %c = ();
	my $posi = $RMout.".tab";
	push(@{$posi_list},$posi);

	open (my $fh, ">>","$posi") or confess "\n      ERROR (Sub get_posi_from_RMout): Failed to create/open $posi $!\n\n";
	#Note that all filters were previously done in the split step	
	for (my $i = 0; $i <= $#RMout; $i++){
		my ($div,$Gname,$Gstart,$Gend,$strand,$Rname,$classfam,$block) 
		   = ($RMout[$i]->[1],$RMout[$i]->[4],$RMout[$i]->[5],$RMout[$i]->[6],$RMout[$i]->[8],$RMout[$i]->[9],$RMout[$i]->[10],$RMout[$i]->[14]);		
 		my ($Rclass,$Rfam,$Rclassfam) = get_Rclass_Rfam($classfam);		
 		
		#count all frgs
		my $posiname = filename($posi);
		($frgs_all->{$posiname})?($frgs_all->{$posiname}+=1):($frgs_all->{$posiname}=1);
		
		#count number of fragments, corrected by blockID. Shared, so can't be complex data...
		my $fullblock = $Gname."_".$block;
		unless ($blockcheck{$Rname.$fullblock}) {
			($frgs_nr->{$posiname})?($frgs_nr->{$posiname}+=1):($frgs_nr->{$posiname}=1);
		}	
		$blockcheck{$Rname.$fullblock}=1; #keep in mind so that counts of fragments are corrected
		
		#now print coordinates for each copy - each file will be accessed by only one thread so no need to use Semaphore	
		my $cc = $Rname."_".$block."_".$Gname;
		($c{$cc})?($c{$cc}++):($c{$cc}=1);
		print $fh "$Gname\t$Gstart\t$Gend\t$strand\t$Rclass\t$Rfam\t$Rclassfam\t$div\t$block\t$c{$cc}\n";
	}	
	close($fh);
	print STDERR "        ..done: get positions from $RMout\n" if ($v);
	return ($posi_list,$posi);
}

#----------------------------------------------------------------------------
# load RM output in array of array
# my @RMout = get_RMout_array($RMout);
#----------------------------------------------------------------------------
sub get_RMout_array {
	my $rm = shift;
	my @array = ();
	open my $rm_fh, "<", $rm or confess "\nERROR (Sub get_RMout_array): can't open $rm $!\n";
	LINE: while(<$rm_fh>) {
		chomp (my $line = $_);
		next LINE if (($line =~ /position|[sS]core|Gname/) || ($line !~ /\w/)); #skip headers and white lines
		$line =~ s/^\s+//; #remove spaces in beginning of lines
		my @line = split(/\s+/,$line);
		push(@array, \@line);
	}	
	return @array;
}

#----------------------------------------------------------------------------
# Get Rclassfam from RMout
# my ($Rclass,$Rfam,$Rclassfam) = get_Rclass_Rfam($classfam);
#----------------------------------------------------------------------------
sub get_Rclass_Rfam {
	my($classfam) = shift;
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
	return ($Rclass,$Rfam,$Rclassfam);
}

#----------------------------------------------------------------------------
# Sub to loop on big hash to extract sequences
# extract_sequences($dir,$faloc,$posi,$flank,$rc,$min_len,$extract,$frgs_all,$v);	
#----------------------------------------------------------------------------
sub extract_sequences {
	my ($dir,$faloc,$posi,$flank,$rc,$min_len,$extract,$frgs_all,$v) = @_;
	print STDERR "        ..in progress: extract sequences from $posi\n" if ($v);

	#get genome (existence already checked)
	my $sp = filename($posi);
	$sp = $1 if $sp =~ /^(.*)\.out/;
	my $fa = $sp.".fa";
	$fa = $$dir."/".$fa;
			
	# connect to the fasta file
	my $db = Bio::DB::Fasta->new($fa, -reindex => 0) or confess "\n      ERROR (Sub extract_sequences): Failed to open Bio::DB::Fasta object from $fa $!\n\n";
		
	#get values from $extract if relevant
	my ($e_type,$e_nb) = split(",",$$extract);
	my $div_nb = int($e_nb /100 * 5) unless ($e_nb eq "na");
	
	#read posifile to extract sequences - one fasta file per element
	# NB: my $posi = $RMout."##".$Rname.".tab"; # $RMout = RMout first provided; otherwise, with $RMout from thread it's $RMout.tab
	my $posiname = filename($posi);
	my ($RMout,$Rname) = ($1,$2) if $posiname =~ /^(.+?\.out.*?)##(.+?)\.tab$/;
	#deal with output fasta file with extracted sequences - 1 file per repeat for each RMout/genome => all threads may want to print in it
	my $extracted = $faloc."/$Rname.extract.fa";
	my $subset = $faloc."/$Rname.extract.$e_nb.fa" unless ($e_nb eq "na");
	
	#Now deal with extraction
	my %check = ();
	my %prevcoords = ();
	my @seqs = ();
	my $s = 0;
	my $counts = 0;
	my @maxdiv = ();
	my @lowdiv = ();
	my $divc = 0;
	open (my $posi_fh, "<", $posi) or confess "\nERROR (Sub extract_sequences): Failed to open $posi $!\n";
	REP: while (<$posi_fh>) {
		chomp (my $l = $_);
		my ($Gname,$Gstart,$Gend,$strand,$Rclass,$Rfam,$Rclassfam,$div,$block,$c) = split ("\t",$l);
		$counts++;		
		#If first copy of any given block is processed that means print what is stored + store
		if ($c == 1) {
			#extract previous unless it is first line and nothing is stored; temp values because $seqs and $div may not be defined.
			my ($seqs,$divtemp) = extract_one_seq(\%prevcoords,$sp,$Rname,$db,$flank,$rc,$min_len,$v) if ($prevcoords{1});
			$seqs[$s]=$seqs if ($seqs);
			#If relevant, deal with the max and low %div sequences (eg if sequence went through filter > 80nt)
			unless (($e_nb eq "na") || (!$divtemp)) {
				push(@maxdiv,$divtemp);
				push(@lowdiv,$divtemp);
				unless ($divc < $div_nb ) {
					#sort the @alldiv array by the values of div, eg ->[0] ($newID is in [1])
					@maxdiv = sort { ($b -> [0] <=> $a -> [0]) } @maxdiv; #sort descending
					@maxdiv = splice(@maxdiv,0,$e_nb); #splice to keep first $e_nb values
					@lowdiv = sort { ($a -> [0] <=> $b -> [0]) } @lowdiv; #sort ascending	
					@lowdiv = splice(@lowdiv,0,$e_nb); #splice to keep first $e_nb values
				}
				$divc++;	
			}	
			#extraction done so increment the array coordinate in $seqs
			$s++ if ($seqs);
			#empty hashes and store this current line as the new previous
			%prevcoords = ();
		} 
		#Now store lines, in hash that has been reset or not (eg $c = 1 others)
		$prevcoords{$c} = "$Gname\t$Gstart\t$Gend\t$strand\t$Rclassfam\t$div";
		
		#Check for last line whatever the number; Ignore $div for the last one, makes it simpler and really no big deal.
		my ($seqs,$divtemp) = extract_one_seq(\%prevcoords,$sp,$Rname,$db,$flank,$rc,$min_len,$v) if ($counts == $frgs_all->{$posiname});
		$seqs[$s]=$seqs if ($seqs);
		
		#print $seqs if more than 500 to avoid big memory usage, unless specified not to
		if (($#seqs > 499) && ($e_type ne "only")) {
			open (my $fh, ">>","$extracted") or confess "\nERROR (Sub extract_sequences): Failed to create/open $extracted $!\n";		
			for (my $i = 0; $i <= $#seqs; $i++) {
				print $fh ">".$seqs[$i]."\n" if ($seqs[$i]);
			} 
			close($fh);
			#reinitialize
			@seqs = ();
			$s = 0;
		}
	}
	close($posi_fh);
	
	#now print whatever is left in @seqs and was not printed, unless specified not to
	my $seqsnb = @seqs;
	if (($seqsnb > 0) && ($e_type ne "only")) {
		open (my $fh, ">>","$extracted") or confess "\nERROR (Sub extract_sequences): Failed to create/open $extracted $!\n";		
		for (my $i = 0; $i <= $#seqs; $i++) {
			print $fh ">".$seqs[$i]."\n" if ($seqs[$i]);
		} 
		close($fh);
	}
	
	print STDERR "        ..done: extract all sequences from $posi\n" if ($v);

	#Now print a subset if relevant, from file with all sequences (unless it's not there obviously).
	unless (($e_nb eq "na") || (! -e $extracted)) {
		print STDERR "        ..in progress: extract a subset of $e_nb sequences from $posi\n" if (($v) && ($e_type ne "compl"));
		print STDERR "        ..in progress: extract a subset of $e_nb sequences + complementary set from $posi\n" if (($v) && ($e_type eq "compl"));
		my $reind;
		my $ind = "$extracted.index";
		(-e $ind)?($reind = 0):($reind = 1);	
		my $extract_db = Bio::DB::Fasta->new($extracted, -reindex=>$reind, -makeid=>\&make_my_id_m) or confess "\nERROR (Sub extract_sequences): Failed to open Bio::DB::Fasta object from $extracted $!\n";
		my @extr_ids = $extract_db->ids();
		my $nb_extr_ids = @extr_ids;	
		#1. if nb asked for is >= nb of sequences for this repeat, then nothing to do, just cp file
		if ($nb_extr_ids <= $e_nb) {
			my $subset = $1.".$e_nb.fa" if $extracted =~ /^(.*)\.fa/;
			`cp $extracted $subset`;
			print STDERR "        ..done: extract a subset of $e_nb sequences from $posi [same file as $extracted, not enough sequences] \n" if (($v) && ($e_type ne "compl"));
			print STDERR "        ..done: extract a subset of $e_nb sequences + complementary set from $posi [no complementary set: no enough sequences]\n" if (($v) && ($e_type eq "compl"));
		} else {
			extract_random($extract_db,$e_type,$e_nb,$extracted,$div_nb,\@maxdiv,\@lowdiv);
			print STDERR "        ..done: extract a subset of $e_nb sequences from $posi\n" if (($v) && ($e_type ne "compl"));
			print STDERR "        ..done: extract a subset of $e_nb sequences + complementary set from $posi\n" if (($v) && ($e_type eq "compl"));
		}	
	}
	return;	
}			
			
#----------------------------------------------------------------------------
# Extract a sequence using coordinates 
# $seqs = extract_one_seq(\%prevcoords,$sp,$Rname,$db,$flank,$rc,$min_len,$v);
#----------------------------------------------------------------------------
sub extract_one_seq {					
	my ($coords,$sp,$Rname,$db,$flank,$rc,$min_len,$v) = @_;
	
	#get list of the ID of the genome file
	my @dbIDs = $db->get_all_ids();
	
	my %check = ();
	my $Gs;
	my $Ge;
	my $len;
	my $nb = 0;	
	my $div_corr;
	my $pondiv;
	my @storediv = ();
	my $str;
	#First round to get start and end, if more than one fragment per block, for the name of the sequence
	foreach my $c (sort keys %{ $coords }) {
		my ($Gname,$Gstart,$Gend,$strand,$Rclassfam,$div) = split ("\t",$coords->{$c});
		my $currlen = $Gend - $Gstart + 1;
		#Correct coords for new name + get real length
		$div_corr = $div;
		if ($c == 1) {
			$Gs = $Gstart;
			$Ge = $Gend; #to intialize
			$len = $currlen;
			$pondiv = $div*$currlen;
		} else {
			$Ge = $Gend; #to correct - Gs remains the same
			$len+=$currlen;
			$pondiv+=$div*$currlen;
		}
		$nb++;
	}
	#correct %div if $nb >1
	$div_corr = $pondiv/$len if ($nb > 1);
	
	#no extraction if too small
	if ($$min_len ne "na") {
		return if ($len < $$min_len);
	}
	
	#Now loop to extract sequence(s)
	my $newId;
	my $descr;
	my $seq = "";
	my $fl5;
	my $fl3;
	foreach my $c (sort keys %{ $coords }) {
		my ($Gname,$Gstart,$Gend,$strand,$Rclassfam,$div) = split ("\t",$coords->{$c});

		#name of the sequence extracted will be only gb when there is gi and gb in name (otherwise too long)
		my $Gnamemod = $Gname;
		$Gnamemod =~ s/^gi\|.*\|gb/gb/;
		
		#prep new sequence name
		$descr = $Rclassfam."__frg=".$nb."__len=".$len."_st=".$strand."_div=".$div_corr."_sp=".$sp;
		$newId = $Rname."__".$nb."_".$Gnamemod."--".$Gs."-".$Ge;
		
		#Check if $Gname is in $db and get $Glen at the same time or exit the sub
		my $Glen;
		if ($db->length($Gname)) {
			$Glen = $db->length($Gname); 
		} else {
			print STDERR "        ERROR (Sub extract_one_seq): $Gname was not found in genome file\n" if ((! $check{$Gname}) && ($v));
			$check{$Gname}=1;
			return;
		}	
						
		#get coords with flanks if relevant, but only if outside fragments 
		unless ($$flank eq "na") {
			if ($c==1) { #first fragment => get flank at start
				$Gstart = $Gstart-$$flank;
				$Gstart = 1 if ($Gstart < 1); #correct to avoid negative start
			} 
			if ($c == $nb) { #last fragment => get flank at end
				$Gend = $Gend+$$flank;
				my $Glen = $db->length($Gname); 
				$Gend = $Glen if ($Gend > $Glen); #correct to avoid getting out of scaffold
			}
		}			
		# now extract target sequence, with or without flankings
		my $subSeq = $db->seq($Gname,$Gstart,$Gend);
		if ($c == 1){
			$seq = $subSeq; #get the sequence of the first fragment
		} else {
			$seq = $seq.$subSeq; # add other piece(s)
		}
		
		# if relevant, get actual flankings
		unless ($$flank eq "na") {
			if ($c == 1) {
				$fl5 = $Gs - $Gstart; #real TE start minus start of flanking => real 5' flanking size
				$fl3 = $Gend - $Ge; #Do to avoid errors Use of uninitialized value $fl3 - will be corrected after
			}
			if ($c == $nb) {
				$fl3 = $Gend - $Ge; #new end after fl determination minus actual end => real 3' flanking size
			}
		}
		$str = $strand;
	}
	
	#now get RC if needed
	my $cseq = $seq;
	if (($$rc eq "y") && ($str ne "+")) {
		my $revseq = revcom($cseq);
		$cseq = $revseq->seq();
	}
	#correct desc if relevant to add flanking length
	$descr = $descr."_5fl=".$fl5."_3fl=".$fl3 unless ($$flank eq "na");	
	#get full sequence now and store it based on %div
	my $fullseq = "$newId\t$descr\n$cseq";
	($storediv[0],$storediv[1]) = ($div_corr,$fullseq);
	#now return these values
	return ($fullseq,\@storediv);	
}

#----------------------------------------------------------------------------
# Modified version of make_my_id for Bio::DB::Fasta
# Because I want to keep description here.
#----------------------------------------------------------------------------
sub make_my_id_m {
	my $line = shift;
	$line =~ /^(.*)$/; #keep description
	return $1;
}

#----------------------------------------------------------------------------
# Extract random set of sequences using a Bio::DB::Fasta object
# extract_random($extract_db,$e_type,$e_nb,$extracted,$div_nb,\@maxdiv,\@lowdiv);
#----------------------------------------------------------------------------
sub extract_random {
	my ($db,$type,$nb,$file,$div_nb,$maxdiv,$lowdiv) = @_;
	my @ids = $db->ids();
	my %ids = ();
	
	#Extract the X seqs, with the maxdiv and lowdiv
	my %already = (); #to keep in mind which ones are already picked, so they are not printed again if randomly selected
	my $subset = $1.".$nb.fa" if $file =~ /^(.*)\.fa/;
	open (my $fh, ">","$subset") or confess "\nERROR (Sub extract_random): Failed to create/open $subset $!\n";
	#first print the low and high div
	for (my $i = 0; $i <= $div_nb; $i++) {
		print $fh ">".$maxdiv->[$i]->[1]."\n";
		print $fh ">".$lowdiv->[$i]->[1]."\n";
		$already{$maxdiv->[$i]} = 1;
		$already{$lowdiv->[$i]} = 1;
		$ids{$maxdiv->[$i]}=1;
		$ids{$lowdiv->[$i]}=1;
	} 
	#Randomize some other seqs now, using ids
	# shuffle array
	my @shuffled = shuffle(@ids);
	# keep only a slice of the array => $nb values of this array
	my @slice = @shuffled[ 0 .. $nb-1 ]; #even if ALL max and low div are in the slice, there will be enough ids to go through
	my $slicenb = @slice;
	# extract the subset of sequences
	my $i = 0;
	RAND: foreach my $id (@slice) {
		my $seq = $db->seq($id);
		if  (! $seq) {
			print "\n     ERROR (sub extract_random): $id not found in $file\n" 
		} else {
			print $fh "$id\n$seq\n" unless ($already{$id});
			$ids{$id}=1;
			$i++;
		}
		last RAND if ($i == $nb); #exit extraction loop if enough sequences
	}
	close $fh;

	#get complementary sequences if relevant
	if ($type eq "compl") {
		my $outc = $1.".".$nb."_compl.fa" if $file =~ /^(.*)\.fa/;
		open (my $outcfh, ">","$outc") or die "\n     ERROR (sub extract_random): Failed to create file $outc $!\n";
		foreach my $id (@ids) {
			my $seq = $db->seq($id);
			if  (! $seq) {
				print "\n     ERROR (sub extract_random): $id not found in $file\n" 
			} elsif (! $ids{$id}) {
				print $outcfh "$id\n$seq\n";
			}
		}
		close $outcfh;
	}
	return;
}

#----------------------------------------------------------------------------
# Cat the outputs
# cat_fa($RMouts,$dir,$v);
#----------------------------------------------------------------------------
sub cat_fa {
	my ($list,$dir,$v) = @_;
	print STDERR "\n --- Concatenating extracted fasta sequences (for each Repeat Masker output)\n" if ($v);
	FILE: foreach my $folder (@$list) {
		chomp $folder;
		#next FILE unless (-e $folder); #for some reason there was an extra empty value in the list. Better check on file existence.
		print "     cat $dir/$folder/*extract.fa > $dir/$folder.extracted.fa\n";
		`cat $dir/$folder/*extract.fa > $dir/$folder.extracted.fa`;
	}
	return;
}	








		
		
		