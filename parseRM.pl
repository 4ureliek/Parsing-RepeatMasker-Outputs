#!/usr/bin/perl -w
##########################################################
# Author  :  Aurelie Kapusta
# version :  see below + see changelog
# email   :  4urelie.k@gmail.com
# PURPOSE :  parsing Repeat Masker outputs, .out and/or .align
##########################################################
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Bio::SeqIO;

my $version = "4.0";

my $changelog = "
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
#            Fixed the -parse option; changes in storing of %div, %del and %ins to remove redundant info and make it run faster
#            Count fragments that are start to end of consensus, +/- 2
#            Count fragments that are > 50% and > 90 % of the consensus (different than above)
#            Get median of fragments length
#	- v4.1 = 09 Sep 2016
#            Fixed bug that would crash when -p is set without -f

\n";
	
my $usage = "\nUsage [v$version]: 
	perl parseRM.pl --in <genome.out/align>
	               [--parse] [--fa <maskedfile.fa>] [--nrem] [-glen <X>] [--rlib <repeat_library.fa>]
	               [--age <file.txt> OR <X,X>] [--te] [--my <file.txt> OR <X,X>] [--land <max,bin>] 
	               [--dir] [--kim] [-simple <type>] [--edit <file>] [--what <type,name>] [--contain] 
	               [--v] [--updates] [--help]
	
	This script will process a RepeatMasker output .out or .align file (or several, see -dir option),
	with 3 non exclusive parsings (they can all be set):
	-p => to get a summary of the masking, as well as amount or DNA, counts of fragments, + several other details
          for each repeat name (all-repeats file), family, class and total amount (summary file)
	-a => to determine the amounts of DNA in a genome that is masked by repeats of different lineages / %divergence categories
	-l => to split the amount of DNA by bins of %div or My, allowing to generate landscape graphs
          for each repeat name, family or class (one output for each)
	
    Note: if all 3 options are set (-age, -TEage and -land) are chosen, there will be an additional output 
          with bins by %div or My, but split by age cagtegories specified in the -age file instead of repeat name/family/class            
	
    MANDATORY ARGUMENTS:	
     -i,--in      => (STRING) 
                     RepeatMasker output .out or RepeatMasker output .align 
                     (use -dir and provide a directory path here to parse both)
                     Use of .align is best for -l, because he %div corrected for higher mutation rate at CpG sites 
                     and with Kimura 2-Parameter divergence metric is in the .align files (graphs are smoother)
                     However, because the script will treat separately repeats that overlap 
                     (which are merged in the .out, with one name kept), meaning that for the -p parsing, 
                     the amount of DNA masked by these repeats will be much higher than if the .out is parsed

  
    OPTIONAL ARGUMENTS RELATED TO --parse
     -p,--parse   => (BOOL)
                     To get the summary of Repeat Masker such as number of counts, etc. 
                     Following options are only relevant if this is chosen
     -f,--fa      => (BOOL)
                     If the file THAT WAS MASKED, typically a genome, is in the same directory as the .out or .align 
                     If not provided or not found, % of genome masked won't be calculated
                     Names should correspond before .align and .fa(sta), and/or .out and .fa(sta)
     -n,--nrem    => (BOOL)
                     To remove Ns from the genome file before getting its length
     -g,--glen    => (INT or STRING)
                     Alternatively to genome file(s), you can provide the value, in nt, that will be used as the 100%
                     If several genomes are being looked at (--dir chosen for example) this can be a file,
                     containing 2 columns: filename \\t X                           
                        with filename = file name of the file to be parsed that is in the directory --in
                        and X = the value
     -r,--rlib    => (STRING)
                     Library of consensus sequences used to mask, to have their lengths included in the output 
                     If several maskings are parsed, you can concatenate the different libraries
                     (same if RM was run with -lib, concat RM library and the user one)   
                     If a repeat could not be found or if -r is not chosen the value in the output will be \"nd\".              
          
    OPTIONAL ARGUMENTS RELATED TO --age    
     -a,--age     => (STRING)
                     Option1: load a file (set --te so the file is leaded properly)
                        File should contain TE age data, tabular delimited. Columns should be:
                        Rname \\t Rclass \\t Rfam \\t Lineage
                        (with Rname an exact match to the repeat names in repeat masker output; class and family can be changed)
                        If -d is chosen, the file must contain info for ALL repeats of all files
                     Option2: 1 or 2 values (INT or FLOAT), separated by a comma if 2 values: 13.5,15 (or 15,13.5)
                        any DNA masked by a TE with %div < the smallest number given (here 13.5) will be placed as lineage specific 
                        any DNA masked by a TE with %div > the largest number given (here 15) will be placed as ancient
                        the rest will be placed as \"nd\". 
                     If -m is set, then numbers here HAVE TO BE in My instead of %divergence
                     If -d is chosen and different values should be used for each input files
                        you can provide here a file containing 2 columns: filename \\t X                           
                        with filename = input file name that is in the directory -in
                        and X = what would be given as Option2 for each file
     -t,--te      => (BOOL)
                     Set this if --age is a file with TE age info
     -m,--my      => (INT or STRING)
                     To output age data in My and not %div (set it only one time if noth -age and -land are chosen)
                     Substitution rates (per My) should be provided here:
                     If -d not chosen, use -m substitution_rate (ex. -m 0.00199)
                     If -d is chosen, then different values can be provided for each input files,
                        if you use -m subst-rates.tab, with subst-rates.tab = a file with 2 columns: 
                        filename_as_in_dir \\t substitution_rate
                         
    OPTIONAL ARGUMENTS RELATED TO --land
     -l,--land    => (STRING)
                     To generate additional outputs that can be used to make landscape graphs (by repeat, by family and by class)
                     2 values should be given, separated by a comma: max,bin
                        with max being the value of the last bin, where everything higher than this %div or My will be put
                        and bin = size of the bins
                     If -m is set, then numbers here HAVE TO BE in My instead of %divergence
                     Typically: -land 50,1 (if %div are corrected for CpGs then values tend to be higher)
                     Or: -land 5,0.25 if you wish to look only at recent stuff in a very dynamic genome
     -m,--my      => (INT or STRING)
                     To output bins in My and not %div: substitution rates need to be provided, see above
                         
    OTHER OPTIONAL ARGUMENTS
     -d,--dir     => (BOOL)
                     Chose this if --in is in fact a directory containing a bunch of .out or .align files                       
                     Note that the directory should contain no other file with .out or .align or, 
                        if genomes are provided, end by .fa (genomes); if so, names should be identical 
                        before .align and .fa(sta), and/or .out and .fa(sta)
                     No technical objection to have a mixture of .out and .align in the folder 
                     (how to read them is based on having .out or .align in their names)
     -k,--kim     => (BOOL)
                     If no .align file available, use this to correct %div with Kimura equation: %div = -300/4*log(1-%div*4/300)
                     Do not chose this if you are parsing any .align files! 
     -s,--simple  => (STRING)
                     If you wish to KEEP (some) nonTE sequences (as far as annotation goes of course)
                     By default, all non TE sequences are filtered out.
                     Chose between all or nonTE
                        all   = KEEP ALL nonTE stuff (Low_complexity, Simple_repeat, Satellites, etc)
                        nonTE = KEEP only non TE sequences when class = nonTE (other stuff will be filtered out)                   
     -t,--te    => (STRING)
                     File with TE information as follow: Rname \\t Rclass \\t Rfam
                     If you need to correct some class and/or families of some repeats.
                     It may complement the -age file. Not all repeats need to be in it, just the ones to correct.
                     The match will be done on Rname, so it has to be the same as in the repeat masker output;
                     An easy way to get it is to run this script a first time and use the output all-repeats.tab
                     However, this can be proven difficult when parsing .align files, with the whole LINE1-3end or LINE_5end renamings.                          
     -w,--what    => (STRING)
                     Typically: -w <type,name> to run the script on only a subset of repeats. Not case sensitive.
                     The type can be: name, class or family and it will be EXACT MATCH unless -c is chosen as well
                     ex: -w name,nhAT1_ML => only fragments corresponding to the repeat named exactly nhAT1_ML will be looked at
                         -w class,DNA => all repeats with class named exactly DNA (as in ...#DNA/hAT or ...#DNA/Tc1)
                         -w family,hAT => all repeats with family named exactly hAT (so NOT ...#DNA/hAT-Charlie for example)
     -c,--contain => (BOOL) 
                     To check if the \"name\" determined with -filter is included in the value in Repeat Masker output, instead of exact match
                     ex: name,HERVK => all fragments containing HERVK in their name
                         family,hAT => all repeats with family containing hAT (...#DNA/hAT, ...#DNA/hAT-Charlie, etc)
     -v,--v       => (BOOL)
                     Verbose mode, make the script talks to you
                     Print the version if only option
     -u,--updates => (BOOL)
                     Print change log (updates)
     -h,--help    => (BOOL)
                     Print this usage\n\n";      


################################################################################
# Get arguments/options, check some of them, print details of the run if verbose chosen
################################################################################
my $nonTE = "no";
my ($in,$p,$aa,$l,$gfile,$glen,$Nrem,$lib,$My,$filter,$contain,$TEage,$TEs,$dir,$k,$help,$v,$chlog);
GetOptions ('in=s'     => \$in, 
            'dir'      => \$dir, 
            'parse'    => \$p, 
            'age=s'    => \$aa, 
            'land=s'   => \$l, 
            'my=s'     => \$My,
            'te'       => \$TEage,  
            'fa'       => \$gfile, 
            'nrem'     => \$Nrem, 
            'glen=s'   => \$glen, 
            'rlib=s'   => \$lib, 
            'kim'      => \$k, 
            'simple=s' => \$nonTE, 
            'edit=s'   => \$TEs, 
            'what=s'   => \$filter, 
            'contain'  => \$contain, 
            'updates'  => \$chlog, 
            'help'     => \$help, 
            'v'        => \$v);

#check step to see if mandatory argument is provided + if help/changelog
die "\n Script parseRM.pl version $version\n\n" if ((! $in) && (! $p) && (! $aa) && (! $l) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ($help);
die $usage if ((! $in) && (! $p) && (! $aa) && (! $l));
die "\n please chose one of --parse, --age or --land (use -h to see usage)\n\n" if ((! $p) && (! $aa) && (! $l));

#avoid / at the end of path in case it's a directory
$in = $1 if (($dir) && ($in =~ /^(.*)\/$/));

#"log"
($k)?($k = "y"):($k = "n");
select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately
print STDERR "\n--------------------------------------------------\n" if ($v);
print STDERR " --- Script parseRM.pl started (v$version)\n" if ($v);
print STDERR "      - Directory containing input files = $in (--dir chosen)\n" if (($dir) && ($v));
print STDERR "      - Input file = $in\n" if ((! $dir) && ($v));
die "\n $in does not exist?\n\n" if (($in) && (! -e $in));
print STDERR "      - All non TE repeats will be filtered out\n" if (($nonTE eq "no") && ($v));
print STDERR "      - Elements with class = nonTE won't be filtered out (-simple nonTE chosen)\n" if (($nonTE eq "nonTE") && ($v));
print STDERR "      - Non TE repeats won't be filtered out (-simple all chosen)\n" if (($nonTE eq "all") && ($v));
print STDERR "      - Repeats will be filtered based on $filter (--what)\n" if (($filter) && ($v));
print STDERR "        regexp will be used and not exact match (--contain)\n" if (($filter) && ($contain) && ($v));
print STDERR "      - %div will be corrected using -300/4*log(1-%div*4/300)\n" if (($k ne "n") && ($v));
print STDERR "      => global parsing will be done, with following options:\n" if (($p) && ($v));
print STDERR "         - genome file(s) will be looked for (--fa chosen)\n" if (($gfile) && ($v));
print STDERR "         - Ns will be removed from the genome before getting its length (--nrem chosen)\n" if (($Nrem) && ($v));
print STDERR "         - length(s) are given instead (--glen $glen)\n" if (($glen) && ($v));
print STDERR "         - library file = $lib\n" if (($lib) && ($v));
print STDERR "           WARN: $lib does not exist? Length of repeats won't be in the output\n" if (($lib) && (! -e $lib) && ($v));
undef $lib if (($lib) && (! -e $lib));
print STDERR "      => age split will be done, with following options:\n" if (($aa) && ($v));
print STDERR "         - Age info will be based on $aa\n" if (($aa) && ($v));
print STDERR "         - Age will be in My, using substitution rates ($My)\n" if (($aa) && ($My) && ($v));
print STDERR "      => landscape data will be generated, with following options:\n" if (($l) && ($v));
print STDERR "         - With max_bin,bin_len = $l\n" if (($l) && ($v));
print STDERR "         - Age will be in My, using substitution rates ($My)\n" if (($l) && ($My) && ($v));
print STDERR "--------------------------------------------------\n" if ($v);

($dir)?($dir = "y"):($dir = "n");
($Nrem)?($Nrem = "y"):($Nrem = "n");
($contain)?($contain = "y"):($contain = "n");
$My = "n" unless ($My);
$filter = "na" unless ($filter);

################################################################################
# MAIN
################################################################################
#Get list of input files if -dir, or just load the one from -in
my @files = ();
($dir eq "y")?(@files = `ls $in`):(push(@files,$in));
my $path = get_path($in);

### Deal with -parse stuff if relevant
print STDERR " --- Prepping steps for --parse...\n" if (($p) && ($v));
#genome lengths
my $genlen;
print STDERR "     - Loading sequence length(s) from fasta file...\n" if (($p) && ($glen) && ($v));
$genlen = load_val($in,$glen,"len") if ($glen);
print STDERR "     - Getting sequence length(s) from $gfile...\n" if (($p) && ($gfile) && ($v));
$genlen = get_length_noBio($in,$dir,$Nrem,$v) if ($gfile); #then should look for genomes

#lib lengths
print STDERR "     - Getting lengths of consensus sequences from $lib...\n" if (($p) && ($lib) && ($v));
my $liblen;
$liblen = get_lengths_noBio($lib,$v) if ($lib);

### Deal with -age stuff if relevant
print STDERR "\n --- Loading age info from $aa\n" if (($aa) && ($v));
my $age;
$age = load_val($in,$aa,"div") if (($aa) && (! $TEage));
$age = load_TEdb($aa) if (($aa) && ($TEage));
($TEage)?($TEage = "y"):($TEage = "n");

#prep cat output if relevant
my $fall = $in.".splitage_all.tab";
if (($aa) && ($dir eq "y")) {	
	`rm $fall` if (-e $fall);
	open my $fhall, ">>", $fall or confess "\nERROR (main): could not open to read $fall $!\n";
	print $fhall "#nr_masked = amount of masked nt to consider\n";
	print $fhall "#tot_masked = total amount of masked nt, including redundant maskins / overlaps\n";
	print $fhall "#double_masked = total amount of nt that were masked by at least 2 repeats = redundant masking / overlaps\n\n";	
	print $fhall "#Input_file\tnr_masked\ttot_masked\tdouble_masked\tAgeCat\tCounts_this_age\tnr_Counts_this_age\tnr_masked_this_age\t%nr_masked_this_age\n\n";	
	close $fhall; 
}

### Prep steps not specific to p, a or l
# Load substitution rates from $My if relevant
my $srates = ();
$srates = load_val($in,$My,"rates") unless ($My eq "n");
#Deal with $TEs if relevant
my $TE = ();
$TE = get_TEs_infos($TEs) if ($TEs);

$p = "n" unless ($p);
$aa = "n" unless ($aa);
$l = "n" unless ($l);

#load all the RM files now and loop
print STDERR "\n --- Now parsing RM output(s)\n" if ($v);
parseRM($p,$aa,$l,$filter,$contain,$genlen,$liblen,$age,$My,$TEage,$srates,$TE,$path,$in,\@files,$k,$dir,$nonTE,$v);

print STDERR "\n --- Script done\n" if ($v);
print STDERR "    -> age split files: *.agesplit.tab for all files: $fall\n" if (($aa ne "n") && ($v));
print STDERR "                        + $fall\n" if (($aa ne "n") && ($dir eq "y") && ($v));
print STDERR "    -> parsed files: *.parseRM.all-repeats.tab\n" if (($p ne "n") && ($v));
print STDERR "                     *.parseRM.summary.tab\n" if (($p ne "n") && ($v));
print STDERR "    -> landscape files: *.landscape.*.Rname.tab\n" if (($l ne "n") && ($v));
print STDERR "                        *.landscape.*.Rfam.tab\n" if (($l ne "n") && ($v));
print STDERR "                        *.landscape.*.Rclass.tab\n" if (($l ne "n") && ($v));
print STDERR "                        *.landscape.*.Rage.tab\n" if (($l ne "n") && ($TEage ne "n") && ($aa ne "n") && ($v));
print STDERR "\n" if ($v);
exit;

##########################################################################################################
# SUBROUTINES
# ##########################################################################################################
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
# my $path = get_path($filename);
#----------------------------------------------------------------------------
sub get_path {
	my($file) = shift;
	($file =~ /\//)?($file =~ s/(.*)\/.*$/$1/):($file = ".");
	return $file;
}

#----------------------------------------------------------------------------
# load age data and subsitution rates in array
# $genlen = load_val($in,$glen,"len") if ($glen);
# $age    = load_val($in,$aa,"div");
# $srates = load_val($in,$My,"rates");
#----------------------------------------------------------------------------
sub load_val {
	my ($in,$data,$type) = @_;
	my %hash = ();	
	if (-e $data) { #then it's a file
		open my $fh, "<", $data or confess "\nERROR (Sub load_val): could not open to read $data $!\n"; 
 		while (<$fh>) {
 			chomp (my $line = $_);
 			my ($f,$d) = split('\s+',$line);
 			my $fname = filename($f); #just in case user put the path instead of the name
 			if ($type eq "rates") {
				$hash{$fname}=$d;
			} else {
				my ($m,$n);
				($d =~ /,/)?(($m,$n) = split(",",$d)):(($m,$n)=($d,$d));
				my @v = (); #values
				($m < $n)?(@v=($m,$n)):(@v=($n,$m));				
				$hash{$fname}=\@v;
			}			
 		}
 		close ($fh);
	} else {
		#Then values should be the variable
		confess "\nERROR (Sub load_val): $data does not exist as a file but also not numerical?\n\n" if ($data !~ /[0-9\.,]+/);
		my $fname = filename($in);
		if ($type eq "rates") {
			$hash{$fname}=$data;
		} else {
			my ($m,$n);
			($data =~ /,/)?(($m,$n) = split(",",$data)):(($m,$n)=($data,$data));		
			my @v = (); 
			($m < $n)?(@v=($m,$n)):(@v=($n,$m));			
			$hash{$fname}=\@v; 
		}	
	}
	return(\%hash);
}

#----------------------------------------------------------------------------
# Get total length of all sequences
# This sub does not use BioPerl - avoid having to index the genome
# $genlen = get_length_noBio($in,$dir,$Nrem,$v) if ($gfile); 
#----------------------------------------------------------------------------
sub get_length_noBio {
	my ($in,$dir,$Nrem,$v) = @_;	
	#Could be several files to deal with
	my @list;	
	if ($dir eq "n") {
		my $fa = $in;
		$fa = $1 if (($in =~ /(.*)\.out/) || ($in =~ /(.*)\.align/));
		$fa =~ $1 if ($fa =~ /(.*\.fa)sta$/);
		if (! -e $fa) {
			print STDERR "       WARN: $fa could not be found, % genome won't be determined\n";
			return;
		} else {
			push(@list,$fa);
		}	
	} else {
		@list = `ls $in`; #do it on all fasta files located in $in basically
	}
	
	my %len = ();	
	FA: foreach my $fa (@list) {
		chomp $fa;		
		next FA unless (($fa =~ /\.fa$/) || ($fa =~ /\.fasta$/) || ($fa =~ /\.fsa$/));
		print STDERR "        -> $fa\n" if ($v);	
		$fa = Nrem($fa,$v) if ($Nrem eq "y");	
		my $do = 0;
		my $totalfile = "$fa.length";	
		if (-e $totalfile) {
			print STDERR "           total length has been previously calculated ($totalfile exists)\n" if ($v);
			open (my $length_fh, "<", $totalfile) or ((warn "           WARN: could not open to read $totalfile, but length will be recalculated ($!)\n") && ($do++));
			unless ($do == 1) {
				while (<$length_fh>) {
					$len{filename($in)}=$_;
				}
				close ($length_fh);
				print STDERR "             => ".$len{filename($in)}." nt\n" if ($v);
				return(\%len);
			}
		}
		if ((! -e $totalfile) || ($do == 1)) {
			print STDERR "           obtaining total length\n" if ($v);
			#looping through fasta file
			my $id = "";
			my $l = 0;
			my $c = 0;			
			open (my $fa_fh, "<", $fa) or warn "           WARN: could not open to read $fa, % genomes won't be calculated ($!)\n";
			while (<$fa_fh>) {
				chomp (my $line = $_);
				if (substr($line,0,1) eq ">") {
					#first get and print unless first header
					unless ($c == 0) {
						$len{filename($in)}+=$l;
					}
					$c=1;
					#store header and reinitialize length
					my @id = split (/\s+/,$line);
					$id = $id[0];
					$id =~ s/>//;
					$l = 0;
				} else {
					#get length; could be more than one line so increment
					$l+=length($line);
				}
			}
			#get and print len last sequence
			$len{filename($in)}+=$l;
			close ($fa_fh);
			open (my $fh, ">", $totalfile) or warn "           WARN: could not open to write $totalfile $!\n";
			print $fh "$len{filename($in)}";
			close $fh;
			print STDERR "           => ".$len{filename($in)}." nt\n" if ($v);
		}
	}	
	return (\%len);
}

#----------------------------------------------------------------------------
# Get lengths of all sequences and store that by file and sequence ID. Note that if some are not unique, it just replaces by last length.
# Specifically for TEs.
# This sub does not use BioPerl - avoid having to index the file
# $liblen = get_lengths_noBio($lib,$v) if ($lib);
#----------------------------------------------------------------------------
sub get_lengths_noBio {
	my ($fa,$v) = @_;	
	print STDERR "        -> $fa\n" if ($v);
	unless (($fa =~ /\.fa$/) || ($fa =~ /\.fasta$/) || ($fa =~ /\.fsa$/)) {
		print STDERR "       WARN: $fa not a fasta file? => skipped\n" if ($v);
		return;
	}	
	if (! -e $fa) {
		print STDERR "       WARN: $fa could not be found, % genome won't be determined\n";
		return;
	}
	my $lengthfile = "$fa.lengths";
	my $do = 0;
	my %len = ();
	if (-e $lengthfile) {
		print STDERR "           lengths have been previously calculated ($lengthfile exists) => extracting\n" if ($v);
		#extract lengths now
		open (my $lengths_fh, "<", $lengthfile) or ((warn "           WARN: could not open to read $lengthfile, but lengths will be recalculated ($!)\n") && ($do++));
		unless ($do == 1) {
		while (<$lengths_fh>) {
				chomp (my $line = $_);
				my ($id,$l) = split(/\s+/,$line);
				my ($Rname,$Rclassfam) = split('#',$id);
				$len{lc($Rname)}=$l;
			}	
			close ($lengths_fh);
			return(\%len);
		}
	}
	if ((! -e $lengthfile) || ($do == 1)) {
		print STDERR "           obtaining lengths\n" if ($v);
		my $id = "";
		my $l = 0;
		my $c = 0;
		open (my $fa_fh, "<", $fa) or warn "           WARN: could not open to read $fa, repeat lengths won't be in the output $!\n";
		open (my $len_fh, ">", $lengthfile) or warn "           WARN: could not open to write $lengthfile, but lengths will be calculated $!\n";
		while (<$fa_fh>) {
			chomp (my $line = $_);
			if (substr($line,0,1) eq ">") {
				#first get and print unless first header
				unless ($c == 0) {
					print $len_fh "$id\t$l\n";
					my ($Rname,$Rclassfam) = split('#',$id);
					$len{lc($Rname)}=$l;
				}
				$c=1;
				#store header and reinitialize length
				my @id = split (/\s+/,$line);
				$id = $id[0];
				$id =~ s/>//;
				$l = 0;
			} else {
				#get length; could be more than one line so increment
				$l+=length($line);
			}
		}
		#get and print len last sequence
		print $len_fh "$id\t$l\n";
		my ($Rname,$Rclassfam) = split('#',$id);
		$len{lc($Rname)}=$l;
		close ($fa_fh);
		close ($len_fh);
		print STDERR "           lengths are now extracted in $lengthfile\n" if ($v);
	}	
	return (\%len);
}

#----------------------------------------------------------------------------
# remove Ns in a genome file
# $fa = Nrem($fa,$v) if ($Nrem eq "y");
#----------------------------------------------------------------------------
sub Nrem {
	my ($fa,$v) = @_;
	my $genometmp = "$fa.Nrem";
	my $genome = "$fa.Nrem.fa";
	
	if (-e $genome) {
		print STDERR "           Ns already removed from $fa ($genome exists) - skipping N removal step\n";
	} else {
		#remove Ns
		print STDERR "           Ns may be not removed from $fa ($genome does not exists) - Removing N...\n";
		my $Gentmp = Bio::SeqIO->new(-file => $fa, -format => "fasta") or die "Failed to create Bio::SeqIO object from $fa $!\n";
		open(my $GminN_fh, ">", $genometmp) or confess "\nERROR (sub Nrem): could not open to write $genometmp $!\n";	
		while( my $seq = $Gentmp->next_seq() ) {
			my $sequence = $seq->seq;
			$sequence =~ s/[Nn]//g;
			my $id = $seq->display_id."\t".$seq->desc;
			print $GminN_fh ">$id\n$sequence\n";
		}
		# Rewrite sequences in fasta format just to be sure
		my $Gen = Bio::SeqIO->new(-file => $genometmp, -format => "fasta") or confess "\nERROR (sub Nrem): could not create Bio::SeqIO object from $genometmp $!\n";
		my $GenMinusN = Bio::SeqIO->new(-file => ">$genome", -format => "fasta") or confess "\nERROR (sub Nrem): could not open to write $genome $!\n";	
		while( my $seq2 = $Gen->next_seq() ) {
			$GenMinusN->write_seq($seq2);		
		}
	}
	unlink($genometmp);
	return $genome;
}

#----------------------------------------------------------------------------
# load TE data if it's a lineage file
# $age = load_TEdb($aa) if (($aa) && ($TEage));
#----------------------------------------------------------------------------
sub load_TEdb {
	my $input = shift; #file name
	my %TEs = ();
	open(my $fh, "<", $input) or confess print STDERR "ERROR (sub get_TEs_infos): could not open $input $!\n";
	LINE: while(<$fh>) {
		chomp (my $line = $_);
		next LINE if ($line !~ /\w/);
		my @TEs = split('\t', $line); #some spaces in some TE names of Jainy... Tab is better
		my $lcRname = lc($TEs[0]);
		$TEs{$lcRname} = \@TEs;
	}	
	close $fh;
	return \%TEs;
}

#----------------------------------------------------------------------------
# get TE infos
# my $TE = get_TEs_infos($TEs) if ($TEs);
#----------------------------------------------------------------------------
sub get_TEs_infos {
	my $input = shift; #file name
	my %TEs = ();
	open(my $fh, "<", $input) or confess "\nERROR (Sub get_TEs_infos): could not open to read $input $!\n"; 
	LINE: while(<$fh>) {
		chomp (my $line = $_);
		next LINE if ($line !~ /\w/);
		my @TEs = split('\t', $line); 
		my $lcRname = lc ($TEs[0]);
		$TEs{$lcRname} = \@TEs;
	}	
	close ($fh);
	return \%TEs;
}

#############################################################################
# THE MAIN SUB THAT DOES THE WORK
# split age amounts [take care of overlaps as well]
# parseRM($p,$aa,$l,$filter,$contain,$genlen,$liblen,$age,$My,$TEdb,$srates,$TE,$path,$in,\@files,$k,$dir,$nonTE,$v);
#############################################################################
sub parseRM {
	my ($p,$aa,$l,$filter,$contain,$genlen,$liblen,$age,$My,$TEdb,$srates,$TE,$path,$in,$files,$k,$dir,$nonTE,$v) = @_;
	FILE: foreach my $f (@{$files}) {
		chomp $f;
		next FILE unless ($f); #ls is weird sometimes, blank values
		$f = $in."/".filename($f) if ($dir eq "y");
		$f = $path."/".filename($f) if ($dir eq "n");
		print STDERR "     -> $f..\n" if ($v);	
		print STDERR "        ..skipped, does not exist?\n" if ((! -e $f) && ($v));
		next FILE unless (-e $f);
				
		#prep landscape stuff if relevant
		my($max,$bin);
		if ($l ne "n") {
			($max,$bin) = split(",",$l);
			prep_landscape($f,$max,$bin,$TEage,$My);
		}
			
		#load it all in an array of array and filter at that point (reduce memory)
		my ($big);
		print STDERR "        ..loading in array..\n" if ($v);
		($big,$TE) = get_RMalign_array($f,$nonTE,$filter,$contain,$TE,$age,$TEage,$v) if ($f=~ /\.align/); 
		($big,$TE) = get_RMout_array($f,$nonTE,$filter,$contain,$k,$TE,$age,$TEage,$v) if ($f=~ /\.out$/);
		print STDERR "          ..done\n" if ($v);
		my $nb = @{$big};
		#loop through the array to store all %div per nt piece so that nt can be split into age category and bins [not the best memory usage wise but whatever]
		print STDERR "        ..Looping through array (parsing)..\n" if ($v);
		print STDERR "          Number of processed Gnames:\n" if ($v);
		my ($alldiv,$len,$landscape,$masked,$counts,$tot,$parsed,$id,$indel);
		$tot->{'nr'}=0;
		$tot->{'double'}=0;
		my $nbscaff = 0;
		my ($div,$del,$ins,$Gname,$Gst,$Gen,$strand,$RfullID,$Rid);
		LINE: for (my $i = 0; $i < $nb; $i++){
			($div,$del,$ins,$Gname,$Gst,$Gen,$strand,$RfullID,$Rid) = ($big->[$i][0],$big->[$i][1],$big->[$i][2],$big->[$i][3],$big->[$i][4],$big->[$i][5],$big->[$i][6],$big->[$i][7],$big->[$i][8]);			
			
			#To reduce memory usage I can loop on positions and store data as soon as it goes to the next Gname
			if (($i != 0) && ($Gname ne $big->[$i-1][3])) {
				($landscape,$masked,$counts,$tot,$parsed) = get_data($len,$alldiv,$TE,$f,$age,$TEage,$id,$indel,$max,$bin,$landscape,$masked,$counts,$tot,$parsed);
				undef $len;
				undef $alldiv;
				undef $id;
				undef $indel;
				$nbscaff++;
				print STDERR "          ..$nbscaff ($Gname - ne $big->[$i-1][3])\n";
# 				print STDERR "          ..$nbscaff done\n" if (($nbscaff <= 10000) && (($nbscaff / 10) =~ /^10*$/)); #should mean that it is a multiple of 10
# 				print STDERR "          ..$nbscaff done\n" if (($nbscaff > 10000) && (($nbscaff / 10) =~ /^[1-9]00000*$/)); #increments of 10,000
			}	
			
			#OK so now, store stuff
			#First, see if %div should be converted in My using substitution rate provided and just replace the value
			if ($My ne "n") {				
				my $subs = $srates->{filename($f)};
				my $myears = $div / 100 / ($subs * 2);
				$div = $myears;
			}
						
			#Now store list of all %div(or My) per position as well as the associated repeat info
			for (my $n = $Gst; $n <=$Gen; $n++) {
				($alldiv->[$n])?($alldiv->[$n]=$alldiv->[$n]."#".$div):($alldiv->[$n]=$div); 
				$id->{$n}{$div}{'f'}=$RfullID;
				$id->{$n}{$div}{'r'}=$Rid;
				unless ($p eq "n") {
					$indel->{$n}{$div}[0]=$del;
					$indel->{$n}{$div}[1]=$ins;
				}	
				$len=$Gen unless ($len);
				$len=$Gen if ($Gen > $len); #increment last position to loop on
			}
		}
		print STDERR "          ..done\n" if ($v);
		
		#Now, print all data contained in $landscape,$masked,$counts,$tot
		print STDERR "        ..printing parsed data..\n" if ($v);
		$tot->{'double_na'} = 0; #value needed for now in print_split
		print_parsed($f,$masked,$counts,$tot,$parsed,$genlen,$liblen,$TE,$v) unless ($p eq "n");
		print_split($f,$dir,$masked,$counts,$tot,$v) unless ($aa eq "n");
		print_landscape($f,$landscape,$max,$bin,$My,$v) unless ($l eq "n");
		print STDERR "          .. done\n" if ($v);
	}
	return;		
}

#############################################################################
# THE SECOND SUB THAT DOES THE WORK
# Parse the stored positional data for each Gname
# ($landscape,$masked,$counts,$tot,$parsed) = get_data($len,$alldiv,$TE,$f,$age,$TEage,$id,$indel,$max,$bin,$landscape,$masked,$counts,$tot,$parsed);
#############################################################################
sub get_data {			
	my ($len,$alldiv,$TE,$f,$age,$TEage,$id,$indel,$max,$bin,$landscape,$masked,$counts,$tot,$parsed) = @_;
	my $nr_check;
	NT: for (my $n = 1; $n <= $len; $n++) {
		next NT unless ($alldiv->[$n]); #if not masking this nt
		
		#Do split age / landscape / parsing stuff:
		my @div = ();					
		($alldiv->[$n]=~/#/)?(@div = split("#",$alldiv->[$n])):($div[0] = $alldiv->[$n]); #get all the %div values			 
		@div = sort @div if ($div[0]); #sort, so that lowest can be kept
		my $fullid = $id->{$n}{$div[0]}{'f'}; #get corresponding repeat
		my $Rid = $id->{$n}{$div[0]}{'r'}; #get corresponding repeat
		my ($name,$rest) = split("-#-",$fullid);
		my $lcname = lc($name);
		my ($class,$fam) = ($TE->{$lcname}[1],$TE->{$lcname}[2]);
		$tot->{'nr'}++;
		$tot->{'double'}++ if ($div[1]); #this meant overlap at that position
				
		#Do parsing spe stuff
		if ($p ne "n") {
			#Now store values, basically increment for each base => nr
			($masked->{'pn'}{$name}{'nr'})?($masked->{'pn'}{$name}{'nr'}++):($masked->{'pn'}{$name}{'nr'}=1);
			($masked->{'pc'}{$class}{'nr'})?($masked->{'pc'}{$class}{'nr'}++):($masked->{'pc'}{$class}{'nr'}=1);
			($masked->{'pf'}{$fam}{'nr'})?($masked->{'pf'}{$fam}{'nr'}++):($masked->{'pf'}{$fam}{'nr'}=1);
			if ($div[1]) { #this meant overlap at that position
				($masked->{'pn'}{$name}{'double'})?($masked->{'pn'}{$name}{'double'}++):($masked->{'pn'}{$name}{'double'}=1);
				($masked->{'pc'}{$class}{'double'})?($masked->{'pc'}{$class}{'double'}++):($masked->{'pc'}{$class}{'double'}=1);
				($masked->{'pf'}{$fam}{'double'})?($masked->{'pf'}{$fam}{'double'}++):($masked->{'pf'}{$fam}{'double'}=1);
			}
			#count fragments using the full ID; note that with .align files it will be a bit off
			unless ($nr_check->{'pnr'}{$fullid}) {
				($counts->{'pn'}{$name}{'nr'})?($counts->{'pn'}{$name}{'nr'}+=1):($counts->{'pn'}{$name}{'nr'}=1);
				($counts->{'pc'}{$class}{'nr'})?($counts->{'pc'}{$class}{'nr'}+=1):($counts->{'pc'}{$class}{'nr'}=1);
				($counts->{'pf'}{$fam}{'nr'})?($counts->{'pf'}{$fam}{'nr'}+=1):($counts->{'pf'}{$fam}{'nr'}=1);
				$nr_check->{'pnr'}{$fullid}=1;
			}
			#get total counts with repeat ID; note that with .align files it will be a bit off			
			unless ($nr_check->{'ptot'}{$Rid}) {
				($counts->{'pn'}{$name}{'tot'})?($counts->{'pn'}{$name}{'tot'}+=1):($counts->{'pn'}{$name}{'tot'}=1);
				($counts->{'pc'}{$class}{'tot'})?($counts->{'pc'}{$class}{'tot'}+=1):($counts->{'pc'}{$class}{'tot'}=1);
				($counts->{'pf'}{$fam}{'tot'})?($counts->{'pf'}{$fam}{'tot'}+=1):($counts->{'pf'}{$fam}{'tot'}=1);
				$nr_check->{'ptot'}{$Rid}=1;
			}
			#make the list of %div/My, %del, %ins to be able to do median for each repeat. Average can be done at the same time.
			$parsed->{$name}{'div'}?($parsed->{$name}{'div'}=$parsed->{$name}{'div'}.",".$div[0]):($parsed->{$name}{'div'}=$div[0]);
			$parsed->{$name}{'del'}?($parsed->{$name}{'del'}=$parsed->{$name}{'del'}.",".$indel->{$n}{$div[0]}[0]):($parsed->{$name}{'del'}=$indel->{$n}{$div[0]}[0]);
			$parsed->{$name}{'ins'}?($parsed->{$name}{'ins'}=$parsed->{$name}{'ins'}.",".$indel->{$n}{$div[0]}[1]):($parsed->{$name}{'ins'}=$indel->{$n}{$div[0]}[1]);
		}	
		
		#Do split age spe stuff
		my $type;
		if ($aa ne "n") {
			$type = det_age_type($div[0],$f,$age,$TEage,$name);
			#Now store values, basically increment for each base => nr
			($masked->{'a'}{$type}{'nr'})?($masked->{'a'}{$type}{'nr'}++):($masked->{'a'}{$type}{'nr'}=1);
			if ($div[1]) { #this meant overlap at that position
				($masked->{'a'}{$type}{'double'})?($masked->{'a'}{$type}{'double'}++):($masked->{'a'}{$type}{'double'}=1);
			}
			#count fragments using the full ID; not really possible to do nr when .align, but check counting full IDs
			unless ($nr_check->{'a'}{$fullid}) {
				($counts->{'a'}{$type}{'nr'})?($counts->{'a'}{$type}{'nr'}+=1):($counts->{'a'}{$type}{'nr'}=1);
				$nr_check->{'a'}{$Rid}=1;
			}
			#get total counts with repeat ID; note that with .align files it will be a bit off			
			unless ($nr_check->{'atot'}{$fullid}) {
				($counts->{'a'}{$type}{'tot'})?($counts->{'a'}{$type}{'tot'}+=1):($counts->{'a'}{$type}{'tot'}=1);
				$nr_check->{'atot'}{$Rid}=1;
			}
		}
		
		#Do landscape -> basically increment+1 in the proper bin
		if ($l ne "n") {
			unless ($div[0] > $max) {
				FINDBIN: for (my $j = $bin; $j <= $max; $j+=$bin) {
					my $coord = $j-$bin; 
					if (($div[0] >= $coord) && ($div[0] < $j)) {
						$landscape->{"Rname"}{"$name\t$class\t$fam"}{$coord}++; #since it's per position here, simple increment
						$landscape->{"Rclass"}{$class}{$coord}++; 
						$landscape->{"Rfam"}{"$class\t$fam"}{$coord}++; 
						$landscape->{"Rage"}{$type}{$coord}++ if (($aa ne "n") && ($TEage eq "y"));
						last FINDBIN;
					}
				}	
			}
		}
	}
	return($landscape,$masked,$counts,$tot,$parsed);
}

#----------------------------------------------------------------------------
# load RM output .align in array of array
# my @big = get_RMalign_array($f,$nonTE,$filter,$contain,$TE,$age,$TEage,$v);
#----------------------------------------------------------------------------
sub get_RMalign_array {
	my ($f,$nonTE,$filter,$contain,$TE,$age,$TEage,$v) = @_;
	my @array = ();
	my $i = 0;
	my $prevskip;
	open my $fh, "<", $f or confess "\nERROR (Sub get_RMalign_array): can't open to read $f $!\n";	
	LINE: while(<$fh>) {
		chomp (my $line = $_);
		next LINE if (($line !~ /\w/) || (substr($line,0,1) =~ /\s/) || (substr($line,0,2) eq "C ") || (substr($line,0,6) eq "Matrix") || (substr($line,0,11) eq "Transitions") || (substr($line,0,8) eq "Gap_init"));		
		if (substr($line,0,1) =~ /[0-9]/) { #kindda the RM.out line, without block info
			my @line = split(/\s+/,$line);
			#773 12.31 10.00 0.00 Scaffold10260 147 276 (1002) Cannrnd-3_family-5#LTR/ERVL 440 582 (2) 1 m_b1349s001i0
			#453 15.99 4.90 1.90 Scaffold10260 1174 1275 (3) C Cannrnd-5_family-2606#LTR/ERV1 (0) 580 476 2 m_b1349s001i1
			#=> apparently there is only a strand column when it is C. Interesting.
			my ($sc,$div,$del,$ins,$Gname,$Gst,$Gen) = ($line[0],$line[1],$line[2],$line[3],$line[4],$line[5],$line[6]);
			my ($strand,$Rfullname);
			if (($line[8] eq "C") || ($line[8] eq "+")) { #just in case it gets changed in a later version
				($strand,$Rfullname) = ($line[8],$line[9]);
			} else {
				($strand,$Rfullname) = ("+",$line[8]);
			}			
			my ($Rname,$classfam)=split("#",$Rfullname);
			my ($Rclass,$Rfam,$Rclassfam) = get_Rclass_Rfam($classfam);
			
			#Correct class anf fam if relevant
			my $lcRname = lc($Rname);
			($Rclass,$Rfam) = ($age->{$lcRname}[1],$age->{$lcRname}[2]) if (($TEage eq "y") && ($age->{$f}{$lcRname}[1]));
			($TE->{$lcRname}[1])?(($Rclass,$Rfam) = ($TE->{$lcRname}[1],$TE->{$lcRname}[2])):(($TE->{$lcRname}[1],$TE->{$lcRname}[2]) = ($Rclass,$Rfam));
			
			#filter stuff (if relevant)
			my $skip = "no";
			$skip = get_RM_array_skip($Rname,$Rclass,$Rfam,$nonTE,$filter,$contain) unless (($nonTE eq "all") && ($filter eq "na"));
			$prevskip = $skip;
			next LINE if ($skip eq "yes");
			
			my $RfullID = $Rname."-#-".$Gname.":".$Gst."_".$Gen."_".$strand;
			my $Rid = $Rname."-#-".$Gname.":".$Gst."-".$Gen;
			my @new = ($div,$del,$ins,$Gname,$Gst,$Gen,$strand,$RfullID,$Rid);
			push(@array,\@new);
			$i++;
		}

		#grab corrected %div when it exists -> it will be value for $i-1 (NB: missing for many - seems to be simple repeats...??)
		if (substr($line,0,6) eq "Kimura") {	
			unless ($prevskip eq "yes") {		
				my ($whatever,$div) = split("=",$line);
				$div =~ s/\s//g;
				$array[$i-1]->[1]=$div;
			}
		}
	}	
	close ($fh);
	return (\@array,$TE);
}

#----------------------------------------------------------------------------
# load RM output .out in array of array
# my @big = get_RMout_array($f,$nonTE,$filter,$contain,$k,$TE,$age,$TEage,$v);
#----------------------------------------------------------------------------
sub get_RMout_array {
	my ($f,$nonTE,$filter,$contain,$k,$TE,$age,$TEage,$v) = @_;
	my @array = ();
	my $skipped = 0;
	open my $fh, "<", $f or confess "\nERROR (Sub get_RMout_array): can't open to read $f $!\n";
	LINE: while(<$fh>) {
		chomp (my $line = $_);		
		next LINE if (($line =~ /position|[sS]core|Gname/) || ($line !~ /\w/)); #skip headers and white lines
		$line =~ s/^\s+//; #remove spaces in beginning of lines
		my @line = split(/\s+/,$line);			
		my ($Gname,$Gst,$Gen,$Rname,$block) = ($line[4],$line[5],$line[6],$line[9],$line[14]);
		$skipped++ unless ($block);
		next LINE unless ($block); #apparently, that happens		
		my ($Rclass,$Rfam,$Rclassfam) = get_Rclass_Rfam($line[10]);
		my $RfullID = $Rname."-#-".$block.":".$Gname;
		my $Rid = $Rname."-#-".$Gname.":".$Gst."-".$Gen;
		
		#Correct class and fam if relevant + complete the $TE table
		my $lcRname = lc($Rname);
		($Rclass,$Rfam) = ($age->{$lcRname}[1],$age->{$lcRname}[2]) if (($TEage eq "y") && ($age->{$f}{$lcRname}[1]));
		($TE->{$lcRname}[1])?(($Rclass,$Rfam) = ($TE->{$lcRname}[1],$TE->{$lcRname}[2])):(($TE->{$lcRname}[1],$TE->{$lcRname}[2]) = ($Rclass,$Rfam));
		
		#filter stuff (if relevant)
		my $skip = "no";
		$skip = get_RM_array_skip($Rname,$Rclass,$Rfam,$nonTE,$filter,$contain) unless (($nonTE eq "all") && ($filter eq "na"));
		next LINE if ($skip eq "yes");
		
		#deal with %div:
		my $div = $line[1];						
		# correct the % of divergence if $k
		$div = -300/4*log(1-$div*4/300) if ($k ne "n"); #note that in perl, log in ln; log10
			
		my @new = ($div,$line[2],$line[3],$Gname,$line[5],$line[6],$line[8],$RfullID,$Rid);	
		push(@array,\@new);
	}		
	close ($fh);
	print STDERR "          WARN: $skipped lines skipped because they had no block info\n" if (($v) && ($skipped > 0));	
	return (\@array,$TE);
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
# Filter out nonTE stuff and TEs based on $filter while loading
# $skip = get_RM_array_skip($Rname,$Rclass,$Rfam,$nonTE,$filter,$contain) unless (($nonTE eq "all") && ($filter eq "na"));
#----------------------------------------------------------------------------
sub get_RM_array_skip {
	my ($Rname,$Rclass,$Rfam,$nonTE,$filter,$contain) = @_;
	my $skip = "no";
	#filter out non TE stuff unless asked not
	unless (($nonTE eq "nonTE") || ($nonTE eq "all")) {
		$skip = "yes" if ($Rclass eq "nonTE");
	}
	unless ($nonTE eq "all") {
		$skip = "yes" if (($Rclass eq "Simple_repeat") || ($Rclass eq "Low_complexity") 
					   || ($Rclass eq "Satellite") || ($Rclass =~ /RNA$/) || ($Rclass =~ /omeric$/) || ($Rclass eq "ARTEFACT"));
	}

	#filter out stuff if relevant		
	unless ($filter eq "na") {
		my ($f_type,$f_name) = split(",",$filter);
		my ($lcRname,$lcRclass,$lcRfam,$lcf_name) = (lc($Rname),lc($Rclass),lc($Rfam),lc($f_name));
		if ($contain eq "y") {
			#check if what's det as the filter is included in the names
			$skip = "yes" unless ((($f_type eq "name") && ($lcRname =~ /$lcf_name/))
						      || (($f_type eq "class") && ($lcRclass =~ /$lcf_name/)) 
						      || (($f_type eq "family") && ($lcRfam =~ /$lcf_name/)));
		} else {
			$skip = "yes" unless ((($f_type eq "name") && ($lcf_name eq $lcRname))
						      || (($f_type eq "class") && ($lcf_name eq $lcRclass)) 
						      || (($f_type eq "family") && ($lcf_name eq $lcRfam)));
		}	
	}
	return ($skip);
}

#----------------------------------------------------------------------------
# prep landscape outputs, prep
# prep_landscape($f,$max_bin,$bin_len,$TEage);
#----------------------------------------------------------------------------
sub prep_landscape {
	my ($f,$max,$bin,$TEage,$My) = @_;
	#4 outputs => 4 files to prep
	my ($outrname,$outfam,$outclass,$outlineage);
	my ($headrname,$headfam,$headclass,$headlineage);
	if ($My ne "n") {
		$outrname = "$f.landscape.My.Rname.tab";
		$outfam = "$f.landscape.My.Rfam.tab";
		$outclass = "$f.landscape.My.Rclass.tab";
		$outlineage = "$f.landscape.My.Rage.tab";
		$headrname = "\t\t\tMillion years bins:\nRname\tRclass\tRfam";
		$headfam = "\t\tMillion years bins:\nRclass\tRfam";
		$headclass = "\tMillion years bins:\nRclass";
		$headlineage = "\tMillion years bins:\nLineage";
	} else {
		$outrname = "$f.landscape.Div.Rname.tab";
		$outfam = "$f.landscape.Div.Rfam.tab";
		$outclass = "$f.landscape.Div.Rclass.tab";
		$outlineage = "$f.landscape.Div.Rage.tab";
		$headrname = "\t\t\t% of divergence bins:\nRname\tRclass\tRfam";
		$headfam = "\t\t% of divergence bins:\nRclass\tRfam";
		$headclass = "\t% of divergence bins:\nRclass";
		$headlineage = "\t% of divergence bins:\nLineage";	
	}
	prep_landscape_print($outrname,$headrname,$max,$bin);
	prep_landscape_print($outfam,$headfam,$max,$bin);
	prep_landscape_print($outclass,$headclass,$max,$bin);
	prep_landscape_print($outlineage,$headlineage,$max,$bin) if ($TEage eq "y");
	return;
}

#----------------------------------------------------------------------------
# prep landscape outputs, print
# prep_landscape_print($outrname,$headrname,$max,$bin);
#----------------------------------------------------------------------------
sub prep_landscape_print {
	my ($out,$head,$max,$bin) = @_;
	open(my $fh, ">", $out) or confess "\nERROR (sub prep_landscape_print): could not open to write $out $!\n";
	print $fh "$head";
	for (my $i = 0; $i < $max; $i+=$bin) {
		my $i2 = $i+$bin;
		print $fh "\t\[$i;$i2\[";
	}
	print $fh "\n\n";
	close $fh;
	return;
}
	
#----------------------------------------------------------------------------
# Det the age type for this $div
# my $type = det_age_type($div,$f,$age,$TEage,$Rname);
#----------------------------------------------------------------------------				
sub det_age_type {
	my ($div,$f,$age,$TEage,$Rname) = @_;	
	my $type = "na";
	if ($TEage eq "y") { #then it was a file, so use $Rname
		my $lcRname = lc($Rname);
		$type = $age->{$lcRname}[3];
	} else {
		my $fname = filename($f);
		my ($min,$max)=($age->{$fname}[0],$age->{$fname}[1]);
		$type = "LS" if ($div <= $min);
		$type = "A" if ($div >= $max);
		$type = "nd" if (($div < $max) && ($div > $min));
	}
	return $type;
}

#----------------------------------------------------------------------------
# print parsed data
# print_parsed($f,$masked,$counts,$tot,$parsed,$genlen,$liblen,$TE,$v) unless ($p eq "n");
#----------------------------------------------------------------------------
sub print_parsed {
	my ($f,$masked,$counts,$tot,$parsed,$genlen,$liblen,$TE,$v) = @_;
	my $fname = filename($f);
	
	#summary file
	my $out = "$f.parseRM.Summary.tab";
	print STDERR "          -> $out\n" if ($v);
	open(my $fh, ">", $out) or confess "\nERROR (sub print_parsed): could not open to write $out $!\n";
	print $fh  "#This file gives summary info for total, by class and by class/fam
#See file $f.parseRM.all-repeats.tab for more info about repeats one by one
#Overlap or double corresponds to DNA fragments that are masked by several elements. These amounts need to be subtracted in order to get more accurate TE amount.
#If a .align file was parsed, these amounts will be much higher than for the associated .out
#Note overlaps may not be parsed correctly if names are not formatted consistently (Rname#Rclass/Rfam)\n\n";
	$genlen->{$fname} = "nd" unless ($genlen->{$fname});
	my $pertotnr = "nd";
	$pertotnr = $tot->{'nr'} / $genlen->{$fname} * 100 if ($genlen->{$fname} ne "nd");
	print $fh "\n#TOTAL (nr)\n";
	print $fh "#nt_total\tnt_masked\t%_masked\tnt_masked_double\n";
	print $fh "$genlen->{$fname}\t$tot->{'nr'}\t$pertotnr\t$tot->{'double'}\n";		
	print $fh "\n#BY CLASS\n";
	print $fh "#class\tnt_masked\t%_masked\tnt_masked_double\n";
	foreach my $key (keys %{$masked->{'pc'}}) {
		my $per = "nd";
		$per = $masked->{'pc'}{$key}{'nr'} / $genlen->{$fname} * 100 if ($genlen->{$fname} ne "nd");
		$masked->{'pc'}{$key}{'double'}=0 unless ($masked->{'pc'}{$key}{'double'});
		print $fh "$key\t$masked->{'pc'}{$key}{'nr'}\t$per\t$masked->{'pc'}{$key}{'double'}\n";	
	}
	print $fh "\n#BY FAMILY\n";
	print $fh "#family\tnt_masked\t%_masked\tnt_masked_double\n";
	foreach my $key (keys %{$masked->{'pf'}}) {
		my $per = "nd";
		$per = $masked->{'pf'}{$key}{'nr'} / $genlen->{$fname} * 100 if ($genlen->{$fname} ne "nd");
		$masked->{'pf'}{$key}{'double'}=0 unless ($masked->{'pf'}{$key}{'double'});
		print $fh "$key\t$masked->{'pf'}{$key}{'nr'}\t$per\t$masked->{'pf'}{$key}{'double'}\n";
	
	}
	close $fh;
	
	#all repeats file
	my $allrep = "$f.parseRM.all-repeats.tab";
	print STDERR "          -> $allrep\n" if ($v);
	open($fh, ">", $allrep) or confess "\nERROR (sub print_parsed): could not open to write $allrep $!\n";
	print $fh "#Rname\tRclass\tRfam\tRlen\tFRG_NB\tFRG_NB_StartToEnd\tNR_FRG_NB\tLEN_MASKED_NR\tAVG_%DIV\tMED_%DIV\tAVG_%DEL\tMED_%DEL\tAVG_%INS\tMED_%INS\tAVG_LEN_MASKED\t%_GENOME\n\n";
	foreach my $name (keys %{$masked->{'pn'}}) {
		my $lc = lc($name);
		my $rlen = "nd";
		$rlen = $liblen->{$lc} if ($liblen->{$lc});
		my ($dela,$delm,$diva,$divm,$insa,$insm) = get_avg_med_from_list($parsed,$name);
		my $len = $masked->{'pn'}{$name}{'nr'};
		my $avglen = $len / $tot->{'nr'};
		my $perlen = "nd";
		$perlen = $len / $genlen->{$fname} if ($genlen->{$fname} ne "nd");
		print $fh "$name\t$TE->{$lc}[1]\t$TE->{$lc}[2]\t$rlen\t$counts->{'pn'}{$name}{'tot'}\tnd\t$counts->{'pn'}{$name}{'nr'}\t$len\t$diva\t$divm\t$dela\t$delm\t$insa\t$insm\t$avglen\t$perlen\n";
	}
	close $fh;
	return;
}

#----------------------------------------------------------------------------
# Get average and median from a list
# Here avg is OK because I have a value per position (not very good memory wise... I should have a list that increment +1 when same % follows)
# But median is not really median of fragments though.
# my ($div,$del,$ins) = get_avg_med_from_list($parsed,$name);
#----------------------------------------------------------------------------
sub get_avg_med_from_list {
	my ($parsed,$name) = @_;
	my @med = ();
	my @avg = ();
	foreach my $t (sort keys %{$parsed->{$name}}) { #del div ins in that order since sorted
		my @list = split (',', $parsed->{$name}{$t});
		push(@med,median(\@list));
		push(@avg,average(\@list));
	}	
	return($avg[0],$med[0],$avg[1],$med[1],$avg[2],$med[2]);
}

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

sub average { 
	my ($array_ref) = @_; 
	my $count = scalar @$array_ref;
	my $total = 0;
	foreach my $v (@$array_ref) {
		$total+=$v;
	}
	return($total/$count);
} 

#----------------------------------------------------------------------------
# print split age data
# print_split($f,$dir,$masked,$counts,$tot,$v) unless ($aa eq "n");
#----------------------------------------------------------------------------
sub print_split {
	my ($f,$dir,$masked,$counts,$tot,$v) = @_;
	my $out = $f.".splitage.tab";
	print STDERR "          -> $out\n" if ($v);	
	my $fall = $path.".splitage_all.tab";
	open my $fhall, ">>", $fall or confess "\nERROR (Sub split_age): could not open to write $fall $!\n" if ($dir eq "y"); 
	open my $fh, ">", $out or confess "\nERROR (Sub print_data): could not open to write $out $!\n";
	print $fh "#nr_masked = amount of masked nt to consider\n";
	print $fh "#tot_masked = total amount of masked nt, including redundant maskins / overlaps\n";
	print $fh "#double_masked = total amount of nt that were masked by at least 2 repeats = redundant masking / overlaps\n\n";	
	print $fh "#Input_file\tnr_masked\ttot_masked\tdouble_masked\tAgeCat\tCounts_this_age\tnr_Counts_this_age\tnr_masked_this_age\t%nr_masked_this_age\n\n";	
	foreach my $type (keys %{$masked->{'a'}}) {
		my $total = $tot->{'nr'} + $tot->{'double'};
		my $nr_per = $masked->{'a'}{$type}{'nr'}/$tot->{'nr'}*100;
		my $tot_c = "na";
		$tot_c = $counts->{$type}{'tot'} if ($counts->{$type}{'tot'});
		print $fh    "$f\t$tot->{'nr'}\t$total\t$tot->{'double'}\t$type\t$tot_c\t$counts->{$type}{'nr'}\t$masked->{'a'}{$type}{'nr'}\t$nr_per\n";
		print $fhall "$f\t$tot->{'nr'}\t$total\t$tot->{'double'}\t$type\t$tot_c\t$counts->{$type}{'nr'}\t$masked->{'a'}{$type}{'nr'}\t$nr_per\n" if ($dir eq "y");
	}	
	close $fh;
	close $fhall if ($dir eq "y"); 
	return;
}

#----------------------------------------------------------------------------
# print landscape data
# print_landscape($f,$landscape,$max,$bin,$My,$v) unless ($l eq "n");
#----------------------------------------------------------------------------
sub print_landscape {
	my ($f,$landscape,$max,$bin,$My,$v) = @_;
	my $n = "Div";
	$n = "My" if ($My ne "n");
	foreach my $type (keys %{$landscape}) {		
		my $out = $f.".landscape.$n.$type.tab";
		print STDERR "          -> $out\n" if ($v);
		open (my $fh, ">>", $out) or confess "\nERROR (sub print_landscape): could not open to write $out $!\n";		
		foreach my $key (keys %{$landscape->{$type}}) {
			print $fh "$key";
			for (my $k=0; $k<$max; $k+=$bin) {
				($landscape->{$type}{$key}{$k})?(print $fh "\t$landscape->{$type}{$key}{$k}"):(print $fh "\t0");
			}
			print $fh "\n";	
		}
		close $fh;
	}	
	return;
}



