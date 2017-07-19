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

my $version = "5.2";

my $changelog;
set_chlog();
sub set_chlog {
	$changelog = "
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
#            Bug fix - one sub was calling itself... It was probably the issue of non progressing
#	- v5.1 = 09 Jul 2017
#            Bug fix - hash per genome sequence were not emptied
#            Speed up the parsing - thanks to GitHub user WuChangCheng (BioWu), who used NYTProf to identify 
#               that the '\@div = sort \@div if (\$div[0])' line was super long, I changed way of doing things           
#	- v5.2 = 19 Jul 2017
#            Bug fix in the 'double_masking' & thus total_masked amounts (nr were OK)
#            Bug fix in loading the kimura corrected %div from .align
#            Check step in case no repeat to parse in a file (for example if only nonTE stuff)
#            Usage update
\n";
	return;
}

my $usage;
set_usage();
sub set_usage {
    $usage = "\nUsage [v$version]: 
    perl parseRM.pl -i <genome.(out|align)>
            [-p] [-f <genome.fa>] [-n] [-g <X>] [-r <repeat_library.fa>]
            [-a <file.txt> OR <X,X>] [-t] [-m <file.txt OR X,X>] [-l <max,bin>] 
            [-d] [-k] [-s <type>] [-e <file>] [-w <type,name>] [-c] [-v] [-u] [-h]
	
    This script will process RepeatMasker outputs .out or .align file(s), 
    with 3 non exclusive parsings that can all be set together:
     -p To get a summary of the masking, as well as amount or DNA, 
          counts of fragments, + several other details for each repeat name 
          (all-repeats file), family, class and total amount (summary file)
          It is very long right now, so for a .out file, used parseRM_simple.pl
     -a To determine the amounts of DNA in a genome that is masked by repeats 
          of different lineages / %divergence categories
     -l To split the amount of DNA by bins of %div or My, allowing to generate 
          landscape graphs for each repeat name, family or class (one output for each)
	
    Note: if all 3 options -a, -t and -l are set, there will be an additional output 
          with bins by %div or My, but by age categories (specified in -a) 
	
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
         than if the .out is parsed.

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
         If not provided/not found, the % of genome masked won't be calculated.
         Names should correspond before .align and .fa(sta), and/or .out and .fa(sta)
     -n,--nrem (BOOL)
         To remove Ns from the genome file before getting its length 
         (to calculate the percentages on the amount of DNA that is not Ns)
     -g,--glen (INT or STRING)
         Alternatively to genome file(s), you can provide the total length of the genome (in nt)
         If several genomes are being looked at (-d chosen for example) this can be a file,
         containing 2 columns: filename \\t X                           
            filename = name of the file to be parsed
            X = the value (in nt)
     -r,--rlib (STRING)
         To add the length of the consensus sequence included in the output,
         set here the library of consensus sequences used to mask the genome 
         If several maskings are parsed, you can concatenate the different libraries
         (same if RM was run with -lib, concatenate the RM library and the user one)   
         If a repeat could not be found, or if -r is not chosen, the value 
         in the output will be \"nd\".              
          
    OPTIONAL ARGUMENTS RELATED TO --age    
     -a,--age (STRING)
         Option1: load a file (need to also set -t)
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
               you can provide here a file containing 2 columns: filename \\t X                           
               filename = name of the file to be parsed
               X = the value (in nt)
     -t,--te (BOOL)
         Set this if -a points to a file
     -m,--my (INT or STRING)
         To output age data in My and not %div. 
         Substitution rates (per My) should be provided, such as:
         if -d not chosen, use -m substitution_rate (ex. -m 0.00199)
         if -d is chosen, then different values can be provided for each input files,
            if you use -m subst-rates.tab, and subst-rates.tab is a file with 2 columns: 
            filename_as_in_dir \\t substitution_rate
                         
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
         Or: -l 5,0.25 if you wish to look only at recent stuff in a very dynamic genome
     -m,--my (INT or STRING)
         To output bins in My and not %div: substitution rates need to be provided, see above
                         
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
         %div are already corrected, with the kimura equation
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
         Print this usage\n\n";      
	return;
}

#-----------------------------------------------------------------------------
#-------------------------- LOAD AND CHECK OPTIONS ---------------------------
#-----------------------------------------------------------------------------
my $nonTE = "no";
my ($in,$p,$aa,$land);
my ($gfile,$glen,$Nrem,$lib);
my ($filter,$contain,$TEage,$TEs);
my ($My,$dir,$k);
my ($help,$v,$chlog);
GetOptions ('in=s'     => \$in, 
            'dir'      => \$dir, 
            'parse'    => \$p, 
            'age=s'    => \$aa, 
            'land=s'   => \$land, 
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

#check steps on the options
check_opt();
#print some log if $v
print_log("1") if ($v);

#make STDERR buffer flush immediately
select((select(STDERR), $|=1)[0]); 

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
#----- Get list of input files if -dir, or just load the one from -in
my @files = ();
($dir)?(@files = `ls $in`):(push(@files,$in));
my $path = get_path($in);

#----- Prep --parse stuff if relevant
my $genlen = ();
my $liblen = ();
prep_parse() if ($p);

#----- Prep --age stuff if relevant
my $age = ();
prep_age() if ($aa);

#----- Prep --land stuff if relevant
my ($max,$bin) = split(",",$land) if ($land);

#----- Prep steps not specific to p, a or l
#Load substitution rates from $My if relevant
my $srates = load_val($My,"rates") if ($My);
#Deal with $TEs if relevant
my $TE = load_TE_info($TEs) if ($TEs);

#----- Now parse all the RM files
print STDERR "\n --- Now parsing RM output(s)\n" if ($v);
my $nb;
my $tot = ();
$tot->{'nr'}= 0;
$tot->{'double'}= 0;
#main parsers, per Gname:
my %lowdiv = ();
my %id = ();
my %indel = ();
my %double = ();
#processed data:
my $parsed = ();
my $masked = ();
my $counts = ();
my $nr_check = ();
my $landscape = ();
#Loading data
my $nbscaff = 0;
my $nbscafftmp = 0;
my $f;
my $big = ();
my ($div,$del,$ins,$Gname,$Gst,$Gen,$strand,$RfullID,$Rid); #'columns' in $big
my ($Rname,$classfam);
my ($block,$Rfullname);
my ($Rclass,$Rfam,$Rclassfam); 
my ($skipped,$prevskip);
my ($fname,$out,$allrep);
FILE: foreach my $file (@files) {
	chomp $file;
	next FILE unless ($file); #ls is weird sometimes, blank values??
	$f = $in."/".filename($file) if ($dir);
	$f = $path."/".filename($file) if (! $dir);
	$fname = filename($f);
	print STDERR "     -> $f..\n" if ($v);
	
	#Check if f can/should be parsed
	if (! -e $f) {
		print STDERR "        ..skipped, does not exist?\n" if ($v);
		next FILE;
	}
	if ($f !~ /\.align/ && $f !~ /\.out/) {
		print STDERR "        ..skipped, not .out or .align?\n" if ($v);
		next FILE;	
	}
	
	#Now parse:
	#filter & load files in an array of array
	print STDERR "        ..loading in array..\n" if ($v);
	$big = (); #reinitialize for each file
	load_RM_in_array();
	print STDERR "          ..done\n" if ($v);

	print STDERR "          WARN: no repeat to parse in $file? Skipping\n" if (! $big && $v);
	next FILE if (! $big);
			
	#loop through the array to store all %div per nt piece so that nt can be split 
	#into age category and bins [not the best memory usage wise]
	print STDERR "        ..Looping through array (parsing)..\n" if ($v);		
	$nb = @{$big};	
	parseRM_table();
		
	#Now, print all data contained in hash tables
	print STDERR "        ..printing parsed data..\n" if ($v);	
	$tot->{'double_na'} = 0; #value needed for now in print_split
	if ($p) {
		$allrep = "$f.parseRM.all-repeats.tab";
		print_parsed_summary();	
		print_parsed_allrep();
	}	
	print_age() if ($aa);
	print_landscape() if ($land);	
	print STDERR "          .. done\n" if ($v);
}
#----- Done - log & exit
print_log("2") if ($v);
exit;

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

#--------------------------------- MAIN SUBS ---------------------------------
#----------------------------------------------------------------------------
sub load_RM_in_array {
	my $i = 0;
	$prevskip = 0;
	$skipped = 0; #for RMout (if no block)
	open (my $fh, "<", $f) or confess "\nERROR (Sub get_RMalign_array): can't open to read $f $!\n";	
	LINE: while(<$fh>) {	
		chomp (my $l = $_);
		my $next = check_line($l);
		next LINE if ($next eq "y");
	
		#Load these values:
		if ($f=~ /\.align/) {
			get_RMalign_val($l,$i);
		} elsif ($f =~ /\.out/) {
			get_RMout_val($l);	
			next LINE if (! $block);
		}
		
		#get, + correct if needed, the Rclass and Rfam
		($Rclass,$Rfam,$Rclassfam) = get_Rclass_Rfam($Rname,$classfam);
		
		#filter stuff (if relevant)
		my $skip = "no";
		$skip = get_RM_array_skip($Rname,$Rclass,$Rfam) unless ($nonTE eq "all" && ! $filter);
		$prevskip = $skip;
		next LINE if ($skip eq "yes");
		
		#Reconstruct large array
		if ($f=~ /\.align/) {
			$RfullID = $Rname."-#-".$Gname.":".$Gst."_".$Gen."_".$strand;
		} else {
			#For .out, use blocks:
			$RfullID = $Rname."-#-".$block.":".$Gname;
		}	
		$Rid = $Rname."-#-".$Gname.":".$Gst."-".$Gen;
		
		#Now load
		if (substr($l,0,1) =~ /[0-9]/) {
			my @new = ($div,$del,$ins,$Gname,$Gst,$Gen,$strand,$RfullID,$Rid);
			push(@{$big},\@new);
			$i++;
		}		
	}		
	close $fh;
	print STDERR "          WARN: $skipped lines skipped because they had no block info\n" if ($v && $skipped > 0);
	return;
}		

#----------------------------------------------------------------------------
sub get_RMalign_val {
	my $l = shift;
	my $i = shift;
	#Load values:
	if (substr($l,0,1) =~ /[0-9]/) {
		#kindda the RM.out line, without block info
		#773 12.31 10.00 0.00 Scaffold10260 147 276 (1002) Cannrnd-3_family-5#LTR/ERVL 440 582 (2) 1 m_b1349s001i0
		#453 15.99 4.90 1.90 Scaffold10260 1174 1275 (3) C Cannrnd-5_family-2606#LTR/ERV1 (0) 580 476 2 m_b1349s001i1
		#=> apparently there is only a strand column when it is C. Interesting.
		my @l = split(/\s+/,$l);
		my $sc;
		($sc,$div,$del,$ins,$Gname,$Gst,$Gen) = ($l[0],$l[1],$l[2],$l[3],$l[4],$l[5],$l[6]);
		if (($l[8] eq "C") || ($l[8] eq "+")) { #just in case it gets changed in a later version
			($strand,$Rfullname) = ($l[8],$l[9]);
		} else {
			($strand,$Rfullname) = ("+",$l[8]);
		}			
		($Rname,$classfam)=split("#",$Rfullname);
	#check the %div:	
	} elsif (substr($l,0,6) eq "Kimura" && $prevskip eq "no") {	
		# (NB: missing for many - seems to be simple repeats...??)
		my ($whatever,$kimdiv) = split("=",$l);
		$kimdiv =~ s/\s//g;
		$big->[$i-1]->[0]=$kimdiv;
	}
	return;
}

#----------------------------------------------------------------------------
sub get_RMout_val {
	my $l = shift;
	$l =~ s/^\s+//; #remove spaces in beginning of lines
	my @l = split(/\s+/,$l);			
	my $sc;
	($sc,$div,$del,$ins) = ($l[0],$l[1],$l[2],$l[3]);
	($Gname,$Gst,$Gen) = ($l[4],$l[5],$l[6]);
	($strand,$Rname,$classfam,$block) = ($l[8],$l[9],$l[10],$l[14]);
	$skipped++ unless ($block);						
	# correct the % of divergence if $k
	$div = -300/4*log(1-$div*4/300) if ($k);
}

#-----------------------------------------------------------------------------
sub parseRM_table {
	my $len = 0;
	LINE: for (my $i = 0; $i < $nb; $i++){
		($div,$del,$ins) = ($big->[$i][0],$big->[$i][1],$big->[$i][2]);
		($Gname,$Gst,$Gen) = ($big->[$i][3],$big->[$i][4],$big->[$i][5]);
		($strand,$RfullID,$Rid) = ($big->[$i][6],$big->[$i][7],$big->[$i][8]);			
		
		#Process data
		#First, see if %div should be converted in My using substitution rate provided and just replace the value
		if ($My) {				
			my $subs = $srates->{filename($f)};
			my $myears = $div / 100 / ($subs * 2);
			$div = $myears;
		}
		#Now store lowest %div(or My) per position as well as the associated repeat info
		for (my $n = $Gst; $n <=$Gen; $n++) {
			$tot->{'tot'}++;
			$double{$n}++ if ($lowdiv{$n});
			$tot->{'double'}++ if ($lowdiv{$n});
			if (! $lowdiv{$n} || ($lowdiv{$n} && $lowdiv{$n} > $div)) {
				$lowdiv{$n} = $div;
				$id{$n}{$div}{'f'}=$RfullID;
				$id{$n}{$div}{'r'}=$Rid;
				if ($p) {
					$indel{$n}{$div}[0]=$del;
					$indel{$n}{$div}[1]=$ins;
				}	
			}		
		}	
					
		#Now parse all this if this is the last line of a Gname, or of file
		if ($i == $nb-1 || $Gname ne $big->[$i+1][3]) {
			my $totdb = keys %double;
			parseRM_table_divlist();
			$nbscaff++;
			$nbscafftmp++;
			#empty these hash for each Gname
			%lowdiv = ();
			%id = (); 
			%indel = ();
			%double = ();
			if ($nbscafftmp == 100) {
				print STDERR "          ..$nbscaff Gname done\n" if ($v);
				$nbscafftmp =0;
			}
		}
	}
	print STDERR "          ..done\n" if ($v);
	return;		
}

#----------------------------------------------------------------------------
sub parseRM_table_divlist {
	NT: foreach my $n (sort keys %lowdiv) {
		my $lowest = $lowdiv{$n};
		$RfullID = $id{$n}{$lowest}{'f'}; #get corresponding repeat
		my $rest;
		($Rname,$rest) = split("-#-",$RfullID);
		$Rid = $id{$n}{$lowest}{'r'}; #get corresponding repeat
		my $lcname = lc($Rname);
		($Rclass,$Rfam) = ($TE->{$lcname}[1],$TE->{$lcname}[2]);
		parse_all_parse($n,$lowest,$double{$n}) if ($p);			
		my $type = "na";
		$type = parse_all_age($lowest,$double{$n}) if ($aa);
		parse_all_land($type,$lowest) if ($land);
	}
	return;
}




#---------------------------------- GENERAL ---------------------------------
#----------------------------------------------------------------------------
sub check_opt {
	die "\n Script parseRM.pl version $version\n\n" if (! $in && ! $p && ! $aa && ! $land && ! $help && ! $chlog && $v);
	die $changelog if ($chlog);
	die $usage if ($help);
	die $usage if (! $in && ! $p && ! $aa && ! $land);
	die "\n please chose one of --parse, --age or --land (use -h to see usage)\n\n" if (! $p && ! $aa && ! $land);
	die "\n $in does not exist?\n\n" if ($in && ! -e $in);
	die "\n $lib does not exist?\n\n" if ($lib && ! -e $lib);
	check_file($glen) if ($glen);
	check_file($aa) if ($aa);
	check_file($My) if ($My);	
	#avoid / at the end of path in case it's a directory
	$in = $1 if (($dir) && ($in =~ /^(.*)\/$/));	
	return;
}

#----------------------------------------------------------------------------
sub check_file {
	my $data = shift;
	die "\n $data does not exist as a file but also not numerical?\n\n" if (! -e $data && $data !~ /[0-9\.,]+/);
	return:
}

#----------------------------------------------------------------------------
sub print_log {	
	my $type = shift;
	if ($type == 1) {
		print STDERR "\n--------------------------------------------------\n";
		print STDERR " --- Script parseRM.pl started (v$version), with:\n";
		print STDERR "      - Directory containing input files = $in (-d chosen)\n" if ($dir);
		print STDERR "      - Input file = $in\n" if (! $dir);	
		print STDERR "      - All non TE repeats will be filtered out\n" if ($nonTE eq "no");
		print STDERR "      - Elements with class = nonTE won't be filtered out (-simple nonTE chosen)\n" if ($nonTE eq "nonTE");
		print STDERR "      - Non TE repeats won't be filtered out (-simple all chosen)\n" if ($nonTE eq "all");
		print STDERR "      - Repeats will be filtered based on $filter (--what)\n" if ($filter);
		print STDERR "        regexp will be used and not exact match (--contain)\n" if ($filter && $contain);
		print STDERR "      - %div will be corrected using -300/4*log(1-%div*4/300)\n" if ($k);
		print STDERR "      => global parsing will be done, with following options:\n" if ($p);
		print STDERR "         - genome file(s) will be looked for (--fa chosen)\n" if ($gfile);
		print STDERR "         - Ns will be removed from the genome before getting its length (--nrem chosen)\n" if ($Nrem);
		print STDERR "         - length(s) are given instead (--glen $glen)\n" if ($glen);
		print STDERR "         - library file = $lib\n" if ($lib);
		print STDERR "         - no options\n" if ($p && ! $gfile && ! $Nrem && ! $glen && ! $lib);
		print STDERR "      => age split will be done, with following options:\n" if ($aa);
		print STDERR "         - Age info will be based on $aa\n" if ($aa);
		print STDERR "         - Age will be in My, using substitution rates ($My)\n" if ($aa && $My);
		print STDERR "      => landscape data will be generated, with following options:\n" if ($land);
		print STDERR "         - With max_bin,bin_len = $land\n" if ($land);
		print STDERR "         - Age will be in My, using substitution rates ($My)\n" if ($land && $My);
		print STDERR "--------------------------------------------------\n";
	} else {
		print STDERR "\n --- Script done\n";
		if ($dir) {
			print STDERR "    -> age split files: *.agesplit.tab for all files\n" if ($aa);
			print STDERR "                        + all in $in.splitage_all.tab\n" if ($aa);
			print STDERR "    -> parsed files: *.parseRM.all-repeats.tab\n" if ($p);
			print STDERR "                     *.parseRM.summary.tab\n" if ($p);
			print STDERR "    -> landscape files: *.landscape.*.Rname.tab\n" if ($land);
			print STDERR "                        *.landscape.*.Rfam.tab\n" if ($land);
			print STDERR "                        *.landscape.*.Rclass.tab\n" if ($land);
			print STDERR "                        *.landscape.*.Rage.tab\n" if ($land && $TEage && $aa);
		}
		print STDERR "--------------------------------------------------\n\n";
	}	
	return;
}

#----------------------------------------------------------------------------
sub load_val {
	my $data = shift;
	my $type = shift;
	my $hash = ();	
	if (-e $data) { #then it's a file
		open (my $fh, "<", $data) or confess "\nERROR (Sub load_val): could not open to read $data $!\n"; 
 		while (<$fh>) {
 			chomp (my $l = $_);
 			my ($fname,$d) = split(/\s+/,$l);
 			if ($type eq "rates") {
				$hash->{$fname}=$d;
			} else {
				$hash = load_val_hash($hash,$fname,$d);
			}
 		}
 		close ($fh);
	} else {
		my $fname = filename($in);
		if ($type eq "rates") {
			$hash->{$fname}=$data;
		} else {
			$hash = load_val_hash($hash,$fname,$data);
		}	
	}
	return($hash);
}

#----------------------------------------------------------------------------
sub load_val_hash {	
	my $hash = shift;
	my $fname = shift;		
	my $d = shift;
	my ($m,$n);
	($d =~ /,/)?(($m,$n) = split(",",$d)):(($m,$n)=($d,$d));
	my @v = (); #values
	($m < $n)?(@v=($m,$n)):(@v=($n,$m));				
	$hash->{$fname}=\@v;
	return $hash;
}	

#----------------------------------------------------------------------------
sub get_Rclass_Rfam {
	my($Rn,$cf) = @_;
	my ($Rc,$Rf);
	if ($cf =~ /\//) {
		($Rc,$Rf) = split(/\//, $cf);
	} else {
		$Rf = $cf;
		$Rf=~ s/^(.*)\..*$/$1/;
		$Rc = $cf;
		$Rc =~ s/^.*\.(.*)$/$1/;
	}
	my $Rcf = "$Rc/$Rf";
	
	#Correct class anf fam if relevant
	my $lcRname = lc($Rn);
	($Rc,$Rf) = ($age->{$lcRname}[1],$age->{$lcRname}[2]) if ($TEage && $age->{$f}{$lcRname}[1]);
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
	if ($f=~ /\.align/) {
		return "y" if ($l !~ /\w/ || substr($l,0,1) =~ /\s/ || substr($l,0,2) eq "C ");
		return "y" if (substr($l,0,6) eq "Matrix" || substr($l,0,11) eq "Transitions" || substr($l,0,8) eq "Gap_init");	
	} else {	
		return "y" if (($l =~ /position|[sS]core|Gname/) || ($l !~ /\w/));
	}
	return "n";
}

#----------------------------------------------------------------------------
sub get_RM_array_skip {
	my ($cRn,$cRc,$cRf) = @_;
	my $skip = "no";
	#filter out non TE stuff unless asked not
	unless (($nonTE eq "nonTE") || ($nonTE eq "all")) {
		$skip = "yes" if ($cRc eq "nonTE");
	}
	unless ($nonTE eq "all") {
		$skip = "yes" if (($cRc eq "Simple_repeat") 
		               || ($cRc eq "Low_complexity") 
					   || ($cRc eq "Satellite") 
					   || ($cRc =~ /RNA$/) 
					   || ($cRc =~ /omeric$/) 
					   || ($cRc eq "ARTEFACT"));
	}
	#filter out stuff if relevant		
	if ($filter) {
		my ($f_type,$f_name) = split(",",$filter);
		my ($lcRname,$lcRclass,$lcRfam,$lcf_name) = (lc($cRn),lc($cRc),lc($cRf),lc($f_name));
		if ($contain) {
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
sub increment_hash {
	my $hash = shift;	
	my $ptype = shift;
	my $data = shift;
	my $type = shift;
	my $dbl = shift;
	my @ntype = $type;
	push (@ntype,'double') if ($dbl);
	foreach my $ntype (@ntype) {
		if ($hash->{$ptype}{$data}{$ntype}) {
			$hash->{$ptype}{$data}{$ntype}++;
		} else {
			$hash->{$ptype}{$data}{$ntype}=1;
		}
	}
	return $hash;
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



#----------------------------- RELATED TO PARSE -----------------------------
#----------------------------------------------------------------------------
sub prep_parse {
	print STDERR " --- Prepping steps for --parse...\n" if ($v);
	#genome length
	if ($glen) {
		print STDERR "     - Loading sequence length(s) from $glen...\n" if ($v);
		$genlen = load_val($glen,"len");
	} elsif ($gfile) {
		print STDERR "     - Getting sequence length(s) from genome file...\n" if ($v);
		$genlen = get_tot_length(); 
	}
	#lib lengths
	if ($lib) {
		print STDERR "     - Getting lengths of consensus sequences from $lib...\n" if ($v);
		$liblen = get_all_lengths();
	}	
	return;
}

#----------------------------------------------------------------------------
sub get_tot_length {
	my @list;	
	if (! $dir) {
		my $fa = $in;
		$fa = $1 if ($in =~ /(.*)\.out/ || $in =~ /(.*)\.align/);
		$fa =~ $1 if ($fa =~ /(.*\.fa)sta$/);
		if (! -e $fa) {
			print STDERR "       WARN: $fa could not be found, % genome won't be determined\n";
			return;
		} else {
			push(@list,$fa);
		}	
	} else {
		@list = `ls $in`;
	}
	#Now get length & print it	
	my %len = ();	
	FA: foreach my $fa (@list) {
		chomp $fa;		
		next FA unless (($fa =~ /\.fa$/) || ($fa =~ /\.fasta$/) || ($fa =~ /\.fsa$/));
		print STDERR "        -> $fa\n" if ($v);	
		$fa = Nrem($fa) if ($Nrem);	
		my $do = 0;
		my $flen = "$fa.length";	
		if (-e $flen) {
			print STDERR "           total length has been previously calculated ($flen exists)\n" if ($v);
			open (my $fh, "<", $flen) or ((warn "           WARN: could not open to read $flen, but length will be recalculated ($!)\n") && ($do++));
			unless ($do == 1) {
				while (<$fh>) {
					$len{filename($in)}=$_;
				}
				close ($fh);
				print STDERR "             => ".$len{filename($in)}." nt\n" if ($v);
				return(\%len);
			}
		}		
		if ((! -e $flen) || ($do == 1)) {
			print STDERR "           obtaining total length\n" if ($v);
			#looping through fasta file
			my $id = "";
			my $ln = 0;
			my $c = 0;			
			open (my $fh, "<", $fa) or warn "           WARN: could not open to read $fa, % genomes won't be calculated ($!)\n";
			while (<$fh>) {
				chomp (my $l = $_);
				if (substr($l,0,1) eq ">") {
					#first get and print unless first header
					$len{filename($in)}+=$ln unless ($c == 0);
					$c=1;
					#store header and reinitialize length
					my @id = split(/\s+/,$l);
					$id = $id[0];
					$id =~ s/>//;
					$ln = 0;
				} else {
					$ln+=length($l);
				}
			}
			#get and print len last sequence
			$len{filename($in)}+=$ln;
			close ($fh);
			open (my $fhl, ">", $flen) or warn "           WARN: could not open to write $flen $!\n";
			print $fhl "$len{filename($in)}";
			close $fhl;
			print STDERR "           => ".$len{filename($in)}." nt\n" if ($v);
		}
	}	
	return (\%len);
}

#----------------------------------------------------------------------------
sub get_all_lengths {
	#Get lengths of all sequences and store that by file and sequence ID. 
	#Note that if some are not unique, it just replaces by last length.
	#Specifically for TEs.
	print STDERR "        -> $lib\n" if ($v);
	unless (($lib =~ /\.fa$/) || ($lib =~ /\.fasta$/) || ($lib =~ /\.fsa$/)) {
		print STDERR "       WARN: $lib not a fasta file? => skipped\n" if ($v);
		return;
	}	
	if (! -e $lib) {
		print STDERR "       WARN: $lib could not be found, % genome won't be determined\n";
		return;
	}
	my $lfile = "$lib.lengths";
	my $do = 0;
	my %len = ();
	if (-e $lfile) {
		print STDERR "           lengths have been previously calculated ($lfile exists) => extracting\n" if ($v);
		#extract lengths now
		open (my $lfh, "<", $lfile) or ((warn "           WARN: could not open to read $lfile, but lengths will be recalculated ($!)\n") && ($do++));
		unless ($do == 1) {
		while (<$lfh>) {
				chomp (my $l = $_);
				my ($id,$ln) = split(/\s+/,$l);
				my ($Rname,$Rclassfam) = split('#',$id);
				$len{lc($Rname)}=$ln;
			}	
			close ($lfh);
			return(\%len);
		}
	}
	if ((! -e $lfile) || ($do == 1)) {
		print STDERR "           obtaining lengths\n" if ($v);
		my $id = "";
		my $ln = 0;
		my $c = 0;
		open (my $fa_fh, "<", $lib) or warn "           WARN: could not open to read $lib, repeat lengths won't be in the output $!\n";
		open (my $len_fh, ">", $lfile) or warn "           WARN: could not open to write $lfile, but lengths will be calculated $!\n";
		while (<$fa_fh>) {
			chomp (my $l = $_);
			if (substr($l,0,1) eq ">") {
				#first get and print unless first header
				unless ($c == 0) {
					print $len_fh "$id\t$ln\n";
					my ($Rn,$Rc) = split('#',$id);
					$len{lc($Rn)}=$ln;
				}
				$c=1;
				#store header and reinitialize length
				my @id = split (/\s+/,$l);
				$id = $id[0];
				$id =~ s/>//;
				$ln = 0;
			} else {
				#get length; could be more than one line so increment
				$l+=length($l);
			}
		}
		#get and print len last sequence
		print $len_fh "$id\t$ln\n";
		my ($Rn,$Rcf) = split('#',$id);
		$len{lc($Rn)}=$ln;
		close ($fa_fh);
		close ($len_fh);
		print STDERR "           lengths are now extracted in $lfile\n" if ($v);
	}	
	return (\%len);
}

#----------------------------------------------------------------------------
sub Nrem {
	my $fa = shift;
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
		my $Gen = Bio::SeqIO->new(-file => $genometmp, -format => "fasta") 
		          or confess "\nERROR (sub Nrem): could not create Bio::SeqIO object from $genometmp $!\n";
		my $GenMinusN = Bio::SeqIO->new(-file => ">$genome", -format => "fasta") 
		                or confess "\nERROR (sub Nrem): could not open to write $genome $!\n";	
		while( my $seq2 = $Gen->next_seq() ) {
			$GenMinusN->write_seq($seq2);		
		}
	}
	unlink($genometmp);
	return $genome;
}

#----------------------------------------------------------------------------				
sub parse_all_parse {
	my $n = shift;
	my $div = shift;	
	my $dbl = shift;	
	
	#Now store values, basically increment for each base => nr
	$masked = increment_hash($masked,'pn',$Rname,'nr',$dbl);
	$masked = increment_hash($masked,'pc',$Rclass,'nr',$dbl);
	$masked = increment_hash($masked,'pf',$Rfam,'nr',$dbl);

	#count fragments using the full ID; note that with .align files it will be a bit off
	if (! $nr_check->{'pnr'}{$RfullID}) {
		$counts = increment_hash($counts,'pn',$Rname,'nr');
		$counts = increment_hash($counts,'pc',$Rclass,'nr');
		$counts = increment_hash($counts,'pf',$Rfam,'nr');
		$nr_check->{'pnr'}{$RfullID}=1;
	}
	if (! $nr_check->{'ptot'}{$Rid}) {
		$counts = increment_hash($counts,'pn',$Rname,'tot');
		$counts = increment_hash($counts,'pc',$Rclass,'tot');
		$counts = increment_hash($counts,'pf',$Rfam,'tot');
		$nr_check->{'ptot'}{$Rid}=1;
	}

	#make the list of %div/My, %del, %ins to be able to do median for each repeat. Average can be done at the same time.
	if ($parsed->{$Rname}{'div'}) {
		$parsed->{$Rname}{'div'}=$parsed->{$Rname}{'div'}.",".$div;
	} else {
		$parsed->{$Rname}{'div'}=$div;
	}
	increment_parsed('div',$div);
	increment_parsed('del',$indel{$n}{$div}[0]);
	increment_parsed('ins',$indel{$n}{$div}[1]);	
	return;	
}

#----------------------------------------------------------------------------				
sub increment_parsed {
	my $type = shift;
	my $val = shift;
	if ($parsed->{$Rname}{$type}) {
		$parsed->{$Rname}{$type}=$parsed->{$Rname}{$type}.",".$val;
	} else {
		$parsed->{$Rname}{$type}=$val;
	}
	return;
}

#----------------------------------------------------------------------------
sub print_parsed_summary {		
	$out = "$f.parseRM.summary.tab";
	open(my $fh, ">", $out) or confess "\nERROR (sub print_parsed_summary): could not open to write $out $!\n";
	print $fh  "#This file gives some summary info about the masking, total, by class and by family\n";
	print $fh  "#(all repeats are detailed in the file $allrep)\n";
	print $fh  "#Overlap or double corresponds to DNA fragments that are masked by several elements.\n";
	print $fh  "#These amounts need to be subtracted in order to get more accurate TE amount.\n";
	print $fh  "#If a .align file was parsed, these amounts will be much higher than for the associated .out\n";
	print $fh  "#Note that overlaps may not be parsed correctly if names are not formatted consistently (Rname#Rclass/Rfam)\n";	
	
	$genlen->{$fname} = "nd" unless ($genlen->{$fname});
	my $pertotnr = "nd";
	$tot->{'nr'} = $tot->{'tot'} - $tot->{'double'};
	$pertotnr = $tot->{'nr'} / $genlen->{$fname} * 100 if ($genlen->{$fname} ne "nd");

	print $fh "\n#TOTAL:\n";
	print $fh "#nt_total_in_genome\tnt_masked-minus-double\t%_masked\tFYI:nt_masked_double\n";
	print $fh "$genlen->{$fname}\t$tot->{'nr'}\t$pertotnr\t$tot->{'double'}\n";	
		
	print $fh "\n#BY CLASS\n";
	print $fh "#class\tnt_masked-minus-double\t%_masked\tFYI:nt_masked_double\n";	
	print_parsed_summary_details('pc',$fh);
	
	print $fh "\n#BY FAMILY\n";
	print $fh "#family\tnt_masked-minus-double\t%_masked\tFYI:nt_masked_double\n";
	print_parsed_summary_details('pf',$fh);
	close $fh;
	print STDERR "          -> $out\n" if ($v);
	return;
}	

#----------------------------------------------------------------------------
sub print_parsed_summary_details {	
	my $type = shift;
	my $fh = shift;
	foreach my $key (keys %{$masked->{$type}}) {
		my $nt = $masked->{$type}{$key}{'nr'};
		my $per = "nd";
		$per = $nt / $genlen->{$fname} * 100 if ($genlen->{$fname} ne "nd");
		my $db = 0;
		$db = $masked->{$type}{$key}{'double'} if ($masked->{$type}{$key}{'double'});
		my $tot = $nt + $db;
		print $fh "$key\t$nt\t$per\t$db\n";	
	}
}	

#----------------------------------------------------------------------------
sub print_parsed_allrep {			
	open(my $fh, ">", $allrep) or confess "\nERROR (sub print_parsed_allrep): could not open to write $allrep $!\n";
	print $fh "#Rname\tRclass\tRfam\tRlen\tFRG_NB_all";
#	print $fh "\tFRG_NB_StartToEnd";
	print $fh "\tFRG_NB_Reconstructed_repeats";
	print $fh "\tLEN_MASKED_NR\tAVG_%DIV\tMED_%DIV\tAVG_%DEL\tMED_%DEL\tAVG_%INS\tMED_%INS\tAVG_LEN_MASKED\t%_GENOME\n\n";
	foreach my $name (keys %{$masked->{'pn'}}) {
		my $lc = lc($name);
		my $rlen = "nd";
		$rlen = $liblen->{$lc} if ($liblen->{$lc});
		my ($dela,$delm,$diva,$divm,$insa,$insm) = get_avg_med_from_list($parsed,$name);
		my $len = $masked->{'pn'}{$name}{'nr'};
		my $avglen = $len / $tot->{'nr'};
		my $perlen = "nd";
		$perlen = $len / $genlen->{$fname} if ($genlen->{$fname} ne "nd");
		print $fh "$name\t$TE->{$lc}[1]\t$TE->{$lc}[2]\t$rlen\t$counts->{'pn'}{$name}{'tot'}";
#		print $fh "\tnd";
		print $fh "\t$counts->{'pn'}{$name}{'nr'}\t$len\t$diva\t$divm\t$dela\t$delm\t$insa\t$insm\t$avglen\t$perlen\n";
	}
	close $fh;
	print STDERR "          -> $allrep\n" if ($v);
	return;
}

#----------------------------------------------------------------------------
sub get_avg_med_from_list {
	my ($parsed,$name) = @_;
	my @med = ('nd','nd','nd');
	my @avg = ();
	foreach my $t (sort keys %{$parsed->{$name}}) { #del div ins in that order since sorted
		my @list = split (',', $parsed->{$name}{$t});
		push(@med,median(\@list)) unless ($f =~ /\.align/);
		push(@avg,average(\@list));
	}	
	return($avg[0],$med[0],$avg[1],$med[1],$avg[2],$med[2]);
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
	foreach my $v (@$array_ref) {
		$total+=$v;
	}
	return($total/$count);
} 



#------------------------------ RELATED TO AGE ------------------------------
#----------------------------------------------------------------------------
sub prep_age {
	if ($TEage) {
		print STDERR "\n --- Loading age info from $aa\n" if ($v);
		$age = load_TE_info($aa);
	} else {
		print STDERR "\n --- Loading age value ($aa)\n" if ($v);
		$age = load_val($aa,"div");
	}
	return;
}

#----------------------------------------------------------------------------
sub load_TE_info {
	my $input = shift;
	my %TEs = ();
	open(my $fh, "<", $input) or confess "\nERROR (Sub get_TEs_infos): could not open to read $input $!\n"; 
	LINE: while(<$fh>) {
		chomp (my $l = $_);
		next LINE if ($l !~ /\w/);
		my @TEs = split(/\t/, $l); 
		my $lcRname = lc($TEs[0]);
		$TEs{$lcRname} = \@TEs;
	}	
	close ($fh);
	return \%TEs;
}

#----------------------------------------------------------------------------				
sub det_age_type {
	my ($div,$Rname) = @_;	
	my $type = "na";
	if ($TEage) { #then it was a file, so use $Rname
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
sub parse_all_age {
	my $div = shift;
	my $div1 = shift;
	my $type = det_age_type($div,$Rname);
	#Now store values, basically increment for each base => nr
	#Now store values, basically increment for each base => nr
	$masked = increment_hash($masked,'a',$type,'nr',$div1);

	#count fragments using the full ID; not really possible to do nr when .align, but check counting full IDs
	if (! $nr_check->{'a'}{$RfullID}) {
		$counts = increment_hash($counts,'a',$type,'nr');
		$nr_check->{'a'}{$RfullID}=1;
	}
	#get total counts with repeat ID; note that with .align files it will be a bit off			
	if (! $nr_check->{'atot'}{$Rid}) {
		$counts = increment_hash($counts,'a',$type,'tot');
		$nr_check->{'atot'}{$Rid}=1;
	}
	return $type;
}

#----------------------------------------------------------------------------
sub print_age {
	my $out = $f.".splitage.tab";
	print STDERR "          -> $out\n" if ($v);	
	my $fall = $in.".splitage_all.tab" if ($dir);

	#Prep outputs if directory	
	open my $fh, ">", $out or confess "\nERROR (Sub print_age): could not open to write $out $!\n";
	prep_age_out_headers($fh);
	open my $fhall, ">", $fall or confess "\nERROR (Sub print_age): could not open to write $fall $!\n" if ($dir);
	prep_age_out_headers($fhall) if ($dir);

	foreach my $type (keys %{$masked->{'a'}}) {
		my $total = $tot->{'nr'} + $tot->{'double'};
		my $nr_per = $masked->{'a'}{$type}{'nr'}/$tot->{'nr'}*100;
		my $nr_c = "na";
		$nr_c = $counts->{'a'}{$type}{'nr'} if ($counts->{'a'}{$type}{'nr'});
		my $tot_c = "na";
		$tot_c = $counts->{'a'}{$type}{'tot'} if ($counts->{'a'}{$type}{'tot'});		
		print $fh    "$f\t$tot->{'nr'}\t$total\t$tot->{'double'}\t$type";
		print $fhall "$f\t$tot->{'nr'}\t$total\t$tot->{'double'}\t$type" if ($dir);
		print $fh    "\t$tot_c\t$nr_c\t$masked->{'a'}{$type}{'nr'}\t$nr_per\n";		
		print $fhall "\t$tot_c\t$nr_c\t$masked->{'a'}{$type}{'nr'}\t$nr_per\n" if ($dir);
	}	
	close $fh;
	close $fhall if ($dir); 
	return;
}

#----------------------------------------------------------------------------
sub prep_age_out_headers {
	my $fh = shift;
	print $fh "#nr_masked = amount of masked nt to consider\n";
	print $fh "#tot_masked = total amount of masked nt, including redundant maskins / overlaps\n";
	print $fh "#double_masked = total amount of nt that were masked by at least 2 repeats = redundant masking / overlaps\n\n";	
	print $fh "#Input_file\tnr_masked\ttot_masked\tdouble_masked\tAgeCat";
	print $fh "\tCounts_this_age\tnr_Counts_this_age\tnr_masked_this_age\t%nr_masked_this_age\n\n";	
	return;
}
	


#----------------------------- RELATED TO LAND ------------------------------
#----------------------------------------------------------------------------
sub parse_all_land {
	my $type = shift;
	my $div = shift;
	unless ($div > $max) {
		FINDBIN: for (my $j = $bin; $j <= $max; $j+=$bin) {
			my $coord = $j-$bin; 
			if ($div >= $coord && $div < $j) {
				$landscape->{"Rname"}{"$Rname\t$Rclass\t$Rfam"}{$coord}++; #since it's per position here, simple increment
				$landscape->{"Rclass"}{$Rclass}{$coord}++; 
				$landscape->{"Rfam"}{"$Rclass\t$Rfam"}{$coord}++; 
				$landscape->{"Rage"}{$type}{$coord}++ if ($aa && $TEage);
				last FINDBIN;
			}
		}	
	}	
	return;
}

#----------------------------------------------------------------------------
sub print_landscape {
	my $n = "Div";
	$n = "My" if ($My);
	foreach my $type (keys %{$landscape}) {		
		my $out = $f.".landscape.$n.$type.tab";
		open (my $fh, ">", $out) or confess "\nERROR (sub print_landscape): could not open to write $out $!\n";
		prep_landscape_out($n,$type,$fh);			
		foreach my $key (keys %{$landscape->{$type}}) {
			print $fh "$key";
			for (my $k=0; $k<$max; $k+=$bin) {
				if ($landscape->{$type}{$key}{$k}) {
					print $fh "\t$landscape->{$type}{$key}{$k}";
				} else {
					print $fh "\t0";
				}	
			}
			print $fh "\n";	
		}
		close $fh;
		print STDERR "          -> $out\n" if ($v);
	}	
	return;
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
	for (my $i = 0; $i < $max; $i+=$bin) {
		my $i2 = $i+$bin;
		print $fh "\t\[$i;$i2\[";
	}
	print $fh "\n\n";
	return;
}





