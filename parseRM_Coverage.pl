#!/usr/bin/perl -w
#######################################################
# Author  :  Aurelie K
# version :  2.0 (see below)
# email   :  4urelie.k@gmail.com  
# Purpose :  It reads a repeat masker output file and gives the coverage per position for each repeat
#				-> output files can be directly used to plot in R
#######################################################
# UPDATES
#	- v1.0 = March 2012
#		Script written with Xiaoyu Zhuo
#	- v2.0 = May 2014
#		Re wrote the script with subroutines
#		-lib or most common length is kept for consensus length
#		Added several options (filter, min_len, min_frg, R etc)
#   - v2.1 = May 2014
#		R plots through Statistics::R module
#
#######################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Array::Transpose::Ragged qw/transpose_ragged/;
use Bio::SeqIO;
use Statistics::R; #required to make plots through R; comment this if not needed

my $version = "2.1";
my $usage = "\nUsage [$version]: 
    perl <scriptname.pl> -RMout <FastaFile.out> [-filter <type,name>] [-min_frg <X>] [-min_len <X>] [-big] [-TEs <TEclass>] [-v] [-R] [-Rfile]
	
	This script will create a file per repeat (Rname) and then use them to calculate coverage
	Output files can be directly used to plot the coverage in R

    MANDATORY ARGUMENT:	
    -RMout <FastaFile.out>
        => <FastaFile.out> = Repeat masker output file (typically, Genome.out).
           If several .out files to parse, concatenate them before using this script, using following command:
           cat file1.out file2.out file3.out > FastaFile.out
           OR, if all .out in current folder are to parse:
           cat *.out > FastaFile.out
	 
    OPTIONAL ARGUMENTS (flagged by brackets [...] around each of them)
    -filter <type,name>
        => run the script on only a subset of repeats. Case does not matter.
           The type can be: name, class or family.
           ex: name,nhAT1_ML => only fragments corresponding to the repeat named exactly nhAT1_ML will in the output 
               class,DNA => all repeats with class named exactly DNA (as in ...#DNA/hAT or ...#DNA/Tc1)
               family,hAT => all repeats CONTAINING hAT in family (as in ...#DNA/hAT or ...#DNA/hAT-Charlie)
    -min_frg <X>
        => filter for a minimum number of fragments (e.g. coverage not so useful if less than 3 fragments for example)
           ex: -min_frg 3   
    -min_len <X>
        => filter for a minimum length of the repeat consensus (in nt)
           ex: -min_len 80              
	-bigtable
		=> add an output: a big tabulated file with all repeats in it, in addition to a file per repeat
    -TEs <TEclass>
        => optional input file, with TE information as follow, 3 first columns mandatory: 
           Rname \\t Rclass \\t Rfam \\t Rclass/Rfam \\t etc (only first 3 columns will be used)
           an easy way to get it is to run my other script parseRM.pl.
           -> copy the first columns of the output all-repeats.tab, modify them accordingly and then copy this in a text file
    -v
        verbose mode, make the script talks to you
     -R
    	To directly get the images of the plots (png)
    -Rfile
    	To print a file with R command lines to plot the coverage graphs
    	Behavior will be different if bigtable or not
    -h|help 
        Print this help\n\n";


################################################################################
# Get arguments/options
################################################################################
my ($RMout,$TEclass,$filter,$min_frg,$min_len,$bigtable,$lib,$help,$R,$Rfile,$verbose);
GetOptions ('RMout=s' => \$RMout, 'filter=s' => \$filter, 'min_frg=s' => \$min_frg, 'min_len=s' => \$min_len, 'big' => \$bigtable, 'TEs=s' => \$TEclass, 'lib=s' => \$lib, 'h' => \$help, 'help' => \$help, 'v' => \$verbose, 'R' => \$R, 'Rfile' => \$Rfile);

#check step to see if mandatory argument is provided + if help
die $usage if ((! $RMout) || ($help));

print "\n --- Script parseRM_Coverage.pl started (v$version)\n" if ($verbose);

#get filter if relevant
my ($f_type,$f_name) = split(",",$filter) if ($filter);

#Get name and path of the RMout file + create folder for outputs
my $path = path($RMout);
my $name = filename($RMout);
print " --- making output directory\n" if ($verbose);
my $f_out = $f_name.".coverage" if ($filter);
$f_out =~ s/[\[\]]//g if ($filter); #problem of brackets in folder name
($filter)?($path = make_outdir($path,$name,$f_out)):($path = make_outdir($path,$name,"coverage"));
print "      -> $path\n" if ($verbose);

#Get consensus lengths if -lib provided
my $replen; #ref of hash that will contain consensus lengths
print " --- getting TE consensus lengths from $lib\n" if (($lib) && ($verbose));
$replen = get_fasta_lengths($lib) if ($lib);

#Get TE infos if provided
my $TE; #will contain TE info
print " --- getting TE infos from $TEclass\n" if (($TEclass) && ($verbose));
$TE = get_TEs_infos($TEclass) if ($TEclass);



################################################################################
# read input file - create files for each TE
################################################################################
print " --- Reading $RMout to get all infos\n" if ($verbose);
print "     With filter = $filter\n" if (($verbose) && ($filter));
print "     Only when repeat has more than $min_frg fragments\n" if (($verbose) && ($min_frg));
print "     Only when repeat consensus is longer that $min_len nt\n" if (($verbose) && ($min_len));

my %check = ();
my $lib_incomplete = 0;
#my %tots = ();
my %RMout = ();
my %which_len;
my %frgs = ();
open RMOUT, "<$RMout" or die print "ERROR (main): can't open $RMout $!\n";
LINE: while (<RMOUT>){
	chomp (my $line = $_);	
	next LINE if (($line =~ /position|[sS]core|Gname/) || ($line !~ /\w/)); #skip headers and white lines
	$line =~ s/^\s+//; #remove spaces in beginning of lines
	my ($Gstart,$Gend,$strand,$Rname,$Rclass,$Rfam,$Rclassfam,$Rstart,$Rend,$Rleft,$Gmasked) = get_info_from_RMout_file($line);
	
	#filter out stuff if relevant
	if ($filter) {
		next LINE unless ((($f_type eq "name") && (lc($f_name) eq lc($Rname))) ||
						 (($f_type eq "class") && (lc($f_name) eq lc($Rclass))) ||
						 (($f_type eq "family") && (lc($f_name) =~ /lc($Rfam)/)));
	}
	
	#check if there was a TEclass or TEage file provided, if not fake the hash
	unless ($TEclass){
		my @TE = ($Rname,$Rclass,$Rfam,"$Rclass/$Rfam");
		$TE->{lc($Rname)} = \@TE;
	}	
	#get class, fam if relevant, e.g. not already defined (case of TEclass provided but not all TEs in it)
	unless($TE->{lc($Rname)}) {
		my @TE = ($Rname,$Rclass,$Rfam,"$Rclass/$Rfam");
		$TE->{lc($Rname)} = \@TE;
		print "      !! $Rname not found in $TEclass => RM output class / fam used instead\n" if (($TEclass) && ($verbose) && (! $check{lc($Rname)}{"class"}));
		$check{lc($Rname)}{"class"}=1;
	}
	
	#store length of consensus unless $repinfo has value for this TE (eg was in the library provided)
	my $Rlen = $Rend + $Rleft; #brackets got replaced
	my $tempRname = $Rname;
	$tempRname =~ s/\[|\]//g;
	unless ($replen->{lc($tempRname)}) {
		print "      !! $Rname not found in $lib => most common reconstructed length will be used\n" if (($lib) && ($verbose) && (! $check{lc($Rname)}{"len"}));
		$check{lc($Rname)}{"len"}=1;
		$lib_incomplete = 1;
		($which_len{$Rname}{$Rlen})?($which_len{$Rname}{$Rlen}++):($which_len{$Rname}{$Rlen}=1);
	}
	
	#count number of fragments
	($frgs{$Rname})?($frgs{$Rname}++):($frgs{$Rname}=1);

	#print in files, with current Rlen
	my $Rfile = "$path/$Rname.txt";
	unless (-e $Rfile) {
		open RFILE, ">$Rfile" or die print "ERROR (main): can't open/create $Rfile $!\n";
		print RFILE "#Rstart\tRend\tRleft\tRlen\tStrand\n";
	} else {
		open RFILE, ">>$Rfile" or die print "ERROR (main): can't open/create $Rfile $!\n";
	}
	print RFILE "$Rstart\t$Rend\t$Rleft\t$Rlen\t$strand\n";
	close RFILE;

	#get counts for summary file
	#($tots{$Rname}->{"counts"})?($tots{$Rname}->{"counts"}++):($tots{$Rname}->{"counts"}=1);
}
close RMOUT;

#Now make decision for $Rlen using the $which_len hash => store in %Rlen
print " --- chose most common consensus lengths\n" if (($verbose) && (! $lib));
print " --- chose most common consensus lengths [not all sequences were in $lib]\n" if (($verbose) && ($lib) && ($lib_incomplete == 1));
my $Rlen = get_Rlen(\%which_len) unless (($lib) && ($lib_incomplete == 0));


################################################################################
# process files and create files with coverage
################################################################################
print " --- Going through files with coordinates to get coverage\n" if ($verbose);
my %big = () if ($bigtable);
my @files = `ls $path`;
my %covfilelist = (); #for $R
FILES: foreach my $file (@files) {
	chomp $file;		
	#get name of repeat studied from file name
	my $Rname = $1 if ($file =~ s/^(.*)\.txt$/$1/);
	
	#skip this one if not enough fragments to see a coverage trend [defined by user]
	next FILES if (($min_frg) && ($frgs{$Rname} < $min_frg));
	
	#skip this one if not long enough
	next FILES if (($min_len) && ($Rlen->{$Rname} < $min_len));
	
	# get consensus length
	my $tempRname = $Rname;
	$tempRname =~ s/\[|\]//g;
	$replen->{lc($tempRname)} = $Rlen->{$Rname} unless ($replen->{lc($tempRname)});
	my $len = $replen->{lc($tempRname)};
	
	#read file and stock values in hash that will be read after to get coverage
	my @coords = ();
	open IN, "<$path/$file.txt" or die print "ERROR (main): can't open $path/$file.txt $!\n";
	LINE: while (<IN>) {
		chomp (my $line = $_);
		next LINE if ($line =~ /^#/);
		my ($Rstart,$Rend,$Rleft,$Rlen,$strand) = split(/\s/,$line);		
		#Store coordinates 
		push(@coords,($Rstart,$Rend));
	}
	close IN;
	
	#initialize coverage array
	my @cov = ();
	for (my $c = 0; $c < $len; $c++){
		push(@cov, 0);
	}
	
	#fill coverage
	for (my $i=0, my $j=1; $i<$#coords+1; $i+=2,$j+=2){ #$j<$#coords+2
		my ($s,$e) = ($coords[$i], $coords[$j]);
		foreach my $k ($s-1..$e-1){
			$cov[$k]++;
		}
	}
	
	my $output = "$path/$Rname.coverage.tab";
	$covfilelist{$Rname} = "$Rname.coverage.tab" if ($R);
	open OUT, ">$output" or die print "ERROR (main): Failed to create $output $!\n";
	for (my $i=0; $i<=$#cov; $i++) { 
		print OUT "$cov[$i]\n";		
	}
	close OUT;

	$big{$Rname}=\@cov; 	#if bigtable option then memorize all coverages
}

#now if bigtable => get all repeats in one file
print " --- Getting bif file with coverage for all repeats\n" if (($verbose) && ($bigtable));
my $bigOut = "no";
if ($bigtable) {
	#get big table with everything
	my @coverage = ();
	foreach my $key (sort keys %big) {
		my @cov = @{ $big{$key} }; #de-reference array that was value
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
	$bigOut = "$path/_ALL.coverage.tab";
	open OUT, ">$bigOut" or die print "ERROR (main): Failed to create $bigOut $!\n";
	for ($line=0; $line<=$maxlength; $line+=1) { #for each line,
		for ($col=0; $col<=$catnb; $col+=1) {				#print all columns
			if (defined ($TR_coverage[$line][$col])) {
				print OUT "$TR_coverage[$line][$col]\t";			
			} else {
				print OUT "0\t";
			}
		}
		print OUT "\n";
	}
	close OUT;
}

#print R command lines if relevant
print " --- Printing R command lines\n" if (($verbose) && ($Rfile)); 
get_Rcmdlines(\%big,$bigOut,$replen,\%covfilelist,$TE) if ($Rfile);

#get plots through R
print " --- Getting R plots\n" if (($verbose) && ($R)); 
get_Rplots(\%big,$replen,\%covfilelist,$TE,$path) if ($R);


print " --- Script done -> see files in $path\n\n" if ($verbose); 
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
# $path = make_outdir($path,$name,"coverage");
#----------------------------------------------------------------------------
sub make_outdir {
	my $path = shift;
	my $name = shift;
	my $type = shift;
	my $c=0;
	my $full;
	($type)?($full = "$path/$name.$type.parsed"):($full = "$path/$name.parsed");
	until (!-d "$full.$c"){
		$c++;
	}
	mkdir ("$full.$c", 0755) or confess "ERROR (sub make_outdir): couldn't mkdir $path $!\n";
	return "$full.$c";
}

#----------------------------------------------------------------------------
# get_fasta_lengths (returns ref of hash, fasta IDs as keys)
# $rep_len = get_fasta_lengths($repfile");
#----------------------------------------------------------------------------
sub get_fasta_lengths {
	my $repfile = shift;
	my $lib = Bio::SeqIO->new(-file => $repfile, -format => "fasta") or confess "ERROR (sub get_fasta_lengths): can't open $repfile $!\n";
	my %replen = ();
	while( my $seq = $lib->next_seq() ) {
		my $id = $seq->display_id;
		my $len = $seq->length;
		$id =~ s/\[|\]//g;
		$id =~ s/^(.*)#/$1/;
		my $lcRname = lc ($id);
		$replen{$lcRname} = $len;
	}
	return \%replen;
}

#----------------------------------------------------------------------------
# get TE infos
#----------------------------------------------------------------------------
sub get_TEs_infos {
	my $input = shift; #file name
	my %TEs = ();
	open(my $input_fh, "<", $input) or confess print "ERROR (sub get_TEs_infos): could not open $input $!\n";
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
# get values from RM output line (simple array > split in subroutine) => while looping on simple array
# => my ($Gstart,$Gend,$strand,$Rname,$Rclass,$Rfam,$Rclassfam,$Rstart,$Rend,$Rleft,$Gmasked) = get_info_from_RMout_file($line);
#----------------------------------------------------------------------------
sub get_info_from_RMout_file {
	my $line = shift;
	my ($score,$div,$del,$ins,$Gname,$Gstart,$Gend,$Gleft,$strand,$Rname,$classfam,$repstart,$Rend,$repleft,$block) = split(/\s+/,$line);
	my ($Rstart,$Rleft);
	unless ($strand eq "+") {
		$strand = "-";
		$Rstart = $repleft;
		$Rleft = $repstart;
	} else {
		$Rstart = $repstart;
		$Rleft = $repleft;
	}
	$Rleft =~ s/[\(\)]//g; #remove parentheses
	my ($Rclass,$Rfam,$Rclassfam) = get_Rclass_Rfam($classfam);
	my $Gmasked = ($Gend - $Gstart + 1);
	
	return ($Gstart,$Gend,$strand,$Rname,$Rclass,$Rfam,$Rclassfam,$Rstart,$Rend,$Rleft,$Gmasked);
}

#----------------------------------------------------------------------------
# Get Rclassfam from RMout
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
# get Rlen from hash, returns ref of hash
# my $Rlen = get_Rlen(\%which_len);
#----------------------------------------------------------------------------
sub get_Rlen {
	my $which_len = shift;
	my %Rlen;
	foreach my $R (sort keys %{which_len}){
		my @lengths = sort {$which_len->{$R}{$b} <=> $which_len->{$R}{$a}} keys %{$which_len->{$R}};
		$Rlen{$R}=$lengths[0]; #most common length
	}
	return \%Rlen;
}

#----------------------------------------------------------------------------
# print R command lines to plot coverages
# get_Rcmdlines($big,$bigOut,$Rlen,$covfilelist);
#----------------------------------------------------------------------------
sub get_Rcmdlines {
	my ($big,$bigOut,$replen,$filelist,$TE) = @_;
	#output file	
	my $out = "$path/_ALL.R.cmdlines.txt";
	open OUTR, ">$out" or confess "ERROR (sub get_Rcmdlines): Failed to create $out $!\n";
	print OUTR "setwd(\"\") #fill with location of files on computer\n\n";
	print OUTR "dat<-read.csv(\"_ALL.R.cmdlines.txt\", sep=\"\t\", header = TRUE)\n" unless ($bigOut eq "no");
	#now loop
	my $nb = 1;
	my $i = 1;
	foreach my $Rname (keys %{ $big }) {
		my $filename = $filelist->{$Rname};
		#deal with consensus length and tick marks
		my $tempRname = $Rname;
		$tempRname =~ s/\[|\]//g;
		my $len = $replen->{lc($tempRname)};
		my $ticks = int( $len/100 );
		$ticks = 10 if ($ticks < 10);
		my $fullname = $Rname."#".$TE->{lc($Rname)}->[1]."/".$TE->{lc($Rname)}->[2];
		#now print
		print OUTR "\npar(mfrow=c(4,2))\n" if ($nb == 1);		
		if ($bigOut eq "no") {
			print OUTR "dat<-read.csv(\"$filename\", header = FALSE)\n";
			print OUTR "plot(dat\$V1, col=\"blue\",pch = 18, xlab = \"position in consensus\", ylab = \"coverage\", main=\"$fullname\", xaxp = c(0,$len,$ticks))\n";
		} else {
			print OUTR "plot(dat\$V$i, col=\"blue\",pch = 18, xlab = \"position in consensus\", ylab = \"coverage\", main=\"$fullname\", xaxp = c(0,$len,$ticks))\n";
		}
		$i++;	
		$nb = 0 if ($nb == 8);
		$nb++;
	}	
	close OUTR;
#Table[!is.na(Table$Col),])
}

#----------------------------------------------------------------------------
# get R plots
# get_Rplots($big,$replen,$filelist,$TE,$path);
#----------------------------------------------------------------------------
sub get_Rplots {
	my ($big,$replen,$filelist,$TE,$path) = @_;

	#Start R bridge
	my $R = Statistics::R->new() ;
	$R->startR ;

	foreach my $Rname (keys %{ $big }) {
		my $filename = $filelist->{$Rname};
		
		#get infos
		my $fullname = $Rname."#".$TE->{lc($Rname)}->[1]."/".$TE->{lc($Rname)}->[2];		
		my $tempRname = $Rname;
		$tempRname =~ s/\[|\]//g;
		my $len = $replen->{lc($tempRname)};
		
		#plot
		$R->send(qq`dat <- read.table("$path/$filename", header = FALSE)`);
		my $out = "$path/$filename.pdf"; 
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

