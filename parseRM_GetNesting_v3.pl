#!/usr/bin/perl -w
#----------------------------------------------------------------------------
# Author  :  Aurelie Kapusta; adapted from another script (Qi Wang, former student in C. Feschotte Lab)
# email   :  4urelie.k@gmail.com
# PURPOSE :  parsing Repeat Masker .out output and figure out nesting info
#----------------------------------------------------------------------------
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Data::Dumper;

my $VERSION = "v2.0";
my $SCRIPTNAME = "parseRM_GetNesting.pl";
my $CHANGELOG;
set_chlog();
sub set_chlog {
	$CHANGELOG = "
# v1.0 = Dec 2012
# v1.1 = Nov 2013
#		 - made loop on columns dynamic to take in account different outputs (with correction of fragmentation or not for ex.)
#		 - updated filtering out headers of RMoutput
# v1.2 = Mar 2014
#		 - changed structure of the code; subs, loops, \$NESTING hash => start reducing code length and variable amounts + make the code more clear
#		 - added option TEclass to correct class/fam
#		 - This was needed to implement counting of nesting cases based on classes
# v1.3 = Apr 2014
#	 	 - Bug fix for TEc (was not defined if there was not a / in the classfam column or RMoutput)
#	 	 - Bug fix 
# v2.0 = Oct 2017
#		 - Rewrite for form (subroutines etc)
#        - works on original RM.out file
#        - remove -a
#        - make use of the block number!
\n";
	return 1;
}

my $USAGE;
set_usage();
sub set_usage {
$USAGE = "
    PLEASE CITE: 
       $SCRIPTNAME, v$VERSION, https://github.com/4ureliek/Parsing-RepeatMasker-Outputs 

    Usage [v$VERSION]:   
       perl $SCRIPTNAME -i <genome.out> [-t <TEinfo.tab>] [-v] [-l] [-h]
	
    MANDATORY ARGUMENTS:
    -i,--in (STRING) 
             File to analyse  = .out file from RepeatMasker
	
    [OPTIONAL ARGUMENTS]:
    -t,--tes (STRING) 
             Tabulated file with TEs information, with columns as follow: 
                Rname \\t Rclass \\t Rfamily
             This is to allow cleaner output for the classes (useful if some repeat names are not formatted 
             like: Rname#Rclass/Rfamily). To build this file, all repeat names + current class / fam can be 
             obtained easily using output of the parseRM.pl script (then checked and corrected)
    -m,--more (BOOL)
             To look for more nesting blocks after they are merged based on the block IDs from RM
             [NOT IMPLEMENTED YET]
    -v,--v (BOOL)
             Verbose mode, make the script talks to you; here will also print a summary
             Print the version if only option
    -l,--log (BOOL)
             Print change log (updates)
    -h,--help (BOOL)
             Print this usage

    WHAT IT DOES: 
    This script reads a Repeat Masker output (.out) and corrects coordinates based on 
    the nesting blocks (last column of RM.out). Note that the TE that suffered 
    the insertion is the nesting TE. The TE that inserted is the nested TE.
			   
    Once this is done, if -m is set the script will then look for more nesting blocks, 
    see method below.
			   
    Two output files:
      - original lines of RM output (TEs only), with annotations of nesting/nested
      - all lines, with edited coordinates: TEs fragmented by nesting events are merged in one line 
        (noted by additional column with number of frags in it). 
        Note that coordinates of the nesting TE will be true, but to calculate
        the repeat length use the additional column with real lenght.
	 
    POST BLOCK ID METHOD, IF -m SET (from Qi Wang, former student of Cedric Feschotte):
    For each TE, the previous and next TE in the file are compared. 
    If the previous and next TEs have:
      - the same repeat name
      - the same strand
      - genomic end of the previous TE is within 50 bp of the genomic Start of the nested TE
      - genomic start of the next TE is within 50 bp of the genomic end of the nested TE
      - repeat end (in consensus) of the previous TE is within +/- 20 bp of the repeat Start (in consensus) of the next TE
    Then the TE is determined to be nested within the previous/next TE.
	 
    Only the following nesting block structures are searched for:
    Frg A in B:                           [BBBBBB][AAAAAA][BBBBBB]
    Frg A in B in C:              [CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC]
    Frg A in B in C in D: [DDDDDD][CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC][DDDDDD]
		
    Two indep. nested frg in C:                      [CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC]
    Two indep. nested frg in C, nested in D: [DDDDDD][CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC][DDDDDD]	
    
    Two frg in C:                      [CCCCCC][BBBBBB][AAAAAA][CCCCCC]
    Two frg in C, nested in D: [DDDDDD][CCCCCC][BBBBBB][AAAAAA][CCCCCC][DDDDDD]
    
    Three frg in C:                      [CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC]
    Three frg in C, nested in D: [DDDDDD][CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC][DDDDDD]
	\n";
	return 1;
}

#-----------------------------------------------------------------------------
#-------------------------- LOAD AND CHECK OPTIONS ---------------------------
#-----------------------------------------------------------------------------
my ($IN,$TES,$MORE);
my ($HELP,$V,$CHLOG);
GetOptions ('in=s'     => \$IN, 
            'tes=s'    => \$TES,
            'more'     => \$MORE,
            'log'      => \$CHLOG, 
            'help'     => \$HELP, 
            'v'        => \$V);

#check steps on the options
check_opt();
sub check_opt {
	die "\n Script $SCRIPTNAME version $VERSION\n\n" if (! $IN && ! $HELP && ! $CHLOG && $V);
	die $CHANGELOG if ($CHLOG);
	die $USAGE if (! $IN);
	die "\n -i $IN does not exist?\n\n" if (! -e $IN);
	die "\n -t $TES does not exist?\n\n" if ($TES && ! -e $TES);
	return 1;
}

#make STDERR buffer flush immediately
select((select(STDERR), $|=1)[0]); 

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
print STDERR "\n --- $SCRIPTNAME v$VERSION started, with:\n" if ($V);
print STDERR "       -i $IN\n" if ($V);
print STDERR "       -t $TES\n" if ($TES && $V);
print STDERR "       -m\n" if ($MORE && $V);

# Get TE infos if provided
my %TE = ();
if ($TES) {
	print STDERR " --- Getting TE info from $TES\n" if ($V);
	get_TEs_infos();
}

# Loop on RMout
print STDERR " --- Loading $IN\n" if ($V);
my %RM = ();
load_RM();

print STDERR " --- Looping through $IN and merge when same block ID\n" if ($V);
my %COUNT = ();
my @NEW = ();
loop_RM_blocks();

my @FINAL = ();
my %NESTING = ();
if ($MORE) {
	print STDERR " --- Looping through again and merge when nesting, as described in the usage\n" if ($V);
	my $NID = 1;
	%NESTING=(
		  '101' => 0,                   #[BBBBBB][AAAAAA][BBBBBB]
		 '21012' => 0,          #[CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC]
		'3210123' => 0, #[DDDDDD][CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC][DDDDDD]

		 '2X2X2' => 0,          #[CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC]
		'32X2X23' => 0, #[DDDDDD][CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC][DDDDDD]

		 '2XXX2' => 0,          #[CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC]
		'32XXX23' => 0, #[DDDDDD][CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC][DDDDDD]

		 '2012' => 0,           #[CCCCCC][AAAAAA][BBBBBB][CCCCCC]
		'320123' => 0,  #[DDDDDD][CCCCCC][AAAAAA][BBBBBB][CCCCCC][DDDDDD]
	);
	loop_RM_nesting();
	sub loop_RM_nesting {
		
		#TO IMPLEMENT
		@FINAL = @NEW;
		unlink @NEW;
		
		return 1;
	}
} else {
	@FINAL = @NEW;
	unlink @NEW;
}

@FINAL = sort { $a->[14] <=> $b->[14] } @FINAL;

print STDERR " --- Now print output files\n" if ($V);
my $CORR = $IN;
$CORR =~ s/\.out$/.Corr.out/;
print_RM();
print_summary() if ($V);
exit;

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#----------------------------------------------------------------------------
sub get_TEs_infos {
	open(my $fh, "<", $TES) or confess "\n ERROR (sub ): could not open to read $TES!\n";
	while(defined(my $l = <$fh>)) {	
		# 0       1       2         3        4       5			6
		#Rname	Rclass	Rfam	Rclassfam	%div	~Age	Acient/LineageSpe
		chomp ($l);
		next if ($l !~ /\w/ || substr($l,0,1) eq "#");
		my @TEs = split(/\t/, $l); 
		my $lcRname = lc ($TEs[0]);
		$TE{$lcRname} = \@TEs;
	}
	close $fh;
	return 1;
}

#----------------------------------------------------------------------------
sub load_RM {
	# Structure =
	# 0		1	2	3	4		5		6		7		8		9		10			11			12		13		14    15
	# Score	Div	Del	Ins	Gname	Gstart	Gend	Gleft	strand	Rname	Rclass/fam	Repstart	Rend	Repleft	ID    *
	open(my $fh, "<", $IN) or confess "\n ERROR (sub ): could not open to read $IN!\n";
	while(defined(my $l = <$fh>)) {	
		chomp($l);	
		next if ($l =~ /position|[sS]core|Gname/ || $l !~ /\w/);
		#remove spaces in beginning of lines
		$l =~ s/^\s+//; 
		$l =~ s/\*$//;
		my @l = split(/\s+/,$l);
		$l[8] = "-" if ($l[8] eq "C");
	
		#In case there are some low complexity stuff, we don't want them (are not TEs)
		my $Rname = lc($l[9]);	
		next if ($TES && $TE{$Rname} =~ "nonTE");
		my $Rclass = $l[10];
		next if ($Rclass eq "Simple_repeat"
		      || $Rclass eq "Low_complexity"
		      || $Rclass eq "Satellite"
		      || $Rclass =~ /RNA$/
		      || $Rclass =~ /omeric$/
		      || $Rclass eq "ARTEFACT");
		
		#correct class/family in @l if relevant
		if ($TES && $TE{$Rname}->[1]) {
			$l[10] = $TE{$Rname}->[1]."/".$TE{$Rname}->[2];
		}
				
		#now, for relevant lines, load in hash based on the blockID
		my $id = $l[14];
		push (@{$RM{$id}}, \@l);
	}
	close $fh;
	return;
}

#----------------------------------------------------------------------------
sub loop_RM_blocks {
	for my $id (keys %RM) {
		#merge in one annotation if relevant
		if ($RM{$id}[1]) {
			my $new = merge_blocks(\@{$RM{$id}});
			$new->[14]=$id;
			push (@NEW,$new);
			$COUNT{$new->[10]}++;
		} else {
			push(@NEW,$RM{$id}[0]);
		}	
	}	
	return 1;
}

#----------------------------------------------------------------------------
sub merge_blocks {
	my $bl = shift;
	my @bl = @{$bl};
	my @new = ();
	my @val = ();
	for (my $i = 0; $i <= $#bl; $i++){
		#to get average of the first values (TO DO: PONDER THEM)
		for (my $j = 0; $j <= 3; $j++){ 
			push(@{$val[$j]},$bl[$i]->[$j]);
		}
		$new[4]=$bl[$i]->[4];
		#Gst:
		$new[5]=$bl[$i]->[5] if (! $new[5] || $new[5] > $bl[$i]->[5]);
		#en:
		$new[6]=$bl[$i]->[6] if (! $new[6] || $new[6] < $bl[$i]->[6]);
		#left (in negative):
		$new[7]="nd";
		for (my $j = 8; $j <= 10; $j++){ 
			$new[$j]=$bl[$i]->[$j];
		}
		for (my $j = 11; $j <= 13; $j++){ 
			$new[$j]="nd";
		}
		my $len = $bl[$i]->[6] - $bl[$i]->[5] +1;
		$new[15]+=$len;
	}
	for (my $j = 0; $j <= 3; $j++){
		next unless ($val[$j]);
		$new[$j] = average(\@{$val[$j]});
	}
	return \@new;
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

#----------------------------------------------------------------------------
sub print_RM {
	# Output files
	my $header = "Score\t%Div\t%Del\t%Ins\tGname\tGstart\tGend\tGleft\tStrand\tRname\tRclass/fam\tRstart(left)\tRend\tRend(left)\tblockID";
	open(my $fhc, ">", $CORR)  or confess "\n ERROR (sub print_RM): could not open to write $CORR!\n";
	print $fhc "$header\tRealMaskedLen\n";
	for (my $i = 0; $i < $#FINAL; $i++){
		my $line = join("\t",@{$FINAL[$i]});
		print $fhc "$line\n";
	}
	close($fhc);
	return 1;
}

#----------------------------------------------------------------------------
sub print_summary {
	print STDERR " --- SUMMARY:\n";
	if ($MORE) {
		print STDERR "        $NESTING{'101'} cases like [BBBBBB][AAAAAA][BBBBBB]
				(TE A nested in TE B)
			$NESTING{'21012'} cases like [CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC]
				(TE A nested in TE B, nested in TE C)
			$NESTING{'3210123'} cases like [DDDDDD][CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC][DDDDDD] 
				(TE A nested in TE B, nested in TE C, nested in TE D)
		
			$NESTING{'2XXX2'} cases like [CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC]
				(TEs A, E and B are nested in TE C - A B E could be pieces of same TE badly masked, since they are close)
			$NESTING{'32XXX23'} cases like [DDDDDD][CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC][DDDDDD]
				(TEs A, E and B are nested in TE C, nested in TE D - A B E could be pieces of same TE badly masked, since they are close)
		
			$NESTING{'2X2X2'} cases like [CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC]
				(TEs X are 2 independent nesting events in TE C)
			$NESTING{'32X2X23'} cases like [DDDDDD][CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC][DDDDDD]
				(TEs X are 2 independent nesting events in TE C, nesting in D)
		
			$NESTING{'2012'} cases like [CCCCCC][AAAAAA][BBBBBB][CCCCCC] 
				(TEs A and B are nested in TE C - A B could be pieces of same TE badly masked, since they are close)
			$NESTING{'320123'} cases like [DDDDDD][CCCCCC][AAAAAA][BBBBBB][CCCCCC][DDDDDD]
				(TEs A and B are nested in TE C, nested in TE D - A B could be pieces of same TE badly masked, since they are close)
			
			
			Note that these could be false if the masking is messy.\n\n";
	}
	print STDERR " --- Number of blocks that were merged (nesting elements) by class/fam:\n";
	foreach my $class (keys %COUNT) {
		print STDERR "     $class\t$COUNT{$class}\n";
	}
	return 1;
}	


