#!/usr/bin/perl -w
#----------------------------------------------------------------------------
# Author  :  Aurelie Kapusta
# email   :  4urelie.k@gmail.com
# PURPOSE :  parsing Repeat Masker .out output and merge interrupted repeats
#----------------------------------------------------------------------------
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Data::Dumper;

my $VERSION = "v3.0";
my $SCRIPTNAME = "parseRM_merge_interrupted";
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
#        - implement -m to check for more interrupted repeats
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
       perl $SCRIPTNAME -i <genome.out> [-t <TEinfo.tab>] [-m] [-n] [-v] [-l] [-h]

    This script reads a Repeat Masker output (.out) and corrects coordinates based on 
    the nesting blocks (last column of RM.out). Once this is done, if -m is set 
    the script will then look for more interrupted repeats.
			   
    The output will contain edited coordinates for nesting TEs 
    (TEs fragmented by nesting events are merged in one line). 
    Note that coordinates of the nesting TE will be true: to calculate
    the repeat length use the additional column with its real lenght.
    The scores, %div %del %ins are averaged (BUT NOT PONDERED FOR NOW). 
    The ( ) around Gleft and Rleft are removed.
	
    MANDATORY ARGUMENTS:
    -i,--in (STRING) 
             File to analyse  = .out file from RepeatMasker
	
    OPTIONAL ARGUMENTS:
    -t,--tes (STRING) 
             Tabulated file with TEs information, with columns as follow: 
                Rname \\t Rclass \\t Rfamily
             This is to allow cleaner output for the classes (useful if some repeat names are not formatted 
             like: Rname#Rclass/Rfamily). To build this file, all repeat names + current class / fam can be 
             obtained easily using output of the parseRM.pl script (then checked and corrected)
    -m,--more (BOOL)
             To look for more interrupted repeats; will be merged if:
              - same sequence name (genomic), strand, and repeat name
              - repeat end and start (in consensus) within +/- 20 bp of each other
              - not farther apart than 50nt if adjacent, and not farther apart than 5kb per block if not
    -n,--nonTE (BOOL)
             To keep the non TE repeats. Default behavior is to skip the repeat when class/fam value is:
             Simple_repeat, Low_complexity, Satellite, *RNA, *omeric, ARTEFACT or nonTE.           
    -v,--v (BOOL)
             Verbose mode, make the script talks to you; here will also print a summary
             Print the version if only option
    -l,--log (BOOL)
             Print change log (updates)
    -h,--help (BOOL)
             Print this usage
\n";
	return 1;
}

#-----------------------------------------------------------------------------
#-------------------------- LOAD AND CHECK OPTIONS ---------------------------
#-----------------------------------------------------------------------------
my ($IN,$TES,$MORE,$NONTE);
my ($HELP,$V,$CHLOG);
GetOptions ('in=s'     => \$IN, 
            'tes=s'    => \$TES,
            'more'     => \$MORE,
            'nonTE'    => \$NONTE,
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
print STDERR "       -n\n" if ($NONTE && $V);

#Get TE infos if provided
my %TE = ();
if ($TES) {
	print STDERR " --- Getting TE info from $TES\n" if ($V);
	get_TEs_infos();
}

#Now load & loop on RM
print STDERR " --- Loading $IN\n" if ($V);
my %RM = ();
load_RM();

print STDERR " --- Looping through $IN and merge when same block ID\n" if ($V);
my %COUNT = ();
my @NEW = ();
loop_RM('b');

%RM = ();
my %CHECK = ();
my $NID = 1; #new ID
if ($MORE) {
	@NEW = sort { $a->[4] cmp $b->[4] || $a->[5] <=> $b->[5] || $a->[6] <=> $b->[6] } @NEW;
	print STDERR " --- Looping through for more interrupted repeats\n" if ($V);
	load_RM_NEW();
	@NEW = ();
	loop_RM('m');
}

print STDERR " --- Printing output file\n" if ($V);
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
		my $Rname = lc($l[9]);
		
		#In case there are some low complexity stuff, skip unless -n set
		if (! $NONTE) {
			next if ($TES && $TE{$Rname} =~ "nonTE");
			my $Rclass = $l[10];
			next if ($Rclass eq "Simple_repeat"
				  || $Rclass eq "Low_complexity"
				  || $Rclass eq "Satellite"
				  || $Rclass =~ /RNA$/
				  || $Rclass =~ /omeric$/
				  || $Rclass eq "ARTEFACT");
		}
				
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
sub loop_RM {
	my $t = shift;
	for my $id (keys %RM) {
		#merge in one annotation if relevant
		if ($RM{$id}[1]) {
			my $new = merge_blocks(\@{$RM{$id}},$t);
			#keep the block id
			$new->[14]=$id;
			#now push
			push (@NEW,$new);
			#count by classfam
			$COUNT{$t}{$new->[10]}++;
		} else {
			if ($t eq "b") {
				#add repeat length in additional column
				$RM{$id}[0]->[15] = $RM{$id}[0]->[6] - $RM{$id}[0]->[5] +1;
				#remove brackets
				$RM{$id}[0]->[7] =~ s/\((.+?)\)/$1/;
				if ($RM{$id}[0]->[8] eq "+") {
					$RM{$id}[0]->[13] =~ s/\((.+?)\)/$1/;
				} else {
					$RM{$id}[0]->[11] =~ s/\((.+?)\)/$1/;
				}	
			}
			push(@NEW,$RM{$id}[0]);
		}	
	}	
	return 1;
}

#----------------------------------------------------------------------------
sub merge_blocks {
	my $bl = shift;
	my $t = shift;
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
		#left:
		$bl[$i]->[7] =~ s/\((.+?)\)/$1/;
		$new[7]=$bl[$i]->[7] if (! $new[7] || $new[7] > $bl[$i]->[7]);
		for (my $j = 8; $j <= 10; $j++){ 
			$new[$j]=$bl[$i]->[$j];
		}
		#Rst & Rleft:
		if ($new[8] eq "+") {
			$bl[$i]->[13] =~ s/\((.+?)\)/$1/;
		} else {
			$bl[$i]->[11] =~ s/\((.+?)\)/$1/;
		}
		#Rst/left
		$new[11]=$bl[$i]->[11] if (! $new[11] || $new[11] > $bl[$i]->[11]);
		#Ren
		$new[12]=$bl[$i]->[12] if (! $new[12] || $new[12] < $bl[$i]->[12]);
		#Rleft/st
		$new[13]=$bl[$i]->[13] if (! $new[13] || $new[13] > $bl[$i]->[13]);
		#save real length
		my $len;
		if ($t eq "b") {
			$len = $bl[$i]->[6] - $bl[$i]->[5] +1;
		} else {
			$len = $bl[$i]->[15];
		}	
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
sub load_RM_NEW {
	#detect interrupted repeats and push them in hash - block id = the one of $j
	#same as load RM but this time checks for interrupted repeats in +3 frame
	LINE: for (my $j = 0; $j < $#NEW; $j++){		
		#next if that line has already been pushed
		next LINE if ($CHECK{$j});
		#so push the current line
		my $id = $NEW[$j]->[14];
		push(@{$RM{$id}}, $NEW[$j]);
		#search next TEs to see if anything is interrupted
		TE: for (my $i = $j+1; $i < $#NEW; $i++){
			#exit the loop if not the same Gname
			last TE if ($NEW[$j]->[4] ne $NEW[$i]->[4]);
			#next if not the same strand
			next TE if ($NEW[$j]->[8] ne $NEW[$i]->[8]);
			#next if not the same Rname
			next TE if ($NEW[$j]->[9] ne $NEW[$i]->[9]);	
			#deal with genomic distance
			if ($i == $j+1) {
				#next if adjacent blocks and not in 50nt in genomic of each other
				next TE if ($NEW[$i]->[5] - $NEW[$j]->[6] > 50);
			} else {
				#not adjacent; last TE if too far away (limit of 5kb per block in between)
				my $diff = ($i - $j)*5000;
				last TE if ($NEW[$i]->[5] - $NEW[$j]->[6] > $diff);
			}
			#next check if in 20nt in consensus
			my $dist;
			if ($NEW[$j]->[8] eq "+") {
				$dist = $NEW[$j+1]->[11] - $NEW[$j]->[12] + 1;
			} else {
				$dist = $NEW[$j+1]->[12] - $NEW[$j]->[13] + 1;
			}	
			next TE if ($dist > 20);
			#OK, went through for $i => push
			push(@{$RM{$id}}, $NEW[$i]);
			$CHECK{$i} = 1;
# 			print STDERR "@{$NEW[$j]}\n";
# 			print STDERR "    => MERGED WITH @{$NEW[$i]}\n";
			$COUNT{'b'}{$NEW[$i]->[10]}++;
		
		}
	}
	return;
}	

#----------------------------------------------------------------------------
sub print_RM {
	@NEW = sort { $a->[4] cmp $b->[4] || $a->[5] <=> $b->[5] || $a->[6] <=> $b->[6] } @NEW;
	#Output file
	my $header = "Score\t%Div\t%Del\t%Ins\tGname\tGstart\tGend\tGleft\tStrand\tRname\tRclass/fam\tRstart(left)\tRend\tRleft(start)\tblockID";
	open(my $fhc, ">", $CORR)  or confess "\n ERROR (sub print_RM): could not open to write $CORR!\n";
	print $fhc "$header\tRealMaskedLen\n";
	for (my $i = 0; $i <= $#NEW; $i++){
		my $line = join("\t",@{$NEW[$i]});
		print $fhc "$line\n";
	}
	close($fhc);
	return 1;
}

#----------------------------------------------------------------------------
sub print_summary {
	print STDERR " --- SUMMARY:\n";
	print STDERR "     Number of TEs merged based on block IDs, by class/family:\n";
	foreach my $cf (keys %{$COUNT{'b'}}) {
		print STDERR "     $cf\t$COUNT{'b'}{$cf}\n";
	}
	if ($MORE) {
		print STDERR "     Number TEs were merged on second round (-m), by class/family:\n";
		foreach my $cf (keys %{$COUNT{'m'}}) {
			print STDERR "     $cf\t$COUNT{'m'}{$cf}\n";
		}
	}	
	return 1;
}	

	