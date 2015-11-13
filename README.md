Parsing-RepeatMasker-Outputs => parseRM
========================================================
Last update: 2015 November 13


Of interest if you are using the software Repeat Masker, for transposable elements (TEs) annotation [see http://www.repeatmasker.org/].
These script may help refining a custom (de novo) TE library, but won't help you to actually make one.

This is a collection of perl scripts I wrote for me and my lab to facilitate TE annotation (see list below)
and the extraction of information from Repeat Masker output file ".out".

Note that for all of the scripts, a complete usage will be obtained by simply launch them without any argument or option -h

========================================================
parseRM.pl

    WHAT IT DOES
    Parse RepeatMasker outputs to get summary info for each repeat as well as masked amounts by class and family.
    This is the first basic script to run - most labs have their own!
    
    For the repeats, it will provide:
     - fragment number (frg nb): 
          all
          frg nb from start to end of consensus
          frg nb corrected for interupted repeats (using the "ID". e.g. from the Repeat Masker .out)
     - %div, ins, del: pondered average, median
     - length masked + %genome by this repeat
     - amount of DNA that is masked several times (usually 2) by this element 
      
    Summary will provide amounts and % masked by various class, families
    Plus the script provides a set of files with overlap info to help giving real amounts/%

========================================================

parseRM_GetLandscape.pl

    NOTE THAT IN THE LATEST REPEAT MASKER RELEASES, AN EQUIVALENT SCRIPT IS PROVIDED:
    -> "the Kimura divergence is now calculated for each alignment and placed in the *.algn files, 
    and we now make available our software for drawing repeat landscapes (util/createRepeateLandscape.pl)" 
    [repeatmasker.org]
    
    WHAT MY SCRIPT DOES: 
        For each fragment (line of repeat masker output), amount of DNA masked is put in a bin, 
        based on the % of divergence to the consenus of this fragment.
        This is done per repeat and per class, as well as per lineage is -TEage is used (see detailed options with the fake -h)
        This script generates outputs that are tabulated text files. 
  	
        If you have .aln files, you should correct the %div with CG correction (Kimura).
        How: check you RM version
  	    If > 4.0.5, it is included in the -aln files
  	    If < 4.0.5, look in the RepeatMasker directory, Utils directory, called calcDivergenceFromAlign.pl and run it
  	    			then run this script on the modified .out

========================================================

parseRM_GetNesting.pl

    WHAT IT DOES: 
        This script reads a Repeat Masker output (.out) and find nested groups
        see below for more details on method and outputs
        Note that the TE that suffered the insertion is the nesting TE. The TE that inserted is the nested TE.
        It provides 3 outputs:
            - only nested/nesting blocks
            - original lines of RM output, with annotations of nesting/nested
            - all lines, but with corrections of coordinates. TEs fragmented by nesting events are merged in one line (noted by additional column with number of frags in it)
            note that coordinates of the nesting TE will be true, but WRONG TO CALCULATE LENGTH MASKED BY IT. For that, use additional column with real lenght
	 
    METHOD (from Qi Wang)
        For each TE, the previous and next TE in the file are compared. If the previous and next TEs have:
            - the same Repeat Name
            - the same strand
            - genomic End of the previous TE is within 50 bp of the genomic Start of the nested TE
            - genomic Start of the next TE is within 50 bp of the genomic End of the nested TE
            - repeat End (in consensus) of the previous TE is within +/- 20 bp of the repeat Start (in consensus) of the next TE
            Then the TE is determined to be nested within the previous/next TE.
	 
    SUMMARY OF THE NESTED/NESTING STRUCTURES FOUND BY THE SCRIPT (frg = fragment):
        Frg A in B:                           [BBBBBB][AAAAAA][BBBBBB]
        Frg A in B in C:              [CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC]
        Frg A in B in C in D: [DDDDDD][CCCCCC][BBBBBB][AAAAAA][BBBBBB][CCCCCC][DDDDDD]
		
        Two indep. nested frg in C:                      [CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC]
        Two indep. nested frg in C, nested in D: [DDDDDD][CCCCCC][XXXXXX][CCCCCC][XXXXXX][CCCCCC][DDDDDD]	
		
        Two frg in C:                      [CCCCCC][BBBBBB][AAAAAA][CCCCCC]
        Two frg in C, nested in D: [DDDDDD][CCCCCC][BBBBBB][AAAAAA][CCCCCC][DDDDDD]
		   
        Three frg in C:                      [CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC]
        Three frg in C, nested in D: [DDDDDD][CCCCCC][AAAAAA][EEEEEE][BBBBBB][CCCCCC][DDDDDD]
		
========================================================

parseRM_Coverage.pl

    WHAT IT DOES: 
        This script will output the coverage of a repeat plotted on the consensus length
        (if the plotting fails, -Rfile will output command lines for R to get the plots) 
        It is identical to the script in the Utils of this repository: https://github.com/4ureliek/TEanalysis
        Here, for a simple parsing of repeat masker output, use -type RMout

        Check the usage for the many options (such as filtering etc)!
        perl parseRM_Coverage.pl -h

    NOTE: 
        Requires the Statistics::R perl module
        
========================================================

parseRM_ExtractSeqs_P.pl

    WHAT IT DOES: 
        This script reads a Repeat Masker output (.out) and extract sequences.
        It will reconstructed interrupted repeats, and allows extraction of a random subset + its complementary
  
        Check the usage for the many options (such as filtering etc)!
        perl parseRM_ExtractSeqs_P.pl -h
 
    NOTE: 
        Requires the forks perl module
