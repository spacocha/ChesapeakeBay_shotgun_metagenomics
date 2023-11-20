#! /usr/bin/perl -w
#
#	Use this program to make tab files for FileMaker
#

	die "Usage: input_file > Redirect output\n" unless (@ARGV);
	
	chomp (@ARGV);
	($file) = (@ARGV);
	$/ = ">";
	open (IN, "<$file") or die "Can't open $file\n";
	while ($line1 = <IN>){	
		chomp ($line1);
		next unless ($line1);
		$sequence = ();
		(@pieces) = split ("\n", $line1);
		($info) = shift (@pieces);
		(@infopieces)=split (" ", $info);
		if ($info=~/[Oo]xidative/){
			($newinfo)=$infopieces[0]."_Oxidative";
		} elsif ($info=~/[Rr]eductive/){
			($newinfo)=$infopieces[0]."_Reductive";
		} else {
			($newinfo)=$infopieces[0]."_Unclassified";
		}
		#die "$newinfo $info\n";
		print ">$newinfo\n";
		$sequence=();
		foreach $piece (@pieces){
			(@letters)=$piece=~/[A-z]+/g;
			foreach $let (@letters){
				if ($sequence){	
					$newsequence="$sequence"."$let";
					$sequence=$newsequence;
				} else {
					$sequence=$let;
				}
			}
		}
		#die "Finished $sequence\n";
		print "$sequence\n";

	}
	
	close (IN);

	
