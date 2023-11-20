#! /usr/bin/perl -w
#
#	Use this program to make tab files for FileMaker
#

	die "Usage: fasta_file tab_file index_for_header >Redirect output\n" unless (@ARGV);
	
	chomp (@ARGV);
	($fafile, $bfile, $index) = (@ARGV);
	$index--;
	open (IN, "<$bfile") or die "Can't open $bfile\n";
	while ($line1 = <IN>){	
		chomp ($line1);
		next unless ($line1);
		(@pieces)=split ("\t", $line1);
		#die "$index $pieces[$index]\n";
		$hash{$pieces[$index]}++;
	}
	
	close (IN);

	$/ = ">";
        open (IN, "<$fafile") or die "Can't open $fafile\n";
        while ($line1 = <IN>){
                chomp ($line1);
                next unless ($line1);
		($longhead, @seqs)=split ("\n", $line1);
		($head, @pieces)=split (" ", $longhead);
		#die "|$head|\t$longhead\n";
		if ($hash{$head}){
				print ">$line1\n";
		}
	}
	close (IN);
