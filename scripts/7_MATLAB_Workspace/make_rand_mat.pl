#! /usr/bin/perl -w

die "Usage: mat > redirect\n" unless (@ARGV);
($mat) = (@ARGV);

chomp ($mat);

use List::Util qw(shuffle);

open (IN, "<$mat" ) or die "Can't open $mat\n";
while ($line =<IN>){
	chomp ($line);
	next unless ($line);
	($KO, @pieces)=split ("\t", $line);
	if (@headers){
		#save data in hash
		$i=0;
		$j=@pieces;
		until($i>=$j){
			$allheaders{$headers[$i]}++;
			$hashHK{$headers[$i]}{$KO}=$pieces[$i];
			$i++;
		}
	} else {	
		(@headers)=@pieces;
	}
}
close(IN);

foreach $header (sort keys %allheaders){
	@randarray=();
	@KOarray=();
	foreach $KO (sort keys %{$hashHK{$header}}){
		#Collect all values within each header
		push (@randarray, $hashHK{$header}{$KO});
		push (@KOarray, $KO)
	}
	@randarray=shuffle(@randarray);
	$i=0;
	$j=@randarray;
	until ($i>=$j){
		#values were shuffled, KOs were not
		$randhash{$KOarray[$i]}{$header}=$randarray[$i];
		$i++;
	}

}

print "KEGGID";
foreach $header (sort keys %allheaders){
	print "\t$header";
}
print "\n";

#make a matrix of values to compare
foreach $KO (sort keys %randhash){
	print "$KO";
	foreach $header (sort keys %allheaders){
			print "\t$randhash{$KO}{$header}";
	}
	print "\n";
}


