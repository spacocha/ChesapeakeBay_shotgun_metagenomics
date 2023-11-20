#! /usr/bin/perl -w
#
#

	die "Usage: table list_of_divided_subtables> Redirect\n" unless (@ARGV);
($list, $table) = (@ARGV);
chomp ($table, $list);

open (IN, "<$list") or die "Can't open $list\n";
while ($file=<IN>){
    chomp ($file);
    next unless ($file);
    $first=1;
    open (IN2, "<$file") or die "Can't open $file\n";
    while ($line=<IN2>){
	    	chomp ($line);
	    	next unless ($line);
		($KO, @pieces)=split("\t", $line);
		if ($first){
			$first=0;
			(@headers)=@pieces;
		} else {
			$i=0;
			$j=@pieces;
			until ($i >=$j){
				if ($KO=~/^K[0-9]{5}\.[12]/){
					#This has the right format
					($replaceKO)=$KO=~/^(K[0-9]{5})\.[12]/;
					$replacehash{$replaceKO}++;
					$mathash{$KO}{$headers[$i]}+=$pieces[$i];
					$allheaders{$headers[$i]}++;
				} else {
					#there will be some NA, keep as NA
					$mathash{$KO}{$headers[$i]}+=$pieces[$i];
                                        $allheaders{$headers[$i]}++;
				}
				$i++;
			}
		}
	}
	close (IN2);
}
close (IN);

$first=1;
open (IN, "<$table") or die "Can't open $table\n";
while ($line=<IN>){
    chomp ($line);
    next unless ($line);
    ($KO, @pieces)=split("\t", $line);
    if ($first){
	(@headers)=@pieces;
        $first=0;
    } else {
	unless ($replacehash{$KO}){
		#don't include any of the replaced KOs
		$i=0;
		$j=@pieces;
		until ($i >=$j){
			$mathash{$KO}{$headers[$i]}+=$pieces[$i];
			$allheaders{$headers[$i]}++;
			$i++;
	    	}
	}
    }
}
close (IN);

print "KO";
foreach $head (sort keys %allheaders){
	print "\t$head";
}
print "\n";

foreach $KO (sort keys %mathash){
	print "$KO";
	foreach $header (sort keys %allheaders){
		if ($mathash{$KO}{$header}){
			print "\t$mathash{$KO}{$header}";
		} else {
			print "\t0";
		}
	}
	print "\n";
}
