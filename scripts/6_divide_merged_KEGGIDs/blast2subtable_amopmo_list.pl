#! /usr/bin/perl -w

die "Usage: blast KO full_table > Redirect\n" unless (@ARGV);
($blastfile, $oKO, $matfile) = (@ARGV);
chomp ($blastfile);
chomp ($oKO);
chomp ($matfile);

open (IN, "<$blastfile" ) or die "Can't open $blastfile\n";
while ($line =<IN>){
    chomp ($line);
    next unless ($line);
    ($node, $blastmatch)=split("\t", $line);
    #Keep nodes in the file and assign divided KOs based on blast
    #The blastmatch will include info about type at end
    #Determine what type based on end 
    (@pieces)=split("_", $blastmatch);
    ($type)=pop(@pieces);
    #die "$type\n";
    if ($type=~/amo/){
	    $dKO="${oKO}.1";
    } elsif ($type=~/pmo/){
    	    $dKO="${oKO}.2";
    } elsif ($type=~/xmo/){
	    $dKO=$oKO;
    } else {
	    $dKO="NA";
    }
    $hash{$node}=$dKO;
}
close(IN);

open (IN, "<$matfile" ) or die "Can't open $matfile\n";
while ($line =<IN>){
    chomp ($line);
    next unless ($line);
    ($node, @pieces) = split ("\t", $line);
    ($tKO)=pop(@pieces);
    if (@headers){
	    if ($hash{$node}){
		#die "Made it\n";
 		#only work with nodes that are in divided database
		$i=0;
		$j=@pieces;
		until ($i >=$j){
			$mathash{$hash{$node}}{$headers[$i]}+=$pieces[$i];
			$allheaders{$headers[$i]}++;
 			$i++;
		}
	    }
    } else {
	    (@headers)=@pieces;
    }
} 
close (IN);

print "KOID";
foreach $header (sort keys %allheaders){
	print "\t$header";
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

