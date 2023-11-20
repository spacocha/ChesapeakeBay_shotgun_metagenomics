#! /usr/bin/perl -w
#
#

	die "Usage: KO_list RPM_mat> Redirect\n" unless (@ARGV);
($kofile, $mat) = (@ARGV);
chomp ($kofile, $mat);

open (IN, "<${kofile}") or die "Can't open $kofile\n";
while ($line=<IN>){
    chomp ($line);
    next unless ($line);
    (@pieces)=split ("\t", $line);
    #which kos should we keep?
    #die "|$pieces[0]|\n";
    $kohash{$pieces[0]}++;
}
close (IN);

#Header contains sample names
$first=1;
open (IN, "<${mat}") or die "Can't open $mat\n";
while ($line=<IN>){
    chomp ($line);
    next unless ($line);
    ($node, @pieces)=split ("\t",  $line);
    #Headers include both headers and KO/taxa info
    #die "$node\n";
    ($Conf)=pop(@pieces);
    ($KGenome)=pop(@pieces);
    ($Genus)=pop(@pieces);
    ($Class)=pop(@pieces);
    ($domain)=pop(@pieces);
    ($KOIDwNA)=pop(@pieces);
    ($KOID)=pop(@pieces);
    #die "$Conf\t$KGenome\tGenus\t$Class\t$domain\t$KOID\n" unless ($first);
    if ($first){
	(@headers)=@pieces;
	$first=0;
    } else {
	$i=0;
    	$j=@pieces;	
	#Add the coverage for each KO that it was assigned	
	until ($i>=$j){
		#die "$Class\n" if ($Class eq "Cyanobacteria");
		if ($kohash{$KOID}){
			 #it's a photosystem gene, so divide Cyano non-cyano
			if ($Class eq "Cyanobacteria"){
				#die "Here\n";
				$newKOID="${KOID}.1";
			} else {
				$newKOID="${KOID}.2";	
			}
			$hash{$newKOID}{$headers[$i]}+=$pieces[$i];
                        $allheaders{$headers[$i]}++;
		}
		$i++;	
    	}
   }
}
close (IN);

print "KO";
foreach $head (sort keys %allheaders){
	print "\t$head";
}
print "\n";

foreach $KO (sort keys %hash){
	print "$KO";
	foreach $value (sort keys %allheaders){
		print "\t$hash{$KO}{$value}";
	}
	print "\n";
}

