#! /usr/bin/perl -w

die "uat: metabolic EOF matrix (38 x 38)
uct: common EOF matrix (38 x 38)
Usage: uat uct output_prefix\n" unless (@ARGV);
($umt, $uct, $prefix) = (@ARGV);
chomp ($umt);
chomp ($uct);
chomp ($prefix);

$lineno=1;
open (IN, "<$umt" ) or die "Can't open $umt\n";
while ($line =<IN>){
    chomp ($line);
    next unless ($line);
    (@pieces)=split(",", $line);
    #die "$pieces[0]\n";
    $i=0;
    $j=@pieces;
    until ($i>=$j){
	#EOF1 is $i = 0, so change index to match EOF num
	$ei=$i+1;
	#lineno corresponds to libraries
	#Those are actually related to number 1...47 by indenv
	$umthash{$ei}{$lineno}=$pieces[$i];
	$allheaders{$lineno}++;
	#die "ei $ei lineno $lineno pieces $pieces[$i]\n" if ($lineno == 2 && $ei == 2); 
	$i++;
    }
    $lineno++;
} 
close (IN);

$lineno=1;
open (IN, "<$uct" ) or die "Can't open $uct\n";
while ($line =<IN>){
    chomp ($line);
    next unless ($line);
    (@pieces)=split(",", $line);
    #die "$pieces[0]\n";
    $i=0;
    $j=@pieces;
    until ($i>=$j){
        #Only use the first 4 EOFs
        #EOF1 is $i = 0, so change index to match EOF num 
        $ei=$i+1;
        #$i is the library, which corresponds to the number in indenv
        #only keep the first four EOFs
        #lineno corresponds to libraries
        #Those are actually related to number 1...47 by indenv
        $ucthash{$ei}{$lineno}=$pieces[$i];
	$allheaders{$lineno}++;
        #die "ei $ei lineno $lineno pieces $pieces[$i]\n" if ($lineno == 2 && $ei == 2);  
        $i++;
    }   
    $lineno++;
} 
close (IN);

open (OUTX, ">${prefix}.X.txt" ) or die "Can't open ${prefix}.X.txt\n";
open (OUTY, ">${prefix}.Y.txt" ) or die "Can't open ${prefix}.Y.txt\n";
open (OUTLOG, ">${prefix}.log.txt" ) or die "Can't open ${prefix}.log.txt\n";
print OUTX "mEOF";
print OUTY "cEOF";
foreach $header (sort {$a <=> $b} keys %allheaders){
	print OUTX "\t$header";
	print OUTY "\t$header";
}
print OUTX "\n";
print OUTY "\n";

#make a matrix of values to compare
foreach $mEOF (sort {$a <=> $b} keys %umthash){
	foreach $cEOF (sort {$a <=> $b} keys %ucthash){
		#don't repeat
		next if ($done{$mEOF}{$cEOF});
		print OUTX "$mEOF";
		print OUTY "$cEOF";
		print OUTLOG "mEOF_${mEOF},cEOF_${cEOF}\t";
		$done{$mEOF}{$cEOF}++;
		foreach $header (sort {$a <=> $b} keys %allheaders){
				print OUTX "\t$umthash{$mEOF}{$header}";
				print OUTY "\t$ucthash{$cEOF}{$header}";
		}
		print OUTX "\n";
		print OUTY "\n";
	}
}

print OUTLOG "\n";

