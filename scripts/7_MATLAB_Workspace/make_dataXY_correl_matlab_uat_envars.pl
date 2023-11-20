#! /usr/bin/perl -w

die "uat: EOF matrix (38 x 38)
MATLAB_metadata_3: export from regression_models_3.m (38 x 27)
output_prefix: The beginning part of the file names
Usage: uat MATLAB_metadata_3 output_prefix\n" unless (@ARGV);
($umt, $metadata, $prefix) = (@ARGV);
chomp ($umt);
chomp ($metadata);
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
	#Only use the first 2 EOFs
	if ($i<=1){
		#EOF1 is $i = 0, so change index to match EOF num
		$ei=$i+1;
		#$i is the library, which corresponds to the number in indenv
		#only keep the first four EOFs
		#lineno corresponds to libraries
		#Those are actually related to number 1...47 by indenv
		$umthash{$ei}{$lineno}=$pieces[$i];
		$allheaders{$lineno}++;
		#die "ei $ei lineno $lineno pieces $pieces[$i]\n" if ($lineno == 2 && $ei == 2); 
	}
	$i++;
    }
    $lineno++;
} 
close (IN);

$lineno=1;
open (IN, "<$metadata" ) or die "Can't open $metadata\n";
while ($line =<IN>){
	chomp ($line);
	next unless ($line);
	(@pieces)=split (",", $line);
	if ($lineno==1){
		(@headers)=@pieces;
		$lineno++;
		next;
	} else {
		#lineno - 1 are libraries
		$si=$lineno-1;
		$i=0;
		$j=@pieces;
		until ($i >=$j){
			$datahash{$headers[$i]}{$si}=$pieces[$i];
			$allheaders{$si}++;
			$i++;
		}
		$lineno++;
	}
}
close(IN);

open (OUTX, ">${prefix}.X.txt" ) or die "Can't open ${prefix}.X.txt\n";
open (OUTY, ">${prefix}.Y.txt" ) or die "Can't open ${prefix}.Y.txt\n";
open (OUTLOG, ">${prefix}.log.txt" ) or die "Can't open ${prefix}.log.txt\n";
print OUTX "Metadata";
print OUTY "EOF";
foreach $header (sort {$a <=> $b} keys %allheaders){
	print OUTX "\t$header";
	print OUTY "\t$header";
}
print OUTX "\n";
print OUTY "\n";

#make a matrix of values to compare
foreach $datapt (sort keys %datahash){
	foreach $umt (sort {$a <=> $b} keys %umthash){
		#don't repeat
		next if ($done{$datapt}{$umt});
		print OUTX "$datapt";
		print OUTY "$umt";
		print OUTLOG "${datapt},${umt}\t";
		$done{$datapt}{$umt}++;
		foreach $header (sort {$a <=> $b} keys %allheaders){
				print OUTX "\t$datahash{$datapt}{$header}";
				print OUTY "\t$umthash{$umt}{$header}";
		}
		print OUTX "\n";
		print OUTY "\n";
	}
}

print OUTLOG "\n";

