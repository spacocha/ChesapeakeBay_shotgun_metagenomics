#! /usr/bin/perl -w

die "lnmetabolicgenesnoatp: gene x sample matrix (103 x 38)
konoatp: KO names (1 x 103)
MATLAB_metadata_3: export from regression_models_3.m (38 x 27)
allphoto: ln of all photosynthetic genes (run regression_model_3.m in MATLAB)
exclude: kovalues to exclude
output_prefix: The beginning part of the file names
Usage: lnmetabolicgenesnoatp allphoto konoatp exclude_list MATLAB_metadata_3 output_prefix\n" unless (@ARGV);
($mat, $allphotofile, $konoatp, $excludelist, $metadata, $prefix) = (@ARGV);
chomp ($mat);
chomp ($allphotofile);
chomp ($konoatp);
chomp ($excludelist);
chomp ($metadata);
chomp ($prefix);

$lineno=1;
open (IN, "<$konoatp") or die "Can't open $konoatp\n";
while ($line = <IN>){
	chomp ($line);
	next unless ($line);
	#ko num corresponds to the order of genes in the mat 
	$kohash{$lineno}=$line;
	$lineno++;
}
close (IN);

open (IN, "<$excludelist") or die "Can't open $excludelist\n";
while ($line = <IN>){
        chomp ($line);
        next unless ($line);
	($name,$exko)=split("\t",$line);
        #ko num corresponds to the order of genes in the mat
        $excludehash{$exko}++;
}
close (IN);

#lnmetabolicgenesnoatp has no header or gene info
#arrayed as gene x library (38 columns for each library
#103 rows for each gene
$lineno=1;
open (IN, "<$mat" ) or die "Can't open $mat\n";
while ($line =<IN>){
    chomp ($line);
    next unless ($line);
    (@pieces)=split(",", $line);
    #die "$pieces[0]\n";
    $i=0;
    $j=@pieces;
    until ($i>=$j){
		#lineno corresponds to gene names in %kohash
		#$i+1 related to number 1...47 by indenv
		$ei=$i+1;
		if ($kohash{$lineno}){
			#It has a translated KO
			#Is it in exclude list
			unless ($excludehash{$kohash{$lineno}}){
				$mathash{$kohash{$lineno}}{$ei}=$pieces[$i];
				$allheaders{$ei}++;
				#die "ei $ei lineno $lineno pieces $pieces[$i]\n" if ($lineno == 2 && $ei == 2); 
			}
		} else {
			die "Missing $lineno $kohash{$lineno}\n";
		}
		$i++;
    }
    $lineno++;
} 
close (IN);

$lineno=1;
#lineno equals library no
open (IN, "<$allphotofile" ) or die "Can't open $allphotofile\n";
while ($line =<IN>){
    chomp ($line);
    next unless ($line);
    (@pieces)=split(",", $line);
    #die "$pieces[0]\n";
    #add this to mathash
    $mathash{'allphoto'}{$lineno}=$pieces[0];
    $allheaders{$lineno}++;
    #die "ei $ei lineno $lineno pieces $pieces[$i]\n" if ($lineno == 2 && $ei == 2);
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
print OUTY "Gene";
foreach $header (sort {$a <=> $b} keys %allheaders){
	print OUTX "\t$header";
	print OUTY "\t$header";
}
print OUTX "\n";
print OUTY "\n";

#make a matrix of values to compare
foreach $datapt (sort keys %datahash){
	foreach $KOID (sort {$a <=> $b} keys %mathash){
		#don't repeat
		next if ($done{$datapt}{$KOID});
		print OUTX "$datapt";
		print OUTY "$KOID";
		print OUTLOG "${datapt},${KOID}\t";
		$done{$datapt}{$KOID}++;
		foreach $header (sort {$a <=> $b} keys %allheaders){
				print OUTX "\t$datahash{$datapt}{$header}";
				print OUTY "\t$mathash{$KOID}{$header}";
		}
		print OUTX "\n";
		print OUTY "\n";
	}
}

print OUTLOG "\n";

