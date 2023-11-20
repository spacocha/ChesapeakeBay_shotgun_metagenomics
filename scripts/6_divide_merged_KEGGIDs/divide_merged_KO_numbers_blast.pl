#! /usr/bin/perl -w
#
#

	die "Usage: table old_KO ox_blast red_blast2> Redirect\n" unless (@ARGV);
($table, $oldKEGG, $oxblast, $redblast) = (@ARGV);
chomp ($table, $oldKEGG, $oxblast, $redblast);

$first=1;
open (IN, "<$oxblast") or die "Can't open $oxblast\n";
while ($line=<IN>){
    chomp ($line);
    next unless ($line);
    (@pieces)=split("\t", $line);
    if ($first){
	$first=0;
    } else {
    	$oxhash{$pieces[0]}++;
    }
}
close (IN);

$first=1;
open (IN, "<$redblast") or die "Can't open $redblast\n";
while ($line=<IN>){
    chomp ($line);
    next unless ($line);
    (@pieces)=split("\t", $line);
    if ($first){
	$first=0;
    } else {
    	if ($oxhash{$pieces[0]}){
		die "In both $pieces[0]\n";
    	} else {
    		$redhash{$pieces[0]}++;
    	}
    }
}
close (IN);

open (IN, "<$table") or die "Can't open $table\n";
while ($line=<IN>){
    chomp ($line);
    next unless ($line);
    ($node, @pieces)=split("\t", $line);
    ($conf) = pop(@pieces);
    ($KEDDtaxID)=pop(@pieces);
    ($KEGGgen)=pop(@pieces);
    ($KEGGphy)=pop(@pieces);
    ($KEGGdom)=pop(@pieces);
    ($KEGGID2)=pop(@pieces);
    ($KEGGID1)=pop(@pieces);
    if (@headers){
     	#die "$node $KEGGID1\n";
	if ($KEGGID1 eq $oldKEGG){
		if ($oxhash{$node}){
			#die "Oxidative $node\n";
			$i=0;
			$j=@pieces;
			until ($i >=$j){
				$oxsum{$headers[$i]}+=$pieces[$i];
				$i++;
			}
		} elsif ($redhash{$node}){
			#die "Reductive $node\n";
                        $i=0;
                        $j=@pieces;
                        until ($i >=$j){
                                $redsum{$headers[$i]}+=$pieces[$i];
                                $i++;
                        }
		} else {
			#die "Neither $node\n";
                        $i=0;
                        $j=@pieces;
                        until ($i >=$j){
                                $nonesum{$headers[$i]}+=$pieces[$i];
                                $i++;
                        }
		}
	}
    } else {
	(@headers)=@pieces;
	foreach $head (@headers){
		$allheaders{$head}++;
	}
    }
}

print "KO";
foreach $head (sort keys %allheaders){
	print "\t$head";
}
print "\n${oldKEGG}.1";
foreach $head (sort keys %allheaders){
	if ($oxsum{$head}){
		print "\t$oxsum{$head}";
	} else {
		print "\t0";
	}
}

print "\n${oldKEGG}.2";
foreach $head (sort keys %allheaders){
	if ($redsum{$head}){ 
		print "\t$redsum{$head}";
	} else {
		print "\t0";
	}
}

print "\nNA";
foreach $head (sort keys %allheaders){
	if ($nonesum{$head}){ 
		print "\t$nonesum{$head}";
	} else {
		print "\t0";
	}
}

print "\n";
