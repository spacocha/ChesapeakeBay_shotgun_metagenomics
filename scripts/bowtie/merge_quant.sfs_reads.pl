#! /usr/bin/perl -w
#
#


	die "Usage: KO list> Redirect\n" unless (@ARGV);
($KOfile, $list) = (@ARGV);
chomp ($KOfile, $list);

open (IN, "<${KOfile}") or die "Can't open $KOfile\n";
while ($line=<IN>){
	chomp ($line);
	next unless ($line);
	($name, $KO, $L1, $L2, $L3)=split ("\t", $line);
	$transhash{$name}{"KO"}=$KO;
	$transhash{$name}{"L1"}=$L1;
	$transhash{$name}{"L2"}=$L2;
	$transhash{$name}{"L3"}=$L3;
}
close (IN);


open (IN, "<${list}") or die "Can't open $list\n";
while ($file=<IN>){
    chomp ($file);
    next unless ($file);
    ($prefix)=$file=~/\/.*(HJKKLBCX2_[12]_[A-Z]{8}~[A-Z]{8})/;
    die "$file $prefix\n" unless ($prefix);
    open (IN2, "<${file}") or die "Can't open $file\n";
    while ($line=<IN2>){
	chomp ($line);
	next unless ($line);
	next if ($line=~/^Name/);
	($name, $length, $efflen, $TPM, $reads)=split("\t", $line);
	#die "$name $length $efflen $TPM $reads\n";
	if ($transhash{$name}{"KO"}){
		$hash{$transhash{$name}{"KO"}}{$prefix}{"reads"}+=$reads;
    	}
    }
    close (IN2);
}

close (IN);
print "Read";
foreach $name (sort keys %hash){
        foreach $prefix (sort keys %{$hash{$name}}){
                print "\t$prefix";
        }
	print "\n";
	last;
}

#Next I have to figure out a way to print out the entire matrix
foreach $name (sort keys %hash){
	print "$name";
	foreach $prefix (sort keys %{$hash{$name}}){
		print "\t$hash{$name}{$prefix}{'reads'}";
	}
	print "\n";
}

