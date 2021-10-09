#! /usr/bin/perl -w
#
#

	die "Usage: KOdb mat> Redirect\n" unless (@ARGV);
($list, $mat) = (@ARGV);
chomp ($list, $mat);

open (IN1, "<$list") or die "Can't open $list\n";
while ($line1=<IN1>){
	chomp($line1);
	next unless ($line1);
	open (IN, "<${line1}") or die "Can't open $line1\n";
		while ($line=<IN>){
			chomp ($line);
			next unless ($line);
			($user, @rest)=split ("\t", $line);
        		($gene)=$user=~/^user:(.+)$/;
        		#die "$user $gene\n";
			$allhash{$gene}=$line;
	}
}
close (IN);

open (IN, "<${mat}") or die "Can't open $mat\n";
while ($line=<IN>){
    chomp ($line);
    next unless ($line);
    ($gene)=split ("\t",  $line);
    print "$line";
    if ($allhash{$gene}){
	($user, @rest)=split ("\t", $allhash{$gene});
	foreach $piece (@rest){
		print "\t$piece";
	}
	print "\n";
   } else {
	print "\tNA\tNA\tNA\tNA\tNA\tNA\n";
   }

}
close (IN);
