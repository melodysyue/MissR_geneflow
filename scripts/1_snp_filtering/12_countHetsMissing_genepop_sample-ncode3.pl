#/usr/bin/perl -w
use strict;

my$header=<>;
print "Sample\tTotalLoci\tMissingGenos\tHeterozygous\tMissingPerc\tGenoRate\tHetPerc\n";
my@loci;
my$locusSection=1;
while(my$line=<>){
	chomp $line;
	if($line=~/^pop/i){
		$locusSection=0;
		next;
	}
	if($locusSection==1){
		push@loci,$line;
	}else{
		my$missingGenos=0;
		my$hetGenos=0;
		my$totalLoci=0;
		my@columns=split " ", $line;
		print "$columns[0]";
		foreach my$i (1..$#columns){
			$totalLoci++;
			my$allele1=substr $columns[$i], 0, 3;
			my$allele2=substr $columns[$i], 3, 3;
			if($allele1 eq "000"){
				$missingGenos++;
			}elsif($allele1 ne $allele2){
				$hetGenos++;
			}
		}
		my$missingPerc=$missingGenos/$totalLoci;
		my$genoRate=1-$missingPerc;
		my$hetPerc=$hetGenos/$totalLoci;
		print "\t$totalLoci\t$missingGenos\t$hetGenos\t$missingPerc\t$genoRate\t$hetPerc\n";
	}
}
