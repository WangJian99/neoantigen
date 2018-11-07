#!/usr/bin/perl
use strict;
use warnings;

open IN,$ARGV[0] or die $!;

my %hash;
#my $head = <IN>;
#print "$head";

while (my $line=<IN>){
	chomp $line;
	next if $line =~ /^ORDER_ID/;
	my @line=split(/\t/,$line);
	my $id=join("_",@line[0..4,7]);
	my $aff=$line[9];
	next if $aff eq "-";
	#print "$id\n";
	$hash{$id}{$aff}=$line;

}
close IN;

foreach my $k (keys %hash){
	foreach my $m (sort {$a <=> $b} keys %{$hash{$k}}){
		print "$hash{$k}{$m}\n";
		last;
	}

}
