#!/usr/bin/perl
use strict;

my %h = ();
my $count = '';

while ( <> ) {
  chomp;
  $h{$_}++;
}

foreach ( sort { $h{$a} <=> $h{$b} } keys %h ) {
	print "$_\t$h{$_}\n";
}


