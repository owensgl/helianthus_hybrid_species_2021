#!/bin/perl
use strict;
use warnings;

while(<STDIN>){
  chomp;
  if ($. == 1){
    print "$_";
    next;
  }
  if ($_ =~ m/^#/){
    print "\n$_";
    next;
  }
  my @a = split(/\t/,$_);
  my $alt = $a[3];
  my $length = length($alt);
  if ($length > 1){
    next;
  }
  print "\n$_";
}
