#!/bin/perl
use strict;
use warnings;

print "#NEXUS\n\n";
print "BEGIN Trees;\n";
#Pipe in a file with tree on every line;
my $counter = 1;
while(<STDIN>){
  chomp;
  print "\nTree gt$counter=$_";
  $counter++;
}
print "\n\nEND;"
