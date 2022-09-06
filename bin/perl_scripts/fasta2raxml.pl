#!/bin/perl
use warnings;
use strict;
#This takes a list of fasta files, and merges them into a raxml file with partitions file, for use in quartet sampling
#Pipe in a list of fasta, all with same samples and same length sequence
my %sequences;
my $random_sample;
my $partition_string;
my $running_basecount = 0;
my $prefix = $ARGV[0];
my $gene_counter = 1;
while(<STDIN>){
  chomp;
  my $file = $_;
  my $counter= 0;
  my $current_name;
  open FILE, $file;
  while(my $line = <FILE>){
    chomp $line;
    if ($counter % 2 == 0){
      $line =~ s/>//g;
      $current_name = $line;
      $random_sample = $current_name;
    }else{
      $sequences{$current_name}.=$line;
    }
    if ($counter == 1){
      my $start = $running_basecount+1;
      my $end = $running_basecount + length($line);
      my $name = "p$gene_counter";
      $gene_counter++;
      $partition_string .= "DNA, $name=$start-$end\n";
      $running_basecount+=length($line);
    }
    $counter++;
  }
}
open(RAXML, '>', "$prefix.phy") or die $!;
open(PARTITIONS, '>', "$prefix.partitions") or die $!;

print PARTITIONS "$partition_string";
my $n_sequences = scalar keys %sequences;
my $sequences_length = length($sequences{$random_sample});
print RAXML " $n_sequences $sequences_length\n";
foreach my $sample (sort keys %sequences){
  print RAXML "$sample $sequences{$sample}\n";
}
