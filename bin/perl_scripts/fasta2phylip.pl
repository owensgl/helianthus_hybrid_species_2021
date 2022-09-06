#!/bin/perl
use strict;
use warnings;

#Pipe in a file with tree on every line;
my $counter = 1;
my %data;
#Pipe in a list of files;
while(<STDIN>){
  chomp;
  my $loci = "loci_$counter";
  $counter++;
  my $filename = $_;
  open (FILE, '<', $filename);
  my $line_counter = 0;
  my $sample_buffer;
  while(my $line = <FILE>){
    chomp $line;
    if ($line_counter % 2 == 0){
      #It's a header line
      my $sample = $line;
      $sample =~ s/\>//;
      $sample_buffer = $sample;
    }else{
      $data{$loci}{$sample_buffer} = $line;
    }
    $line_counter++;
  }
}
#Figure out dimensions of the data
my @loci = sort keys %data;
my @samples = sort keys %{$data{$loci[0]}};
my $n_samples = @samples;
my $total_characters;
foreach my $locus (@loci){
  my $loci_length = length($data{$locus}{$samples[0]});
  $total_characters += $loci_length;
}
foreach my $locus (@loci){
  my $length = length($data{$locus}{$samples[0]});
  print "$n_samples $length\n";
  foreach my $sample (@samples){
    print "$locus^$sample\t";
    print "$data{$locus}{$sample}\n";
  }
}
