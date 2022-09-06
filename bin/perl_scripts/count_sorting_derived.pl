#!/bin/perl
use strict;
use warnings;
#This script takes a vcf and finds the ancestral/derived alleles based on the outgroup. It then finds how many samples share the derived allele. It outputs counts per outgroup

my %outgroups;
$outgroups{"GRO_2043"}++;
$outgroups{"DIV_1956"}++;
$outgroups{"DEC_1895"}++;
$outgroups{"664647_GIG"}++;
my %sample;
my %counts;
my $site_counter = 1;
while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  if ($_ =~ /^##/){next;}
  if ($_ =~ /^#/){
    foreach my $i (9..$#a){
      $sample{$i} = $a[$i];
    }
  }else{

    $site_counter++;
    if ($site_counter % 50000 == 0){
      print STDERR ("Processing $a[0] $a[1]\n");
    }
    my %calls;
    my %all_bases;
    my $ancestors_genotyped = 0;
    #Load in genotypes
    foreach my $i (9..$#a){
      my @field = split(/:/,$a[$i]);
      my @genotype;
      if ($field[0] =~ /\|/){
        @genotype = split(/\|/,$field[0]);
      }else{
        @genotype = split(/\//,$field[0]);
      }
      if ($genotype[0] eq '.'){next;}
      if ($outgroups{$sample{$i}}){
        $ancestors_genotyped++;
      }
      foreach my $j (0..1){
        $calls{$sample{$i}}{$genotype[$j]}++;
        $all_bases{$genotype[$j]}++;
      }
    }
    #Check if there are any ancestors genotyped
    if ($ancestors_genotyped < 1){next;}
#    print ("$ancestors_genotyped genotyped\n");
    my $total_alleles = keys %all_bases;
    #Check which allele is ancestral. 
    my %ancestral_alleles;
    foreach my $ancestor (sort keys %outgroups){
      foreach my $alleles (sort keys %{$calls{$ancestor}}){
        $ancestral_alleles{$alleles}++;
      }
    }
    my $ancestral_alleles_count = keys %ancestral_alleles;
    #Skip sites where the outgroup has all the alleles
    if ($ancestral_alleles_count == $total_alleles){next;}

    #Find derived alleles;
    my %derived_alleles;
    foreach my $allele (sort keys %all_bases){
      unless($ancestral_alleles{$allele}){
        $derived_alleles{$allele}++;
#	print "At $a[0].$a[1] the derived allele is $allele\n";
      }
    }
    foreach my $derived_allele (sort keys %derived_alleles){
      my $sample_string;
      foreach my $i (9..$#a){
        if ($calls{$sample{$i}}{$derived_allele}){
          $sample_string .= "$sample{$i},";
        }
      }
      $counts{$sample_string}++;    
    }
  }
}
foreach my $string (sort keys %counts){
  print "$string\t$counts{$string}\n";
}

