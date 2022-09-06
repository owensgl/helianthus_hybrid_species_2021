#!/bin/perl
use warnings;
use strict;
#Takes a bunch of fasta files and turns them into a VCF. Assumes that VCFs are haploid.
my $gene_location_file = "gene_locations.txt"; 
my $gene_name_suffix = ".sorted.fa.aligned.2.fas.gblocks-gb.fasta.tree.subset.fasta";
my $gene_file_path = "/media/owens/Childs/helianthus_phylogeny/sequence_capture/Helianthus_fasta/";
open GENE, $gene_location_file;

my %genes;
my %gene_chr;
my %gene_start;
while(<GENE>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $start = $a[1];
  my $gene = $a[3];
  $genes{$gene}++;
  $gene_chr{$gene} = $chr;
  $gene_start{$gene} = $start;
}
close GENE;
my %sequence_data;
my %samples;
my %chrs;

foreach my $gene (sort keys %genes){
  my $gene_file = "$gene_file_path$gene$gene_name_suffix";
  open FILE, $gene_file;
print STDERR "Loading $gene\n";
  my $current_sample;
  my $current_chr = $gene_chr{$gene};
  my $original_start = $gene_start{$gene};
  my $current_position = $original_start;
  $chrs{$current_chr}{$original_start}++;
  while(<FILE>){
    chomp;
    if ($_ =~ m/^>/){ #It's a sample name
      my $sample = $_;
      $sample =~ s/>//g;
      $current_sample = $sample;
      $samples{$sample}++;
      $current_position = $original_start;
    }else{
      my @bases = split(//,$_);
      foreach my $base (@bases){
        $sequence_data{$current_chr}{$current_position}{$current_sample} = $base;
        $current_position++;
      }
    }
  }
}
#Now print out a VCF #CRAZY;
print "##fileformat=VCFv?\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
my $last_species;
foreach my $species (sort keys %samples){
  print "\t$species";
  $last_species = $species;
}
foreach my $chr (sort keys %chrs){
  foreach my $start (sort {$a <=> $b} keys %{$chrs{$chr}}){
    my $position = $start;
    while(exists($sequence_data{$chr}{$position})){
      #Count the number of alleles;
      my %alleles;
      foreach my $sample (sort keys %samples){
        unless($sequence_data{$chr}{$position}{$sample}){
	  $sequence_data{$chr}{$position}{$sample} = "-";
        }
        if ($sequence_data{$chr}{$position}{$sample} eq "-"){next;}
        $alleles{$sequence_data{$chr}{$position}{$sample}}++;
      }
      my @alleles = sort keys %alleles;
      my $ref = $alleles[0];;
      
      my $ref = $alleles[0];
      my $alt = ".";
      if ($alleles[1]){
        $alt = $alleles[1];
      }
      if ($alleles[2]){
        $alt .= ",$alleles[2]";
      }
      if ($alleles[3]){
        $alt .= ",$alleles[3]";
      }
      unless($ref){
        $ref = ".";
      }
      print "\n$chr\t$position\t.\t$ref\t$alt\t.\t.\t.\t.";
      foreach my $sample (sort keys %samples){
        if ($sequence_data{$chr}{$position}{$sample} eq "-"){
          print "\t./.";
        }elsif($sequence_data{$chr}{$position}{$sample} eq $ref){
          print "\t0/0";
        }elsif($sequence_data{$chr}{$position}{$sample} eq $alleles[1]){
          print "\t1/1";
        }elsif($sequence_data{$chr}{$position}{$sample} eq $alleles[2]){
          print "\t2/2";
        }elsif($sequence_data{$chr}{$position}{$sample} eq $alleles[3]){
          print "\t3/3";
        }  
      }
      $position++;
    }

  }
}

