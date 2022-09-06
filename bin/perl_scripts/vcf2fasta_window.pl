#!/bin/perl
#This script takes a vcf file and outputs a fasta.
use strict;
use warnings;
my $IUPAC = "TRUE";
my $window_size = 50000;
my $missing_character = "N"; #Normally N, but ? for splitstree
my(%table) = (
        'AC' => 'M',
        'CA' => 'M',
        'AG' => 'R',
        'GA' => 'R',
        'AT' => 'W',
        'TA' => 'W',
        'CG' => 'S',
        'GC' => 'S',
        'CT' => 'Y',
        'TC' => 'Y',
        'TG' => 'K',
        'GT' => 'K',
        'AA' => 'A',
        'CC' => 'C',
        'GG' => 'G',
        'TT' => 'T',
        'NN' => 'N',
	'N' => 'N',
);

my %ind;
my %sites;
my $counter = 0;
my $start_position;
my $final_position;
while(<STDIN>){
  chomp;
  if ($_ =~ m/^##/){next;}
  my @a = split(/\t/,$_);
  if ($_ =~ m/^#/){
    foreach my $i (9..$#a){
      $ind{$i} = $a[$i];
    }
  }else{
    $counter++;
    my $chr = $a[0];
    my $pos = $a[1];
    if ($counter == 1){
      $start_position = "$chr-$pos";
    }
    $final_position = $pos;
    my $ref = $a[3];
    #Only take first base of ref to fix vcf issue
    my @refs = split(//,$ref);
    $ref = $refs[0];
    my @alts = split(/,/,$a[4]);
    foreach my $i(9..$#a){
      my @tmp = split(/:/,$a[$i]);
      my $call = $tmp[0];
      $call =~ s/\|/\//;
      my @bases = split(/\//,$call);
      my $current_call;
      if (($call eq '\.') or ($call eq "\.\/")){
	$current_call = $missing_character;
        next;
      }
      foreach my $j (0..1){
        unless(exists($bases[$j])){
          $current_call .= $missing_character;
	  next;
	}
        if ($bases[$j] eq "0"){
          $current_call .= $ref;
        }elsif($bases[$j] eq "1"){
          $current_call .= $alts[0];
	}elsif($bases[$j] eq "2"){
	  $current_call .= $alts[1];
	}elsif($bases[$j] eq "3"){
	  $current_call .= $alts[2];
        }elsif($bases[$j] eq '.'){
	  $current_call .= $missing_character;
	}
       
      }
      if ($IUPAC eq "TRUE"){
        my $code = $table{$current_call};
        unless($code){
	  print "$current_call\n";
	}
        $sites{$ind{$i}} .= $code;
      }else{
        my @calls = split(//,$current_call);
        my $rand = int(rand(2));
	$sites{$ind{$i}} .= $calls[$rand];
      }
    }
  }
}
my $output_name = "$start_position-$final_position.fasta";
open(OUT, '>', $output_name);
&print_fasta;

sub print_fasta {
  my @a = sort values %ind;
  foreach my $i (0..$#a){
    print  OUT ">$a[$i]\n";
    print  OUT "$sites{$a[$i]}\n";
  }
}

