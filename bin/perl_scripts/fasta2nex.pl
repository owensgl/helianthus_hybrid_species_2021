#!/bin/perl
use strict;
use warnings;

my $phylonet_suffix = "TRUE"; #Turn this off to prevent the extra MCMC_SEQ tag at the end;
my $species_string = "<ann:ANN1029,ANN1283;arg:ARG0143,ARG0295;pet_pet:PET0568,PET0495;pet_fal:PET0765,PET0424;niv:PET0695,PET0662;deb:DEB_1837,DEB_1135;ano:HK24_16_01_W_097,NHKG_16_01_R_196;des:GO7_GB071,GO8_GB179;par:PAR_3,PAR_posas_01;outgroup:DIV_1956,DEC_1895,664647_GIG,GRO_2043>"; #Specific set of samples to species for phylonet;
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
print "#NEXUS\n\n";
print "BEGIN data;\n";
#Figure out dimensions of the data
my @loci = sort keys %data;
my @samples = sort keys %{$data{$loci[0]}};
my $n_samples = @samples;
my $total_characters;
foreach my $locus (@loci){
  my $loci_length = length($data{$locus}{$samples[0]});
  $total_characters += $loci_length;
}
print "Dimensions ntax=$n_samples nchar=$total_characters;\n";
print 'Format datatype=dna symbols="ACTG" missing=N gap=-;';
print "\n";
print "Matrix\n";
foreach my $locus (@loci){
  my $length = length($data{$locus}{$samples[0]});
  print "[$locus, $length]\n";
  foreach my $sample (@samples){
    print "$sample ";
    print "$data{$locus}{$sample}\n";
  }
}
print ";\nEnd;\n";
if ($phylonet_suffix eq "TRUE"){
  print "BEGIN PHYLONET;\n";
  print "MCMC_SEQ -loci (";  
  foreach my $locus (@loci){
    print "$locus";
    if ($locus ne $loci[$#loci]){
      print ",";
    }
  }
  print ")";
  print " -cl 2000000 -bl 200000 -sf 10000 ";
  print "-tm $species_string;\n";
  print "END;";
}
