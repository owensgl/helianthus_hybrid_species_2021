i=$1
bcftools view -r Ha412HOChr$i gowens22.mappable.variant.combined.dp4.missing80.vcf.gz | grep -v "^##" > chromosomes/chr$i.vcf
  cd chromosomes
  tail -n +2 chr$i.vcf | split -l 10000 - split_chr${i}_
  head -n 1 chr$i.vcf > tmp_chr${i}_header;
  echo "Splitting VCF $i"
  for file in split_chr${i}_*; 
  do     
    cat tmp_chr${i}_header "$file" > ${file}_labelled;     
    rm $file;
  done
  echo "Making fasta $i"
  for file in split_chr${i}_*;
  do
    cat "$file" | perl /home/owens/bin/reformat/vcf2fasta_window.pl; 
    rm "$file";
  done
  rm chr$i.vcf;
  cd ..
