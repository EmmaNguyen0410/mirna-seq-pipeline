#!/bin/bash

# Loop over all matching .tabular files
for file in *.fastq.fasta.gz.tabular; do
  # Remove the .fastq.fasta.gz and .tabular parts, add .txt
  base_name=$(basename "$file" .fastq.fasta.gz.tabular)
  new_txt="${base_name}.txt"

  # Copy original to .txt version
  cp "$file" "$new_txt"

  # Run the awk pipeline
  awk -F '[-:]' '{ print $2 $3 }' "$new_txt" | awk '{ print $2 "\t" $1 }' > "$new_txt"
done
