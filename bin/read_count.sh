for file in results/collapsed_reads/*; do
    filename=$(basename "$file")  # Extract file name
    awk -F '[-:]' '{ print $2 $3 }' "$file" | awk '{ print $2 "\t" $1 }' > "results/isomir_input/$filename"
done