#!/bin/bash

listfile="CEEHRC.CEMT0004.txt"
outdir="/path/to/output"

mkdir -p "$outdir"

# 找 Input 文件
input_file=$(grep 'Input.signal_unstranded.bin1000.step1000.bp$' "$listfile")

if [ -z "$input_file" ]; then
    echo "No Input file found!"
    exit 1
fi

# 遍历 signal 文件
grep 'signal_unstranded.bin1000.step1000.bp$' "$listfile" | \
grep -v 'Input.signal_unstranded.bin1000.step1000.bp$' | \
while read signal_file; do

    base=$(basename "$signal_file")

    outfile="${outdir}/${base%.bp}.divInput.bp"

    echo "Processing:"
    echo "  Signal: $signal_file"
    echo "  Input : $input_file"
    echo "  Output: $outfile"

    paste "$signal_file" "$input_file" | \
    awk 'BEGIN{OFS="\t"}
    {
        signal=$4
        input=$8

        if(input==0){
            ratio="NA"
        } else {
            ratio=signal/input
        }

        print $1,$2,$3,ratio
    }' > "$outfile"

done