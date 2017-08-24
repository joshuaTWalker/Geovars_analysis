#!/bin/bash
for f in ./*.bam
do
sed -i -e 's/S\_/S/g' $f
done
