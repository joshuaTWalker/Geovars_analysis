#!/bin/bash
module load picard
ls *.bam > temp.sh
sed -i -e 's/^/I\=/' temp.sh
sed -i '1ipicard MergeSamFiles' temp.sh
echo "O=output_all_bams.bam" >> temp.sh to run
sed -i -e ':a;N;$!ba;s/\n/ /g' temp.sh
bash temp.sh
rm temp.sh temp.she
