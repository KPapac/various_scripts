#ls -d Sample_* | while read line; do 
#echo $line
#cd $line
#ls *_vs_wPau*.bam

bam_file=$1
depth_file=$(basename -s .bam $1)_depth
samtools depth $bam_file > $depth_file
cov_from_depth_any_length.pl $depth_file 100 > ${depth_file}_100
R --slave --args ${depth_file}_100 < /home/kostas/scripts/plotCoverage_mod.R
#cd ..
#done 
