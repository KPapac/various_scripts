set -euox pipefail 

# Chromosomal assemblies
for gbk in annotations_to_submit/{wPauO11,wProPE2,wWilP98}.gbk
do

basename=$(basename -s .gbk $gbk)
chr_file=${basename}_chrfile
chr_name=$(head -n1 $gbk | awk '{print $2}')

# Sets up chromosome file
echo "${chr_name} 1 Circular-Chromosome" > $chr_file

# Sets up manifest and chromosome files
manifest_file="${basename}_manifest"
study="PRJEB76506"
togrep=${basename:4}
sample_code=$(grep "$togrep" sample_IDs | cut -f2 -d',')
assembly_name=${basename}_chromosome
program="SPAdes genome assembler v3.15.5"
flatfile="${basename}.embl"

# Make the EMBL flat file
python3 fix_gbk_for_ENA.py $gbk $flatfile -f circular
gzip $flatfile $chr_file

cat <<EOF > $manifest_file
STUDY           $study
SAMPLE          $sample_code
ASSEMBLYNAME    $assembly_name
ASSEMBLY_TYPE   isolate
COVERAGE        10
PROGRAM         $program
PLATFORM        TODO
MINGAPLENGTH    1
MOLECULETYPE    genomic DNA
FLATFILE        ${flatfile}.gz
CHROMOSOME_LIST ${chr_file}.gz
EOF

done

# Draft genomes now!
for gbk in $(ls annotations_to_submit/ -I "wProPE2*" -I "wWilP98*" -I "wPauO11*")
do

gbk=annotations_to_submit/$gbk
basename=$(basename -s .gbk $gbk)

# Sets up manifest and chromosome files
manifest_file="${basename}_manifest"
study="PRJEB76506"
togrep=${basename:4}
sample_code=$(grep "$togrep" sample_IDs | cut -f2 -d',')
assembly_name=${basename}_chromosome
program="SPAdes genome assembler v3.15.5"
flatfile="${basename}.embl"

# Make the EMBL flat file
python3 fix_gbk_for_ENA.py $gbk $flatfile -f linear 
gzip $flatfile 

cat <<EOF > $manifest_file
STUDY           $study
SAMPLE          $sample_code
ASSEMBLYNAME    $assembly_name
ASSEMBLY_TYPE   isolate
COVERAGE        10
PROGRAM         $program
PLATFORM        TODO
MINGAPLENGTH    1
MOLECULETYPE    genomic DNA
FLATFILE        ${flatfile}.gz
EOF
done
