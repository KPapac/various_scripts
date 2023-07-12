set -eu -o pipefail

parse_arguments() {
  # Default values
  database=""
  query=""

  # Parse command-line arguments
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -db)
        shift
        database="$1"
        ;;
      -query)
        shift
        query="$1"
        ;;
      *)
        echo "Invalid argument: $1"
        exit 1
        ;;
    esac
    shift
  done

  # Check if both arguments were provided
  if [[ -z "$database" || -z "$query" ]]; then
    echo "Both -db and -query arguments are required."
    echo "USAGE: blast_to_art.sh -db <Path_to_database_of_Pseudogenes> -query <Path_to_FASTA>"
    exit 1
  fi
}

# Call the function with the command-line arguments
parse_arguments "$@"

## Running blast. First awk filters the blast file, second awk converts blastout to tab for artemis.
## Matches from blast that have 95% identity + 90% of the pseudogene covered in the alignment + an evalue less than 1e-5.
## Header: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle slen

blastn -db $database -query $query -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle slen" |\
  awk -F"\t" '($3>95 && $4/$14>0.9 && $11<0.00001) {print}' |\
awk -F"\t"\
  '{if ($9<$10) print "FT   misc_feature    "$7".."$8; 
            else print "FT   misc_feature    complement("$7".."$8")"} \
  {print "FT                   /color=6"}
  {print "FT                   /product="$13}
  {print "FT                   /note=PercID:"$3", eval:"$11", PercPseudoCov:"$4/$14}
  '
