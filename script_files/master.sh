#!/bin/bash

#Pipline for generation of Influenza A + B whole genome sequences, 
#including analysis producing mutation calling, vaccine-contamination protocols and more....

#!/bin/bash

#COMMAND LINE OPTIONS
while getopts ":i:h" opt; do
  case ${opt} in
    i )
      run_folder=${OPTARG}
      ;;
    h )
      echo "Usage: ./master.sh -i [input_file]"
      echo ""
      echo "Options:"
      echo "i    Specify the input file (required)"
      echo "h    Display this help message"
      exit 0
      ;;
    \? )
      echo "Invalid option: -$OPTARG" 1>&2
      exit 1
      ;;
    : )
      echo "Option -$OPTARG requires an argument." 1>&2
      exit 1
      ;;
  esac
done

if [ -z "$run_folder" ]; then
  echo "Error: You must specify an input file with the -i option"
  exit 1
fi

echo "The run folder is $run_folder"


#TECHNICAL INFO OF PIPLINE
date=$(date +"%Y-%m-%d_%H-%M-%S")
startdir=$(pwd)
script="$startdir/script_files"
demultiplexed="not determined"
runname=$(basename $startdir)

#FOLDER FOR RESULTS AND SEQUENCES
mkdir results
cd results
mkdir fasta bam stat mutation

cd $startdir

result_folder="$startdir/results"
fasta_folder="$result_folder/fasta"
bam_folder="$result_folder/bam"
stat_folder="$result_folder/stat"
mutation_folder="$result_folder/mutation"
dataset_folder="$startdir/dataset"

#TRANSFER OF RAW SEQUENCINF DATA FROM GRIDION
source "$script"/config.sh
#run_folder=${1}

if [ -d "$run_folder" ]; then
  echo "Data is located, starts analysis"
  cd "$run_folder"
else
  echo "Transfering data"
  mkdir "$run_folder"
  cd "$run_folder"
  rsync -avr --exclude '*.fast5' "grid@$ip_address:/data/$run_folder/*" ./
  echo "Transfer complete, starts analysis"
fi

#CHECKS DEMULTIPLEX STATUS

if [ -d "$(find . -type d -name 'demultiplexed' -print -quit)" ]; then
  echo "Data is demultiplexed"
  demultiplexed="TRUE_DP"
  export demultiplexed_folder=$(find . -type d -name "demultiplexed" -print -quit)
elif [ -n "$(find . -type d -name '*barcode[0-9][0-9]' -print -quit)" ]; then
  echo "Barcode directory found"
  demultiplexed="TRUE_FP"
  export fastq_pass_folder=$(find . -type d -name "fastq_pass" -print -quit)
else
  echo "Data is not demultiplexed"
  demultiplexed="FALSE"
    export fastq_pass_folder=$(find . -type d -name "fastq_pass" -print -quit)
fi


#DEMULIPLEXING IF NECESSARY

# REMOVE PRIMERS


#EPI2ME NEXTFLOW - CONSENSUS AND SUBTYPING ANALYSIS
if [[ "$demultiplexed" == "TRUE_FP" ]]; then
  input_fastq="$fastq_pass_folder"
elif [[ "$demultiplexed" == "TRUE_DP" ]]; then
  input_fastq="$demultiplexed_folder"
else
  echo "Error: Unknown demultiplexed status"
  exit 1
fi

nextflow run epi2me-labs/wf-flu --fastq $input_fastq --out_dir $result_folder/epi2me_wf_flu_output --min_qscore 14  --min_coverage 50 --downsample 600 --reference "$startdir/references/epi2me/reference_epi2me_FULL_NAMES.fasta"

#PLANNED NEXTFLOW IMPLEMENTATION

#GENERATION OF STATISTICS
#add removal of "-" from fasta sequences
cd $startdir
"$script"/technical_stat.sh

#MUTATION GENERATION
cd $startdir
"$script"/mutation_list.sh

#H1N1 resitance
cd $startdir
"$script"/H1N1_resistens_checker.sh

#PRIMER CHECKER
cd $startdir
#"$script"/primer_checker.sh

#NOISE FINDER

#FINALIZIN§G SUMMARY
python3 "$script/group_sample_summary.py" "$stat_folder/"$runname"_long_quality_mutation_summary.csv" "$result_folder"/"$runname"_final.csv

