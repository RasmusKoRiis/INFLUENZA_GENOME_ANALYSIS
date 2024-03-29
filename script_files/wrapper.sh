#!/bin/bash
#Pipline for generation of Influenza A + B whole genome sequences, 
#including analysis producing mutation calling, vaccine-contamination protocols and more....

export INFLUENZA_V1_VERSION="1.0"
export LAB_PROTOCOL=""

#COMMAND LINE OPTIONS
while getopts ":i:dh" opt; do
  case ${opt} in
    i )
      run_folder=${OPTARG}
      ;;
    d )
      demultiplexing=TRUE
      ;;
    h )
      echo "Usage: ./master.sh -i [input_file] [-d]"
      echo ""
      echo "Options:"
      echo "i    Specify the input file (required)"
      echo "d    Specify if samples should be demultiplexed"
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

#SETS BASEDIR
basedir=$(pwd)

#DOWNLOAD SCRIPT
git clone https://github.com/RasmusKoRiis/INFLUENZA_GENOME_ANALYSIS.git

#CHECK INPUT FOLDER
if [ ! -d "$run_folder" ]; then
  echo "$run_folder does not exist, creating it and running rsync"
  ip_address='X'
  mkdir ${run_folder}_data
  cd ${run_folder}_data
  rsync -avr --exclude '*.fast5' grid@${ip_address}:/data/${run_folder}/* ./
  cd $basedir
  find $basedir -type d -name "fastq_pass" -exec cp -R {}/ $basedir \; -exec mv {} $basedir/${run_folder}_fastq \;
  cd $basedir
fi

#CHECK DEMULTIPLEX STATUS
if [ "$demultiplexing" == "TRUE" ]; then
  guppy_barcoder -i ${run_folder}_fastq/fastq_pass -s input_fastq --barcode_kits INFLUENSA --enable_trim_barcodes
  rm -r guppy_barcoder-core-dump-db
  python3 INFLUENZA_GENOME_ANALYSIS/script_files/rename_fastq_folders.py input_fastq *csv
  input_fastq=${basedir}/input_fastq
  echo "Demultiplexing and renaming done"
else 
  input_fastq=${basedir}/${run_folder}
fi

echo "The run folder is $input_fastq"

cp -a $input_fastq INFLUENZA_GENOME_ANALYSIS
cd INFLUENZA_GENOME_ANALYSIS
mv $input_fastq input_fastq_processed
input_fastq_processed=input_fastq_processed


#TECHNICAL INFO OF PIPLINE
date=$(date +"%Y-%m-%d_%H-%M-%S")
startdir=$(pwd)
script="$startdir/script_files"
runname=$(basename $startdir)
script_name=

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
reference="$startdir/references"

cd $startdir


#PRE POCESSING

# Build the QA-Docker image
image_name_qa="new_influensa_pipeline_qa_v1_3"
container_name_qa="influenza_qa_container_v0_1"
docker_file_qa="Dockerfile.QA"  

# Building the Docker image for QA
docker buildx build --platform linux/amd64 -t $image_name_qa -f Dockerfile.QA .

docker run --rm -it --name $container_name_qa \
  -v $startdir/qa_docker:/qa_docker \
  -e RUNNAME=$run_folder -e INFLUENZA_V1_VERSION=INFLUENZA_V1_VERSION \
  $image_name_qa bash -c "script_files/master_QA.sh && cp -r /app /qa_docker"

mkdir "$startdir/output_fastq_cleaned"
fastq_out_cleaned="$startdir/output_fastq_cleaned"

mv $startdir/qa_docker/app/output_fastq_cleaned/*fastq $fastq_out_cleaned

cd $fastq_out_cleaned

for file in *_cleaned.fastq; do
    # Remove '_cleaned' from the filename
    new_file="${file/_cleaned/}"

    # Create a directory named after the new file without '.fastq' extension
    mkdir "${new_file%.fastq}"

    # Rename the file and move it to the new directory
    mv "$file" "${new_file%.fastq}/$new_file"
done

cd $startdir


#CHECK IF IRMA SHOULD RUN

# EPI2ME NEXTFLOW 
nextflow run epi2me-labs/wf-flu -r v0.0.6 --fastq $fastq_out_cleaned/  --out_dir $result_folder/epi2me_wf_flu_output --min_qscore 10  --min_coverage 50 --reference "$startdir/references/epi2me/reference_epi2me_FULL_NAMES.fasta"

cd $startdir

container_name="influenza_container"
image_name="new_influensa_pipeline_v0.3"

docker buildx build --platform linux/amd64 -t $image_name .

# Run the Docker container to execute the rest of the pipeline and copy the results
docker run --rm -it --name $container_name \
  -v $startdir/results_docker:/results_docker \
  -e RUNNAME=$run_folder -e INFLUENZA_V1_VERSION=INFLUENZA_V1_VERSION  \
  $image_name bash -c "script_files/master_NF.sh && cp -r /app/results /results_docker"


#Copy and clean up folders
cd $basedir
cp -r $runname/results_docker/results $basedir
cp results/stats/*_summary.csv $basedir
#rm -r INFLUENZA_GENOME_ANALYSIS






