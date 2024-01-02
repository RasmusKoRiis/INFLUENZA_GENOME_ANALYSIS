# Use an Ubuntu base image
#Build this docker on a mac M1/M2 with:
#docker buildx build --platform linux/amd64 -t new_influensa_pipline_v1 .
# RUN - 
FROM ubuntu:20.04

# Set environment variables for non-interactive installation
ENV DEBIAN_FRONTEND=noninteractive

# Install necessary dependencies
RUN apt-get update && \
    apt-get install -y wget openjdk-11-jre-headless bash python3 python3-pip git build-essential zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev cmake

# Install Python packages
#RUN pip3 install pandas matplotlib seaborn biopython pysam

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh && \
    bash Miniconda3-py39_4.10.3-Linux-x86_64.sh -b -p /miniconda && \
    rm Miniconda3-py39_4.10.3-Linux-x86_64.sh
ENV PATH="/miniconda/bin:${PATH}"
RUN conda update -y conda

# Install Python packages
RUN conda install -y -c conda-forge pandas
RUN conda install -y -c conda-forge matplotlib
RUN conda install -y -c conda-forge seaborn
RUN conda install -y -c conda-forge biopython=1.81
RUN conda install -y -c bioconda pysam

# Install nextalign
RUN conda install -c bioconda -y nextalign

# Install Nextflow
# Download and install the specified version of Nextflow
RUN wget -qO- https://get.nextflow.io/ | bash \
    && mv nextflow /usr/local/bin/


# Install seqkit
RUN wget https://github.com/shenwei356/seqkit/releases/download/v0.16.1/seqkit_linux_amd64.tar.gz && \
    tar -xvzf seqkit_linux_amd64.tar.gz && \
    mv seqkit /usr/local/bin/ && \
    rm seqkit_linux_amd64.tar.gz

# Install weeSAM
RUN git clone https://github.com/centre-for-virus-research/weeSAM.git && \
    chmod +x weeSAM/weeSAM && \
    mv weeSAM/weeSAM /usr/local/bin/ && \
    rm -rf weeSAM

# Install Nextclade
RUN wget https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-unknown-linux-gnu && \
    chmod +x nextclade-x86_64-unknown-linux-gnu && \
    mv nextclade-x86_64-unknown-linux-gnu /usr/local/bin/nextclade

# Install bam-readcount
RUN git clone https://github.com/genome/bam-readcount.git && \
    cd bam-readcount && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    mv bin/bam-readcount /usr/local/bin/ && \
    cd ../.. && \
    rm -rf bin/bam-readcount



# Set the working directory
WORKDIR /app

# Copy the script and reference files into the container
COPY script_files /app/script_files
COPY references /app/references
COPY dataset /app/dataset
COPY results /app/results


# Make the script executable
RUN chmod +x /app/script_files/master_NF.sh 

# Set the default command to run the script
CMD ["/app/script_files/master_NF.sh"]

