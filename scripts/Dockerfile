FROM ubuntu:latest

# Set the maintainer label
LABEL maintainer="manish@south.du.ac.in"

# Update the package list and install dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    python3-numpy \
    python3-pandas \
    python3-matplotlib \
    git \
    build-essential \
    wget \
    unzip \
    bzip2 \
    default-jdk \
    default-jre \
    mash \
    skesa \
    seqtk \
    spades \
    bowtie2 \
    samtools \
    racon \
    pilon \
    libz-dev \
    perl \
    libfindbin-libs-perl \
    openjdk-11-jre-headless \
    ncbi-blast+ \
    bedtools \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    libreadline-dev \
    libsqlite3-dev \
    llvm \
    libncurses5-dev \
    libncursesw5-dev \
    xz-utils \
    tk-dev \
    libffi-dev \
    liblzma-dev \
    curl \
    hmmer \
    augustus \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

#RUN apt-get install -y bbmap
# Download and install BBMap
RUN wget https://sourceforge.net/projects/bbmap/files/latest/download -O bbmap.tar.gz && \
    mkdir /opt/bbmap && \
    tar -xvzf bbmap.tar.gz -C /opt/bbmap --strip-components=1 && \
    rm bbmap.tar.gz
# Set BBMap path
ENV PATH="/opt/bbmap:$PATH"

# Create a virtual environment and install MultiQC within it
RUN python3 -m venv /opt/new_venv && \
    /opt/new_venv/bin/pip install --upgrade pip && \
    /opt/new_venv/bin/pip install multiqc

# Add the virtual environment to the PATH
ENV PATH="/opt/new_venv/bin:$PATH"

# Install FLASH
RUN wget http://ccb.jhu.edu/software/FLASH/FLASH-1.2.11.tar.gz -O /opt/FLASH.tar.gz && \
    tar -xzf /opt/FLASH.tar.gz -C /opt && \
    cd /opt/FLASH-1.2.11 && \
    make && \
    cp flash /usr/local/bin/

# Install Python 3.7 from source # i changed here <- from python3.7 to 3.8
RUN mkdir /opt/python && \
    curl -L https://www.python.org/ftp/python/3.7.12/Python-3.7.12.tgz | tar -xz -C /opt/python --strip-components=1 && \
    cd /opt/python && \
    ./configure --enable-optimizations && \
    make -j$(nproc) && \
    make altinstall

# Verify Python installation and SSL support
RUN /usr/local/bin/python3.7 -c "import ssl; print(ssl.OPENSSL_VERSION)"

# Use ensurepip to install setuptools and pip for Python 3.7
RUN /usr/local/bin/python3.7 -m ensurepip --upgrade

# Install pip for Python 3.7
RUN /usr/local/bin/python3.7 -m pip install --no-cache --upgrade pip setuptools

# Install Prodigal from GitHub
RUN git clone https://github.com/hyattpd/Prodigal.git /opt/prodigal && \
    cd /opt/prodigal && \
    make && \
    cp prodigal /usr/local/bin/prodigal

# Install MetaEuk using precompiled binaries (AVX2 or SSE4.1)
RUN wget https://mmseqs.com/metaeuk/metaeuk-linux-avx2.tar.gz -O /opt/metaeuk.tar.gz && \
    tar -xzvf /opt/metaeuk.tar.gz -C /opt/ && \
    rm /opt/metaeuk.tar.gz && \
    ln -s /opt/metaeuk/bin/metaeuk /usr/local/bin/metaeuk

# Install specific Python libraries
RUN /usr/local/bin/python3.7 -m pip install numpy==1.17.4 pandas==1.3.5 matplotlib==3.5.0 biopython==1.80 requests==2.31.0

# Install ggplot
RUN /usr/local/bin/python3.7 -m pip install ggplot==0.11.5

# Create a virtual environment using Python 3.7 and activate it
RUN /usr/local/bin/python3.7 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Upgrade pip and setuptools inside the virtual environment
RUN pip install --upgrade pip setuptools

# Install FastQC
RUN wget -O /opt/fastqc_v0.11.9.zip https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip /opt/fastqc_v0.11.9.zip -d /opt && \
    chmod +x /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc

# Clone the Unicycler repository and install it
RUN git clone https://github.com/rrwick/Unicycler.git /opt/Unicycler && \
    cd /opt/Unicycler && \
    python3 -m pip install .

# Install Trimmomatic
RUN wget -O /opt/Trimmomatic-0.39.zip http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip /opt/Trimmomatic-0.39.zip -d /opt && \
    chmod +x /opt/Trimmomatic-0.39/trimmomatic-0.39.jar && \
    ln -s /opt/Trimmomatic-0.39/trimmomatic-0.39.jar /usr/local/bin/trimmomatic

# Set the environment variable for Trimmomatic
ENV TRIMMOMATIC_JAR=/opt/Trimmomatic-0.39/trimmomatic-0.39.jar

# Install PlentyofBugs from GitHub
RUN python3 -m venv /opt/plentyofbugs_env && \
    /opt/plentyofbugs_env/bin/pip install --upgrade pip && \
    git clone --depth 1 https://github.com/nickp60/plentyofbugs.git /opt/PlentyofBugs && \
    cd /opt/PlentyofBugs && \
    /opt/plentyofbugs_env/bin/pip install .

# Add PlentyofBugs to the PATH
ENV PATH="/opt/plentyofbugs_env/bin:$PATH"

# Install MUMmer
RUN wget https://downloads.sourceforge.net/project/mummer/mummer/3.23/MUMmer3.23.tar.gz && \
    tar -xzvf MUMmer3.23.tar.gz && \
    cd MUMmer3.23 && \
    make && \
    cp -r * /usr/local/bin && \
    cd .. && \
    rm -rf MUMmer3.23 MUMmer3.23.tar.gz

# Clone the QUAST repository and check out the desired version
RUN git clone https://github.com/ablab/quast.git && \
    cd quast && \
    git checkout tags/quast_5.2.0 -b quast_5.2.0 && \
    ./install_full.sh || (echo "MetaQUAST installation failed, but proceeding with lightweight QUAST installation") && \
    ln -s /quast/quast.py /usr/local/bin/quast

# Install BUSCO
RUN git clone https://gitlab.com/ezlab/busco.git /opt/busco && \
    cd /opt/busco && \
    /usr/local/bin/python3.7 -m pip install numpy pandas matplotlib biopython && \
    /usr/local/bin/python3.7 -m pip install .

# Install necessary dependencies
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    git \
    zlib1g-dev \
    gdb \
    && rm -rf /var/lib/apt/lists/*

# Clone the AlignGraph repository
RUN git clone https://github.com/baoe/AlignGraph.git /opt/AlignGraph

# Set the working directory to AlignGraph
WORKDIR /opt/AlignGraph/AlignGraph

# Add AlignGraph directory to PATH
ENV PATH="/opt/AlignGraph/AlignGraph:$PATH"

# Install blat
# Set the working directory
WORKDIR /opt/

# Download precompiled BLAT binaries
# Download and install precompiled BLAT binary
RUN wget -O /opt/blat \
    https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat \
    && chmod +x /opt/blat \
    && cp /opt/blat /usr/local/bin/blat \
    && /usr/local/bin/blat 2>&1 | head -n 2 || true

#RUN wget --no-check-certificate \
#    -O /opt/blatSuite.34.zip \
#    https://hgwdev.cse.ucsc.edu/~kent/exe/opteron/blatSuite.34.zip

# Extract BLAT binaries
#RUN unzip /opt/blatSuite.34.zip -d /opt/ && \
#    find /opt -type f -name blat -exec cp -v {} /usr/local/bin/blat \; && \
#    chmod +x /usr/local/bin/blat

# Add /usr/local/bin to PATH
ENV PATH="/usr/local/bin:$PATH"

# Install blat
# Download precompiled BLAT binaries
# RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat \
#    && wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/gfClient \
#    && wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/gfServer

# Move BLAT binaries to /usr/local/bin and give them execute permissions
#RUN chmod +x blat gfClient gfServer \
#    && mv blat gfClient gfServer /usr/local/bin/ \
#    && ln -s /usr/local/bin/blat /usr/bin/blat \
#    && ln -s /usr/local/bin/gfClient /usr/bin/gfClient \
#    && ln -s /usr/local/bin/gfServer /usr/bin/gfServer

# Add BLAT to PATH so it can be used by AlignGraph
#ENV PATH="/usr/local/bin:$PATH"

# install blat using bioconda
# Create a virtual environment using Python 3.7
#RUN /usr/bin/python3 -m venv /opt/blat_env

# Set PATH for the virtual environment
#ENV PATH="/opt/blat_env/bin:$PATH"

# Download and install Miniconda
#RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
#    chmod +x Miniconda3-latest-Linux-x86_64.sh && \
#    ./Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
#    rm Miniconda3-latest-Linux-x86_64.sh && \
#    /opt/conda/bin/conda clean -afy

# Install Conda in the newly created virtual environment and install BLAT
#RUN /opt/conda/bin/conda create -n blat_env python=3.8 -y && \
#    /opt/conda/bin/conda run -n blat_env conda install bioconda::blat -y && \
#    ln -s /opt/conda/envs/blat_env/bin/blat /usr/local/bin/blat && \
#    /opt/conda/bin/conda clean -afy

# Update PATH to include Conda binaries
#ENV PATH="/opt/conda/bin:$PATH"

#Install MUMer with compiled dependencies
# Install required packages
#RUN apt-get update && apt-get install -y \
#    build-essential \
#    g++ \
#    make \
#    git \
#    perl \
#    wget \
#    csh \
#    likwid \
#    && apt-get clean \
#    && rm -rf /var/lib/apt/lists/*
#
# Clone the MUMmer repository
#RUN git clone https://github.com/garviz/MUMmer.git /opt/MUMmer

#COPY Makefiles from the local machine to the Docker image
#COPY MUMmer/Makefile /opt/MUMmer/
#COPY MUMmer/src/kurtz/Makefile /opt/MUMmer/src/kurtz/
#COPY MUMmer/src/kurtz/streesrc/Makefile /opt/MUMmer/src/kurtz/streesrc/
#COPY MUMmer/src/kurtz/libbasedir/Makefile /opt/MUMmer/src/kurtz/libbasedir/
#COPY MUMmer/src/tigr/Makefile /opt/MUMmer/src/tigr/

# Install MUMmer
#RUN cd /opt/MUMmer && make check && make clean && make install

# Add MUMmer binary directory to PATH
#ENV PATH="/opt/MUMmer:${PATH}"

# Install biopython
# Download precompiled BLAT binaries
#RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat \
#    && wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/gfClient \
#    && wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/gfServer

# Move BLAT binaries to /usr/local/bin and give them execute permissions
#RUN chmod +x blat gfClient gfServer \
#    && mv blat gfClient gfServer /usr/local/bin/

# Set the working directory
WORKDIR /data

# Copy the genome assembly script into the Docker image
COPY dnrbg.sh /usr/local/bin/
COPY pad_reads_1.py /usr/local/bin/

# Make the scripts executable
RUN chmod +x /usr/local/bin/dnrbg.sh /usr/local/bin/pad_reads_1.py

# Ensure /data is created as a volume mount point
VOLUME ["/data"]

# Set the entrypoint to the script and allow passing of input arguments
ENTRYPOINT ["/usr/local/bin/dnrbg.sh"]
