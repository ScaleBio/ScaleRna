FROM continuumio/miniconda3

RUN apt-get clean && apt-get update && \
 DEBIAN_FRONTEND=noninteractive apt-get -y install \
 build-essential \
 curl libz-dev \
 unzip && \
rm -rf /var/lib/apt/lists/*

# Install the conda environment
COPY scalereport.conda.yml /environment.yml
RUN conda env create --quiet -f /environment.yml && \
 conda clean -a --yes && \
 conda env export --name scalereport > /scalereport.yml
# Instead of 'conda activate'
ENV PATH="/opt/conda/envs/scalereport/bin:${PATH}"
