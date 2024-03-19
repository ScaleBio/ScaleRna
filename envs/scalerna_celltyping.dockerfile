FROM nfcore/base:2.1
# Path to yaml file containing conda dependencies
ARG CONDA_YML_PATH=.
RUN apt-get clean && apt-get update && \
 DEBIAN_FRONTEND=noninteractive apt-get -y install \
 build-essential \
 curl libz-dev \
 unzip && \
rm -rf /var/lib/apt/lists/*

# Install the conda environment
COPY $CONDA_YML_PATH/scalerna_celltyping.conda.yml /environment.yml
RUN conda install -n base conda-libmamba-solver
RUN conda env create --experimental-solver libmamba --quiet -f /environment.yml && \
 conda clean -a --yes && \
 conda env export --name scalerna_celltyping > scalerna_celltyping.yml
# Instead of 'conda activate'
ENV PATH="/opt/conda/envs/scalerna_celltyping/bin:${PATH}"

COPY $CONDA_YML_PATH/installGithubPackages.r ./installGithubPackages.r

RUN ./installGithubPackages.r
