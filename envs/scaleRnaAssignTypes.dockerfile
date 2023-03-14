FROM nfcore/base:2.1

RUN apt-get clean && apt-get update && \
 DEBIAN_FRONTEND=noninteractive apt-get -y install \
 build-essential \
 curl libz-dev \
 unzip libhdf5-dev && \
 rm -rf /var/lib/apt/lists/*

COPY assignTypes.conda.yml /environment.yml

RUN conda env create --quiet -f /environment.yml && \
 conda clean -a --yes && \
 conda env export --name scaleRnaAssignTypesNf > scaleRnaAssignTypesNf.yml

ENV PATH="/opt/conda/envs/nfRnaAnalysis/bin:${PATH}"

COPY installGithubDependenciesForCellTyping.R /installGithubDependencies.R

RUN Rscript /installGithubDependencies.R 
