FROM nfcore/base:2.1

RUN apt-get clean && apt-get update && \
DEBIAN_FRONTEND=noninteractive apt-get -y install \
build-essential \
curl libz-dev \
unzip libhdf5-dev && \
rm -rf /var/lib/apt/lists/*

COPY scrubletEnvironment.yml /environment.yml

RUN conda install -n base conda-libmamba-solver && \
    conda env create --quiet -f /environment.yml --experimental-solver libmamba && \
    conda clean -a --yes && \
    conda env export --name scrublet > scrublet.yml

ENV PATH="/opt/conda/envs/scrublet/bin:${PATH}"
