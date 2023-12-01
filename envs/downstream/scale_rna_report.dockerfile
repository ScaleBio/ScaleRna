FROM rocker/verse

RUN apt-get clean && apt-get update && \
DEBIAN_FRONTEND=noninteractive apt-get -y install \
build-essential \
curl libz-dev \
unzip libhdf5-dev && \
rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c("argparse", "Seurat"))'
