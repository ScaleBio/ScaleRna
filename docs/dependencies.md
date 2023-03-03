# Dependency Management

The ScaleBio RNA workflow requires a number of dependencies to run. These include ScaleBio developed and open-source executables, python libraries, etc. There are three alternative ways to provide these dependencies; select one of these, depending on what is easiest on your system, and follow the instructions below.

## Using Docker or Singularity
If your system supports [docker containers](https://www.docker.com/), this is the recommended way to handle all dependencies for the ScaleBio RNA workflow. We provide pre-build docker containers and the workflow is setup to automatically use them.
This is enabled by adding `-profile docker` to the nextflow command-line.

If your system does not support *docker*, [singularity](https://sylabs.io/docs/) is an alternative that is enabled on many HPC clusters (2.3.x or newer). Setting `-profile docker,singularity` (**no space**) will use the _singularity_ engine for all dependencies. The environment variable `NXF_SINGULARITY_CACHEDIR` can be used to control where singularity images are stored. This should be a writable location that is available on all compute nodes. Similarly `TMPDIR` should be changed from the default `/tmp` to a location writable from the container if necessary.

See [Nextflow Containers](https://www.nextflow.io/docs/latest/container.html) for details and additional configuration options. 

One important point is that all input and output paths need to be available (_bind_) inside the containers. For _docker_, Nextflow will set the relevant options automatically at runtime; for _singularity_ this requires _user mounts_ to be enabled in the system-wide configuration (see the notes in the [Nextflow singularity documentation](https://www.nextflow.io/docs/latest/container.html#singularity)).

## Using Conda
Another option is using the [Conda](https://docs.conda.io/en/latest) package manager. Nextflow can automatically create conda environments with most dependencies. This mode is selected by setting `-profile conda`. In this case the following additional steps need to be completed:
- Install and update conda
    - `conda update -n base -c defaults conda`
- Install [ScaleBio Tools](scaleBioTools.md)
    - `/PATH/TO/ScaleRNA/envs/download-scale-tools.sh`
- If running from a sequencer runFolder (.bcls) Illumina [BCL Convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) is required to be installed (and available on `$PATH`)

See the [Nextflow documentation](https://www.nextflow.io/docs/latest/conda.html) for additional detail of conda support in Nextflow.

## Manual Dependency installation
As a final alternative it is also possible to simply install the required dependencies directly, either by hand or using Conda.
A list of all requirements can be found in `envs/scaleRna.conda.yml` and `envs/scalereport.conda.yml`. All tools need to be available on `$PATH` or in `/PATH/TO/ScaleRNA/bin/`
