# Dependency Management

The ScaleBio Seq Suite: RNA workflow requires a number of dependencies to run. These include ScaleBio developed and open-source executables, python libraries, etc. Nextflow provides different options to automatically provide these dependencies. Pick one of these approaches and select the corresponding Nextflow configuration `-profile`, or install dependencies manually. 

## Using Docker or Singularity
If your system supports [docker containers](https://www.docker.com/), this is the recommended way to handle all dependencies for the ScaleBio RNA workflow. We provide pre-built docker containers and the workflow is setup to automatically use them.
This is enabled by adding `-profile docker` to the nextflow command-line.

If your system does not support *docker*, [singularity](https://sylabs.io/docs/) (or new versions called *Apptainer*) is an alternative that is enabled on many HPC clusters (2.3.x or newer). Setting `-profile singularity` will use the _singularity_ engine for all dependencies. The environment variable `NXF_SINGULARITY_CACHEDIR` can be used to control where singularity images are stored. This should be a writable location that is available on all compute nodes. Similarly `TMPDIR` should be changed from the default `/tmp` to a location writable from the container if necessary.

See [Nextflow Containers](https://www.nextflow.io/docs/latest/container.html) for details and additional configuration options. 

One important point is that all input and output paths need to be available (_bind_) inside the containers. For _docker_, Nextflow will set the relevant options automatically at runtime; for _singularity_ this requires _user mounts_ to be enabled in the system-wide configuration (see the notes in the [Nextflow singularity documentation](https://www.nextflow.io/docs/latest/container.html#singularity)).

## Using Conda
Another option is using the [Conda](https://docs.conda.io/en/latest) package manager. Nextflow can automatically create conda environments with most dependencies. This mode is selected by setting `-profile conda`. 

In this case the following additional steps need to be completed:
- Install and update conda
    - `conda update -n base -c defaults conda`
- Install [ScaleBio Tools](scaleBioTools.md)
    - `/PATH/TO/ScaleRNA/envs/download-scale-tools.sh`
- If running from a sequencer runFolder (.bcls) Illumina [BCL Convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) is required to be installed (and available on `$PATH`)
- If nextflow throws an error installing packages while using `-profile conda`, there is a more verbose yaml file with a comprehensive list of all python packages and their versions, specified in [envs](../envs) ([scalerna_verbose.conda.yml](../envs/scalerna_verbose.conda.yml) and [scalernareport_verbose.conda.yml](../envs/scalernareport_verbose.conda.yml))
    - To run nextflow with these conda files, edit the [nextflow.config](../nextflow.config) and replace the process.conda section of the conda profile with the appropriate yml. So `process.conda = "$projectDir/envs/scalerna.conda.yml"` becomes `process.conda = "$projectDir/envs/scalerna_verbose.conda.yml` and `conda = "$projectDir/envs/scalernareport.conda.yml"` becomes `conda = "$projectDir/envs/scalernareport_verbose.conda.yml"`

See the [Nextflow documentation](https://www.nextflow.io/docs/latest/conda.html) for additional detail of conda support in Nextflow. Generally automatic installation will work best if the `base` conda environment is clean, i.e. does not contain extra channels or complex packages.

## Manual Dependency installation
As a final alternative it is also possible to simply install the required dependencies directly, either by hand or using Conda.
A list of all requirements can be found in `envs/scalerna.conda.yml` and `envs/scalernareport.conda.yml`. All tools need to be available on `$PATH` or in `/PATH/TO/ScaleRNA/bin/`
