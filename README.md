# Nextflow FLAIRR pipeline based on pipeAIRR

This is a re-work of the [the Yaari lab pipeAIRR script](https://github.com/yaarilab/FLAIRRSeq) written by Ayelet Peres. 


## Installation

The pipeline will run on Linux, or Windows Subsystem for Linux (WSL).

Beyond this repo, you will need to install:
- [Nextflow](https://www.nextflow.io/)
- [Docker client](https://www.docker.com/) or [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html)
- Some Python modules (see python/requirements.txt)

The examples assume Docker is used.

## Usage

1. Review the configuration in processed_samples/FLAIRRSeq.config and make any changes to reflect location of your reference sets, number of threads, etc.
2. Review the script in processed_samples/process_samples.sh and make changes to the wildcards and directories, to suit the naming of your samples and FASTQ files.
3. Edit processed_samples/flairr_logs.toml and modify the first line to reflect the absolute location of the processed_samples directory in your file system (this toml file is not used by nextflow, but it is used by subsequent python scripts to produce summary analyses).
4. Run the pipeline with the following command:
```bash
cd processed_samples
./process_samples.sh
```
The script shows both the use of 'standard' germline reference sets and the use of personalised reference sets for each sample, which could, for example,
be taken from genomic sources. See the params at the top of preprocess/main.nf and annotate/main.nf for configuration, eg to annotate different loci. Any params
setting can be over-ridden on the command line, as shown in process_samples.sh.

## Output
Output is created in subdirectories under processed_samples. The script as written retains the nextflow work directories ./.nextflow, in both processed_samples itself and in the subdirectories. These work directories can grow large, and should be deleted unless needed for debugging. 


## Running with WSL

The pipeline uses the Immcantation Docker image. The current image contains an installation of usearch which is not compatible with WSL. To use the pipeline in WSL, you will need to build a new Docker image, as follows:

1. Run the container interactively: 
```bash
docker run -it --cpus 1 immcantation/suite:4.4.0 /bin/bash
```
2. In another terminal, find the container ID:
```bash
docker ps
```
3. In the container, type:
```bash
mkdir muscledir  
cd muscledir  
wget https://www.drive5.com/muscle/muscle_src_3.8.1551.tar.gz  
tar xzvf muscle_src_3.8.1551.tar.gz  
make
```
The last stage of the build (ln) will fail. Repeat the ln command, but without the -static option. Then test that muscle has build correctly:
```bash
./muscle --help
```
Move it to the bin directory and clean up:
```bash
mv muscle /usr/local/bin  
cd ..  
rm -rf muscledir
```
4. Exit the container and commit the changes:
```bash
docker commit <container_id> my_immcantation_suite:4.4.0
```
5. Modify processed_samples/FLAIRRSeq.config to use my_immcantation_suite:4.4.0 rather than immcantation_suite:4.4.0

