# DNA methylation pipelines

This repository contains pipelines to process DNA methylation data for RRBS and WGBS experiments. It does adapter trimming, mapping, methylation calling, and produces other outputs. You can download the latest version from the [releases page](https://github.com/databio/dnameth_pipelines/releases) and a history of version changes is in the [CHANGELOG](CHANGELOG.md).

## Pipeline features at-a-glance

These features are explained in more detail later in this README.

Description pending.

## Installing

**Prerequisite python packages**. This pipeline uses [pypiper](https://github.com/databio/pypiper) to run a single sample, [looper](https://github.com/pepkit/looper) to handle multi-sample projects (for either local or cluster computation), and [pararead](https://github.com/databio/pararead) for parallel processing sequence reads. You can do a user-specific install of these like this:

```
pip install --user https://github.com/databio/pypiper/zipball/master
pip install --user https://github.com/pepkit/looper/zipball/master
pip install --user https://github.com/databio/pararead/zipball/master
```

**Required executables**. You will need some common bioinformatics tools installed. The list is specified in the pipeline configuration files (`.yaml` files in [src/](src/)).

**Genome resources**. This pipeline requires genome assemblies produced by [refgenie](https://github.com/databio/refgenie). You may [download pre-indexed references](http://cloud.databio.org/refgenomes) or you may index your own (see [refgenie instructions](https://github.com/databio/refgenie#indexing-your-own-reference-genome)). Any prealignments you want to do use will also require refgenie assemblies. Some common examples are provided by [ref_decoy](https://github.com/databio/ref_decoy).

**Clone the pipeline**. Clone this repository using one of these methods:
- using SSH: `git clone git@github.com:databio/dnameth_pipelines.git`
- using HTTPS: `git clone https://github.com/databio/dnameth_pipelines.git`

## Configuring

There are two configuration options: You can either set up environment variables to fit the default configuration, or change the configuration file to fit your environment. Choose one:

**Option 1: Default configuration** (recommended; e.g. [src/rrbs.yaml](src/rrbs.yaml)). 
  - Make sure the executable tools (java, samtools, bowtie2, etc.) are in your PATH.
  - Set up environment variables to point to `jar` files for the java tools (`picard` and `trimmomatic`).
  ```
  export PICARD="/path/to/picard.jar"
  export TRIMMOMATIC="/path/to/trimmomatic.jar"
  ```
  
  - Define environment variable `GENOMES` for refgenie genomes. 
  ```
  export GENOMES="/path/to/genomes/folder/"
  ```
  

**Option 2: Custom configuration**. Instead, you can also put absolute paths to each tool or resource in the configuration file to fit your local setup. Just change the pipeline configuration file ([src/rrbs.yaml](src/rrbs.yaml)) appropriately. 


## Running the pipeline

You never need to interface with the pipeline directly, but you can if you want. Just run `python src/rrbs.py -h` to see usage. But the best way to use this pipeline is to run it using looper. You will need to tell looper about your project. Example project data are in the [examples/test_project](examples/test_project) folder. Run the pipeline across all samples in the test project with this command:
```
looper run examples/test_project/test_config.yaml
```

If the looper executable in not your `$PATH`, add the following line to your `.bashrc` or `.profile`:

```
export PATH=$PATH:~/.local/bin
```

Now, adapt the example project to your project. Here's a quick start: You need to build two files for your project (follow examples in the [examples/test_project](examples/test_project/) folder):

- [project config file](examples/test_config.yaml) -- describes output locations, pointers to data, etc.
- [sample annotation file](examples/microtest_annotation.csv) -- comma-separated value (CSV) list of your samples.

Your annotation file must specify these columns:
- sample_name
- library (must be 'RRBS' or 'WGBS')
- organism (may be 'human' or 'mouse')
- read1
- read2
- whatever else you want

Run your project as above, by passing your project config file to `looper run`. More detailed instructions and advanced options for how to define your project are in the [Looper documentation on defining a project](http://looper.readthedocs.io/en/latest/define-your-project.html). Of particular interest may be the section on [using looper derived columns](http://looper.readthedocs.io/en/latest/advanced.html#pointing-to-flexible-data-with-derived-columns).

## Using a cluster

Once you've specified your project to work with this pipeline, you will also inherit all the power of looper for your project.  You can submit these jobs to a cluster with a simple change to your configuration file. Follow instructions in [configuring looper to use a cluster](http://looper.readthedocs.io/en/latest/cluster-computing.html).

Looper can also summarize your results, monitor your runs, clean intermediate files to save disk space, and more. You can find additional details on what you can do with this in the [looper docs](http://looper.readthedocs.io/). 

## Contributing

Pull requests welcome. Active development should occur in a development or feature branch.
