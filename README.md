# fungiflow
Workflow for identifying biosynthetic gene clusters from raw sequence data with minimal inputs.

This repository encompasses an entire workflow from raw sequence data up to and including biosynthetic gene cluster identification.

## Installation
Clone this GitHub repositroy by entering ```git clone ``` in the directory you would like to install in.

This workflow is designed to operate in parallel on an HPC, so a lot of cpus and memory is required. If you are running this on a personal computer, I would suggest a minimum of 8 cpus and 40 GB of memory, as some of the programs in the workflow require a lot of memory.

To make this workflow minimally difficult to operate, all the software needed has been packaged into a series of Singularity images, meaning that the only pre-requisite you need for the workflow to function other than a bash environment is Singularity and its prequisites.

- bash environment (tested with Ubuntu 18.02)
- Singularity >/ 3.5.3
