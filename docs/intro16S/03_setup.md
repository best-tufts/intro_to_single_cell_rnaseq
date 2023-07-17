
!!! example "Prerequisites"
    - [Request an account](http://research.uit.tufts.edu/) on the Tufts HPC Cluster
    - Connect to the [VPN](https://access.tufts.edu/vpn) if off campus

## Navigate To The Cluster

Once you have an account and are connected to the VPN/Tufts Network, navigate to the [OnDemand Website](https://ondemand.pax.tufts.edu/){:target="_blank" rel="noopener"} and log in with your tufts credentials. Once you are logged in you'll notice a few navigation options:

![](images/ondemandLayout.png)

Click on `Interactive Apps > RStudio Pax` and you will see a form to fill out to request compute resources to use RStudio on the Tufts HPC cluster. We will fill out the form with the following entries:

- `Number of hours` : `3`
- `Number of cores` : `1`
- `Amount of memory` : `8GB`
- `R version` : `4.0.0`
- `Reservation for class, training, workshop` : `Bioinformatics Workshop`---> NOTE: This reservation closed on Nov 9, 2022, use Default if running through the materials after that date.
- `Load Supporting Modules`: `boost/1.63.0-python3 java/1.8.0_60 gsl/2.6`

Click `Lauch` and wait until your session is ready. Click `Connect To RStudio Server`, and you will notice a new window will pop up with RStudio. 

## Project Setup

We are going to create a new project to begin:

1. Go to `File` > `New Project`
2. `New Directory`
3. `New Project`
4. Create a name for your project (e.g. `intro-to-16S`)
5. `Create Project`

## File Organization

In our project we will need some folders to contain our scripts, data and results:

- Click the New Folder icon
- Create a folder called data and click ok
- Following the same process, create a scripts folder and a results folder


## Data & Scripts

Today we will be working with data from Rosshart et al. (2107) where wild-type and laboratory strain mouse microbiomes were assessed. To copy over 
this data we will enter the following command into the console:

```R
file.copy(from="/cluster/tufts/bio/tools/training/microbiome16S/raw_fastq/",to="./data/", recursive = TRUE)
file.copy(from="/cluster/tufts/bio/tools/training/microbiome16S/meta/metaData.txt",to="./data/", recursive = TRUE)
file.copy(from="/cluster/tufts/bio/tools/training/microbiome16S/silva/silva_nr99_v138.1_train_set.fa.gz",to="./data/")
file.copy(from="/cluster/tufts/bio/tools/training/microbiome16S/scripts/dada2pipeline.Rmd",to="./scripts/")
```

Now that we have our data and scripts copied, let's navigate to our scripts folder and open up "dada2pipeline.Rmd".

## Libraries

To run a code chunk in this R markdown file, click the play button at the top right hand side of the code chunk. We will practice by running the code chunk that loads the R libraries we will need for this workshop:

```R
# load our libraries
.libPaths(c('/cluster/tufts/hpc/tools/R/4.0.0',.libPaths()))
library(dada2)
library(phyloseq)
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(phangorn)
library(msa)
```
