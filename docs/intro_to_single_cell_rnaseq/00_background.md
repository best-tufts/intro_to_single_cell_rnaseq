## Introduction To Single-Cell RNA-Seq Bioinformatics

**Content developed by [Data Intensive Studies Center](https://disc.tufts.edu/)**

- Rebecca Batorsky, PhD, Data Scientist
- Eric Reed, PhD, Data Scientist 

!!! info ""

## Schedule

**Lecture**

1. [Quality Control, Integration, Clustering](slides/lecture_part_1.pdf)
2. [Cell-type identification and Differential Expression](slides/lecture_part_2.pdf)

**Hands On Activity**

1. [Setup](01_setup.md) 
2. [Quality Control](02_quality_control.md)
3. [Integration](03_integration.md)
4. [Clustering](04_clustering.md)
5. [Cell-type identification](05_cell_type_identification.md)
6. [Differential Expression](06_differential_expression.md)


!!! example "Prerequisites"
    - [Request an account](http://research.uit.tufts.edu/) on the Tufts HPC Cluster
        - Note if you signed up for the Introduction to Single-Cell RNA-Seq workshop this will have been already taken care of for you!
    - Connect to the [VPN](https://access.tufts.edu/vpn) if off campus
    

## Navigate To The Cluster

Once you have an account and are connected to the VPN/Tufts Network, navigate to the [OnDemand Website](https://ondemand.pax.tufts.edu/){:target="_blank" rel="noopener"} and log in with your tufts credentials. Once you are logged in you'll notice a few navigation options:

!!! info "OnDemand Layout"

    ![](images/ondemand_layout_pic.png)

Click on `Interactive Apps > RStudio Pax` and you will see a form to fill out to request compute resources to use RStudio on the Tufts HPC cluster. We will fill out the form with the following entries:

- `Number of hours` : `5`
- `Number of cores` : `1`
- `Amount of memory` : `16GB`
- `R version` : `4.0.0`
- `Reservation for class, training, workshop` : `Bioinformatics Workshops`
    - **NOTE: This reservation will closed on Nov 5th 2023, use `Default` if running through the materials after that date.**

Click `Launch` and wait until your session is ready. Click `Connect To RStudio Server`, and you will notice a new window will pop up with RStudio. 

??? question "Are you connected to RStudio?"
    - Yes (put up a green check mark in zoom)
    - No (raise hand in zoom)


## Data & Scripts

To copy over the data and scripts we will need for the workshop into our home directory, enter the following command into the console:

```R
file.copy(from="/cluster/tufts/bio/tools/training/intro_to_scrnaseq",to="~/", recursive = TRUE)
```

## Project Setup
Now we are going to use this folder to create a new R project. R projects are great for managing analyses in a portable, self-contained folder. To create an R project from within our `intro_to_scrnaseq` directory we will:

- Go to `File` > `New Project`
- `Existing Directory`
- Browse for the `intro_to_scrnaseq` folder
- Click `Create Project`

Let's navigate to our project in our home directory and open up our workshop script:

- Click on the `Files` tab in the lower right hand Rstudio pane
- Click on the `scripts` folder
- Click on the `intro_to_scrnaseq.Rmd` script
