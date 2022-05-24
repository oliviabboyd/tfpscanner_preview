
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tfpscanner : Transmission fitness polymorphism scanner

This package provides an analytical pipeline for rapid variant scanning
based on integration of scalable phylogenetic analysis with non-genetic
epidemiological data streams. The scanner computes a test statistic for
every possible bipartition of the provided phylogeny to identify
relative growth rates and relative evolutionary rates of clades with the
phylogeny. For every clade in a partitioned tree, matched comparison
clades are selected based on time (and if specified, space), for further
statistical analyses.

The main analyses conducted for each clade, or node with included
descendants, are as follows:

  - Molecular clock outlier statistic using root to tip regression
  - Simple logisitic growth rate estimate
  - Generalised additive model (GAM) combined with a Gaussian process
    model to estimate growth over time

Additional options in the scanning tool may be specified, including:

  - (GAM) combined wtih a model of spatial correlation between
    neighbouring regions to estimate frequency of a clade over time and
    space
  - Other covariate analyses, specified by the user, to run a
    conditional logistic regression on a specified covariate of interest
  - Identification of cluster defining mutations, identified as
    mutations present in \> x% of sequences in the specified cluster and
    less than x% of sequences in the comparison sample
  - MLEsky analyses to estime effective population size (Ne) over time

The output from the scanning tool includes a new output directory,
specified as `"tfpscan-{Sys.Date()}"`, along with a rds file containing
output descriptives (size, lineage, sample date) and statistics (clock
outlier, logistic growth rate, GAM) for every node included and a rds
file containing the environment from the scanning tool run. For each
node, logistic growth rate p-values as well as suport for logistic model
vs. GAM are also reported.

Within the output directory, the scanning tool creates a further folder
directory for each node included in analyses. Node speicifc outputs
within individual node directories include seperate csv files with the
following:

  - Summary statistics (clock outlier, logstic growth rate, GAM logistic
    growth rate, support values)
  - Sequences with specific sample times, region, and further
    descriptors reported if included (lineage, mutations, covariate)
  - Regional composition of sequences within individual nodes
  - Lineage composition of sequences within individual nodes
  - Co-circulating lineage composition from comparison matched sample
    nodes

Additional options for node specific outputs may be specified,
including:

  - GAM figure over time, plotting the logistics odds of a sequences
    being from the specified node
  - Map of GAM + neighbour joining model over space and time (in epi
    weeks)
  - Tree figure for the node of interest, colored by lineage
  - CSV file with defining mutations for cluster of interest specified
  - Online tree viewer for the whole phylogeny with linked hover
    function for statistics associated with each node

For a complete description of the statistical methodology underpinning
this package, see our preprint:
<https://www.biorxiv.org/content/10.1101/2021.01.18.427056v1.full>

## Installation

In R, install the `devtools` package and run

    devtools::install_github('mrc-ide/tfpscanner')

## Documentation

### Required

### 

See the vignettes in the R package for examples of how to use
treeviewer() and get\_mlesky\_node() functions.
