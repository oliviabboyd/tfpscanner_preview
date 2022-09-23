

#dplyr version

## scanner testing
library( ape)
library( data.table )
library( lubridate )
library( glue ) 
library( mgcv )
library( ggplot2 )
library( ggtree ) 
library( phangorn )
library( rbenchmark )

#load trees and amds 

tre <- read.tree("C:/Users/oboyd/Documents/COGUK_server/sp_analyses/tree_archive/cog_global_2021-12-31_tree.newick")

amd <- read.csv("C:/Users/oboyd/Documents/amd_phe_merged_2022-03-25.csv")

amd1 <- fread("C:/Users/oboyd/Documents/amd_phe_merged_2022-03-25.csv")


# fread vs benchmark testing ####

#benchmark("readcsv" = {
#  amd <- read.csv("C:/Users/oboyd/Documents/amd_phe_merged_2022-03-25.csv")
#},
#"fread" = {
#  amd1 <- fread("C:/Users/oboyd/Documents/amd_phe_merged_2022-03-25.csv")
#},
#replications = 3,
#columns = c("test", "replications", "elapsed",
#            "relative", "user.self", "sys.self"))

#     test replications elapsed relative user.self sys.self
#2   fread            3   27.26    1.000     36.31     4.06
#1 readcsv            3  680.91   24.978    669.21     3.87

min_descendants = 50
max_descendants = 20e3
min_cluster_age_yrs = 1/12
min_date = "2021-09-01"
max_date = "2021-12-14"
min_blen = 1/30e3/2
ncpu = 1
output_dir = paste0('tfpscan-', Sys.Date())
num_ancestor_comparison = 500
factor_geo_comparison = 5
Tg = 7/365
report_freq = 50 
mutation_cluster_frequency_threshold = 0.75
test_cluster_odds = c() 
test_cluster_odds_value = c() 
root_on_tip = 'Wuhan/WH04/2020'
root_on_tip_sample_time = 2020 
detailed_output = FALSE 
compute_gam = FALSE
compute_geo = FALSE
compute_cluster_muts = FALSE
compute_cluster_tree = FALSE



if (!dir.exists( output_dir ))
  dir.create( output_dir )

max_time <- Inf 
if ( !is.null( max_date )){
  max_time <- decimal_date(ymd( max_date ))
} else{
  max_date <- Sys.Date()
}

min_time <- -Inf 
if (!is.null( min_date ))
  min_time <- decimal_date( ymd( min_date ))

# load algn md 

#verify df character values in correct format
amd1$sequence_name <- as.character(amd1$sequence_name)
amd1$region <- as.character(amd1$region)

#filter NA in required columns
amd1 <- amd1[ !is.na( amd1$sequence_name ) , ]
amd1 <- amd1[ !is.na( amd1$sample_date ) , ]
amd1 <- amd1[ !is.na( amd1$region ) , ]

#verify df date in correct format 
amd1$sample_date <- as.Date( as.character( amd1$sample_date ))
amd1 <- amd1[amd1$sequence_name %in% tre$tip.label,]
stopifnot( all(amd$sequence_name %in% tre$tip.label) )

if ( !('sample_time' %in% colnames(amd1))){
  amd1$sample_time = decimal_date (amd1$sample_date)
}

amd1$sts <- amd1$sample_time 

if (!('lineage' %in% colnames( amd1 ))){
  amd1$lineage <- 'lineage_not_provided'
} else { 
  
  amd1$lineage <- as.character(amd1$lineage)
  
}

if (!('mutations' %in% colnames( amd1 ))){
  amd1$mutations <- 'mutations_not_provided'
} else { 
  
  amd1$mutations = as.character(amd1$mutations)
  
}

# filter by sample time 
amd1 <- amd1 [ (amd1$sample_time >= min_time) & (amd1$sample_time <= max_time) , ] 
sts <- setNames( amd1$sample_time, amd1$sequence_name )

# retain only required variables 
keep_col <- unique( c('sequence_name', 'sample_time', 'sample_date', 'region', 'lineage', 'mutations', test_cluster_odds) )
amd1 <- amd1[ , ..keep_col] 

# prune tree
if ( !is.rooted( tre ) ){
  if ( !( root_on_tip %in% amd1$sequence_name)){
    stopifnot( root_on_tip %in% tre$tip.label )
    
    root_tre <- data.frame(sequence_name = root_on_tip
                           , sample_time = root_on_tip_sample_time
                           , sample_date = as.Date( date_decimal( root_on_tip_sample_time ) ))
    
    root_tre[setdiff(names(amd1), names(root_tre))] <- NA
    amd1 <- rbind(amd1, root_tre)
    rm(root_tre)
    
  }
}




