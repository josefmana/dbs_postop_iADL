# All analyses reported in the article (Beydicek et al., under review) ran in R version 4.0.5 (2021-03-31),
# on x86_64-w64-mingw32/x64 (64-bit) platform under Windows 10 x64 (build 19043).

# I used the following versions of packages employed: dplyr_1.0.7, tidyverse_1.3.0, DiagrammeR_1.0.9, DiagrammeRsvg_0.1,
# rsvg_2.3.1, brms_2.16.3, loo_2.4.1, tidybayes_2.3.1, ggplot2_3.3.3 and  patchwork_1.1.1

# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list required packages
pkgs <- c(
  "dplyr", # for objects manipulation
  "tidyverse", # for pivot_longer
  "DiagrammeR", # for DAG creation
  "DiagrammeRsvg", # for saving plots crated in DiagrammeR
  "rsvg", # for saving plots crated in DiagrammeR
  "brms", # for Bayesian model fitting / interface with Stan
  "tidybayes", # for posteriors manipulations
  "ggplot2", # for plotting
  "patchwork" # for ggplots manipulations
)

# load required packages
# prints NULL if a package is already installed
sapply(
  pkgs, # packages to be loaded/installed
  function(i)
    if ( !require( i , character.only = T ) ){
      # it's important to have the 'character.only = T' command here
      install.packages(i)
      library(i)
    }
)

# set ggplot theme
theme_set( theme_classic(base_size = 25) )

# prepare colors to use in graphs (a colorblind-friendly palette)
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

# create folders "models", "figures" and "tables" to store results in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply(
  c("models", "figures", "tables"), # folders to be created
  function(i)
    if( !dir.exists(i) ) dir.create(i)
)

# Note that although I set a seed for all models, the results are only exactly
# reproducible on the same operating system with the same C++ compiler and version.

# read data sets
d.par <- read.table( "20220224_stim_pars.txt" , header = T ) # stimulation parameters
d0 <- read.csv( "20220207_dbs_iadl_data_full.csv" , sep = "," ) %>%
  # keeping only those variables I will need (PDAQ, ID, DRS-2, BDI-II, LEDD and demographics)
  select(
    id, assessment_type, levodopa_equivalent, drsii_total, bdi, # predictors
    time_after_surg, age_neuropsy, sex, edu_type, years_edu, # demographics
    starts_with("pdaq") # outcome
  ) %>%
  # rename the variables so that I don't go crazy
  rename(
    ass_type = assessment_type , led = levodopa_equivalent , drs = drsii_total
  ) %>%
  # change "post_1y" to "post" because they are all one-year after DBS operation
  mutate(
    ass_type = as.factor(
      ifelse( ass_type == "post_1y" , "1_post" , "0_pre" )
    )
  )

# read item descriptions for plotting
it_desc <- read.csv( "item_nms.csv" , sep = "," , row.names = 1 )


# ----------- heuristic causal model (DAG)  -----------

# will be created via DiagrammeR
# for now, I know I need to adjust for preop. DRS2, BDI-II and LEDD for total effect
# and adjust for postop. DRS-2, BDI-II and LEDD for direct effect


# ----------- data pre-processing  -----------

# check there ain't no NAs in the data set
sapply( names(d0) , function(i) sum( is.na(d0[[i]]) ) ) # seems ok
all( sapply( names(d0) , function(i) sum( is.na(d0[[i]]) ) ) == 0 ) # indeed, it is ok, no NAs

# extract the number of patients
N = length( unique(d0$id) ) # N = 32

# find a reasonable scaling of predictors
# first get an idea about in-sample variability
sapply(
  levels(d0$ass_type) , function(i) # loop through pre/post assessments
    sapply(
      c("drs","bdi","led") , function(j) # loop through all predictors
        round( sd( d0[[j]][ d0$ass_type == i ] , na.rm = T ) , 2 )
    )
)

# will center all predictors of interest on their pre-surgery mean
# extract the means from data and save them to a table
scl <- data.frame(
  # first a column of means to use for scaling
  M = c(
    mean( d0$drs[d0$ass_type == "0_pre"] , na.rm = T ),
    mean( d0$bdi[d0$ass_type == "0_pre"] , na.rm = T ),
    mean( d0$led[d0$ass_type == "0_pre"] , na.rm = T )
  ),
  # then add user(me)-specified scaling factors based on theoretical consideration and data variability
  SD = c(
    3, # 3 point for DRS-2, half-way between max (144/144) and likely MCI (138/144), close to SDs pre and post
    7, # 7 points for BDI/II, equivalent to feeling better/worse in a third of the items (7/21), close to SDs pre and post
    500 # 500 mg for LEDD, equivalent to 2 Isicom 250 or 5 Isicom 100 pills, close to SD post
  ),
  # add names to rows so that I know which row belongs to which variable
  row.names = c( "drs" , "bdi" , "led" )
)

# transform the data from wide to long format over PDAQ items
d1 <- d0 %>% pivot_longer(
  cols = paste0( "pdaq_" , 1:15 ) , names_to = "item" , values_to = "resp"
) %>% mutate(
  # change names of PDAQ items to integers
  # the gsub part removes everything before the underscore (text)
  # so that only integer remains and is easy to transform via as.integer()
  # could have also gone via as.integer( as.factor() )
  item = as.integer( gsub( "^.*\\_" , "" , item ) )
)

# calculate number of PDAQ items (K = 15)
K <- length( unique(d1$item) )
K == max( d1$item ) # A-OK

# check that patient/item/assessment type each combination is there exactly once
table( d1$id , d1$item , d1$ass_type ) # a-ok

# check the contrast for pre/post factor
contrasts( d1$ass_type ) # ok, it's a dummy coding as I wanted
#contrasts( d1$ass_type ) <- -contr.sum(2)/2 # pre = -0.5, post = 0.5

# scale DRS and LEDD such that the model samples more efficiently
# make the outcome ordinal such that brms can take it in
d2 <- d1 %>%
  mutate(
    drs = ( drs - scl["drs","M"] ) / scl["drs","SD"],
    bdi = ( drs - scl["bdi","M"] ) / scl["bdi","SD"],
    led = ( led - scl["led","M"] ) / scl["led","SD"],
    resp = as.ordered( resp )
  )

# plot raw data to have something nice to look at while Stan will be fitting models
table( d2[ , c("item","resp" , "ass_type") ] ) %>%
  as.data.frame  %>%
  mutate(
    item = factor( rep( it_desc$name , 10 ), levels = it_desc$name , ordered = T ),
    Time = factor(
      ifelse( ass_type == "0_pre" , "pre-surgery", "post-surgery"),
      levels = paste0(c("pre","post"),"-surgery"), ordered = T
    ),
    Freq = Freq + 1 # shift such that zero-response frequencies are still visualized
  ) %>%
  # plotting proper
  ggplot( aes(x = resp , y = Freq , fill = Time ) ) +
  geom_bar( stat = "identity" , position = position_dodge( width = .8 ) , width = .7 ) +
  facet_wrap( ~ item  , scales = "free" , nrow = 5 , ncol = 3 ) +
  labs( x = "Response" , y = "Frequency" ) +
  scale_fill_manual( values = cbPal[c(1,6)] ) +
  scale_y_continuous( limits = c(0, 32), breaks = seq(1,31,5) , labels = seq(0,30,5) ) +
  theme_minimal( base_size = 25 ) +
  theme( legend.position = "none" , panel.grid.minor = NULL )

# save the plot as Fis S1
ggsave( "figures/Fig S1 Item-level responses.png" , width = 1.5*10.9, height = 2*11.8, dpi = "retina" )
