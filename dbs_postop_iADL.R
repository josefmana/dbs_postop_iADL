# All analyses reported in the article (Beydicek et al., under review) ran in R version 4.0.5 (2021-03-31),
# on x86_64-w64-mingw32/x64 (64-bit) platform under Windows 10 x64 (build 19043).

# I used the following versions of packages employed: dplyr_1.0.7, tidyverse_1.3.0, DiagrammeR_1.0.9,
# DiagrammeRsvg_0.1, rsvg_2.3.1, brms_2.16.3, loo_2.4.1, tidybayes_2.3.1, ggplot2_3.3.3 and  patchwork_1.1.1

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
  "bayestestR", # for calculation of the probability of direction (pd)
  "ggplot2", # for plotting
  "ggpubr" # for ggplots manipulations
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
theme_set( bayesplot::theme_default( base_size = 20 ) )

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
  # first create a column of means to use for scaling
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

# scale DRS and LEDD such that the model samples more efficiently
# make the outcome ordinal such that brms can take it in
d2 <- d1 %>%
  mutate(
    drs = ( drs - scl["drs","M"] ) / scl["drs","SD"],
    bdi = ( bdi - scl["bdi","M"] ) / scl["bdi","SD"],
    led = ( led - scl["led","M"] ) / scl["led","SD"],
    resp = as.ordered( resp )
  )

# check the contrast for pre/post factor
contrasts( d2$ass_type ) # ok, it's a dummy coding as I wanted
#contrasts( d2$ass_type ) <- -contr.sum(2)/2 # pre = -0.5, post = 0.5

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

# save the plot as Fig S1
ggsave( "figures/Fig S1 Item-level responses.jpg", width = 1.75*10.9, height = 2*11.8, dpi = "retina" )


# ----------- model fitting (direct effect) -----------

# set rstan options
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores
s = 87542 # seed for reproducibility
ch = 4 # number of chains
it = 1500 # iterations per chain
wu = 500 # warm-up iterations, to be discarded
ad = .90 # adapt_delta parameter

# set up some reasonable (regularizing) priors
p <- c(
  prior( normal( 0 , 0.5 ) , class = b ),
  prior( student_t( 3 , 0 , 2.5 ) , class = Intercept ),
  prior( student_t( 3 , 0 , 2.5 ) , class = sd , group = id ),
  prior( student_t( 3 , 0 , 2.5 ) , class = sd , group = item )
)

# fit an ordinal-logit GLMM via brms
m <- brm(
  formula = resp ~ ass_type * led + ass_type * drs + ass_type * bdi + ( 1 | id ) + ( 1 | item ),
  prior = p, family = cumulative( link = "logit" , link_disc = "log" , threshold = "flexible" ), data = d2,
  seed = s, chains = ch, iter = it, warmup = wu, control = list( adapt_delta = ad ),
  file = "models/direct_effect.rds"
)

# check some model diagnostics
max( rhat(m)[ !is.nan(rhat(m)) ] ) # 1.005
max( loo(m)$diagnostics$pareto_k ) # 0.361
plot( loo(m) ) # good


# ----------- model post-processing -----------

# extract names of parameters of m
pars <- get_variables(m)

# create a Table 2 from posterior draws of model m
t2 <- m %>%
  spread_draws( `b_.*` , regex = T ) %>% # extracting only fixed-effects for Table 2
  median_hdi( .width = .95) %>% # extract 95% posterior HDIs
  pivot_longer( cols = pars[1:11], values_to = "b" , names_to = "Parameter" ) # pivot such that rows are predictors 

# create a 95% HDI column in t2
for ( i in 1:nrow(t2) )t2[ i , "95% HDI"] <- paste0(
  "[", sprintf( "%.2f" , round( t2[ i , paste0(t2$Parameter[i], ".lower") ] , 2 ) ), ", ", # HDI lowermost point
  sprintf( "%.2f" , round( t2[ i , paste0(t2$Parameter[i], ".upper") ] , 2 ) ), "]" # HDI uppermost point
)

# add p-direction via bayestestR
add_pd <- function( m = m ) {
  dum_pd <- m %>% bayestestR::p_direction()
  return( sprintf( "%.3f" , round( dum_pd$pd , 3 ) ) )
}

# finish the table by getting rid of dummy columns
t2 <- t2 %>%
  select( -ends_with("upper") , -ends_with("lower") , -.width , -.point , -.interval ) %>%
  # and tidy columns "b" and "Parameter"
  mutate(
    b = sprintf( "%.2f" , round( b , 2 ) ),
    Parameter = recode(Parameter, 
                       "b_Intercept[1]" = "threshold[0-1]",
                       "b_Intercept[2]" = "threshold[1-2]",
                       "b_Intercept[3]" = "threshold[2-3]",
                       "b_Intercept[4]" = "threshold[3-4]",
                       "b_ass_type1_post" = "Time of assessment",
                       "b_led" = "LEDD",
                       "b_drs" = "DRS-2",
                       "b_bdi" = "BDI-II",
                       "b_ass_type1_post:led" = "Time of assessment × LEDD",
                       "b_ass_type1_post:drs" = "Time of assessment × DRS-2",
                       "b_ass_type1_post:bdi" = "Time of assessment × BDI-II"
    ),
    `Pr(b > 0)` = add_pd(m)
  )

# calculate contrast between pre- vs. post-surgery PDAQ while keeping remaining predictors unaltered
# i.e, looking for i in {0, 1, 2, 3, 4} to calculate from posteiro the contrast:
# Pr(answer = i | ass_type = post, drs = 0, led = 0, bdi = 0) -
# Pr(answer = i | ass_type = pre, drs = 0, led = 0, bdi = 0)
ppred <- posterior_epred(
  # the model and  predictors' levels of interest
  m, newdata = data.frame(
    ass_type = c( "0_pre", "1_post" ), led = 0, drs = 0, bdi = 0, id = NA, item = NA
  ),
  # exclude group-level variance, i.e., we will get only expectations base of the "fixed-effect" part of the model
  re_formula = NA
  # this gives me a 4000 (iterations) x 6 (predictor levels) x 5 (possible answers) array
  # the indexing works like this: ppred[ iteration , assessment (pre vs post) , answer ]
)

# create function for extracting median and hdi
extract_md_hdi <- function( ppred = ppred , row ) {
  # get medians for the row in question
  md <- apply( ppred[ , row , ] , 2 , median )
  # compute 95 % HDI using bayestestR hdi() function
  hdi <- apply( ppred[ , row , ] , 2 , bayestestR::hdi , ci = .95 )
  hdi <- do.call( rbind.data.frame , hdi ) # will need this for the next step to work
  # create the outcome
  # showing the results as percentages
  outcome <- paste0(
    sprintf( "%.1f" , round( 100*md , 1 ) ) , " [",
    sprintf( "%.1f" , round( 100*hdi$CI_low , 1 ) ), ", ",
    sprintf( "%.1f" , round( 100*hdi$CI_high , 1 ) ), "]"
  )
  # return the outcome
  return(outcome)
}

# create a function for extracting post-minus-pre contrast for each level of LED
extract_contrast <- function( ppred = ppred ) {
  # subtract predicted post-minus-pre-surgery response probabilities
  contr <- ppred[ , 2 , ] - ppred[ , 1 , ]
  # compute median expectation of the contrast where pre = pre and post = pre+3
  contr.md <- apply( contr , 2 , median )
  # compute 95 % hdi using bayestestR hdi() function (because tidybayes' hdi didn't work for some contrasts)
  contr.hdi <- apply( contr , 2 , bayestestR::hdi , ci = .95 )
  contr.hdi <- do.call( rbind.data.frame , contr.hdi ) # will need this for the next step to work
  # show the results as percentages
  outcome <- paste0(
    sprintf( "%.1f" , round( 100*contr.md , 1 ) ) , " [",
    sprintf( "%.1f" , round( 100*contr.hdi$CI_low , 1 ) ), ", ",
    sprintf( "%.1f" , round( 100*contr.hdi$CI_high , 1 ) ), "]"
  )
  # return the outcome
  return(outcome)
}

# Prepare supplementary table 1 containing expected response probabilities stratified by assessment (pre-, post-),
# as well as post-minus-pre contrast, i.e., it's a table of posterior predictive distribution across assessments
# marginalized over all other predictors (LED, DRS-2, BDI-II) at their mean in-sample level
# in the article, this table is presented in text
t.s1 <- data.frame(
  `Pre-surgery` = extract_md_hdi( ppred , 1 ),
  `Post-surgery` = extract_md_hdi( ppred , 2 ),
  `Post-minus-pre-surgery` = extract_contrast( ppred ),
  row.names = paste0( "Pr(resp = ", 0:4, ")")
) %>% t

# next prepare posterior predictive distribution for expected response probabilities across assessments (pre-, post-)
# and across eleven distinct LEDD levels
ppred2 <- posterior_epred(
  # the model and  predictors' levels of interest
  m, newdata = data.frame(
    ass_type = c( rep("0_pre",11), rep("1_post",11) ) ,
    led = ( rep( seq(0 , 5000 ,500 )  , 2 ) - scl["led","M"] ) / scl["led","SD"],
    # now the remaining predictors to marginalize over
    drs = 0, bdi = 0, id = NA , item = NA
  ),
  # exclude group-level variance, i.e., we will get only expectations base of the "fixed-effect" part of the model
  re_formula = NA
)

# for better readability, present in mean ± SD format, create a function to do this
extract_m_sd <- function( ppred = ppred , row ) {
  # get means and sds for the row in question
  m <- apply( ppred[ , row , ] , 2 , mean )
  sd <- apply( ppred[ , row , ] , 2 , sd )
  # create the outcome
  # showing the results as percentages
  outcome <- paste0(
    sprintf( "%.1f" , round( 100*m , 1 ) ) , " ± ",
    sprintf( "%.1f" , round( 100*sd , 1 ) ) , " %"
  )
  # return the outcome
  return(outcome)
}

# prepare a data.frame for Table 3 to be filled-in by posterior predictions across different levels of LEDD
t3 <- data.frame(
  Assessment = c( rep("Pre-surgery",11), rep("Post-surgery",11) ), # technically splitting to two tables
  # now all different LEDD values I'd like to see
  LEDD = sprintf(
    "%.0f" , round(
      rep( seq(0 , 5000 ,500 )  , 2 ), 0
    )
  ), ans_0 = NA , ans_1 = NA , ans_2 = NA , ans_3 = NA , ans_4 = NA # responses
) %>% rename(
  # change names of Pr(outcome) to a better one
  "Pr(resp = 0)" = "ans_0" , "Pr(resp = 1)" = "ans_1" , "Pr(resp = 2)" = "ans_2" , "Pr(resp = 3)" = "ans_3" , "Pr(resp = 4)" = "ans_4"
)

# fill-in Table 3
for ( i in 1:nrow(t3) ) t3[ i , 3:ncol(t3) ] <- extract_m_sd( ppred2 , i )

# prepare list to hold figures of conditional (marginalized over in-sample means) main effects
f <- list()

# loop through main effects to extract conditional effects via build-in brms plotting function
for ( i in c("ass_type","led","drs","bdi") ) f[[i]] <- conditional_effects( m , categorical = T , effects = i )

# prepare tables for plotting
for ( i in names(f) ) f[[i]][[paste0( i, ":cats__" )]] <- f[[i]][[paste0( i, ":cats__" )]] %>%
  mutate(
    `Response:` = recode(
      cats__ , "4" = '"none"' , "3" = '"a little"' , "2" = '"somewhat"' , "1" = '"a lot"' , "0" = '"cannot do"'
    )
  )

# tidy the "ass_type" part of this figure
f$A <- f$ass_type[["ass_type:cats__"]] %>%
  ggplot( aes(x = ass_type , y = estimate__ , ymin = lower__ , ymax = upper__ , color = `Response:`) ) +
  geom_point( position = position_dodge(.66) , size = 7 ) +
  geom_errorbar( position = position_dodge(.66) , size = 1 , width = .5) +
  scale_x_discrete(
    name = "Time of assessment" , labels = c("Pre-surgery","Post-surgery")
  ) +
  scale_y_continuous(
    name = "Probability of response (%)" , limits = c(0,1) , breaks = seq(0,1,.1) , labels = seq(0,100,10)
  )

# continue with the LEDD figue
f$B <- f$led[["led:cats__"]] %>%
  ggplot( aes(x = led , y = estimate__ , ymin = lower__ , ymax = upper__ , color = `Response:` , fill = `Response:`) ) +
  geom_line( size = 2 ) +
  geom_ribbon( alpha = .1 , linetype = 0 ) +
  scale_x_continuous(
    name = "LEDD (mg)", labels = seq(500,4000,500),
    breaks = (seq(500,4000,500) - scl["led","M"] ) / scl["led","SD"]
  ) +
  scale_y_continuous(
    name = "Probability of response (%)" , limits = c(0,1) , breaks = seq(0,1,.1) , labels = seq(0,100,10)
  )

# do the same for the DRS-2 figure
f$C <- f$drs[["drs:cats__"]] %>%
  ggplot( aes(x = drs , y = estimate__ , ymin = lower__ , ymax = upper__ , color = `Response:` , fill = `Response:`) ) +
  geom_line( size = 2 ) +
  geom_ribbon( alpha = .1 , linetype = 0 ) +
  scale_x_continuous(
    name = "DRS-2", labels = seq(132,144,3),
    breaks = (seq(132,144,3) - scl["drs","M"] ) / scl["drs","SD"]
  ) +
  scale_y_continuous(
    name = "Probability of response (%)" , limits = c(0,1) , breaks = seq(0,1,.1) , labels = seq(0,100,10)
  )

# finish with the BDI-II figure
f$D <- f$bdi[["bdi:cats__"]] %>%
  ggplot( aes(x = bdi , y = estimate__ , ymin = lower__ , ymax = upper__ , color = `Response:` , fill = `Response:`) ) +
  geom_line( size = 2 ) +
  geom_ribbon( alpha = .1 , linetype = 0 ) +
  scale_x_continuous(
    name = "BDI-II", labels = seq(0,27,3),
    breaks = (seq(0,27,3) - scl["bdi","M"] ) / scl["bdi","SD"]
  ) +
  scale_y_continuous(
    name = "Probability of response (%)" , limits = c(0,1) , breaks = seq(0,1,.1) , labels = seq(0,100,10)
  )

# make the legend text bigger
for ( i in c("A","B","C", "D") ) f[[i]] <- f[[i]] + theme( legend.text = element_text( size = 18) )

# put them all together to create Figure 2
f2 <- ggarrange(
  f$A, f$B, f$C, f$D, nrow = 2, ncol = 2,
  labels = c( "A", "B", "C", "D" ), font.label = list( size = 24 ),
  common.legend = T, legend = "bottom"
)

# save Fig.2
ggsave( "figures/Fig 2 Marginal effects.jpg", width = 1.5 * 9.09, height = 2.5 * 4.31, dpi = "retina" )
