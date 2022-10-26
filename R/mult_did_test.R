##
## Dynamic Treatment Effects via Difference-in-Differences
##  - multiple period DiD
##  - staggered DiD
##  - - vs. static DiD
##
## Problem that BSTS can solve: 
##  How to bin event study design
##


library(did)

sp <- reset.sim()

set.seed(1814)

time.periods <- 4

## cannot reproduce results at 
## https://cran.r-project.org/web/packages/did/vignettes/did-basics.html
## sp$te.e <- 1:time.periods
sp$te.e <- (1:time.periods) - 1
 
dta <- build_sim_dataset(sp)

nrow(dta)

head(dta)

example_attgt <- att_gt(yname = "Y",
                        tname = "period",
                        idname = "id",
                        gname = "G",
                        xformla = ~X,
                        data = dta
)

summary(example_attgt)

ggdid(example_attgt)


##------------------------
## AGGREGATION
##------------------------


agg.simple <- aggte(example_attgt, type='simple')
summary(agg.simple)

## DYNAMIC EFFECTS AND EVENT STUDIES
agg.es <- aggte(example_attgt, type = "dynamic")
summary(agg.es)
ggdid(agg.es)

## Group-Specific Effects
agg.gs <- aggte(example_attgt, type = "group")
summary(agg.gs)
ggdid(agg.gs)


## Calendar Time Effects
agg.ct <- aggte(example_attgt, type = "calendar")
summary(agg.ct)
ggdid(agg.ct)



##-------------------------
## Alternative Control Groups
##-------------------------
## Not yet treated
example_attgt_altcontrol <- att_gt(yname = "Y",
                                   tname = "period",
                                   idname = "id",
                                   gname = "G",
                                   xformla = ~X,
                                   data = dta,
                                   control_group = "notyettreated"          
)
summary(example_attgt_altcontrol)




##------------------------------
## EXAMPLE DATA 
##  - Minimum Wage effect on county-level teen employment rates
##------------------------------
data(mpdta)
head(mpdta)

# estimate group-time average treatment effects without covariates
mw.attgt <- att_gt(yname = "lemp",
                   gname = "first.treat",
                   idname = "countyreal",
                   tname = "year",
                   xformla = ~1,
                   data = mpdta,
)

# summarize the results
summary(mw.attgt)

ggdid(mw.attgt, ylim = c(-.3,.3))

# aggregate the group-time average treatment effects
mw.dyn <- aggte(mw.attgt, type = "dynamic")
summary(mw.dyn)
ggdid(mw.dyn, ylim = c(-.3,.3))


## Balance across group exposure levels
## (avoid confounding dynamics and selective treatment timing among different groups)
mw.dyn.balance <- aggte(mw.attgt, type = "dynamic", balance_e=1)
summary(mw.dyn.balance)














