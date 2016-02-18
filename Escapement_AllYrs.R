# Author: Kevin See
# Purpose: this R script estimates adult Chinook and steelhead LGR passage by week based on both window counts and counts in a fish trap in the ladder
# Created: 2/19/2015
# Last Modified: 2/17/2016
# Notes: It is run on a weekly time-step

library(R2jags)

#-----------------------------------------------------------------
# load the data for all years
load(file = 'DataPrepped.rda')
# load the data list prepped for JAGS
load(file = 'JAGS_DataList.rda')

# create folder to store results
if(!file.exists('ModelFits')) dir.create(file.path(getwd(), 'ModelFits'))

#-----------------------------------------------------------------
# set mcmc parameters
#-----------------------------------------------------------------
# number of total samples
mcmc.chainLength = as.integer(40000) 

# number of burn-in samples
mcmc.burn = as.integer (10000)

# thinning interval
mcmc.thin = 30

# number of MCMC chains
mcmc.chains = 4

# how many poseterior samples will we have?
data.frame(Each.Chain = (mcmc.chainLength - mcmc.burn) / mcmc.thin, All.Chains = (mcmc.chainLength - mcmc.burn) / mcmc.thin * mcmc.chains)


#-----------------------------------------------------------------
# run model for several years
#-----------------------------------------------------------------
model.loc = 'ModelFiles/LGR_TotalEscape_JAGS.txt'

jags.params = c('X.log.all', 'X.all', 'X.day', 'X.night', 'X.reasc', 'X.all.wild.morph', 'X.new.wild.morph', 'X.all.wild.pbt', 'X.new.wild.pbt', 'X.reasc.wild.morph', 'X.reasc.wild.pbt', 'X.night.wild.morph', 'X.night.wild.pbt', 'X.tot.all', 'X.tot.day', 'X.tot.night', 'X.tot.reasc', 'X.tot.all.wild.morph','X.tot.new.wild.morph', 'X.tot.night.wild.morph', 'X.tot.reasc.wild.morph', 'X.tot.all.wild.pbt', 'X.tot.new.wild.pbt', 'X.tot.night.wild.pbt', 'X.tot.reasc.wild.pbt', 'prop.tagged.morph', 'prop.tagged.pbt', 'X.sigma', 'true.prop', 'win.prop.avg', 'win.prop.true', 'win.prop.sigma', 'hist.prop', 'reasc.avg', 'reasc.true', 'reasc.sigma', 'acf', 'wnc.avg', 'wnc.true', 'wnc.sigma', 'trap.rate.true', 'r')

# if using fixed, known trap rate
# jags.params = c(jags.params[-match('trap.rate.true', jags.params)], 'trap.bump')


# set initial values
jags.inits = function(){
  return(list('X.log.all' = with(jags.data, log(Y.window + 1)),
              'trap.rate.true' = with(jags.data, trap.rate + 0.0001)))
}

set.seed(17)
for(i in 1:length(jags_data_list)) {
  cat(paste('Starting', names(jags_data_list)[i], 'run \n'))
  
  # pull out data for JAGS, and data for plotting later
  jags.data = jags_data_list[[i]]

  ptm = proc.time()
  adult.pass.mod = try(jags(data = jags.data, 
                            inits = jags.inits, 
                            parameters.to.save = jags.params, 
                            model.file = model.loc, 
                            n.chains = mcmc.chains, 
                            n.burnin = mcmc.burn, 
                            n.thin = mcmc.thin, 
                            n.iter = mcmc.chainLength, 
                            DIC = FALSE))
  proc.time() - ptm # returns the CPU time used
  if(class(adult.pass.mod)=='try-error') {
    cat(paste('Error with', names(jags_data_list)[i], 'run. \n'))
    rm(adult.pass.mod, jags.data, ptm)
    next
  }
  
  # save results
  save(adult.pass.mod, jags.data, file = paste0('ModelFits/LGD_', gsub('\\.', '_', names(jags_data_list)[i]), '.rda'))

  # remove some objects before running next species/year
  rm(adult.pass.mod, jags.data, ptm)
}



