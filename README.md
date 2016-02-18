# LGR_TotalEscapement
Estimate total adult Chinook and steelhead escapement past Lower Granite Dam.

## File descriptions
- DataInstructions.txt: Where to go to download data from DART
- DataPrep.R: pull together window counts, trap rates, data from trap, etc. and summarize on weekly basis
- Escapement_AllYrs.R: Run JAGS model of total escapement, save results
- Summarize_Results: Pull together posteriors across all years and species. Contains scripts for MCMC diagnostic plots, as well as plots of results.

- ModelFiles: Contains JAGS model files
- ModelFits: Contains results (e.g. JAGS model objects), one for each species/year combination
