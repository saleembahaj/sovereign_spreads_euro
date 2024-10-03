This folder contains:
1) the narrative of events (event spreadsheet.xlsx) 
2) the daily and monthly event reaction series ("country"_instrument.xlsx)
3) the data used in the VAR (VAR_Database.xlsx)
3) Replication code for the VAR estimates (Figures 4-7 in the paper) (run runner.m in matlab).

Interested users may also want to review the matlab function for the Gibbs sampler: PBVARX_HIERARCHICAL_FUN_COMP.m
The sub functions in runner.m : data_read.m and specification.m build the sample and set up the sampler.
None of the other matlab functions are of particular interest.

Note that the specification is set to run a shorter MCMC chain than in the main specification in the text. This is to save on computational time for those who just want to hit run.

Saleem Bahaj
28th December 2018
Contact: saleembahaj@googlemail.com


   