# Sovereign Spreads in the Euro Area: Cross Border Transmission and Macroeconomic Implications

# Summary
Large movements in sovereign spreads were at the heart of the euro crisis. This paper builds a new high-frequency narrative dataset of country specific events in the crisis period to identify shocks to sovereign spreads that are orthogonal to the economy. It finds that an increase in sovereign spreads has a contractionary macroeconomic impact with transmission running through a deterioration in private financial conditions. Moreover, the market reactions to foreign events explains a meaningful share of the variation in a sovereign’s cost of borrowing during the crisis.

---

# Authors and Reference:
[Sovereign Spreads in the Euro Area](https://www.sciencedirect.com/science/article/pii/S0304393219300066) Previously circulated under the title “Systemic Sovereign Risk: Macroeconomic Implications in the Euro Area”. 

- [Saleem Bahaj](https://sites.google.com/site/saleembahaj/home)
- Acknowledgements: Aidan Saggers and Sian Besley provided able research assistance

---
# Full Dataset
- See below for the daily and monthly event reaction series  
  [Italy](./Italy_Instrumets.xlsx)  
  [Spain](./Spain_Instrumets.xlsx)  
  [Ireland](./Ireland_Instrumets.xlsx)  
  [Portugal](./Portugal_Instrumets.xlsx)

- See below for the summary dataset for the country specific events  
  [Event](./EventSpreadsheet.xlsx)

---
# Codes
- See below for the matlab function for the Gibbs Sampler  
  [Gibbs Sampler](./PBVARX_HIERARCHICAL_FUN_COMP.m)
- See below for the main function: runner.m  
  [Runner](./Runner.m)
- The sub functions in runner.m: data_read.m and specification.m build the sample and set up the sampler. See below for these functions  
  [data_read](./data_read.m)  
  [specification](./specification.m)

  

# Figures

Below are some key figures visualizing the dataset:

---
- [Fig. 1: cumulative reaction to foreign events](./Fig1.jpg)
- [Fig.7a: historical decomposition](./fig7a.jpg)


