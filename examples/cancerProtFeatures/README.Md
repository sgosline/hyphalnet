## Pan-Cancer Feature Extraction

For this example we focus on using protein abundances to build patient-specific networks that can then be compared across tumor type to derive specific communities

### Getting started
Currently we have built a suite of tools to evaluate the feasibility of this approach. To do this we employ the `hyphalnet` package to build networks of patients and communities across those patients.

``` python
python3 panCancerTest.py --enrich --getDistances

```

This will run the network analysis on 4 test datasets from the Proteomic data commons and compute the Gene Ontology enrichment and also output the data for the next analysis.

``` R
Rscript panCancerEnrichment.R
Rscript netStats.R
Rscript communityStats.R
```
This will compute the plotting and analysis of the ability of Forests to properly summarize the data, and the ability of the communities to represent some sort of summary.

### Assessing community detection
While we have the framework to detect communities we still have to ascertain how to tune this properly. This really relies on the algorithmic develeopment to improve this.

### Open Issues
1- Currently I'm downloading the files manually - how do we automate this?
2- Heavy reliance on `OmicsIntegrator` makes it difficult to expand
3-
