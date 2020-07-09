## Pan-Cancer Feature Extraction

For this example we focus on using protein abundances to build patient-specific networks that can then be compared across tumor type to derive specific communities

### Getting started
Currently we have built a suite of tools to evaluate the feasibility of this approach. To do this we employ the `hyphalnet` package to build networks of patients and communities across those patients. You will first need to downlooad a network file and build the graphs.

#### Files needed

Make sure you have downloaded and built the files according to the [main README file](../../README.Md).

_GO enrichment files_
You will also need to create a `db` directory and download these files:
1- [GO OBO file](http://purl.obolibrary.org/obo/go/go-basic.obo)
2- [Gene GO mapping](ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz)


#### Creating HyphalNetworks

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

I'm currently working on the following:

1 - Currently I'm downloading the files manually - how do we automate this?
2 - How do we assess which communities are meaningful?
