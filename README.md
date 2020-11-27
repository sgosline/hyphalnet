# HyphalNet

This project represents the implementation of community detection algorithms in biological networks. Our basic approach is shown below. It requires using an underlying protein-rotein interaciton network to derive sample-specific networks, that can then be compressed to a multigraph and reduced to functionally-specific communities.

![Community Detection](/img/community_detection.jpg)

We are currently evaluating this approach in three separate domains.

## Install
To install the package, clone the repository.
1. ``` git clone https://github.com/sgosline/hyphalnet.git```
2. Then run: ```sudo python setup.py install```

This will install the package and its files.

## Getting Started
To run hyphalNet you need to download an interactome and create an `igraph` Graph and compress it using `pickle`. This script will do it for you.
1. Download the latest [STRING DB](https://stringdb-static.org/download/protein.links.detailed.v11.0/9606.protein.links.detailed.v11.0.txt.gz) interactome
2. Download the [mapping file](https://stringdb-static.org/download/protein.info.v11.0/9606.protein.info.v11.0.txt.gz)
3. Run:
```
python ./bin/graphTools.py --graphFile [stringfile] --graphSource=string --nodeMapping=[mappingfile] --dest=./data
```

This will load up the interactomes in the correct format (for now we are working across 3 formats)


## Cancer Signatures
We have built a script that downloads data from the [CPTAC]() repositories using the Python package `cptac`. To construct cancer signatures from the available data.

To build these scripts, run the script

``` python
python bin/buildCancerSigs.py --help
```
This will show the available options to build new cancer signatures.

The output of this will be a hyphalNetwork object in `pkl` format that can be used to map to novel datasets.

## Signature statistics
To collect information about the signatures, particularly ahead of paper publication, we run the following script to collect gene enrichment statistics, and distance metrics within the cancer signature data itself. Some of these are more useful than others.

### Network Mutual Information
We wanted to show that the signatures we were getting were unique to the data collected. To do that we evaluated the newtork mutual information statistic. This probably does not need to be run more generally.

``` python
python bin/getSigNMI.py

```

### Signature gene enrichment
We also wanted to capture tio biological activity of the underlying cancer signatures. This statistic _is_ of interest for others wanting to do more analysis downstream.

``` python
python bin/hyphEnrich.py
```

## Applying signature to new data
Lastly we need to apply the signatures to novel datasets. This script expects a csv of proteomics measurements in tidied format.

``` python
python bin/mapSigsToData.py
```
