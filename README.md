# HyphalNet

This project represents the implementation of community detection algorithms in biological networks. Our basic approach is shown below. It requires using an underlying protein-rotein interaciton network to derive sample-specific networks, that can then be compressed to a multigraph and reduced to functionally-specific communities.

![Community Detection](/img/community_detection.jpg)

We are currently evaluating this approach in three separate domains.

## Getting Started
To run hyphalNet you need to download an interactome.

1 - Download the latest [STRING DB](https://stringdb-static.org/download/protein.links.detailed.v11.0/9606.protein.links.detailed.v11.0.txt.gz) interactome
2 - Download the [mapping file](https://stringdb-static.org/download/protein.info.v11.0/9606.protein.info.v11.0.txt.gz)
3 - run:
```
python ./bin/graphTools.py --graphFile [stringfile] --graphSource=string --nodeMapping=[mappingfile] --dest=./data
```

This will load up the interactomes in the correct format (for now we are working across 3 formats)


## Cancer Proteomics
We hypothesize that we can identify underlying patterns across patient cohorts in patient datasets from the [Proteomic Data Commons](https://pdc.cancer.gov/).

This approach can be described in more detail [here](examples/cancerProtFeatures) along with test cases that work on local files.

## Mutational networks
We are also evaluating the ability to identify mutational networks that are shared across patients or clinical samples. This approach will be described [here](examples/mutationDrugResponse).

## Proteomics drug response
Building networks out of proteomics data, we also hope to identify specific communties that give rise to changes in drug response.
