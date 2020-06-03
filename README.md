# HyphalNet

This project represents the implementation of community detection algorithms in biological networks. Our basic approach is shown below. It requires using an underlying protein-rotein interaciton network to derive sample-specific networks, that can then be compressed to a multigraph and reduced to functionally-specific communities.

We are currently evaluating this approach in three separate domains.

## Cancer Proteomics
We hypothesize that we can identify underlying patterns across patient cohorts in patient datasets from the [Proteomic Data Commons](https://pdc.cancer.gov/).

This approach can be described in more detail [here](examples/cancerProtFeatures).

## Mutational networks
We are also evaluating the ability to identify mutational networks that are shared across patients or clinical samples. This approach will be described [here](examples/mutationDrugResponse).

## Proteomics drug response
Building networks out of proteomics data, we also hope to identify specific communties that give rise to changes in drug response.
