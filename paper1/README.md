# The resurgence of dengue in Brazil, 2018-2019
This repository contains the data, results, and code for the manuscript of Brito et al 2020.

## data

The directory `data` contains sequences, metadata and epidemiological data used to generate the figures shown in the paper. In that folder, one can find:

```
data/
│
├── denv1/
│   ├── genomics		→ DENV-1 genome and envelope sequences, and their metadata
│   └── phylogenetics 		→ DENV-1 maximum likelihood and Bayesian phylogenies
│
├── denv2/
│   ├── genomics		→ DENV-2 genome and envelope sequences, and their metadata
│   └── phylogenetics 		→ DENV-2 maximum likelihood and Bayesian phylogenies
│
├── epidemiology 		→ Dengue, Zika, and Chikungunya case data, per country, region and state
│
└── mosquito-suitability	→ Weekly estimates of mosquito suitability

code/
```

## Nextstrain

The directory `auspice`, in the main folder of 'DENV-genomics' contains the json files produced using `augur` from the nextstrain pipeline. These results can be visualized using `auspice`, accessing the links below:

* [Spatiotemporal dispersal of DENV-1 in the Americas (using Envelope sequences)](https://nextstrain.org/community/grubaughlab/DENV-genomics/DENV1-Brazil)
* [Spatiotemporal dispersal of DENV-2 in the Americas (using Envelope sequences)](https://nextstrain.org/community/grubaughlab/DENV-genomics/DENV2-Brazil)

---

**Grubaugh Lab** | Yale School of Public Health (YSPH) | [https://grubaughlab.com/](https://grubaughlab.com/)
