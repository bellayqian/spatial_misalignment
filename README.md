# Spatial Misalignment Correction Methods

This repository contains R scripts for addressing spatial misalignment in health data analysis, with a focus on comparing dasymetric mapping (a naive approach) with advanced Bayesian methods.

## Overview

Spatial misalignment occurs when data are collected at different geographic scales or boundaries than those required for analysis. This is a common challenge in spatial epidemiology and health geography when working with administrative data.

## Methods Implemented

### Dasymetric Mapping

`Dasymetric_mapping.R` implements a traditional approach that redistributes data from one set of geographic boundaries to another using ancillary information (typically population distribution).

The process involves:
1. Starting with source data (e.g., county-level health statistics)
2. Using ancillary data (e.g., census block populations) to determine redistribution weights
3. Creating weights based on the ancillary data
4. Applying these weights to redistribute the data to the target geography

**Limitations:**
- Assumes uniform distribution within the smallest unit (census blocks)
- Edge effects can create artifacts at boundaries
- Quality depends heavily on the ancillary data's accuracy
- Temporal misalignment between different data sources can introduce error

### Atom-based Regression Model (ABRM)

`abrm_sim.R` and `nimble_abrm_0219.R` implement a fully model-based approach to addressing spatial misalignment that does not rely on assumptions of spatial smoothness like other methods.

ABRM uses "atoms" (the areas created by intersecting all sets of units on which variables of interest are measured) as the units of analysis. This approach treats atom-level unmeasured values of covariates or outcomes as latent variables and builds models allowing them to be sampled conditional on observed values.

Benefits of the ABRM approach:
- Provides a fully model-based approach without requiring spatial alignment prior to analysis
- Properly represents uncertainties arising from all analytic procedures
- Enables more robust confounding adjustment through higher resolution analysis
- Provides more realistic uncertainty estimates compared to deterministic reapportionment approaches
- Does not rely on assumptions of smoothness in the underlying spatial risk surface
- May significantly reduce bias in effect estimates (research shows it can reduce bias by up to 85% compared to dasymetric methods)

Note: While ABRM has existed in statistical literature for over two decades, it has rarely been used in practice due to implementation challenges and computational demands. Our implementation in NIMBLE makes these methods more accessible.

## Scripts

- **Dasymetric_mapping.R**: Implementation of the dasymetric mapping approach
- **abrm_sim.R**: Simulation studies for the Aggregate Bayesian Regression Model
- **nimble_abrm_0219.R**: NIMBLE implementation of the Aggregate Bayesian Regression Model
- **comparison_script.R**: Code to compare the performance of different spatial misalignment methods
- **function_0221.R**: Utility functions used by the other scripts

## Usage

Each script is documented with comments explaining the required inputs and expected outputs. Generally, you will need:

1. Source geography data (shapefiles or spatial objects)
2. Target geography data (shapefiles or spatial objects)
3. Health or demographic data associated with the source geography
4. Ancillary data (typically population counts at a fine resolution)

## Requirements

- R (>= 4.0.0)
- Required packages:
  - sf
  - sp
  - nimble
  - raster
  - dplyr
  - ggplot2
  - rgdal (if working with older shapefiles)

## Citation

If you use these methods in your research, please cite this repository and the relevant methodological papers.
TBD

## Contact

For questions or collaborations, please open an issue or contact the repository owner. TBD