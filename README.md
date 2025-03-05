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

### Aggregate Bayesian Regression Model (ABRM)

`abrm_sim.R` and `nimble_abrm_0219.R` implement Bayesian approaches to spatial misalignment that account for the underlying spatial structure and uncertainty in the data.

Benefits of the Bayesian approach:
- Properly accounts for uncertainty in the redistribution process
- Can incorporate spatial correlation structures
- Allows for the integration of multiple data sources
- Provides full posterior distributions of estimates

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

## Contact

For questions or collaborations, please open an issue or contact the repository owner.