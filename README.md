# Spatial Misalignment Analysis with ABRM

This repository contains code and data for analyzing spatially misaligned data using Atom-Based Regression Models (ABRM) and comparing them with traditional dasymetric mapping approaches.

## Interactive Dashboard

An interactive dashboard visualizing the atom-based regression models is available at:
[https://bellayqian.shinyapps.io/abrm-dashboard/](https://bellayqian.shinyapps.io/abrm-dashboard/)

The dashboard allows you to explore:
- County boundaries
- Congressional district boundaries
- Atom counts and distributions
- Data tables for counties and atoms
- Method comparison results for Poisson, Binomial, and Normal distributions

## Repository Structure

- **poisson_analysis**: R scripts for statistical analysis comparing ABRM with dasymetric mapping
- **dashboard**: Interactive Shiny dashboard for visualizing counties, congressional districts, and atoms
- **data_preprocessing**: Scripts for preparing and preprocessing data

## Funding and Project Information

This work was funded by the Robert Wood Johnson Foundation, Grant 81746. Project details are provided below.

**Project Title:** Aligning spatially misaligned data for health equity analysis, action, and accountability

**Principal Investigators:** Dr. Nancy Krieger (PI) and Dr. Rachel Nethery (co-PI)

**Start Date:** July 2024

**Project Team and Collaborators:**
- Yunzhe Qian (Bella), MS (Research Assistant, Dept of Biostatistics, HSPH)
- Rachel Nethery, PhD (Assistant Professor, Dept of Biostatistics, HSPH)
- Nancy Krieger, PhD (Professor, Department of Social and Behavioral Sciences (SBS), HSPH)
- Nykesha Johnson, MPH (Statistical Data Analyst/Data Manager, SBS, HSPH)

## About

This work is an extension of:

Nethery, Rachel C., et al. "Addressing spatial misalignment in population health research: a case study of US congressional district political metrics and county health data." *MedRxiv* (2023).

Spatial misalignment—which occurs when data on multiple variables are collected using mismatched geographic boundary definitions—is a longstanding challenge in public health research. For instance, congressional districts can cut across multiple counties, and environmental hazard zones may cross census tract boundaries, in both cases creating intersecting areas that complicate efforts to study the relationships between health outcomes and their social, political, and environmental determinants.

Atom-based regression models (ABRM) offer a promising alternative by using atoms, the intersecting areas of all relevant units, as the fundamental units of analysis. By preserving the original spatial resolution of the data, ABRM account for uncertainty in statistical relationships while offering a robust method for handling misaligned data.

## Project Goal

To address these challenges, our work focuses on creating a comprehensive computational framework for ABRM, including developing open-source tools in R, to make ABRM more accessible to government agencies, health institutions, academic public health researchers, policy makers, and community-based organizations that produce, analyze, and/or use data from spatially misaligned geographic areas to inform their work for health equity.

## Citation

If you use this code or methodology in your research, please consider citing the paper above.

## Contact

**Yunzhe Qian (Bella):** Research Assistant, Dept. of Biostatistics, Harvard T.H. Chan School of Public Health

**Rachel Nethery (co-PI):** Assistant Professor of Biostatistics, Harvard T.H. Chan School of Public Health

**Nancy Krieger (PI):** Professor, Department of Social and Behavioral Sciences, Harvard T.H. Chan School of Public Health