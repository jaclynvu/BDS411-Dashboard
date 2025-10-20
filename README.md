# BDS411-Dashboard

Co-Authors: Morgan K Anderson, Connor Michael Eck, Daniel Hickey, Hung Phan, & Jaclyn Vu

This is a project from our Biological Data Science course: Analysis of Biological Data: Case Studies

We worked as a team to build a ShinyApp to analyze plant single cell data and describe its function in a report which includes your data visualizations. 

This repository will contain our app and my individual report. 

Co-Authors: Morgan K. Anderson, Connor Michael Eck, Daniel Hickey, Hung Phan, & Jaclyn Vu

Course: Biological Data Science (BDS 411) – Analysis of Biological Data: Case Studies

# Overview

This repository contains a ShinyApp and accompanying report developed for the course Analysis of Biological Data: Case Studies. The project focuses on visualizing and comparing single-cell RNA sequencing (scRNA-seq) data from Arabidopsis thaliana root tissue.

Our team created an interactive R Shiny dashboard that allows users to explore UMAP and PCA projections, examine correlations between clusters, and compare known versus predicted biomarkers.

# Objectives

The client sought an application to:

- Visualize UMAP and PCA projections of single-cell data.

- Quantitatively and visually compare similarities among cell clusters.

- Display cell counts per cluster and Pearson correlations.

- Explore predicted and known biomarkers interactively.

# Data

Five .csv files were provided:

File	Description
- AT_root_scRNAseq_UMAP.csv	Cell barcodes with UMAP coordinates
- AT_root_scRNAseq_PCA.csv	Cell barcodes with 50 principal components
- AT_root_scRNAseq_clustermap.csv	Cluster IDs for 3,355 cells (0–20)
- AT_root_scRNAseq_putativemarkers.csv	Known gene markers
- AT_root_scRNAseq_anticipatedmarkers.csv	Predicted gene markers

# Features

UMAP vs PCA Tab:

Compare UMAP and PCA side-by-side with consistent color palettes. Select which clusters to display and choose PCA components (PC1–PC50).

Bar Chart & Heatmap Tab:

View the number of cells per cluster and a Pearson correlation heatmap of average PCA values.

Marker Expression Tab:

Select predicted or known markers to visualize expression across clusters on the UMAP projection.

# Implementation

The app uses R Shiny with ggplot2, dplyr, reshape2, and corrplot.
Each team member implemented a functional component—Jaclyn Vu led the UMAP and PCA visualization modules and bar chart development.

To run locally:

git clone https://github.com/<your-username>/BDS411-Dashboard.git

cd BDS411-Dashboard

In R:

install.packages(c("shiny", "ggplot2", "dplyr", "corrplot", "reshape2"))

shiny::runApp("app.R")


Ensure the CSV files are in the same directory as app.R.

# Future Work

- Add zooming and rotation for interactive plots.

- Simplify PCA component selection (e.g., sliders).

- Enable user-uploaded scRNA-seq datasets for custom exploration.

# References

Jovic et al., 2022. Single-cell RNA sequencing technologies and applications. PMC8964935

Shahan et al., 2022. A single-cell Arabidopsis root atlas reveals developmental trajectories. ScienceDirect

National Cancer Institute, 2023. Definition of biomarker. NCI Dictionary
