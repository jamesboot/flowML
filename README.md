# flowML
README last updated: 11/02/2024
## Description
Analysis of flow cytometry data and application of ML algorithms to classify events.
## 🎯 Aim
- Analyse flow cytometry data in a high throughput manner
- Explore flow cytometry data using exploratory techniques such as tSNE and cluster analysis
- Apply ML algorithms to create a classifier of flow cytometry events
## 🖼️ Background
- Flow cytometry data of platelets treated with different agonists
- Can exploratory techniques such as tSNE and cluster analysis distinguish treatment groups?
- Can ML algorithms be trained to identify platelets treated with agonists vs untreated?
## 📝 Getting started
### Dependencies
```
# Load packages
library('stringr')
library('ggplot2')
library('dplyr')
library('limma')
library('reshape2')
library('tidyverse')
library('corrplot')
library('cluster') 
library('factoextra')
library('flowCore')
library('FlowSOM')
library('data.table')
library('Spectre')
```
- All code and functions used for analysis are contained in `Analysis.R`
