---
title: "Deep diversification of an AAV capsid protein by machine learning"
date: 2021-09-10
lastmod: 2021-09-22
draft: false
tags: ["AAV capsid protein variants"]
categories: ["AAV capsid protein variants","Machine learning"]
author: "Xiaoni"

weight: 1

mathjax: true

# menu:
#   main:
#     parent: "docs"
#     weight: 1
---

Modern experimental technologies can assay large numbers of biological sequences, but engineered protein libraries rarely exceed the sequence diversity of natural protein families. Machine learning (ML) models trained directly on experimental data without biophysical modeling provide one route to accessing the full potential diversity of engineered proteins.

<!--more-->

[![fig1](fig1.png)](https://www.nature.com/articles/s41587-020-00793-4)

## Abstract

Apply deep learning to design highly diverse AAV2 capsid protein variants that remain viable for packaging of a DNA payload. Focusing on **a 28-AAA segment**, the present study generated **201,426** variants of the AAV2 WT sequence yielding **110,689** viable engineered capsids, **57,348** of which surpass the average diversity of natural AAV serotype sequences, with 12–29 mutations across this region. Even when trained on limited data, **deep neural network models** accurately predict capsid viability across diverse variants. 

## Significance

***This approach unlocks vast areas of functional but previously unreachable sequence space, with many potential applications for the generation of improved viral vectors and protein therapeutics.***

## Data availability

Experimental data for all three experiments are available at NCBI SRA [#PRJNA673640](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA673640/).

## Code availability

- The **TensorFlow 1.3**  API was used to implement and train all models.
- To construct the A39 training data, synthesize, process and analyze the experimental data with code [@github](https://github.com/churchlab/Deep_diversification_AAV).
- To reproduce the analysis figures from the main text with [ipython notebooks](https://github.com/google-research/google-research/tree/master/aav).



