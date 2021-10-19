# CODA: Calibrated Optimal Decision Making with Multiple Data Sources and Limited Outcome

This repository is the official implementation of [CODA: Calibrated Optimal Decision Making with Multiple Data Sources and Limited Outcome](https://arxiv.org/abs/2104.10554).

## Introduction

We consider the optimal decision-making problem in a primary sample of interest with multiple auxiliary sources available. The outcome of interest is limited in the sense that it is only observed in the primary sample. In reality, such multiple data sources may belong to heterogeneous studies and thus cannot be combined directly. This paper proposes a new framework to handle heterogeneous studies and address the limited outcome simultaneously through a novel calibrated optimal decision making (CODA) method, by leveraging the common intermediate outcomes in multiple data sources. Specifically, CODA allows the baseline covariates across different samples to have either homogeneous or heterogeneous distributions. Under a mild and testable assumption that the conditional means of intermediate outcomes in different samples are equal given baseline covariates and the treatment information, we show that the proposed CODA estimator of the conditional mean outcome is asymptotically normal and more efficient than using the primary sample solely. In addition, the variance of the CODA estimator can be easily obtained using the simple plug-in method due to the rate double robustness. 

## Requirements

 - R 3.6
 - `foreach`
 - `doParallel`
 - `policytree`
 - `randomForestSRC` 
 - `MASS` 

## Contents

  1. `README.txt`: implementation details of source code and contents 

  2. Source codes of CODA and data generation environment

     a). `CODA-HO.R`: main function for CODA with Homogeneous Baseline Covariates (CODA-HO);

     b). `CODA-HE.R`: main function for CODA with Heterogeneous Baseline Covariates (CODA-HE).
 
