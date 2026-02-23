# Dihadron

## Description
The basic post-processing framework for TPC-TPC di-hadron correlations, written by Zhiyong (zhiyong.lu@cern.ch).

This framework works for the output from the di-hadron-cor task (O2Physics/PWGCF/TwoParticleCorrelations/Tasks/diHadronCor.cxx)

## Usage
The sequence is: 
- remember to change the input directory of AnalysisResults.root in Process_dPhidEta.cxx
- root -l Process_dPhidEta.cxx
- root -l Process_CreateBootstrapSample.cxx
- root -l Process_TemplateFit.cxx

Note: This repo focuses on eta-dependent dihadron analysis.
