#!/bin/bash

root -q -l -b Process_dPhidEta.cxx
root -q -l -b Process_CreateBootstrapSample.cxx
root -q -l -b Process_TemplateFit.cxx
root -q -l -b Process_Uncorrected_FourierFit.cxx
root -q -l -b Process_Bootstrap_FourierFit.cxx
root -q -l -b Process_CompareVnMethods.cxx
root -q -l -b Process_ReconstructV2.cxx
root -q -l -b Process_ReconstructV2_TPC_B.cxx
root -q -l -b Process_3x2PC.cxx
root -q -l -b Process_CompareAll3Methods.cxx
root -q -l -b Process_FlowDecorrelation.cxx
