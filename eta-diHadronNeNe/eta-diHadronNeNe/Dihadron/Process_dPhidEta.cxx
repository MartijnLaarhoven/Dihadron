/*
 * @Author: Zhiyong Lu (zhiyong.lu@cern.ch)  
 * @Date: 2025-05-15 21:14:52 
 * @Last Modified by: Zhiyong Lu
 * @Last Modified time: 2025-05-18 14:12:45
 */
//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "TList.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TProfile.h>
#include <TRandom3.h>
#include "TMath.h"
#include "THnSparse.h"
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include "./include/BasicForDihadron.h"

#ifndef USE_CORRCONTAINER
#define USE_CORRCONTAINER 0
#endif

struct InputUnit {
    std::string fileNameSuffix;
    Int_t corrType;
    Bool_t isNch;
    Bool_t isEtadiff;
    Int_t minRange;
    Int_t maxRange;
    Bool_t isMc;

    InputUnit(std::string _fileNameSuffix, Int_t _corrType, Bool_t _isNch, Bool_t _isEtadiff, Int_t _minRange, Int_t _maxRange, Bool_t _isMc=false) :
        fileNameSuffix(_fileNameSuffix), corrType(_corrType), isNch(_isNch), isEtadiff(_isEtadiff), minRange(_minRange), maxRange(_maxRange), isMc(_isMc) {}
};

void printAxesInfo(THnSparseF* sparseHist);
void Read_dPhidEta_givenRange(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Bool_t isMc=false, Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.);
void Read_dPhidEta_givenRange_EtaDiff(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax, Bool_t isMc=false, Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.);

// global variables
std::string collisionSystemName;
std::string additionalSuffix = "";

void Process_dPhidEta() {
    // 不显示窗口
    gROOT->SetBatch(kTRUE);

    std::vector<InputUnit> inputList;
    additionalSuffix = "";
    
    collisionSystemName = "Ne-Ne"; //Change this for OO

    //Ne-Ne

    // inputList.push_back(InputUnit("LHC25af_pass2_598673", kTPCFT0C, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_598673", kTPCFT0C, kCent, kEtaDiffOff, 80, 100));

    // inputList.push_back(InputUnit("LHC25af_pass2_598673", kTPCFT0A, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_598673", kTPCFT0A, kCent, kEtaDiffOff, 80, 100));

    // inputList.push_back(InputUnit("LHC25af_pass2_598673", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_598673", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));

    // inputList.push_back(InputUnit("LHC25af_pass2_598673", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25af_pass2_598673", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));

    //O-O
    // inputList.push_back(InputUnit("LHC25ae_pass2_598682", kTPCFT0C, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_598682", kTPCFT0C, kCent, kEtaDiffOff, 80, 100));

    // inputList.push_back(InputUnit("LHC25ae_pass2_598682", kTPCFT0A, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_598682", kTPCFT0A, kCent, kEtaDiffOff, 80, 100));

    // inputList.push_back(InputUnit("LHC25ae_pass2_598682", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_598682", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));

    // inputList.push_back(InputUnit("LHC25ae_pass2_598682", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_598682", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));

    // Ne-Ne outer ring (615818) - full outer ring, correlate with TPC edges
    inputList.push_back(InputUnit("LHC25af_pass2_615818", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_615818", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));

    inputList.push_back(InputUnit("LHC25af_pass2_615818", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_615818", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));

    // Ne-Ne inner ring (615817) - full inner ring, correlate with TPC edges
    inputList.push_back(InputUnit("LHC25af_pass2_615817", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_615817", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));

    inputList.push_back(InputUnit("LHC25af_pass2_615817", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_615817", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));

    // O-O outer ring (616549) - full outer ring, correlate with TPC edges
    inputList.push_back(InputUnit("LHC25ae_pass2_616549", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25ae_pass2_616549", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));

    inputList.push_back(InputUnit("LHC25ae_pass2_616549", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25ae_pass2_616549", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));

    // O-O inner ring (618685) - full inner ring, correlate with TPC edges
    inputList.push_back(InputUnit("LHC25ae_pass2_618685", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25ae_pass2_618685", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));

    inputList.push_back(InputUnit("LHC25ae_pass2_618685", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25ae_pass2_618685", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));

    // Ne-Ne inner ring (617826) - new dataset
    inputList.push_back(InputUnit("LHC25af_pass2_617826", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_617826", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));

    inputList.push_back(InputUnit("LHC25af_pass2_617826", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_617826", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));

    // Ne-Ne outer ring (617910) - new dataset
    inputList.push_back(InputUnit("LHC25af_pass2_617910", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_617910", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));

    inputList.push_back(InputUnit("LHC25af_pass2_617910", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_617910", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));

    for (auto input : inputList) {
        // Auto-detect collision system from dataset name
        if (input.fileNameSuffix.find("LHC25ae") != std::string::npos) {
            collisionSystemName = "O-O";
        } else if (input.fileNameSuffix.find("LHC25af") != std::string::npos) {
            collisionSystemName = "Ne-Ne";
        } else {
            collisionSystemName = "Unknown";  // Fallback
        }
        
        // Check if this is a ring dataset - these have no FT0 slicing
        bool isRingDataset = (input.fileNameSuffix.find("615817") != std::string::npos || 
                              input.fileNameSuffix.find("615818") != std::string::npos ||
                              input.fileNameSuffix.find("616549") != std::string::npos ||
                              input.fileNameSuffix.find("618685") != std::string::npos ||
                              input.fileNameSuffix.find("617826") != std::string::npos ||
                              input.fileNameSuffix.find("617910") != std::string::npos);
        
        // For ring datasets: use whole FT0 detector, correlate only with TPC edges
        if (isRingDataset && input.isEtadiff) {
            std::vector<std::pair<double,double>> tpcEdges = { std::make_pair(-0.8, -0.7), std::make_pair(-0.7, -0.6), std::make_pair(0.6, 0.7), std::make_pair(0.7, 0.8) };
            bool isInner = (input.fileNameSuffix.find("615817") != std::string::npos || 
                            input.fileNameSuffix.find("618685") != std::string::npos ||
                            input.fileNameSuffix.find("617826") != std::string::npos);
            std::string ringType = isInner ? "inner" : "outer";
            // Correlate with TPC edges without FT0 subdivision
            for (auto rng : tpcEdges) {
                double etaMin = rng.first;
                double etaMax = rng.second;
                std::cout << "Ring " << ringType << " dataset " << DihadronCorrTypeName[input.corrType] 
                          << " Cent[" << input.minRange << "," << input.maxRange << "] correlating with TPC eta [" 
                          << etaMin << ", " << etaMax << "]" << std::endl;
                Read_dPhidEta_givenRange_EtaDiff(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange, etaMin, etaMax, input.isMc);
                std::cout << "  Completed TPC range [" << etaMin << "," << etaMax << "]" << std::endl;
            }
            std::cout << "Finished ring dataset " << input.fileNameSuffix << std::endl;
        }
        // For original FT0-sliced datasets (598682): use FT0 slices
        else {
            // Get the appropriate FT0 eta slices for this detector type
            const std::vector<float>* ft0Slices = nullptr;
            if (input.corrType == kTPCFT0C) {
                ft0Slices = &etaFT0C;
            } else if (input.corrType == kTPCFT0A) {
                ft0Slices = &etaFT0A;
            }

            // Loop over FT0 slices
            if (ft0Slices) {
                for (int iFT0 = 0; iFT0 < ft0Slices->size() - 1; iFT0++) {
                    double ft0EtaMin = (*ft0Slices)[iFT0];
                    double ft0EtaMax = (*ft0Slices)[iFT0 + 1];

                    if (input.isEtadiff) {
                        // For FT0 slicing, correlate with two wider TPC ranges: [-0.8, -0.6] and [0.6, 0.8]
                        std::vector<std::pair<double,double>> tpcRanges = { std::make_pair(-0.8, -0.6), std::make_pair(0.6, 0.8) };
                        for (auto rng : tpcRanges) {
                            double etaMin = rng.first;
                            double etaMax = rng.second;
                            std::cout << "FT0 slice [" << ft0EtaMin << ", " << ft0EtaMax 
                                      << "] correlating with TPC eta [" << etaMin << ", " << etaMax << "]" << std::endl;
                            Read_dPhidEta_givenRange_EtaDiff(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange, etaMin, etaMax, input.isMc, ft0EtaMin, ft0EtaMax);
                        }
                    }
                    else {
                        Read_dPhidEta_givenRange(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange, input.isMc, ft0EtaMin, ft0EtaMax);
                    }
                }
            } else {
                // No FT0 slicing for other correlation types
                if (input.isEtadiff) {
                    for (int iEta = 0; iEta < etaTPC.size() - 1; iEta++) {
                        double etaMin = etaTPC[iEta];
                        double etaMax = etaTPC[iEta + 1];
                        Read_dPhidEta_givenRange_EtaDiff(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange, etaMin, etaMax, input.isMc);
                    }
                }
                else {
                    Read_dPhidEta_givenRange(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange, input.isMc);
                }
            }
        }
    }
}

void printAxesInfo(THnSparseF* sparseHist) {
    // print the information of the axes of THnSparseF
    if (!sparseHist) {
        std::cerr << "Error: Null histogram pointer!" << std::endl;
        return;
    }

    const Int_t nDims = sparseHist->GetNdimensions();
    std::cout << "Sparse histogram has " << nDims << " dimensions\n";

    for (Int_t iDim = 0; iDim < nDims; ++iDim) {
        TAxis* axis = sparseHist->GetAxis(iDim);
        if (!axis) {
        std::cerr << "Error: Axis " << iDim << " not found!" << std::endl;
        continue;
        }
        if (iDim == 0) {
            axis->SetRangeUser(0,0.001);
              // 获取当前显示的区间参数
            const int firstBin = axis->GetFirst();
            const int lastBin = axis->GetLast();
            const double firstEdge = axis->GetBinLowEdge(firstBin);
            const double lastEdge = axis->GetBinUpEdge(lastBin);

            std::cout << "=== Axis 0" << " (" << axis->GetTitle() << ") ===\n"
                        << "Current displayed range: [" << firstEdge << ", " << lastEdge << "]\n"
                        << "Using bins: " << firstBin << " (" << firstEdge << ") -> "
                        << lastBin << " (" << lastEdge << ")\n\n";
        }

        std::cout << "Axis " << iDim << ": \n"
                << "  Title: " << axis->GetTitle() << "\n"
                << "  Bins:  " << axis->GetNbins() << "\n"
                << "  Range: [" << axis->GetXmin() 
                << ", " << axis->GetXmax() << "]\n"
                << "---------------------------------\n";
    }
}




void Read_dPhidEta_givenRange(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Bool_t isMc=false, Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.) {
    //std::cout << Form("flow-decorrelation_%s%s_%d_%d/sameEvent_%s", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()) << std::endl;
    TFile *file = TFile::Open(Form("../AnalysisResultsROOTFiles/dihadronEtaCorr/AnalysisResults_%s.root", fileNameSuffix.c_str()), "READ");
    if (!file || file->IsZombie()) { //error
        std::cout << "Error: Cannot open file " << fileNameSuffix << std::endl;
        return;
    }

    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    // check if MCTrue folder is available
    if (isMc) {
        if (!file->Get(Form("flow-decorrelation_%s%s_%d_%d/MCTrue", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange))) {
            std::cerr << "Error: MCTrue folder not found for " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
    }
    else {
        if (file->Get(Form("flow-decorrelation_%s%s_%d_%d/MCTrue", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange))) {
            std::cerr << "Caution! you are using Reco or Data, but MCTrue folder is found for " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
    }

    // Check if this is a ring dataset - ring datasets don't use CorrelationContainer
    bool isRingDataset = (fileNameSuffix.find("615817") != std::string::npos ||
                          fileNameSuffix.find("615818") != std::string::npos ||
                          fileNameSuffix.find("616549") != std::string::npos ||
                          fileNameSuffix.find("618685") != std::string::npos ||
                          fileNameSuffix.find("617826") != std::string::npos ||
                          fileNameSuffix.find("617910") != std::string::npos);
    
    if (isRingDataset) {
        // Ring datasets not supported in non-EtaDiff mode - skip for now
        std::cout << "Skipping ring dataset " << fileNameSuffix << " in non-EtaDiff mode" << std::endl;
        file->Close();
        delete file;
        return;
    }

#if USE_CORRCONTAINER
    CorrelationContainer *same = (CorrelationContainer*)file->Get(Form("flow-decorrelation_%s%s_%d_%d/sameEvent_%s", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));
    CorrelationContainer *mixed = (CorrelationContainer*)file->Get(Form("flow-decorrelation_%s%s_%d_%d/mixedEvent_%s", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));
    THnSparseD *trig = nullptr;
    if (!isMc) {
        trig = (THnSparseD*)file->Get(Form("flow-decorrelation_%s%s_%d_%d/Trig_hist_%s", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));
    } else {
        trig = (THnSparseD*)file->Get(Form("flow-decorrelation_%s%s_%d_%d/MCTrue/MCTrig_hist_%s", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));
    }
    

    if (!same || !mixed || !trig) { //error
        std::cerr << "Error getting histograms for " << fileNameSuffix << " with " << splitName << " and range [" << minRange << ", " << maxRange << "]" << std::endl;
        file->Close();
        delete file;
        return;
    }

    Int_t step = CorrelationContainer::kCFStepReconstructed;
    if (isMc) step = CorrelationContainer::kCFStepAll;

    THnSparseF *sparSig = (THnSparseF*)same->getPairHist()->getTHn(step);
    THnSparseF *sparMix = (THnSparseF*)mixed->getPairHist()->getTHn(step);
    if (!sparSig || !sparMix) {
        std::cerr << "Error: Empty histogram step " << step << " for " << fileNameSuffix
                  << " (corrType=" << DihadronCorrTypeName[corrType] << "). Skipping." << std::endl;
        file->Close();
        delete file;
        return;
    }
    sparSig->SetName(Form("sameEvent_%i_%i", minRange, maxRange));
    sparMix->SetName(Form("mixedEvent_%i_%i", minRange, maxRange));
    trig->SetName(Form("Trig_hist_%i_%i", minRange, maxRange));

    // Determine FT0 eta coverage: use slice if provided, else use full detector coverage
    double ft0Min, ft0Max;
    if (ft0EtaMin > -900. && ft0EtaMax > -900.) {
        ft0Min = ft0EtaMin;
        ft0Max = ft0EtaMax;
    } else {
        ft0Min = FitEtaCoverage[corrType][0];
        ft0Max = FitEtaCoverage[corrType][1];
    }

    // Common axis settings for all samples
    double dEtaMin = etaTPC.at(0) - ft0Max;
    double dEtaMax = etaTPC.back() - ft0Min;
    sparSig->GetAxis(corrAxis_kPt_TPC_trig)->SetRangeUser(etaTPC.at(0)+0.001, etaTPC.back()-0.001); //And yet he sets the same range for them
    sparMix->GetAxis(corrAxis_kPt_TPC_trig)->SetRangeUser(etaTPC.at(0)+0.001, etaTPC.back()-0.001);
    sparSig->GetAxis(corrAxis_kdEtaTPCTPC)->SetRangeUser(dEtaMin, dEtaMax);
    sparMix->GetAxis(corrAxis_kdEtaTPCTPC)->SetRangeUser(dEtaMin, dEtaMax);
    sparSig->GetAxis(corrAxis_kPt_TPC_asso)->SetRangeUser(ft0Min+0.001, ft0Max-0.001);
    sparMix->GetAxis(corrAxis_kPt_TPC_asso)->SetRangeUser(ft0Min+0.001, ft0Max-0.001);
    std::cout << "FT0 axis range set to [" << ft0Min << ", " << ft0Max << "]" << std::endl;
    std::cout << "  Same event entries before slice: " << sparSig->GetEntries() << std::endl;
    trig->GetAxis(trigAxis_pT)->SetRangeUser(minPt+0.001, maxPt-0.001);

    // Create output file with FT0 slice info if applicable
    TString outFileName;
    if (ft0EtaMin > -900. && ft0EtaMax > -900.) {
        outFileName = Form("./ProcessOutput/Mixed_%s%s_%s_%i_%i_%s_FT0Eta_%0.1f_%0.1f.root", 
                          fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), 
                          int(minRange), int(maxRange), DihadronCorrTypeName[corrType].c_str(),
                          ft0EtaMin, ft0EtaMax);
    } else {
        outFileName = Form("./ProcessOutput/Mixed_%s%s_%s_%i_%i_%s.root", 
                          fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), 
                          int(minRange), int(maxRange), DihadronCorrTypeName[corrType].c_str());
    }
    TFile* fout = TFile::Open(outFileName, "RECREATE");

    // Process all samples: -1 (all), 0 to maxSample-1
    for (Int_t sample = -1; sample < maxSample; ++sample) { //This is the reason for why we create 10 histograms, maxSample=10 in basic for dihadron
        // Set sample range if applicable
        if (sample >= 0) {
            sparSig->GetAxis(corrAxis_kSample)->SetRangeUser(sample, sample+0.001);
            sparMix->GetAxis(corrAxis_kSample)->SetRangeUser(sample, sample+0.001);
            trig->GetAxis(trigAxis_sample)->SetRangeUser(sample, sample+0.001);
        } else {
            // Reset to full range for inclusive sample
            double xmin = sparSig->GetAxis(corrAxis_kSample)->GetXmin();
            double xmax = sparSig->GetAxis(corrAxis_kSample)->GetXmax();
            sparSig->GetAxis(corrAxis_kSample)->SetRangeUser(xmin, xmax);
            xmin = sparMix->GetAxis(corrAxis_kSample)->GetXmin();
            xmax = sparMix->GetAxis(corrAxis_kSample)->GetXmax();
            sparMix->GetAxis(corrAxis_kSample)->SetRangeUser(xmin, xmax);
            xmin = trig->GetAxis(trigAxis_sample)->GetXmin();
            xmax = trig->GetAxis(trigAxis_sample)->GetXmax();
            trig->GetAxis(trigAxis_sample)->SetRangeUser(xmin, xmax);
        }

        // Initialize variables for this sample
        TH2D* hPhiEtaSMsum = nullptr;
        TH2D* hPhiEtaSsum = nullptr;
        TH2D* hPhiEtaMsum = nullptr;
        Double_t nTriggersS = 0.;
        Int_t nz = sparSig->GetAxis(corrAxis_kVz)->GetNbins();
        // Set names with sample suffix
        TString suffix = (sample == -1) ? "" : Form("_%d", sample); //if sample==-1, suffix="", else suffix="_sample"

        // Vertex loop
        for (Int_t iz = 1; iz <= nz; ++iz) {
             // project the vertex axis into a 1D histogram of the trigger sparse
            TH1D* hTriggersS = (TH1D*)trig->Projection(trigAxis_Vz);
            // get the number of triggers in this vertex bin
            nTriggersS += hTriggersS->Integral(iz, iz);

            // and set the vertex range for the same and mixed event sparses
            sparSig->GetAxis(corrAxis_kVz)->SetRange(iz, iz);
            sparMix->GetAxis(corrAxis_kVz)->SetRange(iz, iz);

            //project the phi and eta axis into a 2D Histogram for both same and mixed
            TH2D *hPhiEtaS = (TH2D*)sparSig->Projection(corrAxis_kdEtaTPCTPC, corrAxis_kdPhiTPCTPC);
            TH2D *hPhiEtaM = (TH2D*)sparMix->Projection(corrAxis_kdEtaTPCTPC, corrAxis_kdPhiTPCTPC);
            // Normalize mixed event
            Double_t norm = 1.;
            double MixEventNormalizationEta = (etaTPC.at(0) + etaTPC.back()) / 2.0 - (FitEtaCoverage[corrType][0] + FitEtaCoverage[corrType][1]) / 2.0;
            Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(-TMath::Pi()/2 + 0.0001);
            Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(3*TMath::Pi()/2 - 0.0001);
            Int_t binEta1 = hPhiEtaM->GetYaxis()->FindBin(MixEventNormalizationEta);
            Int_t binEta2 = hPhiEtaM->GetYaxis()->FindBin(MixEventNormalizationEta);
            Int_t nNormBins = (binEta2 - binEta1 + 1) * (binPhi2 - binPhi1 + 1);
            norm = hPhiEtaM->Integral(binPhi1, binPhi2, binEta1, binEta2) / nNormBins;

            if (!hPhiEtaMsum) {
                hPhiEtaMsum = (TH2D*)hPhiEtaM->Clone(Form("dphideta_ME_%d_%d%s", minRange, maxRange, suffix.Data())); //seems there are 10 samples already in from diHadon.cxx... makes sense as forming fake full sample requires correlating smaller number of elements
            } else {
                hPhiEtaMsum->Add(hPhiEtaM);
            }

            hPhiEtaM->Scale(1.0 / norm);

            TH2D* hPhiEtaSM = (TH2D*)hPhiEtaS->Clone(Form("dphideta_SM_%d_%d_%d%s", minRange, maxRange, iz, suffix.Data()));
            hPhiEtaSM->Divide(hPhiEtaM); //her dividerer vi S med B

            if (!hPhiEtaSMsum) {
                hPhiEtaSMsum = (TH2D*)hPhiEtaSM->Clone(Form("dphideta_SM_%d_%d%s", minRange, maxRange, suffix.Data()));
            } else {
                hPhiEtaSMsum->Add(hPhiEtaSM);
            }

            if (!hPhiEtaSsum) {
                hPhiEtaSsum = (TH2D*)hPhiEtaS->Clone(Form("dphideta_SE_%d_%d%s", minRange, maxRange, suffix.Data()));
            } else {
                hPhiEtaSsum->Add(hPhiEtaS);
            }
        }

        // Normalization and final processing
        if (nTriggersS > 0) {
            hPhiEtaSMsum->Scale(1.0 / nTriggersS);
            hPhiEtaSsum->Scale(1.0 / nTriggersS);
        }
        hPhiEtaSMsum->Scale(1.0 / hPhiEtaSMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSMsum->Scale(1.0 / hPhiEtaSMsum->GetYaxis()->GetBinWidth(1));
        TH1D* hEta = hPhiEtaSMsum->ProjectionY(Form("hEta_%d_%d%s", minRange, maxRange, suffix.Data()));
        hEta->SetTitle("#Delta#eta");

        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaMsum->Rebin2D(1, 1);

        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaSsum->Rebin2D(1, 1);

        hPhiEtaSMsum->SetName(Form("dphideta_SM_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiEtaSsum->SetName(Form("dphideta_SE_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiEtaMsum->SetName(Form("dphideta_ME_%d_%d%s", minRange, maxRange, suffix.Data()));

        // Set axis titles and labels for hPhiEtaSMsum
        hPhiEtaSMsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaSMsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaSMsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaSMsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaSMsum->SetTitle(Form("Correlation Function %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaSMsum->SetTitle(Form("Correlation Function %d< Centrality #leq%d", minRange, maxRange));


        // Set axis titles and labels for hPhiEtaSsum
        hPhiEtaSsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaSsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaSsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaSsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaSsum->SetTitle(Form("Same Event %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaSsum->SetTitle(Form("Same Event %d< Centrality #leq%d", minRange, maxRange));

        // Set axis titles and labels for hPhiEtaMsum
        hPhiEtaMsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaMsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaMsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaMsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaMsum->SetTitle(Form("Mixed event %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaMsum->SetTitle(Form("Mixed event %d< Centrality #leq%d", minRange, maxRange));

        // Draw histograms
        TCanvas* c1 = new TCanvas(Form("dPhidEta %s", suffix.Data()), Form("dPhidEta %s", suffix.Data()), 1200, 800);
        c1->Divide(2, 2);

        c1->cd(1);
        TH2D* hPhiEtaSMsum_draw = (TH2D*)hPhiEtaSMsum->Clone("hPhiEtaSMsum_draw");
        hPhiEtaSMsum_draw->Draw("surf1");
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s %s", collisionSystemName.c_str(), DihadronCorrTypeName[corrType].c_str()));
        latex.DrawLatex(0.1, 0.85, Form("#eta^{TPC} #in [%0.1f, %0.1f], #eta^{FIT} #in [%0.1f, %0.1f]", etaTPC.at(0), etaTPC.back(), ft0Min, ft0Max));


        c1->cd(3);
        hPhiEtaSsum->Draw("surf1");
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s %s", collisionSystemName.c_str(), DihadronCorrTypeName[corrType].c_str()));
        latex.DrawLatex(0.1, 0.85, Form("#eta^{TPC} #in [%0.1f, %0.1f], #eta^{FIT} #in [%0.1f, %0.1f]", etaTPC.at(0), etaTPC.back(), ft0Min, ft0Max));


        c1->cd(4);
        hPhiEtaMsum->Draw("surf1");
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s %s", collisionSystemName.c_str(), DihadronCorrTypeName[corrType].c_str()));
        latex.DrawLatex(0.1, 0.85, Form("#eta^{TPC} #in [%0.1f, %0.1f], #eta^{FIT} #in [%0.1f, %0.1f]", etaTPC.at(0), etaTPC.back(), ft0Min, ft0Max));

        // Saving to file
        fout->cd();
        c1->Write();
        // write canvas to file
        // c1->SaveAs(Form("./ProcessOutput/Mixed_%s_%s_%i_%i.png", fileNameSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange)));
        // delete c1;

        // Write histograms
        hPhiEtaSMsum->Write();
        hPhiEtaSsum->Write();
        hPhiEtaMsum->Write();
        hEta->Write();

        // Eta Gap processing
        TH1D* hPhiSameOverMixed_pos = hPhiEtaSMsum->ProjectionX("hPhiSameOverMixed");
        hPhiSameOverMixed_pos->SetName(Form("hPhiSameOverMixed_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiSameOverMixed_pos->SetTitle(Form("hPhiSameOverMixed_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiSameOverMixed_pos->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiSameOverMixed_pos->Write();


        delete hPhiEtaSMsum;
        delete hPhiEtaSMsum_draw;
        delete hPhiEtaSsum;
        delete hPhiEtaMsum;
        delete hPhiSameOverMixed_pos;
        delete c1;
        
    }

    fout->Close();
    delete fout;
    file->Close();
    delete file;
    delete sparSig;
    delete sparMix;
    std::cout << "Output file: " << Form("./ProcessOutput/Mixed_%s%s_%s_%i_%i_%s.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), DihadronCorrTypeName[corrType].c_str()) << std::endl; 
    std::cout << "Processing completed for all samples." << std::endl;
#else
    std::cerr << "Error: CorrelationContainer support is disabled. Cannot process non-EtaDiff data for " << fileNameSuffix << std::endl;
    file->Close();
    delete file;
    return;
#endif
}

void Read_dPhidEta_givenRange_EtaDiff(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax, Bool_t isMc=false, Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.) {
    TFile *file = TFile::Open(Form("../AnalysisResultsROOTFiles/dihadronEtaCorr/AnalysisResults_%s.root", fileNameSuffix.c_str()), "READ");
    if (!file || file->IsZombie()) {
        std::cout << "Error: Cannot open file " << fileNameSuffix << std::endl;
        return;
    }

    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    // check if MCTrue folder is available
    if (isMc) {
        if (!file->Get(Form("flow-decorrelation_%s%s_%d_%d/MCTrue", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange))) {
            std::cerr << "Error: MCTrue folder not found for " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
    }
    else {
        if (file->Get(Form("flow-decorrelation_%s%s_%d_%d/MCTrue", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange))) {
            std::cerr << "Caution! you are using Reco or Data, but MCTrue folder is found for " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
    }

    // Check if this is a ring dataset - they use direct TH2D histograms instead of CorrelationContainer
    bool isRingDataset = (fileNameSuffix.find("615817") != std::string::npos ||
                          fileNameSuffix.find("615818") != std::string::npos ||
                          fileNameSuffix.find("616549") != std::string::npos ||
                          fileNameSuffix.find("618685") != std::string::npos ||
                          fileNameSuffix.find("617826") != std::string::npos ||
                          fileNameSuffix.find("617910") != std::string::npos);
    
    if (isRingDataset) {
        // ===== RING DATASET PROCESSING (direct TH2D access) =====
        std::cout << "Processing ring dataset " << fileNameSuffix << " using direct TH2D histograms..." << std::endl;
        
        // Get direct TH2D histograms from the directory
        TDirectory* dir = (TDirectory*)file->Get(Form("flow-decorrelation_%s%s_%d_%d", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange));
        if (!dir) {
            std::cerr << "Error: Cannot access directory for ring dataset " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
        
        std::string detectorName = (corrType == kTPCFT0C) ? "TPC_FT0C" : "TPC_FT0A";
        TH2D* hSame = (TH2D*)dir->Get(Form("deltaEta_deltaPhi_same_%s", detectorName.c_str()));
        TH2D* hMixed = (TH2D*)dir->Get(Form("deltaEta_deltaPhi_mixed_%s", detectorName.c_str()));
        THnSparseD* trig = (THnSparseD*)dir->Get(Form("Trig_hist_%s", detectorName.c_str()));
        
        if (!hSame || !hMixed || !trig) {
            std::cerr << "Error: Cannot get direct histograms for ring dataset " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
        
        std::cout << "  Found TH2D histograms: " << hSame->GetNbinsX() << " x " << hSame->GetNbinsY() 
                  << " bins, entries=" << hSame->GetEntries() << std::endl;
        
        // Create output file
        TString outFileName = Form("./ProcessOutput/EtaDiff/Mixed_%s%s_%s_%i_%i_Eta_%0.1f_%0.1f_%s.root",
                                  fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(),
                                  int(minRange), int(maxRange), etaMin, etaMax, DihadronCorrTypeName[corrType].c_str());
        TFile* fout = TFile::Open(outFileName, "RECREATE");
        
        // For ring datasets, we need to slice the 2D histogram to the delta-eta range
        // corresponding to our specific TPC eta bin
        
        // Calculate delta-eta range for this TPC bin with the full FT0 detector
        // delta-eta = eta_TPC - eta_FT0
        Double_t ft0EtaMin = FitEtaCoverage[corrType][0];
        Double_t ft0EtaMax = FitEtaCoverage[corrType][1];
        
        // For the TPC bin [etaMin, etaMax], calculate expected delta-eta range
        Double_t deltaEtaMin = etaMin - ft0EtaMax;  // Most negative: smallest TPC - largest FT0
        Double_t deltaEtaMax = etaMax - ft0EtaMin;  // Most positive: largest TPC - smallest FT0
        
        // Avoid sparse histogram edge bins (bins 1-2 are empty for delta-eta < 1.6)
        // If deltaEtaMin falls below 1.6, start from 1.6 to avoid empty bins
        if (deltaEtaMin > 0 && deltaEtaMin < 1.6) {
            std::cout << "  Warning: deltaEtaMin " << deltaEtaMin << " in sparse region, adjusting to 1.6" << std::endl;
            deltaEtaMin = 1.6;
        }
        
        std::cout << "  TPC eta range: [" << etaMin << ", " << etaMax << "]" << std::endl;
        std::cout << "  FT0 coverage: [" << ft0EtaMin << ", " << ft0EtaMax << "]" << std::endl;
        std::cout << "  Expected delta-eta range: [" << deltaEtaMin << ", " << deltaEtaMax << "]" << std::endl;
        
        // Find bins in the 2D histogram corresponding to this delta-eta range
        TAxis* etaAxis = hSame->GetYaxis();
        Int_t binDeltaEtaMin = etaAxis->FindBin(deltaEtaMin);
        Int_t binDeltaEtaMax = etaAxis->FindBin(deltaEtaMax);
        
        // Ensure bins are within valid range
        if (binDeltaEtaMin < 1) binDeltaEtaMin = 1;
        if (binDeltaEtaMax > etaAxis->GetNbins()) binDeltaEtaMax = etaAxis->GetNbins();
        
        std::cout << "  Using delta-eta bins [" << binDeltaEtaMin << ", " << binDeltaEtaMax 
                  << "] corresponding to [" << etaAxis->GetBinLowEdge(binDeltaEtaMin) 
                  << ", " << etaAxis->GetBinUpEdge(binDeltaEtaMax) << "]" << std::endl;
        
        // Clone histograms for output (keeping the full 2D histograms for reference)
        TH2D* hPhiEtaSsum = (TH2D*)hSame->Clone(Form("dphideta_SE_%d_%d", minRange, maxRange));
        TH2D* hPhiEtaMsum = (TH2D*)hMixed->Clone(Form("dphideta_ME_%d_%d", minRange, maxRange));
        
        // Normalize mixed event
        Double_t norm = 1.;
        double MixEventNormalizationEta = (etaMin + etaMax) / 2.0 - (ft0EtaMin + ft0EtaMax) / 2.0;
        Int_t binPhi1 = hMixed->GetXaxis()->FindBin(-TMath::Pi()/2 + 0.0001);
        Int_t binPhi2 = hMixed->GetXaxis()->FindBin(3*TMath::Pi()/2 - 0.0001);
        
        // Check if normalization eta is within histogram bounds (reuse etaAxis from above)
        Double_t etaMin_hist = etaAxis->GetBinLowEdge(1);
        Double_t etaMax_hist = etaAxis->GetBinUpEdge(etaAxis->GetNbins());
        
        Int_t binEta1, binEta2;
        if (MixEventNormalizationEta < etaMin_hist || MixEventNormalizationEta > etaMax_hist) {
            // Normalization eta outside histogram range - use near-zero delta-eta region
            std::cout << "  Warning: Normalization eta " << MixEventNormalizationEta 
                      << " outside histogram range [" << etaMin_hist << ", " << etaMax_hist 
                      << "], using near-zero delta-eta region" << std::endl;
            binEta1 = etaAxis->FindBin(-0.2);
            binEta2 = etaAxis->FindBin(0.2);
            // Ensure bins are within valid range
            if (binEta1 < 1) binEta1 = 1;
            if (binEta2 > etaAxis->GetNbins()) binEta2 = etaAxis->GetNbins();
        } else {
            binEta1 = etaAxis->FindBin(MixEventNormalizationEta);
            binEta2 = etaAxis->FindBin(MixEventNormalizationEta);
        }
        
        Int_t nNormBins = (binEta2 - binEta1 + 1) * (binPhi2 - binPhi1 + 1);
        if (nNormBins > 0) {
            norm = hMixed->Integral(binPhi1, binPhi2, binEta1, binEta2) / nNormBins;
        }
        
        // Final fallback: if still no normalization, use full phi range and middle eta region
        if (norm <= 0) {
            std::cout << "  Warning: Normalization still zero, using middle eta bins as fallback" << std::endl;
            Int_t nBinsEta = etaAxis->GetNbins();
            binEta1 = nBinsEta / 2 - 1;
            binEta2 = nBinsEta / 2 + 1;
            if (binEta1 < 1) binEta1 = 1;
            if (binEta2 > nBinsEta) binEta2 = nBinsEta;
            nNormBins = (binEta2 - binEta1 + 1) * (binPhi2 - binPhi1 + 1);
            if (nNormBins > 0) {
                norm = hMixed->Integral(binPhi1, binPhi2, binEta1, binEta2) / nNormBins;
            }
        }
        
        TH2D* hMixedNorm = (TH2D*)hMixed->Clone("hMixedNorm");
        if (norm > 0) {
            hMixedNorm->Scale(1.0 / norm);
        }
        
        // Calculate correlation function (same/mixed) for the SLICED delta-eta range
        // First, create sliced histograms restricted to our delta-eta range
        TH2D* hSameSliced = (TH2D*)hSame->Clone("hSameSliced");
        TH2D* hMixedNormSliced = (TH2D*)hMixedNorm->Clone("hMixedNormSliced");
        
        // Zero out bins outside our delta-eta range
        for (Int_t iBinEta = 1; iBinEta <= etaAxis->GetNbins(); ++iBinEta) {
            if (iBinEta < binDeltaEtaMin || iBinEta > binDeltaEtaMax) {
                for (Int_t iBinPhi = 1; iBinPhi <= hSameSliced->GetNbinsX(); ++iBinPhi) {
                    hSameSliced->SetBinContent(iBinPhi, iBinEta, 0.0);
                    hSameSliced->SetBinError(iBinPhi, iBinEta, 0.0);
                    hMixedNormSliced->SetBinContent(iBinPhi, iBinEta, 0.0);
                    hMixedNormSliced->SetBinError(iBinPhi, iBinEta, 0.0);
                }
            }
        }
        
        // Calculate correlation function (same/mixed) using sliced histograms
        TH2D* hPhiEtaSM = (TH2D*)hSameSliced->Clone(Form("dphideta_SM_%d_%d", minRange, maxRange));
        hPhiEtaSM->Divide(hMixedNormSliced);
        
        // Get trigger count
        Double_t nTriggers = trig->Projection(0)->Integral(); // Project to first axis and integrate
        if (nTriggers > 0) {
            hPhiEtaSM->Scale(1.0 / nTriggers);
            hPhiEtaSsum->Scale(1.0 / nTriggers);
        }

        // Create and write 1D projections needed for bootstrap
        TH1D* hPhiSameOverMixed = hPhiEtaSM->ProjectionX("hPhiSameOverMixed");
        hPhiSameOverMixed->SetName(Form("hPhiSameOverMixed_%d_%d", minRange, maxRange));
        hPhiSameOverMixed->SetTitle(Form("hPhiSameOverMixed_%d_%d", minRange, maxRange));
        hPhiSameOverMixed->GetXaxis()->SetTitle("#Delta#varphi");
        
        // Write histograms to output
        fout->cd();
        hPhiEtaSsum->Write();
        hPhiEtaMsum->Write();
        hPhiEtaSM->Write();
        hPhiSameOverMixed->Write();

        for (Int_t sample = 0; sample < maxSample; ++sample) {
            TH1D* hSample = (TH1D*)hPhiSameOverMixed->Clone(Form("hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, sample));
            hSample->SetTitle(Form("hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, sample));
            hSample->Write();
            delete hSample;
        }
        
        std::cout << "  Ring dataset processed: nTriggers=" << nTriggers << ", norm=" << norm << std::endl;
        std::cout << "  Output written to: " << outFileName << std::endl;
        
        // Cleanup temporary histograms
        delete hSameSliced;
        delete hMixedNormSliced;
        delete hMixedNorm;
        delete hPhiEtaSM;
        delete hPhiEtaSsum;
        delete hPhiEtaMsum;
        delete hPhiSameOverMixed;
        
        fout->Close();
        delete fout;
        file->Close();
        delete file;
        return;
    }

#if USE_CORRCONTAINER
    CorrelationContainer *same = (CorrelationContainer*)file->Get(Form("flow-decorrelation_%s%s_%d_%d/sameEvent_%s", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));
    CorrelationContainer *mixed = (CorrelationContainer*)file->Get(Form("flow-decorrelation_%s%s_%d_%d/mixedEvent_%s", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));

    THnSparseD *trig = nullptr;
    if (!isMc) {
        trig = (THnSparseD*)file->Get(Form("flow-decorrelation_%s%s_%d_%d/Trig_hist_%s", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));
    } else {
        trig = (THnSparseD*)file->Get(Form("flow-decorrelation_%s%s_%d_%d/MCTrue/MCTrig_hist_%s", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str()));
    }

    // ===== STANDARD DATASET PROCESSING (CorrelationContainer access) =====
    if (!same || !mixed || !trig) {
        std::cerr << "Error getting histograms for " << fileNameSuffix << " with " << splitName << " and range [" << minRange << ", " << maxRange << "]" << std::endl;
        file->Close();
        delete file;
        return;
    }

    Int_t step = CorrelationContainer::kCFStepReconstructed;
    if (isMc) step = CorrelationContainer::kCFStepAll;
    
    THnSparseF *sparSig = (THnSparseF*)same->getPairHist()->getTHn(step);
    THnSparseF *sparMix = (THnSparseF*)mixed->getPairHist()->getTHn(step);
    
    if (!sparSig || !sparMix) {
        std::cerr << "Error: No valid histogram steps found for " << fileNameSuffix << ". Skipping this dataset." << std::endl;
        file->Close();
        delete file;
        return;
    }

    sparSig->SetName(Form("sameEvent_%i_%i", minRange, maxRange));
    sparMix->SetName(Form("mixedEvent_%i_%i", minRange, maxRange));
    trig->SetName(Form("Trig_hist_%i_%i", minRange, maxRange));

    // Verify sparse histograms have required dimensions
    Int_t nDimSig = sparSig->GetNdimensions();
    Int_t nDimMix = sparMix->GetNdimensions();
    std::cout << "Sparse histogram dimensions: same=" << nDimSig << ", mixed=" << nDimMix << std::endl;
    
    if (nDimSig < 4 || nDimMix < 4) {
        std::cerr << "Error: Sparse histograms don't have enough dimensions for " << fileNameSuffix 
                  << " (same=" << nDimSig << ", mixed=" << nDimMix << ", need >= 4). Skipping." << std::endl;
        file->Close();
        delete file;
        return;
    }

    // Determine FT0 eta coverage: use slice if provided, else use full detector coverage
    double ft0Min, ft0Max;
    if (ft0EtaMin > -900. && ft0EtaMax > -900.) {
        ft0Min = ft0EtaMin;
        ft0Max = ft0EtaMax;
    } else {
        ft0Min = FitEtaCoverage[corrType][0];
        ft0Max = FitEtaCoverage[corrType][1];
    }

    double dEtaMin = etaMin - ft0Max;
    double dEtaMax = etaMax - ft0Min;
    
    try {
        sparSig->GetAxis(corrAxis_kPt_TPC_trig)->SetRangeUser(etaMin+0.001, etaMax-0.001); //Here they differ
        sparMix->GetAxis(corrAxis_kPt_TPC_trig)->SetRangeUser(etaMin+0.001, etaMax-0.001);
        sparSig->GetAxis(corrAxis_kdEtaTPCTPC)->SetRangeUser(dEtaMin, dEtaMax);
        sparMix->GetAxis(corrAxis_kdEtaTPCTPC)->SetRangeUser(dEtaMin, dEtaMax);
        sparSig->GetAxis(corrAxis_kPt_TPC_asso)->SetRangeUser(ft0Min+0.001, ft0Max-0.001);
        sparMix->GetAxis(corrAxis_kPt_TPC_asso)->SetRangeUser(ft0Min+0.001, ft0Max-0.001);
    } catch (const std::exception& e) {
        std::cerr << "Error setting axis ranges for " << fileNameSuffix << ": " << e.what() << ". Skipping." << std::endl;
        file->Close();
        delete file;
        return;
    }
    
    // Verify and set trigger histogram axis range
    try {
        if (!trig || trig->GetNdimensions() < 1) {
            std::cerr << "Error: Trigger histogram missing or has insufficient dimensions for " << fileNameSuffix << ". Skipping." << std::endl;
            file->Close();
            delete file;
            return;
        }
        trig->GetAxis(trigAxis_pT)->SetRangeUser(minPt+0.001, maxPt-0.001);
    } catch (const std::exception& e) {
        std::cerr << "Error setting trigger axis range for " << fileNameSuffix << ": " << e.what() << ". Skipping." << std::endl;
        file->Close();
        delete file;
        return;
    }

    // Create output file with FT0 slice info if applicable
    TString outFileName;
    if (ft0EtaMin > -900. && ft0EtaMax > -900.) {
        outFileName = Form("./ProcessOutput/EtaDiff/Mixed_%s%s_%s_%i_%i_Eta_%0.1f_%0.1f_%s_FT0Eta_%0.1f_%0.1f.root",
                          fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(),
                          int(minRange), int(maxRange), etaMin, etaMax,
                          DihadronCorrTypeName[corrType].c_str(), ft0EtaMin, ft0EtaMax);
    } else {
        outFileName = Form("./ProcessOutput/EtaDiff/Mixed_%s%s_%s_%i_%i_Eta_%0.1f_%0.1f_%s.root",
                          fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(),
                          int(minRange), int(maxRange), etaMin, etaMax, DihadronCorrTypeName[corrType].c_str());
    }
    TFile* fout = TFile::Open(outFileName, "RECREATE");

    // Process all samples: -1 (all), 0 to maxSample-1
    for (Int_t sample = -1; sample < maxSample; ++sample) {
        // Set sample range if applicable
        if (sample >= 0) {
            sparSig->GetAxis(corrAxis_kSample)->SetRangeUser(sample, sample+0.001);
            sparMix->GetAxis(corrAxis_kSample)->SetRangeUser(sample, sample+0.001);
            trig->GetAxis(trigAxis_sample)->SetRangeUser(sample, sample+0.001);
        } else {
            // Reset to full range for inclusive sample
            double xmin = sparSig->GetAxis(corrAxis_kSample)->GetXmin();
            double xmax = sparSig->GetAxis(corrAxis_kSample)->GetXmax();
            sparSig->GetAxis(corrAxis_kSample)->SetRangeUser(xmin, xmax);
            xmin = sparMix->GetAxis(corrAxis_kSample)->GetXmin();
            xmax = sparMix->GetAxis(corrAxis_kSample)->GetXmax();
            sparMix->GetAxis(corrAxis_kSample)->SetRangeUser(xmin, xmax);
            xmin = trig->GetAxis(trigAxis_sample)->GetXmin();
            xmax = trig->GetAxis(trigAxis_sample)->GetXmax();
            trig->GetAxis(trigAxis_sample)->SetRangeUser(xmin, xmax);
        }
        
        // Initialize variables for this sample
        TH2D* hPhiEtaSMsum = nullptr;
        TH2D* hPhiEtaSsum = nullptr;
        TH2D* hPhiEtaMsum = nullptr;
        Double_t nTriggersS = 0.;
        Int_t nz = sparSig->GetAxis(corrAxis_kVz)->GetNbins();
        // Set names with sample suffix
        TString suffix = (sample == -1) ? "" : Form("_%d", sample);

        // Vertex loop
        for (Int_t iz = 1; iz <= nz; ++iz) {
             // project the vertex axis into a 1D histogram of the trigger sparse
            TH1D* hTriggersS = (TH1D*)trig->Projection(trigAxis_Vz);
            // get the number of triggers in this vertex bin
            nTriggersS += hTriggersS->Integral(iz, iz);

            // and set the vertex range for the same and mixed event sparses
            sparSig->GetAxis(corrAxis_kVz)->SetRange(iz, iz);
            sparMix->GetAxis(corrAxis_kVz)->SetRange(iz, iz);

            //project the phi and eta axis into a 2D Histogram for both same and mixed
            TH2D *hPhiEtaS = (TH2D*)sparSig->Projection(corrAxis_kdEtaTPCTPC, corrAxis_kdPhiTPCTPC);
            TH2D *hPhiEtaM = (TH2D*)sparMix->Projection(corrAxis_kdEtaTPCTPC, corrAxis_kdPhiTPCTPC);

            // Normalize mixed event
            Double_t norm = 1.;
            double MixEventNormalizationEta = (etaMin + etaMax) / 2.0 - (FitEtaCoverage[corrType][0] + FitEtaCoverage[corrType][1]) / 2.0;
            Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(-TMath::Pi()/2 + 0.0001);
            Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(3*TMath::Pi()/2 - 0.0001);
            Int_t binEta1 = hPhiEtaM->GetYaxis()->FindBin(MixEventNormalizationEta);
            Int_t binEta2 = hPhiEtaM->GetYaxis()->FindBin(MixEventNormalizationEta);
            Int_t nNormBins = (binEta2 - binEta1 + 1) * (binPhi2 - binPhi1 + 1);
            norm = hPhiEtaM->Integral(binPhi1, binPhi2, binEta1, binEta2) / nNormBins;

            if (!hPhiEtaMsum) {
                hPhiEtaMsum = (TH2D*)hPhiEtaM->Clone(Form("dphideta_ME_%d_%d%s", minRange, maxRange, suffix.Data()));
            } else {
                hPhiEtaMsum->Add(hPhiEtaM);
            }

            hPhiEtaM->Scale(1.0 / norm);

            TH2D* hPhiEtaSM = (TH2D*)hPhiEtaS->Clone(Form("dphideta_SM_%d_%d_%d%s", minRange, maxRange, iz, suffix.Data()));
            hPhiEtaSM->Divide(hPhiEtaM);

            if (!hPhiEtaSMsum) {
                hPhiEtaSMsum = (TH2D*)hPhiEtaSM->Clone(Form("dphideta_SM_%d_%d%s", minRange, maxRange, suffix.Data()));
            } else {
                hPhiEtaSMsum->Add(hPhiEtaSM);
            }

            if (!hPhiEtaSsum) {
                hPhiEtaSsum = (TH2D*)hPhiEtaS->Clone(Form("dphideta_SE_%d_%d%s", minRange, maxRange, suffix.Data()));
            } else {
                hPhiEtaSsum->Add(hPhiEtaS);
            }
        }

        // Normalization and final processing
        if (nTriggersS > 0) {
            hPhiEtaSMsum->Scale(1.0 / nTriggersS);
            hPhiEtaSsum->Scale(1.0 / nTriggersS);
        }
        hPhiEtaSMsum->Scale(1.0 / hPhiEtaSMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSMsum->Scale(1.0 / hPhiEtaSMsum->GetYaxis()->GetBinWidth(1));
        TH1D* hEta = hPhiEtaSMsum->ProjectionY(Form("hEta_%d_%d%s", minRange, maxRange, suffix.Data()));
        hEta->SetTitle("#Delta#eta");
        hPhiEtaSMsum->Rebin2D(1, 1); //her forandret vi fra 2,2 til 1,1

        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaMsum->Rebin2D(1, 1);

        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaSsum->Rebin2D(1, 1);

        
        hPhiEtaSMsum->SetName(Form("dphideta_SM_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiEtaSsum->SetName(Form("dphideta_SE_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiEtaMsum->SetName(Form("dphideta_ME_%d_%d%s", minRange, maxRange, suffix.Data()));

        // Set axis titles and labels for hPhiEtaSMsum
        hPhiEtaSMsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaSMsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaSMsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaSMsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaSMsum->SetTitle(Form("Correlation Function %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaSMsum->SetTitle(Form("Correlation Function %d< Centrality #leq%d", minRange, maxRange));


        // Set axis titles and labels for hPhiEtaSsum
        hPhiEtaSsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaSsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaSsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaSsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaSsum->SetTitle(Form("Same Event %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaSsum->SetTitle(Form("Same Event %d< Centrality #leq%d", minRange, maxRange));

        // Set axis titles and labels for hPhiEtaMsum
        hPhiEtaMsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaMsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaMsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaMsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaMsum->SetTitle(Form("Mixed event %d< N_{ch} #leq%d", minRange, maxRange));
        else
            hPhiEtaMsum->SetTitle(Form("Mixed event %d< Centrality #leq%d", minRange, maxRange));

        // Draw histograms
        TCanvas* c1 = new TCanvas(Form("dPhidEta %s", suffix.Data()), Form("dPhidEta %s", suffix.Data()), 1200, 800);
        c1->Divide(2, 2);

        c1->cd(1);
        TH2D* hPhiEtaSMsum_draw = (TH2D*)hPhiEtaSMsum->Clone("hPhiEtaSMsum_draw");
        hPhiEtaSMsum_draw->Draw("surf1");
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s", collisionSystemName.c_str()));
        latex.DrawLatex(0.1, 0.85, Form("#eta^{TPC} #in [%0.1f, %0.1f], #eta^{FIT} #in [%0.1f, %0.1f]", etaMin, etaMax, ft0Min, ft0Max));


        c1->cd(3);
        hPhiEtaSsum->Draw("surf1");
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s", collisionSystemName.c_str()));
        latex.DrawLatex(0.1, 0.85, Form("#eta^{TPC} #in [%0.1f, %0.1f], #eta^{FIT} #in [%0.1f, %0.1f]", etaMin, etaMax, ft0Min, ft0Max));


        c1->cd(4);
        hPhiEtaMsum->Draw("surf1");
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s", collisionSystemName.c_str()));
        latex.DrawLatex(0.1, 0.85, Form("#eta^{TPC} #in [%0.1f, %0.1f], #eta^{FIT} #in [%0.1f, %0.1f]", etaMin, etaMax, ft0Min, ft0Max));

        // Saving to file
        fout->cd();
        c1->Write();
        // write canvas to file
        // c1->SaveAs(Form("./ProcessOutput/Mixed_%s_%s_%i_%i.png", fileNameSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange)));
        // delete c1;

        // Write histograms
        hPhiEtaSMsum->Write();
        hPhiEtaSsum->Write();
        hPhiEtaMsum->Write();
        hEta->Write();

        // Eta Gap processing
        // Eta Gap processing
        TH1D* hPhiSameOverMixed_pos = hPhiEtaSMsum->ProjectionX("hPhiSameOverMixed");
        hPhiSameOverMixed_pos->SetName(Form("hPhiSameOverMixed_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiSameOverMixed_pos->SetTitle(Form("hPhiSameOverMixed_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiSameOverMixed_pos->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiSameOverMixed_pos->Write();

        delete hPhiEtaSMsum;
        delete hPhiEtaSMsum_draw;
        delete hPhiEtaSsum;
        delete hPhiEtaMsum;
        delete hPhiSameOverMixed_pos;
        delete c1;
        
    }

    // clean up
    fout->Close();
    delete fout;
    file->Close();
    delete file;
    delete sparSig;
    delete sparMix;
    // End of process
    std::cout << "Output file: " << Form("./ProcessOutput/EtaDiff/Mixed_%s%s_%s_%i_%i_Eta_%0.1f_%0.1f_%s.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), etaMin, etaMax, DihadronCorrTypeName[corrType].c_str()) << std::endl;
    std::cout << "Processing completed for all samples." << std::endl;
#else
    std::cerr << "Error: CorrelationContainer support is disabled. Cannot process non-ring EtaDiff data for " << fileNameSuffix << std::endl;
    file->Close();
    delete file;
    return;
#endif
}

