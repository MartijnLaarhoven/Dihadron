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
#include "TSystem.h"
#include <TProfile.h>
#include <TRandom3.h>
#include "TMath.h"
#include "THnSparse.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "./include/BasicForDihadron.h"

struct InputUnit {
    std::string fileNameSuffix;
    Bool_t isNch;
    Bool_t isEtadiff;
    Int_t minRange;
    Int_t maxRange;
    Bool_t isMc;

    InputUnit(std::string _fileNameSuffix, Bool_t _isNch, Bool_t _isEtadiff, Int_t _minRange, Int_t _maxRange, Bool_t _isMc=false) :
        fileNameSuffix(_fileNameSuffix), isNch(_isNch), isEtadiff(_isEtadiff), minRange(_minRange), maxRange(_maxRange), isMc(_isMc) {}
};

void printAxesInfo(THnSparseF* sparseHist);
TH1D* ProcessEtaGap(TH2D* hPhiEtaSMsum, Int_t minRange, Int_t maxRange, Int_t sample);
void Read_dPhidEta_givenRange(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Bool_t isMc);
void Read_dPhidEta_givenRange_EtaDiff(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaTrigMin, Double_t etaTrigMax, Double_t etaAssoMin, Double_t etaAssoMax, Bool_t isMc);
void Read_dPhidEta_givenRange_EtaDiff_Reversed(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaTrigMin, Double_t etaTrigMax, Bool_t isMc);

// global variables
std::string collisionSystemName = "peripheral NeNe";
std::string additionalSuffix = "";

// Helper function to set collision system name based on filename
std::string GetCollisionSystemName(const std::string& filename) {
    if (filename.find("ae") != std::string::npos) {
        return "O-O";
    } else if (filename.find("af") != std::string::npos) {
        return "Ne-Ne";
    }
    return "Unknown";
}

void Process_dPhidEta() {
    // 不显示窗口
    gROOT->SetBatch(kTRUE);

    std::vector<InputUnit> inputList;
    additionalSuffix = "";

    // Set collision system name based on filename
    std::string inputFileName = "LHC25ae_pass2_604826";
    collisionSystemName = GetCollisionSystemName(inputFileName);

    inputList.push_back(InputUnit(inputFileName, kCent, kEtaDiffOff, 0, 20));
    inputList.push_back(InputUnit(inputFileName, kCent, kEtaDiffOn, 0, 20));

    inputList.push_back(InputUnit(inputFileName, kCent, kEtaDiffOff, 80, 100));    
    inputList.push_back(InputUnit(inputFileName, kCent, kEtaDiffOn, 80, 100));    

    // Process standard configuration
    for (auto input : inputList) {
        if (input.isEtadiff) {
            for (int iEta = 0; iEta < etaBins.size() - 1; iEta++) {
                double etaMin = etaBins[iEta];
                double etaMax = etaBins[iEta + 1];
                // Negative trigger eta correlated with positive associated eta
                // Split [0.7, 0.8] into two 0.05-width bins for finer dEta resolution
                std::vector<std::pair<double,double>> assoEtaBins = {{0.7, 0.75}, {0.75, 0.8}};
                for (auto& assoEta : assoEtaBins) {
                    Read_dPhidEta_givenRange_EtaDiff(input.fileNameSuffix, input.isNch, input.minRange, input.maxRange, etaMin, etaMax, assoEta.first, assoEta.second, input.isMc);
                }
            }
        }
        else {
            Read_dPhidEta_givenRange(input.fileNameSuffix, input.isNch, input.minRange, input.maxRange, input.isMc);
        }
    }

    // Process reversed eta configuration with new datasets (same file naming, no extra suffix)
    std::cout << "\n========== PROCESSING REVERSED ETA CONFIGURATION (O-O) ==========" << std::endl;
    std::string reversedFileName = "LHC25ae_pass2_604830";
    collisionSystemName = GetCollisionSystemName(reversedFileName);
    additionalSuffix = "";

    std::vector<InputUnit> reversedList;
    reversedList.push_back(InputUnit(reversedFileName, kCent, kEtaDiffOff, 0, 20));
    reversedList.push_back(InputUnit(reversedFileName, kCent, kEtaDiffOff, 80, 100));
    reversedList.push_back(InputUnit(reversedFileName, kCent, kEtaDiffOn, 0, 20));
    reversedList.push_back(InputUnit(reversedFileName, kCent, kEtaDiffOn, 80, 100));

    // Reversed eta bins: trigger positive [0.0, 0.8] in 0.1 steps
    std::vector<float> reversedEtaBins = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

    for (auto input : reversedList) {
        if (input.isEtadiff) {
            for (int iEta = 0; iEta < reversedEtaBins.size() - 1; iEta++) {
                double etaMin = reversedEtaBins[iEta];
                double etaMax = reversedEtaBins[iEta + 1];
                // Positive trigger eta correlated with negative associated eta
                // Split [-0.8, -0.7] into two 0.05-width bins for finer dEta resolution
                std::vector<std::pair<double,double>> assoEtaBins = {{-0.8, -0.75}, {-0.75, -0.7}};
                for (auto& assoEta : assoEtaBins) {
                    Read_dPhidEta_givenRange_EtaDiff(input.fileNameSuffix, input.isNch, input.minRange, input.maxRange, etaMin, etaMax, assoEta.first, assoEta.second, input.isMc);
                }
            }
        } else {
            Read_dPhidEta_givenRange(input.fileNameSuffix, input.isNch, input.minRange, input.maxRange, input.isMc);
        }
    }

    // Process Ne-Ne collision system (same pattern)
    std::cout << "\n========== PROCESSING STANDARD CONFIGURATION (Ne-Ne) ==========" << std::endl;
    std::string neInputFileName = "LHC25af_pass2_611697";
    collisionSystemName = GetCollisionSystemName(neInputFileName);

    std::vector<InputUnit> neInputList;
    neInputList.push_back(InputUnit(neInputFileName, kCent, kEtaDiffOff, 0, 20));
    neInputList.push_back(InputUnit(neInputFileName, kCent, kEtaDiffOn, 0, 20));

    neInputList.push_back(InputUnit(neInputFileName, kCent, kEtaDiffOff, 80, 100));    
    neInputList.push_back(InputUnit(neInputFileName, kCent, kEtaDiffOn, 80, 100));    

    for (auto input : neInputList) {
        if (input.isEtadiff) {
            for (int iEta = 0; iEta < etaBins.size() - 1; iEta++) {
                double etaMin = etaBins[iEta];
                double etaMax = etaBins[iEta + 1];
                // Negative trigger eta correlated with positive associated eta
                // Split [0.7, 0.8] into two 0.05-width bins for finer dEta resolution
                std::vector<std::pair<double,double>> assoEtaBins = {{0.7, 0.75}, {0.75, 0.8}};
                for (auto& assoEta : assoEtaBins) {
                    Read_dPhidEta_givenRange_EtaDiff(input.fileNameSuffix, input.isNch, input.minRange, input.maxRange, etaMin, etaMax, assoEta.first, assoEta.second, input.isMc);
                }
            }
        }
        else {
            Read_dPhidEta_givenRange(input.fileNameSuffix, input.isNch, input.minRange, input.maxRange, input.isMc);
        }
    }

    // Process reversed eta configuration for Ne-Ne
    std::cout << "\n========== PROCESSING REVERSED ETA CONFIGURATION (Ne-Ne) ==========" << std::endl;
    std::string neReversedFileName = "LHC25af_pass2_604820";
    collisionSystemName = GetCollisionSystemName(neReversedFileName);
    additionalSuffix = "";

    std::vector<InputUnit> neReversedList;
    neReversedList.push_back(InputUnit(neReversedFileName, kCent, kEtaDiffOff, 0, 20));
    neReversedList.push_back(InputUnit(neReversedFileName, kCent, kEtaDiffOff, 80, 100));
    neReversedList.push_back(InputUnit(neReversedFileName, kCent, kEtaDiffOn, 0, 20));
    neReversedList.push_back(InputUnit(neReversedFileName, kCent, kEtaDiffOn, 80, 100));

    for (auto input : neReversedList) {
        if (input.isEtadiff) {
            for (int iEta = 0; iEta < reversedEtaBins.size() - 1; iEta++) {
                double etaMin = reversedEtaBins[iEta];
                double etaMax = reversedEtaBins[iEta + 1];
                // Positive trigger eta correlated with negative associated eta
                // Split [-0.8, -0.7] into two 0.05-width bins for finer dEta resolution
                std::vector<std::pair<double,double>> assoEtaBins = {{-0.8, -0.75}, {-0.75, -0.7}};
                for (auto& assoEta : assoEtaBins) {
                    Read_dPhidEta_givenRange_EtaDiff(input.fileNameSuffix, input.isNch, input.minRange, input.maxRange, etaMin, etaMax, assoEta.first, assoEta.second, input.isMc);
                }
            }
        } else {
            Read_dPhidEta_givenRange(input.fileNameSuffix, input.isNch, input.minRange, input.maxRange, input.isMc);
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

// Helper function for Eta Gap processing
TH1D* ProcessEtaGap(TH2D* hPhiEtaSMsum, Int_t minRange, Int_t maxRange, Int_t sample) {
    double dEta = 1.0;
    // when validate with cumulant method, use dEta = 0.
    // double dEta = 0.0;
    double etaMin = -1.6;
    double etaMax = 1.6;

    Int_t minEtaBin = hPhiEtaSMsum->GetYaxis()->FindBin(dEta + 0.0001);
    Int_t maxEtaBin = hPhiEtaSMsum->GetYaxis()->FindBin(etaMax - 0.0001);
    TH1D* pos = hPhiEtaSMsum->ProjectionX("hPhiSameOverMixed_pos", minEtaBin, maxEtaBin);

    minEtaBin = hPhiEtaSMsum->GetYaxis()->FindBin(etaMin + 0.0001);
    maxEtaBin = hPhiEtaSMsum->GetYaxis()->FindBin(-dEta - 0.0001);
    TH1D* neg = hPhiEtaSMsum->ProjectionX("hPhiSameOverMixed_neg", minEtaBin, maxEtaBin);

    pos->Add(neg);
    pos->Scale(0.5);
    // Set names with sample suffix
    TString suffix = (sample == -1) ? "" : Form("_%d", sample);
    pos->SetName(Form("hPhiSameOverMixed_%d_%d%s", minRange, maxRange, suffix.Data()));
    pos->SetTitle(Form("hPhiSameOverMixed_%d_%d%s", minRange, maxRange, suffix.Data()));
    pos->GetXaxis()->SetTitle("#Delta#varphi");

    return pos;
}

void Read_dPhidEta_givenRange(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Bool_t isMc=false) {
    TFile *file = TFile::Open(Form("../AnalysisResultsROOTFiles/dihadronCorr/AnalysisResults_%s.root", fileNameSuffix.c_str()), "READ");
    if (!file || file->IsZombie()) {
        std::cout << "Error: Cannot open file " << fileNameSuffix << std::endl;
        return;
    }

    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    // Debug: Print what's being searched for and list directory contents
    std::cout << "Looking for: eta-dihadron_" << additionalSuffix << splitName << "_" << minRange << "_" << maxRange << std::endl;
    std::cout << "File contents:" << std::endl;
    file->ls();
    std::cout << "========================" << std::endl;

    // check if MCTrue folder is available
    if (isMc) {
        if (!file->Get(Form("eta-dihadron_%s%s_%d_%d/MCTrue", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange))) {
            std::cerr << "Error: MCTrue folder not found for " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
    }
    else {
        if (file->Get(Form("eta-dihadron_%s%s_%d_%d/MCTrue", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange))) {
            std::cerr << "Caution! you are using Reco or Data, but MCTrue folder is found for " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
    }


    CorrelationContainer *same = (CorrelationContainer*)file->Get(Form("eta-dihadron_%s%s_%d_%d/sameEvent", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange));
    CorrelationContainer *mixed = (CorrelationContainer*)file->Get(Form("eta-dihadron_%s%s_%d_%d/mixedEvent", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange));
    THnSparseD *trig = nullptr;
    if (!isMc) {
        trig = (THnSparseD*)file->Get(Form("eta-dihadron_%s%s_%d_%d/Trig_hist", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange));
    } else {
        trig = (THnSparseD*)file->Get(Form("eta-dihadron_%s%s_%d_%d/MCTrue/MCTrig_hist", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange));
    }
    

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

    sparSig->SetName(Form("sameEvent_%i_%i", minRange, maxRange));
    sparMix->SetName(Form("mixedEvent_%i_%i", minRange, maxRange));
    trig->SetName(Form("Trig_hist_%i_%i", minRange, maxRange));

    // Debug: Print sparse structure and entries
    std::cout << "Same event sparse entries: " << sparSig->GetEntries() << std::endl;
    std::cout << "Mixed event sparse entries: " << sparMix->GetEntries() << std::endl;
    std::cout << "Trigger sparse entries: " << trig->GetEntries() << std::endl;
    std::cout << "Number of dimensions in same sparse: " << sparSig->GetNdimensions() << std::endl;
    std::cout << "Number of dimensions in mixed sparse: " << sparMix->GetNdimensions() << std::endl;
    std::cout << "Number of dimensions in trig sparse: " << trig->GetNdimensions() << std::endl;
    
    // Print axis info
    for (int i = 0; i < sparSig->GetNdimensions(); i++) {
        std::cout << "Axis " << i << ": " << sparSig->GetAxis(i)->GetTitle() << " (bins: " << sparSig->GetAxis(i)->GetNbins() << ")" << std::endl;
    }

    // Common axis settings for all samples
    // NOTE: For eta-differential analysis, axes 2 and 3 are eta (not pT)
    // So we should NOT set pT ranges - they don't exist in this sparse!
    // sparSig->GetAxis(corrAxis_kPt_TPC_trig)->SetRangeUser(minPt+0.001, maxPt-0.001);
    // sparMix->GetAxis(corrAxis_kPt_TPC_trig)->SetRangeUser(minPt+0.001, maxPt-0.001);
    sparSig->GetAxis(corrAxis_kdEtaTPCTPC)->SetRangeUser(-1.6, 1.6);
    sparMix->GetAxis(corrAxis_kdEtaTPCTPC)->SetRangeUser(-1.6, 1.6);
    // sparSig->GetAxis(corrAxis_kPt_TPC_asso)->SetRangeUser(minPt+0.001, maxPt-0.001);
    // sparMix->GetAxis(corrAxis_kPt_TPC_asso)->SetRangeUser(minPt+0.001, maxPt-0.001);
    // trig->GetAxis(trigAxis_pT)->SetRangeUser(minPt+0.001, maxPt-0.001);
    
    // Debug: print dEta axis binning information
    TAxis* detaAxis = sparSig->GetAxis(corrAxis_kdEtaTPCTPC);
    std::cout << "dEta axis info - Nbins: " << detaAxis->GetNbins() << ", Xmin: " << detaAxis->GetXmin() << ", Xmax: " << detaAxis->GetXmax() << ", BinWidth: " << detaAxis->GetBinWidth(1) << std::endl;
    
    std::cout << "After range setting - Same entries: " << sparSig->GetEntries() 
              << ", Mixed entries: " << sparMix->GetEntries() << std::endl;

    // Create output file
    TFile* fout = TFile::Open(Form("./ProcessOutput/Mixed_%s%s_%s_%i_%i.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange)), "RECREATE");

    // Check if sample axis exists
    Bool_t hasSampleAxis = (sparSig->GetNdimensions() > corrAxis_kSample);
    Int_t maxSampleToProcess = hasSampleAxis ? maxSample : 0;
    std::cout << "Has sample axis: " << hasSampleAxis << ", will process " << (hasSampleAxis ? maxSampleToProcess + 1 : 1) << " sample(s)" << std::endl;

    // Process all samples: -1 (all), 0 to maxSample-1
    for (Int_t sample = -1; sample < maxSampleToProcess; ++sample) {
        // Set sample range if applicable
        if (sample >= 0 && hasSampleAxis) {
            sparSig->GetAxis(corrAxis_kSample)->SetRangeUser(sample, sample+0.001);
            sparMix->GetAxis(corrAxis_kSample)->SetRangeUser(sample, sample+0.001);
            if (trig->GetNdimensions() > trigAxis_sample) {
                trig->GetAxis(trigAxis_sample)->SetRangeUser(sample, sample+0.001);
            }
        } else if (hasSampleAxis) {
            // Reset to full range for inclusive sample
            double xmin = sparSig->GetAxis(corrAxis_kSample)->GetXmin();
            double xmax = sparSig->GetAxis(corrAxis_kSample)->GetXmax();
            sparSig->GetAxis(corrAxis_kSample)->SetRangeUser(xmin, xmax);
            xmin = sparMix->GetAxis(corrAxis_kSample)->GetXmin();
            xmax = sparMix->GetAxis(corrAxis_kSample)->GetXmax();
            sparMix->GetAxis(corrAxis_kSample)->SetRangeUser(xmin, xmax);
            if (trig->GetNdimensions() > trigAxis_sample) {
                xmin = trig->GetAxis(trigAxis_sample)->GetXmin();
                xmax = trig->GetAxis(trigAxis_sample)->GetXmax();
                trig->GetAxis(trigAxis_sample)->SetRangeUser(xmin, xmax);
            }
        }

        // Initialize variables for this sample
        TH2D* hPhiEtaSMsum = nullptr;
        TH2D* hPhiEtaSsum = nullptr;
        TH2D* hPhiEtaMsum = nullptr;
        Double_t nTriggersS = 0.;
        Int_t nz = sparSig->GetAxis(corrAxis_kVz)->GetNbins();
        // Set names with sample suffix
        TString suffix = (sample == -1) ? "" : Form("_%d", sample);

        std::cout << "==> Processing sample " << sample << ", nz = " << nz << std::endl;

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
            Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(-TMath::Pi()/2 + 0.0001);
            Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(3*TMath::Pi()/2 - 0.0001);
            Int_t binEta1 = hPhiEtaM->GetYaxis()->FindBin(-0.0001);
            Int_t binEta2 = hPhiEtaM->GetYaxis()->FindBin(0.0001);
            Int_t nNormBins = (binEta2 - binEta1 + 1) * (binPhi2 - binPhi1 + 1);
            norm = hPhiEtaM->Integral(binPhi1, binPhi2, binEta1, binEta2) / nNormBins;
            if (norm <= 0) {
                // Fallback: use full Δη range if Δη≈0 window is empty (e.g., eta-diff bins)
                Int_t binEtaAll1 = 1;
                Int_t binEtaAll2 = hPhiEtaM->GetYaxis()->GetNbins();
                Int_t binPhiAll1 = 1;
                Int_t binPhiAll2 = hPhiEtaM->GetXaxis()->GetNbins();
                Int_t nAllBins = (binEtaAll2 - binEtaAll1 + 1) * (binPhiAll2 - binPhiAll1 + 1);
                norm = hPhiEtaM->Integral(binPhiAll1, binPhiAll2, binEtaAll1, binEtaAll2) / nAllBins;
            }

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
        // No rebinning - keep original fine binning from sparse histogram
        // hPhiEtaSMsum->Rebin2D(1, 1);

        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetYaxis()->GetBinWidth(1));
        // hPhiEtaMsum->Rebin2D(1, 1);

        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetYaxis()->GetBinWidth(1));
        // hPhiEtaSsum->Rebin2D(1, 1);

        
        hPhiEtaSMsum->SetName(Form("dphideta_SM_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiEtaSsum->SetName(Form("dphideta_SE_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiEtaMsum->SetName(Form("dphideta_ME_%d_%d%s", minRange, maxRange, suffix.Data()));

        // set the ranges
        hPhiEtaSMsum->GetYaxis()->SetRangeUser(-1.6, 1.6);
        hPhiEtaSsum->GetYaxis()->SetRangeUser(-1.6, 1.6);
        hPhiEtaMsum->GetYaxis()->SetRangeUser(-1.6, 1.6);

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
        hPhiEtaSMsum_draw->GetYaxis()->SetRangeUser(-1.5, 1.5);
        hPhiEtaSMsum_draw->Draw("surf1");
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s", collisionSystemName.c_str()));
        latex.DrawLatex(0.1, 0.85, Form("p_{T}^{trig} #in [%0.1f, %0.1f] GeV/c, p_{T}^{asso} #in [%0.1f, %0.1f] GeV/c", minPt, maxPt, minPt, maxPt));


        c1->cd(3);
        TH2D* hPhiEtaSsum_draw = (TH2D*)hPhiEtaSsum->Clone("hPhiEtaSsum_draw");
        hPhiEtaSsum_draw->GetYaxis()->SetRangeUser(-1.5, 1.5);
        hPhiEtaSsum_draw->Draw("surf1");
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s", collisionSystemName.c_str()));
        latex.DrawLatex(0.1, 0.85, Form("p_{T}^{trig} #in [%0.1f, %0.1f] GeV/c, p_{T}^{asso} #in [%0.1f, %0.1f] GeV/c", minPt, maxPt, minPt, maxPt));


        c1->cd(4);
        TH2D* hPhiEtaMsum_draw = (TH2D*)hPhiEtaMsum->Clone("hPhiEtaMsum_draw");
        hPhiEtaMsum_draw->GetYaxis()->SetRangeUser(-1.5, 1.5);
        hPhiEtaMsum_draw->Draw("surf1");
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s", collisionSystemName.c_str()));
        latex.DrawLatex(0.1, 0.85, Form("p_{T}^{trig} #in [%0.1f, %0.1f] GeV/c, p_{T}^{asso} #in [%0.1f, %0.1f] GeV/c", minPt, maxPt, minPt, maxPt));

        // Saving to file
        fout->cd();
        c1->Write();
        // Write histograms
        hPhiEtaSMsum->Write();
        hPhiEtaSsum->Write();
        hPhiEtaMsum->Write();
        hEta->Write();

        // Eta Gap processing
        TH1D* hPhiSameOverMixed_pos = ProcessEtaGap(hPhiEtaSMsum, minRange, maxRange, sample);
        hPhiSameOverMixed_pos->Write();

        delete hPhiEtaSMsum;
        delete hPhiEtaSMsum_draw;
        delete hPhiEtaSsum_draw;
        delete hPhiEtaSsum;
        delete hPhiEtaMsum_draw;
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
    std::cout << "Output file: " << Form("./ProcessOutput/Mixed_%s%s_%s_%i_%i.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange)) << std::endl; 
    std::cout << "Processing completed for all samples." << std::endl;
}

void Read_dPhidEta_givenRange_EtaDiff(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax, Double_t etaAssoMin, Double_t etaAssoMax, Bool_t isMc=false) {
    TFile *file = TFile::Open(Form("../AnalysisResultsROOTFiles/dihadronCorr/AnalysisResults_%s.root", fileNameSuffix.c_str()), "READ");
    if (!file || file->IsZombie()) {
        std::cout << "Error: Cannot open file " << fileNameSuffix << std::endl;
        return;
    }

    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    // Debug: Print what's being searched for
    std::cout << "Looking for: eta-dihadron_" << additionalSuffix << splitName << "_" << minRange << "_" << maxRange << " (Eta range: " << etaMin << " to " << etaMax << ")" << std::endl;

    // check if MCTrue folder is available
    if (isMc) {
        if (!file->Get(Form("eta-dihadron_%s%s_%d_%d/MCTrue", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange))) {
            std::cerr << "Error: MCTrue folder not found for " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
    }
    else {
        if (file->Get(Form("eta-dihadron_%s%s_%d_%d/MCTrue", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange))) {
            std::cerr << "Caution! you are using Reco or Data, but MCTrue folder is found for " << fileNameSuffix << std::endl;
            file->Close();
            delete file;
            return;
        }
    }


    CorrelationContainer *same = (CorrelationContainer*)file->Get(Form("eta-dihadron_%s%s_%d_%d/sameEvent", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange));
    CorrelationContainer *mixed = (CorrelationContainer*)file->Get(Form("eta-dihadron_%s%s_%d_%d/mixedEvent", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange));
    THnSparseD *trig = nullptr;
    if (!isMc) {
        trig = (THnSparseD*)file->Get(Form("eta-dihadron_%s%s_%d_%d/Trig_hist", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange));
    } else {
        trig = (THnSparseD*)file->Get(Form("eta-dihadron_%s%s_%d_%d/MCTrue/MCTrig_hist", additionalSuffix.c_str(), splitName.c_str(), minRange, maxRange));
    }
    

    if (!same || !mixed || !trig) {
        std::cerr << "Error getting histograms for " << fileNameSuffix << " with " << splitName << " and range [" << minRange << ", " << maxRange << "]" << std::endl;
        file->Close();
        delete file;
        return;
    }

    // Try to get the CF step with data - use kCFStepAll (step 6) for eta-diff
    Int_t step = 6;  // kCFStepAll = 6, which typically has all entries
    THnSparseF *sparSig = nullptr;
    THnSparseF *sparMix = nullptr;
    
    // First try step 6 (kCFStepAll)
    std::cout << "DEBUG: Trying CF step " << step << std::endl;
    sparSig = (THnSparseF*)same->getPairHist()->getTHn(step);
    sparMix = (THnSparseF*)mixed->getPairHist()->getTHn(step);
    
    if (!sparSig || !sparMix) {
        std::cerr << "ERROR: Could not get sparse histograms from step " << step << std::endl;
        file->Close();
        delete file;
        return;
    }
    
    // Check actual bin content, not GetEntries() which can be misleading
    Long64_t nSameNonZero = 0;
    Long64_t nMixedNonZero = 0;
    for (Long64_t i = 0; i < sparSig->GetNbins(); ++i) {
        if (sparSig->GetBinContent(i) > 0) nSameNonZero++;
        if (sparMix->GetBinContent(i) > 0) nMixedNonZero++;
        if (nSameNonZero > 0 && nMixedNonZero > 0) break; // Early exit once we know there's data
    }
    
    std::cout << "DEBUG: Step " << step << " - Same non-zero bins: " << nSameNonZero 
              << ", Mixed non-zero bins: " << nMixedNonZero << std::endl;

    sparSig->SetName(Form("sameEvent_%i_%i", minRange, maxRange));
    sparMix->SetName(Form("mixedEvent_%i_%i", minRange, maxRange));
    trig->SetName(Form("Trig_hist_%i_%i", minRange, maxRange));

    // Common axis settings for all samples
    // In eta-differential analysis: axis 2 = trigger eta, axis 3 = associated eta
    // Select the trigger eta bin using bin indices to avoid edge precision issues
    auto* axisTrig = sparSig->GetAxis(corrAxis_kPt_TPC_trig);
    auto* axisTrigMix = sparMix->GetAxis(corrAxis_kPt_TPC_trig);
    auto* axisTrigSparse = trig->GetAxis(trigAxis_pT);
    
    // Debug: print axis structure
    std::cout << "DEBUG: Selecting eta bin [" << etaMin << ", " << etaMax << "]" << std::endl;
    std::cout << "DEBUG: Axis 2 range: [" << axisTrig->GetXmin() << ", " << axisTrig->GetXmax() 
              << "], nbins: " << axisTrig->GetNbins() << std::endl;
    for (Int_t ib = 1; ib <= axisTrig->GetNbins(); ++ib) {
        std::cout << "  Bin " << ib << ": [" << axisTrig->GetBinLowEdge(ib) 
                  << ", " << axisTrig->GetBinUpEdge(ib) << "]" << std::endl;
    }
    
    Int_t binMin = axisTrig->FindBin(etaMin + 1e-6);
    Int_t binMax = axisTrig->FindBin(etaMax - 1e-6);
    std::cout << "DEBUG: FindBin(" << etaMin + 1e-6 << ") = " << binMin << std::endl;
    std::cout << "DEBUG: FindBin(" << etaMax - 1e-6 << ") = " << binMax << std::endl;
    
    if (binMin > binMax) std::swap(binMin, binMax);
    axisTrig->SetRange(binMin, binMax);
    axisTrigMix->SetRange(binMin, binMax);
    Int_t binMinTrig = axisTrigSparse->FindBin(etaMin + 1e-6);
    Int_t binMaxTrig = axisTrigSparse->FindBin(etaMax - 1e-6);
    if (binMinTrig > binMaxTrig) std::swap(binMinTrig, binMaxTrig);
    axisTrigSparse->SetRange(binMinTrig, binMaxTrig);
    
    std::cout << "DEBUG: After SetRange - Same entries: " << sparSig->GetEntries() 
              << ", Mixed entries: " << sparMix->GetEntries() 
              << ", Trig entries: " << trig->GetEntries() << std::endl;
    
    // Check if we have enough statistics
    Long64_t nSameEntries = 0;
    Long64_t nMixedEntries = 0;
    for (Long64_t i = 0; i < sparSig->GetNbins(); ++i) {
        if (sparSig->GetBinContent(i) > 0) nSameEntries++;
        if (sparMix->GetBinContent(i) > 0) nMixedEntries++;
    }
    std::cout << "DEBUG: Non-zero bins - Same: " << nSameEntries << ", Mixed: " << nMixedEntries << std::endl;
    
    if (nSameEntries == 0 || nMixedEntries == 0) {
        std::cerr << "WARNING: Eta bin [" << etaMin << ", " << etaMax << "] has no statistics. Skipping." << std::endl;
        file->Close();
        delete file;
        return;
    }
    
    // Calculate expected dEta range based on trigger and associated eta
    // Associated eta ranges are passed as parameters to handle different configurations
    // dEta = eta_asso - eta_trig, so:
    // dEtaMin = etaAssoMin - etaMax (most negative)
    // dEtaMax = etaAssoMax - etaMin (most positive)
    Double_t expectedDEtaMin = etaAssoMin - etaMax;
    Double_t expectedDEtaMax = etaAssoMax - etaMin;
    std::cout << "Expected dEta range for trigger eta [" << etaMin << ", " << etaMax 
              << "] and asso eta [" << etaAssoMin << ", " << etaAssoMax << "]: [" 
              << expectedDEtaMin << ", " << expectedDEtaMax << "]" << std::endl;
    
    // Do NOT restrict DeltaEta range here - let the projection get all data
    // We'll restrict the display range on the final histograms instead

    // Create output file
    TFile* fout = TFile::Open(Form("./ProcessOutput/EtaDiff/Mixed_%s%s_%s_%i_%i_Eta_%0.1f_%0.1f.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), etaMin, etaMax), "RECREATE");

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
            Int_t binPhi1 = hPhiEtaM->GetXaxis()->FindBin(-TMath::Pi()/2 + 0.0001);
            Int_t binPhi2 = hPhiEtaM->GetXaxis()->FindBin(3*TMath::Pi()/2 - 0.0001);
            Int_t binEta1 = hPhiEtaM->GetYaxis()->FindBin(-0.0001);
            Int_t binEta2 = hPhiEtaM->GetYaxis()->FindBin(0.0001);
            Int_t nNormBins = (binEta2 - binEta1 + 1) * (binPhi2 - binPhi1 + 1);
            norm = hPhiEtaM->Integral(binPhi1, binPhi2, binEta1, binEta2) / nNormBins;
            if (norm <= 0) {
                // Fallback: use full Δη range if Δη≈0 window is empty (e.g., eta-diff bins)
                Int_t binEtaAll1 = 1;
                Int_t binEtaAll2 = hPhiEtaM->GetYaxis()->GetNbins();
                Int_t binPhiAll1 = 1;
                Int_t binPhiAll2 = hPhiEtaM->GetXaxis()->GetNbins();
                Int_t nAllBins = (binEtaAll2 - binEtaAll1 + 1) * (binPhiAll2 - binPhiAll1 + 1);
                norm = hPhiEtaM->Integral(binPhiAll1, binPhiAll2, binEtaAll1, binEtaAll2) / nAllBins;
            }

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
        hPhiEtaSMsum->Rebin2D(2, 2);

        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaMsum->Scale(1.0 / hPhiEtaMsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaMsum->Rebin2D(1, 1);

        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetXaxis()->GetBinWidth(1));
        hPhiEtaSsum->Scale(1.0 / hPhiEtaSsum->GetYaxis()->GetBinWidth(1));
        hPhiEtaSsum->Rebin2D(1, 1);

        
        hPhiEtaSMsum->SetName(Form("dphideta_SM_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiEtaSsum->SetName(Form("dphideta_SE_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiEtaMsum->SetName(Form("dphideta_ME_%d_%d%s", minRange, maxRange, suffix.Data()));

        // Generate the 1D projections BEFORE setting axis ranges (so they include all dEta data)
        TH1D* hPhiSameOverMixed_temp = (TH1D*)hPhiEtaSMsum->ProjectionX(Form("hPhiSameOverMixed_%d_%d%s_temp", minRange, maxRange, suffix.Data()));
        
        // Rebin: 2 in phi direction, 1 in dEta direction (keeps original dEta binning)
        hPhiEtaSMsum->Rebin2D(2, 1);
        // Don't rebin mixed and single - keep them as is
        // hPhiEtaMsum->Rebin2D(1, 1);
        // hPhiEtaSsum->Rebin2D(1, 1);
        Int_t firstBinWithData = -1;
        Int_t lastBinWithData = -1;
        Int_t nYbins = hPhiEtaSMsum->GetYaxis()->GetNbins();
        
        // Scan Y-axis (dEta) to find first and last bins with non-zero content
        for (Int_t iy = 1; iy <= nYbins; ++iy) {
            // Sum all X bins (dPhi) for this Y bin
            Double_t sumContent = 0;
            for (Int_t ix = 1; ix <= hPhiEtaSMsum->GetXaxis()->GetNbins(); ++ix) {
                sumContent += hPhiEtaSMsum->GetBinContent(ix, iy);
            }
            if (sumContent > 0) {
                if (firstBinWithData == -1) firstBinWithData = iy;
                lastBinWithData = iy;
            }
        }
        
        // Calculate the actual dEta range from the bin edges
        Double_t actualDEtaMin, actualDEtaMax;
        if (firstBinWithData > 0 && lastBinWithData > 0) {
            actualDEtaMin = hPhiEtaSMsum->GetYaxis()->GetBinLowEdge(firstBinWithData);
            actualDEtaMax = hPhiEtaSMsum->GetYaxis()->GetBinUpEdge(lastBinWithData);
            std::cout << "Actual dEta range found: [" << actualDEtaMin << ", " << actualDEtaMax << "]" << std::endl;
        } else {
            // Fallback to calculated range if no data found
            actualDEtaMin = expectedDEtaMin;
            actualDEtaMax = expectedDEtaMax;
            std::cout << "No data found in Y-axis, using calculated range: [" << actualDEtaMin << ", " << actualDEtaMax << "]" << std::endl;
        }
        
        // DO NOT set Y-axis range on original histograms - keep them with full data for file writing
        // We'll only restrict the display range on the clones used for drawing

        // Set axis titles and labels for hPhiEtaSMsum
        hPhiEtaSMsum->GetXaxis()->SetTitle("#Delta#varphi");
        hPhiEtaSMsum->GetYaxis()->SetTitle("#Delta#eta");
        hPhiEtaSMsum->GetXaxis()->SetTitleSize(0.05);
        hPhiEtaSMsum->GetYaxis()->SetTitleSize(0.05);
        if (isNch)
            hPhiEtaSMsum->SetTitle(Form("Correlation Function %d< N_{ch} #leq%d, #eta_{trig}#in[%0.1f,%0.1f]", minRange, maxRange, etaMin, etaMax));
        else
            hPhiEtaSMsum->SetTitle(Form("Correlation Function %d< Centrality #leq%d, #eta_{trig}#in[%0.1f,%0.1f]", minRange, maxRange, etaMin, etaMax));

        // Draw histograms - with restricted dEta range for clarity
        TCanvas* c1 = new TCanvas(Form("dPhidEta %s", suffix.Data()), Form("dPhidEta %s", suffix.Data()), 1200, 800);
        c1->Divide(2, 2);

        c1->cd(1);
        TH2D* hPhiEtaSMsum_draw = (TH2D*)hPhiEtaSMsum->Clone("hPhiEtaSMsum_draw");
        hPhiEtaSMsum_draw->GetYaxis()->SetRangeUser(actualDEtaMin, actualDEtaMax);
        hPhiEtaSMsum_draw->Draw("surf1");
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s", collisionSystemName.c_str()));
        latex.DrawLatex(0.1, 0.85, Form("#eta_{trig}#in[%0.1f,%0.1f], #eta_{asso}#in[%0.1f,%0.1f]", etaMin, etaMax, etaAssoMin, etaAssoMax));
        latex.DrawLatex(0.1, 0.80, Form("p_{T}^{trig}#in[%0.1f,%0.1f] GeV/c, p_{T}^{asso}#in[%0.1f,%0.1f] GeV/c", pTMin, pTMax, minPt, maxPt));

        c1->cd(3);
        TH2D* hPhiEtaSsum_draw = (TH2D*)hPhiEtaSsum->Clone("hPhiEtaSsum_draw");
        hPhiEtaSsum_draw->GetYaxis()->SetRangeUser(actualDEtaMin, actualDEtaMax);
        hPhiEtaSsum_draw->Draw("surf1");
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s", collisionSystemName.c_str()));
        latex.DrawLatex(0.1, 0.85, Form("#eta_{trig}#in[%0.1f,%0.1f], #eta_{asso}#in[%0.1f,%0.1f]", etaMin, etaMax, etaAssoMin, etaAssoMax));
        latex.DrawLatex(0.1, 0.80, Form("p_{T}^{trig}#in[%0.1f,%0.1f] GeV/c, p_{T}^{asso}#in[%0.1f,%0.1f] GeV/c", pTMin, pTMax, minPt, maxPt));

        c1->cd(4);
        TH2D* hPhiEtaMsum_draw = (TH2D*)hPhiEtaMsum->Clone("hPhiEtaMsum_draw");
        hPhiEtaMsum_draw->GetYaxis()->SetRangeUser(actualDEtaMin, actualDEtaMax);
        hPhiEtaMsum_draw->Draw("surf1");
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.1, 0.90, Form("ALICE %s", collisionSystemName.c_str()));
        latex.DrawLatex(0.1, 0.85, Form("#eta_{trig}#in[%0.1f,%0.1f], #eta_{asso}#in[%0.1f,%0.1f]", etaMin, etaMax, etaAssoMin, etaAssoMax));
        latex.DrawLatex(0.1, 0.80, Form("p_{T}^{trig}#in[%0.1f,%0.1f] GeV/c, p_{T}^{asso}#in[%0.1f,%0.1f] GeV/c", pTMin, pTMax, minPt, maxPt));

        // Saving to file
        fout->cd();
        c1->Write();
        // Write histograms
        hPhiEtaSMsum->Write();
        hPhiEtaSsum->Write();
        hPhiEtaMsum->Write();
        hEta->Write();

        // For eta-differential analysis, project over full delta eta range (not eta gap)
        // Use the pre-generated projection (created before axis range was set)
        hPhiSameOverMixed_temp->SetName(Form("hPhiSameOverMixed_%d_%d%s", minRange, maxRange, suffix.Data()));
        hPhiSameOverMixed_temp->SetTitle(Form("Same/Mixed #Delta#phi (all #Delta#eta); #Delta#phi [rad]; S/M"));
        hPhiSameOverMixed_temp->Write();

        delete hPhiEtaSMsum;
        delete hPhiEtaSMsum_draw;
        delete hPhiEtaSsum_draw;
        delete hPhiEtaSsum;
        delete hPhiEtaMsum_draw;
        delete hPhiEtaMsum;
        delete hPhiSameOverMixed_temp;
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
    std::cout << "Output file: " << Form("./ProcessOutput/EtaDiff/Mixed_%s%s_%s_%i_%i_Eta_%0.1f_%0.1f.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), etaMin, etaMax) << std::endl;
    std::cout << "Processing completed for all samples." << std::endl;
}

// Reversed eta configuration: trigger positive [0, 0.8], associated negative [-0.8, -0.7]
void Read_dPhidEta_givenRange_EtaDiff_Reversed(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaTrigMin, Double_t etaTrigMax, Bool_t isMc) {
    TFile *file = TFile::Open(Form("../AnalysisResultsROOTFiles/dihadronCorr/AnalysisResults_%s.root", fileNameSuffix.c_str()), "READ");
    if (!file || file->IsZombie()) {
        std::cout << "Error: Cannot open file " << fileNameSuffix << std::endl;
        return;
    }

    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    // Get the 2D histograms directly (simplified approach without CorrelationContainer)
    TH2D *h2DSame = (TH2D*)file->Get(Form("eta-dihadron_%s_%d_%d/deltaEta_deltaPhi_same", splitName.c_str(), minRange, maxRange));
    TH2D *h2DMixed = (TH2D*)file->Get(Form("eta-dihadron_%s_%d_%d/deltaEta_deltaPhi_mixed", splitName.c_str(), minRange, maxRange));

    if (!h2DSame || !h2DMixed) {
        std::cerr << "Error: Could not find 2D histograms for reversed eta!" << std::endl;
        file->Close();
        delete file;
        return;
    }

    // REVERSED: Associated eta range [-0.8, -0.7] 
    double etaAssocMin = -0.8;
    double etaAssocMax = -0.7;

    std::cout << "\nProcessing reversed eta configuration:" << std::endl;
    std::cout << "  Trigger eta: [" << etaTrigMin << ", " << etaTrigMax << "]" << std::endl;
    std::cout << "  Associated eta: [" << etaAssocMin << ", " << etaAssocMax << "]" << std::endl;
    std::cout << "  Centrality: [" << minRange << ", " << maxRange << "%]" << std::endl;

    TAxis* axisEtaSame = h2DSame->GetXaxis();
    TAxis* axisPhiSame = h2DSame->GetYaxis();
    TAxis* axisEtaMixed = h2DMixed->GetXaxis();
    TAxis* axisPhiMixed = h2DMixed->GetYaxis();

    // Select the associated eta range
    Int_t binEtaMin = axisEtaSame->FindBin(etaAssocMin);
    Int_t binEtaMax = axisEtaSame->FindBin(etaAssocMax - 0.001);
    
    if (binEtaMin > binEtaMax) std::swap(binEtaMin, binEtaMax);

    // Project onto Delta Phi axis for the selected eta range
    TH1D *hPhiSame = (TH1D*)h2DSame->ProjectionY("hPhiSame", binEtaMin, binEtaMax);
    TH1D *hPhiMixed = (TH1D*)h2DMixed->ProjectionY("hPhiMixed", binEtaMin, binEtaMax);

    if (!hPhiSame || !hPhiMixed) {
        std::cerr << "Error: Could not project histograms!" << std::endl;
        file->Close();
        delete file;
        return;
    }

    // Normalize mixed event
    Int_t binPhi1 = hPhiMixed->GetXaxis()->FindBin(-TMath::Pi()/2 + 0.0001);
    Int_t binPhi2 = hPhiMixed->GetXaxis()->FindBin(3*TMath::Pi()/2 - 0.0001);
    Double_t norm = hPhiMixed->Integral(binPhi1, binPhi2) / (double)(binPhi2 - binPhi1 + 1);
    if (norm > 0) {
        hPhiMixed->Scale(1.0 / norm);
    }

    // Divide same by mixed
    hPhiSame->Divide(hPhiMixed);

    hPhiSame->SetName(Form("hPhiSameOverMixed_%d_%d", minRange, maxRange));
    hPhiSame->SetTitle(Form("Correlation (Reversed Eta) %d-%d%%", minRange, maxRange));

    // Save to file
    gSystem->mkdir("./ProcessOutput/EtaDiff", kTRUE);
    TFile* fout = TFile::Open(Form("./ProcessOutput/EtaDiff/Mixed_%s%s_%s_%i_%i_Eta_%0.1f_%0.1f.root", 
                                   fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), etaTrigMin, etaTrigMax), "RECREATE");
    hPhiSame->Write();

    fout->Close();
    delete fout;
    file->Close();
    delete file;

    std::cout << "Output file: " << Form("./ProcessOutput/EtaDiff/Mixed_%s%s_%s_%i_%i_Eta_%0.1f_%0.1f.root", 
                                         fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), etaTrigMin, etaTrigMax) << std::endl;
}