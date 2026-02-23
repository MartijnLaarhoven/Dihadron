/*
 * @Author: Martijn Laarhoven (martijn.laarhoven@cern.ch)
 * @Date: 2026-02-11
 * @Last Modified by: Martijn Laarhoven
 * @Last Modified time: 2026-02-11
 */
// Uncorrected v_n^Δ extraction using Fourier fit on raw S/M correlations
// No peripheral subtraction, no ZYAM baseline - just pure correlation Fourier decomposition
// Used as baseline to see impact of non-flow corrections

#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "./include/BasicForDihadron.h"

struct InputUnit {
    std::string fileNameSuffix;
    Bool_t isNch;
    Bool_t isEtadiff;
    Int_t minRange;
    Int_t maxRange;

    InputUnit(std::string _fileNameSuffix, Bool_t _isNch, Bool_t _isEtadiff, Int_t _minRange, Int_t _maxRange) :
        fileNameSuffix(_fileNameSuffix), isNch(_isNch), isEtadiff(_isEtadiff), minRange(_minRange), maxRange(_maxRange) {}
};

struct VnUnit {
    Double_t v2, v2_err;
    Double_t v3, v3_err;
    Double_t v4, v4_err;
    
    VnUnit() : v2(-1), v2_err(10), v3(-1), v3_err(10), v4(-1), v4_err(10) {}
};

void ProcessConfig_Uncorrected(Bool_t isNch, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName);
void ProcessConfig_Uncorrected_EtaDiff(Bool_t isNch, Bool_t isEtaDiff, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName, const std::vector<float>& etaBins);
TGraphErrors* FourierFit_Uncorrected_GetGraph(Bool_t isNch, InputUnit data, Double_t etaMin, Double_t etaMax);
void CombineEtaDiffV2Plots_Uncorrected(const std::string& negFile, const std::string& posFile, const std::string& outTag, const std::string& splitName);
void CombineSystemPairs_Uncorrected(const std::string& system1File, const std::string& system2File, const std::string& outFile, const std::string& graphName);

void Process_Uncorrected_FourierFit() {
    gROOT->SetBatch(kTRUE);

    // Define eta bins for both configurations
    std::vector<float> etaBinsNeg = {-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0};
    std::vector<float> etaBinsPos = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

    std::vector<InputUnit> inputList;

    // O-O collisions
    std::string inputFileName = "LHC25ae_pass2_604826";
    std::vector<InputUnit> dataListOO;
    dataListOO.push_back(InputUnit(inputFileName, kCent, kEtaDiffOff, 0, 20));
    
    InputUnit templOO(inputFileName, kCent, kEtaDiffOff, 0, 20);
    ProcessConfig_Uncorrected(kCent, templOO, dataListOO, inputFileName + "_Uncorrected");
    ProcessConfig_Uncorrected_EtaDiff(kCent, kEtaDiffOn, templOO, dataListOO, inputFileName + "_Uncorrected", etaBinsNeg);

    // O-O reversed
    std::string inputFileNameNew = "LHC25ae_pass2_604830";
    std::vector<InputUnit> dataListOONew;
    dataListOONew.push_back(InputUnit(inputFileNameNew, kCent, kEtaDiffOff, 0, 20));
    
    InputUnit templOONew(inputFileNameNew, kCent, kEtaDiffOff, 0, 20);
    ProcessConfig_Uncorrected(kCent, templOONew, dataListOONew, inputFileNameNew + "_Uncorrected");
    ProcessConfig_Uncorrected_EtaDiff(kCent, kEtaDiffOn, templOONew, dataListOONew, inputFileNameNew + "_Uncorrected", etaBinsPos);

    // Ne-Ne collisions
    std::string inputFileNameNeStd = "LHC25af_pass2_611697";
    std::vector<InputUnit> dataListNeStd;
    dataListNeStd.push_back(InputUnit(inputFileNameNeStd, kCent, kEtaDiffOff, 0, 20));
    
    InputUnit templNeStd(inputFileNameNeStd, kCent, kEtaDiffOff, 0, 20);
    ProcessConfig_Uncorrected(kCent, templNeStd, dataListNeStd, inputFileNameNeStd + "_Uncorrected");
    ProcessConfig_Uncorrected_EtaDiff(kCent, kEtaDiffOn, templNeStd, dataListNeStd, inputFileNameNeStd + "_Uncorrected", etaBinsNeg);

    // Ne-Ne reversed
    std::string inputFileNameNeRev = "LHC25af_pass2_604820";
    std::vector<InputUnit> dataListNeRev;
    dataListNeRev.push_back(InputUnit(inputFileNameNeRev, kCent, kEtaDiffOff, 0, 20));
    
    InputUnit templNeRev(inputFileNameNeRev, kCent, kEtaDiffOff, 0, 20);
    ProcessConfig_Uncorrected(kCent, templNeRev, dataListNeRev, inputFileNameNeRev + "_Uncorrected");
    ProcessConfig_Uncorrected_EtaDiff(kCent, kEtaDiffOn, templNeRev, dataListNeRev, inputFileNameNeRev + "_Uncorrected", etaBinsPos);

    // Combine O-O datasets (neg+pos)
    std::cout << "\n========== Combining O-O Datasets ==========" << std::endl;
    std::string oo_neg_file = "./TemplateFit/Uncorrected_FourierFit/EtaDiff/Vn_Combined_LHC25ae_pass2_604826_Cent_0_20_Cent.root";
    std::string oo_pos_file = "./TemplateFit/Uncorrected_FourierFit/EtaDiff/Vn_Combined_LHC25ae_pass2_604830_Cent_0_20_Cent.root";
    std::string oo_combined_file = "./TemplateFit/Uncorrected_FourierFit/EtaDiff/Vn_Combined_LHC25ae_pass2_604826_vs_LHC25ae_pass2_604830_Cent.root";
    CombineSystemPairs_Uncorrected(oo_neg_file, oo_pos_file, oo_combined_file, "gV2_Uncorrected_Combined");

    // Combine Ne-Ne datasets (std+rev)
    std::cout << "\n========== Combining Ne-Ne Datasets ==========" << std::endl;
    std::string ne_std_file = "./TemplateFit/Uncorrected_FourierFit/EtaDiff/Vn_Combined_LHC25af_pass2_611697_Cent_0_20_Cent.root";
    std::string ne_rev_file = "./TemplateFit/Uncorrected_FourierFit/EtaDiff/Vn_Combined_LHC25af_pass2_604820_Cent_0_20_Cent.root";
    std::string ne_combined_file = "./TemplateFit/Uncorrected_FourierFit/EtaDiff/Vn_Combined_LHC25af_pass2_611697_vs_LHC25af_pass2_604820_Cent.root";
    CombineSystemPairs_Uncorrected(ne_std_file, ne_rev_file, ne_combined_file, "gV2_Uncorrected_Combined");

    std::cout << "========================================" << std::endl;
    std::cout << "Uncorrected Fourier Fit Analysis Complete" << std::endl;
    std::cout << "Output: TemplateFit/Uncorrected_FourierFit/EtaDiff/" << std::endl;
    std::cout << "========================================" << std::endl;
}

void ProcessConfig_Uncorrected(Bool_t isNch, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName) {
    std::string splitName = isNch ? "Mult" : "Cent";
    std::cout << "\n========== Processing Uncorrected (Integrated) ==========" << std::endl;
    std::cout << "Dataset: " << templ.fileNameSuffix << " (" << splitName << " " << templ.minRange << "-" << templ.maxRange << "%)" << std::endl;
}

void ProcessConfig_Uncorrected_EtaDiff(Bool_t isNch, Bool_t isEtaDiff, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName, const std::vector<float>& etaBins) {
    std::string splitName = isNch ? "Mult" : "Cent";
    std::cout << "\n========== Processing Uncorrected (EtaDiff) ==========" << std::endl;
    std::cout << "Dataset: " << templ.fileNameSuffix << " (" << splitName << " " << templ.minRange << "-" << templ.maxRange << "%)" << std::endl;

    gSystem->mkdir("./TemplateFit/Uncorrected_FourierFit/EtaDiff", kTRUE);

    // Prepare output files
    std::string negFile = Form("./TemplateFit/Uncorrected_FourierFit/EtaDiff/Vn_%s_%s_%d_%d_Eta_Neg.root", 
                               templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange);
    std::string posFile = Form("./TemplateFit/Uncorrected_FourierFit/EtaDiff/Vn_%s_%s_%d_%d_Eta_Pos.root", 
                               templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange);
    
    TFile* fNeg = new TFile(negFile.c_str(), "RECREATE");
    TFile* fPos = new TFile(posFile.c_str(), "RECREATE");

    // Process each eta bin and write directly to appropriate file
    for (size_t iEta = 0; iEta < etaBins.size() - 1; iEta++) {
        Double_t etaMin = etaBins[iEta];
        Double_t etaMax = etaBins[iEta + 1];

        TGraphErrors* gV2 = FourierFit_Uncorrected_GetGraph(isNch, templ, etaMin, etaMax);
        if (gV2) {
            if (etaMax <= 0) {
                fNeg->cd();
                gV2->Write();
            } else {
                fPos->cd();
                gV2->Write();
            }
            delete gV2;
        }
    }
    
    fNeg->Close();
    fPos->Close();
    delete fNeg;
    delete fPos;

    std::cout << "  ✓ Saved negative eta bins to: " << negFile << std::endl;
    std::cout << "  ✓ Saved positive eta bins to: " << posFile << std::endl;

    // Combine negative and positive eta bins
    std::string outTag = Form("./TemplateFit/Uncorrected_FourierFit/EtaDiff/Vn_Combined_%s_%s_%d_%d", 
                              templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange);
    
    CombineEtaDiffV2Plots_Uncorrected(negFile, posFile, outTag, splitName);
}

TGraphErrors* FourierFit_Uncorrected_GetGraph(Bool_t isNch, InputUnit data, Double_t etaMin, Double_t etaMax) {
    std::string splitName = isNch ? "Mult" : "Cent";
    const Bool_t useDiff = !(etaMin == 0 && etaMax == 0);

    if (useDiff) {
        std::cout << "    Processing eta bin [" << etaMin << ", " << etaMax << "]" << std::endl;
    }

    // Read raw S/M histogram from bootstrap samples (use regular bootstrap, not peripheral subtracted)
    TFile* datafile = nullptr;
    std::string filePath;
    
    if (useDiff) {
        filePath = Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%d_%d_Eta_%0.1f_%0.1f.root", 
                        data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, etaMin, etaMax);
    } else {
        filePath = Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d.root", 
                        data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange);
    }

    datafile = new TFile(filePath.c_str(), "READ");
    if (!datafile || !datafile->IsOpen()) {
        std::cerr << "ERROR: Cannot open file: " << filePath << std::endl;
        if (datafile) delete datafile;
        return nullptr;
    }

    // Collect bootstrap samples and compute v2 for each using max-min method
    std::vector<double> v2Values;
    double sumV2 = 0.0, sumV2Sq = 0.0;
    Int_t nValidSamples = 0;

    for (Int_t sample = 0; sample < maxSample * maxSample; ++sample) {
        TH1D* h = (TH1D*)datafile->Get(Form("bsSample_hPhiSameOverMixed_%d_%d_%d", data.minRange, data.maxRange, sample));
        if (!h) continue;

        // Simple max-min method for uncorrected v2
        double hmax = h->GetMaximum();
        double hmin = h->GetMinimum();
        double denom = hmax + hmin;
        
        if (denom > 0.0 && hmin > 0.0) {
            double v2_sample = (hmax - hmin) / denom / 2.0;
            v2Values.push_back(v2_sample);
            sumV2 += v2_sample;
            sumV2Sq += v2_sample * v2_sample;
            nValidSamples++;
        }
    }

    datafile->Close();
    delete datafile;

    if (nValidSamples == 0) {
        std::cerr << "No valid samples found for " << data.fileNameSuffix << " eta [" << etaMin << ", " << etaMax << "]" << std::endl;
        return nullptr;
    }

    // Compute mean and error from bootstrap
    double meanV2 = sumV2 / nValidSamples;
    double errV2 = (nValidSamples > 1) ? TMath::Sqrt((sumV2Sq / nValidSamples) - (meanV2 * meanV2)) : 0;

    std::cout << "    v2 (uncorrected, max-min) = " << meanV2 << " ± " << errV2 
              << " (from " << nValidSamples << " samples)" << std::endl;

    // Create and return graph
    TGraphErrors* gV2 = new TGraphErrors();
    gV2->SetName(Form("gV2_Uncorrected_%0.1f_%0.1f", etaMin, etaMax));
    gV2->SetMarkerStyle(23);
    gV2->SetMarkerColor(kMagenta);
    gV2->SetLineColor(kWhite);
    gV2->SetLineWidth(0);
    gV2->SetPoint(0, (etaMin + etaMax) / 2.0, meanV2);
    gV2->SetPointError(0, 0, errV2);
    
    return gV2;
}

void CombineEtaDiffV2Plots_Uncorrected(const std::string& negFile, const std::string& posFile, const std::string& outTag, const std::string& splitName) {
    TFile* fNeg = TFile::Open(negFile.c_str(), "READ");
    TFile* fPos = TFile::Open(posFile.c_str(), "READ");
    
    if (!fNeg || !fNeg->IsOpen()) {
        std::cerr << "Cannot open negative eta file: " << negFile << std::endl;
        if (fNeg) delete fNeg;
        if (fPos) { fPos->Close(); delete fPos; }
        return;
    }
    if (!fPos || !fPos->IsOpen()) {
        std::cerr << "Cannot open positive eta file: " << posFile << std::endl;
        if (fPos) delete fPos;
        if (fNeg) { fNeg->Close(); delete fNeg; }
        return;
    }

    std::cout << "  Combining uncorrected v2 from neg and pos eta files..." << std::endl;

    TGraphErrors* gCombined = new TGraphErrors();
    gCombined->SetName("gV2_Uncorrected_Combined");
    
    Int_t point = 0;
    
    // Add negative eta points
    TIter nextNeg(fNeg->GetListOfKeys());
    TKey* keyNeg;
    while ((keyNeg = (TKey*)nextNeg())) {
        if (TString(keyNeg->GetName()).Contains("gV2_Uncorrected")) {
            TGraphErrors* g = (TGraphErrors*)fNeg->Get(keyNeg->GetName());
            if (g && g->GetN() > 0) {
                Double_t x, y;
                g->GetPoint(0, x, y);
                gCombined->SetPoint(point, x, y);
                gCombined->SetPointError(point, 0, g->GetErrorY(0));
                std::cout << "    Added negative eta point " << point << ": eta=" << x << ", v2=" << y << std::endl;
                point++;
            }
        }
    }
    
    // Add positive eta points
    TIter nextPos(fPos->GetListOfKeys());
    TKey* keyPos;
    while ((keyPos = (TKey*)nextPos())) {
        if (TString(keyPos->GetName()).Contains("gV2_Uncorrected")) {
            TGraphErrors* g = (TGraphErrors*)fPos->Get(keyPos->GetName());
            if (g && g->GetN() > 0) {
                Double_t x, y;
                g->GetPoint(0, x, y);
                gCombined->SetPoint(point, x, y);
                gCombined->SetPointError(point, 0, g->GetErrorY(0));
                std::cout << "    Added positive eta point " << point << ": eta=" << x << ", v2=" << y << std::endl;
                point++;
            }
        }
    }

    // Sort by x-coordinate
    gCombined->Sort();
    
    // Set line properties for combined graph - NO LINES
    gCombined->SetLineColor(kWhite);
    gCombined->SetLineWidth(0);
    gCombined->SetMarkerStyle(23);
    gCombined->SetMarkerColor(kMagenta);
    gCombined->SetMarkerSize(1.2);
    
    std::string outPath = Form("%s_Cent.root", outTag.c_str());
    TFile* outfile = new TFile(outPath.c_str(), "RECREATE");
    gCombined->Write();
    outfile->Close();
    delete outfile;
    
    std::cout << "  ✓ Combined " << point << " eta points → " << outPath << std::endl;

    fNeg->Close();
    fPos->Close();
    delete fNeg;
    delete fPos;
}

void CombineSystemPairs_Uncorrected(const std::string& system1File, const std::string& system2File, const std::string& outFile, const std::string& graphName) {
    std::cout << "  Combining " << system1File << " + " << system2File << std::endl;
    
    TFile* f1 = TFile::Open(system1File.c_str(), "READ");
    TFile* f2 = TFile::Open(system2File.c_str(), "READ");
    
    if (!f1 || !f1->IsOpen() || !f2 || !f2->IsOpen()) {
        std::cerr << "  ✗ Error: Cannot open system files" << std::endl;
        if (f1) { f1->Close(); delete f1; }
        if (f2) { f2->Close(); delete f2; }
        return;
    }
    
    TGraphErrors* g1 = (TGraphErrors*)f1->Get(graphName.c_str());
    TGraphErrors* g2 = (TGraphErrors*)f2->Get(graphName.c_str());
    
    if (!g1 || !g2 || g1->GetN() == 0 || g2->GetN() == 0) {
        std::cerr << "  ✗ Error: Cannot find or empty graph in files" << std::endl;
        if (f1) { f1->Close(); delete f1; }
        if (f2) { f2->Close(); delete f2; }
        return;
    }
    
    // Create combined graph with all points
    TGraphErrors* gCombined = new TGraphErrors();
    gCombined->SetName(graphName.c_str());
    
    Int_t point = 0;
    
    // Add all points from first graph
    for (Int_t i = 0; i < g1->GetN(); ++i) {
        Double_t x, y;
        g1->GetPoint(i, x, y);
        gCombined->SetPoint(point, x, y);
        gCombined->SetPointError(point, 0, g1->GetErrorY(i));
        point++;
    }
    
    // Add all points from second graph
    for (Int_t i = 0; i < g2->GetN(); ++i) {
        Double_t x, y;
        g2->GetPoint(i, x, y);
        gCombined->SetPoint(point, x, y);
        gCombined->SetPointError(point, 0, g2->GetErrorY(i));
        point++;
    }
    
    // Sort by x-coordinate
    gCombined->Sort();
    
    // Set line properties
    gCombined->SetLineColor(kWhite);
    gCombined->SetLineWidth(0);
    gCombined->SetMarkerStyle(23);
    gCombined->SetMarkerColor(kMagenta);
    
    // Write to file
    TFile* outfile = TFile::Open(outFile.c_str(), "RECREATE");
    outfile->cd();
    gCombined->Write();
    outfile->Close();
    delete outfile;
    
    f1->Close();
    f2->Close();
    delete f1;
    delete f2;
    
    std::cout << "  ✓ Combined " << point << " total eta points → " << outFile << std::endl;
}
