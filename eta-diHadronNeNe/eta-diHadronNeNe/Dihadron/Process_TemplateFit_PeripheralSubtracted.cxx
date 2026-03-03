/*
 * @Author: Martijn Laarhoven (martijn.laarhoven@cern.ch)
 * @Date: 2026-02-06
 * @Last Modified by: Martijn Laarhoven
 * @Last Modified time: 2026-02-06
 */
//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

// CRITICAL: Set batch mode BEFORE including ROOT GUI headers
#include "TROOT.h"
#include "TApplication.h"

// Initialize batch mode early to prevent GUI initialization
namespace {
    struct BatchInit {
        BatchInit() {
            gROOT->SetBatch(kTRUE);
            gROOT->SetStyle("Plain");
        }
    } batchInit;
}

#include "TFile.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include "./include/BasicForDihadron.h"
#include "./include/Bootstrap.h"

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

// Forward declarations
void ProcessConfig_PeripheralSubtracted(Bool_t isNch, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName);
void ProcessConfig_PtDiff_PeripheralSubtracted(Bool_t isNch, Bool_t isEtaDiff, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName, const std::vector<float>& etaBins);
VnUnit* TemplateFit_PeripheralSubtracted(Bool_t isNch, InputUnit templ, InputUnit data, Bool_t cn2Tovn2, Double_t etaMin=0, Double_t etaMax=0);
void PlotPeripheralSubtractedFit(TH1D* hAvg, Bool_t isNch, const std::string& fileSuffix, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax, const VnUnit* vn, double baseline);
void CreateMultiPanelEtaDiffPlots(Bool_t isNch, const InputUnit& data, const std::vector<float>& etaBins, const std::string& outputFileTag);
void CombineEtaDiffV2Plots(const std::string& negFile, const std::string& posFile, const std::string& outTag, const std::string& splitName);
void CompareCollisionSystems(const std::string& ooFile, const std::string& neFile, const std::string& splitName);
void CombineEtaDiffV2Plots(const std::string& negFile, const std::string& posFile, const std::string& outTag, const std::string& splitName);

// Helper function to set collision system name based on filename
std::string GetCollisionSystemName(const std::string& filename) {
    if (filename.find("ae") != std::string::npos) {
        return "O-O";
    } else if (filename.find("af") != std::string::npos) {
        return "Ne-Ne";
    }
    return "Unknown";
}

std::string collisionSystemName = "peripheral NeNe";

void Process_TemplateFit_PeripheralSubtracted() {
    // Batch mode already set at global initialization
    gStyle->SetOptStat(0);
    
    // Define eta bins for both configurations
    std::vector<float> etaBinsNeg = {-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0};
    std::vector<float> etaBinsPos = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    
    std::vector<InputUnit> inputList;
    std::string inputFileNameNew = "LHC25ae_pass2_604830";
    std::string inputFileNameNeStd = "LHC25af_pass2_611697";
    std::string inputFileNameNeRev = "LHC25af_pass2_604820";

    // O-O standard configuration
    std::string inputFileName = "LHC25ae_pass2_604826";
    collisionSystemName = GetCollisionSystemName(inputFileName);
    
    InputUnit templ(inputFileName, kCent, kEtaDiffOff, 80, 100);
    std::vector<InputUnit> dataList;
    dataList.push_back(InputUnit(inputFileName, kCent, kEtaDiffOff, 0, 20));

    ProcessConfig_PeripheralSubtracted(kCent, templ, dataList, inputFileName + "_PeripheralSubtracted");
    
    // Eta-diff version
    dataList.clear();
    dataList.push_back(InputUnit(inputFileName, kCent, kEtaDiffOn, 0, 20));
    ProcessConfig_PtDiff_PeripheralSubtracted(kCent, kEtaDiffOn, templ, dataList, inputFileName + "_PeripheralSubtracted", etaBinsNeg);

    // O-O reversed configuration
    collisionSystemName = GetCollisionSystemName(inputFileNameNew);

    InputUnit templNew(inputFileNameNew, kCent, kEtaDiffOff, 80, 100);
    dataList.clear();
    dataList.push_back(InputUnit(inputFileNameNew, kCent, kEtaDiffOff, 0, 20));
    ProcessConfig_PeripheralSubtracted(kCent, templNew, dataList, inputFileNameNew + "_PeripheralSubtracted");

    dataList.clear();
    dataList.push_back(InputUnit(inputFileNameNew, kCent, kEtaDiffOn, 0, 20));
    ProcessConfig_PtDiff_PeripheralSubtracted(kCent, kEtaDiffOn, templNew, dataList, inputFileNameNew + "_PeripheralSubtracted", etaBinsPos);

    // Ne-Ne standard configuration
    collisionSystemName = GetCollisionSystemName(inputFileNameNeStd);

    InputUnit templNeStd(inputFileNameNeStd, kCent, kEtaDiffOff, 80, 100);
    dataList.clear();
    dataList.push_back(InputUnit(inputFileNameNeStd, kCent, kEtaDiffOff, 0, 20));
    ProcessConfig_PeripheralSubtracted(kCent, templNeStd, dataList, inputFileNameNeStd + "_PeripheralSubtracted");

    dataList.clear();
    dataList.push_back(InputUnit(inputFileNameNeStd, kCent, kEtaDiffOn, 0, 20));
    ProcessConfig_PtDiff_PeripheralSubtracted(kCent, kEtaDiffOn, templNeStd, dataList, inputFileNameNeStd + "_PeripheralSubtracted", etaBinsNeg);

    // Ne-Ne reversed configuration
    collisionSystemName = GetCollisionSystemName(inputFileNameNeRev);

    InputUnit templNeRev(inputFileNameNeRev, kCent, kEtaDiffOff, 80, 100);
    dataList.clear();
    dataList.push_back(InputUnit(inputFileNameNeRev, kCent, kEtaDiffOff, 0, 20));
    ProcessConfig_PeripheralSubtracted(kCent, templNeRev, dataList, inputFileNameNeRev + "_PeripheralSubtracted");

    dataList.clear();
    dataList.push_back(InputUnit(inputFileNameNeRev, kCent, kEtaDiffOn, 0, 20));
    ProcessConfig_PtDiff_PeripheralSubtracted(kCent, kEtaDiffOn, templNeRev, dataList, inputFileNameNeRev + "_PeripheralSubtracted", etaBinsPos);

    // Generate combined plots for each collision system
    const std::string splitName = "Cent";
    
    // Plot 1: O-O combined (negative + positive eta)
    std::cout << "\n========== CREATING O-O COMBINED PLOT ==========" << std::endl;
    CombineEtaDiffV2Plots(
        Form("./TemplateFit/PeripheralSubtracted/EtaDiff/Vn_%s_PeripheralSubtracted_%s.root", inputFileName.c_str(), splitName.c_str()),
        Form("./TemplateFit/PeripheralSubtracted/EtaDiff/Vn_%s_PeripheralSubtracted_%s.root", inputFileNameNew.c_str(), splitName.c_str()),
        Form("%s_vs_%s", inputFileName.c_str(), inputFileNameNew.c_str()),
        splitName
    );

    // Plot 2: Ne-Ne combined (negative + positive eta) - only if files exist
    std::cout << "\n========== CREATING Ne-Ne COMBINED PLOT ==========" << std::endl;
    TFile* checkNeStd = TFile::Open(Form("./TemplateFit/PeripheralSubtracted/EtaDiff/Vn_%s_PeripheralSubtracted_%s.root", inputFileNameNeStd.c_str(), splitName.c_str()), "READ");
    TFile* checkNeRev = TFile::Open(Form("./TemplateFit/PeripheralSubtracted/EtaDiff/Vn_%s_PeripheralSubtracted_%s.root", inputFileNameNeRev.c_str(), splitName.c_str()), "READ");
    
    if (checkNeStd && checkNeStd->IsOpen() && checkNeRev && checkNeRev->IsOpen()) {
        if (checkNeStd) { checkNeStd->Close(); delete checkNeStd; }
        if (checkNeRev) { checkNeRev->Close(); delete checkNeRev; }
        
        CombineEtaDiffV2Plots(
            Form("./TemplateFit/PeripheralSubtracted/EtaDiff/Vn_%s_PeripheralSubtracted_%s.root", inputFileNameNeStd.c_str(), splitName.c_str()),
            Form("./TemplateFit/PeripheralSubtracted/EtaDiff/Vn_%s_PeripheralSubtracted_%s.root", inputFileNameNeRev.c_str(), splitName.c_str()),
            Form("%s_vs_%s", inputFileNameNeStd.c_str(), inputFileNameNeRev.c_str()),
            splitName
        );

        // Plot 3: O-O vs Ne-Ne comparison - only if both combined plots exist
        std::cout << "\n========== CREATING O-O vs Ne-Ne COMPARISON PLOT ==========" << std::endl;
        CompareCollisionSystems(
            Form("./TemplateFit/PeripheralSubtracted/EtaDiff/Vn_Combined_%s_vs_%s_%s.root", inputFileName.c_str(), inputFileNameNew.c_str(), splitName.c_str()),
            Form("./TemplateFit/PeripheralSubtracted/EtaDiff/Vn_Combined_%s_vs_%s_%s.root", inputFileNameNeStd.c_str(), inputFileNameNeRev.c_str(), splitName.c_str()),
            splitName
        );
    } else {
        std::cout << "Ne-Ne standard configuration files not found. Skipping Ne-Ne combined and comparison plots." << std::endl;
        if (checkNeStd) { checkNeStd->Close(); delete checkNeStd; }
        if (checkNeRev) { checkNeRev->Close(); delete checkNeRev; }
    }
}

void ProcessConfig_PeripheralSubtracted(Bool_t isNch, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName) {
    std::sort(dataList.begin(), dataList.end(), [](const InputUnit& a, const InputUnit& b) {
        return a.minRange < b.minRange;
    });
    
    std::cout << "Data list (PeripheralSubtracted): " << std::endl;
    for (auto input : dataList) {
        std::cout << "  " << input.fileNameSuffix << " " << input.minRange << " " << input.maxRange << std::endl;
    }

    std::vector<VnUnit*> vnResults;
    std::vector<double> xCenters;
    std::vector<double> xErrors;
    
    for (auto& data : dataList) {
        VnUnit* res = TemplateFit_PeripheralSubtracted(isNch, templ, data, kTRUE);
        if (!res) {
            std::cerr << "Skip bin " << data.minRange << " to " << data.maxRange << " (fit failed or missing files)." << std::endl;
            continue;
        }
        vnResults.push_back(res);
        xCenters.push_back(0.5 * (data.minRange + data.maxRange));
        xErrors.push_back(0.5 * (data.maxRange - data.minRange));
    }

    if (vnResults.empty()) {
        std::cerr << "No valid results!" << std::endl;
        return;
    }

    std::string splitName = isNch ? "Mult" : "Cent";
    
    // Create output histograms
    Int_t nBins = vnResults.size();
    TH1D* hV2 = new TH1D("hV2", Form("v_{2};Centrality;v_{2}"), nBins, 0, nBins);
    TH1D* hV3 = new TH1D("hV3", Form("v_{3};Centrality;v_{3}"), nBins, 0, nBins);
    TH1D* hV4 = new TH1D("hV4", Form("v_{4};Centrality;v_{4}"), nBins, 0, nBins);
    
    for (Int_t i = 0; i < (Int_t)vnResults.size(); ++i) {
        hV2->SetBinContent(i + 1, vnResults[i]->v2);
        hV2->SetBinError(i + 1, vnResults[i]->v2_err);
        hV3->SetBinContent(i + 1, vnResults[i]->v3);
        hV3->SetBinError(i + 1, vnResults[i]->v3_err);
        hV4->SetBinContent(i + 1, vnResults[i]->v4);
        hV4->SetBinError(i + 1, vnResults[i]->v4_err);
    }

    // Create V2Δ graph
    Int_t nGraphBins = static_cast<Int_t>(vnResults.size());
    if (nGraphBins <= 0) {
        std::cerr << "No bins for graph!" << std::endl;
        delete hV2;
        delete hV3;
        delete hV4;
        for (auto res : vnResults) delete res;
        return;
    }

    TGraphErrors* gV2 = new TGraphErrors(nGraphBins);
    gV2->SetName("gV2Delta");
    gV2->SetTitle(Form("v_{2#Delta};Centrality;v_{2#Delta}"));
    for (Int_t i = 0; i < nGraphBins; ++i) {
        gV2->SetPoint(i, xCenters[i], vnResults[i]->v2);
        gV2->SetPointError(i, xErrors[i], vnResults[i]->v2_err);
    }

    // Write to file
    gSystem->mkdir("./TemplateFit/PeripheralSubtracted", kTRUE);
    std::string outPath = Form("./TemplateFit/PeripheralSubtracted/Vn_%s_%s.root",
                               outputFileName.c_str(), splitName.c_str());
    TFile* outputFile = TFile::Open(outPath.c_str(), "RECREATE");
    if (!outputFile || !outputFile->IsOpen()) {
        std::cerr << "Cannot open output file: " << outPath << std::endl;
        delete gV2;
        delete hV2;
        delete hV3;
        delete hV4;
        for (auto res : vnResults) delete res;
        return;
    }
    
    hV2->SetDirectory(nullptr);
    hV3->SetDirectory(nullptr);
    hV4->SetDirectory(nullptr);
    
    hV2->Write();
    hV3->Write();
    hV4->Write();
    gV2->Write();

    std::cout << "Output file: " << outPath << std::endl;
    outputFile->Close();
    delete outputFile;
    
    // Save V2Δ plot
    gSystem->mkdir("./TemplateFit/PeripheralSubtracted/PDFs", kTRUE);
    TCanvas* cV2 = new TCanvas("cV2Delta_PeripheralSub", "cV2Delta_PeripheralSub", 800, 600);
    gV2->SetMarkerStyle(20);
    gV2->SetMarkerColor(kBlack);
    gV2->SetLineColor(kBlack);
    gV2->Draw("AP");
    cV2->SaveAs(Form("./TemplateFit/PeripheralSubtracted/PDFs/V2Delta_%s_%s.pdf",
                     outputFileName.c_str(), splitName.c_str()));
    cV2->Close();
    delete cV2;

    // Cleanup
    delete gV2;
    delete hV2;
    delete hV3;
    delete hV4;
    for (auto res : vnResults) {
        delete res;
    }
}

void CreateMultiPanelEtaDiffPlots(Bool_t isNch, const InputUnit& data, const std::vector<float>& etaBins, const std::string& outputFileTag) {
    std::string splitName = isNch ? "Mult" : "Cent";
    const std::string systemName = GetCollisionSystemName(data.fileNameSuffix);
    Int_t nEtaBins = etaBins.size() - 1;
    
    // Determine grid layout (prefer 4x4, 3x3, etc.)
    Int_t nCols = (Int_t)std::ceil(std::sqrt(nEtaBins));
    Int_t nRows = (Int_t)std::ceil((double)nEtaBins / nCols);
    
    // Create large canvas for all eta bins
    TCanvas* cAll = new TCanvas("cAllEtaBins_ZYAM", "All Eta Bins ZYAM", 1600, 1200);
    cAll->Divide(nCols, nRows, 0.005, 0.005);
    
    gStyle->SetOptStat(0);
    
    for (Int_t iEta = 0; iEta < nEtaBins; ++iEta) {
        double etaMin = etaBins[iEta];
        double etaMax = etaBins[iEta + 1];
        
        // Load peripheral-subtracted histogram
        std::string filePath = Form("./ProcessOutput/PeripheralSubtraction/EtaDiff/BootstrapSample_%s_%s_%d_%d_Eta_%0.1f_%0.1f_PeripheralSub.root",
                                     data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, etaMin, etaMax);
        TFile* datafile = TFile::Open(filePath.c_str(), "READ");
        if (!datafile || !datafile->IsOpen()) {
            std::cerr << "Cannot open: " << filePath << std::endl;
            continue;
        }
        
        // Get average histogram
        TH1D* hAvg = nullptr;
        Int_t nSamples = 0;
        for (Int_t sample = 0; sample < maxSample * maxSample; ++sample) {
            TH1D* h = (TH1D*)datafile->Get(Form("hPhiSameOverMixed_%d_%d_%d", data.minRange, data.maxRange, sample));
            if (!h) continue;
            
            if (!hAvg) {
                hAvg = (TH1D*)h->Clone(Form("hAvg_eta%d", iEta));
                hAvg->SetDirectory(nullptr);
                hAvg->Sumw2();
            } else {
                hAvg->Add(h);
            }
            nSamples++;
        }
        
        if (hAvg && nSamples > 0) {
            hAvg->Scale(1.0 / nSamples);
        }
        
        datafile->Close();
        delete datafile;
        
        if (!hAvg || hAvg->GetEntries() == 0) {
            std::cerr << "No valid data for eta [" << etaMin << ", " << etaMax << "]" << std::endl;
            continue;
        }
        
        // Create fit
        TF1* fit = new TF1(Form("fit_eta%d", iEta),
            "[0]*(1 + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x))",
            -TMath::Pi()/2.0, 1.5*TMath::Pi());
        
        double hmin = hAvg->GetMinimum();
        double hmax = hAvg->GetMaximum();
        double baselineGuess = (hmin + hmax) / 2.0;
        double v2Guess = (hmax - hmin) / (hmax + hmin);
        
        fit->SetParameters(baselineGuess, v2Guess, 0.0, 0.0);
        fit->SetParLimits(0, hmin * 0.5, hmax * 1.5);
        fit->SetParLimits(1, -0.3, 0.4);
        fit->SetParLimits(2, -0.15, 0.15);
        fit->SetParLimits(3, -0.15, 0.15);
        fit->SetRange(-TMath::Pi()/2.0, 1.5*TMath::Pi());
        hAvg->Fit(fit, "RQ0");
        
        // Draw in panel
        cAll->cd(iEta + 1);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.12);
        gPad->SetTopMargin(0.08);
        gPad->SetTicks(1,1);
        
        hAvg->SetMarkerStyle(20);
        hAvg->SetMarkerColor(kBlack);
        hAvg->SetLineColor(kBlack);
        hAvg->SetMarkerSize(0.6);
        hAvg->GetXaxis()->SetTitle("#Delta#phi");
        hAvg->GetYaxis()->SetTitle("C(#Delta#phi)");
        hAvg->GetXaxis()->SetTitleSize(0.06);
        hAvg->GetYaxis()->SetTitleSize(0.06);
        hAvg->GetXaxis()->SetLabelSize(0.05);
        hAvg->GetYaxis()->SetLabelSize(0.05);
        hAvg->Draw("P");
        
        fit->SetLineColor(kRed+1);
        fit->SetLineWidth(2);
        fit->Draw("same");
        
        // Draw baseline
        TF1* baselineFunc = new TF1(Form("baselineFunc_eta%d", iEta), "[0]", -TMath::Pi()/2.0, 1.5*TMath::Pi());
        baselineFunc->SetParameter(0, fit->GetParameter(0));
        baselineFunc->SetLineColor(kBlue+1);
        baselineFunc->SetLineWidth(1);
        baselineFunc->SetLineStyle(2);
        baselineFunc->Draw("same");
        
        // Add eta range label
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(43);
        latex.SetTextSize(14);
        latex.DrawLatex(0.20, 0.85, Form("%.1f<#eta<%.1f", etaMin, etaMax));
        latex.DrawLatex(0.20, 0.75, Form("v_{2}=%.3f", fit->GetParameter(1)));
        
        delete baselineFunc;
        delete fit;
        delete hAvg;
    }
    
    // Save multi-panel plot
    gSystem->mkdir("./TemplateFit/PeripheralSubtracted/EtaDiff/PDFs", kTRUE);
    cAll->SaveAs(Form("./TemplateFit/PeripheralSubtracted/EtaDiff/PDFs/AllEtaBins_%s_%s_%d_%d.pdf",
                      data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange));
    
    // Properly clean up canvas and pads
    cAll->SetBit(kCannotPick);
    cAll->Close();
    delete cAll;
    gROOT->GetListOfCanvases()->Delete();
    
    std::cout << "Created multi-panel plot: AllEtaBins_" << data.fileNameSuffix << "_" << splitName 
              << "_" << data.minRange << "_" << data.maxRange << ".pdf" << std::endl;
}

void ProcessConfig_PtDiff_PeripheralSubtracted(Bool_t isNch, Bool_t isEtaDiff, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName, const std::vector<float>& etaBins) {
    std::sort(dataList.begin(), dataList.end(), [](const InputUnit& a, const InputUnit& b) {
        return a.minRange < b.minRange;
    });
    
    std::cout << "Data list (PeripheralSubtracted, EtaDiff): " << std::endl;
    for (auto input : dataList) {
        std::cout << "  " << input.fileNameSuffix << " " << input.minRange << " " << input.maxRange << std::endl;
    }

    std::vector<VnUnit*> vnResults;
    std::vector<double> xCenters;
    std::vector<double> xErrors;
    
    // Store canvases for multi-panel plot
    std::vector<TCanvas*> allCanvases;

    // Process each eta bin
    for (Int_t iEta = 0; iEta < (Int_t)etaBins.size() - 1; iEta++) {
        double etaMin = etaBins[iEta];
        double etaMax = etaBins[iEta + 1];

        for (auto& data : dataList) {
            VnUnit* res = TemplateFit_PeripheralSubtracted(isNch, templ, data, kTRUE, etaMin, etaMax);
            if (!res) {
                std::cerr << "Skip bin " << etaMin << " to " << etaMax << " (fit failed or missing files)." << std::endl;
                continue;
            }
            vnResults.push_back(res);
            xCenters.push_back(0.5 * (etaMin + etaMax));
            xErrors.push_back(0.5 * (etaMax - etaMin));
        }
    }

    if (vnResults.empty()) {
        std::cerr << "No valid results!" << std::endl;
        return;
    }

    std::string splitName = isNch ? "Mult" : "Cent";
    
    // Create output histograms
    Int_t nBins = vnResults.size();
    TH1D* hV2 = new TH1D("hV2", Form("v_{2};#eta_{trig};v_{2}"), nBins, 0, nBins);
    TH1D* hV3 = new TH1D("hV3", Form("v_{3};#eta_{trig};v_{3}"), nBins, 0, nBins);
    TH1D* hV4 = new TH1D("hV4", Form("v_{4};#eta_{trig};v_{4}"), nBins, 0, nBins);

    for (Int_t i = 0; i < nBins; ++i) {
        hV2->SetBinContent(i + 1, vnResults[i]->v2);
        hV2->SetBinError(i + 1, vnResults[i]->v2_err);
        hV3->SetBinContent(i + 1, vnResults[i]->v3);
        hV3->SetBinError(i + 1, vnResults[i]->v3_err);
        hV4->SetBinContent(i + 1, vnResults[i]->v4);
        hV4->SetBinError(i + 1, vnResults[i]->v4_err);
    }

    // Create V2Δ graph
    TGraphErrors* gV2 = new TGraphErrors(nBins);
    gV2->SetName("gV2Delta");
    gV2->SetTitle(Form("v_{2#Delta};#eta_{trig};v_{2#Delta}"));
    for (Int_t i = 0; i < nBins; ++i) {
        gV2->SetPoint(i, xCenters[i], vnResults[i]->v2);
        gV2->SetPointError(i, xErrors[i], vnResults[i]->v2_err);
    }

    // Write to file
    gSystem->mkdir("./TemplateFit/PeripheralSubtracted/EtaDiff", kTRUE);
    TFile outputFile(Form("./TemplateFit/PeripheralSubtracted/EtaDiff/Vn_%s_%s.root",
                          outputFileName.c_str(), splitName.c_str()), "RECREATE");
    hV2->Write();
    hV3->Write();
    hV4->Write();
    gV2->Write();

    std::cout << "Output file: ./TemplateFit/PeripheralSubtracted/EtaDiff/Vn_" << outputFileName << "_" << splitName << ".root" << std::endl;
    
    // Save V2Δ plot
    gSystem->mkdir("./TemplateFit/PeripheralSubtracted/EtaDiff/PDFs", kTRUE);
    TCanvas* cV2 = new TCanvas("cV2Delta_PeripheralSub", "cV2Delta_PeripheralSub", 800, 600);
    gV2->SetMarkerStyle(20);
    gV2->SetMarkerColor(kBlack);
    gV2->SetLineColor(kBlack);
    gV2->Draw("AP");
    cV2->SaveAs(Form("./TemplateFit/PeripheralSubtracted/EtaDiff/PDFs/V2Delta_%s_%s.pdf",
                     outputFileName.c_str(), splitName.c_str()));
    delete cV2;
    outputFile.Close();

    // Multi-panel overview plot removed - individual eta bin plots are more informative
    // for (auto& data : dataList) {
    //     CreateMultiPanelEtaDiffPlots(isNch, data, etaBins, outputFileName);
    // }

    // Cleanup
    delete hV2;
    delete hV3;
    delete hV4;
    delete gV2;
    for (auto res : vnResults) delete res;
}

void CombineEtaDiffV2Plots(const std::string& negFile, const std::string& posFile, const std::string& outTag, const std::string& splitName) {
    TFile* fNeg = TFile::Open(negFile.c_str(), "READ");
    TFile* fPos = TFile::Open(posFile.c_str(), "READ");

    if (!fNeg || !fNeg->IsOpen() || !fPos || !fPos->IsOpen()) {
        std::cerr << "Cannot open eta-diff files for combined plot: " << negFile << " or " << posFile << std::endl;
        if (fNeg) { fNeg->Close(); delete fNeg; }
        if (fPos) { fPos->Close(); delete fPos; }
        return;
    }

    TGraphErrors* gNeg = (TGraphErrors*)fNeg->Get("gV2Delta");
    TGraphErrors* gPos = (TGraphErrors*)fPos->Get("gV2Delta");
    if (!gNeg || !gPos) {
        std::cerr << "Missing gV2Delta in one of the eta-diff files." << std::endl;
        fNeg->Close();
        fPos->Close();
        delete fNeg;
        delete fPos;
        return;
    }

    struct Point { double x, y, ex, ey; };
    std::vector<Point> points;

    const Int_t nNeg = gNeg->GetN();
    for (Int_t i = 0; i < nNeg; ++i) {
        double x, y;
        gNeg->GetPoint(i, x, y);
        points.push_back({x, y, gNeg->GetErrorX(i), gNeg->GetErrorY(i)});
    }

    const Int_t nPos = gPos->GetN();
    for (Int_t i = 0; i < nPos; ++i) {
        double x, y;
        gPos->GetPoint(i, x, y);
        points.push_back({x, y, gPos->GetErrorX(i), gPos->GetErrorY(i)});
    }

    std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) { return a.x < b.x; });

    TGraphErrors* gCombined = new TGraphErrors(points.size());
    gCombined->SetName("gV2Delta_Combined");
    gCombined->SetTitle("v_{2#Delta};#eta_{trig};v_{2#Delta}");

    for (size_t i = 0; i < points.size(); ++i) {
        gCombined->SetPoint(i, points[i].x, points[i].y);
        gCombined->SetPointError(i, points[i].ex, points[i].ey);
    }

    gSystem->mkdir("./TemplateFit/PeripheralSubtracted/EtaDiff/PDFs", kTRUE);
    gSystem->mkdir("./TemplateFit/PeripheralSubtracted/EtaDiff", kTRUE);

    TFile outFile(Form("./TemplateFit/PeripheralSubtracted/EtaDiff/Vn_Combined_%s_%s.root", outTag.c_str(), splitName.c_str()), "RECREATE");
    gCombined->Write();
    outFile.Close();

    TCanvas* c = new TCanvas("cV2Delta_Combined", "cV2Delta_Combined", 800, 600);
    gCombined->SetMarkerStyle(20);
    gCombined->SetMarkerColor(kBlack);
    gCombined->SetLineColor(kBlack);
    gCombined->Draw("AP");
    gCombined->GetXaxis()->SetLimits(-0.8, 0.8);
    gCombined->GetXaxis()->SetRangeUser(-0.8, 0.8);
    c->SaveAs(Form("./TemplateFit/PeripheralSubtracted/EtaDiff/PDFs/V2Delta_Combined_%s_%s.pdf", outTag.c_str(), splitName.c_str()));

    // Properly clean up canvas
    c->SetBit(kCannotPick);
    c->Close();
    delete c;
    gROOT->GetListOfCanvases()->Delete();
    
    delete gCombined;
    fNeg->Close();
    fPos->Close();
    delete fNeg;
    delete fPos;
}

void CompareCollisionSystems(const std::string& ooFile, const std::string& neFile, const std::string& splitName) {
    TFile* fOO = TFile::Open(ooFile.c_str(), "READ");
    TFile* fNe = TFile::Open(neFile.c_str(), "READ");

    if (!fOO || !fOO->IsOpen() || !fNe || !fNe->IsOpen()) {
        std::cerr << "Cannot open collision system files for comparison: " << ooFile << " or " << neFile << std::endl;
        if (fOO) { fOO->Close(); delete fOO; }
        if (fNe) { fNe->Close(); delete fNe; }
        return;
    }

    TGraphErrors* gOO = (TGraphErrors*)fOO->Get("gV2Delta_Combined");
    TGraphErrors* gNe = (TGraphErrors*)fNe->Get("gV2Delta_Combined");
    if (!gOO || !gNe) {
        std::cerr << "Missing gV2Delta_Combined in one of the collision system files." << std::endl;
        fOO->Close();
        fNe->Close();
        delete fOO;
        delete fNe;
        return;
    }

    // Create comparison plot
    TCanvas* c = new TCanvas("cCompare", "cCompare", 900, 600);
    
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();

    auto updateMinMax = [&](TGraphErrors* g) {
        const Int_t n = g->GetN();
        for (Int_t i = 0; i < n; ++i) {
            double x, y;
            g->GetPoint(i, x, y);
            double ey = g->GetErrorY(i);
            minY = std::min(minY, y - ey);
            maxY = std::max(maxY, y + ey);
        }
    };

    updateMinMax(gOO);
    updateMinMax(gNe);

    if (minY == std::numeric_limits<double>::max()) {
        minY = 0.0;
        maxY = 0.01;
    }

    double margin = 0.1 * (maxY - minY);
    if (margin <= 0) margin = 0.001;

    gOO->SetTitle(Form("v_{2#Delta};#eta_{trig};v_{2#Delta}"));
    gOO->GetXaxis()->SetTitle("#eta_{trig}");
    gOO->GetYaxis()->SetTitle("v_{2#Delta}");
    gOO->GetXaxis()->SetLimits(-0.8, 0.8);
    gOO->GetXaxis()->SetRangeUser(-0.8, 0.8);
    gOO->GetYaxis()->SetRangeUser(minY - margin, maxY + margin);
    gOO->Draw("AP");

    gNe->SetMarkerStyle(21);
    gNe->SetMarkerColor(kBlue);
    gNe->SetLineColor(kBlue);
    gNe->Draw("P SAME");

    // Add legend
    TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(gOO, "O-O", "p");
    leg->AddEntry(gNe, "Ne-Ne", "p");
    leg->Draw();

    gSystem->mkdir("./TemplateFit/PeripheralSubtracted/EtaDiff/PDFs", kTRUE);
    c->SaveAs(Form("./TemplateFit/PeripheralSubtracted/EtaDiff/PDFs/V2Delta_Comparison_OO_vs_NeNe_%s.pdf", splitName.c_str()));
    
    std::cout << "Comparison plot created: V2Delta_Comparison_OO_vs_NeNe_" << splitName << ".pdf" << std::endl;

    delete leg;
    
    // Properly clean up canvas
    c->SetBit(kCannotPick);
    c->Close();
    delete c;
    gROOT->GetListOfCanvases()->Delete();
    
    fOO->Close();
    fNe->Close();
    delete fOO;
    delete fNe;
}

void PlotPeripheralSubtractedFit(TH1D* hAvg, Bool_t isNch, const std::string& fileSuffix, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax, const VnUnit* vn, double baseline) {
    if (!hAvg || !vn) return;
    gStyle->SetOptStat(0);

    std::string splitName = isNch ? "Mult" : "Cent";
    const bool useDiff = !(etaMin == 0 && etaMax == 0);
    const std::string systemName = GetCollisionSystemName(fileSuffix);

    // Create canvas with two pads (data + residuals)
    TCanvas* c = new TCanvas("cZYAMTemplateFit", "cZYAMTemplateFit", 800, 600);
    c->Range(0,0,1,1);
    
    // Top pad for data and fit
    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
    pad1->SetBorderMode(0);
    pad1->SetBorderSize(2);
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.05);
    pad1->SetTopMargin(0.05);
    pad1->SetBottomMargin(0.13);
    pad1->SetTicks(1,1);
    pad1->Draw();
    
    // Bottom pad for residuals
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.25);
    pad2->SetBorderMode(0);
    pad2->SetBorderSize(2);
    pad2->SetLeftMargin(0.12);
    pad2->SetRightMargin(0.05);
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.25);
    pad2->SetTicks(1,1);
    pad2->Draw();

    // Draw in top pad
    pad1->cd();
    
    // Create background histogram for axis
    double ymin = hAvg->GetMinimum();
    double ymax = hAvg->GetMaximum();
    double ydelta = ymax - ymin;
    
    TH1D* hbkg = new TH1D("hbkg_zyam", "hbkg_zyam", 1, -TMath::Pi()/2.0, 3*TMath::Pi()/2.0);
    hbkg->SetStats(0);
    hbkg->GetXaxis()->SetTitle("#Delta#phi [rad]");
    hbkg->GetYaxis()->SetTitle("C(#Delta#phi)");
    hbkg->GetXaxis()->SetTitleSize(0.045);
    hbkg->GetXaxis()->SetLabelSize(0.04);
    hbkg->GetYaxis()->SetTitleSize(0.045);
    hbkg->GetYaxis()->SetLabelSize(0.04);
    hbkg->GetXaxis()->SetTitleOffset(1.2);
    hbkg->GetYaxis()->SetTitleOffset(1.2);
    hbkg->GetYaxis()->SetRangeUser(ymin - ydelta*0.2, ymax + ydelta*0.7);
    hbkg->Draw();

    // Plot data
    hAvg->SetMarkerStyle(20);
    hAvg->SetMarkerColor(kBlack);
    hAvg->SetLineColor(kBlack);
    hAvg->SetMarkerSize(1.0);
    hAvg->Draw("same p");

    // Create fit function: C(Δφ) = B*(1 + 2*v2*cos(2Δφ) + 2*v3*cos(3Δφ) + 2*v4*cos(4Δφ))
    TF1* fit = new TF1("fitZYAM",
        "[0]*(1 + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x))",
        -TMath::Pi()/2.0, 1.5*TMath::Pi());
    
    // Set starting parameters from data
    double hmin = hAvg->GetMinimum();
    double hmax = hAvg->GetMaximum();
    double baselineGuess = (hmin + hmax) / 2.0;
    double v2Guess = (hmax - hmin) / (hmax + hmin);
    
    fit->SetParameters(baselineGuess, v2Guess, 0.0, 0.0);
    
    // Constrain baseline to reasonable range around data mean
    fit->SetParLimits(0, hmin * 0.5, hmax * 1.5);
    
    // Constrain vn to physical ranges (slightly relaxed for better fits)
    fit->SetParLimits(1, -0.3, 0.4);  // v2delta typically 0-0.2, allow some flexibility
    fit->SetParLimits(2, -0.15, 0.15);  // v3delta small
    fit->SetParLimits(3, -0.15, 0.15);  // v4delta small
    
    // Use improved fitting options:
    // R = use range, M = MINUIT fitter, E = Minos errors, N = no draw, 0 = don't add to hist list
    fit->SetRange(-TMath::Pi()/2.0, 1.5*TMath::Pi());
    hAvg->Fit(fit, "RQ0");  // Stable fit for bootstrap outputs
    
    fit->SetLineColor(kRed+1);
    fit->SetLineWidth(2);
    fit->Draw("same");
    
    // Draw baseline as horizontal dashed line
    TF1* baselineFunc = new TF1("baselineFunc", "[0]", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    baselineFunc->SetParameter(0, fit->GetParameter(0));
    baselineFunc->SetLineColor(kBlue+1);
    baselineFunc->SetLineWidth(2);
    baselineFunc->SetLineStyle(2);
    baselineFunc->Draw("same");

    // Add legend
    TLegend* leg = new TLegend(0.5, 0.60, 0.9, 0.90);
    leg->SetBorderSize(0);
    leg->AddEntry(hAvg, "Peripheral-subtracted data", "lep");
    leg->AddEntry(fit, "Fit: B(1+2v_{n}cos(n#Delta#phi))", "l");
    leg->AddEntry(baselineFunc, Form("Normalization B = %.4f", fit->GetParameter(0)), "l");
    leg->Draw();

    // Add text labels
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(43);
    latex.SetTextSize(20);
    latex.DrawLatex(0.15, 0.88, Form("ALICE %s", systemName.c_str()));
    if (isNch) {
        latex.DrawLatex(0.15, 0.83, Form("%d < N_{ch} < %d", minRange, maxRange));
    } else {
        latex.DrawLatex(0.15, 0.83, Form("%d < Cent < %d", minRange, maxRange));
    }
    if (useDiff) {
        latex.DrawLatex(0.15, 0.78, Form("%.1f < #eta_{trig} < %.1f", etaMin, etaMax));
    }
    
    // Show fitted parameters
    latex.DrawLatex(0.15, 0.73, Form("v_{2#Delta} = %.4f #pm %.4f", fit->GetParameter(1), fit->GetParError(1)));
    latex.DrawLatex(0.15, 0.69, Form("v_{3#Delta} = %.4f #pm %.4f", fit->GetParameter(2), fit->GetParError(2)));
    latex.DrawLatex(0.15, 0.65, Form("v_{4#Delta} = %.4f #pm %.4f", fit->GetParameter(3), fit->GetParError(3)));
    
    // Calculate chi²/ndf for top pad display
    double chi2 = fit->GetChisquare();
    int ndf = fit->GetNDF();
    double chi2ndf = (ndf > 0) ? chi2 / ndf : 0;
    TLatex chi2Label;
    chi2Label.SetNDC();
    chi2Label.SetTextFont(43);
    chi2Label.SetTextSize(20);
    chi2Label.DrawLatex(0.50, 0.60, Form("#chi^{2}/ndf = %.1f/%d = %.2f", chi2, ndf, chi2ndf));

    // Draw residuals in bottom pad
    pad2->cd();
    
    // Create residual histogram
    TH1D* hResidual = (TH1D*)hAvg->Clone("hResidual");
    for (Int_t i = 1; i <= hResidual->GetNbinsX(); ++i) {
        double x = hResidual->GetBinCenter(i);
        double fitVal = fit->Eval(x);
        double resVal = hResidual->GetBinContent(i) - fitVal;
        hResidual->SetBinContent(i, resVal);
    }
    
    // Background for residuals
    TH1D* hbkg2 = new TH1D("hbkg2_zyam", "hbkg2_zyam", 1, -TMath::Pi()/2.0, 3*TMath::Pi()/2.0);
    hbkg2->SetStats(0);
    hbkg2->GetXaxis()->SetTitle("#Delta#phi [rad]");
    hbkg2->GetYaxis()->SetTitle("Residuals");
    hbkg2->GetXaxis()->SetTitleSize(0.12);
    hbkg2->GetXaxis()->SetLabelSize(0.10);
    hbkg2->GetYaxis()->SetTitleSize(0.12);
    hbkg2->GetYaxis()->SetLabelSize(0.10);
    hbkg2->GetXaxis()->SetTitleOffset(1.2);
    hbkg2->GetYaxis()->SetTitleOffset(0.5);
    
    double resymax = hResidual->GetMaximum();
    double resymin = hResidual->GetMinimum();
    double resydelta = resymax - resymin;
    hbkg2->GetYaxis()->SetRangeUser(resymin - resydelta*0.3, resymax + resydelta*0.3);
    hbkg2->Draw();
    
    // Plot residuals
    hResidual->SetMarkerStyle(20);
    hResidual->SetMarkerColor(kBlack);
    hResidual->SetLineColor(kBlack);
    hResidual->SetMarkerSize(1.0);
    hResidual->Draw("same p");
    
    // Zero line
    TLine* zeroLine = new TLine(-TMath::Pi()/2.0, 0, 3*TMath::Pi()/2.0, 0);
    zeroLine->SetLineStyle(2);
    zeroLine->SetLineColor(kGray+1);
    zeroLine->Draw();

    if (useDiff) {
        gSystem->mkdir("./TemplateFit/PeripheralSubtracted/EtaDiff/PDFs", kTRUE);
        c->SaveAs(Form("./TemplateFit/PeripheralSubtracted/EtaDiff/PDFs/TemplateFit_%s_%s_%d_%d_Eta_%0.1f_%0.1f.pdf",
                       fileSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax));
    } else {
        gSystem->mkdir("./TemplateFit/PeripheralSubtracted/PDFs", kTRUE);
        c->SaveAs(Form("./TemplateFit/PeripheralSubtracted/PDFs/TemplateFit_%s_%s_%d_%d.pdf",
                       fileSuffix.c_str(), splitName.c_str(), minRange, maxRange));
    }

    delete hbkg;
    delete hbkg2;
    delete hResidual;
    delete zeroLine;
    delete baselineFunc;
    delete fit;
    delete leg;
    
    // Properly clean up pads before canvas to avoid segfault during ROOT cleanup
    pad1->SetBit(kCannotPick);
    pad2->SetBit(kCannotPick);
    pad1->Close();
    pad2->Close();
    delete pad1;
    delete pad2;
    
    // Finally delete canvas
    c->SetBit(kCannotPick);
    c->Close();
    delete c;
    
    gROOT->GetListOfCanvases()->Delete();
}

VnUnit* TemplateFit_PeripheralSubtracted(Bool_t isNch, InputUnit templ, InputUnit data, Bool_t cn2Tovn2, Double_t etaMin, Double_t etaMax) {
    std::string splitName = isNch ? "Mult" : "Cent";
    // useDiff is true if eta parameters are specified (non-zero values passed)
    const Bool_t useDiff = !(etaMin == 0 && etaMax == 0);
    
    if (useDiff) {
        std::cout << "    Processing eta bin [" << etaMin << ", " << etaMax << "]" << std::endl;
    }

    // Open peripheral-subtracted data file
    TFile* datafile = nullptr;
    std::string filePath;
    if (useDiff) {
        filePath = Form("./ProcessOutput/PeripheralSubtraction/EtaDiff/BootstrapSample_%s_%s_%d_%d_Eta_%0.1f_%0.1f_PeripheralSub.root", 
                                 data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, etaMin, etaMax);
        datafile = new TFile(filePath.c_str(), "READ");
    } else {
        filePath = Form("./ProcessOutput/PeripheralSubtraction/BootstrapSample_%s_%s_%d_%d_PeripheralSub.root", 
                                 data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange);
        datafile = new TFile(filePath.c_str(), "READ");
    }

    if (!datafile || !datafile->IsOpen()) {
        std::cerr << "ERROR: Cannot open file: " << filePath << std::endl;
        if (datafile) delete datafile;
        return nullptr;
    }

    // Prepare bootstrap value arrays (like in regular template fit)
    int Nobs = 3;  // v2, v3, v4
    int NofSample = maxSample * maxSample;
    int Nbin = 1;
    
    std::vector<std::vector<std::vector<double>>> ValueArray(Nobs, std::vector<std::vector<double>>(NofSample, std::vector<double>(Nbin, 0)));
    std::vector<std::vector<std::vector<double>>> ValueErrorArray(Nobs, std::vector<std::vector<double>>(NofSample, std::vector<double>(Nbin, 0)));
    std::vector<std::vector<double>> ErrorArray(Nobs, std::vector<double>(Nbin, 0));  // 2D, not 3D!
    
    TH1D* hAvg = nullptr;
    // Store all histogram bins from all samples to compute bootstrap errors
    std::vector<std::vector<double>> hBinValues;  // [bin][sample] = value
    Int_t nSamples = 0;
    
    // Collect all bootstrap samples for averaging
    for (Int_t sample = 0; sample < NofSample; ++sample) {
        TH1D* h = (TH1D*)datafile->Get(Form("hPhiSameOverMixed_%d_%d_%d", data.minRange, data.maxRange, sample));
        if (!h) {
            ValueArray[0][sample][0] = -999;
            ValueArray[1][sample][0] = -999;
            ValueArray[2][sample][0] = -999;
            continue;
        }
        
        nSamples++;
        
        // Ensure this histogram has proper error tracking
        if (h->GetSumw2N() == 0) {
            h->Sumw2();
        }
        
        // Store bin values for bootstrap error calculation
        if (hBinValues.empty()) {
            hBinValues.resize(h->GetNbinsX());
        }
        for (Int_t b = 0; b < h->GetNbinsX(); ++b) {
            hBinValues[b].push_back(h->GetBinContent(b + 1));
        }
        
        // Accumulate average for plotting
        if (!hAvg) {
            hAvg = (TH1D*)h->Clone("hAvg_temp");
            hAvg->SetDirectory(nullptr);
            hAvg->Sumw2();
        } else {
            hAvg->Add(h);
        }
        
        // Fit this bootstrap sample
        TF1* fitSample = new TF1(Form("fitSample_%d", sample),
            "[0]*(1 + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x))",
            -TMath::Pi()/2.0, 1.5*TMath::Pi());
        
        double hmin = h->GetMinimum();
        double hmax = h->GetMaximum();
        double baselineGuess = (hmin + hmax) / 2.0;
        double v2Guess = (hmax - hmin) / (hmax + hmin) / 2.0;
        
        // CRITICAL: Ensure non-zero bin errors to avoid fitting issues
        // This is essential for proper chi-square calculation and fit convergence
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            double err = h->GetBinError(i);
            double val = h->GetBinContent(i);
            if (!(err > 0.0) || std::isnan(err) || std::isinf(err)) {
                // Set error to sqrt(content) or minimum value
                double fallback = std::sqrt(std::abs(val));
                if (!(fallback > 0.0)) fallback = 1e-6;
                h->SetBinError(i, fallback);
            }
        }
        
        fitSample->SetParameters(baselineGuess, v2Guess, 0.0, 0.0);
        fitSample->SetParLimits(0, hmin * 0.7, hmax * 1.3);
        fitSample->SetParLimits(1, -0.2, 0.2);
        fitSample->SetParLimits(2, -0.1, 0.1);
        fitSample->SetParLimits(3, -0.1, 0.1);
        fitSample->SetRange(-TMath::Pi()/2.0, 1.5*TMath::Pi());
        
        Int_t fitStatus = h->Fit(fitSample, "RQ0");
        
        // Extract parameters only if fit succeeded
        if (fitStatus == 0) {
            ValueArray[0][sample][0] = fitSample->GetParameter(1);
            ValueErrorArray[0][sample][0] = fitSample->GetParError(1);
            ValueArray[1][sample][0] = fitSample->GetParameter(2);
            ValueErrorArray[1][sample][0] = fitSample->GetParError(2);
            ValueArray[2][sample][0] = fitSample->GetParameter(3);
            ValueErrorArray[2][sample][0] = fitSample->GetParError(3);
        } else {
            ValueArray[0][sample][0] = -999;
            ValueArray[1][sample][0] = -999;
            ValueArray[2][sample][0] = -999;
        }
        
        delete fitSample;
    }

    if (!hAvg || hAvg->GetEntries() == 0) {
        std::cerr << "No valid data histograms found." << std::endl;
        datafile->Close();
        delete datafile;
        if (hAvg) delete hAvg;
        return nullptr;
    }

    if (nSamples > 0) {
        // Average the histogram normally first
        hAvg->Scale(1.0 / nSamples);
        
        // CRITICAL: Use bootstrap RMS as bin errors instead of statistical errors
        // This gives proper error estimates from the scatter of 900 samples
        for (Int_t bin = 0; bin < (Int_t)hBinValues.size(); ++bin) {
            // Calculate mean and RMS of this bin across all samples
            double sum = 0.0, sum2 = 0.0;
            int validCount = static_cast<int>(hBinValues[bin].size());
            
            for (double val : hBinValues[bin]) {
                sum += val;
                sum2 += val * val;
            }
            
            if (validCount > 1) {
                double mean = sum / validCount;
                double variance = (sum2 / validCount) - (mean * mean);
                double rms = std::sqrt(std::max(0.0, variance));  // Protect against negative variance
                // Set bin error to RMS of bootstrap samples (represents intrinsic scatter)
                hAvg->SetBinError(bin + 1, rms);
            } else {
                hAvg->SetBinError(bin + 1, 1e-6);  // Fallback for single sample
            }
        }
    }
    
    // Count successful bootstrap samples
    int successfulSamples = 0;
    for(int sample = 0; sample < NofSample; sample++) {
        if (ValueArray[0][sample][0] > -900) {
            successfulSamples++;
        }
    }
    
    if (useDiff) {
        std::cout << "    [eta " << etaMin << "," << etaMax << "] Bootstrap fits: " << successfulSamples << "/" << NofSample << " successful" << std::endl;
    }
    
    // Require minimum successful samples
    int minSamples = NofSample / 50;
    if (minSamples < 5) minSamples = 5;
    if (successfulSamples < minSamples) {
        std::cerr << "Too few successful bootstrap fits (" << successfulSamples << "/" << NofSample << ")" << std::endl;
        datafile->Close();
        delete datafile;
        delete hAvg;
        return nullptr;
    }
    
    // Calculate bootstrap errors
    for(int iobs = 0; iobs < Nobs; iobs++){
        CalculateBootstrapError(ValueArray[iobs], ValueErrorArray[iobs], ErrorArray[iobs], 1.);
    }

    // Fit the averaged histogram for plotting purposes
    // C(Δφ) = B * (1 + 2*v2*cos(2Δφ) + 2*v3*cos(3Δφ) + 2*v4*cos(4Δφ))
    
    double hmax = hAvg->GetMaximum();
    double hmin = hAvg->GetMinimum();
    double integral = hAvg->Integral();
    double mean = (integral > 0) ? integral / hAvg->GetNbinsX() : 0;
    
    // Create fit function for plotting
    TF1* fitFunc = new TF1("fitZYAM_plot",
        "[0]*(1 + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x))",
        -TMath::Pi()/2.0, 1.5*TMath::Pi());
    
    // Set starting parameters
    double baselineGuess = (hmin + hmax) / 2.0;
    double v2Guess = (hmax - hmin) / (hmax + hmin) / 2.0;
    fitFunc->SetParameters(baselineGuess, v2Guess, 0.0, 0.0);
    
    // Constrain parameters
    fitFunc->SetParLimits(0, hmin * 0.7, hmax * 1.3);
    fitFunc->SetParLimits(1, -0.2, 0.2);
    fitFunc->SetParLimits(2, -0.1, 0.1);
    fitFunc->SetParLimits(3, -0.1, 0.1);
    
    fitFunc->SetRange(-TMath::Pi()/2.0, 1.5*TMath::Pi());
    hAvg->Fit(fitFunc, "RQ0");
    
    VnUnit* vnResult = new VnUnit();
    
    // Use bootstrap-averaged parameters
    vnResult->v2 = ValueArray[0][0][0];  // Mean from CalculateBootstrapError
    vnResult->v2_err = ErrorArray[0][0];
    vnResult->v3 = ValueArray[1][0][0];
    vnResult->v3_err = ErrorArray[1][0];
    vnResult->v4 = ValueArray[2][0][0];
    vnResult->v4_err = ErrorArray[2][0];
    
    delete fitFunc;

    if (useDiff) {
        std::cout << "  Eta [" << etaMin << ", " << etaMax << "]: v2=" << vnResult->v2 << " (max/min=" << (hmax/hmin) << ")" << std::endl;
    } else {
        std::cout << "Peripheral-subtracted v2 = " << vnResult->v2 << " (averaged over " << nSamples << " samples, max/min=" << (hmax/hmin) << ")" << std::endl;
    }

    PlotPeripheralSubtractedFit(hAvg, isNch, data.fileNameSuffix, data.minRange, data.maxRange, etaMin, etaMax, vnResult, mean);

    delete hAvg;
    datafile->Close();
    delete datafile;

    return vnResult;
}
