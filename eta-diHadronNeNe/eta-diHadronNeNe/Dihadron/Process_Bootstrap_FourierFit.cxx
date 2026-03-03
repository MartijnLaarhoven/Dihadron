/*
 * @Author: Martijn Laarhoven (martijn.laarhoven@cern.ch) with AI assistance
 * @Date: 2026-03-03
 * Bootstrap Fourier fit for v_n extraction with proper uncertainty quantification
 * Applies Fourier fit to bootstrap resampled data for robust error estimation
 */

#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLine.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "./include/BasicForDihadron.h"

const int numBootstrap = 100; // Number of bootstrap samples

struct VnData {
    Double_t etaMin, etaMax;
    Double_t vn, vn_err;
    VnData(Double_t etaMin, Double_t etaMax, Double_t vn, Double_t vn_err) :
        etaMin(etaMin), etaMax(etaMax), vn(vn), vn_err(vn_err) {}
};

struct FitResult {
    Double_t A;
    Double_t v1;
    Double_t v1_err;
    Double_t v2;
    Double_t v2_err;
    Double_t v3;
    Double_t v3_err;
    Double_t v4;
    Double_t v4_err;
    Double_t chi2ndf;
    Bool_t success;
};

bool FourierFitCorrelation(TH1D* hist, FitResult& fitResult, const std::string& fitName) {
    fitResult = {0, 0, 0.01, 0, 0.01, 0, 0.01, 0, 0.01, -1.0, kFALSE};
    if (!hist || hist->GetEntries() == 0) {
        return false;
    }

    // Clone histogram to avoid modifying original
    TH1D* hFit = (TH1D*)hist->Clone(Form("hFit_%s", fitName.c_str()));
    hFit->SetDirectory(nullptr);

    // Add systematic uncertainty to histogram errors (mirrors TemplateFit approach)
    // Systematic fraction is applied in two places: histogram errors AND chi2 calculation
    const double systematic_fraction = 0.0002; // 0.02% systematic uncertainty (same as TemplateFit)
    for (int i = 1; i <= hFit->GetNbinsX(); i++) {
        double content = hFit->GetBinContent(i);
        double stat_error = hFit->GetBinError(i);
        double syst_error = systematic_fraction * TMath::Abs(content);
        double total_error = TMath::Sqrt(stat_error * stat_error + syst_error * syst_error);
        hFit->SetBinError(i, total_error);
    }

    TF1* fitFunc = new TF1(fitName.c_str(),
        "[0]*(1 + 2*[1]*cos(x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x))",
        -TMath::Pi()/2.0, 1.5*TMath::Pi());

    fitFunc->SetParameter(0, hFit->GetMaximum());
    fitFunc->SetParameter(1, -0.001);
    fitFunc->SetParameter(2, 0.003);
    fitFunc->SetParameter(3, 0.0003);
    fitFunc->SetParameter(4, 0.00005);

    // Use chi-square fit with bin errors for proper error propagation
    int fitStatus = hFit->Fit(fitFunc, "RWQN0E");
    if (fitStatus != 0) {
        delete fitFunc;
        delete hFit;
        return false;
    }

    // Manually calculate chi2 using the histogram with inflated errors
    // Matches TemplateFit's findChi2ndf method exactly
    double chi2 = 0.0;
    int nBins = hFit->GetNbinsX();
    int nParams = 5; // A, v1, v2, v3, v4
    int ndf = nBins - nParams;
    
    for (int i = 1; i <= nBins; i++) {
        double content = hFit->GetBinContent(i);
        double stat_error = hFit->GetBinError(i);
        if (stat_error > 0) {
            double xCenter = hFit->GetBinCenter(i);
            double fitted = fitFunc->Eval(xCenter);
            
            // Apply systematic uncertainty again in chi2 calculation (same as TemplateFit)
            double syst_error = systematic_fraction * TMath::Abs(content);
            double error_total = TMath::Sqrt(stat_error * stat_error + syst_error * syst_error);
            
            double residual = content - fitted;
            chi2 += (residual / error_total) * (residual / error_total);
        }
    }
    if (ndf <= 0) ndf = 1;
    
    double chi2ndf_manual = chi2 / ndf;

    fitResult.A = fitFunc->GetParameter(0);
    fitResult.v1 = fitFunc->GetParameter(1);
    fitResult.v1_err = fitFunc->GetParError(1);
    fitResult.v2 = fitFunc->GetParameter(2);
    fitResult.v2_err = fitFunc->GetParError(2);
    fitResult.v3 = fitFunc->GetParameter(3);
    fitResult.v3_err = fitFunc->GetParError(3);
    fitResult.v4 = fitFunc->GetParameter(4);
    fitResult.v4_err = fitFunc->GetParError(4);
    fitResult.chi2ndf = chi2ndf_manual;
    fitResult.success = kTRUE;

    delete fitFunc;
    delete hFit;
    return true;
}

void Process_Bootstrap_FourierFit() {
    gROOT->SetBatch(kTRUE);
    gSystem->mkdir("./TemplateFit/EtaDiff/Bootstrap_FourierFit", kTRUE);

    // Define ring datasets for eta-diHadronNeNe long-range correlations
    std::vector<std::string> datasets = {
        "LHC25af_pass2_615818",  // Ne-Ne outer
        "LHC25af_pass2_615817",  // Ne-Ne inner
        "LHC25ae_pass2_616549",  // O-O outer
        "LHC25ae_pass2_618685",  // O-O inner
        "LHC25af_pass2_617826",  // Ne-Ne inner variant
        "LHC25af_pass2_617910"   // Ne-Ne outer variant
    };

    std::cout << "\n========================================" << std::endl;
    std::cout << "Bootstrap Fourier Fit for Vn Extraction" << std::endl;
    std::cout << "========================================\n" << std::endl;

    // Process each dataset
    for (const auto& dataset : datasets) {
        std::cout << "Processing dataset: " << dataset << std::endl;

        // Define eta bins
        std::vector<std::pair<Double_t, Double_t>> etaBins = {
            {-0.8, -0.7},
            {0.7, 0.8}
        };

        // Create output file for this dataset
        std::string outFileName = Form("./TemplateFit/EtaDiff/Bootstrap_FourierFit/Vn_Bootstrap_FourierFit_%s_Cent_0_20.root",
                                       dataset.c_str());
        TFile* fOut = new TFile(outFileName.c_str(), "RECREATE");

        // Vectors to store Vn vs eta data
        std::vector<VnData> v1Data, v2Data, v3Data, v4Data;

        // Process each eta bin
        for (const auto& etaBin : etaBins) {
            Double_t etaMin = etaBin.first;
            Double_t etaMax = etaBin.second;
            
            // Read bootstrap samples from ProcessOutput
            std::string inputFile = Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_Cent_0_20_Eta_%.1f_%.1f.root",
                                        dataset.c_str(), etaMin, etaMax);

            TFile* fIn = TFile::Open(inputFile.c_str(), "READ");
            if (!fIn || fIn->IsZombie()) {
                std::cout << "  Warning: Could not open " << inputFile << std::endl;
                continue;
            }

            // Fit each bootstrap sample
            std::vector<Double_t> v1vals, v2vals, v3vals, v4vals;
            int successCount = 0;

            for (int iBootstrap = 0; iBootstrap < numBootstrap; iBootstrap++) {
                TH1D* hCorr = (TH1D*)fIn->Get(Form("bsSample_hPhiSameOverMixed_0_20_%d", iBootstrap));
                
                if (!hCorr) {
                    continue;
                }

                FitResult fitResult;
                std::string fitName = Form("bsFit_%s_%.1f_%.1f_%d", dataset.c_str(), etaMin, etaMax, iBootstrap);
                bool fitOk = FourierFitCorrelation(hCorr, fitResult, fitName);
                
                if (fitOk) {
                    v1vals.push_back(fitResult.v1);
                    v2vals.push_back(fitResult.v2);
                    v3vals.push_back(fitResult.v3);
                    v4vals.push_back(fitResult.v4);
                    successCount++;
                }
            }

            fIn->Close();

            if (successCount > 0) {
                // Calculate mean and standard deviation
                auto calcStats = [](const std::vector<Double_t>& vals) -> std::pair<Double_t, Double_t> {
                    Double_t sum = 0.0, sum2 = 0.0;
                    for (auto v : vals) {
                        sum += v;
                        sum2 += v * v;
                    }
                    Double_t mean = sum / vals.size();
                    Double_t variance = (sum2 / vals.size()) - (mean * mean);
                    Double_t stddev = (variance > 0) ? TMath::Sqrt(variance) : 0.0;
                    return {mean, stddev};
                };

                auto v1Stats = calcStats(v1vals);
                auto v2Stats = calcStats(v2vals);
                auto v3Stats = calcStats(v3vals);
                auto v4Stats = calcStats(v4vals);

                v1Data.push_back(VnData(etaMin, etaMax, v1Stats.first, v1Stats.second));
                v2Data.push_back(VnData(etaMin, etaMax, v2Stats.first, v2Stats.second));
                v3Data.push_back(VnData(etaMin, etaMax, v3Stats.first, v3Stats.second));
                v4Data.push_back(VnData(etaMin, etaMax, v4Stats.first, v4Stats.second));

                std::cout << "    eta [" << etaMin << ", " << etaMax << "]: "
                          << successCount << "/" << numBootstrap << " fits successful" << std::endl;
                std::cout << "      V1 = " << v1Stats.first << " +/- " << v1Stats.second << std::endl;
                std::cout << "      V2 = " << v2Stats.first << " +/- " << v2Stats.second << std::endl;
                std::cout << "      V3 = " << v3Stats.first << " +/- " << v3Stats.second << std::endl;
                std::cout << "      V4 = " << v4Stats.first << " +/- " << v4Stats.second << std::endl;
            } else {
                std::cout << "  Warning: All bootstrap fits failed for eta [" << etaMin << ", " << etaMax << "]" << std::endl;
            }
        }

        // Create and write Vn vs eta graphs
        if (!v1Data.empty()) {
            // V1
            TGraphErrors* gV1vsEta = new TGraphErrors(v1Data.size());
            for (size_t i = 0; i < v1Data.size(); i++) {
                Double_t etaMid = (v1Data[i].etaMin + v1Data[i].etaMax) / 2.0;
                gV1vsEta->SetPoint(i, etaMid, v1Data[i].vn);
                gV1vsEta->SetPointError(i, 0.05, v1Data[i].vn_err);
            }
            gV1vsEta->SetName("gV1_Bootstrap_FourierFit");
            gV1vsEta->SetTitle(Form("Bootstrap FourierFit V1#Delta vs TPC Eta - %s;TPC #eta;V_{1}^{#Delta}", dataset.c_str()));
            gV1vsEta->SetMarkerStyle(20);
            gV1vsEta->SetMarkerSize(0.8);
            gV1vsEta->SetLineColor(kRed);
            gV1vsEta->SetMarkerColor(kRed);

            fOut->cd();
            gV1vsEta->Write();

            // V2
            TGraphErrors* gV2vsEta = new TGraphErrors(v2Data.size());
            for (size_t i = 0; i < v2Data.size(); i++) {
                Double_t etaMid = (v2Data[i].etaMin + v2Data[i].etaMax) / 2.0;
                gV2vsEta->SetPoint(i, etaMid, v2Data[i].vn);
                gV2vsEta->SetPointError(i, 0.05, v2Data[i].vn_err);
            }
            gV2vsEta->SetName("gV2_Bootstrap_FourierFit");
            gV2vsEta->SetTitle(Form("Bootstrap FourierFit V2#Delta vs TPC Eta - %s;TPC #eta;V_{2}^{#Delta}", dataset.c_str()));
            gV2vsEta->SetMarkerStyle(21);
            gV2vsEta->SetMarkerSize(0.8);
            gV2vsEta->SetLineColor(kBlue);
            gV2vsEta->SetMarkerColor(kBlue);

            fOut->cd();
            gV2vsEta->Write();

            // V3
            TGraphErrors* gV3vsEta = new TGraphErrors(v3Data.size());
            for (size_t i = 0; i < v3Data.size(); i++) {
                Double_t etaMid = (v3Data[i].etaMin + v3Data[i].etaMax) / 2.0;
                gV3vsEta->SetPoint(i, etaMid, v3Data[i].vn);
                gV3vsEta->SetPointError(i, 0.05, v3Data[i].vn_err);
            }
            gV3vsEta->SetName("gV3_Bootstrap_FourierFit");
            gV3vsEta->SetTitle(Form("Bootstrap FourierFit V3#Delta vs TPC Eta - %s;TPC #eta;V_{3}^{#Delta}", dataset.c_str()));
            gV3vsEta->SetMarkerStyle(22);
            gV3vsEta->SetMarkerSize(0.8);
            gV3vsEta->SetLineColor(kGreen+2);
            gV3vsEta->SetMarkerColor(kGreen+2);

            fOut->cd();
            gV3vsEta->Write();

            // V4
            TGraphErrors* gV4vsEta = new TGraphErrors(v4Data.size());
            for (size_t i = 0; i < v4Data.size(); i++) {
                Double_t etaMid = (v4Data[i].etaMin + v4Data[i].etaMax) / 2.0;
                gV4vsEta->SetPoint(i, etaMid, v4Data[i].vn);
                gV4vsEta->SetPointError(i, 0.05, v4Data[i].vn_err);
            }
            gV4vsEta->SetName("gV4_Bootstrap_FourierFit");
            gV4vsEta->SetTitle(Form("Bootstrap FourierFit V4#Delta vs TPC Eta - %s;TPC #eta;V_{4}^{#Delta}", dataset.c_str()));
            gV4vsEta->SetMarkerStyle(23);
            gV4vsEta->SetMarkerSize(0.8);
            gV4vsEta->SetLineColor(kMagenta);
            gV4vsEta->SetMarkerColor(kMagenta);

            fOut->cd();
            gV4vsEta->Write();

            // Create combined plot for V1, V2, V3, V4
            std::string pdfName = Form("./TemplateFit/EtaDiff/Bootstrap_FourierFit/Vn_Bootstrap_FourierFit_%s_Cent_0_20.pdf",
                                      dataset.c_str());
            TCanvas* c = new TCanvas("cVn", "Vn Bootstrap FourierFit", 1600, 400);
            c->Divide(4, 1);

            c->cd(1);
            gPad->SetLeftMargin(0.15);
            gV1vsEta->Draw("AP");
            gV1vsEta->GetYaxis()->SetRangeUser(-0.01, 0.01);
            gV1vsEta->GetXaxis()->SetLimits(-1.0, 1.0);

            c->cd(2);
            gPad->SetLeftMargin(0.15);
            gV2vsEta->Draw("AP");
            gV2vsEta->GetYaxis()->SetRangeUser(-0.005, 0.015);
            gV2vsEta->GetXaxis()->SetLimits(-1.0, 1.0);

            c->cd(3);
            gPad->SetLeftMargin(0.15);
            gV3vsEta->Draw("AP");
            gV3vsEta->GetYaxis()->SetRangeUser(-0.005, 0.005);
            gV3vsEta->GetXaxis()->SetLimits(-1.0, 1.0);

            c->cd(4);
            gPad->SetLeftMargin(0.15);
            gV4vsEta->Draw("AP");
            gV4vsEta->GetYaxis()->SetRangeUser(-0.002, 0.002);
            gV4vsEta->GetXaxis()->SetLimits(-1.0, 1.0);

            c->Print(pdfName.c_str());
            std::cout << "  Created plot: " << pdfName << std::endl;

            delete c;
            delete gV1vsEta;
            delete gV2vsEta;
            delete gV3vsEta;
            delete gV4vsEta;
        }

        fOut->Close();
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "Bootstrap Fourier Fit Complete" << std::endl;
    std::cout << "Bootstrap fits: " << numBootstrap << "/sample" << std::endl;
    std::cout << "Output: ./TemplateFit/EtaDiff/Bootstrap_FourierFit/" << std::endl;
    std::cout << "========================================\n" << std::endl;
}
