/*
 * @Author: Martijn Laarhoven (martijn.laarhoven@cern.ch)
 * @Date: 2026-03-03
 * Uncorrected v_n extraction using Fourier fit on Mixed/Same-Event correlations
 * Extracts V1Delta, V2Delta, V3Delta, V4Delta as baseline for eta-dihadron long-range correlations
 * Generates plots for Vn vs eta
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

    TF1* fitFunc = new TF1(fitName.c_str(),
        "[0]*(1 + 2*[1]*cos(x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x))",
        -TMath::Pi()/2.0, 1.5*TMath::Pi());

    fitFunc->SetParameter(0, hist->GetMaximum());
    fitFunc->SetParameter(1, -0.001);
    fitFunc->SetParameter(2, 0.003);
    fitFunc->SetParameter(3, 0.0003);
    fitFunc->SetParameter(4, 0.00005);

    // Use chi-square fit with bin errors for proper error propagation
    // "W" option uses bin errors as weights, "E" improves error estimates
    int fitStatus = hist->Fit(fitFunc, "RWQN0E");
    if (fitStatus != 0) {
        delete fitFunc;
        return false;
    }

    fitResult.A = fitFunc->GetParameter(0);
    fitResult.v1 = fitFunc->GetParameter(1);
    fitResult.v1_err = fitFunc->GetParError(1);
    fitResult.v2 = fitFunc->GetParameter(2);
    fitResult.v2_err = fitFunc->GetParError(2);
    fitResult.v3 = fitFunc->GetParameter(3);
    fitResult.v3_err = fitFunc->GetParError(3);
    fitResult.v4 = fitFunc->GetParameter(4);
    fitResult.v4_err = fitFunc->GetParError(4);
    fitResult.chi2ndf = (fitFunc->GetNDF() > 0) ? fitFunc->GetChisquare() / fitFunc->GetNDF() : -1.0;
    fitResult.success = kTRUE;

    delete fitFunc;
    return true;
}

void PlotFourierFitStyle(TH1D* hInput,
                         const FitResult& fitResult,
                         const std::string& dataset,
                         Double_t etaMin,
                         Double_t etaMax,
                         const std::string& outPdf) {
    if (!hInput) {
        return;
    }

    TH1D* hCorr = (TH1D*)hInput->Clone(Form("hCorr_%s_%.1f_%.1f", dataset.c_str(), etaMin, etaMax));
    hCorr->SetDirectory(nullptr);

    const double xMin = -TMath::Pi()/2.0;
    const double xMax = 1.5*TMath::Pi();

    TF1* fTotal = new TF1(Form("fTotal_%s_%.1f_%.1f", dataset.c_str(), etaMin, etaMax),
        "[0]*(1 + 2*[1]*cos(x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x))", xMin, xMax);
    fTotal->SetParameters(fitResult.A, fitResult.v1, fitResult.v2, fitResult.v3, fitResult.v4);
    fTotal->SetLineColor(kRed + 1);
    fTotal->SetLineWidth(3);

    TF1* fV1 = new TF1(Form("fV1_%s_%.1f_%.1f", dataset.c_str(), etaMin, etaMax),
        "[0]*(1 + 2*[1]*cos(x))", xMin, xMax);
    TF1* fV2 = new TF1(Form("fV2_%s_%.1f_%.1f", dataset.c_str(), etaMin, etaMax),
        "[0]*(1 + 2*[1]*cos(2*x))", xMin, xMax);
    TF1* fV3 = new TF1(Form("fV3_%s_%.1f_%.1f", dataset.c_str(), etaMin, etaMax),
        "[0]*(1 + 2*[1]*cos(3*x))", xMin, xMax);
    TF1* fV4 = new TF1(Form("fV4_%s_%.1f_%.1f", dataset.c_str(), etaMin, etaMax),
        "[0]*(1 + 2*[1]*cos(4*x))", xMin, xMax);

    fV1->SetParameters(fitResult.A, fitResult.v1);
    fV2->SetParameters(fitResult.A, fitResult.v2);
    fV3->SetParameters(fitResult.A, fitResult.v3);
    fV4->SetParameters(fitResult.A, fitResult.v4);

    fV1->SetLineColor(kOrange + 7);
    fV2->SetLineColor(kBlue + 1);
    fV3->SetLineColor(kGreen + 2);
    fV4->SetLineColor(kMagenta + 1);
    fV1->SetLineStyle(2);
    fV2->SetLineStyle(2);
    fV3->SetLineStyle(2);
    fV4->SetLineStyle(2);
    fV1->SetLineWidth(2);
    fV2->SetLineWidth(2);
    fV3->SetLineWidth(2);
    fV4->SetLineWidth(2);

    TCanvas* c = new TCanvas(Form("cFit_%s_%.1f_%.1f", dataset.c_str(), etaMin, etaMax), "Uncorrected Fourier Fit", 900, 700);
    TPad* padTop = new TPad("padTop", "padTop", 0.0, 0.28, 1.0, 1.0);
    TPad* padBottom = new TPad("padBottom", "padBottom", 0.0, 0.0, 1.0, 0.28);
    padTop->SetBottomMargin(0.02);
    padTop->SetLeftMargin(0.12);
    padBottom->SetTopMargin(0.02);
    padBottom->SetBottomMargin(0.30);
    padBottom->SetLeftMargin(0.12);
    padTop->Draw();
    padBottom->Draw();

    padTop->cd();
    hCorr->SetStats(0);
    hCorr->SetTitle(Form("Uncorrected Fourier Fit %s, %.1f < #eta_{TPC} < %.1f", dataset.c_str(), etaMin, etaMax));
    hCorr->GetYaxis()->SetTitle("SE/ME");
    hCorr->GetXaxis()->SetLabelSize(0);
    hCorr->SetMarkerStyle(20);
    hCorr->SetMarkerSize(0.8);
    
    // Auto-scale Y-axis based on data range with margins
    Double_t dataMin = hCorr->GetMinimum();
    Double_t dataMax = hCorr->GetMaximum();
    Double_t dataRange = dataMax - dataMin;
    Double_t yMin = dataMin - 0.15 * dataRange;  // 15% margin below
    Double_t yMax = dataMax + 0.25 * dataRange;  // 25% margin above for legend
    hCorr->GetYaxis()->SetRangeUser(yMin, yMax);
    
    hCorr->Draw("E1");

    fTotal->Draw("same");
    fV1->Draw("same");
    fV2->Draw("same");
    fV3->Draw("same");
    fV4->Draw("same");

    TLegend* legend = new TLegend(0.52, 0.55, 0.92, 0.90);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.035);
    legend->AddEntry(hCorr, "Data (SE/ME)", "lep");
    legend->AddEntry(fTotal, "Total Fourier fit", "l");
    legend->AddEntry(fV1, Form("V_{1#Delta}=%.5f #pm %.2e", fitResult.v1, fitResult.v1_err), "l");
    legend->AddEntry(fV2, Form("V_{2#Delta}=%.5f #pm %.2e", fitResult.v2, fitResult.v2_err), "l");
    legend->AddEntry(fV3, Form("V_{3#Delta}=%.5f #pm %.2e", fitResult.v3, fitResult.v3_err), "l");
    legend->AddEntry(fV4, Form("V_{4#Delta}=%.5f #pm %.2e", fitResult.v4, fitResult.v4_err), "l");
    // Show chi2/ndf in appropriate format
    if (fitResult.chi2ndf < 0.01) {
        legend->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.2e", fitResult.chi2ndf), "");
    } else {
        legend->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.2f", fitResult.chi2ndf), "");
    }
    legend->Draw();

    padBottom->cd();
    TH1D* hRatio = (TH1D*)hCorr->Clone(Form("hRatio_%s_%.1f_%.1f", dataset.c_str(), etaMin, etaMax));
    hRatio->SetDirectory(nullptr);
    
    // Calculate ratio and find min/max for auto-scaling
    Double_t ratioMin = 1e10;
    Double_t ratioMax = -1e10;
    for (int i = 1; i <= hRatio->GetNbinsX(); ++i) {
        const double x = hRatio->GetBinCenter(i);
        const double y = hRatio->GetBinContent(i);
        const double e = hRatio->GetBinError(i);
        const double f = fTotal->Eval(x);
        if (f != 0.0 && y > 0.0) {
            double ratio = y / f;
            hRatio->SetBinContent(i, ratio);
            hRatio->SetBinError(i, e / f);
            if (ratio < ratioMin) ratioMin = ratio;
            if (ratio > ratioMax) ratioMax = ratio;
        } else {
            hRatio->SetBinContent(i, 1.0);
            hRatio->SetBinError(i, 0.0);
        }
    }
    
    // Auto-scale ratio plot - always zoom to show variation
    Double_t ratioRange = ratioMax - ratioMin;
    if (ratioRange < 0.02) {
        // Very tight fit - zoom in to show structure
        Double_t ratioMid = (ratioMin + ratioMax) / 2.0;
        ratioRange = 0.02;  // Force minimum visible range of 2%
        ratioMin = ratioMid - ratioRange / 2.0;
        ratioMax = ratioMid + ratioRange / 2.0;
    }
    Double_t ratioYmin = ratioMin - 0.4 * ratioRange;
    Double_t ratioYmax = ratioMax + 0.4 * ratioRange;

    hRatio->SetStats(0);
    hRatio->SetTitle("");
    hRatio->GetYaxis()->SetTitle("Data/Fit");
    hRatio->GetXaxis()->SetTitle("#Delta#phi [rad]");
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetTitleSize(0.10);
    hRatio->GetYaxis()->SetLabelSize(0.09);
    hRatio->GetYaxis()->SetTitleOffset(0.45);
    hRatio->GetXaxis()->SetTitleSize(0.11);
    hRatio->GetXaxis()->SetLabelSize(0.10);
    hRatio->GetXaxis()->SetTitleOffset(1.0);
    hRatio->GetYaxis()->SetRangeUser(ratioYmin, ratioYmax);
    hRatio->Draw("E1");

    TLine* unity = new TLine(xMin, 1.0, xMax, 1.0);
    unity->SetLineStyle(2);
    unity->Draw("same");

    c->SaveAs(outPdf.c_str());

    delete unity;
    delete hRatio;
    delete legend;
    delete padTop;
    delete padBottom;
    delete c;
    delete fV1;
    delete fV2;
    delete fV3;
    delete fV4;
    delete fTotal;
    delete hCorr;
}

void Process_Uncorrected_FourierFit() {
    gROOT->SetBatch(kTRUE);
    gSystem->mkdir("./TemplateFit/EtaDiff/Uncorrected", kTRUE);

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
    std::cout << "Extracting Uncorrected Vn from Correlations" << std::endl;
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
        std::string outFileName = Form("./TemplateFit/EtaDiff/Uncorrected/Vn_Uncorrected_%s_Cent_0_20.root",
                                       dataset.c_str());
        TFile* fOut = new TFile(outFileName.c_str(), "RECREATE");

        // Vectors to store Vn vs eta data
        std::vector<VnData> v1Data, v2Data, v3Data, v4Data;

        // Read Vn values from Mixed files
        for (const auto& etaBin : etaBins) {
            Double_t etaMin = etaBin.first;
            Double_t etaMax = etaBin.second;
            // Read from ProcessOutput Mixed files
            std::string inputFile = Form("./ProcessOutput/EtaDiff/Mixed_%s_Cent_0_20_Eta_%.1f_%.1f_TPC_FT0C.root",
                                        dataset.c_str(), etaMin, etaMax);

            TFile* fIn = TFile::Open(inputFile.c_str(), "READ");
            if (!fIn || fIn->IsZombie()) {
                std::cout << "  Warning: Could not open " << inputFile << std::endl;
                continue;
            }

            // Try to read the SE/ME ratio histogram
            TH1D* hCorr = nullptr;
            hCorr = (TH1D*)fIn->Get("hPhiSameOverMixed_0_20");
            
            if (!hCorr) {
                // Try alternative naming
                hCorr = (TH1D*)fIn->Get("hPhiSameOverMixed");
            }

            if (hCorr) {
                FitResult fitResult;
                std::string fitName = Form("fourierFit_%s_%.1f_%.1f", dataset.c_str(), etaMin, etaMax);
                bool fitOk = FourierFitCorrelation(hCorr, fitResult, fitName);
                if (fitOk) {
                    v1Data.push_back(VnData(etaMin, etaMax, fitResult.v1, fitResult.v1_err));
                    v2Data.push_back(VnData(etaMin, etaMax, fitResult.v2, fitResult.v2_err));
                    v3Data.push_back(VnData(etaMin, etaMax, fitResult.v3, fitResult.v3_err));
                    v4Data.push_back(VnData(etaMin, etaMax, fitResult.v4, fitResult.v4_err));

                    std::cout << "    eta [" << etaMin << ", " << etaMax << "]: "
                              << "V1 = " << fitResult.v1 << " +/- " << fitResult.v1_err << ", "
                              << "V2 = " << fitResult.v2 << " +/- " << fitResult.v2_err << ", "
                              << "V3 = " << fitResult.v3 << " +/- " << fitResult.v3_err << ", "
                              << "V4 = " << fitResult.v4 << " +/- " << fitResult.v4_err
                              << ", chi2/ndf = " << fitResult.chi2ndf
                              << std::endl;

                    std::string fitPdfName = Form("./TemplateFit/EtaDiff/Uncorrected/FourierFit_%s_Cent_0_20_Eta_%.1f_%.1f.pdf",
                                                  dataset.c_str(), etaMin, etaMax);
                    PlotFourierFitStyle(hCorr, fitResult, dataset, etaMin, etaMax, fitPdfName);
                    std::cout << "      Created fit-shape plot: " << fitPdfName << std::endl;
                } else {
                    std::cout << "  Warning: Fourier fit failed for " << inputFile << std::endl;
                }
            } else {
                std::cout << "  Warning: Could not find dPhi histogram in " << inputFile << std::endl;
            }

            fIn->Close();
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
            gV1vsEta->SetName("gV1_Uncorrected");
            gV1vsEta->SetTitle(Form("Uncorrected V1#Delta vs TPC Eta - %s;TPC #eta;V_{1}^{#Delta}", dataset.c_str()));
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
            gV2vsEta->SetName("gV2_Uncorrected");
            gV2vsEta->SetTitle(Form("Uncorrected V2#Delta vs TPC Eta - %s;TPC #eta;V_{2}^{#Delta}", dataset.c_str()));
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
            gV3vsEta->SetName("gV3_Uncorrected");
            gV3vsEta->SetTitle(Form("Uncorrected V3#Delta vs TPC Eta - %s;TPC #eta;V_{3}^{#Delta}", dataset.c_str()));
            gV3vsEta->SetMarkerStyle(22);
            gV3vsEta->SetMarkerSize(0.8);
            gV3vsEta->SetLineColor(kGreen+2);
            gV3vsEta->SetMarkerColor(kGreen+2);

            fOut->cd();
            gV3vsEta->Write();

            // Create combined plot for V1, V2, V3
            std::string pdfName = Form("./TemplateFit/EtaDiff/Uncorrected/Vn_Uncorrected_%s_Cent_0_20.pdf",
                                      dataset.c_str());
            TCanvas* c = new TCanvas("cVn", "Vn Uncorrected", 1200, 400);
            c->Divide(3, 1);

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

            c->Print(pdfName.c_str());
            std::cout << "  Created plot: " << pdfName << std::endl;

            delete c;
            delete gV1vsEta;
            delete gV2vsEta;
            delete gV3vsEta;
        }

        fOut->Close();
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "Uncorrected Fourier Fit Complete" << std::endl;
    std::cout << "Output: ./TemplateFit/EtaDiff/Uncorrected/" << std::endl;
    std::cout << "========================================\n" << std::endl;
}

void CompareUncorrectedVsTemplateFit() {
    gROOT->SetBatch(kTRUE);
    gSystem->mkdir("./TemplateFit/EtaDiff/MethodComparison", kTRUE);

    std::cout << "\n========================================" << std::endl;
    std::cout << "Comparing Uncorrected vs TemplateFit Methods" << std::endl;
    std::cout << "========================================\n" << std::endl;

    std::vector<std::string> datasets = {
        "LHC25af_pass2_615818",  // Ne-Ne outer
        "LHC25af_pass2_615817",  // Ne-Ne inner
        "LHC25ae_pass2_616549",  // O-O outer
        "LHC25ae_pass2_618685",  // O-O inner
        "LHC25af_pass2_617826",  // Ne-Ne inner variant
        "LHC25af_pass2_617910"   // Ne-Ne outer variant
    };

    std::vector<std::pair<Double_t, Double_t>> etaBins = {
        {-0.8, -0.7},
        {0.7, 0.8}
    };

    for (const auto& dataset : datasets) {
        std::cout << "Processing dataset: " << dataset << std::endl;

        // Storage for both eta bins
        std::vector<Double_t> etaPoints;
        std::vector<Double_t> v2_uncorr, v2_uncorr_err, v2_template, v2_template_err;
        std::vector<Double_t> v3_uncorr, v3_uncorr_err, v3_template, v3_template_err;
        std::vector<Double_t> v4_uncorr, v4_uncorr_err, v4_template, v4_template_err;

        for (const auto& [etaMin, etaMax] : etaBins) {
            Double_t etaMid = (etaMin + etaMax) / 2.0;
            etaPoints.push_back(etaMid);

            // Read Uncorrected FourierFit results
            std::string uncorrFile = Form("./TemplateFit/EtaDiff/Uncorrected/Vn_Uncorrected_%s_Cent_0_20.root", dataset.c_str());
            TFile* fUncorr = TFile::Open(uncorrFile.c_str(), "READ");
            if (fUncorr && fUncorr->IsOpen()) {
                TGraphErrors* gV2 = (TGraphErrors*)fUncorr->Get("gV2_Uncorrected");
                TGraphErrors* gV3 = (TGraphErrors*)fUncorr->Get("gV3_Uncorrected");
                TGraphErrors* gV4 = (TGraphErrors*)fUncorr->Get("gV4_Uncorrected");
                
                if (gV2 && gV3 && gV4) {
                    // Find point matching this eta
                    for (int i = 0; i < gV2->GetN(); i++) {
                        Double_t x, y;
                        gV2->GetPoint(i, x, y);
                        if (TMath::Abs(x - etaMid) < 0.06) {
                            v2_uncorr.push_back(y);
                            v2_uncorr_err.push_back(gV2->GetErrorY(i));
                            gV3->GetPoint(i, x, y);
                            v3_uncorr.push_back(y);
                            v3_uncorr_err.push_back(gV3->GetErrorY(i));
                            gV4->GetPoint(i, x, y);
                            v4_uncorr.push_back(y);
                            v4_uncorr_err.push_back(gV4->GetErrorY(i));
                            break;
                        }
                    }
                }
                fUncorr->Close();
            } else {
                std::cout << "  Warning: Cannot open " << uncorrFile << std::endl;
                v2_uncorr.push_back(0); v2_uncorr_err.push_back(0);
                v3_uncorr.push_back(0); v3_uncorr_err.push_back(0);
                v4_uncorr.push_back(0); v4_uncorr_err.push_back(0);
            }

            // Read TemplateFit results
            std::string suffix = (etaMin < 0) ? Form("_%.1f_%.1f", etaMin, etaMax) : Form("_%.1f_%.1f", etaMin, etaMax);
            std::string templateFile = Form("./TemplateFit/EtaDiff/Vn_%s_TPC_FT0C_TPCEta%s_Cent_0_20.root", 
                                           dataset.c_str(), suffix.c_str());
            TFile* fTemplate = TFile::Open(templateFile.c_str(), "READ");
            if (fTemplate && fTemplate->IsOpen()) {
                TH1D* hV2 = (TH1D*)fTemplate->Get("hV2");
                TH1D* hV3 = (TH1D*)fTemplate->Get("hV3");
                TH1D* hV4 = (TH1D*)fTemplate->Get("hV4");
                
                if (hV2 && hV3 && hV4) {
                    v2_template.push_back(hV2->GetBinContent(1));
                    v2_template_err.push_back(hV2->GetBinError(1));
                    v3_template.push_back(hV3->GetBinContent(1));
                    v3_template_err.push_back(hV3->GetBinError(1));
                    v4_template.push_back(hV4->GetBinContent(1));
                    v4_template_err.push_back(hV4->GetBinError(1));
                }
                fTemplate->Close();
            } else {
                std::cout << "  Warning: Cannot open " << templateFile << std::endl;
                v2_template.push_back(0); v2_template_err.push_back(0);
                v3_template.push_back(0); v3_template_err.push_back(0);
                v4_template.push_back(0); v4_template_err.push_back(0);
            }
        }

        // Create comparison graphs
        TGraphErrors* gV2_uncorr = new TGraphErrors(etaPoints.size());
        TGraphErrors* gV2_template = new TGraphErrors(etaPoints.size());
        TGraphErrors* gV3_uncorr = new TGraphErrors(etaPoints.size());
        TGraphErrors* gV3_template = new TGraphErrors(etaPoints.size());
        TGraphErrors* gV4_uncorr = new TGraphErrors(etaPoints.size());
        TGraphErrors* gV4_template = new TGraphErrors(etaPoints.size());

        for (size_t i = 0; i < etaPoints.size(); i++) {
            gV2_uncorr->SetPoint(i, etaPoints[i], v2_uncorr[i]);
            gV2_uncorr->SetPointError(i, 0.05, v2_uncorr_err[i]);
            gV2_template->SetPoint(i, etaPoints[i], v2_template[i]);
            gV2_template->SetPointError(i, 0.05, v2_template_err[i]);

            gV3_uncorr->SetPoint(i, etaPoints[i], v3_uncorr[i]);
            gV3_uncorr->SetPointError(i, 0.05, v3_uncorr_err[i]);
            gV3_template->SetPoint(i, etaPoints[i], v3_template[i]);
            gV3_template->SetPointError(i, 0.05, v3_template_err[i]);

            gV4_uncorr->SetPoint(i, etaPoints[i], v4_uncorr[i]);
            gV4_uncorr->SetPointError(i, 0.05, v4_uncorr_err[i]);
            gV4_template->SetPoint(i, etaPoints[i], v4_template[i]);
            gV4_template->SetPointError(i, 0.05, v4_template_err[i]);
        }

        // Style settings
        gV2_uncorr->SetMarkerStyle(20); gV2_uncorr->SetMarkerColor(kBlue); gV2_uncorr->SetLineColor(kBlue);
        gV2_template->SetMarkerStyle(24); gV2_template->SetMarkerColor(kRed+1); gV2_template->SetLineColor(kRed+1);
        gV3_uncorr->SetMarkerStyle(21); gV3_uncorr->SetMarkerColor(kGreen+2); gV3_uncorr->SetLineColor(kGreen+2);
        gV3_template->SetMarkerStyle(25); gV3_template->SetMarkerColor(kMagenta+1); gV3_template->SetLineColor(kMagenta+1);
        gV4_uncorr->SetMarkerStyle(22); gV4_uncorr->SetMarkerColor(kOrange+7); gV4_uncorr->SetLineColor(kOrange+7);
        gV4_template->SetMarkerStyle(26); gV4_template->SetMarkerColor(kCyan+2); gV4_template->SetLineColor(kCyan+2);

        // Create comparison plot
        TCanvas* c = new TCanvas("cComp", "Method Comparison", 1400, 500);
        c->Divide(3, 1);

        // V2 panel
        c->cd(1);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gV2_uncorr->SetTitle(Form("V_{2}^{#Delta} Comparison - %s;TPC #eta;V_{2}^{#Delta}", dataset.c_str()));
        gV2_uncorr->GetYaxis()->SetRangeUser(0.0, 0.004);
        gV2_uncorr->GetXaxis()->SetLimits(-1.0, 1.0);
        gV2_uncorr->Draw("AP");
        gV2_template->Draw("P SAME");
        TLegend* leg1 = new TLegend(0.18, 0.75, 0.55, 0.90);
        leg1->SetBorderSize(0); leg1->SetFillStyle(0);
        leg1->AddEntry(gV2_uncorr, "Uncorrected FourierFit", "lep");
        leg1->AddEntry(gV2_template, "TemplateFit", "lep");
        leg1->Draw();

        // V3 panel
        c->cd(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gV3_uncorr->SetTitle(Form("V_{3}^{#Delta} Comparison - %s;TPC #eta;V_{3}^{#Delta}", dataset.c_str()));
        gV3_uncorr->GetYaxis()->SetRangeUser(0.0, 0.0006);
        gV3_uncorr->GetXaxis()->SetLimits(-1.0, 1.0);
        gV3_uncorr->Draw("AP");
        gV3_template->Draw("P SAME");
        TLegend* leg2 = new TLegend(0.18, 0.75, 0.55, 0.90);
        leg2->SetBorderSize(0); leg2->SetFillStyle(0);
        leg2->AddEntry(gV3_uncorr, "Uncorrected FourierFit", "lep");
        leg2->AddEntry(gV3_template, "TemplateFit", "lep");
        leg2->Draw();

        // V4 panel
        c->cd(3);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gV4_uncorr->SetTitle(Form("V_{4}^{#Delta} Comparison - %s;TPC #eta;V_{4}^{#Delta}", dataset.c_str()));
        gV4_uncorr->GetYaxis()->SetRangeUser(-0.00005, 0.00015);
        gV4_uncorr->GetXaxis()->SetLimits(-1.0, 1.0);
        gV4_uncorr->Draw("AP");
        gV4_template->Draw("P SAME");
        TLegend* leg3 = new TLegend(0.18, 0.75, 0.55, 0.90);
        leg3->SetBorderSize(0); leg3->SetFillStyle(0);
        leg3->AddEntry(gV4_uncorr, "Uncorrected FourierFit", "lep");
        leg3->AddEntry(gV4_template, "TemplateFit", "lep");
        leg3->Draw();

        std::string pdfName = Form("./TemplateFit/EtaDiff/MethodComparison/VnComparison_%s_Cent_0_20.pdf", dataset.c_str());
        c->Print(pdfName.c_str());
        std::cout << "  Created comparison: " << pdfName << std::endl;

        delete leg1; delete leg2; delete leg3;
        delete c;
        delete gV2_uncorr; delete gV2_template;
        delete gV3_uncorr; delete gV3_template;
        delete gV4_uncorr; delete gV4_template;
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "Method Comparison Complete" << std::endl;
    std::cout << "Output: ./TemplateFit/EtaDiff/MethodComparison/" << std::endl;
    std::cout << "========================================\n" << std::endl;
}
