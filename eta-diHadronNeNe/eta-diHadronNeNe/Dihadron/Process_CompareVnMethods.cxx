/*
 * @Author: Martijn Laarhoven (martijn.laarhoven@cern.ch)
 * @Date: 2026-03-03
 * Compare V2Delta, V3Delta, V4Delta from Uncorrected FourierFit vs TemplateFit methods
 * Creates side-by-side comparison plots for each dataset
 */

#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TMath.h"
#include <iostream>
#include <string>
#include <vector>
#include "./include/BasicForDihadron.h"

void Process_CompareVnMethods() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gSystem->mkdir("./TemplateFit/EtaDiff/MethodComparison", kTRUE);

    std::cout << "\n========================================" << std::endl;
    std::cout << "Comparing V2Δ, V3Δ, V4Δ Methods" << std::endl;
    std::cout << "========================================\n" << std::endl;

    std::vector<std::string> datasets = {
        "LHC25af_pass2_615818",  // Ne-Ne outer
        "LHC25af_pass2_615817",  // Ne-Ne inner
        "LHC25ae_pass2_616549",  // O-O outer
        "LHC25ae_pass2_618685",  // O-O inner
        "LHC25af_pass2_617826",  // Ne-Ne inner variant
        "LHC25af_pass2_617910"   // Ne-Ne outer variant
    };

    std::vector<std::string> datasetLabels = {
        "Ne-Ne (Outer)",
        "Ne-Ne (Inner)",
        "O-O (Outer)",
        "O-O (Inner)",
        "Ne-Ne (Inner-2)",
        "Ne-Ne (Outer-2)"
    };

    std::vector<std::pair<Double_t, Double_t>> etaBins = {
        {-0.8, -0.7},
        {0.7, 0.8}
    };

    for (size_t iDataset = 0; iDataset < datasets.size(); iDataset++) {
        const auto& dataset = datasets[iDataset];
        std::cout << "Processing: " << dataset << " (" << datasetLabels[iDataset] << ")" << std::endl;

        // Storage for both eta bins
        std::vector<Double_t> etaPoints;
        std::vector<Double_t> v2_uncorr, v2_uncorr_err, v2_template, v2_template_err;
        std::vector<Double_t> v3_uncorr, v3_uncorr_err, v3_template, v3_template_err;
        std::vector<Double_t> v4_uncorr, v4_uncorr_err, v4_template, v4_template_err;

        for (const auto& [etaMin, etaMax] : etaBins) {
            Double_t etaMid = (etaMin + etaMax) / 2.0;
            etaPoints.push_back(etaMid);

            // Read Uncorrected FourierFit results
            TString uncorrFile = Form("./TemplateFit/EtaDiff/Uncorrected/Vn_Uncorrected_%s_Cent_0_20.root", dataset.c_str());
            TFile* fUncorr = TFile::Open(uncorrFile.Data(), "READ");
            bool foundUncorr = false;
            
            if (fUncorr && !fUncorr->IsZombie()) {
                TGraphErrors* gV2 = (TGraphErrors*)fUncorr->Get("gV2_Uncorrected");
                TGraphErrors* gV3 = (TGraphErrors*)fUncorr->Get("gV3_Uncorrected");
                TGraphErrors* gV4 = (TGraphErrors*)fUncorr->Get("gV4_Uncorrected");  // May not exist
                
                if (gV2 && gV3) {  // V4 is optional
                    // Find point matching this eta (check all points)
                    for (int i = 0; i < gV2->GetN(); i++) {
                        Double_t x, y;
                        gV2->GetPoint(i, x, y);
                        // Printf("  Checking: file eta=%.4f vs target eta=%.4f (diff=%.4f)", x, etaMid, TMath::Abs(x - etaMid));
                        if (TMath::Abs(x - etaMid) < 0.01) {  // Tighter tolerance
                            v2_uncorr.push_back(y);
                            v2_uncorr_err.push_back(gV2->GetErrorY(i));
                            gV3->GetPoint(i, x, y);
                            v3_uncorr.push_back(y);
                            v3_uncorr_err.push_back(gV3->GetErrorY(i));
                            if (gV4) {
                                gV4->GetPoint(i, x, y);
                                v4_uncorr.push_back(y);
                                v4_uncorr_err.push_back(gV4->GetErrorY(i));
                            } else {
                                v4_uncorr.push_back(0);
                                v4_uncorr_err.push_back(0);
                            }
                            foundUncorr = true;
                            break;
                        }
                    }
                } else {
                    Printf("  Error: Cannot get V2/V3 TGraphErrors from file!");
                }
                fUncorr->Close();
            }
            
            if (!foundUncorr) {
                Printf("  Warning: Cannot read Uncorrected data for eta=%.2f from %s", etaMid, uncorrFile.Data());
                v2_uncorr.push_back(0); v2_uncorr_err.push_back(0);
                v3_uncorr.push_back(0); v3_uncorr_err.push_back(0);
                v4_uncorr.push_back(0); v4_uncorr_err.push_back(0);
            }

            // Read TemplateFit results
            TString suffix = (etaMin < 0) ? Form("_%.1f_%.1f", etaMin, etaMax) : Form("_%.1f_%.1f", etaMin, etaMax);
            TString templateFile = Form("./TemplateFit/EtaDiff/Vn_%s_TPC_FT0C_TPCEta%s_Cent_0_20.root", 
                                       dataset.c_str(), suffix.Data());
            TFile* fTemplate = TFile::Open(templateFile.Data(), "READ");
            bool foundTemplate = false;
            
            if (fTemplate && !fTemplate->IsZombie()) {
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
                    foundTemplate = true;
                }
                fTemplate->Close();
            }
            
            if (!foundTemplate) {
                Printf("  Warning: Cannot read TemplateFit data from %s", templateFile.Data());
                v2_template.push_back(0); v2_template_err.push_back(0);
                v3_template.push_back(0); v3_template_err.push_back(0);
                v4_template.push_back(0); v4_template_err.push_back(0);
            }
        }

        // Print comparison values
        if (v2_uncorr.size() == etaPoints.size() && v2_template.size() == etaPoints.size()) {
            std::cout << "  V2Δ comparison:" << std::endl;
            for (size_t i = 0; i < etaPoints.size(); i++) {
                Printf("    η=%.2f: Uncorr=%.5f±%.2e, Template=%.5f±%.2e (diff=%.5f)", 
                       etaPoints[i], v2_uncorr[i], v2_uncorr_err[i],
                       v2_template[i], v2_template_err[i], 
                       (v2_uncorr[i] - v2_template[i]));
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
            gV2_template->SetPoint(i, etaPoints[i] + 0.02, v2_template[i]);  // Slight offset for visibility
            gV2_template->SetPointError(i, 0.05, v2_template_err[i]);

            gV3_uncorr->SetPoint(i, etaPoints[i], v3_uncorr[i]);
            gV3_uncorr->SetPointError(i, 0.05, v3_uncorr_err[i]);
            gV3_template->SetPoint(i, etaPoints[i] + 0.02, v3_template[i]);
            gV3_template->SetPointError(i, 0.05, v3_template_err[i]);

            gV4_uncorr->SetPoint(i, etaPoints[i], v4_uncorr[i]);
            gV4_uncorr->SetPointError(i, 0.05, v4_uncorr_err[i]);
            gV4_template->SetPoint(i, etaPoints[i] + 0.02, v4_template[i]);
            gV4_template->SetPointError(i, 0.05, v4_template_err[i]);
        }

        // Style settings
        gV2_uncorr->SetMarkerStyle(20); gV2_uncorr->SetMarkerColor(kBlue); gV2_uncorr->SetLineColor(kBlue);
        gV2_uncorr->SetMarkerSize(1.2);
        gV2_template->SetMarkerStyle(24); gV2_template->SetMarkerColor(kRed+1); gV2_template->SetLineColor(kRed+1);
        gV2_template->SetMarkerSize(1.2);
        
        gV3_uncorr->SetMarkerStyle(21); gV3_uncorr->SetMarkerColor(kBlue); gV3_uncorr->SetLineColor(kBlue);
        gV3_uncorr->SetMarkerSize(1.2);
        gV3_template->SetMarkerStyle(25); gV3_template->SetMarkerColor(kRed+1); gV3_template->SetLineColor(kRed+1);
        gV3_template->SetMarkerSize(1.2);
        
        gV4_uncorr->SetMarkerStyle(22); gV4_uncorr->SetMarkerColor(kBlue); gV4_uncorr->SetLineColor(kBlue);
        gV4_uncorr->SetMarkerSize(1.2);
        gV4_template->SetMarkerStyle(26); gV4_template->SetMarkerColor(kRed+1); gV4_template->SetLineColor(kRed+1);
        gV4_template->SetMarkerSize(1.2);

        // Create comparison plot
        TCanvas* c = new TCanvas("cComp", "Method Comparison", 1600, 600);
        c->Divide(3, 1);

        // V2 panel
        c->cd(1);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.10);
        gPad->SetBottomMargin(0.12);
        gPad->SetGridy();
        
        gV2_uncorr->SetTitle("");
        gV2_uncorr->GetYaxis()->SetTitle("V_{2}^{#Delta}");
        gV2_uncorr->GetXaxis()->SetTitle("TPC #eta");
        gV2_uncorr->GetYaxis()->SetTitleSize(0.05);
        gV2_uncorr->GetXaxis()->SetTitleSize(0.05);
        gV2_uncorr->GetYaxis()->SetRangeUser(0.0015, 0.0035);
        gV2_uncorr->GetXaxis()->SetLimits(-1.0, 1.0);
        gV2_uncorr->Draw("AP");
        gV2_template->Draw("P SAME");
        
        TLatex* lat1 = new TLatex();
        lat1->SetNDC();
        lat1->SetTextSize(0.045);
        lat1->DrawLatex(0.18, 0.92, Form("%s", datasetLabels[iDataset].c_str()));
        
        TLegend* leg1 = new TLegend(0.18, 0.70, 0.65, 0.88);
        leg1->SetBorderSize(0); leg1->SetFillStyle(0);
        leg1->SetTextSize(0.04);
        leg1->AddEntry(gV2_uncorr, "Uncorrected FourierFit", "lep");
        leg1->AddEntry(gV2_template, "TemplateFit (Bootstrap)", "lep");
        leg1->Draw();

        // V3 panel
        c->cd(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.10);
        gPad->SetBottomMargin(0.12);
        gPad->SetGridy();
        
        gV3_uncorr->SetTitle("");
        gV3_uncorr->GetYaxis()->SetTitle("V_{3}^{#Delta}");
        gV3_uncorr->GetXaxis()->SetTitle("TPC #eta");
        gV3_uncorr->GetYaxis()->SetTitleSize(0.05);
        gV3_uncorr->GetXaxis()->SetTitleSize(0.05);
        gV3_uncorr->GetYaxis()->SetRangeUser(0.0002, 0.00045);
        gV3_uncorr->GetXaxis()->SetLimits(-1.0, 1.0);
        gV3_uncorr->Draw("AP");
        gV3_template->Draw("P SAME");
        
        TLatex* lat2 = new TLatex();
        lat2->SetNDC();
        lat2->SetTextSize(0.045);
        lat2->DrawLatex(0.18, 0.92, Form("%s", datasetLabels[iDataset].c_str()));
        
        TLegend* leg2 = new TLegend(0.18, 0.70, 0.65, 0.88);
        leg2->SetBorderSize(0); leg2->SetFillStyle(0);
        leg2->SetTextSize(0.04);
        leg2->AddEntry(gV3_uncorr, "Uncorrected FourierFit", "lep");
        leg2->AddEntry(gV3_template, "TemplateFit (Bootstrap)", "lep");
        leg2->Draw();

        // V4 panel
        c->cd(3);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.05);
        gPad->SetTopMargin(0.10);
        gPad->SetBottomMargin(0.12);
        gPad->SetGridy();
        
        gV4_uncorr->SetTitle("");
        gV4_uncorr->GetYaxis()->SetTitle("V_{4}^{#Delta}");
        gV4_uncorr->GetXaxis()->SetTitle("TPC #eta");
        gV4_uncorr->GetYaxis()->SetTitleSize(0.05);
        gV4_uncorr->GetXaxis()->SetTitleSize(0.05);
        gV4_uncorr->GetYaxis()->SetRangeUser(0.00, 0.00012);
        gV4_uncorr->GetXaxis()->SetLimits(-1.0, 1.0);
        gV4_uncorr->Draw("AP");
        gV4_template->Draw("P SAME");
        
        TLatex* lat3 = new TLatex();
        lat3->SetNDC();
        lat3->SetTextSize(0.045);
        lat3->DrawLatex(0.18, 0.92, Form("%s", datasetLabels[iDataset].c_str()));
        
        TLegend* leg3 = new TLegend(0.18, 0.70, 0.65, 0.88);
        leg3->SetBorderSize(0); leg3->SetFillStyle(0);
        leg3->SetTextSize(0.04);
        leg3->AddEntry(gV4_uncorr, "Uncorrected FourierFit", "lep");
        leg3->AddEntry(gV4_template, "TemplateFit (Bootstrap)", "lep");
        leg3->Draw();

        std::string pdfName = Form("./TemplateFit/EtaDiff/MethodComparison/VnComparison_%s_Cent_0_20.pdf", dataset.c_str());
        c->SaveAs(pdfName.c_str());
        std::cout << "  Created: " << pdfName << std::endl;

        delete lat1; delete lat2; delete lat3;
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
