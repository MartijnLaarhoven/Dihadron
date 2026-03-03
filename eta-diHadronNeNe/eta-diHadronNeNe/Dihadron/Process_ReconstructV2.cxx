#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLine.h"
#include "TPaveText.h"

//==============================================================
// ReconstructV2 from V2Delta values
// Formula: V_n^recon(η^FT0) = V_nΔ(η^FT0, η^TPC-A) / sqrt(V_nΔ(η^TPC-B, η^TPC-A))
// where TPC-A is η_TPC = -0.8 to -0.7 (closest to FT0C)
// and TPC-B is η_TPC = 0.7 to 0.8
//==============================================================

struct VnResult {
    Double_t v2, v2_err;
    Double_t v3, v3_err;
    Double_t v4, v4_err;
};

VnResult ReadVnFromFile(TString filename) {
    VnResult result = {0, 0, 0, 0, 0, 0};
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return result;
    }
    
    // Read histograms instead of tree
    TH1D* hV2 = (TH1D*)file->Get("hV2");
    TH1D* hV3 = (TH1D*)file->Get("hV3");
    TH1D* hV4 = (TH1D*)file->Get("hV4");
    
    if (hV2 && hV2->GetNbinsX() > 0) {
        result.v2 = hV2->GetBinContent(1);
        result.v2_err = hV2->GetBinError(1);
    }
    if (hV3 && hV3->GetNbinsX() > 0) {
        result.v3 = hV3->GetBinContent(1);
        result.v3_err = hV3->GetBinError(1);
    }
    if (hV4 && hV4->GetNbinsX() > 0) {
        result.v4 = hV4->GetBinContent(1);
        result.v4_err = hV4->GetBinError(1);
    }
    
    file->Close();
    delete file;
    return result;
}

void ReconstructV2() {
    // Compare Ne-Ne and O-O systems
    // Using TPC-TPC correlations from 611697 (for short-range) 
    // and FT0-TPC correlations from the respective datasets (for long-range with FT0)
    std::vector<std::string> ft0_datasets = {"615818", "615817", "616549", "618685"};
    std::vector<std::string> datasetLabels = {"615818\nNe-Ne Outer", "615817\nNe-Ne Inner", 
                                               "616549\nO-O Outer", "618685\nO-O Inner"};
    std::vector<std::string> detectors = {"FT0C"};
    
    std::cout << "\n=== Reconstructing V2 from V2Delta ===" << std::endl;
    std::cout << "Formula: V2_recon = V2Delta(FT0, TPC-A) / sqrt(V2Delta(TPC-B, TPC-A))" << std::endl;
    std::cout << "where TPC-A = eta [-0.8, -0.7], TPC-B = eta [0.7, 0.8]" << std::endl;
    std::cout << "Numerator from: respective datasets (LR with FT0)" << std::endl;
    std::cout << "Denominator from: 611697 for Ne-Ne, 604830/604826 for O-O (SR TPC-TPC)\n" << std::endl;
    
    // Create output directory for reconstructed V2
    gSystem->Exec("mkdir -p ./TemplateFit/EtaDiff/ReconstructedV2");
    gSystem->Exec("mkdir -p ./TemplateFit/EtaDiff/ReconstructedV2/Plots");
    
    std::vector<Double_t> v2_values;
    std::vector<Double_t> v2_errors;
    
    for (size_t idx = 0; idx < ft0_datasets.size(); idx++) {
        const auto& dataset = ft0_datasets[idx];
        const auto& detector = detectors[0];
        const bool is_oo = (dataset == "616549" || dataset == "618685");
        const std::string file_prefix = is_oo ? "LHC25ae_pass2" : "LHC25af_pass2";
        std::string tpc_tpc_dataset = "611697";
        if (dataset == "616549") {
            tpc_tpc_dataset = "604830";
        } else if (dataset == "618685") {
            tpc_tpc_dataset = "604826";
        }
        const std::string tpc_tpc_prefix = is_oo ? "LHC25ae_pass2" : "LHC25af_pass2";
        
        // Read TPC-TPC correlation (denominator)
        TString tpc_tpc_file = Form("/home/martijn-laarhoven/Work/dihadronanalysis-master/Dihadron/TemplateFit/EtaDiff/Vn_%s_%s_Cent_0_20.root",
                                    tpc_tpc_prefix.c_str(), tpc_tpc_dataset.c_str());
        TFile* tpc_file = TFile::Open(tpc_tpc_file);
        if (!tpc_file || tpc_file->IsZombie()) {
            std::cerr << "Error: Could not open " << tpc_tpc_file << std::endl;
            continue;
        }
        
        TGraphErrors* g_tpc_tpc = (TGraphErrors*)tpc_file->Get("gV2Delta");
        if (!g_tpc_tpc) {
            std::cerr << "Error: Could not find gV2Delta in " << tpc_tpc_file << std::endl;
            tpc_file->Close();
            continue;
        }
        
        // Extract V2Delta at |eta|=0.75 (TPC-A to TPC-B correlation)
        const Double_t target_eta = 0.75;
        Double_t denom_v2delta = 0, denom_v2delta_err = 0;
        Double_t best_diff = 1e9;
        Double_t best_x = 0.0;
        for (Int_t i = 0; i < g_tpc_tpc->GetN(); i++) {
            Double_t x, y;
            g_tpc_tpc->GetPoint(i, x, y);
            Double_t diff = TMath::Abs(TMath::Abs(x) - target_eta);
            if (diff < best_diff) {
                best_diff = diff;
                best_x = x;
                denom_v2delta = y;
                denom_v2delta_err = g_tpc_tpc->GetErrorY(i);
            }
        }
        
        if (denom_v2delta == 0 || best_diff > 0.2) {
            std::cerr << "Error: Could not find |eta|=0.75 point in TPC-TPC correlation" << std::endl;
            tpc_file->Close();
            continue;
        }
        
        printf("TPC-TPC %s: V2Delta(|eta|=%.2f, x=%.2f) = %.6f +/- %.6e\n\n",
               tpc_tpc_dataset.c_str(), target_eta, best_x, denom_v2delta, denom_v2delta_err);
        
        // Build file path for numerator (FT0-TPC-A)
        TString numerator_file = Form("./TemplateFit/EtaDiff/Vn_%s_%s_TPC_%s_TPCEta_-0.8_-0.7_Cent_0_20.root", 
                          file_prefix.c_str(), dataset.c_str(), detector.c_str());
        
        // Read numerator (FT0-TPC correlation at -0.8 to -0.7)
        VnResult numerator = ReadVnFromFile(numerator_file);
        
        if (numerator.v2 == 0) {
            std::cerr << "Warning: Could not read file for dataset " << dataset << " detector " << detector << std::endl;
            tpc_file->Close();
            continue;
        }
        
        // Print debugging info
        printf("Numerator V2_Delta(FT0, TPC-A): %.6f +/- %.6e\n", numerator.v2, numerator.v2_err);
        printf("Denominator V2_Delta(TPC-B, TPC-A): %.6f +/- %.6e (from %s)\n", denom_v2delta, denom_v2delta_err, tpc_tpc_dataset.c_str());
        
        // Calculate V2_recon using formula: V2_recon = V2Delta(FT0, TPC-A) / sqrt(V2Delta(TPC-B, TPC-A))
        Double_t sqrt_denom = TMath::Sqrt(denom_v2delta);
        Double_t v2_recon = numerator.v2 / sqrt_denom;
        
        // Error propagation for division by sqrt:
        // V = A / sqrt(B)
        // (dV/V)^2 = (dA/A)^2 + (dB/(2B))^2
        Double_t rel_err_num = numerator.v2_err / numerator.v2;
        Double_t rel_err_sqrt_term = 0.5 * (denom_v2delta_err / denom_v2delta);
        Double_t rel_err_recon = TMath::Sqrt(rel_err_num * rel_err_num + rel_err_sqrt_term * rel_err_sqrt_term);
        Double_t v2_recon_err = v2_recon * rel_err_recon;
        
        printf("Numerator relative error: %.4f%%\n", rel_err_num * 100);
        printf("Denominator relative error: %.4f%%\n", denom_v2delta_err / denom_v2delta * 100);
        printf("sqrt(Denominator) relative error: %.4f%%\n", rel_err_sqrt_term * 100);
        printf("V2_recon relative error: %.4f%%\n", rel_err_recon*100);
        printf("V2_recon: %.6f +/- %.6e\n\n", v2_recon, v2_recon_err);
        
        // Store for plotting
        v2_values.push_back(v2_recon);
        v2_errors.push_back(v2_recon_err);
        
        // Store results in output file
        TString output_file = Form("./TemplateFit/EtaDiff/ReconstructedV2/ReconstructedV2_%s_%s_TPC_%s_Cent_0_20.root", 
                       file_prefix.c_str(), dataset.c_str(), detector.c_str());
        
        TFile* out = new TFile(output_file, "RECREATE");
        
        // Create histogram to store result
        TH1D* h_v2 = new TH1D("hV2Reconstructed", "Reconstructed V2", 1, 0, 1);
        h_v2->SetBinContent(1, v2_recon);
        h_v2->SetBinError(1, v2_recon_err);
        
        // Also create tree for consistency with other outputs
        TTree* tree = new TTree("V2", "Reconstructed V2");
        Double_t v2 = v2_recon;
        Double_t v2_err = v2_recon_err;
        tree->Branch("v2", &v2, "v2/D");
        tree->Branch("v2_err", &v2_err, "v2_err/D");
        tree->Fill();
        
        h_v2->Write();
        tree->Write();
        out->Close();
        delete out;
        
        tpc_file->Close();
        
        printf("Dataset %s, %s: V2_recon = %.6f +/- %.6e\n", 
               dataset.c_str(), detector.c_str(), v2_recon, v2_recon_err);
    }
    
    // Create comparison plot
    if (!v2_values.empty() && v2_values.size() >= 2) {
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);
        
        TCanvas canvas("c_v2_comparison", "Reconstructed V2 Comparison", 1200, 900);
        canvas.SetLeftMargin(0.15);
        canvas.SetRightMargin(0.05);
        canvas.SetBottomMargin(0.15);
        canvas.SetTopMargin(0.12);
        
        // Create TGraphErrors with proper error bars
        TGraphErrors* graph = new TGraphErrors(v2_values.size());
        
        for (size_t i = 0; i < v2_values.size(); i++) {
            graph->SetPoint(i, (Double_t)i, v2_values[i]);
            graph->SetPointError(i, 0.0, v2_errors[i]);
        }
        
        graph->SetTitle("");
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(3.5);
        graph->SetMarkerColor(kRed+1);
        graph->SetLineColor(kRed+1);
        graph->SetLineWidth(3);
        
        // Set up axes with tight zoom to show error bars
        graph->GetXaxis()->SetLimits(-0.5, 3.5);
        graph->GetXaxis()->SetNdivisions(505);
        graph->GetXaxis()->SetLabelSize(0.06);
        graph->GetXaxis()->SetTitleSize(0.06);
        graph->GetXaxis()->SetLabelFont(42);
        graph->GetXaxis()->SetTitleFont(42);
        graph->GetXaxis()->SetTitle("Dataset");
        
        // Zoom in on y-axis to make error bars visible
        Double_t min_val = *std::min_element(v2_values.begin(), v2_values.end());
        Double_t max_val = *std::max_element(v2_values.begin(), v2_values.end());
        Double_t center = (min_val + max_val) / 2.0;
        Double_t range = (max_val - min_val);
        Double_t y_min = center - range * 3.0;
        Double_t y_max = center + range * 3.0;
        
        graph->GetYaxis()->SetRangeUser(y_min, y_max);
        graph->GetYaxis()->SetLabelSize(0.06);
        graph->GetYaxis()->SetTitleSize(0.06);
        graph->GetYaxis()->SetTitleOffset(1.3);
        graph->GetYaxis()->SetLabelFont(42);
        graph->GetYaxis()->SetTitleFont(42);
        graph->GetYaxis()->SetTitle("V_{2}^{recon}");
        graph->GetYaxis()->SetNdivisions(508);
        
        // Draw with error bars
        graph->Draw("APE");
        
        // Draw a line connecting the points
        graph->SetLineColor(kGray+2);
        graph->SetLineStyle(2);
        graph->SetLineWidth(2);
        graph->Draw("PE same");
        
        // Add custom x-axis labels
        TLatex xlab;
        xlab.SetTextSize(0.045);
        xlab.SetTextFont(42);
        xlab.SetTextAlign(22);
        Double_t label_y = y_min - (y_max - y_min) * 0.08;
        xlab.DrawLatex(0.0, label_y, "615818 Ne-Ne Outer");
        xlab.DrawLatex(1.0, label_y, "615817 Ne-Ne Inner");
        xlab.DrawLatex(2.0, label_y, "616549 O-O Outer");
        xlab.DrawLatex(3.0, label_y, "618685 O-O Inner");
        
        // Add title and information
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        
        latex.SetTextSize(0.05);
        latex.DrawLatex(0.18, 0.95, "ALICE Ne-Ne and O-O");
        
        latex.SetTextSize(0.045);
        latex.DrawLatex(0.18, 0.89, "Centrality: 0-20%");
        
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.18, 0.83, "V_{2}^{recon} = V_{2#Delta}(FT0,TPC) / #sqrt{V_{2#Delta}(TPC,TPC)}");
        
        // Add a box with values
        TPaveText *pt = new TPaveText(0.55, 0.60, 0.92, 0.88, "NDC");
        pt->SetFillColor(kWhite);
        pt->SetBorderSize(1);
        pt->SetTextAlign(12);
        pt->SetTextFont(42);
        pt->SetTextSize(0.035);
        for (size_t i = 0; i < v2_values.size(); i++) {
            pt->AddText(Form("%s: %.5f #pm %.5f", ft0_datasets[i].c_str(), v2_values[i], v2_errors[i]));
        }
        pt->Draw();
        
        canvas.SaveAs("./TemplateFit/EtaDiff/ReconstructedV2/Plots/ReconstructedV2_Comparison.pdf");
        
        printf("\nPlot saved: ./TemplateFit/EtaDiff/ReconstructedV2/Plots/ReconstructedV2_Comparison.pdf\n");
        printf("Y-axis range: %.6f to %.6f\n", y_min, y_max);
        printf("Difference: %.5f\n", v2_values[0] - v2_values[1]);
        
        delete pt;
        delete graph;
    }
    
    // Create individual plots for each dataset
    for (size_t idx = 0; idx < ft0_datasets.size() && idx < v2_values.size(); idx++) {
        TCanvas c("c_individual", Form("V2_recon_%s", ft0_datasets[idx].c_str()), 800, 600);
        c.SetLeftMargin(0.15);
        c.SetRightMargin(0.1);
        c.SetBottomMargin(0.15);
        c.SetTopMargin(0.1);
        
        TH1D* h_single = new TH1D("h_single", "", 1, 0, 1);
        h_single->SetBinContent(1, v2_values[idx]);
        h_single->SetBinError(1, v2_errors[idx]);
        
        h_single->SetTitle(Form("Reconstructed V_{2}^{recon} - Dataset %s (0-20%% Centrality)", ft0_datasets[idx].c_str()));
        h_single->SetYTitle("V_{2}^{recon}");
        h_single->SetXTitle("");
        h_single->SetMarkerStyle(20);
        h_single->SetMarkerSize(2.0);
        h_single->SetMarkerColor(kBlue);
        h_single->SetLineColor(kBlue);
        h_single->SetLineWidth(2);
        h_single->GetXaxis()->SetLabelSize(0);
        h_single->GetYaxis()->SetLabelSize(0.04);
        h_single->GetYaxis()->SetTitleOffset(1.5);
        h_single->GetYaxis()->SetRangeUser(0, v2_values[idx] * 1.5);
        
        h_single->Draw("E1");
        
        // Add value as text
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.06);
        latex.SetTextFont(62);
        latex.DrawLatex(0.3, 0.7, Form("V_{2}^{recon} = %.6f", v2_values[idx]));
        
        latex.SetTextSize(0.05);
        latex.SetTextFont(42);
        latex.DrawLatex(0.3, 0.6, Form("Error = %.6e", v2_errors[idx]));
        
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.15, 0.95, "ALICE Ne-Ne #sqrt{s_{NN}} = 17.3 GeV, 0-20% Centrality");
        latex.DrawLatex(0.15, 0.90, "V_{2}^{recon} = V_{2#Delta}(FT0,TPC-A) / #sqrt{V_{2#Delta}(TPC-B,TPC-A)}");
        
        c.SaveAs(Form("./TemplateFit/EtaDiff/ReconstructedV2/Plots/ReconstructedV2_%s.pdf", ft0_datasets[idx].c_str()));
        
        printf("Plot saved: ./TemplateFit/EtaDiff/ReconstructedV2/Plots/ReconstructedV2_%s.pdf\n", ft0_datasets[idx].c_str());
        
        delete h_single;
    }
    
    std::cout << "\n=== Reconstruction Complete ===" << std::endl;
}

// Execute when macro is run
void Process_ReconstructV2() {
    ReconstructV2();
}
