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
#include "TLegend.h"

//==============================================================
// Compare ReconstructV2 using TPC-A vs TPC-B in numerator
// Method A: V2 = V2Δ(FT0, TPC-A) / sqrt(V2Δ(TPC-A, TPC-B))
// Method B: V2 = V2Δ(FT0, TPC-B) / sqrt(V2Δ(TPC-A, TPC-B))
// where TPC-A is η ∈ [-0.8, -0.7], TPC-B is η ∈ [0.7, 0.8]
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
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return result;
    }
    
    TH1D* hV2 = (TH1D*)file->Get("hV2");
    TH1D* hV3 = (TH1D*)file->Get("hV3");
    TH1D* hV4 = (TH1D*)file->Get("hV4");
    
    if (hV2) {
        result.v2 = hV2->GetBinContent(1);
        result.v2_err = hV2->GetBinError(1);
    }
    if (hV3) {
        result.v3 = hV3->GetBinContent(1);
        result.v3_err = hV3->GetBinError(1);
    }
    if (hV4) {
        result.v4 = hV4->GetBinContent(1);
        result.v4_err = hV4->GetBinError(1);
    }
    
    file->Close();
    delete file;
    return result;
}

void Process_ReconstructV2_TPC_B() {
    // Datasets to compare: Ne-Ne and O-O systems
    std::vector<std::string> ft0_datasets = {"615818", "615817", "616549", "618685"};
    std::vector<std::string> datasetLabels = {"615818 Ne-Ne Outer", "615817 Ne-Ne Inner", 
                                               "616549 O-O Outer", "618685 O-O Inner"};
    std::vector<std::string> detectors = {"FT0C"};
    
    std::cout << "\n=== Comparing ReconstructV2: TPC-A vs TPC-B in Numerator ===" << std::endl;
    std::cout << "Method A: V2 = V2Δ(FT0,TPC-A) / sqrt(V2Δ(TPC,TPC))" << std::endl;
    std::cout << "Method B: V2 = V2Δ(FT0,TPC-B) / sqrt(V2Δ(TPC,TPC))" << std::endl;
    std::cout << "where TPC-A = eta [-0.8, -0.7], TPC-B = eta [0.7, 0.8]" << std::endl;
    std::cout << "Denominator from: 611697 for Ne-Ne, 604830/604826 for O-O\n" << std::endl;
    
    // Create output directory
    gSystem->Exec("mkdir -p ./TemplateFit/EtaDiff/ReconstructedV2_Comparison");
    gSystem->Exec("mkdir -p ./TemplateFit/EtaDiff/ReconstructedV2_Comparison/Plots");
    
    std::vector<Double_t> v2_tpcA_vals, v2_tpcA_errs;
    std::vector<Double_t> v2_tpcB_vals, v2_tpcB_errs;
    
    // Loop over datasets
    for (size_t idx = 0; idx < ft0_datasets.size(); idx++) {
        const auto& dataset = ft0_datasets[idx];
        const auto& detector = detectors[0];
        
        // Determine which TPC-TPC dataset to use based on collision system
        const bool is_oo = (dataset == "616549" || dataset == "618685");
        const std::string file_prefix = is_oo ? "LHC25ae_pass2" : "LHC25af_pass2";
        std::string tpc_tpc_dataset = "611697";
        if (dataset == "616549") {
            tpc_tpc_dataset = "604830";
        } else if (dataset == "618685") {
            tpc_tpc_dataset = "604826";
        }
        const std::string tpc_tpc_prefix = is_oo ? "LHC25ae_pass2" : "LHC25af_pass2";
        const std::string collision_system = is_oo ? "O-O" : "Ne-Ne";
        
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
            delete tpc_file;
            continue;
        }
        
        printf("\n=== Dataset %s (%s) ===\n", dataset.c_str(), collision_system.c_str());
        printf("TPC-TPC %s: V2Delta(|eta|=%.2f, x=%.2f) = %.6f +/- %.6e\n",
               tpc_tpc_dataset.c_str(), target_eta, best_x, denom_v2delta, denom_v2delta_err);
        
        // Build file paths
        TString ft0_tpcA_file = Form("./TemplateFit/EtaDiff/Vn_%s_%s_TPC_%s_TPCEta_-0.8_-0.7_Cent_0_20.root", 
                         file_prefix.c_str(), dataset.c_str(), detector.c_str());
        TString ft0_tpcB_file = Form("./TemplateFit/EtaDiff/Vn_%s_%s_TPC_%s_TPCEta_0.7_0.8_Cent_0_20.root", 
                         file_prefix.c_str(), dataset.c_str(), detector.c_str());
        
        // Read FT0-TPC correlations
        VnResult ft0_tpcA = ReadVnFromFile(ft0_tpcA_file);
        VnResult ft0_tpcB = ReadVnFromFile(ft0_tpcB_file);
        
        if (ft0_tpcA.v2 == 0 || ft0_tpcB.v2 == 0) {
            std::cerr << "Warning: Could not read files for dataset " << dataset << std::endl;
            tpc_file->Close();
            delete tpc_file;
            continue;
        }
        
        printf("V2Delta(FT0C, TPC-A): %.6f +/- %.6e\n", ft0_tpcA.v2, ft0_tpcA.v2_err);
        printf("V2Delta(FT0C, TPC-B): %.6f +/- %.6e\n", ft0_tpcB.v2, ft0_tpcB.v2_err);
        printf("V2Delta(TPC-A, TPC-B): %.6f +/- %.6e\n", denom_v2delta, denom_v2delta_err);
        
        // Calculate using Method A (TPC-A in numerator)
        Double_t sqrt_denom = TMath::Sqrt(denom_v2delta);
        Double_t v2_methodA = ft0_tpcA.v2 / sqrt_denom;
        
        Double_t rel_err_numA = ft0_tpcA.v2_err / ft0_tpcA.v2;
        Double_t rel_err_sqrt = 0.5 * (denom_v2delta_err / denom_v2delta);
        Double_t rel_err_A = TMath::Sqrt(rel_err_numA * rel_err_numA + rel_err_sqrt * rel_err_sqrt);
        Double_t v2_methodA_err = v2_methodA * rel_err_A;
        
        // Calculate using Method B (TPC-B in numerator)
        Double_t v2_methodB = ft0_tpcB.v2 / sqrt_denom;
        
        Double_t rel_err_numB = ft0_tpcB.v2_err / ft0_tpcB.v2;
        Double_t rel_err_B = TMath::Sqrt(rel_err_numB * rel_err_numB + rel_err_sqrt * rel_err_sqrt);
        Double_t v2_methodB_err = v2_methodB * rel_err_B;
        
        printf("\nMethod A (using TPC-A): %.6f +/- %.6e\n", v2_methodA, v2_methodA_err);
        printf("Method B (using TPC-B): %.6f +/- %.6e\n", v2_methodB, v2_methodB_err);
        printf("Difference: %.6f (%.2f%%)\n\n", v2_methodA - v2_methodB, 
               100.0 * (v2_methodA - v2_methodB) / v2_methodB);
        
        // Store for plotting
        v2_tpcA_vals.push_back(v2_methodA);
        v2_tpcA_errs.push_back(v2_methodA_err);
        v2_tpcB_vals.push_back(v2_methodB);
        v2_tpcB_errs.push_back(v2_methodB_err);
        
        // Store Method B results
        TString output_file = Form("./TemplateFit/EtaDiff/ReconstructedV2_Comparison/ReconstructedV2_MethodB_%s_%s_TPC_%s_Cent_0_20.root", 
                       file_prefix.c_str(), dataset.c_str(), detector.c_str());
        
        TFile* out = new TFile(output_file, "RECREATE");
        
        TH1D* h_v2 = new TH1D("hV2Reconstructed_MethodB", "V_{2} Reconstructed (Method B)", 1, 0, 1);
        h_v2->SetBinContent(1, v2_methodB);
        h_v2->SetBinError(1, v2_methodB_err);
        
        TTree* tree = new TTree("V2", "Reconstructed V2 (Method B)");
        Double_t v2 = v2_methodB;
        Double_t v2_err = v2_methodB_err;
        tree->Branch("v2", &v2, "v2/D");
        tree->Branch("v2_err", &v2_err, "v2_err/D");
        tree->Fill();
        
        h_v2->Write();
        tree->Write();
        out->Close();
        delete out;
        
        tpc_file->Close();
        delete tpc_file;
    }
    
    // Create comparison plot
    if (!v2_tpcA_vals.empty() && v2_tpcA_vals.size() >= 2) {
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);
        
        TCanvas canvas("c_compare_tpc", "ReconstructV2: TPC-A vs TPC-B", 1400, 900);
        canvas.SetLeftMargin(0.12);
        canvas.SetRightMargin(0.05);
        canvas.SetBottomMargin(0.15);
        canvas.SetTopMargin(0.10);
        
        // Create graphs
        TGraphErrors* g_tpcA = new TGraphErrors(v2_tpcA_vals.size());
        TGraphErrors* g_tpcB = new TGraphErrors(v2_tpcB_vals.size());
        
        for (size_t i = 0; i < v2_tpcA_vals.size(); i++) {
            Double_t x_A = (Double_t)i - 0.1;
            Double_t x_B = (Double_t)i + 0.1;
            
            g_tpcA->SetPoint(i, x_A, v2_tpcA_vals[i]);
            g_tpcA->SetPointError(i, 0.0, v2_tpcA_errs[i]);
            
            g_tpcB->SetPoint(i, x_B, v2_tpcB_vals[i]);
            g_tpcB->SetPointError(i, 0.0, v2_tpcB_errs[i]);
        }
        
        // Style
        g_tpcA->SetMarkerStyle(20);
        g_tpcA->SetMarkerSize(3.0);
        g_tpcA->SetMarkerColor(kRed+1);
        g_tpcA->SetLineColor(kRed+1);
        g_tpcA->SetLineWidth(3);
        
        g_tpcB->SetMarkerStyle(21);
        g_tpcB->SetMarkerSize(3.0);
        g_tpcB->SetMarkerColor(kBlue+1);
        g_tpcB->SetLineColor(kBlue+1);
        g_tpcB->SetLineWidth(3);
        
        // Axes
        g_tpcA->SetTitle("");
        g_tpcA->GetXaxis()->SetLimits(-0.5, 3.5);
        g_tpcA->GetXaxis()->SetNdivisions(505);
        g_tpcA->GetXaxis()->SetLabelSize(0.06);
        g_tpcA->GetXaxis()->SetTitleSize(0.06);
        g_tpcA->GetXaxis()->SetTitle("Dataset");
        
        std::vector<Double_t> all_vals;
        for (size_t i = 0; i < v2_tpcA_vals.size(); ++i) {
            all_vals.push_back(v2_tpcA_vals[i]);
            all_vals.push_back(v2_tpcB_vals[i]);
        }
        Double_t min_val = *std::min_element(all_vals.begin(), all_vals.end());
        Double_t max_val = *std::max_element(all_vals.begin(), all_vals.end());
        Double_t center = (min_val + max_val) / 2.0;
        Double_t range = (max_val - min_val);
        Double_t y_min = center - range * 2.5;
        Double_t y_max = center + range * 2.5;
        
        g_tpcA->GetYaxis()->SetRangeUser(y_min, y_max);
        g_tpcA->GetYaxis()->SetLabelSize(0.06);
        g_tpcA->GetYaxis()->SetTitleSize(0.06);
        g_tpcA->GetYaxis()->SetTitleOffset(1.0);
        g_tpcA->GetYaxis()->SetTitle("V_{2}^{recon}");
        g_tpcA->GetYaxis()->SetNdivisions(508);
        
        g_tpcA->Draw("APE");
        g_tpcB->Draw("PE same");
        
        // Labels
        TLatex xlab;
        xlab.SetTextSize(0.045);
        xlab.SetTextFont(42);
        xlab.SetTextAlign(22);
        Double_t label_y = y_min - (y_max - y_min) * 0.08;
        xlab.DrawLatex(0.0, label_y, "615818 Ne-Ne Outer");
        xlab.DrawLatex(1.0, label_y, "615817 Ne-Ne Inner");
        xlab.DrawLatex(2.0, label_y, "616549 O-O Outer");
        xlab.DrawLatex(3.0, label_y, "618685 O-O Inner");
        
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        
        latex.SetTextSize(0.05);
        latex.DrawLatex(0.14, 0.95, "ALICE Ne-Ne and O-O #sqrt{s_{NN}} = 17.3 GeV, 0-20% Centrality");
        
        latex.SetTextSize(0.045);
        latex.DrawLatex(0.14, 0.89, "ReconstructV2: TPC-A vs TPC-B in Numerator");
        
        // Legend
        TLegend* leg = new TLegend(0.14, 0.75, 0.50, 0.86);
        leg->SetBorderSize(1);
        leg->SetFillColor(kWhite);
        leg->SetTextFont(42);
        leg->SetTextSize(0.04);
        leg->AddEntry(g_tpcA, "Method A: V_{2#Delta}(FT0,TPC-A)", "pe");
        leg->AddEntry(g_tpcB, "Method B: V_{2#Delta}(FT0,TPC-B)", "pe");
        leg->Draw();
        
        // Info box
        TPaveText *pt = new TPaveText(0.55, 0.55, 0.93, 0.86, "NDC");
        pt->SetFillColor(kWhite);
        pt->SetBorderSize(1);
        pt->SetTextAlign(12);
        pt->SetTextFont(42);
        pt->SetTextSize(0.03);
        pt->AddText("Method A (TPC-A):");
        for (size_t i = 0; i < v2_tpcA_vals.size(); ++i) {
            pt->AddText(Form("  %s: %.5f #pm %.5f", datasetLabels[i].c_str(), v2_tpcA_vals[i], v2_tpcA_errs[i]));
        }
        pt->AddText("");
        pt->AddText("Method B (TPC-B):");
        for (size_t i = 0; i < v2_tpcB_vals.size(); ++i) {
            pt->AddText(Form("  %s: %.5f #pm %.5f", datasetLabels[i].c_str(), v2_tpcB_vals[i], v2_tpcB_errs[i]));
        }
        pt->Draw();
        
        canvas.SaveAs("./TemplateFit/EtaDiff/ReconstructedV2_Comparison/Plots/ReconstructV2_TPC_A_vs_B.pdf");
        
        std::cout << "\n=== Comparison Plot Created ===" << std::endl;
        std::cout << "Saved: ./TemplateFit/EtaDiff/ReconstructedV2_Comparison/Plots/ReconstructV2_TPC_A_vs_B.pdf" << std::endl;
        
        delete pt;
        delete leg;
        delete g_tpcA;
        delete g_tpcB;
    }
    
    std::cout << "\n=== Comparison Complete ===" << std::endl;
}
