#include <iostream>
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLine.h"
#include "TLegend.h"

//==============================================================
// Calculate Flow Vector Decorrelations in eta
// r_n|n(η) = V_nΔ(FT0, TPC at η) / V_nΔ(FT0, TPC at -η)
// This measures how flow decorrelates with pseudorapidity separation
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

void Process_FlowDecorrelation() {
    
    std::cout << "\n=== Flow Vector Decorrelations in Pseudorapidity ===" << std::endl;
    std::cout << "r_n|n(η) = V_nΔ(FT0, TPC at η) / V_nΔ(FT0, TPC at -η)" << std::endl;
    std::cout << "where η refers to TPC pseudorapidity bins\n" << std::endl;
    
    // Create output directory
    gSystem->Exec("mkdir -p ./TemplateFit/EtaDiff/FlowDecorrelation");
    
    std::vector<std::string> datasets = {"615818", "615817", "616549", "618685"};
    std::vector<std::string> datasetLabels = {"Dataset 615818 Ne-Ne (Outer)", "Dataset 615817 Ne-Ne (Inner)",
                                               "Dataset 616549 O-O (Outer)", "Dataset 618685 O-O (Inner)"};
    
    // Define symmetric TPC eta bins
    // Positive eta bins (TPC-B side)
    std::vector<std::pair<Double_t, Double_t>> eta_bins_pos = {
        {0.7, 0.8}
    };
    
    // Corresponding negative eta bins (TPC-A side)
    std::vector<std::pair<Double_t, Double_t>> eta_bins_neg = {
        {-0.8, -0.7}
    };
    
    std::vector<Double_t> eta_centers = {0.75}; // Center of symmetric eta bin
    
    for (size_t ds_idx = 0; ds_idx < datasets.size(); ds_idx++) {
        const auto& dataset = datasets[ds_idx];
        const bool is_oo = (dataset == "616549" || dataset == "618685");
        const std::string file_prefix = is_oo ? "LHC25ae_pass2" : "LHC25af_pass2";
        
        std::cout << "\n=== " << datasetLabels[ds_idx] << " ===" << std::endl;
        
        std::vector<Double_t> r2_vals, r2_errs;
        std::vector<Double_t> r3_vals, r3_errs;
        std::vector<Double_t> r4_vals, r4_errs;
        
        for (size_t i = 0; i < eta_bins_pos.size(); i++) {
            Double_t eta_low_pos = eta_bins_pos[i].first;
            Double_t eta_high_pos = eta_bins_pos[i].second;
            Double_t eta_low_neg = eta_bins_neg[i].first;
            Double_t eta_high_neg = eta_bins_neg[i].second;
            
            // Read V_nΔ(FT0, TPC at +η)
            TString file_pos = Form("./TemplateFit/EtaDiff/Vn_%s_%s_TPC_FT0C_TPCEta_%.1f_%.1f_Cent_0_20.root",
                                   file_prefix.c_str(), dataset.c_str(), eta_low_pos, eta_high_pos);
            
            // Read V_nΔ(FT0, TPC at -η)
            TString file_neg = Form("./TemplateFit/EtaDiff/Vn_%s_%s_TPC_FT0C_TPCEta_%.1f_%.1f_Cent_0_20.root",
                                   file_prefix.c_str(), dataset.c_str(), eta_low_neg, eta_high_neg);
            
            VnResult vn_pos = ReadVnFromFile(file_pos);
            VnResult vn_neg = ReadVnFromFile(file_neg);
            
            if (vn_pos.v2 == 0 || vn_neg.v2 == 0) {
                std::cerr << "Warning: Could not read files for eta bins [" 
                         << eta_low_pos << ", " << eta_high_pos << "] and ["
                         << eta_low_neg << ", " << eta_high_neg << "]" << std::endl;
                continue;
            }
            
            // Calculate r_n|n = V_nΔ(+η) / V_nΔ(-η)
            Double_t r2 = vn_neg.v2 / vn_pos.v2;  // Note: neg/pos because we want r(η^a, η^b)/r(-η^a, η^b)
            Double_t r3 = vn_neg.v3 / vn_pos.v3;
            Double_t r4 = vn_neg.v4 / vn_pos.v4;
            
            // Error propagation for ratio: (da/a)^2 + (db/b)^2
            Double_t r2_err = r2 * TMath::Sqrt(TMath::Power(vn_neg.v2_err/vn_neg.v2, 2) + 
                                                TMath::Power(vn_pos.v2_err/vn_pos.v2, 2));
            Double_t r3_err = r3 * TMath::Sqrt(TMath::Power(vn_neg.v3_err/vn_neg.v3, 2) + 
                                                TMath::Power(vn_pos.v3_err/vn_pos.v3, 2));
            Double_t r4_err = r4 * TMath::Sqrt(TMath::Power(vn_neg.v4_err/vn_neg.v4, 2) + 
                                                TMath::Power(vn_pos.v4_err/vn_pos.v4, 2));
            
            std::cout << "\nη = ±" << eta_centers[i] << ":" << std::endl;
            std::cout << "  V₂Δ(FT0, -η): " << vn_neg.v2 << " ± " << vn_neg.v2_err << std::endl;
            std::cout << "  V₂Δ(FT0, +η): " << vn_pos.v2 << " ± " << vn_pos.v2_err << std::endl;
            std::cout << "  r₂|₂: " << r2 << " ± " << r2_err << std::endl;
            
            std::cout << "  V₃Δ(FT0, -η): " << vn_neg.v3 << " ± " << vn_neg.v3_err << std::endl;
            std::cout << "  V₃Δ(FT0, +η): " << vn_pos.v3 << " ± " << vn_pos.v3_err << std::endl;
            std::cout << "  r₃|₃: " << r3 << " ± " << r3_err << std::endl;
            
            std::cout << "  V₄Δ(FT0, -η): " << vn_neg.v4 << " ± " << vn_neg.v4_err << std::endl;
            std::cout << "  V₄Δ(FT0, +η): " << vn_pos.v4 << " ± " << vn_pos.v4_err << std::endl;
            std::cout << "  r₄|₄: " << r4 << " ± " << r4_err << std::endl;
            
            r2_vals.push_back(r2);
            r2_errs.push_back(r2_err);
            r3_vals.push_back(r3);
            r3_errs.push_back(r3_err);
            r4_vals.push_back(r4);
            r4_errs.push_back(r4_err);
        }
        
        // Save results to ROOT file
        if (r2_vals.size() != eta_centers.size()) {
            std::cerr << "Warning: Skipping output for " << dataset << " due to missing inputs" << std::endl;
            continue;
        }
        
        TString output_file = Form("./TemplateFit/EtaDiff/FlowDecorrelation/FlowDecorrelation_%s_%s_Cent_0_20.root",
                                  file_prefix.c_str(), dataset.c_str());
        
        TFile* out = new TFile(output_file, "RECREATE");
        
        TGraphErrors* g_r2 = new TGraphErrors(eta_centers.size());
        TGraphErrors* g_r3 = new TGraphErrors(eta_centers.size());
        TGraphErrors* g_r4 = new TGraphErrors(eta_centers.size());
        
        for (size_t i = 0; i < eta_centers.size(); i++) {
            g_r2->SetPoint(i, eta_centers[i], r2_vals[i]);
            g_r2->SetPointError(i, 0.0, r2_errs[i]);
            
            g_r3->SetPoint(i, eta_centers[i], r3_vals[i]);
            g_r3->SetPointError(i, 0.0, r3_errs[i]);
            
            g_r4->SetPoint(i, eta_centers[i], r4_vals[i]);
            g_r4->SetPointError(i, 0.0, r4_errs[i]);
        }
        
        g_r2->SetName("gR2");
        g_r2->SetTitle("r_{2|2} vs |#eta|");
        g_r3->SetName("gR3");
        g_r3->SetTitle("r_{3|3} vs |#eta|");
        g_r4->SetName("gR4");
        g_r4->SetTitle("r_{4|4} vs |#eta|");
        
        g_r2->Write();
        g_r3->Write();
        g_r4->Write();
        
        out->Close();
        delete out;
        
        std::cout << "\nSaved: " << output_file << std::endl;
    }
    
    // Create comparison plot for both datasets
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
    
    TCanvas canvas("c_decorr", "Flow Decorrelations", 1400, 900);
    canvas.SetLeftMargin(0.12);
    canvas.SetRightMargin(0.05);
    canvas.SetBottomMargin(0.12);
    canvas.SetTopMargin(0.10);
    
    // Read back the graphs for all datasets
    std::vector<TFile*> files;
    std::vector<TGraphErrors*> graphs;
    
    for (size_t i = 0; i < datasets.size(); ++i) {
        const bool is_oo = (datasets[i] == "616549" || datasets[i] == "618685");
        const std::string file_prefix = is_oo ? "LHC25ae_pass2" : "LHC25af_pass2";
        TFile* f = TFile::Open(Form("./TemplateFit/EtaDiff/FlowDecorrelation/FlowDecorrelation_%s_%s_Cent_0_20.root", 
                         file_prefix.c_str(), datasets[i].c_str()));
        if (!f) {
            std::cerr << "Error opening output file for " << datasets[i] << std::endl;
            continue;
        }
        
        TGraphErrors* g = (TGraphErrors*)f->Get("gR2");
        if (!g) {
            std::cerr << "Error reading graph from file for " << datasets[i] << std::endl;
            f->Close();
            continue;
        }
        
        files.push_back(f);
        graphs.push_back(g);
    }
    
    if (graphs.empty()) {
        std::cerr << "No graphs loaded for comparison" << std::endl;
        return;
    }
    
    // Define colors and markers for each dataset
    Int_t colors[] = {kRed+1, kBlue+1, kGreen+2, kMagenta+2};
    Int_t markers[] = {20, 21, 22, 23};
    Double_t offsets[] = {-0.03, -0.01, 0.01, 0.03};
    
    // Apply styling and offsets
    for (size_t i = 0; i < graphs.size(); ++i) {
        TGraphErrors* g = graphs[i];
        
        // Offset x-positions for visibility
        for (Int_t j = 0; j < g->GetN(); j++) {
            Double_t x, y;
            g->GetPoint(j, x, y);
            g->SetPoint(j, x + offsets[i], y);
        }
        
        g->SetMarkerStyle(markers[i]);
        g->SetMarkerSize(3.5);
        g->SetMarkerColor(colors[i]);
        g->SetLineColor(colors[i]);
        g->SetLineWidth(3);
        
        if (i == 0) {
            g->SetTitle("");
            g->GetXaxis()->SetTitle("|#eta_{TPC}|");
            g->GetXaxis()->SetLabelSize(0.05);
            g->GetXaxis()->SetTitleSize(0.05);
            g->GetXaxis()->SetTitleOffset(1.1);
            
            g->GetYaxis()->SetTitle("r_{2|2} = V_{2#Delta}(FT0,-#eta) / V_{2#Delta}(FT0,+#eta)");
            g->GetYaxis()->SetLabelSize(0.05);
            g->GetYaxis()->SetTitleSize(0.045);
            g->GetYaxis()->SetTitleOffset(1.3);
            g->GetYaxis()->SetRangeUser(0.85, 1.15);
            
            g->Draw("APE");
        } else {
            g->Draw("PE same");
        }
    }
    
    // Unity line
    TLine* unity = new TLine(graphs[0]->GetXaxis()->GetXmin(), 1.0, 
                             graphs[0]->GetXaxis()->GetXmax(), 1.0);
    unity->SetLineStyle(2);
    unity->SetLineWidth(2);
    unity->SetLineColor(kGray+2);
    unity->Draw();
    
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.14, 0.95, "ALICE Ne-Ne and O-O #sqrt{s_{NN}} = 17.3 GeV, 0-20% Centrality");
    
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.14, 0.83, "Flow Vector Decorrelation");
    
    TLegend* leg = new TLegend(0.55, 0.60, 0.93, 0.86);
    leg->SetBorderSize(1);
    leg->SetFillColor(kWhite);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    for (size_t i = 0; i < graphs.size(); ++i) {
        leg->AddEntry(graphs[i], datasetLabels[i].c_str(), "pe");
    }
    leg->Draw();
    
    canvas.SaveAs("./TemplateFit/EtaDiff/FlowDecorrelation/FlowDecorrelation_Comparison.pdf");
    
    std::cout << "\n=== Comparison Plot Created ===" << std::endl;
    std::cout << "Saved: ./TemplateFit/EtaDiff/FlowDecorrelation/FlowDecorrelation_Comparison.pdf" << std::endl;
    
    // Cleanup
    for (auto* f : files) {
        f->Close();
        delete f;
    }
    delete unity;
    delete leg;
    
    std::cout << "\n=== Flow Decorrelation Analysis Complete ===" << std::endl;
}
