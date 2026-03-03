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
// Calculate V2 using 3x2PC method
// Formula: V2_3x2PC(FT0) = sqrt[(V2Δ(FT0,TPC-A) × V2Δ(FT0,TPC-B)) / V2Δ(TPC-A,TPC-B)]
// where TPC-A is η_TPC = -0.8 to -0.7
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

void Process_3x2PC() {
    // Datasets to compare: Ne-Ne and O-O systems
    std::vector<std::string> ft0_datasets = {"615818", "615817", "616549", "618685"};
    std::vector<std::string> datasetLabels = {"615818 Ne-Ne Outer", "615817 Ne-Ne Inner", 
                                               "616549 O-O Outer", "618685 O-O Inner"};
    std::vector<std::string> detectors = {"FT0C"};
    
    std::cout << "\n=== Calculating V2 using 3x2PC Method ===" << std::endl;
    std::cout << "Formula: V2_3x2PC = sqrt[(V2Δ(FT0,TPC-A) × V2Δ(FT0,TPC-B)) / V2Δ(TPC-A,TPC-B)]" << std::endl;
    std::cout << "where TPC-A = eta [-0.8, -0.7], TPC-B = eta [0.7, 0.8]" << std::endl;
    std::cout << "FT0-TPC correlations from: respective datasets (LR)" << std::endl;
    std::cout << "TPC-TPC correlation from: 611697 for Ne-Ne, 604830/604826 for O-O (SR)\n" << std::endl;
    
    // Create output directory
    gSystem->Exec("mkdir -p ./TemplateFit/EtaDiff/V2_3x2PC");
    gSystem->Exec("mkdir -p ./TemplateFit/EtaDiff/V2_3x2PC/Plots");
    
    std::vector<Double_t> v2_values;
    std::vector<Double_t> v2_errors;
    
    // Loop over datasets and calculate V2 using 3x2PC
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
        
        // Extract V2Delta at |eta|=0.75 (representing TPC-A to TPC-B correlation)
        const Double_t target_eta = 0.75;
        Double_t v2delta_tpc_tpc = 0, v2delta_tpc_tpc_err = 0;
        Double_t best_diff = 1e9;
        Double_t best_x = 0.0;
        for (Int_t i = 0; i < g_tpc_tpc->GetN(); i++) {
            Double_t x, y;
            g_tpc_tpc->GetPoint(i, x, y);
            Double_t diff = TMath::Abs(TMath::Abs(x) - target_eta);
            if (diff < best_diff) {
                best_diff = diff;
                best_x = x;
                v2delta_tpc_tpc = y;
                v2delta_tpc_tpc_err = g_tpc_tpc->GetErrorY(i);
            }
        }
        
        if (v2delta_tpc_tpc == 0 || best_diff > 0.2) {
            std::cerr << "Error: Could not find |eta|=0.75 point in TPC-TPC correlation" << std::endl;
            tpc_file->Close();
            delete tpc_file;
            continue;
        }
        
        printf("\n=== Dataset %s (%s) ===\n", dataset.c_str(), collision_system.c_str());
        printf("TPC-TPC correlation (%s): V2Delta(|eta|=%.2f, x=%.2f) = %.6f +/- %.6e\n", 
               tpc_tpc_dataset.c_str(), target_eta, best_x, v2delta_tpc_tpc, v2delta_tpc_tpc_err);
        
        // Build file paths for FT0-TPC-A and FT0-TPC-B correlations
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
        
        // Print input values
        printf("V2Delta(FT0C, TPC-A): %.6f +/- %.6e\n", ft0_tpcA.v2, ft0_tpcA.v2_err);
        printf("V2Delta(FT0C, TPC-B): %.6f +/- %.6e\n", ft0_tpcB.v2, ft0_tpcB.v2_err);
        printf("V2Delta(TPC-A, TPC-B): %.6f +/- %.6e\n", v2delta_tpc_tpc, v2delta_tpc_tpc_err);
        
        // Calculate V2 using 3x2PC formula:
        // V2 = sqrt[(V2Δ(FT0,TPC-A) × V2Δ(FT0,TPC-B)) / V2Δ(TPC-A,TPC-B)]
        Double_t numerator = ft0_tpcA.v2 * ft0_tpcB.v2;
        Double_t denominator = v2delta_tpc_tpc;
        Double_t v2_3x2pc = TMath::Sqrt(numerator / denominator);
        
        // Error propagation for: V = sqrt[(A × B) / C]
        // (dV/V)^2 = (dA/A)^2 + (dB/B)^2 + (dC/C)^2 / 4
        Double_t rel_err_A = ft0_tpcA.v2_err / ft0_tpcA.v2;
        Double_t rel_err_B = ft0_tpcB.v2_err / ft0_tpcB.v2;
        Double_t rel_err_C = v2delta_tpc_tpc_err / v2delta_tpc_tpc;
        
        Double_t rel_err_product = TMath::Sqrt(rel_err_A * rel_err_A + rel_err_B * rel_err_B);
        Double_t rel_err_ratio = TMath::Sqrt(rel_err_product * rel_err_product + rel_err_C * rel_err_C);
        Double_t rel_err_sqrt = 0.5 * rel_err_ratio;
        Double_t v2_3x2pc_err = v2_3x2pc * rel_err_sqrt;
        
        printf("\nError breakdown:\n");
        printf("  FT0-TPC-A relative error: %.4f%%\n", rel_err_A * 100);
        printf("  FT0-TPC-B relative error: %.4f%%\n", rel_err_B * 100);
        printf("  TPC-TPC relative error: %.4f%%\n", rel_err_C * 100);
        printf("  Product relative error: %.4f%%\n", rel_err_product * 100);
        printf("  Final relative error: %.4f%%\n", rel_err_sqrt * 100);
        
        printf("\nV2_3x2PC: %.6f +/- %.6e\n\n", v2_3x2pc, v2_3x2pc_err);
        
        // Store for plotting
        v2_values.push_back(v2_3x2pc);
        v2_errors.push_back(v2_3x2pc_err);
        
        // Store results in output file
        TString output_file = Form("./TemplateFit/EtaDiff/V2_3x2PC/V2_3x2PC_%s_%s_TPC_%s_Cent_0_20.root", 
                       file_prefix.c_str(), dataset.c_str(), detector.c_str());
        
        TFile* out = new TFile(output_file, "RECREATE");
        
        TH1D* h_v2 = new TH1D("hV2_3x2PC", "V_{2} from 3x2PC", 1, 0, 1);
        h_v2->SetBinContent(1, v2_3x2pc);
        h_v2->SetBinError(1, v2_3x2pc_err);
        
        TTree* tree = new TTree("V2_3x2PC", "V2 from 3x2PC");
        Double_t v2 = v2_3x2pc;
        Double_t v2_err = v2_3x2pc_err;
        tree->Branch("v2", &v2, "v2/D");
        tree->Branch("v2_err", &v2_err, "v2_err/D");
        tree->Fill();
        
        h_v2->Write();
        tree->Write();
        out->Close();
        delete out;
        
        tpc_file->Close();
        delete tpc_file;
        
        printf("Output saved: %s\n\n", output_file.Data());
    }
    
    // Create comparison plot
    if (!v2_values.empty() && v2_values.size() >= 2) {
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);
        
        TCanvas canvas("c_v2_3x2pc", "V2 from 3x2PC Comparison", 1200, 900);
        canvas.SetLeftMargin(0.15);
        canvas.SetRightMargin(0.05);
        canvas.SetBottomMargin(0.15);
        canvas.SetTopMargin(0.12);
        
        TGraphErrors* graph = new TGraphErrors(v2_values.size());
        
        for (size_t i = 0; i < v2_values.size(); i++) {
            graph->SetPoint(i, (Double_t)i, v2_values[i]);
            graph->SetPointError(i, 0.0, v2_errors[i]);
        }
        
        graph->SetTitle("");
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(3.5);
        graph->SetMarkerColor(kBlue+1);
        graph->SetLineColor(kBlue+1);
        graph->SetLineWidth(3);
        
        // Set up axes
        graph->GetXaxis()->SetLimits(-0.5, 3.5);
        graph->GetXaxis()->SetNdivisions(505);
        graph->GetXaxis()->SetLabelSize(0.06);
        graph->GetXaxis()->SetTitleSize(0.06);
        graph->GetXaxis()->SetLabelFont(42);
        graph->GetXaxis()->SetTitleFont(42);
        graph->GetXaxis()->SetTitle("Dataset");
        
        // Zoom y-axis
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
        graph->GetYaxis()->SetTitle("V_{2}^{3x2PC}");
        graph->GetYaxis()->SetNdivisions(508);
        
        graph->Draw("APE");
        
        // Draw connecting line
        graph->SetLineColor(kGray+2);
        graph->SetLineStyle(2);
        graph->SetLineWidth(2);
        graph->Draw("PE same");
        
        // Add x-axis labels
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
        
        latex.SetTextSize(0.038);
        latex.DrawLatex(0.18, 0.83, "V_{2}^{3x2PC} = #sqrt{[V_{2#Delta}(FT0,TPC-A) #times V_{2#Delta}(FT0,TPC-B)] / V_{2#Delta}(TPC,TPC)}");
        
        // Add value box
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
        
        canvas.SaveAs("./TemplateFit/EtaDiff/V2_3x2PC/Plots/V2_3x2PC_Comparison.pdf");
        
        printf("\n=== Comparison Plot ===\n");
        printf("Plot saved: ./TemplateFit/EtaDiff/V2_3x2PC/Plots/V2_3x2PC_Comparison.pdf\n");
        printf("Y-axis range: %.6f to %.6f\n", y_min, y_max);
        printf("Difference: %.5f\n", v2_values[0] - v2_values[1]);
        
        delete pt;
        delete graph;
    }
    
    std::cout << "\n=== 3x2PC Calculation Complete ===" << std::endl;
}
