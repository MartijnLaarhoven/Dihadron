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
// Compare all 3 V2 calculation methods:
// Method 1 (Recon-A): V2 = V2Δ(FT0,TPC-A) / sqrt(V2Δ(TPC,TPC))
// Method 2 (Recon-B): V2 = V2Δ(FT0,TPC-B) / sqrt(V2Δ(TPC,TPC))
// Method 3 (3x2PC): V2 = sqrt[(V2Δ(FT0,TPC-A) × V2Δ(FT0,TPC-B)) / V2Δ(TPC,TPC)]
//==============================================================

struct V2Result {
    Double_t v2;
    Double_t v2_err;
};

V2Result ReadV2FromFile(TString filename, const char* histname = "hV2") {
    V2Result result = {0.0, 0.0};
    
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return result;
    }
    
    // Try multiple histogram names
    TH1D* h = (TH1D*)file->Get(histname);
    if (!h) h = (TH1D*)file->Get("hV2Reconstructed");
    if (!h) h = (TH1D*)file->Get("hV2Reconstructed_MethodB");
    if (!h) h = (TH1D*)file->Get("hV2_3x2PC");
    
    if (!h) {
        std::cerr << "Error: Could not find V2 histogram in " << filename << std::endl;
        file->Close();
        return result;
    }
    
    result.v2 = h->GetBinContent(1);
    result.v2_err = h->GetBinError(1);
    
    file->Close();
    delete file;
    return result;
}

void Process_CompareAll3Methods() {
    
    std::vector<std::string> datasets = {"615818", "615817", "616549", "618685"};
    std::vector<std::string> datasetLabels = {"615818 Ne-Ne (Outer)", "615817 Ne-Ne (Inner)", 
                                               "616549 O-O (Outer)", "618685 O-O (Inner)"};
    
    std::cout << "\n=== Comparing All 3 V2 Methods ===" << std::endl;
    std::cout << "Method 1 (Recon-A): V2 = V2Δ(FT0,TPC-A) / sqrt(V2Δ(TPC,TPC))" << std::endl;
    std::cout << "Method 2 (Recon-B): V2 = V2Δ(FT0,TPC-B) / sqrt(V2Δ(TPC,TPC))" << std::endl;
    std::cout << "Method 3 (3x2PC):   V2 = sqrt[(V2Δ(FT0,TPC-A) × V2Δ(FT0,TPC-B)) / V2Δ(TPC,TPC)]" << std::endl;
    std::cout << "\nwhere TPC-A = eta [-0.8, -0.7], TPC-B = eta [0.7, 0.8]\n" << std::endl;
    
    // Create output directory
    gSystem->Exec("mkdir -p ./TemplateFit/EtaDiff/All3Methods");
    
    std::vector<Double_t> v2_reconA_vals, v2_reconA_errs;
    std::vector<Double_t> v2_reconB_vals, v2_reconB_errs;
    std::vector<Double_t> v2_3x2pc_vals, v2_3x2pc_errs;
    
    // Read results from all methods
    for (size_t idx = 0; idx < datasets.size(); idx++) {
        const auto& dataset = datasets[idx];
        
        const bool is_oo = (dataset == "616549" || dataset == "618685");
        const std::string file_prefix = is_oo ? "LHC25ae_pass2" : "LHC25af_pass2";
        
        std::cout << "=== Dataset " << dataset << " ===" << std::endl;
        
        // Read Method 1 (Recon-A)
        TString reconA_file = Form("./TemplateFit/EtaDiff/ReconstructedV2/ReconstructedV2_%s_%s_TPC_FT0C_Cent_0_20.root", 
                       file_prefix.c_str(), dataset.c_str());
        V2Result reconA = ReadV2FromFile(reconA_file, "hV2Reconstructed");
        
        // Read Method 2 (Recon-B)
        TString reconB_file = Form("./TemplateFit/EtaDiff/ReconstructedV2_Comparison/ReconstructedV2_MethodB_%s_%s_TPC_FT0C_Cent_0_20.root", 
                       file_prefix.c_str(), dataset.c_str());
        V2Result reconB = ReadV2FromFile(reconB_file, "hV2Reconstructed_MethodB");
        
        // Read Method 3 (3x2PC)
        TString pc3x2_file = Form("./TemplateFit/EtaDiff/V2_3x2PC/V2_3x2PC_%s_%s_TPC_FT0C_Cent_0_20.root", 
                      file_prefix.c_str(), dataset.c_str());
        V2Result pc3x2 = ReadV2FromFile(pc3x2_file, "hV2_3x2PC");
        
        printf("  Recon-A (TPC-A): %.6f +/- %.6e\n", reconA.v2, reconA.v2_err);
        printf("  Recon-B (TPC-B): %.6f +/- %.6e\n", reconB.v2, reconB.v2_err);
        printf("  3x2PC:           %.6f +/- %.6e\n", pc3x2.v2, pc3x2.v2_err);
        
        // Calculate differences
        Double_t diff_A_3x2pc = reconA.v2 - pc3x2.v2;
        Double_t diff_B_3x2pc = reconB.v2 - pc3x2.v2;
        Double_t diff_A_B = reconA.v2 - reconB.v2;
        
        printf("\n  Recon-A vs 3x2PC:  %.6f (%.2f%%)\n", diff_A_3x2pc, 100.0 * diff_A_3x2pc / pc3x2.v2);
        printf("  Recon-B vs 3x2PC:  %.6f (%.2f%%)\n", diff_B_3x2pc, 100.0 * diff_B_3x2pc / pc3x2.v2);
        printf("  Recon-A vs Recon-B: %.6f (%.2f%%)\n\n", diff_A_B, 100.0 * diff_A_B / reconB.v2);
        
        v2_reconA_vals.push_back(reconA.v2);
        v2_reconA_errs.push_back(reconA.v2_err);
        v2_reconB_vals.push_back(reconB.v2);
        v2_reconB_errs.push_back(reconB.v2_err);
        v2_3x2pc_vals.push_back(pc3x2.v2);
        v2_3x2pc_errs.push_back(pc3x2.v2_err);
    }
    
    // Create comparison plot
    if (v2_reconA_vals.size() >= 2) {
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);
        
        TCanvas canvas("c_all3methods", "All 3 V2 Methods Comparison", 1400, 900);
        canvas.SetLeftMargin(0.12);
        canvas.SetRightMargin(0.05);
        canvas.SetBottomMargin(0.15);
        canvas.SetTopMargin(0.10);
        
        // Create graphs with x-offset for visibility
        TGraphErrors* g_reconA = new TGraphErrors(v2_reconA_vals.size());
        TGraphErrors* g_reconB = new TGraphErrors(v2_reconB_vals.size());
        TGraphErrors* g_3x2pc = new TGraphErrors(v2_3x2pc_vals.size());
        
        for (size_t i = 0; i < v2_reconA_vals.size(); i++) {
            Double_t x_A = (Double_t)i - 0.15;
            Double_t x_B = (Double_t)i;
            Double_t x_3x2pc = (Double_t)i + 0.15;
            
            g_reconA->SetPoint(i, x_A, v2_reconA_vals[i]);
            g_reconA->SetPointError(i, 0.0, v2_reconA_errs[i]);
            
            g_reconB->SetPoint(i, x_B, v2_reconB_vals[i]);
            g_reconB->SetPointError(i, 0.0, v2_reconB_errs[i]);
            
            g_3x2pc->SetPoint(i, x_3x2pc, v2_3x2pc_vals[i]);
            g_3x2pc->SetPointError(i, 0.0, v2_3x2pc_errs[i]);
        }
        
        // Style for Recon-A (red)
        g_reconA->SetMarkerStyle(20);
        g_reconA->SetMarkerSize(3.0);
        g_reconA->SetMarkerColor(kRed+1);
        g_reconA->SetLineColor(kRed+1);
        g_reconA->SetLineWidth(3);
        
        // Style for Recon-B (orange)
        g_reconB->SetMarkerStyle(21);
        g_reconB->SetMarkerSize(3.0);
        g_reconB->SetMarkerColor(kOrange+1);
        g_reconB->SetLineColor(kOrange+1);
        g_reconB->SetLineWidth(3);
        
        // Style for 3x2PC (blue)
        g_3x2pc->SetMarkerStyle(22);
        g_3x2pc->SetMarkerSize(3.3);
        g_3x2pc->SetMarkerColor(kBlue+1);
        g_3x2pc->SetLineColor(kBlue+1);
        g_3x2pc->SetLineWidth(3);
        
        // Axes
        g_reconA->SetTitle("");
        g_reconA->GetXaxis()->SetLimits(-0.5, 3.5);
        g_reconA->GetXaxis()->SetNdivisions(505);
        g_reconA->GetXaxis()->SetLabelSize(0.06);
        g_reconA->GetXaxis()->SetTitleSize(0.06);
        g_reconA->GetXaxis()->SetTitle("Dataset");
        
        // Y-axis range
        std::vector<Double_t> all_vals;
        all_vals.insert(all_vals.end(), v2_reconA_vals.begin(), v2_reconA_vals.end());
        all_vals.insert(all_vals.end(), v2_reconB_vals.begin(), v2_reconB_vals.end());
        all_vals.insert(all_vals.end(), v2_3x2pc_vals.begin(), v2_3x2pc_vals.end());
        
        Double_t min_val = *std::min_element(all_vals.begin(), all_vals.end());
        Double_t max_val = *std::max_element(all_vals.begin(), all_vals.end());
        Double_t center = (min_val + max_val) / 2.0;
        Double_t range = (max_val - min_val);
        Double_t y_min = center - range * 2.5;
        Double_t y_max = center + range * 2.5;
        
        g_reconA->GetYaxis()->SetRangeUser(y_min, y_max);
        g_reconA->GetYaxis()->SetLabelSize(0.06);
        g_reconA->GetYaxis()->SetTitleSize(0.06);
        g_reconA->GetYaxis()->SetTitleOffset(1.0);
        g_reconA->GetYaxis()->SetTitle("V_{2}");
        g_reconA->GetYaxis()->SetNdivisions(508);
        
        g_reconA->Draw("APE");
        g_reconB->Draw("PE same");
        g_3x2pc->Draw("PE same");
        
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
        latex.DrawLatex(0.14, 0.95, "ALICE Ne-Ne and O-O, 0-20% Centrality");
        
        latex.SetTextSize(0.045);
        latex.DrawLatex(0.14, 0.83, "Comparison of Three V_{2} Calculation Methods");
        
        // Legend on top - wider and slightly lower
        TLegend* leg = new TLegend(0.14, 0.64, 0.85, 0.78);
        leg->SetBorderSize(1);
        leg->SetFillColor(kWhite);
        leg->SetTextFont(42);
        leg->SetTextSize(0.032);
        leg->AddEntry(g_reconA, "Recon-A: V_{2#Delta}(FT0,TPC-A) / #sqrt{V_{2#Delta}(TPC,TPC)}", "pe");
        leg->AddEntry(g_reconB, "Recon-B: V_{2#Delta}(FT0,TPC-B) / #sqrt{V_{2#Delta}(TPC,TPC)}", "pe");
        leg->AddEntry(g_3x2pc, "3#times2PC: #sqrt{[V_{2#Delta}(FT0,TPC-A) #times V_{2#Delta}(FT0,TPC-B)] / V_{2#Delta}(TPC,TPC)}", "pe");
        leg->Draw();
        
        canvas.SaveAs("./TemplateFit/EtaDiff/All3Methods/V2_All3Methods_Comparison.pdf");
        
        std::cout << "\n=== Comparison Plot Created ===" << std::endl;
        std::cout << "Saved: ./TemplateFit/EtaDiff/All3Methods/V2_All3Methods_Comparison.pdf" << std::endl;
        
        // Create separate info plot with values table
        TCanvas canvas_info("c_info", "V2 Methods Values", 1200, 800);
        canvas_info.SetLeftMargin(0.05);
        canvas_info.SetRightMargin(0.05);
        canvas_info.SetBottomMargin(0.05);
        canvas_info.SetTopMargin(0.12);
        
        TPaveText *pt_info = new TPaveText(0.1, 0.1, 0.9, 0.85, "NDC");
        pt_info->SetFillColor(kWhite);
        pt_info->SetBorderSize(1);
        pt_info->SetTextAlign(12);
        pt_info->SetTextFont(42);
        pt_info->SetTextSize(0.03);
        
        pt_info->AddText("");
        pt_info->AddText("V_{2} Values Comparison");
        pt_info->AddText("");
        pt_info->AddText("========================================");
        pt_info->AddText("");
        
        // Add entries for all datasets
        for (size_t i = 0; i < datasets.size(); i++) {
            pt_info->AddText(Form("%s:", datasetLabels[i].c_str()));
            pt_info->AddText("");
            pt_info->AddText(Form("  Recon-A (TPC-A):  V_{2} = %.6f #pm %.6f", v2_reconA_vals[i], v2_reconA_errs[i]));
            pt_info->AddText(Form("  Recon-B (TPC-B):  V_{2} = %.6f #pm %.6f", v2_reconB_vals[i], v2_reconB_errs[i]));
            pt_info->AddText(Form("  3#times2PC:           V_{2} = %.6f #pm %.6f", v2_3x2pc_vals[i], v2_3x2pc_errs[i]));
            pt_info->AddText("");
            pt_info->AddText("========================================");
            pt_info->AddText("");
        }
        
        // Calculate average differences across all datasets
        Double_t avg_A_3x2pc = 0.0;
        Double_t avg_B_3x2pc = 0.0;
        Double_t avg_A_B = 0.0;
        
        for (size_t i = 0; i < v2_reconA_vals.size(); i++) {
            avg_A_3x2pc += 100.0 * (v2_reconA_vals[i] - v2_3x2pc_vals[i]) / v2_3x2pc_vals[i];
            avg_B_3x2pc += 100.0 * (v2_reconB_vals[i] - v2_3x2pc_vals[i]) / v2_3x2pc_vals[i];
            avg_A_B += 100.0 * (v2_reconA_vals[i] - v2_reconB_vals[i]) / v2_reconB_vals[i];
        }
        avg_A_3x2pc /= v2_reconA_vals.size();
        avg_B_3x2pc /= v2_reconB_vals.size();
        avg_A_B /= v2_reconA_vals.size();
        
        pt_info->AddText("Average Differences:");
        pt_info->AddText("");
        pt_info->AddText(Form("  Recon-A vs 3#times2PC:  %.2f%%", avg_A_3x2pc));
        pt_info->AddText(Form("  Recon-B vs 3#times2PC:  %.2f%%", avg_B_3x2pc));
        pt_info->AddText(Form("  Recon-A vs Recon-B:  %.2f%%", avg_A_B));
        pt_info->AddText("");
        
        pt_info->Draw();
        
        TLatex latex_info;
        latex_info.SetNDC();
        latex_info.SetTextFont(42);
        latex_info.SetTextSize(0.04);
        latex_info.DrawLatex(0.1, 0.95, "ALICE Ne-Ne and O-O, 0-20% Centrality");
        
        canvas_info.SaveAs("./TemplateFit/EtaDiff/All3Methods/V2_All3Methods_Values.pdf");
        
        std::cout << "Saved: ./TemplateFit/EtaDiff/All3Methods/V2_All3Methods_Values.pdf" << std::endl;
        
        delete pt_info;
        delete leg;
        
        // Create ratio plots
        TCanvas canvas2("c_ratios", "Method Ratios", 1400, 900);
        canvas2.SetLeftMargin(0.12);
        canvas2.SetRightMargin(0.05);
        canvas2.SetBottomMargin(0.15);
        canvas2.SetTopMargin(0.10);
        
        TGraphErrors* g_ratio_A_3x2pc = new TGraphErrors(v2_reconA_vals.size());
        TGraphErrors* g_ratio_B_3x2pc = new TGraphErrors(v2_reconB_vals.size());
        TGraphErrors* g_ratio_A_B = new TGraphErrors(v2_reconA_vals.size());
        
        for (size_t i = 0; i < v2_reconA_vals.size(); i++) {
            Double_t ratio_A_3x2pc = v2_reconA_vals[i] / v2_3x2pc_vals[i];
            Double_t ratio_B_3x2pc = v2_reconB_vals[i] / v2_3x2pc_vals[i];
            Double_t ratio_A_B = v2_reconA_vals[i] / v2_reconB_vals[i];
            
            // Simplified error propagation (relative errors add in quadrature)
            Double_t err_A_3x2pc = ratio_A_3x2pc * TMath::Sqrt(TMath::Power(v2_reconA_errs[i]/v2_reconA_vals[i], 2) + 
                                                                 TMath::Power(v2_3x2pc_errs[i]/v2_3x2pc_vals[i], 2));
            Double_t err_B_3x2pc = ratio_B_3x2pc * TMath::Sqrt(TMath::Power(v2_reconB_errs[i]/v2_reconB_vals[i], 2) + 
                                                                 TMath::Power(v2_3x2pc_errs[i]/v2_3x2pc_vals[i], 2));
            Double_t err_A_B = ratio_A_B * TMath::Sqrt(TMath::Power(v2_reconA_errs[i]/v2_reconA_vals[i], 2) + 
                                                         TMath::Power(v2_reconB_errs[i]/v2_reconB_vals[i], 2));
            
            g_ratio_A_3x2pc->SetPoint(i, (Double_t)i - 0.15, ratio_A_3x2pc);
            g_ratio_A_3x2pc->SetPointError(i, 0.0, err_A_3x2pc);
            
            g_ratio_B_3x2pc->SetPoint(i, (Double_t)i, ratio_B_3x2pc);
            g_ratio_B_3x2pc->SetPointError(i, 0.0, err_B_3x2pc);
            
            g_ratio_A_B->SetPoint(i, (Double_t)i + 0.15, ratio_A_B);
            g_ratio_A_B->SetPointError(i, 0.0, err_A_B);
        }
        
        // Style for ratios
        g_ratio_A_3x2pc->SetMarkerStyle(20);
        g_ratio_A_3x2pc->SetMarkerSize(3.0);
        g_ratio_A_3x2pc->SetMarkerColor(kViolet+1);
        g_ratio_A_3x2pc->SetLineColor(kViolet+1);
        g_ratio_A_3x2pc->SetLineWidth(3);
        
        g_ratio_B_3x2pc->SetMarkerStyle(21);
        g_ratio_B_3x2pc->SetMarkerSize(3.0);
        g_ratio_B_3x2pc->SetMarkerColor(kGreen+2);
        g_ratio_B_3x2pc->SetLineColor(kGreen+2);
        g_ratio_B_3x2pc->SetLineWidth(3);
        
        g_ratio_A_B->SetMarkerStyle(22);
        g_ratio_A_B->SetMarkerSize(3.3);
        g_ratio_A_B->SetMarkerColor(kAzure+2);
        g_ratio_A_B->SetLineColor(kAzure+2);
        g_ratio_A_B->SetLineWidth(3);
        
        g_ratio_A_3x2pc->SetTitle("");
        g_ratio_A_3x2pc->GetXaxis()->SetLimits(-0.5, 3.5);
        g_ratio_A_3x2pc->GetXaxis()->SetNdivisions(505);
        g_ratio_A_3x2pc->GetXaxis()->SetLabelSize(0.06);
        g_ratio_A_3x2pc->GetXaxis()->SetTitleSize(0.06);
        g_ratio_A_3x2pc->GetXaxis()->SetTitle("Dataset");
        
        g_ratio_A_3x2pc->GetYaxis()->SetRangeUser(0.95, 1.15);
        g_ratio_A_3x2pc->GetYaxis()->SetLabelSize(0.06);
        g_ratio_A_3x2pc->GetYaxis()->SetTitleSize(0.06);
        g_ratio_A_3x2pc->GetYaxis()->SetTitleOffset(1.0);
        g_ratio_A_3x2pc->GetYaxis()->SetTitle("Ratio");
        g_ratio_A_3x2pc->GetYaxis()->SetNdivisions(508);
        
        g_ratio_A_3x2pc->Draw("APE");
        g_ratio_B_3x2pc->Draw("PE same");
        g_ratio_A_B->Draw("PE same");
        
        // Unity line
        TLine* unity = new TLine(-0.5, 1.0, 3.5, 1.0);
        unity->SetLineStyle(2);
        unity->SetLineWidth(2);
        unity->SetLineColor(kGray+2);
        unity->Draw();
        
        xlab.DrawLatex(0.0, 0.93, "615818 Ne-Ne Outer");
        xlab.DrawLatex(1.0, 0.93, "615817 Ne-Ne Inner");
        xlab.DrawLatex(2.0, 0.93, "616549 O-O Outer");
        xlab.DrawLatex(3.0, 0.93, "618685 O-O Inner");
        
        latex.DrawLatex(0.14, 0.95, "ALICE Ne-Ne and O-O, 0-20% Centrality");
        latex.DrawLatex(0.14, 0.89, "Method Ratios");
        
        TLegend* leg2 = new TLegend(0.55, 0.73, 0.91, 0.87);
        leg2->SetBorderSize(1);
        leg2->SetFillColor(kWhite);
        leg2->SetTextFont(42);
        leg2->SetTextSize(0.04);
        leg2->AddEntry(g_ratio_A_3x2pc, "Recon-A / 3#times2PC", "pe");
        leg2->AddEntry(g_ratio_B_3x2pc, "Recon-B / 3#times2PC", "pe");
        leg2->AddEntry(g_ratio_A_B, "Recon-A / Recon-B", "pe");
        leg2->Draw();
        
        canvas2.SaveAs("./TemplateFit/EtaDiff/All3Methods/V2_All3Methods_Ratios.pdf");
        
        std::cout << "Saved: ./TemplateFit/EtaDiff/All3Methods/V2_All3Methods_Ratios.pdf" << std::endl;
        
        delete unity;
        delete leg2;
        delete g_ratio_A_3x2pc;
        delete g_ratio_B_3x2pc;
        delete g_ratio_A_B;
        delete g_reconA;
        delete g_reconB;
        delete g_3x2pc;
    }
    
    std::cout << "\n=== All 3 Methods Comparison Complete ===" << std::endl;
}
