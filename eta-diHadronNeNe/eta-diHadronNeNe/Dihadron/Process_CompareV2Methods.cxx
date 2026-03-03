#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLine.h"

//==============================================================
// Compare V2 from different methods:
// 1. ReconstructV2: V2 = V2Δ(FT0,TPC-A) / sqrt(V2Δ(TPC-A,TPC-B))
// 2. 3x2PC: V2 = sqrt[(V2Δ(FT0,TPC-A) × V2Δ(FT0,TPC-B)) / V2Δ(TPC-A,TPC-B)]
//==============================================================

struct V2Result {
    Double_t v2;
    Double_t v2_err;
    std::string label;
};

V2Result ReadV2FromFile(TString filename, TString histname = "hV2") {
    V2Result result = {0, 0, ""};
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Warning: Could not open file " << filename << std::endl;
        return result;
    }
    
    TH1D* h = (TH1D*)file->Get(histname);
    if (!h) {
        // Try alternative histogram names
        h = (TH1D*)file->Get("hV2Reconstructed");
    }
    if (!h) {
        h = (TH1D*)file->Get("hV2_3x2PC");
    }
    
    if (h) {
        result.v2 = h->GetBinContent(1);
        result.v2_err = h->GetBinError(1);
    } else {
        std::cerr << "Warning: Could not find histogram in " << filename << std::endl;
    }
    
    file->Close();
    delete file;
    return result;
}

void Process_CompareV2Methods() {
    std::vector<std::string> datasets = {"615818", "615817"};
    std::vector<std::string> datasetLabels = {"615818 (Outer)", "615817 (Inner)"};
    
    std::cout << "\n=== Comparing V2 Methods ===" << std::endl;
    
    // Storage for results
    std::vector<Double_t> v2_recon_vals, v2_recon_errs;
    std::vector<Double_t> v2_3x2pc_vals, v2_3x2pc_errs;
    
    // Read results from both methods
    for (const auto& dataset : datasets) {
        // Read ReconstructV2 result
        TString recon_file = Form("./TemplateFit/EtaDiff/ReconstructedV2/ReconstructedV2_LHC25af_pass2_%s_TPC_FT0C_Cent_0_20.root", 
                                  dataset.c_str());
        V2Result recon = ReadV2FromFile(recon_file, "hV2Reconstructed");
        
        // Read 3x2PC result
        TString pc3x2_file = Form("./TemplateFit/EtaDiff/V2_3x2PC/V2_3x2PC_LHC25af_pass2_%s_TPC_FT0C_Cent_0_20.root", 
                                  dataset.c_str());
        V2Result pc3x2 = ReadV2FromFile(pc3x2_file, "hV2_3x2PC");
        
        v2_recon_vals.push_back(recon.v2);
        v2_recon_errs.push_back(recon.v2_err);
        v2_3x2pc_vals.push_back(pc3x2.v2);
        v2_3x2pc_errs.push_back(pc3x2.v2_err);
        
        printf("\nDataset %s:\n", dataset.c_str());
        printf("  ReconstructV2: %.6f +/- %.6e\n", recon.v2, recon.v2_err);
        printf("  3x2PC:         %.6f +/- %.6e\n", pc3x2.v2, pc3x2.v2_err);
        printf("  Difference:    %.6f (%.2f%%)\n", 
               recon.v2 - pc3x2.v2, 
               100.0 * (recon.v2 - pc3x2.v2) / pc3x2.v2);
    }
    
    // Create comparison plot
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
    
    TCanvas canvas("c_compare", "V2 Method Comparison", 1400, 900);
    canvas.SetLeftMargin(0.12);
    canvas.SetRightMargin(0.05);
    canvas.SetBottomMargin(0.15);
    canvas.SetTopMargin(0.10);
    
    // Create graphs for each method
    TGraphErrors* g_recon = new TGraphErrors(datasets.size());
    TGraphErrors* g_3x2pc = new TGraphErrors(datasets.size());
    
    for (size_t i = 0; i < datasets.size(); i++) {
        // Offset x-positions slightly for visibility
        Double_t x_recon = (Double_t)i - 0.1;
        Double_t x_3x2pc = (Double_t)i + 0.1;
        
        g_recon->SetPoint(i, x_recon, v2_recon_vals[i]);
        g_recon->SetPointError(i, 0.0, v2_recon_errs[i]);
        
        g_3x2pc->SetPoint(i, x_3x2pc, v2_3x2pc_vals[i]);
        g_3x2pc->SetPointError(i, 0.0, v2_3x2pc_errs[i]);
    }
    
    // Style for ReconstructV2
    g_recon->SetMarkerStyle(20);
    g_recon->SetMarkerSize(3.0);
    g_recon->SetMarkerColor(kRed+1);
    g_recon->SetLineColor(kRed+1);
    g_recon->SetLineWidth(3);
    
    // Style for 3x2PC
    g_3x2pc->SetMarkerStyle(21);
    g_3x2pc->SetMarkerSize(3.0);
    g_3x2pc->SetMarkerColor(kBlue+1);
    g_3x2pc->SetLineColor(kBlue+1);
    g_3x2pc->SetLineWidth(3);
    
    // Set up axes
    g_recon->SetTitle("");
    g_recon->GetXaxis()->SetLimits(-0.5, 1.5);
    g_recon->GetXaxis()->SetNdivisions(505);
    g_recon->GetXaxis()->SetLabelSize(0.06);
    g_recon->GetXaxis()->SetTitleSize(0.06);
    g_recon->GetXaxis()->SetLabelFont(42);
    g_recon->GetXaxis()->SetTitleFont(42);
    g_recon->GetXaxis()->SetTitle("Dataset");
    
    // Set y-axis range to show both methods
    Double_t all_vals[] = {v2_recon_vals[0], v2_recon_vals[1], v2_3x2pc_vals[0], v2_3x2pc_vals[1]};
    Double_t min_val = *std::min_element(all_vals, all_vals + 4);
    Double_t max_val = *std::max_element(all_vals, all_vals + 4);
    Double_t center = (min_val + max_val) / 2.0;
    Double_t range = (max_val - min_val);
    Double_t y_min = center - range * 2.5;
    Double_t y_max = center + range * 2.5;
    
    g_recon->GetYaxis()->SetRangeUser(y_min, y_max);
    g_recon->GetYaxis()->SetLabelSize(0.06);
    g_recon->GetYaxis()->SetTitleSize(0.06);
    g_recon->GetYaxis()->SetTitleOffset(1.0);
    g_recon->GetYaxis()->SetLabelFont(42);
    g_recon->GetYaxis()->SetTitleFont(42);
    g_recon->GetYaxis()->SetTitle("V_{2}");
    g_recon->GetYaxis()->SetNdivisions(508);
    
    // Draw
    g_recon->Draw("APE");
    g_3x2pc->Draw("PE same");
    
    // Add x-axis labels
    TLatex xlab;
    xlab.SetTextSize(0.055);
    xlab.SetTextFont(42);
    xlab.SetTextAlign(22);
    Double_t label_y = y_min - (y_max - y_min) * 0.08;
    xlab.DrawLatex(0.0, label_y, "615818 Outer");
    xlab.DrawLatex(1.0, label_y, "615817 Inner");
    
    // Add title
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    
    latex.SetTextSize(0.05);
    latex.DrawLatex(0.14, 0.95, "ALICE Ne-Ne #sqrt{s_{NN}} = 17.3 GeV, Centrality 0-20%");
    
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.14, 0.89, "Comparison of V_{2} Calculation Methods");
    
    // Add legend
    TLegend* leg = new TLegend(0.14, 0.75, 0.50, 0.86);
    leg->SetBorderSize(1);
    leg->SetFillColor(kWhite);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->AddEntry(g_recon, "ReconstructV2 Method", "pe");
    leg->AddEntry(g_3x2pc, "3x2PC Method", "pe");
    leg->Draw();
    
    // Add information box
    TPaveText *pt = new TPaveText(0.55, 0.65, 0.93, 0.86, "NDC");
    pt->SetFillColor(kWhite);
    pt->SetBorderSize(1);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.035);
    pt->AddText("ReconstructV2:");
    pt->AddText(Form("  Outer: %.5f #pm %.5f", v2_recon_vals[0], v2_recon_errs[0]));
    pt->AddText(Form("  Inner: %.5f #pm %.5f", v2_recon_vals[1], v2_recon_errs[1]));
    pt->AddText("");
    pt->AddText("3x2PC:");
    pt->AddText(Form("  Outer: %.5f #pm %.5f", v2_3x2pc_vals[0], v2_3x2pc_errs[0]));
    pt->AddText(Form("  Inner: %.5f #pm %.5f", v2_3x2pc_vals[1], v2_3x2pc_errs[1]));
    pt->AddText("");
    Double_t avg_diff = 100.0 * ((v2_recon_vals[0] - v2_3x2pc_vals[0]) / v2_3x2pc_vals[0] + 
                                  (v2_recon_vals[1] - v2_3x2pc_vals[1]) / v2_3x2pc_vals[1]) / 2.0;
    pt->AddText(Form("Avg. Difference: %.1f%%", avg_diff));
    pt->Draw();
    
    // Save plot
    gSystem->Exec("mkdir -p ./TemplateFit/EtaDiff/MethodComparison");
    canvas.SaveAs("./TemplateFit/EtaDiff/MethodComparison/V2_Methods_Comparison.pdf");
    
    std::cout << "\n=== Comparison Plot Created ===" << std::endl;
    std::cout << "Saved: ./TemplateFit/EtaDiff/MethodComparison/V2_Methods_Comparison.pdf" << std::endl;
    std::cout << "Average difference between methods: " << avg_diff << "%" << std::endl;
    
    // Create ratio plot
    TCanvas c_ratio("c_ratio", "Method Ratio", 1200, 900);
    c_ratio.SetLeftMargin(0.15);
    c_ratio.SetRightMargin(0.05);
    c_ratio.SetBottomMargin(0.15);
    c_ratio.SetTopMargin(0.10);
    
    TGraphErrors* g_ratio = new TGraphErrors(datasets.size());
    
    for (size_t i = 0; i < datasets.size(); i++) {
        Double_t ratio = v2_recon_vals[i] / v2_3x2pc_vals[i];
        Double_t ratio_err = ratio * TMath::Sqrt(
            TMath::Power(v2_recon_errs[i] / v2_recon_vals[i], 2) +
            TMath::Power(v2_3x2pc_errs[i] / v2_3x2pc_vals[i], 2)
        );
        
        g_ratio->SetPoint(i, (Double_t)i, ratio);
        g_ratio->SetPointError(i, 0.0, ratio_err);
    }
    
    g_ratio->SetTitle("");
    g_ratio->SetMarkerStyle(20);
    g_ratio->SetMarkerSize(3.5);
    g_ratio->SetMarkerColor(kViolet+1);
    g_ratio->SetLineColor(kViolet+1);
    g_ratio->SetLineWidth(3);
    
    g_ratio->GetXaxis()->SetLimits(-0.5, 1.5);
    g_ratio->GetXaxis()->SetNdivisions(505);
    g_ratio->GetXaxis()->SetLabelSize(0.06);
    g_ratio->GetXaxis()->SetTitleSize(0.06);
    g_ratio->GetXaxis()->SetTitle("Dataset");
    
    g_ratio->GetYaxis()->SetRangeUser(0.98, 1.12);
    g_ratio->GetYaxis()->SetLabelSize(0.06);
    g_ratio->GetYaxis()->SetTitleSize(0.06);
    g_ratio->GetYaxis()->SetTitleOffset(1.3);
    g_ratio->GetYaxis()->SetTitle("V_{2}^{Recon} / V_{2}^{3x2PC}");
    g_ratio->GetYaxis()->SetNdivisions(508);
    
    g_ratio->Draw("APE");
    
    // Draw unity line
    TLine unity(-0.5, 1.0, 1.5, 1.0);
    unity.SetLineStyle(2);
    unity.SetLineColor(kGray+2);
    unity.SetLineWidth(2);
    unity.Draw("same");
    g_ratio->Draw("PE same");
    
    // Add labels
    xlab.DrawLatex(0.0, 0.97, "615818 Outer");
    xlab.DrawLatex(1.0, 0.97, "615817 Inner");
    
    latex.DrawLatex(0.18, 0.95, "ALICE Ne-Ne #sqrt{s_{NN}} = 17.3 GeV, 0-20% Centrality");
    
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.18, 0.88, "Ratio: ReconstructV2 / 3x2PC");
    
    // Add value box
    TPaveText *pt_ratio = new TPaveText(0.55, 0.72, 0.92, 0.84, "NDC");
    pt_ratio->SetFillColor(kWhite);
    pt_ratio->SetBorderSize(1);
    pt_ratio->SetTextAlign(12);
    pt_ratio->SetTextFont(42);
    pt_ratio->SetTextSize(0.04);
    
    Double_t ratio_outer, ratio_inner;
    g_ratio->GetPoint(0, ratio_outer, ratio_outer);
    g_ratio->GetPoint(1, ratio_inner, ratio_inner);
    
    pt_ratio->AddText(Form("Outer (615818): %.4f", ratio_outer));
    pt_ratio->AddText(Form("Inner (615817): %.4f", ratio_inner));
    pt_ratio->Draw();
    
    c_ratio.SaveAs("./TemplateFit/EtaDiff/MethodComparison/V2_Methods_Ratio.pdf");
    
    std::cout << "Saved: ./TemplateFit/EtaDiff/MethodComparison/V2_Methods_Ratio.pdf" << std::endl;
    std::cout << "\n=== Method Comparison Complete ===" << std::endl;
    
    delete pt;
    delete pt_ratio;
    delete leg;
    delete g_recon;
    delete g_3x2pc;
    delete g_ratio;
}
