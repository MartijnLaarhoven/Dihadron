#include "TROOT.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TMultiGraph.h"
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <limits>

// Create comprehensive comparison plots showing all three methods
void Plot_AllMethods_Comparison() {
    gROOT->SetBatch(kTRUE);
    gSystem->mkdir("./TemplateFit/Comparisons", kTRUE);
    gSystem->mkdir("./TemplateFit/Comparisons/PDFs", kTRUE);

    // Dataset names
    std::string ooNeg = "LHC25ae_pass2_604826";
    std::string ooPos = "LHC25ae_pass2_604830";
    std::string neNeg = "LHC25af_pass2_611697";
    std::string nePos = "LHC25af_pass2_604820";

    std::cout << "========================================" << std::endl;
    std::cout << "CREATING COMPREHENSIVE METHOD COMPARISONS" << std::endl;
    std::cout << "========================================" << std::endl;

    // ============================================================
    // PLOT 1: O-O Combined - All Three Methods
    // ============================================================
    std::cout << "\n1. Creating O-O combined comparison (all methods)..." << std::endl;
    
    // Load Template Fit O-O
    TFile* fTF_OO = TFile::Open(Form("./TemplateFit/EtaDiff/Vn_Combined_%s_vs_%s_Cent.root", 
                                      ooNeg.c_str(), ooPos.c_str()), "READ");
    TGraphErrors* gTF_OO = nullptr;
    if (fTF_OO && fTF_OO->IsOpen()) {
        gTF_OO = (TGraphErrors*)fTF_OO->Get("gV2Delta_Combined");
        if (gTF_OO) {
            gTF_OO = (TGraphErrors*)gTF_OO->Clone("gTF_OO");
            std::cout << "  ✓ Template Fit O-O: " << gTF_OO->GetN() << " points" << std::endl;
        }
    } else {
        std::cerr << "  ✗ Template Fit O-O: File not found or cannot open" << std::endl;
    }
    
    // Load ZYAM O-O
    TFile* fZYAM_OO = TFile::Open(Form("./TemplateFit/PeripheralSubtracted/EtaDiff/Vn_Combined_%s_vs_%s_Cent.root",
                                        ooNeg.c_str(), ooPos.c_str()), "READ");
    TGraphErrors* gZYAM_OO = nullptr;
    if (fZYAM_OO && fZYAM_OO->IsOpen()) {
        gZYAM_OO = (TGraphErrors*)fZYAM_OO->Get("gV2Delta_Combined");
        if (gZYAM_OO) {
            gZYAM_OO = (TGraphErrors*)gZYAM_OO->Clone("gZYAM_OO");
            std::cout << "  ✓ ZYAM O-O: " << gZYAM_OO->GetN() << " points" << std::endl;
        }
    }
    
    // NOTE: LM-HM removed from O-O comparison plots for clarity
    
    // Load Uncorrected O-O
    TFile* fUncorr_OO = TFile::Open(Form("./TemplateFit/Uncorrected_FourierFit/EtaDiff/Vn_Combined_%s_vs_%s_Cent.root",
                                         ooNeg.c_str(), ooPos.c_str()), "READ");
    TGraphErrors* gUncorr_OO = nullptr;
    if (fUncorr_OO && fUncorr_OO->IsOpen()) {
        gUncorr_OO = (TGraphErrors*)fUncorr_OO->Get("gV2_Uncorrected_Combined");
        if (gUncorr_OO) {
            gUncorr_OO = (TGraphErrors*)gUncorr_OO->Clone("gUncorr_OO");
            std::cout << "  ✓ Uncorrected O-O: " << gUncorr_OO->GetN() << " points" << std::endl;
        }
    }
    
    if (gTF_OO && gZYAM_OO && gUncorr_OO) {
        TCanvas* cOO = new TCanvas("cOO_AllMethods", "O-O All Methods", 1000, 700);
        
        // Style settings
        gTF_OO->SetMarkerStyle(20);
        gTF_OO->SetMarkerColor(kRed);
        gTF_OO->SetMarkerSize(1.2);
        
        gZYAM_OO->SetMarkerStyle(21);
        gZYAM_OO->SetMarkerColor(kBlue);
        gZYAM_OO->SetMarkerSize(1.2);
        
        gUncorr_OO->SetMarkerStyle(23);
        gUncorr_OO->SetMarkerColor(kMagenta);
        gUncorr_OO->SetMarkerSize(1.2);
        
        // Draw - start with graph with most points to set proper axis ranges
        gZYAM_OO->GetXaxis()->SetTitle("#eta_{trig}");
        gZYAM_OO->GetYaxis()->SetTitle("v_{2}");
        gZYAM_OO->GetXaxis()->SetLimits(-0.8, 0.8);
        gZYAM_OO->GetYaxis()->SetRangeUser(0.0, 0.012);
        
        // Ensure absolutely no lines are drawn - set BEFORE drawing
        gZYAM_OO->SetLineColor(kWhite);
        gZYAM_OO->SetLineWidth(0);
        gTF_OO->SetLineColor(kWhite);
        gTF_OO->SetLineWidth(0);
        gUncorr_OO->SetLineColor(kWhite);
        gUncorr_OO->SetLineWidth(0);
        
        gZYAM_OO->Draw("AP");
        gTF_OO->Draw("P SAME");
        gUncorr_OO->Draw("P SAME");
        
        TLegend* legOO = new TLegend(0.15, 0.65, 0.50, 0.88);
        legOO->SetHeader("O-O Collisions", "C");
        legOO->AddEntry(gUncorr_OO, "Uncorrected (Raw)", "p");
        legOO->AddEntry(gTF_OO, "Template Fit", "p");
        legOO->AddEntry(gZYAM_OO, "ZYAM (Peripheral subtraction)", "p");
        legOO->SetBorderSize(0);
        legOO->SetFillStyle(0);
        legOO->Draw();
        
        cOO->SaveAs("./TemplateFit/Comparisons/PDFs/AllMethods_OO_Comparison.pdf");
        std::cout << "  → Saved: AllMethods_OO_Comparison.pdf" << std::endl;
        delete cOO;
        delete legOO;
    } else {
        std::cout << "  ✗ Missing some O-O data" << std::endl;
    }
    
    // ============================================================
    // PLOT 2: Ne-Ne Combined - All Three Methods
    // ============================================================
    std::cout << "\n2. Creating Ne-Ne combined comparison (all methods)..." << std::endl;
    
    // Load Template Fit Ne-Ne
    TFile* fTF_Ne = TFile::Open(Form("./TemplateFit/EtaDiff/Vn_Combined_%s_vs_%s_Cent.root", 
                                      neNeg.c_str(), nePos.c_str()), "READ");
    TGraphErrors* gTF_Ne = nullptr;
    if (fTF_Ne && fTF_Ne->IsOpen()) {
        gTF_Ne = (TGraphErrors*)fTF_Ne->Get("gV2Delta_Combined");
        if (gTF_Ne) {
            gTF_Ne = (TGraphErrors*)gTF_Ne->Clone("gTF_Ne");
            std::cout << "  ✓ Template Fit Ne-Ne: " << gTF_Ne->GetN() << " points" << std::endl;
        }
    }
    
    // Load ZYAM Ne-Ne
    TFile* fZYAM_Ne = TFile::Open(Form("./TemplateFit/PeripheralSubtracted/EtaDiff/Vn_Combined_%s_vs_%s_Cent.root",
                                        neNeg.c_str(), nePos.c_str()), "READ");
    TGraphErrors* gZYAM_Ne = nullptr;
    if (fZYAM_Ne && fZYAM_Ne->IsOpen()) {
        gZYAM_Ne = (TGraphErrors*)fZYAM_Ne->Get("gV2Delta_Combined");
        if (gZYAM_Ne) {
            gZYAM_Ne = (TGraphErrors*)gZYAM_Ne->Clone("gZYAM_Ne");
            std::cout << "  ✓ ZYAM Ne-Ne: " << gZYAM_Ne->GetN() << " points" << std::endl;
        }
    }
    
    // Load Uncorrected Ne-Ne
    TFile* fUncorr_Ne = TFile::Open(Form("./TemplateFit/Uncorrected_FourierFit/EtaDiff/Vn_Combined_%s_vs_%s_Cent.root",
                                         neNeg.c_str(), nePos.c_str()), "READ");
    TGraphErrors* gUncorr_Ne = nullptr;
    if (fUncorr_Ne && fUncorr_Ne->IsOpen()) {
        gUncorr_Ne = (TGraphErrors*)fUncorr_Ne->Get("gV2_Uncorrected_Combined");
        if (gUncorr_Ne) {
            gUncorr_Ne = (TGraphErrors*)gUncorr_Ne->Clone("gUncorr_Ne");
            std::cout << "  ✓ Uncorrected Ne-Ne: " << gUncorr_Ne->GetN() << " points" << std::endl;
        }
    }
    
    if (gTF_Ne && gZYAM_Ne && gUncorr_Ne) {
        TCanvas* cNe = new TCanvas("cNe_AllMethods", "Ne-Ne All Methods", 1000, 700);
        
        // Style settings
        gTF_Ne->SetMarkerStyle(20);
        gTF_Ne->SetMarkerColor(kRed);
        gTF_Ne->SetMarkerSize(1.2);
        
        gZYAM_Ne->SetMarkerStyle(21);
        gZYAM_Ne->SetMarkerColor(kBlue);
        gZYAM_Ne->SetMarkerSize(1.2);
        
        gUncorr_Ne->SetMarkerStyle(23);
        gUncorr_Ne->SetMarkerColor(kMagenta);
        gUncorr_Ne->SetMarkerSize(1.2);
        
        // Draw - start with graph with most points to set proper axis ranges
        gZYAM_Ne->GetXaxis()->SetTitle("#eta_{trig}");
        gZYAM_Ne->GetYaxis()->SetTitle("v_{2}");
        gZYAM_Ne->GetXaxis()->SetLimits(-0.8, 0.8);
        gZYAM_Ne->GetYaxis()->SetRangeUser(0.0, 0.012);
        
        // Ensure absolutely no lines are drawn - set BEFORE drawing
        gZYAM_Ne->SetLineColor(kWhite);
        gZYAM_Ne->SetLineWidth(0);
        gTF_Ne->SetLineColor(kWhite);
        gTF_Ne->SetLineWidth(0);
        gUncorr_Ne->SetLineColor(kWhite);
        gUncorr_Ne->SetLineWidth(0);
        
        gZYAM_Ne->Draw("AP");
        gTF_Ne->Draw("P SAME");
        gUncorr_Ne->Draw("P SAME");
        
        TLegend* legNe = new TLegend(0.15, 0.65, 0.50, 0.88);
        legNe->SetHeader("Ne-Ne Collisions", "C");
        legNe->AddEntry(gUncorr_Ne, "Uncorrected (Raw)", "p");
        legNe->AddEntry(gTF_Ne, "Template Fit", "p");
        legNe->AddEntry(gZYAM_Ne, "ZYAM (Peripheral subtraction)", "p");
        legNe->SetBorderSize(0);
        legNe->SetFillStyle(0);
        legNe->Draw();
        
        cNe->SaveAs("./TemplateFit/Comparisons/PDFs/AllMethods_NeNe_Comparison.pdf");
        std::cout << "  → Saved: AllMethods_NeNe_Comparison.pdf" << std::endl;
        delete cNe;
        delete legNe;
    } else {
        std::cout << "  ✗ Missing some Ne-Ne data" << std::endl;
    }
    
    // ============================================================
    // PLOT 3: O-O vs Ne-Ne - Template Fit Only
    // ============================================================
    std::cout << "\n3. Creating O-O vs Ne-Ne (Template Fit)..." << std::endl;
    if (gTF_OO && gTF_Ne) {
        TCanvas* cTF = new TCanvas("cTF_Systems", "Template Fit Systems", 1000, 700);
        
        gTF_OO->SetMarkerColor(kRed);
        gTF_OO->SetLineColor(kRed);
        gTF_Ne->SetMarkerColor(kBlue);
        gTF_Ne->SetLineColor(kBlue);
        
        gTF_OO->GetXaxis()->SetTitle("#eta_{trig}");
        gTF_OO->GetYaxis()->SetTitle("v_{2#Delta}");
        gTF_OO->GetXaxis()->SetLimits(-0.8, 0.8);
        gTF_OO->GetYaxis()->SetRangeUser(0.0, 0.010);
        gTF_OO->Draw("AP");
        gTF_Ne->Draw("P SAME");
        
        TLegend* legTF = new TLegend(0.6, 0.75, 0.88, 0.88);
        legTF->SetHeader("Template Fit", "C");
        legTF->AddEntry(gTF_OO, "O-O", "p");
        legTF->AddEntry(gTF_Ne, "Ne-Ne", "p");
        legTF->SetBorderSize(0);
        legTF->SetFillStyle(0);
        legTF->Draw();
        
        cTF->SaveAs("./TemplateFit/Comparisons/PDFs/TemplateFit_OO_vs_NeNe.pdf");
        std::cout << "  → Saved: TemplateFit_OO_vs_NeNe.pdf" << std::endl;
        delete cTF;
        delete legTF;
    }
    
    // ============================================================
    // PLOT 4: O-O vs Ne-Ne - ZYAM Only
    // ============================================================
    std::cout << "\n4. Creating O-O vs Ne-Ne (ZYAM)..." << std::endl;
    if (gZYAM_OO && gZYAM_Ne) {
        TCanvas* cZYAM = new TCanvas("cZYAM_Systems", "ZYAM Systems", 1000, 700);
        
        gZYAM_OO->SetMarkerColor(kRed);
        gZYAM_OO->SetLineColor(kRed);
        gZYAM_Ne->SetMarkerColor(kBlue);
        gZYAM_Ne->SetLineColor(kBlue);
        
        gZYAM_OO->GetXaxis()->SetTitle("#eta_{trig}");
        gZYAM_OO->GetYaxis()->SetTitle("v_{2#Delta}");
        gZYAM_OO->GetXaxis()->SetLimits(-0.8, 0.8);
        gZYAM_OO->GetYaxis()->SetRangeUser(0.0, 0.010);
        gZYAM_OO->Draw("AP");
        gZYAM_Ne->Draw("P SAME");
        
        TLegend* legZYAM = new TLegend(0.6, 0.75, 0.88, 0.88);
        legZYAM->SetHeader("ZYAM Method", "C");
        legZYAM->AddEntry(gZYAM_OO, "O-O", "p");
        legZYAM->AddEntry(gZYAM_Ne, "Ne-Ne", "p");
        legZYAM->SetBorderSize(0);
        legZYAM->SetFillStyle(0);
        legZYAM->Draw();
        
        cZYAM->SaveAs("./TemplateFit/Comparisons/PDFs/ZYAM_OO_vs_NeNe.pdf");
        std::cout << "  → Saved: ZYAM_OO_vs_NeNe.pdf" << std::endl;
        delete cZYAM;
        delete legZYAM;
    }
    
    // Cleanup
    if (fTF_OO) { fTF_OO->Close(); delete fTF_OO; }
    if (fZYAM_OO) { fZYAM_OO->Close(); delete fZYAM_OO; }
    if (fTF_Ne) { fTF_Ne->Close(); delete fTF_Ne; }
    if (fZYAM_Ne) { fZYAM_Ne->Close(); delete fZYAM_Ne; }
    if (fUncorr_OO) { fUncorr_OO->Close(); delete fUncorr_OO; }
    if (fUncorr_Ne) { fUncorr_Ne->Close(); delete fUncorr_Ne; }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "ALL COMPARISON PLOTS COMPLETE" << std::endl;
    std::cout << "Output directory: ./TemplateFit/Comparisons/PDFs/" << std::endl;
    std::cout << "========================================" << std::endl;
}
