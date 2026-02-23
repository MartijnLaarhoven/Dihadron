/*
 * Compare Template Fit vs Peripheral Subtraction Methods
 */
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TROOT.h"
#include <iostream>
#include <string>
#include <vector>

void CompareMethodsForDataset(const std::string& dataset, const std::string& label, const std::string& splitName);
void CompareMethods_OO();
void CompareMethods_OO_Combined();

void Compare_Methods() {
    gROOT->SetBatch(kTRUE);
    
    CompareMethods_OO();
    CompareMethods_OO_Combined();
    
    std::cout << "\nComparison plots created in ./TemplateFit/Comparisons/PDFs/" << std::endl;
}

void CompareMethods_OO() {
    const std::string splitName = "Cent";
    
    // O-O standard (negative eta)
    CompareMethodsForDataset("LHC25ae_pass2_604826", "O-O (negative eta)", splitName);
    
    // O-O reversed (positive eta)
    CompareMethodsForDataset("LHC25ae_pass2_604830", "O-O (positive eta)", splitName);
    
    // Ne-Ne reversed (if available)
    CompareMethodsForDataset("LHC25af_pass2_604820", "Ne-Ne (positive eta)", splitName);
}

void CompareMethodsForDataset(const std::string& dataset, const std::string& label, const std::string& splitName) {
    std::cout << "\n========== Comparing methods for " << label << " ==========" << std::endl;
    
    // Try to open template fit results
    std::string tfFile = Form("./TemplateFit/EtaDiff/Vn_%s_%s_0_20.root", dataset.c_str(), splitName.c_str());
    TFile* fTF = TFile::Open(tfFile.c_str(), "READ");
    
    // Try to open peripheral subtraction results
    std::string psFile = Form("./TemplateFit/PeripheralSubtracted/EtaDiff/Vn_%s_PeripheralSubtracted_%s.root", 
                              dataset.c_str(), splitName.c_str());
    TFile* fPS = TFile::Open(psFile.c_str(), "READ");
    
    if (!fTF || !fTF->IsOpen()) {
        std::cerr << "Template fit file not found: " << tfFile << std::endl;
        if (fTF) { fTF->Close(); delete fTF; }
        if (fPS && fPS->IsOpen()) { fPS->Close(); delete fPS; }
        return;
    }
    
    if (!fPS || !fPS->IsOpen()) {
        std::cerr << "Peripheral subtraction file not found: " << psFile << std::endl;
        if (fTF) { fTF->Close(); delete fTF; }
        if (fPS) { fPS->Close(); delete fPS; }
        return;
    }
    
    // Get graphs
    TGraphErrors* gTF = (TGraphErrors*)fTF->Get("gV2Delta");
    TGraphErrors* gPS = (TGraphErrors*)fPS->Get("gV2Delta");
    
    if (!gTF || !gPS) {
        std::cerr << "Missing graphs in one of the files." << std::endl;
        fTF->Close();
        fPS->Close();
        delete fTF;
        delete fPS;
        return;
    }
    
    // Clone graphs to avoid deletion issues
    TGraphErrors* gTF_copy = (TGraphErrors*)gTF->Clone("gTF_copy");
    TGraphErrors* gPS_copy = (TGraphErrors*)gPS->Clone("gPS_copy");
    
    // Find y-axis range (include negative values)
    double minY = 1e10;
    double maxY = -1e10;
    
    auto updateMinMax = [&](TGraphErrors* g) {
        for (int i = 0; i < g->GetN(); i++) {
            double x, y;
            g->GetPoint(i, x, y);
            double ey = g->GetErrorY(i);
            minY = std::min(minY, y - ey);
            maxY = std::max(maxY, y + ey);
        }
    };

    updateMinMax(gTF_copy);
    updateMinMax(gPS_copy);

    if (minY >= 1e10) {
        minY = 0.0;
        maxY = 0.01;
    }
    
    double margin = 0.2 * (maxY - minY);
    if (margin <= 0) margin = 0.001;
    
    // Create canvas with better formatting
    TCanvas* c = new TCanvas("cCompare", "cCompare", 900, 700);
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetTicks(1, 1);
    
    // Style template fit (circles)
    gTF_copy->SetMarkerStyle(20);
    gTF_copy->SetMarkerColor(kRed+1);
    gTF_copy->SetLineColor(kRed+1);
    gTF_copy->SetMarkerSize(1.3);
    gTF_copy->SetTitle("");
    gTF_copy->GetXaxis()->SetTitle("#eta_{trig}");
    gTF_copy->GetYaxis()->SetTitle("v_{2#Delta}");
    gTF_copy->GetXaxis()->SetLimits(-0.8, 0.8);
    gTF_copy->GetXaxis()->SetRangeUser(-0.8, 0.8);
    gTF_copy->GetYaxis()->SetRangeUser(minY - margin, maxY + margin);
    gTF_copy->GetXaxis()->SetTitleSize(0.05);
    gTF_copy->GetYaxis()->SetTitleSize(0.05);
    gTF_copy->GetXaxis()->SetLabelSize(0.045);
    gTF_copy->GetYaxis()->SetLabelSize(0.045);
    gTF_copy->GetXaxis()->SetTitleOffset(1.1);
    gTF_copy->GetYaxis()->SetTitleOffset(1.3);
    
    // Style peripheral subtraction (open squares)
    gPS_copy->SetMarkerStyle(25);
    gPS_copy->SetMarkerColor(kBlue+1);
    gPS_copy->SetLineColor(kBlue+1);
    gPS_copy->SetMarkerSize(1.3);
    
    // Draw (points only, no lines)
    gTF_copy->Draw("AP");
    gPS_copy->Draw("P SAME");
    
    // Add zero reference line
    TLine* zeroLine = new TLine(-0.8, 0.0, 0.8, 0.0);
    zeroLine->SetLineStyle(kDashed);
    zeroLine->SetLineColor(kGray+1);
    zeroLine->Draw("same");
    
    // Add legend with clearer positioning
    TLegend* leg = new TLegend(0.18, 0.72, 0.55, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.038);
    leg->SetHeader(label.c_str(), "C");
    leg->AddEntry(gTF_copy, "Template Fit (Improved)", "p");
    leg->AddEntry(gPS_copy, "Peripheral Subtraction (ZYAM)", "p");
    leg->Draw();
    
    // Save
    gSystem->mkdir("./TemplateFit/Comparisons/PDFs", kTRUE);
    std::string outName = Form("./TemplateFit/Comparisons/PDFs/MethodComparison_%s_%s.pdf", 
                                dataset.c_str(), splitName.c_str());
    c->SaveAs(outName.c_str());
    
    std::cout << "Comparison plot saved: " << outName << std::endl;
    
    // Cleanup
    delete zeroLine;
    delete leg;
    delete c;
    delete gTF_copy;
    delete gPS_copy;
    fTF->Close();
    fPS->Close();
    delete fTF;
    delete fPS;
}
void CompareMethods_OO_Combined() {
    std::cout << "\n========== Creating combined O-O comparison plot ==========" << std::endl;
    
    const std::string splitName = "Cent";
    
    // Open all four files (TF and PS for both negative and positive eta)
    std::string tfNegFile = Form("./TemplateFit/EtaDiff/Vn_LHC25ae_pass2_604826_%s_0_20.root", splitName.c_str());
    std::string tfPosFile = Form("./TemplateFit/EtaDiff/Vn_LHC25ae_pass2_604830_%s_0_20.root", splitName.c_str());
    std::string psNegFile = "./TemplateFit/PeripheralSubtracted/EtaDiff/Vn_LHC25ae_pass2_604826_PeripheralSubtracted_Cent.root";
    std::string psPosFile = "./TemplateFit/PeripheralSubtracted/EtaDiff/Vn_LHC25ae_pass2_604830_PeripheralSubtracted_Cent.root";
    
    TFile* fTFNeg = TFile::Open(tfNegFile.c_str(), "READ");
    TFile* fTFPos = TFile::Open(tfPosFile.c_str(), "READ");
    TFile* fPSNeg = TFile::Open(psNegFile.c_str(), "READ");
    TFile* fPSPos = TFile::Open(psPosFile.c_str(), "READ");
    
    if (!fTFNeg || !fTFNeg->IsOpen() || !fTFPos || !fTFPos->IsOpen() ||
        !fPSNeg || !fPSNeg->IsOpen() || !fPSPos || !fPSPos->IsOpen()) {
        std::cerr << "Cannot open one or more files for combined O-O comparison" << std::endl;
        if (fTFNeg) { fTFNeg->Close(); delete fTFNeg; }
        if (fTFPos) { fTFPos->Close(); delete fTFPos; }
        if (fPSNeg) { fPSNeg->Close(); delete fPSNeg; }
        if (fPSPos) { fPSPos->Close(); delete fPSPos; }
        return;
    }
    
    // Get all graphs
    TGraphErrors* gTFNeg = (TGraphErrors*)fTFNeg->Get("gV2Delta");
    TGraphErrors* gTFPos = (TGraphErrors*)fTFPos->Get("gV2Delta");
    TGraphErrors* gPSNeg = (TGraphErrors*)fPSNeg->Get("gV2Delta");
    TGraphErrors* gPSPos = (TGraphErrors*)fPSPos->Get("gV2Delta");
    
    if (!gTFNeg || !gTFPos || !gPSNeg || !gPSPos) {
        std::cerr << "Missing gV2Delta in one of the files" << std::endl;
        fTFNeg->Close(); fTFPos->Close(); fPSNeg->Close(); fPSPos->Close();
        delete fTFNeg; delete fTFPos; delete fPSNeg; delete fPSPos;
        return;
    }
    
    // Clone graphs
    TGraphErrors* gTFNeg_copy = (TGraphErrors*)gTFNeg->Clone("gTFNeg_copy");
    TGraphErrors* gTFPos_copy = (TGraphErrors*)gTFPos->Clone("gTFPos_copy");
    TGraphErrors* gPSNeg_copy = (TGraphErrors*)gPSNeg->Clone("gPSNeg_copy");
    TGraphErrors* gPSPos_copy = (TGraphErrors*)gPSPos->Clone("gPSPos_copy");
    
    // Find y-axis range across all graphs
    double minY = 1e10, maxY = -1e10;
    auto updateMinMax = [&](TGraphErrors* g) {
        for (int i = 0; i < g->GetN(); i++) {
            double x, y;
            g->GetPoint(i, x, y);
            double ey = g->GetErrorY(i);
            minY = std::min(minY, y - ey);
            maxY = std::max(maxY, y + ey);
        }
    };
    
    updateMinMax(gTFNeg_copy);
    updateMinMax(gTFPos_copy);
    updateMinMax(gPSNeg_copy);
    updateMinMax(gPSPos_copy);
    
    if (minY >= 1e10) {
        minY = 0.0;
        maxY = 0.01;
    }
    
    double margin = 0.2 * (maxY - minY);
    if (margin <= 0) margin = 0.001;
    
    // Create canvas
    TCanvas* c = new TCanvas("cCompare_OO_Combined", "cCompare_OO_Combined", 900, 700);
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    c->SetTicks(1, 1);
    
    // Style Template Fit - Negative eta (red circles)
    gTFNeg_copy->SetMarkerStyle(20);
    gTFNeg_copy->SetMarkerColor(kRed+1);
    gTFNeg_copy->SetLineColor(kRed+1);
    gTFNeg_copy->SetMarkerSize(1.3);
    gTFNeg_copy->SetTitle("");
    gTFNeg_copy->GetXaxis()->SetTitle("#eta_{trig}");
    gTFNeg_copy->GetYaxis()->SetTitle("v_{2#Delta}");
    gTFNeg_copy->GetXaxis()->SetLimits(-0.8, 0.8);
    gTFNeg_copy->GetXaxis()->SetRangeUser(-0.8, 0.8);
    gTFNeg_copy->GetYaxis()->SetRangeUser(minY - margin, maxY + margin);
    gTFNeg_copy->GetXaxis()->SetTitleSize(0.05);
    gTFNeg_copy->GetYaxis()->SetTitleSize(0.05);
    gTFNeg_copy->GetXaxis()->SetLabelSize(0.045);
    gTFNeg_copy->GetYaxis()->SetLabelSize(0.045);
    gTFNeg_copy->GetXaxis()->SetTitleOffset(1.1);
    gTFNeg_copy->GetYaxis()->SetTitleOffset(1.3);
    
    // Style Template Fit - Positive eta (red open circles)
    gTFPos_copy->SetMarkerStyle(24);
    gTFPos_copy->SetMarkerColor(kRed+1);
    gTFPos_copy->SetLineColor(kRed+1);
    gTFPos_copy->SetMarkerSize(1.3);
    
    // Style Peripheral Subtraction - Negative eta (blue squares)
    gPSNeg_copy->SetMarkerStyle(21);
    gPSNeg_copy->SetMarkerColor(kBlue+1);
    gPSNeg_copy->SetLineColor(kBlue+1);
    gPSNeg_copy->SetMarkerSize(1.3);
    
    // Style Peripheral Subtraction - Positive eta (blue open squares)
    gPSPos_copy->SetMarkerStyle(25);
    gPSPos_copy->SetMarkerColor(kBlue+1);
    gPSPos_copy->SetLineColor(kBlue+1);
    gPSPos_copy->SetMarkerSize(1.3);
    
    // Draw all graphs (points only, no lines)
    gTFNeg_copy->Draw("AP");
    gTFPos_copy->Draw("P SAME");
    gPSNeg_copy->Draw("P SAME");
    gPSPos_copy->Draw("P SAME");
    
    // Add zero reference line
    TLine* zeroLine = new TLine(-0.8, 0.0, 0.8, 0.0);
    zeroLine->SetLineStyle(kDashed);
    zeroLine->SetLineColor(kGray+1);
    zeroLine->Draw("same");
    
    // Add legend
    TLegend* leg = new TLegend(0.18, 0.65, 0.58, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.036);
    leg->SetHeader("O-O collisions, 0-20%", "C");
    leg->AddEntry(gTFNeg_copy, "Template Fit (Improved) - #eta < 0", "p");
    leg->AddEntry(gTFPos_copy, "Template Fit (Improved) - #eta > 0", "p");
    leg->AddEntry(gPSNeg_copy, "Peripheral Sub. (ZYAM) - #eta < 0", "p");
    leg->AddEntry(gPSPos_copy, "Peripheral Sub. (ZYAM) - #eta > 0", "p");
    leg->Draw();
    
    // Save
    gSystem->mkdir("./TemplateFit/Comparisons/PDFs", kTRUE);
    std::string outName = "./TemplateFit/Comparisons/PDFs/MethodComparison_OO_Combined_Cent.pdf";
    c->SaveAs(outName.c_str());
    
    std::cout << "Combined O-O comparison plot saved: " << outName << std::endl;
    
    // Cleanup
    delete zeroLine;
    delete leg;
    delete c;
    delete gTFNeg_copy;
    delete gTFPos_copy;
    delete gPSNeg_copy;
    delete gPSPos_copy;
    fTFNeg->Close(); fTFPos->Close(); fPSNeg->Close(); fPSPos->Close();
    delete fTFNeg; delete fTFPos; delete fPSNeg; delete fPSPos;
}