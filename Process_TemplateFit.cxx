/*
 * @Author: Zhiyong Lu (zhiyong.lu@cern.ch)  
 * @Date: 2025-05-18 13:51:01 
 * @Last Modified by: Zhiyong Lu
 * @Last Modified time: 2025-05-18 19:15:28
 */
//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "TList.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TProfile.h>
#include <TRandom3.h>
#include "TMath.h"
#include "THnSparse.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TSystem.h"
#include "TMultiGraph.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "./include/BasicForDihadron.h"
#include "./include/Bootstrap.h"
#include "./include/TemplateFitter.cxx"
#include "./include/TemplateFunction.C"
#include "./include/ErrorPropagation.h"
#include "./include/plotting.h"

// define struct
struct InputUnit {
    std::string fileNameSuffix;
    Int_t minRange;
    Int_t maxRange;
    InputUnit(std::string _fileNameSuffix, Int_t _minRange, Int_t _maxRange) :
        fileNameSuffix(_fileNameSuffix), minRange(_minRange), maxRange(_maxRange) {}
};

struct ConfigUnit {
    Bool_t isNch;
    Bool_t isEtaDiff;
    InputUnit templ;
    std::vector<InputUnit> dataList;
    std::string outputFileName;
    ConfigUnit(Bool_t _isNch, Bool_t _isEtaDiff, InputUnit _template, std::vector<InputUnit> _dataList, std::string _outputFileName) :
        isNch(_isNch), isEtaDiff(_isEtaDiff), templ(_template), dataList(_dataList), outputFileName(_outputFileName) {}
};

struct VnUnit {
    Double_t v1;
    Double_t v1_err;
    Double_t v2;
    Double_t v2_err;
    Double_t v3;
    Double_t v3_err;
    Double_t v4;
    Double_t v4_err;
    VnUnit(Double_t _v1, Double_t _v1_err, Double_t _v2, Double_t _v2_err, Double_t _v3, Double_t _v3_err, Double_t _v4, Double_t _v4_err) :
        v1(_v1), v1_err(_v1_err), v2(_v2), v2_err(_v2_err), v3(_v3), v3_err(_v3_err), v4(_v4), v4_err(_v4_err) {}
};

// declare functions
void ProcessConfig(Bool_t isNch, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName);
void ProcessConfig_PtDiff(Bool_t isNch, Bool_t isEtaDiff, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName);
VnUnit* TemplateFit(Bool_t isNch, InputUnit templ, InputUnit data, Bool_t cn2Tovn2, Double_t pTMin=0, Double_t pTMax=0, const std::string& systemName = "");
std::vector<Int_t> CheckAndMergeRanges(const std::vector<InputUnit>& inputUnits);
VnUnit* fitSample(Bool_t isNch, TFile* templatefile, InputUnit templ, TFile* datafile, InputUnit data, int sample = -1, Double_t pTMin=0, Double_t pTMax=0, const std::string& systemName = "");
void RooTempFitter(TH1 *lm, TH1 *hm, std::vector<Double_t>& fParamVal, std::vector<Double_t>& fParamErr, Bool_t kRefit);
void PlotFitting(TH1 *lm, TH1 *hm, Bool_t isNch, std::string fileSuffix, Int_t minRange, Int_t maxRange, const std::vector<Double_t>& par, const std::vector<Double_t>& parerr, Double_t pTMin=0, Double_t pTMax=0, const std::string& systemName = "");
void CombineEtaDiffV2Plots_TemplateFit(const std::string& negFile, const std::string& posFile, const std::string& outTag, const std::string& splitName, const std::string& systemName = "");
void CombineEtaDiffV1Plots_TemplateFit(const std::string& negFile, const std::string& posFile, const std::string& outTag, const std::string& splitName, const std::string& systemName = "");

// global variables
std::string collisionSystemName = "peripheral NeNe";
Bool_t gIsEtaDiff = kFALSE;

// Helper function to set collision system name based on filename
std::string GetCollisionSystemName(const std::string& filename) {
    if (filename.find("ae") != std::string::npos) {
        return "O-O";
    } else if (filename.find("af") != std::string::npos) {
        return "Ne-Ne";
    }
    return "Unknown";
}

//==============================================================
void CombineEtaDiffV2Plots_TemplateFit(const std::string& negFile, const std::string& posFile, const std::string& outTag, const std::string& splitName, const std::string& systemName = "") {
    TFile* fNeg = TFile::Open(negFile.c_str(), "READ");
    TFile* fPos = TFile::Open(posFile.c_str(), "READ");

    if (!fNeg || !fNeg->IsOpen() || !fPos || !fPos->IsOpen()) {
        std::cerr << "Cannot open eta-diff files for combined plot: " << negFile << " or " << posFile << std::endl;
        if (fNeg) { fNeg->Close(); delete fNeg; }
        if (fPos) { fPos->Close(); delete fPos; }
        return;
    }

    TGraphErrors* gNeg = (TGraphErrors*)fNeg->Get("gV2Delta");
    TGraphErrors* gPos = (TGraphErrors*)fPos->Get("gV2Delta");
    if (!gNeg || !gPos) {
        std::cerr << "Missing gV2Delta in one of the eta-diff files." << std::endl;
        fNeg->Close();
        fPos->Close();
        delete fNeg;
        delete fPos;
        return;
    }

    struct Point { double x, y, ex, ey; };
    std::vector<Point> points;

    const Int_t nNeg = gNeg->GetN();
    for (Int_t i = 0; i < nNeg; ++i) {
        double x, y;
        gNeg->GetPoint(i, x, y);
        points.push_back({x, y, gNeg->GetErrorX(i), gNeg->GetErrorY(i)});
    }

    const Int_t nPos = gPos->GetN();
    for (Int_t i = 0; i < nPos; ++i) {
        double x, y;
        gPos->GetPoint(i, x, y);
        points.push_back({x, y, gPos->GetErrorX(i), gPos->GetErrorY(i)});
    }

    std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) { return a.x < b.x; });

    TGraphErrors* gCombined = new TGraphErrors(points.size());
    gCombined->SetName("gV2Delta_Combined");
    gCombined->SetTitle("v_{2#Delta} (TemplateFit);#eta_{trig};v_{2#Delta}");

    for (size_t i = 0; i < points.size(); ++i) {
        gCombined->SetPoint(i, points[i].x, points[i].y);
        gCombined->SetPointError(i, points[i].ex, points[i].ey);
    }

    gSystem->mkdir("./TemplateFit/EtaDiff/PDFs", kTRUE);
    gSystem->mkdir("./TemplateFit/EtaDiff", kTRUE);

    TFile outFile(Form("./TemplateFit/EtaDiff/Vn_Combined_%s_%s.root", outTag.c_str(), splitName.c_str()), "RECREATE");
    gCombined->Write();
    outFile.Close();

    // Find y-axis range
    double minY = 1e10, maxY = -1e10;
    for (size_t i = 0; i < points.size(); ++i) {
        minY = std::min(minY, points[i].y - points[i].ey);
        maxY = std::max(maxY, points[i].y + points[i].ey);
    }
    double margin = 0.2 * (maxY - minY);
    if (margin <= 0) margin = 0.001;

    TCanvas* c = new TCanvas("cV2Delta_Combined_TemplateFit", "cV2Delta_Combined_TemplateFit", 900, 700);
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    
    gCombined->SetMarkerStyle(20);
    gCombined->SetMarkerColor(kBlack);
    gCombined->SetLineColor(kBlack);
    gCombined->SetMarkerSize(1.2);
    gCombined->SetTitle("");
    gCombined->Draw("AP");  // "A" = draw axes, "P" = draw points only (no line)
    
    gCombined->GetXaxis()->SetLimits(-0.8, 0.8);
    gCombined->GetXaxis()->SetRangeUser(-0.8, 0.8);
    gCombined->GetYaxis()->SetRangeUser(minY - margin, maxY + margin);
    gCombined->GetXaxis()->SetTitle("#eta_{trig}");
    gCombined->GetYaxis()->SetTitle("v_{2#Delta}");
    gCombined->GetXaxis()->SetTitleSize(0.05);
    gCombined->GetYaxis()->SetTitleSize(0.05);
    gCombined->GetXaxis()->SetLabelSize(0.045);
    gCombined->GetYaxis()->SetLabelSize(0.045);
    gCombined->GetXaxis()->SetTitleOffset(1.1);
    gCombined->GetYaxis()->SetTitleOffset(1.3);
    
    // Add a zero line for reference
    TLine* zeroLine = new TLine(-0.8, 0.0, 0.8, 0.0);
    zeroLine->SetLineStyle(kDashed);
    zeroLine->SetLineColor(kGray+1);
    zeroLine->Draw("same");
    
    c->SaveAs(Form("./TemplateFit/EtaDiff/PDFs/V2Delta_Combined_%s_%s.pdf", outTag.c_str(), splitName.c_str()));

    delete zeroLine;
    delete c;
    delete gCombined;
    fNeg->Close();
    fPos->Close();
    delete fNeg;
    delete fPos;
}

//==============================================================
void CombineEtaDiffV1Plots_TemplateFit(const std::string& negFile, const std::string& posFile, const std::string& outTag, const std::string& splitName, const std::string& systemName = "") {
    TFile* fNeg = TFile::Open(negFile.c_str(), "READ");
    TFile* fPos = TFile::Open(posFile.c_str(), "READ");

    if (!fNeg || !fNeg->IsOpen() || !fPos || !fPos->IsOpen()) {
        std::cerr << "Cannot open eta-diff files for combined V1 plot: " << negFile << " or " << posFile << std::endl;
        if (fNeg) { fNeg->Close(); delete fNeg; }
        if (fPos) { fPos->Close(); delete fPos; }
        return;
    }

    TGraphErrors* gNeg = (TGraphErrors*)fNeg->Get("gV1Delta");
    TGraphErrors* gPos = (TGraphErrors*)fPos->Get("gV1Delta");
    if (!gNeg || !gPos) {
        std::cerr << "Missing gV1Delta in one of the eta-diff files." << std::endl;
        fNeg->Close();
        fPos->Close();
        delete fNeg;
        delete fPos;
        return;
    }

    struct Point { double x, y, ex, ey; };
    std::vector<Point> points;

    const Int_t nNeg = gNeg->GetN();
    for (Int_t i = 0; i < nNeg; ++i) {
        double x, y;
        gNeg->GetPoint(i, x, y);
        points.push_back({x, y, gNeg->GetErrorX(i), gNeg->GetErrorY(i)});
    }

    const Int_t nPos = gPos->GetN();
    for (Int_t i = 0; i < nPos; ++i) {
        double x, y;
        gPos->GetPoint(i, x, y);
        points.push_back({x, y, gPos->GetErrorX(i), gPos->GetErrorY(i)});
    }

    std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) { return a.x < b.x; });

    TGraphErrors* gCombined = new TGraphErrors(points.size());
    gCombined->SetName("gV1Delta_Combined");
    gCombined->SetTitle("v_{1#Delta} (TemplateFit);#eta_{trig};v_{1#Delta}");

    for (size_t i = 0; i < points.size(); ++i) {
        gCombined->SetPoint(i, points[i].x, points[i].y);
        gCombined->SetPointError(i, points[i].ex, points[i].ey);
    }

    gSystem->mkdir("./TemplateFit/EtaDiff/PDFs", kTRUE);
    gSystem->mkdir("./TemplateFit/EtaDiff", kTRUE);

    TFile outFile(Form("./TemplateFit/EtaDiff/Vn_Combined_V1_%s_%s.root", outTag.c_str(), splitName.c_str()), "RECREATE");
    gCombined->Write();
    outFile.Close();

    // Find y-axis range
    double minY = 1e10, maxY = -1e10;
    for (size_t i = 0; i < points.size(); ++i) {
        minY = std::min(minY, points[i].y - points[i].ey);
        maxY = std::max(maxY, points[i].y + points[i].ey);
    }
    double margin = 0.2 * (maxY - minY);
    if (margin <= 0) margin = 0.001;

    TCanvas* c = new TCanvas("cV1Delta_Combined_TemplateFit", "cV1Delta_Combined_TemplateFit", 900, 700);
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);
    
    gCombined->SetMarkerStyle(20);
    gCombined->SetMarkerColor(kRed);
    gCombined->SetLineColor(kRed);
    gCombined->SetMarkerSize(1.2);
    gCombined->SetTitle("");
    gCombined->Draw("AP");
    
    gCombined->GetXaxis()->SetLimits(-0.8, 0.8);
    gCombined->GetXaxis()->SetRangeUser(-0.8, 0.8);
    gCombined->GetYaxis()->SetRangeUser(minY - margin, maxY + margin);
    gCombined->GetXaxis()->SetTitle("#eta_{trig}");
    gCombined->GetYaxis()->SetTitle("v_{1#Delta}");
    gCombined->GetXaxis()->SetTitleSize(0.05);
    gCombined->GetYaxis()->SetTitleSize(0.05);
    gCombined->GetXaxis()->SetLabelSize(0.045);
    gCombined->GetYaxis()->SetLabelSize(0.045);
    gCombined->GetXaxis()->SetTitleOffset(1.1);
    gCombined->GetYaxis()->SetTitleOffset(1.3);
    
    // Add a zero line for reference
    TLine* zeroLine = new TLine(-0.8, 0.0, 0.8, 0.0);
    zeroLine->SetLineStyle(kDashed);
    zeroLine->SetLineColor(kGray+1);
    zeroLine->Draw("same");
    
    c->SaveAs(Form("./TemplateFit/EtaDiff/PDFs/V1Delta_Combined_%s_%s.pdf", outTag.c_str(), splitName.c_str()));

    delete zeroLine;
    delete c;
    delete gCombined;
    fNeg->Close();
    fPos->Close();
    delete fNeg;
    delete fPos;
}

//==============================================================
void Process_TemplateFit() {
    // 不显示窗口
    gROOT->SetBatch(kTRUE);
    std::vector<ConfigUnit> configList;

    // O-O standard (negative eta)
    std::string inputFileName = "LHC25ae_pass2_604826";
    collisionSystemName = GetCollisionSystemName(inputFileName);
    configList.push_back(ConfigUnit(kCent, kEtaDiffOn, InputUnit(inputFileName, 80, 100), 
    {InputUnit(inputFileName, 0, 20)}, 
    inputFileName));

    // O-O reversed (positive eta)
    std::string inputFileNameNew = "LHC25ae_pass2_604830";
    collisionSystemName = GetCollisionSystemName(inputFileNameNew);
    configList.push_back(ConfigUnit(kCent, kEtaDiffOn, InputUnit(inputFileNameNew, 80, 100), 
    {InputUnit(inputFileNameNew, 0, 20)}, 
    inputFileNameNew));

    // Ne-Ne standard (negative eta)
    std::string inputFileNameNeStd = "LHC25af_pass2_611697";
    collisionSystemName = GetCollisionSystemName(inputFileNameNeStd);
    configList.push_back(ConfigUnit(kCent, kEtaDiffOn, InputUnit(inputFileNameNeStd, 80, 100), 
    {InputUnit(inputFileNameNeStd, 0, 20)}, 
    inputFileNameNeStd));

    // Ne-Ne reversed (positive eta)
    std::string inputFileNameNeRev = "LHC25af_pass2_604820";
    collisionSystemName = GetCollisionSystemName(inputFileNameNeRev);
    configList.push_back(ConfigUnit(kCent, kEtaDiffOn, InputUnit(inputFileNameNeRev, 80, 100), 
    {InputUnit(inputFileNameNeRev, 0, 20)}, 
    inputFileNameNeRev));

    for (auto config : configList) {
        if (config.isEtaDiff) {
            ProcessConfig_PtDiff(config.isNch, config.isEtaDiff, config.templ, config.dataList, config.outputFileName);
        } else {
            ProcessConfig(config.isNch, config.templ, config.dataList, config.outputFileName);
        }
    }

    // Combine negative and positive eta sides for template-fit results
    const std::string splitName = "Cent";
    std::string ooNeg = Form("./TemplateFit/EtaDiff/Vn_%s_%s_0_20.root", inputFileName.c_str(), splitName.c_str());
    std::string ooPos = Form("./TemplateFit/EtaDiff/Vn_%s_%s_0_20.root", inputFileNameNew.c_str(), splitName.c_str());

    CombineEtaDiffV2Plots_TemplateFit(ooNeg, ooPos, Form("%s_vs_%s", inputFileName.c_str(), inputFileNameNew.c_str()), splitName, "O-O");
    CombineEtaDiffV1Plots_TemplateFit(ooNeg, ooPos, Form("%s_vs_%s", inputFileName.c_str(), inputFileNameNew.c_str()), splitName, "O-O");
    std::cout << "Created O-O combined plots (V1 and V2) with system name: O-O" << std::endl;

    // Ne-Ne combined
    std::string neNeg = Form("./TemplateFit/EtaDiff/Vn_%s_%s_0_20.root", inputFileNameNeStd.c_str(), splitName.c_str());
    std::string nePos = Form("./TemplateFit/EtaDiff/Vn_%s_%s_0_20.root", inputFileNameNeRev.c_str(), splitName.c_str());
    TFile* checkNeNeg = TFile::Open(neNeg.c_str(), "READ");
    TFile* checkNePos = TFile::Open(nePos.c_str(), "READ");
    
    if (checkNeNeg && checkNeNeg->IsOpen() && checkNePos && checkNePos->IsOpen()) {
        if (checkNeNeg) { checkNeNeg->Close(); delete checkNeNeg; }
        if (checkNePos) { checkNePos->Close(); delete checkNePos; }
        
        CombineEtaDiffV2Plots_TemplateFit(neNeg, nePos, Form("%s_vs_%s", inputFileNameNeStd.c_str(), inputFileNameNeRev.c_str()), splitName, "Ne-Ne");
        CombineEtaDiffV1Plots_TemplateFit(neNeg, nePos, Form("%s_vs_%s", inputFileNameNeStd.c_str(), inputFileNameNeRev.c_str()), splitName, "Ne-Ne");
        
        std::cout << "Created Ne-Ne combined plots (V1 and V2) with system name: Ne-Ne" << std::endl;
    } else {
        std::cout << "Ne-Ne files not available for combining" << std::endl;
        if (checkNeNeg) { checkNeNeg->Close(); delete checkNeNeg; }
        if (checkNePos) { checkNePos->Close(); delete checkNePos; }
    }

    // Ensure ROOT macro exits cleanly in batch mode
    gSystem->Exit(0);

}

//==============================================================
void ProcessConfig(Bool_t isNch, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName) {
    // 按 minRange 排序输入列表
    std::sort(dataList.begin(), dataList.end(), [](const InputUnit& a, const InputUnit& b) {
        return a.minRange < b.minRange;
    });
    
    // print datalist
    std::cout << "Data list: " << std::endl;
    for (const auto& data : dataList) {
        std::cout << "[" << data.minRange << ", " << data.maxRange << "] " << std::endl;
    }

    // 检查范围是否连续
    std::vector<Int_t> mergedRanges = CheckAndMergeRanges(dataList);
    Bool_t isContinuous = !mergedRanges.empty();

    // 执行模板拟合获取所有结果
    std::vector<VnUnit*> vnResults;
    for (const auto& data : dataList) {
        vnResults.push_back(TemplateFit(isNch, templ, data, kTRUE));
    }

    // 创建输出文件
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    TFile outputFile(Form("./TemplateFit/Vn_%s_%s.root", outputFileName.c_str(), splitName.c_str()), "RECREATE");

    if (isContinuous) {
        // 创建可变bin宽度的TH1D
        std::vector<Double_t> binEdges;
        for (auto val : mergedRanges) {
            binEdges.push_back(static_cast<Double_t>(val));
        }

        // 初始化直方图
        TH1D* hV2 = new TH1D("hV2", "v_{2};Centrality;v_{2}", 
                            mergedRanges.size()-1, binEdges.data());
        TH1D* hV3 = new TH1D("hV3", "v_{3};Centrality;v_{3}", 
                            mergedRanges.size()-1, binEdges.data());
        TH1D* hV4 = new TH1D("hV4", "v_{4};Centrality;v_{4}", 
                            mergedRanges.size()-1, binEdges.data());

        // 填充数据
        for (size_t i = 0; i < vnResults.size(); ++i) {
            hV2->SetBinContent(i+1, vnResults[i]->v2);
            hV2->SetBinError(i+1, vnResults[i]->v2_err);
            
            hV3->SetBinContent(i+1, vnResults[i]->v3);
            hV3->SetBinError(i+1, vnResults[i]->v3_err);
            
            hV4->SetBinContent(i+1, vnResults[i]->v4);
            hV4->SetBinError(i+1, vnResults[i]->v4_err);
        }

        // 写入文件
        hV2->Write();
        hV3->Write();
        hV4->Write();
    } else {
        // 创建TGraphErrors
        Int_t nPoints = dataList.size();
        TGraphErrors* gV2 = new TGraphErrors(nPoints);
        TGraphErrors* gV3 = new TGraphErrors(nPoints);
        TGraphErrors* gV4 = new TGraphErrors(nPoints);

        // 设置标题
        std::string Xtitle = "Centrality (%)";
        if (isNch) Xtitle = "Nch";
        gV2->SetNameTitle("gV2", Form("v_{2};%s;v_{2}", Xtitle.c_str()));
        gV3->SetNameTitle("gV3", Form("v_{3};%s;v_{3}", Xtitle.c_str()));
        gV4->SetNameTitle("gV4", Form("v_{4};%s;v_{4}", Xtitle.c_str()));

        // 填充数据
        for (Int_t i = 0; i < nPoints; ++i) {
            Double_t xCenter = 0.5*(dataList[i].minRange + dataList[i].maxRange);
            Double_t xError = 0.5*(dataList[i].maxRange - dataList[i].minRange);
            
            gV2->SetPoint(i, xCenter, vnResults[i]->v2);
            gV2->SetPointError(i, xError, vnResults[i]->v2_err);
            
            gV3->SetPoint(i, xCenter, vnResults[i]->v3);
            gV3->SetPointError(i, xError, vnResults[i]->v3_err);
            
            gV4->SetPoint(i, xCenter, vnResults[i]->v4);
            gV4->SetPointError(i, xError, vnResults[i]->v4_err);
        }

        // 写入文件
        gV2->Write();
        gV3->Write();
        gV4->Write();
    }

    std::cout << "Output file: " << Form("./TemplateFit/Vn_%s_%s.root", outputFileName.c_str(), splitName.c_str()) << std::endl;
    outputFile.Close();
}

//==============================================================
void ProcessConfig_PtDiff(Bool_t isNch, Bool_t isEtaDiff, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName) {
    gIsEtaDiff = isEtaDiff;
    
    // Determine collision system name
    std::string systemName = GetCollisionSystemName(outputFileName);
    
    // Determine which eta bins to use based on dataset
    std::vector<float> etaBinsNeg = {-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0};
    std::vector<float> etaBinsPos = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    
    // Check if this is a reversed dataset
    bool isReversed = (outputFileName.find("604830") != std::string::npos || 
                      outputFileName.find("604820") != std::string::npos);
    std::vector<float>& etaBinsToUse = isReversed ? etaBinsPos : etaBinsNeg;
    
    // Just looping the list - store graphs for combined plots
    std::vector<TGraphErrors*> v1Graphs;
    std::vector<std::string> v1Labels;
    std::vector<TGraphErrors*> v2Graphs;
    std::vector<std::string> v2Labels;
    for (const auto& data : dataList) {
        // 执行模板拟合获取所有结果
        std::vector<VnUnit*> vnResults;
        std::vector<double> xCenters;
        std::vector<double> xErrors;
        for (Int_t iPt = 0; iPt < etaBinsToUse.size() - 1; iPt++) {
            double pTMin = etaBinsToUse[iPt];
            double pTMax = etaBinsToUse[iPt + 1];
            VnUnit* res = TemplateFit(isNch, templ, data, kFALSE, pTMin, pTMax, systemName);
            if (!res) {
                std::cerr << "Fit failed for bin " << pTMin << " to " << pTMax << ", keeping placeholder point with large errors." << std::endl;
                res = new VnUnit(0.0, 10.0, 0.0, 10.0, 0.0, 10.0, 0.0, 10.0);
            }
            if (!std::isfinite(res->v2) || !std::isfinite(res->v2_err)) {
                std::cerr << "Skip bin " << pTMin << " to " << pTMax << " (invalid error/fit: v2="
                          << res->v2 << ", err=" << res->v2_err << ")" << std::endl;
                continue;
            }
            vnResults.push_back(res);
            xCenters.push_back(0.5 * (pTMin + pTMax));
            xErrors.push_back(0.5 * (pTMax - pTMin));
        }

        // 创建输出文件
        std::string splitName = "Mult";
        if (!isNch) splitName = "Cent";
        if (gIsEtaDiff) {
            gSystem->mkdir("./TemplateFit/EtaDiff", kTRUE);
        } else {
            gSystem->mkdir("./TemplateFit/PtDiff", kTRUE);
        }
        TFile outputFile(Form("./TemplateFit/%s/Vn_%s_%s_%i_%i.root",
                              gIsEtaDiff ? "EtaDiff" : "PtDiff",
                              outputFileName.c_str(), splitName.c_str(), data.minRange, data.maxRange), "RECREATE");

        // 初始化直方图
        const char* xTitle = gIsEtaDiff ? "#eta_{trig}" : "p_{T}";
        TH1D* hV1 = nullptr;
        if (gIsEtaDiff) {
            hV1 = new TH1D("hV1", Form("v_{1};%s;v_{1}", xTitle),
                        etaBinsToUse.size()-1, etaBinsToUse.data());
        }
        TH1D* hV2 = new TH1D("hV2", Form("v_{2};%s;v_{2}", xTitle), 
                    etaBinsToUse.size()-1, etaBinsToUse.data());
        TH1D* hV3 = new TH1D("hV3", Form("v_{3};%s;v_{3}", xTitle), 
                    etaBinsToUse.size()-1, etaBinsToUse.data());
        TH1D* hV4 = new TH1D("hV4", Form("v_{4};%s;v_{4}", xTitle), 
                    etaBinsToUse.size()-1, etaBinsToUse.data());

        // 填充数据 (only valid bins)
        for (size_t i = 0; i < vnResults.size(); ++i) {
            int bin = hV2->FindBin(xCenters[i]);
            if (gIsEtaDiff) {
                hV1->SetBinContent(bin, vnResults[i]->v1);
                hV1->SetBinError(bin, vnResults[i]->v1_err);
            }
            hV2->SetBinContent(bin, vnResults[i]->v2);
            hV2->SetBinError(bin, vnResults[i]->v2_err);

            hV3->SetBinContent(bin, vnResults[i]->v3);
            hV3->SetBinError(bin, vnResults[i]->v3_err);

            hV4->SetBinContent(bin, vnResults[i]->v4);
            hV4->SetBinError(bin, vnResults[i]->v4_err);
        }

        // Create V1Δ graph (eta-diff only)
        TGraphErrors* gV1 = nullptr;
        if (gIsEtaDiff) {
            Int_t nBins = static_cast<Int_t>(vnResults.size());
            gV1 = new TGraphErrors(nBins);
            gV1->SetName("gV1Delta");
            gV1->SetTitle(Form("v_{1#Delta};%s;v_{1#Delta}", xTitle));
            for (Int_t i = 0; i < nBins; ++i) {
                gV1->SetPoint(i, xCenters[i], vnResults[i]->v1);
                gV1->SetPointError(i, xErrors[i], vnResults[i]->v1_err);
            }
        }

        // Create V2Δ graph
        Int_t nBins = static_cast<Int_t>(vnResults.size());
        TGraphErrors* gV2 = new TGraphErrors(nBins);
        gV2->SetName("gV2Delta");
        gV2->SetTitle(Form("v_{2#Delta};%s;v_{2#Delta}", xTitle));
        for (Int_t i = 0; i < nBins; ++i) {
            gV2->SetPoint(i, xCenters[i], vnResults[i]->v2);
            gV2->SetPointError(i, xErrors[i], vnResults[i]->v2_err);
        }

        // Write to file
        if (gIsEtaDiff && hV1) {
            hV1->Write();
        }
        hV2->Write();
        hV3->Write();
        hV4->Write();
        if (gIsEtaDiff && gV1) {
            gV1->Write();
        }
        gV2->Write();

        std::cout << "Output file: " << Form("./TemplateFit/%s/Vn_%s_%s_%i_%i.root",
                                              gIsEtaDiff ? "EtaDiff" : "PtDiff",
                                              outputFileName.c_str(), splitName.c_str(), data.minRange, data.maxRange) << std::endl;
        
        // Save per-correlation V1Δ plot (eta-diff only)
        if (gIsEtaDiff && gV1) {
            gSystem->mkdir("./TemplateFit/EtaDiff/PDFs", kTRUE);
            TCanvas* cV1 = new TCanvas("cV1Delta", "cV1Delta", 800, 600);
            gV1->SetMarkerStyle(20);
            gV1->SetMarkerColor(kRed);
            gV1->SetLineColor(kRed);
            gV1->Draw("AP");
            cV1->SaveAs(Form("./TemplateFit/EtaDiff/PDFs/V1Delta_%s_%s_%d_%d.pdf",
                             outputFileName.c_str(), splitName.c_str(), data.minRange, data.maxRange));
            delete cV1;
        }
        
        // Save per-correlation V2Δ plot
        if (gIsEtaDiff) {
            gSystem->mkdir("./TemplateFit/EtaDiff/PDFs", kTRUE);
        } else {
            gSystem->mkdir("./TemplateFit/PtDiff/PDFs", kTRUE);
        }
        TCanvas* cV2 = new TCanvas("cV2Delta", "cV2Delta", 800, 600);
        gV2->SetMarkerStyle(20);
        gV2->SetMarkerColor(kBlack);
        gV2->SetLineColor(kBlack);
        gV2->Draw("AP");
        cV2->SaveAs(Form("./TemplateFit/%s/PDFs/V2Delta_%s_%s_%d_%d.pdf",
                         gIsEtaDiff ? "EtaDiff" : "PtDiff",
                         outputFileName.c_str(), splitName.c_str(), data.minRange, data.maxRange));
        delete cV2;
        outputFile.Close();

        if (nBins > 0) {
            v2Graphs.push_back(gV2);
            v2Labels.push_back(Form("%d-%d", data.minRange, data.maxRange));
            if (gIsEtaDiff && gV1) {
                v1Graphs.push_back(gV1);
                v1Labels.push_back(Form("%d-%d", data.minRange, data.maxRange));
            }
        }

    }

    // Combined overlay of V1Δ for all correlations (eta-diff only)
    if (gIsEtaDiff && v1Graphs.size() > 1) {
        std::string splitName = "Mult";
        if (!isNch) splitName = "Cent";
        gSystem->mkdir("./TemplateFit/EtaDiff/PDFs", kTRUE);
        TCanvas* cAllV1 = new TCanvas("cV1DeltaAll", "cV1DeltaAll", 900, 700);
        TMultiGraph* mgV1 = new TMultiGraph();
        TLegend* legV1 = new TLegend(0.65, 0.70, 0.88, 0.88);
        legV1->SetBorderSize(0);
        for (size_t i = 0; i < v1Graphs.size(); ++i) {
            v1Graphs[i]->SetMarkerStyle(20 + i);
            v1Graphs[i]->SetLineColor(colors[i % 8]);
            v1Graphs[i]->SetMarkerColor(colors[i % 8]);
            mgV1->Add(v1Graphs[i], "LP");
            legV1->AddEntry(v1Graphs[i], v1Labels[i].c_str(), "lp");
        }
        mgV1->SetTitle("v_{1#Delta} (all);#eta_{trig};v_{1#Delta}");
        mgV1->Draw("A");
        legV1->Draw();
        cAllV1->SaveAs(Form("./TemplateFit/EtaDiff/PDFs/V1Delta_All_%s_%s.pdf",
                          outputFileName.c_str(), splitName.c_str()));
        delete legV1;
        delete mgV1;
        delete cAllV1;
    }

    // Combined overlay of V2Δ for all correlations
    if (v2Graphs.size() > 1) {
        std::string splitName = "Mult";
        if (!isNch) splitName = "Cent";
        if (gIsEtaDiff) {
            gSystem->mkdir("./TemplateFit/EtaDiff/PDFs", kTRUE);
        } else {
            gSystem->mkdir("./TemplateFit/PtDiff/PDFs", kTRUE);
        }
        TCanvas* cAll = new TCanvas("cV2DeltaAll", "cV2DeltaAll", 900, 700);
        TMultiGraph* mg = new TMultiGraph();
        TLegend* leg = new TLegend(0.65, 0.70, 0.88, 0.88);
        leg->SetBorderSize(0);
        for (size_t i = 0; i < v2Graphs.size(); ++i) {
            v2Graphs[i]->SetMarkerStyle(20 + i);
            v2Graphs[i]->SetLineColor(colors[i % 8]);
            v2Graphs[i]->SetMarkerColor(colors[i % 8]);
            mg->Add(v2Graphs[i], "LP");
            leg->AddEntry(v2Graphs[i], v2Labels[i].c_str(), "lp");
        }
        const char* xTitle = gIsEtaDiff ? "#eta_{trig}" : "p_{T}";
        mg->SetTitle(Form("v_{2#Delta} (all);%s;v_{2#Delta}", xTitle));
        mg->Draw("A");
        leg->Draw();
        cAll->SaveAs(Form("./TemplateFit/%s/PDFs/V2Delta_All_%s_%s.pdf",
                          gIsEtaDiff ? "EtaDiff" : "PtDiff",
                          outputFileName.c_str(), splitName.c_str()));
        delete leg;
        delete mg;
        delete cAll;
    }

}

void VnPtDiff(VnUnit* VnResult_PtDiff, VnUnit* VnResult_Ref) {
    if (!VnResult_PtDiff || !VnResult_Ref) {
        std::cerr << "Invalid VnUnit pointers provided." << std::endl;
        return;
    }
    // assuming V_{n\Delta}
    if (VnResult_Ref->v2 > 0) {
        double v2_2PC = VnResult_PtDiff->v2 / TMath::Sqrt(VnResult_Ref->v2);
        // synthetic error
        double v2_2PC_err = Error_Ratio_sqrtY(VnResult_PtDiff->v2, VnResult_PtDiff->v2_err, VnResult_Ref->v2, VnResult_Ref->v2_err);
        VnResult_PtDiff->v2 = v2_2PC;
        VnResult_PtDiff->v2_err = v2_2PC_err;
    } else {
        VnResult_PtDiff->v2 = -1.;
        VnResult_PtDiff->v2_err = 10.;
    }
    if (VnResult_Ref->v3 > 0) {
        double v3_2PC = VnResult_PtDiff->v3 / TMath::Sqrt(VnResult_Ref->v3);
        // synthetic error
        double v3_2PC_err = Error_Ratio_sqrtY(VnResult_PtDiff->v3, VnResult_PtDiff->v3_err, VnResult_Ref->v3, VnResult_Ref->v3_err);
        VnResult_PtDiff->v3 = v3_2PC;
        VnResult_PtDiff->v3_err = v3_2PC_err;
    } else {
        VnResult_PtDiff->v3 = -1.;
        VnResult_PtDiff->v3_err = 10.;
    }
    if (VnResult_Ref->v4 > 0) {
        double v4_2PC = VnResult_PtDiff->v4 / TMath::Sqrt(VnResult_Ref->v4);
        // synthetic error
        double v4_2PC_err = Error_Ratio_sqrtY(VnResult_PtDiff->v4, VnResult_PtDiff->v4_err, VnResult_Ref->v4, VnResult_Ref->v4_err);
        VnResult_PtDiff->v4 = v4_2PC;
        VnResult_PtDiff->v4_err = v4_2PC_err;
    } else {
        VnResult_PtDiff->v4 = -1.;
        VnResult_PtDiff->v4_err = 10.;
    }
    return;
}

//==============================================================
VnUnit* TemplateFit(Bool_t isNch, InputUnit templ, InputUnit data, Bool_t cn2Tovn2, Double_t pTMin=0, Double_t pTMax=0, const std::string& systemName = "") {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    const Bool_t useDiff = gIsEtaDiff || (pTMin > 0 && pTMax > 0);

    TFile* templatefile = new TFile(Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d.root", templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange), "READ");
    if (!templatefile || !templatefile->IsOpen()) {
        std::cerr << "Cannot open template file: " << Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d.root", templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange) << std::endl;
        return nullptr;
    }
    // For eta/pT-differential analysis, template is the same (not binned in eta/pT)
    // So templatefile_PtDiff should point to the same integrated template file
    TFile* templatefile_PtDiff = templatefile;

    TFile* datafile = new TFile(Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange), "READ");
    if (!datafile || !datafile->IsOpen()) {
        std::cerr << "Cannot open input file: " << Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange) << std::endl;
        return nullptr;
    }
    TFile* datafile_PtDiff = nullptr;
    if (useDiff) {
        // First try EtaDiff folder (for eta-differential analysis)
        datafile_PtDiff = new TFile(Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%i_%i_Eta_%0.1f_%0.1f.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, pTMin, pTMax), "READ");
        if (!datafile_PtDiff || !datafile_PtDiff->IsOpen()) {
            // If not found, try PtDiff folder (for pT-differential analysis)
            datafile_PtDiff = new TFile(Form("./ProcessOutput/PtDiff/BootstrapSample_%s_%s_%i_%i_Pt_%0.1f_%0.1f.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, pTMin, pTMax), "READ");
            if (!datafile_PtDiff || !datafile_PtDiff->IsOpen()) {
                std::cerr << "Cannot open data diff file (EtaDiff/PtDiff) for bin " << pTMin << " to " << pTMax << std::endl;
                return nullptr;
            }
        }
    }

    VnUnit* vnResult = fitSample(isNch, templatefile, templ, datafile, data, -1, 0, 0, systemName);
    if (!vnResult) {
        std::cerr << "Cannot fit sample: " << data.fileNameSuffix << std::endl;
        return nullptr;
    }
    VnUnit* vnResult_PtDiff = nullptr;
    if (useDiff) {
        vnResult_PtDiff = fitSample(isNch, templatefile_PtDiff, templ, datafile_PtDiff, data, -1, pTMin, pTMax, systemName);
        if (!vnResult_PtDiff) {
            std::cerr << "Cannot fit diff sample: " << data.fileNameSuffix << std::endl;
            return nullptr;
        }
    }

    std::vector<std::vector<std::vector<double>>> ValueArray;
    std::vector<std::vector<std::vector<double>>> ValueErrorArray;
    std::vector<std::vector<double>> ErrorArray;
    int Nobs=4;//v1,v2,v3,v4
    int NofSample = maxSample*maxSample;
    int Nbin = 1;
    ResizeValueArray(ValueArray,ValueErrorArray,ErrorArray,Nobs,NofSample,Nbin);

    for(int sample=0;sample<NofSample;sample++) {
        VnUnit* vnTemp = fitSample(isNch, templatefile, templ, datafile, data, sample, 0, 0, systemName);
        if (!vnTemp) {
            // Skip this failed bootstrap sample
            ValueArray[0][sample][0] = -999;  // Mark as invalid
            ValueErrorArray[0][sample][0] = 0;
            ValueArray[1][sample][0] = -999;
            ValueErrorArray[1][sample][0] = 0;
            ValueArray[2][sample][0] = -999;
            ValueErrorArray[2][sample][0] = 0;
            ValueArray[3][sample][0] = -999;
            ValueErrorArray[3][sample][0] = 0;
            continue;
        }
        ValueArray[0][sample][0] = vnTemp->v1;
        ValueErrorArray[0][sample][0] = vnTemp->v1_err;
        ValueArray[1][sample][0] = vnTemp->v2;
        ValueErrorArray[1][sample][0] = vnTemp->v2_err;
        ValueArray[2][sample][0] = vnTemp->v3;
        ValueErrorArray[2][sample][0] = vnTemp->v3_err;
        ValueArray[3][sample][0] = vnTemp->v4;
        ValueErrorArray[3][sample][0] = vnTemp->v4_err;
        if (useDiff) {
            VnUnit* vnTemp_PtDiff = fitSample(isNch, templatefile_PtDiff, templ, datafile_PtDiff, data, sample, pTMin, pTMax, systemName);
            if (!vnTemp_PtDiff) {
                // Skip this failed bootstrap sample
                ValueArray[0][sample][0] = -999;
                ValueErrorArray[0][sample][0] = 0;
                ValueArray[1][sample][0] = -999;
                ValueErrorArray[1][sample][0] = 0;
                ValueArray[2][sample][0] = -999;
                ValueErrorArray[2][sample][0] = 0;
                ValueArray[3][sample][0] = -999;
                ValueErrorArray[3][sample][0] = 0;
                delete vnTemp;
                continue;
            }
            VnPtDiff(vnTemp_PtDiff, vnTemp);
            ValueArray[0][sample][0] = vnTemp_PtDiff->v1;
            ValueErrorArray[0][sample][0] = vnTemp_PtDiff->v1_err;
            ValueArray[1][sample][0] = vnTemp_PtDiff->v2;
            ValueErrorArray[1][sample][0] = vnTemp_PtDiff->v2_err;
            ValueArray[2][sample][0] = vnTemp_PtDiff->v3;
            ValueErrorArray[2][sample][0] = vnTemp_PtDiff->v3_err;
            ValueArray[3][sample][0] = vnTemp_PtDiff->v4;
            ValueErrorArray[3][sample][0] = vnTemp_PtDiff->v4_err;
            delete vnTemp_PtDiff;
        }
        delete vnTemp;
    }
    // for(int sample=0;sample<NofSample;sample++) {
    //     std::cout << "sample: " << sample << " v2^2: " << ValueArray[0][sample][0]  << std::endl;
    //     std::cout << "sample: " << sample << " v3^2: " << ValueArray[1][sample][0]  << std::endl;
    //     std::cout << "sample: " << sample << " v4^2: " << ValueArray[2][sample][0] << std::endl << std::endl;
    // }
    
    // Count successful bootstrap samples
    int successfulSamples = 0;
    for(int sample=0;sample<NofSample;sample++) {
        if (ValueArray[0][sample][0] > -900) {
            successfulSamples++;
        }
    }
    
    // Require a minimum number of successful bootstrap samples (looser for eta/pT-diff)
    int minSamples = NofSample / 50; // 2%
    if (minSamples < 5) minSamples = 5;
    if (successfulSamples < minSamples) {
        std::cerr << "Too few successful bootstrap fits (" << successfulSamples << "/" << NofSample << ") for " 
                  << data.fileNameSuffix << " bin " << pTMin << " to " << pTMax << std::endl;
        return nullptr;
    }
    
    std::cout << "Bootstrap fits: " << successfulSamples << "/" << NofSample << " successful" << std::endl;
    
    for(int iobs = 0;iobs < Nobs;iobs++){
        CalculateBootstrapError(ValueArray[iobs],ValueErrorArray[iobs],ErrorArray[iobs],1.);
    }

    vnResult->v1_err = ErrorArray[0][0];
    vnResult->v2_err = ErrorArray[1][0];
    vnResult->v3_err = ErrorArray[2][0];
    vnResult->v4_err = ErrorArray[3][0];
    if (useDiff) {
        vnResult_PtDiff->v1_err = ErrorArray[0][0];
        vnResult_PtDiff->v2_err = ErrorArray[1][0];
        vnResult_PtDiff->v3_err = ErrorArray[2][0];
        vnResult_PtDiff->v4_err = ErrorArray[3][0];
    }

    if (cn2Tovn2) {
        if (vnResult->v2 > 0.) {
            vnResult->v2_err = vnResult->v2_err / (2 * sqrt(vnResult->v2));
            vnResult->v2 = sqrt(vnResult->v2);
        }
        else {
            vnResult->v2 = -1;
            vnResult->v2_err = 10.;
        }
        
        if (vnResult->v3 > 0.) {
            vnResult->v3_err = vnResult->v3_err / (2 * sqrt(vnResult->v3));
            vnResult->v3 = sqrt(vnResult->v3);
        }
        else {
            vnResult->v3 = -1;
            vnResult->v3_err = 10.;
        }
        
        if (vnResult->v4 > 0.) {
            vnResult->v4_err = vnResult->v4_err / (2 * sqrt(vnResult->v4));
            vnResult->v4 = sqrt(vnResult->v4);
        }
        else {
            vnResult->v4 = -1;
            vnResult->v4_err = 10.;
        }
        
        // Also apply cn2 to vn2 conversion for differential result
        if (useDiff) {
            if (vnResult_PtDiff->v2 > 0.) {
                vnResult_PtDiff->v2_err = vnResult_PtDiff->v2_err / (2 * sqrt(vnResult_PtDiff->v2));
                vnResult_PtDiff->v2 = sqrt(vnResult_PtDiff->v2);
            }
            else {
                vnResult_PtDiff->v2 = -1;
                vnResult_PtDiff->v2_err = 10.;
            }
            
            if (vnResult_PtDiff->v3 > 0.) {
                vnResult_PtDiff->v3_err = vnResult_PtDiff->v3_err / (2 * sqrt(vnResult_PtDiff->v3));
                vnResult_PtDiff->v3 = sqrt(vnResult_PtDiff->v3);
            }
            else {
                vnResult_PtDiff->v3 = -1;
                vnResult_PtDiff->v3_err = 10.;
            }
            
            if (vnResult_PtDiff->v4 > 0.) {
                vnResult_PtDiff->v4_err = vnResult_PtDiff->v4_err / (2 * sqrt(vnResult_PtDiff->v4));
                vnResult_PtDiff->v4 = sqrt(vnResult_PtDiff->v4);
            }
            else {
                vnResult_PtDiff->v4 = -1;
                vnResult_PtDiff->v4_err = 10.;
            }
            
            // Now that both have been sqrt'd, compute V2Delta
            VnPtDiff(vnResult_PtDiff, vnResult);
        }
    }

    // print result
    if (useDiff) {
        std::cout << "print result: " << data.fileNameSuffix << " pT-diff" << std::endl;
        std::cout << "v1: " << vnResult_PtDiff->v1 << " +/- " << vnResult_PtDiff->v1_err << std::endl;
        std::cout << "v2: " << vnResult_PtDiff->v2 << " +/- " << vnResult_PtDiff->v2_err << std::endl;
        std::cout << "v3: " << vnResult_PtDiff->v3 << " +/- " << vnResult_PtDiff->v3_err << std::endl;
        std::cout << "v4: " << vnResult_PtDiff->v4 << " +/- " << vnResult_PtDiff->v4_err << std::endl;
        return vnResult_PtDiff;
    }
    std::cout << "print result: " << data.fileNameSuffix << std::endl;
    std::cout << "v1: " << vnResult->v1 << " +/- " << vnResult->v1_err << std::endl;
    std::cout << "v2: " << vnResult->v2 << " +/- " << vnResult->v2_err << std::endl;
    std::cout << "v3: " << vnResult->v3 << " +/- " << vnResult->v3_err << std::endl;
    std::cout << "v4: " << vnResult->v4 << " +/- " << vnResult->v4_err << std::endl;

    return vnResult;
}


//==============================================================
std::vector<Int_t> CheckAndMergeRanges(const std::vector<InputUnit>& inputUnits) {
    std::vector<InputUnit> sortedUnits = inputUnits;
    // 按 minRange 升序排序
    std::sort(sortedUnits.begin(), sortedUnits.end(),
              [](const InputUnit& a, const InputUnit& b) {
                  return a.minRange < b.minRange;
              });

    // 检查连续性
    bool isContinuous = true;
    for (size_t i = 0; i < sortedUnits.size() - 1; ++i) {
        if (sortedUnits[i].maxRange != sortedUnits[i + 1].minRange) {
            isContinuous = false;
            break;
        }
    }

    // 如果连续则生成合并后的范围
    std::vector<Int_t> mergedRanges;
    if (isContinuous && !sortedUnits.empty()) {
        mergedRanges.reserve(sortedUnits.size() + 1);
        mergedRanges.push_back(sortedUnits[0].minRange);
        for (const auto& unit : sortedUnits) {
            mergedRanges.push_back(unit.maxRange);
        }
    }
    return mergedRanges;
}

//==============================================================
VnUnit* fitSample(Bool_t isNch, TFile* templatefile, InputUnit templ, TFile* datafile, InputUnit data, int sample = -1, Double_t pTMin=0, Double_t pTMax=0, const std::string& systemName = "") {
    std::vector<Double_t> fParamVal;
    std::vector<Double_t> fParamErr;
    TH1D* lm=0;
    TH1D* hm=0;
    TString suffix = (sample == -1) ? "" : Form("_%d", sample);
    lm = (TH1D*)templatefile->Get(Form("bsSample_hPhiSameOverMixed_%d_%d%s", templ.minRange, templ.maxRange, suffix.Data()));
    hm = (TH1D*)datafile->Get(Form("bsSample_hPhiSameOverMixed_%d_%d%s", data.minRange, data.maxRange, suffix.Data()));
    if (!lm) {
        std::cerr << "Cannot find histogram: " << Form("bsSample_hPhiSameOverMixed_%d_%d%s", templ.minRange, templ.maxRange, suffix.Data()) << std::endl;
        return 0;
    }
    if (!hm) {
        std::cerr << "Cannot find histogram: " << Form("bsSample_hPhiSameOverMixed_%d_%d%s", data.minRange, data.maxRange, suffix.Data()) << std::endl;
        return 0;
    }
    
    // Check template histogram quality
    if (lm->GetEntries() == 0 || lm->Integral() <= 0) {
        std::cerr << "Empty template histogram for " << Form("templ %d-%d%s", templ.minRange, templ.maxRange, suffix.Data()) << std::endl;
        return 0;
    }
    
    if (lm->Integral() <= 0 || hm->Integral() <= 0) {
        std::cerr << "Empty histogram for fit: " << Form("templ %d-%d, data %d-%d%s", templ.minRange, templ.maxRange, data.minRange, data.maxRange, suffix.Data()) << std::endl;
        return 0;
    }

    // Ensure non-zero bin errors to avoid RooChi2Var infinities
    for (int i = 1; i <= hm->GetNbinsX(); ++i) {
        double err = hm->GetBinError(i);
        double val = hm->GetBinContent(i);
        if (!(err > 0.0) || std::isnan(err) || std::isinf(err)) {
            double fallback = std::sqrt(std::abs(val));
            if (!(fallback > 0.0)) fallback = 1e-6;
            hm->SetBinError(i, fallback);
        }
        // Also check template histogram
        double templ_err = lm->GetBinError(i);
        double templ_val = lm->GetBinContent(i);
        if (!(templ_err > 0.0) || std::isnan(templ_err) || std::isinf(templ_err)) {
            double templ_fallback = std::sqrt(std::abs(templ_val));
            if (!(templ_fallback > 0.0)) templ_fallback = 1e-6;
            lm->SetBinError(i, templ_fallback);
        }
    }
    
    // Check for sufficient statistics (use a lower threshold for eta/pT-diff)
    const bool isDiff = gIsEtaDiff || (pTMin != 0.0 || pTMax != 0.0);
    const double minIntegral = isDiff ? 0.1 : 10.0;
    if (hm->Integral() < minIntegral) {
        std::cerr << "Insufficient statistics for fit (integral=" << hm->Integral() << "): " 
                  << Form("data %d-%d%s", data.minRange, data.maxRange, suffix.Data()) << std::endl;
        return 0;
    }
    
    RooTempFitter(lm, hm, fParamVal, fParamErr, kFALSE);
    for (size_t i = 0; i < fParamVal.size(); ++i) {
        if (std::isnan(fParamVal[i]) || std::isinf(fParamVal[i]) ||
            std::isnan(fParamErr[i]) || std::isinf(fParamErr[i])) {
            std::cerr << "Invalid fit parameters for " << Form("%d-%d%s", data.minRange, data.maxRange, suffix.Data()) << std::endl;
            return 0;
        }
    }
    // Guard against silent fit failures (all-zero parameters/errors or unreasonable values)
    bool allZero = true;
    for (size_t i = 0; i < fParamVal.size(); ++i) {
        if (fParamVal[i] != 0.0 || fParamErr[i] != 0.0) {
            allZero = false;
            break;
        }
    }
    if (allZero) {
        std::cerr << "Fit returned all-zero parameters for " << Form("%d-%d%s", data.minRange, data.maxRange, suffix.Data()) << std::endl;
        return 0;
    }
    
    // Check if fit failed (RooFit sets default values F=0, G=0, V2=-1 with err=10 for failed fits)
    // Check both F and V2 to detect Hessian errors
    if ((fParamVal[3] == 0.0 && fParamErr[3] == 10.0) || 
        (fParamVal[0] == -1.0 && fParamErr[0] == 10.0)) {
        // Fit failed - skip silently for bootstrap samples
        return 0;
    }
    
    if (sample == -1) {
        PlotFitting(lm, hm, isNch, data.fileNameSuffix, data.minRange, data.maxRange, fParamVal, fParamErr, pTMin, pTMax, systemName);
    }
    VnUnit* vnResult = new VnUnit(fParamVal[5], fParamErr[5], fParamVal[0], fParamErr[0], fParamVal[1], fParamErr[1], fParamVal[2], fParamErr[2]);
    return vnResult;
}

//==============================================================
void RooTempFitter(TH1 *lm, TH1 *hm, std::vector<Double_t>& fParamVal, std::vector<Double_t>& fParamErr, Bool_t kRefit) {
    fParamVal.clear();
    fParamErr.clear();
    fParamVal.resize(6);
    fParamErr.resize(6);
    if (!lm ||!hm) {
        std::cerr << "Null pointer to histogram" << std::endl;
        return;
    }
    //Initialize fitter with given projections
    TemplateFitter *ft = new TemplateFitter(hm);
    //Setting up variable ( = delta phi, or just "x"):
    ft->AddVariable("x", "x", -TMath::Pi()/2.0, 1.5*TMath::Pi());

    if (!kRefit){
        // Pb-Pb initial value
        ft->AddParameter("Fa","Fa",4.5,0,30);
        ft->AddParameter("Ga","Ga",15,0,1000);
        ft->AddParameter("v1","v1",0.0,-1.0,1.0);
        ft->AddParameter("v2","v2",4e-3,-1.0,1.0);
        ft->AddParameter("v3","v3",6e-4,-1.0,1.0);
        ft->AddParameter("v4","v4",1.8e-4,-1.0,1.0);

        // pp initial value
        // ft->AddParameter("Fa","Fa",0.1,0,30);
        // ft->AddParameter("Ga","Ga",1,0,30);
        // ft->AddParameter("v2","v2",4e-3,-1.0,1.0);
        // ft->AddParameter("v3","v3",6e-4,-1.0,1.0);
        // ft->AddParameter("v4","v4",1.8e-4,-1.0,1.0);
        // // ft->AddParameter("v5","v5",1e-4,0,1.0);
    }
    else{
        ft->AddParameter("Fa","Fa",fParamVal[3],fParamVal[3]-1,fParamVal[3]+1);
        ft->AddParameter("Ga","Ga",fParamVal[4],fParamVal[4]-5,fParamVal[4]+5);
        ft->AddParameter("v1","v1",fParamVal[5],fParamVal[5]-0.0002,fParamVal[5]+0.0002);
        ft->AddParameter("v2","v2",fParamVal[0],fParamVal[0]-0.0002,fParamVal[0]+0.0002);
        ft->AddParameter("v3","v3",fParamVal[1],fParamVal[1]-0.0002,fParamVal[1]+0.0002);
        ft->AddParameter("v4","v4",fParamVal[2],fParamVal[2]-0.0002,fParamVal[2]+0.0002);
    }

    //Construct fit function
    FunctionObject *fobj = new TemplateFunction(lm);
    ft->SetFitFunction(fobj);
    //Perform fit:
    printf("About to fit\n");
    Int_t dummy = ft->Fit(0); //Do not draw performance at this point. Return value of Fit() is false if no base is set.
    if(!dummy) return;

    Double_t F   = ft->getVal(0); //0 for F, 1 for G, 2 for v2, 3 for v3, 4 for v4
    Double_t Fe  = ft->getErr(0);
    Double_t G   = ft->getVal(1); //0 for F, 1 for G, 2 for v2, 3 for v3, 4 for v4
    Double_t Ge  = ft->getErr(1);
    Double_t v11 = ft->getVal(2); //0 for F, 1 for G, 2 for v1, 3 for v2, 4 for v3, 5 for v4
    Double_t v11e= ft->getErr(2);
    Double_t v21 = ft->getVal(3); //0 for F, 1 for G, 2 for v1, 3 for v2, 4 for v3, 5 for v4
    Double_t v21e= ft->getErr(3);
    Double_t v31 = ft->getVal(4); //0 for F, 1 for G, 2 for v1, 3 for v2, 4 for v3, 5 for v4
    Double_t v31e= ft->getErr(4);
    Double_t v41 = ft->getVal(5); //0 for F, 1 for G, 2 for v1, 3 for v2, 4 for v3, 5 for v4
    Double_t v41e= ft->getErr(5);

    printf("Values from fit:\n");
    printf("F  = %f +- %f\n",F,Fe);
    printf("G  = %f +- %f\n",G,Ge);
    printf("V1 = %f +- %f\n",v11,v11e);
    printf("V2 = %f +- %f\n",v21,v21e);
    printf("V3 = %f +- %f\n",v31,v31e);
    printf("V4 = %f +- %f\n",v41,v41e);

    fParamVal[0] = v21; fParamErr[0] = v21e;
    fParamVal[1] = v31; fParamErr[1] = v31e;
    fParamVal[2] = v41; fParamErr[2] = v41e;
    fParamVal[3] = F;   fParamErr[3] = Fe;
    fParamVal[4] = G;   fParamErr[4] = Ge;
    fParamVal[5] = v11; fParamErr[5] = v11e;
}

void DrawText(double xmin, double ymin, double textSize, TString text)
{

	TLatex *textPreliminary = new TLatex(xmin, ymin, Form("%s", text.Data()));
	textPreliminary->SetNDC();
	textPreliminary->SetTextFont(43);
	textPreliminary->SetTextSize(textSize);
	textPreliminary->Draw();

}

void PlotFitting(TH1 *lm, TH1 *hm, Bool_t isNch, std::string fileSuffix, Int_t minRange, Int_t maxRange, const std::vector<Double_t>& par, const std::vector<Double_t>& parerr, Double_t pTMin=0, Double_t pTMax=0, const std::string& systemName = "") {
    gStyle->SetOptStat(0); 
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    
    double v21 = par[0];
    double v31 = par[1];
    double v41 = par[2];
    double F =   par[3];
    double G =   par[4];
    double v11 = par[5];

    double v21e = parerr[0];
    double v31e = parerr[1];
    double v41e = parerr[2];
    double v11e = parerr[5];


    TCanvas* canvas = new TCanvas(Form("Fit"), "Fit", 800, 600);
    canvas->Range(0,0,1,1);
    
    // 创建上下面板
    TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
    pad1->SetBorderMode(0);
    pad1->SetBorderSize(2);
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.05);
    pad1->SetTopMargin(0.05);
    pad1->SetBottomMargin(0.13);
    pad1->SetTicks(1,1);
    pad1->Draw();
    
    TPad* pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.25);
    pad2->SetBorderMode(0);
    pad2->SetBorderSize(2);
    pad2->SetLeftMargin(0.12);
    pad2->SetRightMargin(0.05);
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.25);
    pad2->SetTicks(1,1);
    pad2->Draw();

    // 绘制主图
    pad1->cd();
    
    // 创建背景直方图
    TH1D* hbkg1 = new TH1D("dPhi", "dPhi", 1, -TMath::Pi()/2.0, 3*TMath::Pi()/2.0);
    hbkg1->SetStats(0);
    hbkg1->GetXaxis()->SetTitle("#Delta#phi [rad]");
    hbkg1->GetYaxis()->SetTitle("Arbitrary scaled");
    hbkg1->GetXaxis()->SetTitleSize(0.045);
    hbkg1->GetXaxis()->SetLabelSize(0.04);
    hbkg1->GetYaxis()->SetTitleSize(0.045);
    hbkg1->GetYaxis()->SetLabelSize(0.04);
    hbkg1->GetXaxis()->SetTitleOffset(1.2);
    hbkg1->GetYaxis()->SetTitleOffset(1.2);
    
    // 设置Y轴范围
    double ymax = hm->GetMaximum();
    double ymin = hm->GetMinimum();
    double ydelta = ymax - ymin;
    hbkg1->GetYaxis()->SetRangeUser(ymin-ydelta*0.2, ymax+ydelta*0.7);
    hbkg1->Draw();

    // 绘制原始数据
    hm->SetMarkerStyle(20);
    hm->SetMarkerColor(kBlack);
    hm->SetLineColor(kBlack);
    hm->SetMarkerSize(1.0);
    hm->Draw("same p");

    // 创建拟合曲线
    const int pointBin = (int)hm->GetNbinsX();
    Double_t CopyPointX[pointBin];
    Double_t CopyPointY[pointBin];
    for (int i=0; i<pointBin; ++i){
      CopyPointX[i] = hm->GetBinCenter(i+1);
      double x = hm->GetBinCenter(i+1);
      // CopyPointY[i] = F*lm->GetBinContent(i+1)+G*(1+2*v21*cos(2*x)+2*v31*cos(3*x));
      CopyPointY[i] = F*lm->GetBinContent(i+1)+G*(1+2*v11*cos(x)+2*v21*cos(2*x)+2*v31*cos(3*x)+2*v41*cos(4*x));
    };
    TGraph* gCopy = new TGraph(pointBin,CopyPointX,CopyPointY);

    Double_t PeriPointY[pointBin];
    for (int i=0; i<pointBin; ++i){
      PeriPointY[i] = F*lm->GetBinContent(i+1) + G;
    };
    TGraph* gPeri = new TGraph(pointBin,CopyPointX,PeriPointY);

    Double_t Y0Position = lm->GetXaxis()->FindBin(0.0);
    Double_t Y0 = lm->GetBinContent(Y0Position);


    gCopy->SetLineColor(colors[0]);
    gCopy->SetLineWidth(2);
    gCopy->Draw("same");

    gPeri->SetLineColor(colors[1]);
    gPeri->SetLineWidth(2);
    gPeri->SetLineStyle(5);
    gPeri->Draw("same");


    TF1* fit_p1 = new TF1("fit_p1","[0]*[1] + [2]*(1 + 2*[3]*cos(x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    TF1* fit_p2 = new TF1("fit_p2","[0]*[1] + [2]*(1 + 2*[3]*cos(2*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    TF1* fit_p3 = new TF1("fit_p3","[0]*[1] + [2]*(1 + 2*[3]*cos(3*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    TF1* fit_p4 = new TF1("fit_p4","[0]*[1] + [2]*(1 + 2*[3]*cos(4*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    fit_p1->SetParameters(F,Y0,G,v11);
    fit_p2->SetParameters(F,Y0,G,v21);
    fit_p3->SetParameters(F,Y0,G,v31);
    fit_p4->SetParameters(F,Y0,G,v41);
    fit_p1->SetLineColor(colors[1]);
    fit_p1->SetLineWidth(2);
    fit_p1->SetLineStyle(1);
    fit_p2->SetLineColor(colors[2]);
    fit_p2->SetLineWidth(2);
    fit_p2->SetLineStyle(2);
    fit_p3->SetLineColor(colors[3]);
    fit_p3->SetLineWidth(2);
    fit_p3->SetLineStyle(3);
    fit_p4->SetLineColor(colors[4]);
    fit_p4->SetLineWidth(2);
    fit_p4->SetLineStyle(4);
    fit_p1->Draw("same");
    fit_p2->Draw("same");
    fit_p3->Draw("same");
    fit_p4->Draw("same");

    // 添加图例
    TLegend* leg = new TLegend(0.5, 0.65, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->AddEntry(hm, "Data", "lep");
    leg->AddEntry(gCopy, "FY(#Delta#phi)^{peri} + G(1+#Sigma_{n=1}^{4}2V_{n#Delta}cos(n#Delta#phi))", "l");
    leg->AddEntry(gPeri, "FY(#Delta#phi)^{peri} + G", "l");
    leg->AddEntry(fit_p1, "FY(0)^{peri} + G(1+2V_{1#Delta}cos(#Delta#phi))", "l");
    leg->AddEntry(fit_p2, "FY(0)^{peri} + G(1+2V_{2#Delta}cos(2#Delta#phi))", "l");
    leg->AddEntry(fit_p3, "FY(0)^{peri} + G(1+2V_{3#Delta}cos(3#Delta#phi))", "l");
    leg->AddEntry(fit_p4, "FY(0)^{peri} + G(1+2V_{4#Delta}cos(4#Delta#phi))", "l");
    leg->Draw();

    // 添加文本标签
    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(43);
    tex->SetTextSize(20);
    std::string displaySystemName = systemName.empty() ? collisionSystemName : systemName;
    DrawText(0.2, 0.85, 20, Form("ALICE %s", displaySystemName.c_str()));
    if (!gIsEtaDiff) {
        DrawText(0.2, 0.80, 20, Form("%.1f < |#Delta#eta| < %.1f", 1.0, 1.6));
    }
    if (isNch) {
        DrawText(0.2, 0.75, 20, Form("%d < N_{ch} < %d", minRange, maxRange));
    }
    else {
        DrawText(0.2, 0.75, 20, Form("%d < Cent < %d", minRange, maxRange));
    }

    if (gIsEtaDiff || (pTMin > 0 && pTMax > 0)) {
        if (gIsEtaDiff) {
            DrawText(0.2, 0.70, 20, Form("%.1f < #eta_{trig} < %.1f", pTMin, pTMax));
        } else {
            DrawText(0.2, 0.70, 20, Form("%.1f < p_{T} < %.1f", pTMin, pTMax));
        }
    }
    DrawText(0.2, 0.65, 20, Form("V_{1#Delta} = %.5f #pm %.5f", v11, v11e));
    DrawText(0.2, 0.60, 20, Form("V_{2#Delta} = %.5f #pm %.5f", v21, v21e));
    DrawText(0.2, 0.55, 20, Form("V_{3#Delta} = %.5f #pm %.5f", v31, v31e));

    // 绘制底部残差图
    pad2->cd();

    TF1* fit_p1m = new TF1("fit_p1m","[0]*(1 + 2*[1]*cos(x) + 2*[2]*cos(2*x) + 2*[3]*cos(3*x) + 2*[4]*cos(4*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    fit_p1m->SetParameters(G,v11,v21,v31,v41);

    
    // 创建残差直方图
    TH1D* hsubtract = (TH1D*)hm->Clone(Form("subtract"));
    hsubtract->Add(lm, -F);

    // =============== 新增：计算chi2/ndf ===============
    double chi2 = 0.0;
    int nBins = hsubtract->GetNbinsX();
    int nParams = 6; // 参数个数: F,G,v11,v21,v31,v41
    int ndf = nBins - nParams;
    
    for (int i = 1; i <= nBins; i++) {
        double data = hsubtract->GetBinContent(i);
        double error = hsubtract->GetBinError(i);
        double x = hsubtract->GetBinCenter(i);
        double fit = fit_p1m->Eval(x);
        
        if (error > 0) { // 忽略误差为0的bin
            double residual = data - fit;
            chi2 += (residual * residual) / (error * error);
        }
    }
    double chi2ndf = (ndf > 0) ? chi2 / ndf : 0;
    // =============== 新增：添加chi2/ndf标签 ===============
    pad1->cd();
    TLatex* chi2Label = new TLatex();
    chi2Label->SetNDC();
    chi2Label->SetTextFont(43);
    chi2Label->SetTextSize(20);
    chi2Label->DrawLatex(0.50, 0.60, Form("#chi^{2}/ndf = %.1f/%d = %.2f", chi2, ndf, chi2ndf));

    pad2->cd();
    // TH1D* hResidual = (TH1D*)hsubtract->Clone("residual");
    // for(int i=1; i<=hResidual->GetNbinsX(); ++i){
    //   double data = hm->GetBinContent(i);
    //   double fit = fit_p1m->Eval(hm->GetBinCenter(i));
    //   hResidual->SetBinContent(i, data/fit);
    //   hResidual->SetBinError(i, hm->GetBinError(i)/fit);
    // }
    TH1D* hResidual = (TH1D*)hsubtract->Clone(Form("pull"));
    double ymax_pull = 1, ymin_pull = 1;
    for (int i = 1; i <= hResidual->GetXaxis()->GetNbins(); i++) {
      double dat = hResidual->GetBinContent(i);
      double err = hResidual->GetBinError(i);
      // double lin = gCopy->GetPointY(i-1);
      double lin = fit_p1m->Eval(hResidual->GetBinCenter(i));
      hResidual->SetBinContent(i, dat/lin);
      hResidual->SetBinError(i, err/lin);
      ymin_pull = ymin_pull>(dat/lin)?(dat/lin):ymin_pull;
      ymax_pull = ymax_pull<(dat/lin)?(dat/lin):ymax_pull;
    }
    ymin_pull = 1 - (1-ymin_pull) * 2;
    ymax_pull = 1 + (ymax_pull-1) * 2;
    hResidual->GetYaxis()->SetRangeUser(ymin_pull,ymax_pull);

    // 配置残差图
    hResidual->SetMarkerStyle(20);
    hResidual->SetMarkerColor(kBlack);
    hResidual->SetLineColor(kBlack);
    hResidual->GetYaxis()->SetTitle("Data/Fit");
    hResidual->GetYaxis()->SetTitleSize(0.12);
    hResidual->GetYaxis()->SetLabelSize(0.10);
    hResidual->GetYaxis()->SetTitleOffset(0.35);
    hResidual->GetXaxis()->SetTitleSize(0.12);
    hResidual->GetXaxis()->SetTitle("#Delta#phi [rad]");
    hResidual->GetXaxis()->SetLabelSize(0.10);
    hResidual->GetXaxis()->SetTitleOffset(0.9);
    hResidual->Draw("ep");

    // 添加参考线
    TLine* line = new TLine(-TMath::Pi()/2.0, 1.0, 1.5*TMath::Pi(), 1.0);
    line->SetLineColor(kBlack);
    line->SetLineStyle(kDashed);
    line->Draw();

    // 保存结果
    if (gIsEtaDiff || (pTMin > 0 && pTMax > 0)) {
        if (gIsEtaDiff) {
            gSystem->mkdir("./TemplateFit/EtaDiff/PDFs", kTRUE);
            canvas->SaveAs(Form("./TemplateFit/EtaDiff/PDFs/TemplateFit_%s_%s_%d_%d_Eta_%0.1f_%0.1f.pdf", fileSuffix.c_str(), splitName.c_str(), minRange, maxRange, pTMin, pTMax));
        } else {
            gSystem->mkdir("./TemplateFit/PtDiff/PDFs", kTRUE);
            canvas->SaveAs(Form("./TemplateFit/PtDiff/PDFs/TemplateFit_%s_%s_%d_%d_Pt_%0.1f_%0.1f.pdf", fileSuffix.c_str(), splitName.c_str(), minRange, maxRange, pTMin, pTMax));
        }
    } else {
        gSystem->mkdir("./TemplateFit/PDFs", kTRUE);
        canvas->SaveAs(Form("./TemplateFit/PDFs/TemplateFit_%s_%s_%d_%d.pdf", fileSuffix.c_str(), splitName.c_str(), minRange, maxRange));
    }
    
    // 清理内存
    delete hbkg1;
    delete gCopy;
    delete gPeri;
    delete fit_p2;
    delete fit_p3;
    delete fit_p4;
    delete leg;
    delete tex;
    delete hResidual;
    delete line;
  
}
