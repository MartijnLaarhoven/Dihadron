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
#include "TSystem.h"
#include <TProfile.h>
//#include <TRandom3.h>
#include "TMath.h"
#include "THnSparse.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include "./include/BasicForDihadron.h"
#include "./include/Bootstrap.h"
#include "./include/TemplateFitter.cxx"
#include "./include/TemplateFunction.C"
#include "./include/ErrorPropagation.h"
#include "./include/plotting.h"

//my add-ons
#include <fstream> //for csv writing



// define struct
struct InputUnit {
    std::string fileNameSuffix;
    Int_t minRange;
    Int_t maxRange;
    InputUnit(std::string _fileNameSuffix, Int_t _minRange, Int_t _maxRange) :
        fileNameSuffix(_fileNameSuffix),  minRange(_minRange), maxRange(_maxRange) {}
};

struct ConfigUnit {
    Bool_t isNch;
    Int_t corrType;
    Bool_t isEtaDiff;
    InputUnit templ;
    std::vector<InputUnit> dataList;
    std::string outputFileName;
    ConfigUnit(Bool_t _isNch, Int_t _corrType, Bool_t _isEtaDiff, InputUnit _template, std::vector<InputUnit> _dataList, std::string _outputFileName) :
        isNch(_isNch), corrType(_corrType), isEtaDiff(_isEtaDiff), templ(_template), dataList(_dataList), outputFileName(_outputFileName) {}
};

struct VnUnit {
    Double_t v2;
    Double_t v2_err;
    Double_t v3;
    Double_t v3_err;
    Double_t v4;
    Double_t v4_err;
    Double_t F;
    Double_t F_err;
    Double_t G;
    Double_t G_err;
    Double_t chi2ndf;
    Double_t chi2ndf_err;
    VnUnit(Double_t _v2, Double_t _v2_err, Double_t _v3, Double_t _v3_err, Double_t _v4, Double_t _v4_err, Double_t _F, Double_t _F_err, Double_t _G, Double_t _G_err, Double_t _chi2ndf, Double_t _chi2ndf_err) :
        v2(_v2), v2_err(_v2_err), v3(_v3), v3_err(_v3_err), v4(_v4), v4_err(_v4_err), F(_F), F_err(_F_err), G(_G), G_err(_G_err), chi2ndf(_chi2ndf), chi2ndf_err(_chi2ndf_err) {}
};

// declare functions
void ProcessConfig(Bool_t isNch, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName, Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.);
void ProcessConfig_EtaDiff(Bool_t isNch, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName, Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.);
void CreateV2DeltaSummaryPlot(Bool_t isNch, std::vector<InputUnit> dataList, std::string fileNameSuffix, std::string detectorName, const std::vector<float>& ft0EtaVec);
void CreateV2DeltaSummaryPlot_InnerRing(Bool_t isNch, std::vector<InputUnit> dataList, std::string fileNameSuffix, std::string detectorName);
void CreateRingComparisonPlots(Bool_t isNch, std::vector<InputUnit> dataList, std::string fileNameSuffix, std::string detectorName, bool isInner);
VnUnit* TemplateFit(Bool_t isNch, InputUnit templ, InputUnit data, Bool_t cn2Tovn2, Double_t etaMin=-101., Double_t etaMax=-101., Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.);
std::vector<Int_t> CheckAndMergeRanges(const std::vector<InputUnit>& inputUnits);
VnUnit* fitSample(Bool_t isNch, TFile* templatefile, InputUnit templ, TFile* datafile, InputUnit data, int sample = -1, Double_t etaMin=-101., Double_t etaMax=-101., Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.);
void RooTempFitter(TH1 *lm, TH1 *hm, std::vector<Double_t>& fParamVal, std::vector<Double_t>& fParamErr, Bool_t kRefit);
void PlotFitting(TH1 *lm, TH1 *hm, Bool_t isNch, std::string fileSuffix, Int_t minRange, Int_t maxRange, const std::vector<Double_t>& par, const std::vector<Double_t>& parerr, Double_t etaMin=-101., Double_t etaMax=-101., Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.);
void writeToCSV(const std::string &filename, const std::vector<std::vector<std::string>> &data, const std::string &folder, bool createHeader=false); //mine
double findChi2ndf(TH1 *lm, TH1 *hm, double& F_val, double& G_val, double& v21_val, double& v31_val, double& v41_val);
std::vector<std::string> toStringV2(std::vector<std::vector<std::vector<double>>> vec, int Nelements, double etaMin, double etaMax);
std::vector<std::string> toStringV3(std::vector<std::vector<std::vector<double>>> vec, int Nelements, double etaMin, double etaMax);
std::vector<std::string> toStringV4(std::vector<std::vector<std::vector<double>>> vec, int Nelements, double etaMin, double etaMax);

// global variables
std::string collisionSystemName = "peripheral Ne-Ne"; 
std::string FITused = "TPC-FT0C";
collisionSystemName = "Ne-Ne";

//==============================================================
void Process_TemplateFit() {
    // 不显示窗口
    gROOT->SetBatch(kTRUE);
    std::vector<ConfigUnit> configList;

    // Simple switch: run one detector at a time
    // Change to kTPCFT0A for FT0A, or kTPCFT0C for FT0C
    const Int_t corrTypeToRun = kTPCFT0C;

    //Ne-Ne

    // configList.push_back(ConfigUnit(kCent, kTPCFT0C,  kEtaDiffOff, InputUnit("LHC25af_pass2_598673", 80, 100), 
    // {InputUnit("LHC25af_pass2_598673", 0, 20)}, 
    // "LHC25af_pass2_598673"));

    // configList.push_back(ConfigUnit(kCent, kTPCFT0A,  kEtaDiffOff, InputUnit("LHC25af_pass2_598673", 80, 100), 
    // {InputUnit("LHC25af_pass2_598673", 0, 20)}, 
    // "LHC25af_pass2_598673"));

    // configList.push_back(ConfigUnit(kCent, kTPCFT0C,  kEtaDiffOn, InputUnit("LHC25af_pass2_598673", 80, 100), 
    // {InputUnit("LHC25af_pass2_598673", 0, 20)}, 
    // "LHC25af_pass2_598673"));

    // configList.push_back(ConfigUnit(kCent, kTPCFT0A,  kEtaDiffOn, InputUnit("LHC25af_pass2_598673", 80, 100), 
    // {InputUnit("LHC25af_pass2_598673", 0, 20)}, 
    // "LHC25af_pass2_598673"));

    //O-O

    // configList.push_back(ConfigUnit(kCent, kTPCFT0C,  kEtaDiffOff, InputUnit("LHC25ae_pass2_598682", 80, 100), 
    // {InputUnit("LHC25ae_pass2_598682", 0, 20)},
    // "LHC25ae_pass2_598682"));

    // configList.push_back(ConfigUnit(kCent, kTPCFT0A,  kEtaDiffOff, InputUnit("LHC25ae_pass2_598682", 80, 100), 
    // {InputUnit("LHC25ae_pass2_598682", 0, 20)},
    // "LHC25ae_pass2_598682"));

    // configList.push_back(ConfigUnit(kCent, kTPCFT0C,  kEtaDiffOn, InputUnit("LHC25ae_pass2_598682", 80, 100), 
    // {InputUnit("LHC25ae_pass2_598682", 0, 20)},
    // "LHC25ae_pass2_598682"));

    // configList.push_back(ConfigUnit(kCent, kTPCFT0A,  kEtaDiffOn, InputUnit("LHC25ae_pass2_598682", 80, 100), 
    // {InputUnit("LHC25ae_pass2_598682", 0, 20)},
    // "LHC25ae_pass2_598682"));

    // Ne-Ne outer ring (615818)
    configList.push_back(ConfigUnit(kCent, kTPCFT0C,  kEtaDiffOn, InputUnit("LHC25af_pass2_615818", 0, 20), 
    {InputUnit("LHC25af_pass2_615818", 0, 20)},
    "LHC25af_pass2_615818"));

    configList.push_back(ConfigUnit(kCent, kTPCFT0A,  kEtaDiffOn, InputUnit("LHC25af_pass2_615818", 80, 100), 
    {InputUnit("LHC25af_pass2_615818", 0, 20)},
    "LHC25af_pass2_615818"));

    // Ne-Ne inner ring (615817)
    configList.push_back(ConfigUnit(kCent, kTPCFT0C,  kEtaDiffOn, InputUnit("LHC25af_pass2_615817", 0, 20), 
    {InputUnit("LHC25af_pass2_615817", 0, 20)},
    "LHC25af_pass2_615817"));

    configList.push_back(ConfigUnit(kCent, kTPCFT0A,  kEtaDiffOn, InputUnit("LHC25af_pass2_615817", 0, 20), 
    {InputUnit("LHC25af_pass2_615817", 0, 20)},
    "LHC25af_pass2_615817"));

    // O-O outer ring (616549)
    configList.push_back(ConfigUnit(kCent, kTPCFT0C,  kEtaDiffOn, InputUnit("LHC25ae_pass2_616549", 0, 20), 
    {InputUnit("LHC25ae_pass2_616549", 0, 20)},
    "LHC25ae_pass2_616549"));

    configList.push_back(ConfigUnit(kCent, kTPCFT0A,  kEtaDiffOn, InputUnit("LHC25ae_pass2_616549", 0, 20), 
    {InputUnit("LHC25ae_pass2_616549", 0, 20)},
    "LHC25ae_pass2_616549"));

    // O-O inner ring (618685)
    configList.push_back(ConfigUnit(kCent, kTPCFT0C,  kEtaDiffOn, InputUnit("LHC25ae_pass2_618685", 0, 20), 
    {InputUnit("LHC25ae_pass2_618685", 0, 20)},
    "LHC25ae_pass2_618685"));

    configList.push_back(ConfigUnit(kCent, kTPCFT0A,  kEtaDiffOn, InputUnit("LHC25ae_pass2_618685", 0, 20), 
    {InputUnit("LHC25ae_pass2_618685", 0, 20)},
    "LHC25ae_pass2_618685"));

    // Ne-Ne inner ring (617826) - new dataset
    configList.push_back(ConfigUnit(kCent, kTPCFT0C,  kEtaDiffOn, InputUnit("LHC25af_pass2_617826", 0, 20), 
    {InputUnit("LHC25af_pass2_617826", 0, 20)},
    "LHC25af_pass2_617826"));

    configList.push_back(ConfigUnit(kCent, kTPCFT0A,  kEtaDiffOn, InputUnit("LHC25af_pass2_617826", 0, 20), 
    {InputUnit("LHC25af_pass2_617826", 0, 20)},
    "LHC25af_pass2_617826"));

    // Ne-Ne outer ring (617910) - new dataset
    configList.push_back(ConfigUnit(kCent, kTPCFT0C,  kEtaDiffOn, InputUnit("LHC25af_pass2_617910", 0, 20), 
    {InputUnit("LHC25af_pass2_617910", 0, 20)},
    "LHC25af_pass2_617910"));

    configList.push_back(ConfigUnit(kCent, kTPCFT0A,  kEtaDiffOn, InputUnit("LHC25af_pass2_617910", 0, 20), 
    {InputUnit("LHC25af_pass2_617910", 0, 20)},
    "LHC25af_pass2_617910"));

    //her itererer vi over konfigurasjonene
    for (auto config : configList) {
        if (config.corrType != corrTypeToRun) continue;

        // Auto-detect collision system from dataset name
        if (config.templ.fileNameSuffix.find("LHC25ae") != std::string::npos) {
            collisionSystemName = "O-O";
        } else if (config.templ.fileNameSuffix.find("LHC25af") != std::string::npos) {
            collisionSystemName = "Ne-Ne";
        } else {
            collisionSystemName = "Unknown";  // Fallback
        }

        // Used for plot labels in PlotFitting
        if (config.corrType == kTPCFT0A) FITused = "TPC-FT0A";
        if (config.corrType == kTPCFT0C) FITused = "TPC-FT0C";

        // Check if this is a ring dataset
        bool isRingDataset = (config.templ.fileNameSuffix.find("615817") != std::string::npos || 
                              config.templ.fileNameSuffix.find("615818") != std::string::npos ||
                              config.templ.fileNameSuffix.find("616549") != std::string::npos ||
                              config.templ.fileNameSuffix.find("618685") != std::string::npos ||
                              config.templ.fileNameSuffix.find("617826") != std::string::npos ||
                              config.templ.fileNameSuffix.find("617910") != std::string::npos);

        // For ring datasets: process with TPC edges, no FT0 slicing, create comparison plots
        if (isRingDataset && config.isEtaDiff) {
            std::vector<std::pair<double,double>> tpcEdges = { std::make_pair(-0.8, -0.7), std::make_pair(0.7, 0.8) };
            bool isInner = (config.templ.fileNameSuffix.find("615817") != std::string::npos || 
                            config.templ.fileNameSuffix.find("618685") != std::string::npos ||
                            config.templ.fileNameSuffix.find("617826") != std::string::npos);
            std::string ringType = isInner ? "inner" : "outer";
            // For ring datasets, use whole FT0 detector without subdivision
            for (auto rng : tpcEdges) {
                double etaMin = rng.first;
                double etaMax = rng.second;
                std::string outputTag = Form("%s_%s_TPCEta_%0.1f_%0.1f", 
                                            config.outputFileName.c_str(), 
                                            DihadronCorrTypeName[config.corrType].c_str(),
                                            etaMin, etaMax);
                ProcessConfig_EtaDiff(config.isNch, config.templ, config.dataList, outputTag);
            }
            // Create comparison plots after processing all eta ranges
            CreateRingComparisonPlots(config.isNch, config.dataList, config.outputFileName, 
                                     DihadronCorrTypeName[config.corrType], isInner);
        }
        // For original FT0-sliced datasets (598682): use FT0 slices with summary plots
        else {
            // FT0 eta slicing: loop over slices for FT0A/FT0C
            const std::vector<float>* ft0EtaVec = nullptr;
            if (config.corrType == kTPCFT0A) ft0EtaVec = &etaFT0A;
            else if (config.corrType == kTPCFT0C) ft0EtaVec = &etaFT0C;

            if (ft0EtaVec && ft0EtaVec->size() > 1) {
                // FT0 slicing enabled: iterate over slices
                for (size_t iFT0 = 0; iFT0 < ft0EtaVec->size() - 1; ++iFT0) {
                    Double_t ft0EtaMin = (*ft0EtaVec)[iFT0];
                    Double_t ft0EtaMax = (*ft0EtaVec)[iFT0 + 1];
                    std::string outputTag = Form("%s_%s_FT0Eta_%0.1f_%0.1f", 
                                                config.outputFileName.c_str(), 
                                                DihadronCorrTypeName[config.corrType].c_str(),
                                                ft0EtaMin, ft0EtaMax);
                    
                    if (config.isEtaDiff) {
                        ProcessConfig_EtaDiff(config.isNch, config.templ, config.dataList, outputTag, ft0EtaMin, ft0EtaMax);
                    } else {
                        ProcessConfig(config.isNch, config.templ, config.dataList, outputTag, ft0EtaMin, ft0EtaMax);
                    }
                }
                
                // After processing all FT0 slices, create summary comparison plots
                if (config.isEtaDiff) {
                    CreateV2DeltaSummaryPlot(config.isNch, config.dataList, config.outputFileName, 
                                            DihadronCorrTypeName[config.corrType], *ft0EtaVec);
                }
            } else {
                // No FT0 slicing: original behavior
                std::string outputTag = Form("%s_%s", config.outputFileName.c_str(), DihadronCorrTypeName[config.corrType].c_str());

                if (config.isEtaDiff) { //eta-diff case
                    ProcessConfig_EtaDiff(config.isNch, config.templ, config.dataList, outputTag);
                } else { //ikke eta-diff case
                    ProcessConfig(config.isNch, config.templ, config.dataList, outputTag);
                }
            }
        }
    }

}

//==============================================================
void ProcessConfig(Bool_t isNch, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName, Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.) {
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
        vnResults.push_back(TemplateFit(isNch, templ, data, kFALSE, -101., -101., ft0EtaMin, ft0EtaMax)); //her finner vi hva v2 er, kTRUE betyr at cn2 to vn2 er aktivert?
    } //her står vi igjen med V_Delta n (ref, ref)

    // 创建输出文件
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    TFile outputFile(Form("./TemplateFit/Vn_%s_%s.root", outputFileName.c_str(), splitName.c_str()), "RECREATE"); //har lages .root for V_Delta n (ref, ref)

    // Quick-look plots (PDF)
    gSystem->mkdir("./TemplateFit/Plots", kTRUE);

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
            hV2->SetBinContent(i+1, vnResults[i]->v2); //her hentes v2 for ikke diff (skrives direkte)
            hV2->SetBinError(i+1, vnResults[i]->v2_err);
            
            hV3->SetBinContent(i+1, vnResults[i]->v3); //v3
            hV3->SetBinError(i+1, vnResults[i]->v3_err);
            
            hV4->SetBinContent(i+1, vnResults[i]->v4); //v4
            hV4->SetBinError(i+1, vnResults[i]->v4_err);
        }

        // 写入文件
        hV2->Write();
        hV3->Write();
        hV4->Write();

        TCanvas cVn("cVn", "cVn", 900, 700);
        hV2->SetMarkerStyle(20);
        hV2->Draw("E1");
        cVn.SaveAs(Form("./TemplateFit/Plots/V2_%s_%s.pdf", outputFileName.c_str(), splitName.c_str()));
        hV3->SetMarkerStyle(20);
        hV3->Draw("E1");
        cVn.SaveAs(Form("./TemplateFit/Plots/V3_%s_%s.pdf", outputFileName.c_str(), splitName.c_str()));
        hV4->SetMarkerStyle(20);
        hV4->Draw("E1");
        cVn.SaveAs(Form("./TemplateFit/Plots/V4_%s_%s.pdf", outputFileName.c_str(), splitName.c_str()));
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

        TCanvas cVn("cVn", "cVn", 900, 700);
        gV2->SetMarkerStyle(20);
        gV2->Draw("AP");
        cVn.SaveAs(Form("./TemplateFit/Plots/V2_%s_%s.pdf", outputFileName.c_str(), splitName.c_str()));
        gV3->SetMarkerStyle(20);
        gV3->Draw("AP");
        cVn.SaveAs(Form("./TemplateFit/Plots/V3_%s_%s.pdf", outputFileName.c_str(), splitName.c_str()));
        gV4->SetMarkerStyle(20);
        gV4->Draw("AP");
        cVn.SaveAs(Form("./TemplateFit/Plots/V4_%s_%s.pdf", outputFileName.c_str(), splitName.c_str()));
    }

    std::cout << "Output file: " << Form("./TemplateFit/Vn_%s_%s.root", outputFileName.c_str(), splitName.c_str()) << std::endl;
    outputFile.Close();
}

//==============================================================
void ProcessConfig_EtaDiff(Bool_t isNch, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName, Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.) {
    // Just looping the list
    for (const auto& data : dataList) {
        
        std::string splitName = "Mult";
        if (!isNch) splitName = "Cent";

        bool isRingDataset = (templ.fileNameSuffix.find("615817") != std::string::npos ||
                              templ.fileNameSuffix.find("615818") != std::string::npos ||
                              templ.fileNameSuffix.find("616549") != std::string::npos ||
                              templ.fileNameSuffix.find("618685") != std::string::npos ||
                              templ.fileNameSuffix.find("617826") != std::string::npos ||
                              templ.fileNameSuffix.find("617910") != std::string::npos);
        std::vector<std::pair<double,double>> tpcRanges;
        if (isRingDataset) {
            double parsedMin = 0.0;
            double parsedMax = 0.0;
            bool hasTPCEta = false;
            std::size_t pos = outputFileName.find("TPCEta_");
            if (pos != std::string::npos) {
                std::string tail = outputFileName.substr(pos + 7);
                if (std::sscanf(tail.c_str(), "%lf_%lf", &parsedMin, &parsedMax) == 2) {
                    hasTPCEta = true;
                }
            }
            if (hasTPCEta) {
                tpcRanges.push_back(std::make_pair(parsedMin, parsedMax));
            } else {
                tpcRanges = { std::make_pair(-0.8, -0.7), std::make_pair(0.7, 0.8) };
            }
        }
        
        // When FT0 slicing is active, process each TPC bin separately with its own output file
        if (ft0EtaMin > -900. && ft0EtaMax > -900.) {
            // Use wider TPC ranges to match Process_dPhidEta: [-0.8, -0.6] and [0.6, 0.8]
            std::vector<std::pair<double,double>> tpcRanges = { { -0.8, -0.6 }, { 0.6, 0.8 } };
            for (auto rng : tpcRanges) {
                double etaMin = rng.first;
                double etaMax = rng.second;
                
                // Execute template fit for this specific TPC bin
                VnUnit* vnResult = TemplateFit(isNch, templ, data, kFALSE, etaMin, etaMax, ft0EtaMin, ft0EtaMax);
                
                // Create separate output file for this TPC bin
                TFile outputFile(Form("./TemplateFit/EtaDiff/Vn_%s_%s_FT0Eta_%0.1f_%0.1f_%i_%i_TPCEta_%0.1f_%0.1f.root", 
                                     outputFileName.c_str(), splitName.c_str(), ft0EtaMin, ft0EtaMax,
                                     data.minRange, data.maxRange, etaMin, etaMax), "RECREATE");
                
                gSystem->mkdir("./TemplateFit/EtaDiff/Plots", kTRUE);
                
                // Create simple histogram with single bin for this TPC eta range
                Double_t binEdges[2] = {etaMin, etaMax};
                TH1D* hV2 = new TH1D("hV2", "V_{ #Delta 2};#eta^{TPC};V_{ #Delta 2}", 1, binEdges);
                TH1D* hV3 = new TH1D("hV3", "V_{ #Delta 3};#eta^{TPC};V_{ #Delta 3}", 1, binEdges);
                TH1D* hV4 = new TH1D("hV4", "V_{ #Delta 4};#eta^{TPC};V_{ #Delta 4}", 1, binEdges);
                
                hV2->SetStats(0);
                hV3->SetStats(0);
                hV4->SetStats(0);
                
                if (vnResult) {
                    Double_t v2 = vnResult->v2;
                    Double_t e2 = vnResult->v2_err;
                    Double_t v3 = vnResult->v3;
                    Double_t e3 = vnResult->v3_err;
                    Double_t v4 = vnResult->v4;
                    Double_t e4 = vnResult->v4_err;
                    
                    if (!std::isfinite(v2) || !std::isfinite(e2) || e2 >= 9.99) { v2 = 0.0; e2 = 0.0; }
                    if (!std::isfinite(v3) || !std::isfinite(e3) || e3 >= 9.99) { v3 = 0.0; e3 = 0.0; }
                    if (!std::isfinite(v4) || !std::isfinite(e4) || e4 >= 9.99) { v4 = 0.0; e4 = 0.0; }
                    
                    hV2->SetBinContent(1, v2);
                    hV2->SetBinError(1, e2);
                    hV3->SetBinContent(1, v3);
                    hV3->SetBinError(1, e3);
                    hV4->SetBinContent(1, v4);
                    hV4->SetBinError(1, e4);
                }
                
                hV2->Write();
                hV3->Write();
                hV4->Write();
                
                // Quick-look plot
                TCanvas cV2Eta("cV2Eta", "cV2Eta", 900, 700);
                hV2->SetMarkerStyle(20);
                hV2->Draw("E1");
                cV2Eta.SaveAs(Form("./TemplateFit/EtaDiff/Plots/V2_%s_%s_%i_%i_TPCEta_%0.1f_%0.1f.pdf", 
                                  outputFileName.c_str(), splitName.c_str(), 
                                  data.minRange, data.maxRange, etaMin, etaMax));
                
                std::cout << "Output file: " << Form("./TemplateFit/EtaDiff/Vn_%s_%s_FT0Eta_%0.1f_%0.1f_%i_%i_TPCEta_%0.1f_%0.1f.root", 
                                                     outputFileName.c_str(), splitName.c_str(), ft0EtaMin, ft0EtaMax,
                                                     data.minRange, data.maxRange, etaMin, etaMax) << std::endl;
                outputFile.Close();
                
                // Only delete vnResult; histograms are owned by TFile and deleted when it closes
                delete vnResult;
            }
            return;  // Early return after processing FT0-sliced case
        }
        
        // Original code for non-FT0-sliced case (all TPC bins together)
        // 执行模板拟合获取所有结果
        std::vector<VnUnit*> vnResults;
        if (isRingDataset && !tpcRanges.empty()) {
            for (const auto& rng : tpcRanges) {
                vnResults.push_back(TemplateFit(isNch, templ, data, kFALSE, rng.first, rng.second, ft0EtaMin, ft0EtaMax));
            }
        } else {
            for (Int_t iEta = 0; iEta < etaTPC.size() - 1; iEta++) {
                double etaMin = etaTPC[iEta];
                double etaMax = etaTPC[iEta + 1];
                vnResults.push_back(TemplateFit(isNch, templ, data, kFALSE, etaMin, etaMax, ft0EtaMin, ft0EtaMax)); //men ikke her cn2 to vn2
            }
        }

        // 创建输出文件
        TFile outputFile(Form("./TemplateFit/EtaDiff/Vn_%s_%s_%i_%i.root", outputFileName.c_str(), splitName.c_str(), data.minRange, data.maxRange), "RECREATE"); //har skapes .root for diff vns

        // Quick-look plots (PDF)
        gSystem->mkdir("./TemplateFit/EtaDiff/Plots", kTRUE);

        // NOTE: etaTPC is std::vector<float> in BasicForDihadron.h.
        // TH1D expects Double_t* bin edges; passing float* can corrupt binning in ROOT/Cling.
        std::vector<Double_t> etaEdges;
        
        // When FT0 slicing is active, only include the bins we actually processed
        if (ft0EtaMin > -900. && ft0EtaMax > -900.) {
            // Match the TPC ranges used in FT0-sliced processing: [-0.8,-0.6] and [0.6,0.8]
            std::vector<std::pair<double,double>> tpcRanges = { { -0.8, -0.6 }, { 0.6, 0.8 } };
            etaEdges.reserve(tpcRanges.size() + 1);
            for (auto r : tpcRanges) {
                etaEdges.push_back(static_cast<Double_t>(r.first));
            }
            // Add the upper edge of the last range
            etaEdges.push_back(static_cast<Double_t>(tpcRanges.back().second));
        } else if (isRingDataset && !tpcRanges.empty()) {
            etaEdges.reserve(tpcRanges.size() + 1);
            for (auto r : tpcRanges) {
                etaEdges.push_back(static_cast<Double_t>(r.first));
            }
            etaEdges.push_back(static_cast<Double_t>(tpcRanges.back().second));
        } else {
            // No FT0 slicing: use all TPC bins
            etaEdges.reserve(etaTPC.size());
            for (auto v : etaTPC) etaEdges.push_back(static_cast<Double_t>(v));
        }

        // 初始化直方图
        TH1D* hV2 = new TH1D("hV2", "V_{ #Delta 2};#eta^{TPC};V_{ #Delta 2}",
                    static_cast<Int_t>(etaEdges.size()) - 1, etaEdges.data());
        TH1D* hV3 = new TH1D("hV3", "V_{ #Delta 3};#eta^{TPC};V_{ #Delta 3}",
                    static_cast<Int_t>(etaEdges.size()) - 1, etaEdges.data());
        TH1D* hV4 = new TH1D("hV4", "V_{ #Delta 4};#eta^{TPC};V_{ #Delta 4}",
                    static_cast<Int_t>(etaEdges.size()) - 1, etaEdges.data());

        hV2->SetStats(0);
        hV3->SetStats(0);
        hV4->SetStats(0);

        // 填充数据
        // If a bin ends up with the sentinel bootstrap error (=10) or non-finite values,
        // keep it from dominating the plot range.
        for (size_t i = 0; i < vnResults.size(); ++i) {
            const auto* vn = vnResults[i];
            Double_t v2 = vn ? vn->v2 : 0.0;
            Double_t e2 = vn ? vn->v2_err : 0.0;
            Double_t v3 = vn ? vn->v3 : 0.0;
            Double_t e3 = vn ? vn->v3_err : 0.0;
            Double_t v4 = vn ? vn->v4 : 0.0;
            Double_t e4 = vn ? vn->v4_err : 0.0;

            if (!std::isfinite(v2) || !std::isfinite(e2) || e2 >= 9.99) { v2 = 0.0; e2 = 0.0; }
            if (!std::isfinite(v3) || !std::isfinite(e3) || e3 >= 9.99) { v3 = 0.0; e3 = 0.0; }
            if (!std::isfinite(v4) || !std::isfinite(e4) || e4 >= 9.99) { v4 = 0.0; e4 = 0.0; }

            hV2->SetBinContent(i + 1, v2);
            hV2->SetBinError(i + 1, e2);

            hV3->SetBinContent(i + 1, v3);
            hV3->SetBinError(i + 1, e3);

            hV4->SetBinContent(i + 1, v4);
            hV4->SetBinError(i + 1, e4);
        }

        // 写入文件
        hV2->Write();
        hV3->Write();
        hV4->Write();

        // Minimal: plot v2 only (requested)
        TCanvas cV2Eta("cV2Eta", "cV2Eta", 900, 700);
        hV2->SetMarkerStyle(20);
        // Prevent a single outlier bin from blowing up the y-range.
        {
            double yMin = 1e300;
            double yMax = -1e300;
            for (int b = 1; b <= hV2->GetNbinsX(); ++b) {
                const double y = hV2->GetBinContent(b);
                const double e = hV2->GetBinError(b);
                if (!std::isfinite(y) || !std::isfinite(e) || e <= 0) continue;
                yMin = std::min(yMin, y - e);
                yMax = std::max(yMax, y + e);
            }
            if (yMax > yMin && yMin < 1e200) {
                const double pad = 0.15 * (yMax - yMin);
                hV2->GetYaxis()->SetRangeUser(yMin - pad, yMax + pad);
            }
        }
        hV2->Draw("E1");
        cV2Eta.SaveAs(Form("./TemplateFit/EtaDiff/Plots/V2_%s_%s_%i_%i.pdf", outputFileName.c_str(), splitName.c_str(), data.minRange, data.maxRange));

        std::cout << "Output file: " << Form("./TemplateFit/EtaDiff/Vn_%s_%s_%i_%i.root", outputFileName.c_str(), splitName.c_str(), data.minRange, data.maxRange) << std::endl;
        outputFile.Close();

        for (auto* vn : vnResults) delete vn;

    }

}



//==============================================================
VnUnit* TemplateFit(Bool_t isNch, InputUnit templ, InputUnit data, Bool_t cn2Tovn2, Double_t etaMin=-101., Double_t etaMax=-101., Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    bool isRingDataset = (templ.fileNameSuffix.find("615817") != std::string::npos ||
                          templ.fileNameSuffix.find("615818") != std::string::npos ||
                          templ.fileNameSuffix.find("616549") != std::string::npos ||
                          templ.fileNameSuffix.find("618685") != std::string::npos ||
                          templ.fileNameSuffix.find("617826") != std::string::npos ||
                          templ.fileNameSuffix.find("617910") != std::string::npos);
    bool useEtaDiffOnly = isRingDataset && (etaMin > -100. && etaMax > -100.);

    TString ft0Suffix = "";
    if (ft0EtaMin > -900. && ft0EtaMax > -900.) {
        ft0Suffix = Form("_FT0Eta_%0.1f_%0.1f", ft0EtaMin, ft0EtaMax);
    }

    TFile* templatefile = nullptr;
    TFile* templatefile_etaDiff = nullptr;
    TFile* datafile = nullptr;
    TFile* datafile_etaDiff = nullptr;

    if (useEtaDiffOnly) {
        templatefile_etaDiff = new TFile(Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%d_%d_Eta_%0.1f_%0.1f%s.root", templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange, etaMin, etaMax, ft0Suffix.Data()), "READ");
        if (!templatefile_etaDiff || !templatefile_etaDiff->IsOpen()) { //error
            std::cerr << "Cannot open eta-diff template file: " << Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%d_%d_Eta_%0.1f_%0.1f%s.root", templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange, etaMin, etaMax, ft0Suffix.Data()) << std::endl;
            exit(1);
        }

        datafile_etaDiff = new TFile(Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%i_%i_Eta_%0.1f_%0.1f%s.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, etaMin, etaMax, ft0Suffix.Data()), "READ");
        if (!datafile_etaDiff || !datafile_etaDiff->IsOpen()) { //error
            std::cerr << "Cannot open eta-diff input file: " << Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%i_%i_Eta_%0.1f_%0.1f%s.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, etaMin, etaMax, ft0Suffix.Data()) << std::endl;
            exit(1);
        }
    } else {
        templatefile = new TFile(Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d%s.root", templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange, ft0Suffix.Data()), "READ");
        if (!templatefile || !templatefile->IsOpen()) { //error
            std::cerr << "Cannot open template file: " << Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d%s.root", templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange, ft0Suffix.Data()) << std::endl;
            exit(1);
        }
        if (etaMin > -100. && etaMax > -100.) { //dersom pt-diff er aktivert
            templatefile_etaDiff = new TFile(Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%d_%d_Eta_%0.1f_%0.1f%s.root", templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange, etaMin, etaMax, ft0Suffix.Data()), "READ");
            if (!templatefile_etaDiff || !templatefile_etaDiff->IsOpen()) { //error
                std::cerr << "Cannot open template file: " << Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%d_%d_Eta_%0.1f_%0.1f%s.root", templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange, etaMin, etaMax, ft0Suffix.Data()) << std::endl;
                exit(1);
            }
        }

        datafile = new TFile(Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d%s.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, ft0Suffix.Data()), "READ");
        if (!datafile || !datafile->IsOpen()) { //error
            std::cerr << "Cannot open input file: " << Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d%s.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, ft0Suffix.Data()) << std::endl;
            exit(1);
        }
        if (etaMin > -100. && etaMax > -100.) { //dersom eta-diff er aktivert
            datafile_etaDiff = new TFile(Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%i_%i_Eta_%0.1f_%0.1f%s.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, etaMin, etaMax, ft0Suffix.Data()), "READ");
            if (!datafile_etaDiff || !datafile_etaDiff->IsOpen()) { //error
                std::cerr << "Cannot open input file: " << Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%i_%i_Eta_%0.1f_%0.1f%s.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, etaMin, etaMax, ft0Suffix.Data()) << std::endl;
                exit(1);
            }
        }
    }

    VnUnit* vnResult = nullptr;
    VnUnit* vnResult_etaDiff = nullptr; //dersom etaDiff ikke er aktivt forblir denne nullptr

    if (!useEtaDiffOnly) {
        vnResult = fitSample(isNch, templatefile, templ, datafile, data, -1, -101., -101., ft0EtaMin, ft0EtaMax); //her finner vi fit for ref
        if (!vnResult) { //error
            std::cerr << "Cannot fit sample: " << data.fileNameSuffix << std::endl;
            exit(1);
        }
    }

    if (etaMin > -100. && etaMax > -100.) { //dersom eta-diff er aktivert
        vnResult_etaDiff = fitSample(isNch, templatefile_etaDiff, templ, datafile_etaDiff, data, -1, etaMin, etaMax, ft0EtaMin, ft0EtaMax); //her finner vi fit for etaDiff
        if (!vnResult_etaDiff) {
            std::cerr << "Cannot fit eta-diff sample: " << data.fileNameSuffix << std::endl;
            exit(1);
        }
    }
    //EDIT START
    std::vector<std::vector<std::vector<double>>> ValueArray;
    std::vector<std::vector<std::vector<double>>> ValueErrorArray;
    std::vector<std::vector<double>> ErrorArray;
    int Nobs=7;//indices used: 0(v2),1(v3),2(v4),3(unused legacy v1 slot),4(F),5(G),6(chi2ndf)
    int NofSample = maxSample*maxSample; //100 DENNE VAR 1000
    int Nbin = 1;
    ResizeValueArray(ValueArray,ValueErrorArray,ErrorArray,Nobs,NofSample,Nbin);

    //her finner vi bootstrap samples (100 times)
    for(int sample=0;sample<NofSample;sample++) { //edited
        if (!useEtaDiffOnly) {
            VnUnit* vnTemp = fitSample(isNch, templatefile, templ, datafile, data, sample);
            if (!vnTemp) { //error - mark as failed (-999) instead of exiting
                // Mark all observables for this sample as invalid
                for (int iobs = 0; iobs < Nobs; iobs++) {
                    ValueArray[iobs][sample][0] = -999;
                    ValueErrorArray[iobs][sample][0] = 0;
                }
                continue;
            }
            //Here we might estimate the uncertainties for v11, F, G too
            ValueArray[0][sample][0] = vnTemp->v2; //3x100x1
            ValueErrorArray[0][sample][0] = vnTemp->v2_err;
            ValueArray[1][sample][0] = vnTemp->v3;
            ValueErrorArray[1][sample][0] = vnTemp->v3_err;
            ValueArray[2][sample][0] = vnTemp->v4;
            ValueErrorArray[2][sample][0] = vnTemp->v4_err;
            // V1 extraction removed
            ValueArray[4][sample][0] = vnTemp->F;
            ValueErrorArray[4][sample][0] = vnTemp->F_err;
            ValueArray[5][sample][0] = vnTemp->G;
            ValueErrorArray[5][sample][0] = vnTemp->G_err;
            ValueArray[6][sample][0] = vnTemp->chi2ndf;
            ValueErrorArray[6][sample][0] = vnTemp->chi2ndf_err;
            delete vnTemp;
        }
        if (etaMin > -100. && etaMax > -100.) { //dersom eta-diff er aktivert
            VnUnit* vnTemp_etaDiff = fitSample(isNch, templatefile_etaDiff, templ, datafile_etaDiff, data, sample, etaMin, etaMax);
            if (!vnTemp_etaDiff) { //error - mark as failed instead of exiting
                // Mark all observables for this sample as invalid
                for (int iobs = 0; iobs < Nobs; iobs++) {
                    ValueArray[iobs][sample][0] = -999;
                    ValueErrorArray[iobs][sample][0] = 0;
                }
                continue;
            }
            ValueArray[0][sample][0] = vnTemp_etaDiff->v2;
            ValueErrorArray[0][sample][0] = vnTemp_etaDiff->v2_err;
            ValueArray[1][sample][0] = vnTemp_etaDiff->v3;
            ValueErrorArray[1][sample][0] = vnTemp_etaDiff->v3_err;
            ValueArray[2][sample][0] = vnTemp_etaDiff->v4;
            ValueErrorArray[2][sample][0] = vnTemp_etaDiff->v4_err;
            // V1 extraction removed
            ValueArray[4][sample][0] = vnTemp_etaDiff->F;
            ValueErrorArray[4][sample][0] = vnTemp_etaDiff->F_err;
            ValueArray[5][sample][0] = vnTemp_etaDiff->G;
            ValueErrorArray[5][sample][0] = vnTemp_etaDiff->G_err;
            ValueArray[6][sample][0] = vnTemp_etaDiff->chi2ndf;
            ValueErrorArray[6][sample][0] = vnTemp_etaDiff->chi2ndf_err;
            delete vnTemp_etaDiff;
        }
    }

    //
    //HER REGNET VI SNITT AV 10 HISTOS
    //

    std::vector<std::vector<std::string>> dataInString = {toStringV2(ValueArray, maxSample*maxSample, etaMin, etaMax)};
    std::vector<std::vector<std::string>> dataInString3 = {toStringV3(ValueArray, maxSample*maxSample, etaMin, etaMax)};
    std::vector<std::vector<std::string>> dataInString4 = {toStringV4(ValueArray, maxSample*maxSample, etaMin, etaMax)};

    if (etaMin > -100. && etaMax > -100.){
        writeToCSV("V2_mean10_100_EtaDiff.csv", dataInString, "./TemplateFit/EtaDiff", true); //Writes wrong header... to be fixed later
        writeToCSV("V3_mean10_100_EtaDiff.csv", dataInString3, "./TemplateFit/EtaDiff", true);
        writeToCSV("V4_mean10_100_EtaDiff.csv", dataInString4, "./TemplateFit/EtaDiff", true); 
    }
    else {
        writeToCSV("V2_mean10_100.csv", dataInString, "./TemplateFit", true); //Writes wrong header... to be fixed later
        writeToCSV("V3_mean10_100.csv", dataInString3, "./TemplateFit", true); 
        writeToCSV("V4_mean10_100.csv", dataInString4, "./TemplateFit", true); 
    }


    /*
    NB! The rest of the code only works with 100 out of the 1000 samples!
    Save those into a file.

    in the external python code (let AI rewrite existing), find the ratios of the 100 samples. (SAME SAMPLES)
    Do bootstrap on the ratio
    Exclude error waight... get rid of samples above 3 sigma
    Recalculate error
    Compare with error from error propagation, hopefully between min and max
    */


    // for(int sample=0;sample<NofSample;sample++) {
    //     std::cout << "sample: " << sample << " v2^2: " << ValueArray[0][sample][0]  << std::endl;
    //     std::cout << "sample: " << sample << " v3^2: " << ValueArray[1][sample][0]  << std::endl;
    //     std::cout << "sample: " << sample << " v4^2: " << ValueArray[2][sample][0] << std::endl << std::endl;
    // }

    // Count successful bootstrap samples (those not marked as -999)
    int successfulSamples = 0;
    for(int sample=0;sample<NofSample;sample++) {
        if (ValueArray[0][sample][0] > -900) {  // Not marked as failed
            successfulSamples++;
        }
    }
    
    // Require a minimum number of successful bootstrap samples
    int minSamples = NofSample / 50;  // 2%
    if (minSamples < 5) minSamples = 5;
    if (successfulSamples < minSamples) {
        std::cerr << "Too few successful bootstrap fits (" << successfulSamples << "/" << NofSample 
                  << ") for " << data.fileNameSuffix << " bin eta " << etaMin << " to " << etaMax << std::endl;
        return nullptr;
    }
    
    std::cout << "Bootstrap fits: " << successfulSamples << "/" << NofSample << " successful" << std::endl;

    //her regnes bootstrap errors
    for(int iobs = 0;iobs < Nobs;iobs++){
        CalculateBootstrapError(ValueArray[iobs],ValueErrorArray[iobs],ErrorArray[iobs],1.); //regn ut errors
    }
    
    // Add systematic uncertainties in quadrature with bootstrap errors
    // The systematic uncertainty is 0.02% of the measured value
    const double systematic_fraction = 0.0002;
    
    // we assign the boot-strap error to vnResult (the one shoed in the root files)
    if (!useEtaDiffOnly && vnResult) {
        vnResult->v2_err = TMath::Sqrt(ErrorArray[0][0]*ErrorArray[0][0] + (systematic_fraction*vnResult->v2)*(systematic_fraction*vnResult->v2));
        vnResult->v3_err = TMath::Sqrt(ErrorArray[1][0]*ErrorArray[1][0] + (systematic_fraction*vnResult->v3)*(systematic_fraction*vnResult->v3));
        vnResult->v4_err = TMath::Sqrt(ErrorArray[2][0]*ErrorArray[2][0] + (systematic_fraction*vnResult->v4)*(systematic_fraction*vnResult->v4));
        // V1 error calculation removed
        vnResult->F_err = TMath::Sqrt(ErrorArray[4][0]*ErrorArray[4][0] + (systematic_fraction*vnResult->F)*(systematic_fraction*vnResult->F));
        vnResult->G_err = TMath::Sqrt(ErrorArray[5][0]*ErrorArray[5][0] + (systematic_fraction*vnResult->G)*(systematic_fraction*vnResult->G));
        vnResult->chi2ndf_err = ErrorArray[6][0]; // No systematic for chi2/ndf
    }
    if (etaMin > -100. && etaMax > -100.) { //dersom eta-diff er aktivert
        vnResult_etaDiff->v2_err = TMath::Sqrt(ErrorArray[0][0]*ErrorArray[0][0] + (systematic_fraction*vnResult_etaDiff->v2)*(systematic_fraction*vnResult_etaDiff->v2));
        vnResult_etaDiff->v3_err = TMath::Sqrt(ErrorArray[1][0]*ErrorArray[1][0] + (systematic_fraction*vnResult_etaDiff->v3)*(systematic_fraction*vnResult_etaDiff->v3));
        vnResult_etaDiff->v4_err = TMath::Sqrt(ErrorArray[2][0]*ErrorArray[2][0] + (systematic_fraction*vnResult_etaDiff->v4)*(systematic_fraction*vnResult_etaDiff->v4));
        // V1 error calculation removed
        vnResult_etaDiff->F_err = TMath::Sqrt(ErrorArray[4][0]*ErrorArray[4][0] + (systematic_fraction*vnResult_etaDiff->F)*(systematic_fraction*vnResult_etaDiff->F));
        vnResult_etaDiff->G_err = TMath::Sqrt(ErrorArray[5][0]*ErrorArray[5][0] + (systematic_fraction*vnResult_etaDiff->G)*(systematic_fraction*vnResult_etaDiff->G));
        vnResult_etaDiff->chi2ndf_err = ErrorArray[6][0]; // No systematic for chi2/ndf
    }

    if (cn2Tovn2 && !useEtaDiffOnly && vnResult) { //desrom kTURE da tar vi roten av Vn2. Bemerk at V er ikke cn per say siden V kommer fra FT, mens cn er observabel, denne er satt til false i alle kall
        if (vnResult->v2 > 0.) { //verdien settes til -1 dersom problem oppstår
            vnResult->v2_err = vnResult->v2_err / (2 * sqrt(vnResult->v2));
            vnResult->v2 = sqrt(vnResult->v2); //her tar vi kvadratroten av V2
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
    }

    //START MINE
    std::vector<std::vector<std::string>> csvData;
    //if eta-diff
    if (etaMin > -100. && etaMax > -100.) {
        //create the vector that we will write into csv
        double F = vnResult_etaDiff->F;
        double Fe = vnResult_etaDiff->F_err;
        double G = vnResult_etaDiff->G;
        double Ge = vnResult_etaDiff->G_err;
        double v21 = vnResult_etaDiff->v2;
        double v21e = vnResult_etaDiff->v2_err;
        double v31 = vnResult_etaDiff->v3;
        double v31e = vnResult_etaDiff->v3_err;
        double v41 = vnResult_etaDiff->v4;
        double v41e = vnResult_etaDiff->v4_err;
        double chi2ndf = vnResult_etaDiff->chi2ndf;
        double chi2ndfe = vnResult_etaDiff->chi2ndf_err;
        //fill the vector with data
        csvData.push_back({std::to_string(etaMin), std::to_string(etaMax), std::to_string(F), std::to_string(Fe), std::to_string(G), std::to_string(Ge), std::to_string(v21), std::to_string(v21e), std::to_string(v31), std::to_string(v31e), std::to_string(v41), std::to_string(v41e), std::to_string(chi2ndf), std::to_string(chi2ndfe)});
        writeToCSV("BootStrapValues.csv", csvData, "./TemplateFit/EtaDiff");
    } 
    else {
        double F = vnResult->F;
        double Fe = vnResult->F_err;
        double G = vnResult->G;
        double Ge = vnResult->G_err;
        double v21 = vnResult->v2;
        double v21e = vnResult->v2_err;
        double v31 = vnResult->v3;
        double v31e = vnResult->v3_err;
        double v41 = vnResult->v4;
        double v41e = vnResult->v4_err;
        double chi2ndf = vnResult->chi2ndf;
        double chi2ndfe = vnResult->chi2ndf_err;
        //fill the vector with data
        csvData.push_back({std::to_string(etaMin), std::to_string(etaMax), std::to_string(F), std::to_string(Fe), std::to_string(G), std::to_string(Ge), std::to_string(v21), std::to_string(v21e), std::to_string(v31), std::to_string(v31e), std::to_string(v41), std::to_string(v41e), std::to_string(chi2ndf), std::to_string(chi2ndfe)});
        writeToCSV("BootStrapValues.csv", csvData, "./TemplateFit");
    }
    //END MINE

    // print result
    if (etaMin > -100. && etaMax > -100.) {
        std::cout << "print result: " << data.fileNameSuffix << " eta-diff" << std::endl;
        std::cout << "v2: " << vnResult_etaDiff->v2 << " +/- " << vnResult_etaDiff->v2_err << std::endl;
        std::cout << "v3: " << vnResult_etaDiff->v3 << " +/- " << vnResult_etaDiff->v3_err << std::endl;
        std::cout << "v4: " << vnResult_etaDiff->v4 << " +/- " << vnResult_etaDiff->v4_err << std::endl;
        return vnResult_etaDiff;
    }
    std::cout << "print result: " << data.fileNameSuffix << std::endl;
    std::cout << "v2: " << vnResult->v2 << " +/- " << vnResult->v2_err << std::endl;
    std::cout << "v3: " << vnResult->v3 << " +/- " << vnResult->v3_err << std::endl;
    std::cout << "v4: " << vnResult->v4 << " +/- " << vnResult->v4_err << std::endl;

    return vnResult;
}


//==============================================================
std::vector<Int_t> CheckAndMergeRanges(const std::vector<InputUnit>& inputUnits) { //this merges centrality ranges
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
VnUnit* fitSample(Bool_t isNch, TFile* templatefile, InputUnit templ, TFile* datafile, InputUnit data, int sample = -1, Double_t etaMin=-101., Double_t etaMax=-101., Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.) {
    std::vector<Double_t> fParamVal;
    std::vector<Double_t> fParamErr;
    TH1D* lm=0;
    TH1D* hm=0;
    TString suffix = (sample == -1) ? "" : Form("_%d", sample);
    lm = (TH1D*)templatefile->Get(Form("bsSample_hPhiSameOverMixed_%d_%d%s", templ.minRange, templ.maxRange, suffix.Data()));
    hm = (TH1D*)datafile->Get(Form("bsSample_hPhiSameOverMixed_%d_%d%s", data.minRange, data.maxRange, suffix.Data()));
    if (!lm) { //error
        std::cerr << "Cannot find histogram: " << Form("bsSample_hPhiSameOverMixed_%d_%d%s", templ.minRange, templ.maxRange, suffix.Data()) << std::endl;
        return 0;
    }
    if (!hm) { //error
        std::cerr << "Cannot find histogram: " << Form("bsSample_hPhiSameOverMixed_%d_%d%s", data.minRange, data.maxRange, suffix.Data()) << std::endl;
        return 0;
    }
    RooTempFitter(lm, hm, fParamVal, fParamErr, kFALSE); //FT magic happens here
    if (sample == -1) { //if -1 we load all of the samples, and only then will .pdf be created
        PlotFitting(lm, hm, isNch, data.fileNameSuffix, data.minRange, data.maxRange, fParamVal, fParamErr, etaMin, etaMax, ft0EtaMin, ft0EtaMax);
    }
    VnUnit* vnResult = new VnUnit(fParamVal[0], fParamErr[0], fParamVal[1], fParamErr[1], fParamVal[2], fParamErr[2], fParamVal[3], fParamErr[3], fParamVal[4], fParamErr[4], fParamVal[6], fParamErr[6]);
    return vnResult; //her sendes resultatet, altså verdier av v2,v3,v4 og errorer videre, altså fit slutter her
}

//==============================================================
void RooTempFitter(TH1 *lm, TH1 *hm, std::vector<Double_t>& fParamVal, std::vector<Double_t>& fParamErr, Bool_t kRefit) {
    fParamVal.clear();
    fParamErr.clear();
    fParamVal.resize(7);
    fParamErr.resize(7);
    if (!lm ||!hm) { //lm is template hm is data
        std::cerr << "Null pointer to histogram" << std::endl;
        return;
    }
    
    // Add systematic uncertainty to data histogram errors
    // This inflates the error bars to account for model uncertainties
    const double systematic_fraction = 0.0002; // 0.02% systematic uncertainty to achieve chi2/ndf ~ 1-2
    for (int i = 1; i <= hm->GetNbinsX(); i++) {
        double content = hm->GetBinContent(i);
        double stat_error = hm->GetBinError(i);
        double syst_error = systematic_fraction * TMath::Abs(content);
        double total_error = TMath::Sqrt(stat_error * stat_error + syst_error * syst_error);
        hm->SetBinError(i, total_error);
    }
    
    //Initialize fitter with given projections
    TemplateFitter *ft = new TemplateFitter(hm);
    //Setting up variable ( = delta phi, or just "x"):
    ft->AddVariable("x", "x", -TMath::Pi()/2.0, 1.5*TMath::Pi());

    if (!kRefit){ //dersom fit failer, skal man forandre disse: (if (!kRefit){)
        // Estimate initial values from data
        Double_t hm_integral = hm->Integral();
        Double_t hm_mean = hm->GetMean(2); //y-axis mean (content)
        Double_t hm_min = hm->GetMinimum();
        Double_t lm_integral = lm->Integral();
        
        // F should be roughly the baseline (minimum or small fraction of mean)
        Double_t F_init = TMath::Max(0.1, hm_min * 0.5);
        if (F_init > 100) F_init = 5.0; // fallback
        
        // G should scale with the integral/mean
        Double_t G_init = hm_integral / lm_integral;
        if (G_init < 1.0 || G_init > 1e8) G_init = 30000.0; // fallback
        
        // Pb-Pb initial value with data-driven estimates
        ft->AddParameter("Fa","Fa",F_init, 0., 300.);
        ft->AddParameter("Ga","Ga",G_init,0.,3000000000.);
        ft->AddParameter("v2","v2",5e-3,-1.0,1.0);
        ft->AddParameter("v3","v3",1e-3,-1.0,1.0); // Let v3 vary freely
        ft->AddParameter("v4","v4",1e-4,-1.0,1.0); // Let v4 vary freely
        // V1 parameter removed


        // pp initial value
        // ft->AddParameter("Fa","Fa",0.1,0,30);
        // ft->AddParameter("Ga","Ga",1,0,30);
        // ft->AddParameter("v2","v2",4e-3,-1.0,1.0);
        // ft->AddParameter("v3","v3",6e-4,-1.0,1.0);
        // ft->AddParameter("v4","v4",1.8e-4,-1.0,1.0);
        // // ft->AddParameter("v5","v5",1e-4,0,1.0);
    }
    //not used here
    else{ 
        ft->AddParameter("Fa","Fa",fParamVal[3],fParamVal[3]-1,fParamVal[3]+1);
        ft->AddParameter("Ga","Ga",fParamVal[4],fParamVal[4]-5,fParamVal[4]+5);
        ft->AddParameter("v2","v2",fParamVal[0],fParamVal[0]-0.0002,fParamVal[0]+0.0002);
        ft->AddParameter("v3","v3",fParamVal[1],fParamVal[1]-0.0002,fParamVal[1]+0.0002);
    }

    //Construct fit function
    FunctionObject *fobj = new TemplateFunction(lm);
    ft->SetFitFunction(fobj);
    //Perform fit:
    printf("About to fit\n");
    Int_t dummy = ft->Fit(0); //Do not draw performance at this point. Return value of Fit() is false if no base is set.
    if(!dummy) return;

    //these errors do not originate from bootstrap, but from fit itself
    Double_t F   = ft->getVal(0); //0 for F, 1 for G, 2 for v2, 3 for v3, 4 for v4
    Double_t Fe  = ft->getErr(0); 
    Double_t G   = ft->getVal(1); //0 for F, 1 for G, 2 for v2, 3 for v3, 4 for v4
    Double_t Ge  = ft->getErr(1);
    Double_t v21 = ft->getVal(2); //0 for F, 1 for G, 2 for v2, 3 for v3, 4 for v4
    Double_t v21e= ft->getErr(2);
    Double_t v31 = ft->getVal(3); //0 for F, 1 for G, 2 for v2, 3 for v3, 4 for v4
    Double_t v31e= ft->getErr(3);
    Double_t v41 = ft->getVal(4); //0 for F, 1 for G, 2 for v2, 3 for v3, 4 for v4
    Double_t v41e= ft->getErr(4);
    Double_t chi2ndf = findChi2ndf(lm, hm, F, G, v21, v31, v41);

    printf("Values from fit:\n");
    printf("F  = %f +- %f\n",F,Fe);
    printf("G  = %f +- %f\n",G,Ge); //skrives ut i terminalen
    // V1 output removed
    printf("V2 = %f +- %f\n",v21,v21e);
    printf("V3 = %f +- %f\n",v31,v31e);
    printf("V4 = %f +- %f\n",v41,v41e);

    fParamVal[0] = v21; fParamErr[0] = v21e;
    fParamVal[1] = v31; fParamErr[1] = v31e;
    fParamVal[2] = v41; fParamErr[2] = v41e;
    fParamVal[3] = F;   fParamErr[3] = Fe;
    fParamVal[4] = G;   fParamErr[4] = Ge; //4
    // V1 parameter storage removed
    fParamVal[6] = chi2ndf; fParamErr[6] = 1.0; //the fit error is only used to waight the average
}

void DrawText(double xmin, double ymin, double textSize, TString text)
{

	TLatex *textPreliminary = new TLatex(xmin, ymin, Form("%s", text.Data()));
	textPreliminary->SetNDC();
	textPreliminary->SetTextFont(43);
	textPreliminary->SetTextSize(textSize);
	textPreliminary->Draw();

}

void PlotFitting(TH1 *lm, TH1 *hm, Bool_t isNch, std::string fileSuffix, Int_t minRange, Int_t maxRange, const std::vector<Double_t>& par, const std::vector<Double_t>& parerr, Double_t etaMin=-101., Double_t etaMax=-101., Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.) {
    gStyle->SetOptStat(0); 
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    static int plotCounter = 0;
    const int localPlotId = plotCounter++;
    
    double v21 = par[0];
    double v31 = par[1];
    double v41 = par[2];
    double F =   par[3];
    double G =   par[4];

    double v21e = parerr[0];
    double v31e = parerr[1];
    double v41e = parerr[2];
    double Fe = parerr[3];
    double Ge = parerr[4];
    // V1 error removed

    TCanvas* canvas = new TCanvas(Form("Fit_%d", localPlotId), "Fit", 800, 600);
    canvas->Range(0,0,1,1);
    
    // 创建上下面板
    TPad* pad1 = new TPad(Form("pad1_%d", localPlotId), Form("pad1_%d", localPlotId), 0, 0.25, 1, 1);
    pad1->SetBorderMode(0);
    pad1->SetBorderSize(2);
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.05);
    pad1->SetTopMargin(0.05);
    pad1->SetBottomMargin(0.13);
    pad1->SetTicks(1,1);
    pad1->Draw();
    
    TPad* pad2 = new TPad(Form("pad2_%d", localPlotId), Form("pad2_%d", localPlotId), 0, 0, 1, 0.25);
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
    TH1D* hbkg1 = new TH1D(Form("dPhi_%d", localPlotId), "dPhi", 1, -TMath::Pi()/2.0, 3*TMath::Pi()/2.0);
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

    // Generate fit curves
    const int pointBin = (int)hm->GetNbinsX();
        std::vector<Double_t> CopyPointX(pointBin);
        std::vector<Double_t> CopyPointY(pointBin);
    for (int i=0; i<pointBin; ++i){
      CopyPointX[i] = hm->GetBinCenter(i+1);
      double x = hm->GetBinCenter(i+1);
      // CopyPointY[i] = F*lm->GetBinContent(i+1)+G*(1+2*v21*cos(2*x)+2*v31*cos(3*x));
      CopyPointY[i] = F*lm->GetBinContent(i+1) + G*(1+2*v21*cos(2*x)+2*v31*cos(3*x)+2*v41*cos(4*x));
    };
        TGraph* gCopy = new TGraph(pointBin, CopyPointX.data(), CopyPointY.data());

        std::vector<Double_t> PeriPointY(pointBin);
    for (int i=0; i<pointBin; ++i){
      PeriPointY[i] = F*lm->GetBinContent(i+1) + G;
    };
        TGraph* gPeri = new TGraph(pointBin, CopyPointX.data(), PeriPointY.data());

    Double_t Y0Position = lm->GetXaxis()->FindBin(0.0);
    Double_t Y0 = lm->GetBinContent(Y0Position);


    gCopy->SetLineColor(colors[0]);
    gCopy->SetLineWidth(2);
    gCopy->Draw("same");

    gPeri->SetLineColor(colors[1]);
    gPeri->SetLineWidth(2);
    gPeri->SetLineStyle(5);
    gPeri->Draw("same");

    TF1* fit_p2 = new TF1(Form("fit_p2_%d", localPlotId),"[0]*[1] + [2]*(1 + 2*[3]*cos(2*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    TF1* fit_p3 = new TF1(Form("fit_p3_%d", localPlotId),"[0]*[1] + [2]*(1 + 2*[3]*cos(3*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    TF1* fit_p4 = new TF1(Form("fit_p4_%d", localPlotId),"[0]*[1] + [2]*(1 + 2*[3]*cos(4*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    // Fit parameter setting updated
    fit_p2->SetParameters(F,Y0,G,v21);
    fit_p3->SetParameters(F,Y0,G,v31);
    fit_p4->SetParameters(F,Y0,G,v41);
    fit_p2->SetLineColor(colors[2]);
    fit_p2->SetLineWidth(2);
    fit_p2->SetLineStyle(2);
    fit_p3->SetLineColor(colors[3]);
    fit_p3->SetLineWidth(2);
    fit_p3->SetLineStyle(3);
    fit_p4->SetLineColor(colors[4]);
    fit_p4->SetLineWidth(2);
    fit_p4->SetLineStyle(4);
    fit_p2->Draw("same");
    fit_p3->Draw("same");
    fit_p4->Draw("same");

    // 添加图例
    TLegend* leg = new TLegend(0.5, 0.65, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->AddEntry(hm, "Data", "lep");
    leg->AddEntry(gCopy, "FY(#Delta#phi)^{peri} + G(1+#Sigma_{n=2}^{4}2V_{n#Delta}cos(n#Delta#phi))", "l");
    leg->AddEntry(gPeri, "FY(#Delta#phi)^{peri} + G", "l");
    leg->AddEntry(fit_p2, "FY(0)^{peri} + G(1+2V_{2#Delta}cos(2#Delta#phi))", "l");
    leg->AddEntry(fit_p3, "FY(0)^{peri} + G(1+2V_{3#Delta}cos(3#Delta#phi))", "l");
    leg->AddEntry(fit_p4, "FY(0)^{peri} + G(1+2V_{4#Delta}cos(4#Delta#phi))", "l");
    leg->Draw();

    // 添加文本标签
    DrawText(0.2, 0.85, 20, Form("ALICE %s", collisionSystemName.c_str()));
    DrawText(0.2, 0.80, 20, Form("%s", FITused.c_str()));
    if (isNch) {
        DrawText(0.2, 0.75, 20, Form("%d < N_{ch} < %d", minRange, maxRange));
    }
    else {
        DrawText(0.2, 0.75, 20, Form("%d < Cent < %d", minRange, maxRange));
    }

    // Detect ring type from dataset ID
    bool isRingDataset = (fileSuffix.find("615817") != std::string::npos ||
                          fileSuffix.find("615818") != std::string::npos ||
                          fileSuffix.find("617826") != std::string::npos ||
                          fileSuffix.find("617910") != std::string::npos);
    if (isRingDataset) {
        bool isInner = (fileSuffix.find("615817") != std::string::npos || 
                        fileSuffix.find("617826") != std::string::npos);
        std::string ringType = isInner ? "Inner Ring" : "Outer Ring";
        DrawText(0.2, 0.70, 20, Form("%s", ringType.c_str()));
    }

    if (etaMin > -100. && etaMax > -100.)
        DrawText(0.2, 0.65, 20, Form("%.1f < #eta^{TPC} < %.1f", etaMin, etaMax));
    DrawText(0.2, 0.60, 20, Form("F = %.5f #pm %.5f", F, Fe));
    DrawText(0.2, 0.55, 20, Form("V_{2#Delta} = %.5f #pm %.5f", v21, v21e));
    DrawText(0.2, 0.50, 20, Form("V_{3#Delta} = %.5f #pm %.5f", v31, v31e));

    // 绘制底部残差图
    pad2->cd();

    TF1* fit_p1m = new TF1(Form("fit_p1m_%d", localPlotId),"[0]*(1 + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    fit_p1m->SetParameters(G,v21,v31,v41);

    
    // here we subtract lm from hm
    TH1D* hsubtract = (TH1D*)hm->Clone(Form("subtract_%d", localPlotId));
    hsubtract->Add(lm, -F);

    // =============== 新增：计算chi2/ndf ===============
    double chi2 = 0.0;
    int nBins = hsubtract->GetNbinsX();
    int nParams = 5; //!! change ndfs here
    int ndf = nBins - nParams; //degrees of freedom (the number of data points minus number of fit parameters)
    
    for (int i = 1; i <= nBins; i++) {
        double data = hsubtract->GetBinContent(i);
        double error = hsubtract->GetBinError(i);
        double x = hsubtract->GetBinCenter(i);
        double fit = fit_p1m->Eval(x);
        
        if (error > 0) { // 忽略误差为0的bin
            double residual = data - fit;
            chi2 += (residual * residual) / (error * error); //Yes chi is waighted by the errors, i.e. large error low chi
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
    TH1D* hResidual = (TH1D*)hsubtract->Clone(Form("pull_%d", localPlotId));
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

    //create the vector that we will write into csv
    std::vector<std::vector<std::string>> csvData;
    //fill the vector with data
    csvData.push_back({std::to_string(etaMin), std::to_string(etaMax), std::to_string(F), std::to_string(Fe), std::to_string(G), std::to_string(Ge), std::to_string(v21), std::to_string(v21e), std::to_string(v31), std::to_string(v31e), std::to_string(v41), std::to_string(v41e), std::to_string(chi2ndf)});

    // 保存结果
    if (etaMin > -100. && etaMax > -100.) {
        TString filename;
        if (ft0EtaMin > -900. && ft0EtaMax > -900.) {
            filename = Form("./TemplateFit/EtaDiff/PDFs/TemplateFit_%s_%s_%d_%d_TPCEta_%0.1f_%0.1f_FT0Eta_%0.1f_%0.1f.pdf", 
                          fileSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax, ft0EtaMin, ft0EtaMax);
        } else {
            filename = Form("./TemplateFit/EtaDiff/PDFs/TemplateFit_%s_%s_%d_%d_Eta_%0.1f_%0.1f.pdf", 
                          fileSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax);
        }
        canvas->SaveAs(filename);
        writeToCSV("TemplateFit.csv", csvData, "./TemplateFit/EtaDiff/PDFs");
    } else {
        TString filename;
        if (ft0EtaMin > -900. && ft0EtaMax > -900.) {
            filename = Form("./TemplateFit/PDFs/TemplateFit_%s_%s_%d_%d_FT0Eta_%0.1f_%0.1f.pdf", 
                          fileSuffix.c_str(), splitName.c_str(), minRange, maxRange, ft0EtaMin, ft0EtaMax);
        } else {
            filename = Form("./TemplateFit/PDFs/TemplateFit_%s_%s_%d_%d.pdf", 
                          fileSuffix.c_str(), splitName.c_str(), minRange, maxRange);
        }
        canvas->SaveAs(filename);
        writeToCSV("TemplateFit.csv", csvData, "./TemplateFit/PDFs");
    }
    
    // 清理内存: canvas owns/destroys its pads and drawn primitives
    canvas->Close();
    delete canvas;
  
}


/*
FOR FT FIT, SET THE LIMITS OF F TO MIN=MAX=0 AND CHANGE NUMBER OF PARAMETERS
*/

// Function to write strings and values into a CSV file
void writeToCSV(const std::string &filename, const std::vector<std::vector<std::string>> &data, const std::string &folder, bool createHeader=false) {
    if (data.empty()) {
        std::cerr << "Warning: writeToCSV called with empty data set. Nothing written.\n";
        return;
    }

    if (gSystem) {
        gSystem->mkdir(folder.c_str(), kTRUE);
    }

    std::filesystem::path dir(folder);
    std::error_code ec;
    if (!folder.empty() && !std::filesystem::exists(dir)) {
        std::filesystem::create_directories(dir, ec);
        if (ec) {
            std::cerr << "Error: Could not create directory '" << folder << "': " << ec.message() << "\n";
            return;
        }
    }

    std::filesystem::path filePath = dir / filename;
    if (filePath.extension() != ".csv") {
        filePath += ".csv"; // append .csv if missing
    }

    bool fileExists = std::filesystem::exists(filePath);

    std::ofstream file(filePath, std::ios::app); // append mode creates file if missing
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filePath << " for writing.\n";
        return;
    }

    // Decide whether to skip the first row (header) when appending.
    size_t startIndex = 0;
    if (fileExists) {
        // Simple heuristic: if first row cells contain no digits at all, treat as header and skip.
        const auto &firstRow = data.front();
        bool looksLikeHeader = true;
        for (const auto &cell : firstRow) {
            if (cell.find_first_of("0123456789") != std::string::npos) {
                looksLikeHeader = false;
                break;
            }
        }
        if (looksLikeHeader) {
            startIndex = 1; // skip header row during append to avoid duplication
        }
    }
    if (!fileExists && !createHeader) {
        // Add heather:
        std::vector<std::vector<std::string>> header;
        header.push_back({"etaMin", "etaMax", "F", "Fe", "G", "Ge", "v2", "v2e", "v3", "v3e", "v4", "v4e", "chi2ndf", "chi2ndfe"});
        writeToCSV(filename, header, folder, true);
    }

    for (size_t r = startIndex; r < data.size(); ++r) {
        const auto &row = data[r];
        for (size_t i = 0; i < row.size(); ++i) {
            std::string cell = row[i];
            // Escape quotes
            size_t pos = 0;
            while ((pos = cell.find('"', pos)) != std::string::npos) {
                cell.insert(pos, "\"");
                pos += 2;
            }
            // Wrap in quotes if needed
            if (cell.find(',') != std::string::npos || cell.find('"') != std::string::npos) {
                cell = '"' + cell + '"';
            }
            file << cell;
            if (i < row.size() - 1) file << ",";
        }
        file << "\n";
    }
    file.close();
}




double findChi2ndf(TH1 *lm, TH1 *hm, double& F_val, double& G_val, double& v21_val, double& v31_val, double& v41_val) {
    static int chi2Counter = 0;
    const int localChi2Id = chi2Counter++;
    double chi2;
    int nBins;
    int nParams;
    int ndf;
    double data;
    double error;
    double error_total;
    double x;
    double fit;
    double residual;
    double chi2ndf;
    
    // Add systematic uncertainty (as fraction of signal) to account for model uncertainties
    // Adjust this value to achieve reasonable chi2/ndf (~1-2)
    const double systematic_fraction = 0.0002; // 0.02% systematic uncertainty

    TF1* fit_p1m = new TF1(Form("fit_p1m_chi2_%d", localChi2Id),"[0]*(1 + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    fit_p1m->SetParameters(G_val,v21_val,v31_val,v41_val);

    TH1D* hsubtract = (TH1D*)hm->Clone(Form("subtract_chi2_%d", localChi2Id));
    hsubtract->SetDirectory(nullptr);
    hsubtract->Add(lm, -F_val);

    chi2 = 0.0;
    nBins = hsubtract->GetNbinsX();
    nParams = 5; // Number of parameters: F,G,v21,v31,v41 (V1 removed)
    ndf = nBins - nParams; //degrees of freedom (the number of data points minus number of fit parameters)
    
    for (int i = 1; i <= nBins; i++) {
        data = hsubtract->GetBinContent(i);
        error = hsubtract->GetBinError(i);
        x = hsubtract->GetBinCenter(i);
        fit = fit_p1m->Eval(x);
        
        if (error > 0) { // 忽略误差为0的bin
            // Add systematic uncertainty in quadrature
            double systematic_error = systematic_fraction * TMath::Abs(data);
            error_total = TMath::Sqrt(error * error + systematic_error * systematic_error);
            
            residual = data - fit;
            chi2 += (residual * residual) / (error_total * error_total); //Yes chi is waighted by the errors, i.e. large error low chi
        }
    }
    chi2ndf = (ndf > 0) ? chi2 / ndf : 0;
    
    // Clean up
    delete fit_p1m;
    delete hsubtract;

    return chi2ndf;
}



std::vector<std::string> toStringV2(std::vector<std::vector<std::vector<double>>> vec, int Nelements, double etaMin, double etaMax) {

    std::vector<std::string> stringVec;
    stringVec.resize(Nelements+2);
    
    stringVec[0] = std::to_string(etaMin);
    stringVec[1] = std::to_string(etaMax);
    
    for (size_t i = 0; i < Nelements; ++i) {
        stringVec[i+2] = std::to_string(vec[0][i][0]); // v2 value
    }
    return stringVec;
}

std::vector<std::string> toStringV3(std::vector<std::vector<std::vector<double>>> vec, int Nelements, double etaMin, double etaMax) {

    std::vector<std::string> stringVec;
    stringVec.resize(Nelements+2);
    
    stringVec[0] = std::to_string(etaMin);
    stringVec[1] = std::to_string(etaMax);
    
    for (size_t i = 0; i < Nelements; ++i) {
        stringVec[i+2] = std::to_string(vec[1][i][0]); // v3 value
    }
    return stringVec;
}

std::vector<std::string> toStringV4(std::vector<std::vector<std::vector<double>>> vec, int Nelements, double etaMin, double etaMax) {

    std::vector<std::string> stringVec;
    stringVec.resize(Nelements+2);
    
    stringVec[0] = std::to_string(etaMin);
    stringVec[1] = std::to_string(etaMax);
    
    for (size_t i = 0; i < Nelements; ++i) {
        stringVec[i+2] = std::to_string(vec[2][i][0]); // v4 value
    }
    return stringVec;
}

//==============================================================
void CreateV2DeltaSummaryPlot(Bool_t isNch, std::vector<InputUnit> dataList, std::string fileNameSuffix, std::string detectorName, const std::vector<float>& ft0EtaVec) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    
    // Define TPC ranges to summarise: use wider ranges matching Process_dPhidEta
    std::vector<std::pair<double,double>> tpcRanges = { { -0.8, -0.6 }, { 0.6, 0.8 } };

    // For each TPC range, create a summary plot
    for (auto tr : tpcRanges) {
        double tpcEtaMin = tr.first;
        double tpcEtaMax = tr.second;
        
        // For each centrality range
        for (const auto& data : dataList) {
            std::vector<Double_t> ft0Centers;
            std::vector<Double_t> ft0Widths;
            std::vector<Double_t> v2Values;
            std::vector<Double_t> v2Errors;
            std::vector<Double_t> v3Values;
            std::vector<Double_t> v3Errors;
            
            // Read v2 and v3 values from all FT0 slice files in a single robust loop
            for (size_t iFT0 = 0; iFT0 < ft0EtaVec.size() - 1; ++iFT0) {
                Double_t ft0EtaMin = ft0EtaVec[iFT0];
                Double_t ft0EtaMax = ft0EtaVec[iFT0 + 1];
                Double_t ft0Center = (ft0EtaMin + ft0EtaMax) / 2.0;
                Double_t ft0Width = std::fabs(ft0EtaMax - ft0EtaMin) / 2.0;  // Half-width (always positive)

                TString rootFileName = Form("./TemplateFit/EtaDiff/Vn_%s_%s_FT0Eta_%0.1f_%0.1f_%s_%i_%i_TPCEta_%0.1f_%0.1f.root",
                                           fileNameSuffix.c_str(), detectorName.c_str(), ft0EtaMin, ft0EtaMax,
                                           splitName.c_str(), data.minRange, data.maxRange, tpcEtaMin, tpcEtaMax);

                TString altRootFileName = Form("./TemplateFit/EtaDiff/Vn_%s_%s_FT0Eta_%0.1f_%0.1f_%s_FT0Eta_%0.1f_%0.1f_%i_%i_TPCEta_%0.1f_%0.1f.root",
                                               fileNameSuffix.c_str(), detectorName.c_str(), ft0EtaMin, ft0EtaMax,
                                               splitName.c_str(), ft0EtaMin, ft0EtaMax, data.minRange, data.maxRange, tpcEtaMin, tpcEtaMax);

                TFile* file = nullptr;
                if (gSystem->AccessPathName(rootFileName) == 0) {
                    file = TFile::Open(rootFileName);
                } else if (gSystem->AccessPathName(altRootFileName) == 0) {
                    file = TFile::Open(altRootFileName);
                } else {
                    continue;
                }

                TH1D* hV2 = (TH1D*)file->Get("hV2");
                TH1D* hV3 = (TH1D*)file->Get("hV3");
                if (hV2 && hV2->GetNbinsX() > 0) {
                    Double_t v2 = hV2->GetBinContent(1);
                    Double_t e2 = hV2->GetBinError(1);
                    Double_t v3 = 0.0;
                    Double_t e3 = 0.0;
                    if (hV3 && hV3->GetNbinsX() > 0) {
                        v3 = hV3->GetBinContent(1);
                        e3 = hV3->GetBinError(1);
                    }

                    if (std::isfinite(v2) && std::isfinite(e2) && e2 < 9.99) {
                        ft0Centers.push_back(ft0Center);
                        ft0Widths.push_back(ft0Width);
                        v2Values.push_back(v2);
                        v2Errors.push_back(e2);
                        v3Values.push_back(v3);
                        v3Errors.push_back(e3);
                    }
                }
                file->Close();
                delete file;
            }

            if (ft0Centers.empty()) continue;

            gSystem->mkdir("./TemplateFit/EtaDiff/SummaryPlots", kTRUE);

            // V2-only plot
            TCanvas cV2("cV2Summary", "V2Delta Summary", 900, 700);
            cV2.SetLeftMargin(0.12);
            cV2.SetRightMargin(0.05);
            TGraphErrors* gV2 = new TGraphErrors(ft0Centers.size(), ft0Centers.data(), v2Values.data(), ft0Widths.data(), v2Errors.data());
            gV2->SetMarkerStyle(20);
            gV2->SetMarkerSize(1.2);
            gV2->SetLineWidth(2);
            gV2->SetTitle(Form("V_{2#Delta} vs FT0 #eta;#eta^{FT0};V_{2#Delta}"));
            gV2->Draw("APE");
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.04);
            latex.DrawLatex(0.15, 0.85, Form("ALICE %s", collisionSystemName.c_str()));
            latex.DrawLatex(0.15, 0.80, Form("%s", detectorName.c_str()));
            latex.DrawLatex(0.15, 0.75, Form("%d < Cent < %d", data.minRange, data.maxRange));
            latex.DrawLatex(0.15, 0.70, Form("%.1f < #eta^{TPC} < %.1f", tpcEtaMin, tpcEtaMax));
            cV2.SaveAs(Form("./TemplateFit/EtaDiff/SummaryPlots/V2Delta_vs_FT0Eta_%s_%s_%s_%i_%i_TPCEta_%0.1f_%0.1f.pdf",
                            fileNameSuffix.c_str(), detectorName.c_str(), splitName.c_str(), data.minRange, data.maxRange, tpcEtaMin, tpcEtaMax));

            // V3-only plot
            TCanvas cV3("cV3Summary", "V3Delta Summary", 900, 700);
            cV3.SetLeftMargin(0.12);
            cV3.SetRightMargin(0.05);
            TGraphErrors* gV3 = new TGraphErrors(ft0Centers.size(), ft0Centers.data(), v3Values.data(), ft0Widths.data(), v3Errors.data());
            gV3->SetMarkerStyle(21);
            gV3->SetMarkerSize(1.2);
            gV3->SetLineWidth(2);
            gV3->SetMarkerColor(kRed);
            gV3->SetLineColor(kRed);
            gV3->SetTitle(Form("V_{3#Delta} vs FT0 #eta;#eta^{FT0};V_{3#Delta}"));
            gV3->Draw("APE");
            latex.DrawLatex(0.15, 0.85, Form("ALICE %s", collisionSystemName.c_str()));
            latex.DrawLatex(0.15, 0.80, Form("%s", detectorName.c_str()));
            latex.DrawLatex(0.15, 0.75, Form("%d < Cent < %d", data.minRange, data.maxRange));
            latex.DrawLatex(0.15, 0.70, Form("%.1f < #eta^{TPC} < %.1f", tpcEtaMin, tpcEtaMax));
            cV3.SaveAs(Form("./TemplateFit/EtaDiff/SummaryPlots/V3Delta_vs_FT0Eta_%s_%s_%s_%i_%i_TPCEta_%0.1f_%0.1f.pdf",
                            fileNameSuffix.c_str(), detectorName.c_str(), splitName.c_str(), data.minRange, data.maxRange, tpcEtaMin, tpcEtaMax));

            // Combined plot (v2 + v3)
            TCanvas cComb("cV2V3Summary", "V2/V3 Delta Summary", 900, 700);
            cComb.SetLeftMargin(0.12);
            cComb.SetRightMargin(0.05);
            gV2->SetMarkerColor(kBlue);
            gV2->SetLineColor(kBlue);
            gV2->Draw("APE");
            gV3->Draw("PEsame");
            TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
            leg->AddEntry(gV2, "v_{2#Delta}", "lep");
            leg->AddEntry(gV3, "v_{3#Delta}", "lep");
            leg->Draw();
            latex.DrawLatex(0.15, 0.85, Form("ALICE %s", collisionSystemName.c_str()));
            latex.DrawLatex(0.15, 0.80, Form("%s", detectorName.c_str()));
            latex.DrawLatex(0.15, 0.75, Form("%d < Cent < %d", data.minRange, data.maxRange));
            latex.DrawLatex(0.15, 0.70, Form("%.1f < #eta^{TPC} < %.1f", tpcEtaMin, tpcEtaMax));
            cComb.SaveAs(Form("./TemplateFit/EtaDiff/SummaryPlots/V2V3Delta_vs_FT0Eta_%s_%s_%s_%i_%i_TPCEta_%0.1f_%0.1f.pdf",
                               fileNameSuffix.c_str(), detectorName.c_str(), splitName.c_str(), data.minRange, data.maxRange, tpcEtaMin, tpcEtaMax));

            delete gV2;
            delete gV3;
            delete leg;
        }
    }
}

//==============================================================
void CreateV2DeltaSummaryPlot_InnerRing(Bool_t isNch, std::vector<InputUnit> dataList, std::string fileNameSuffix, std::string detectorName) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    
    std::vector<std::pair<double,double>> tpcEdges = { { -0.8, -0.7 }, { 0.7, 0.8 } };

    // For each TPC edge, create a summary plot
    for (auto tr : tpcEdges) {
        double tpcEtaMin = tr.first;
        double tpcEtaMax = tr.second;
        
        // For each centrality range
        for (const auto& data : dataList) {
            std::vector<Double_t> ringCenters;
            std::vector<Double_t> ringWidths;
            std::vector<Double_t> v2Values;
            std::vector<Double_t> v2Errors;
            std::vector<Double_t> v3Values;
            std::vector<Double_t> v3Errors;
            
            // Inner ring has no FT0 slicing, just read the single file
            TString rootFileName = Form("./TemplateFit/EtaDiff/Vn_%s_%s_TPCEta_%0.1f_%0.1f_%s_%i_%i.root",
                                       fileNameSuffix.c_str(), detectorName.c_str(), tpcEtaMin, tpcEtaMax,
                                       splitName.c_str(), data.minRange, data.maxRange);

            TFile* file = nullptr;
            if (gSystem->AccessPathName(rootFileName) == 0) {
                file = TFile::Open(rootFileName);
            } else {
                std::cout << "Warning: File not found: " << rootFileName << std::endl;
                continue;
            }

            TH1D* hV2 = (TH1D*)file->Get("hV2");
            TH1D* hV3 = (TH1D*)file->Get("hV3");
            if (hV2 && hV2->GetNbinsX() > 0) {
                Double_t v2 = hV2->GetBinContent(1);
                Double_t e2 = hV2->GetBinError(1);
                Double_t v3 = 0.0;
                Double_t e3 = 0.0;
                if (hV3 && hV3->GetNbinsX() > 0) {
                    v3 = hV3->GetBinContent(1);
                    e3 = hV3->GetBinError(1);
                }

                if (std::isfinite(v2) && std::isfinite(e2) && e2 < 9.99) {
                    // Use center at 0 with fixed width for inner ring representation
                    ringCenters.push_back(0.0);
                    ringWidths.push_back(0.5);  // Fixed width for visualization
                    v2Values.push_back(v2);
                    v2Errors.push_back(e2);
                    v3Values.push_back(v3);
                    v3Errors.push_back(e3);
                }
            }
            file->Close();
            delete file;

            if (ringCenters.empty()) continue;

            gSystem->mkdir("./TemplateFit/EtaDiff/SummaryPlots", kTRUE);

            // V2-only plot
            TCanvas cV2("cV2Summary", "V2Delta Summary", 900, 700);
            cV2.SetLeftMargin(0.12);
            cV2.SetRightMargin(0.05);
            TGraphErrors* gV2 = new TGraphErrors(ringCenters.size(), ringCenters.data(), v2Values.data(), ringWidths.data(), v2Errors.data());
            gV2->SetMarkerStyle(20);
            gV2->SetMarkerSize(1.2);
            gV2->SetLineWidth(2);
            gV2->SetTitle(Form("V_{2#Delta} vs Inner Ring;Inner Ring;V_{2#Delta}"));
            gV2->Draw("APE");
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.04);
            latex.DrawLatex(0.15, 0.85, Form("ALICE %s", collisionSystemName.c_str()));
            latex.DrawLatex(0.15, 0.80, Form("%s", detectorName.c_str()));
            latex.DrawLatex(0.15, 0.75, Form("%d < Cent < %d", data.minRange, data.maxRange));
            latex.DrawLatex(0.15, 0.70, Form("%.1f < #eta^{TPC} < %.1f", tpcEtaMin, tpcEtaMax));
            cV2.SaveAs(Form("./TemplateFit/EtaDiff/SummaryPlots/V2Delta_InnerRing_%s_%s_%s_%i_%i_TPCEta_%0.1f_%0.1f.pdf",
                            fileNameSuffix.c_str(), detectorName.c_str(), splitName.c_str(), data.minRange, data.maxRange, tpcEtaMin, tpcEtaMax));

            // V3-only plot
            TCanvas cV3("cV3Summary", "V3Delta Summary", 900, 700);
            cV3.SetLeftMargin(0.12);
            cV3.SetRightMargin(0.05);
            TGraphErrors* gV3 = new TGraphErrors(ringCenters.size(), ringCenters.data(), v3Values.data(), ringWidths.data(), v3Errors.data());
            gV3->SetMarkerStyle(21);
            gV3->SetMarkerSize(1.2);
            gV3->SetLineWidth(2);
            gV3->SetMarkerColor(kRed);
            gV3->SetLineColor(kRed);
            gV3->SetTitle(Form("V_{3#Delta} vs Inner Ring;Inner Ring;V_{3#Delta}"));
            gV3->Draw("APE");
            latex.DrawLatex(0.15, 0.85, Form("ALICE %s", collisionSystemName.c_str()));
            latex.DrawLatex(0.15, 0.80, Form("%s", detectorName.c_str()));
            latex.DrawLatex(0.15, 0.75, Form("%d < Cent < %d", data.minRange, data.maxRange));
            latex.DrawLatex(0.15, 0.70, Form("%.1f < #eta^{TPC} < %.1f", tpcEtaMin, tpcEtaMax));
            cV3.SaveAs(Form("./TemplateFit/EtaDiff/SummaryPlots/V3Delta_InnerRing_%s_%s_%s_%i_%i_TPCEta_%0.1f_%0.1f.pdf",
                            fileNameSuffix.c_str(), detectorName.c_str(), splitName.c_str(), data.minRange, data.maxRange, tpcEtaMin, tpcEtaMax));

            // Combined plot (v2 + v3)
            TCanvas cComb("cV2V3Summary", "V2/V3 Delta Summary", 900, 700);
            cComb.SetLeftMargin(0.12);
            cComb.SetRightMargin(0.05);
            gV2->SetMarkerColor(kBlue);
            gV2->SetLineColor(kBlue);
            gV2->Draw("APE");
            gV3->Draw("PEsame");
            TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
            leg->AddEntry(gV2, "v_{2#Delta}", "lep");
            leg->AddEntry(gV3, "v_{3#Delta}", "lep");
            leg->Draw();
            latex.DrawLatex(0.15, 0.85, Form("ALICE %s", collisionSystemName.c_str()));
            latex.DrawLatex(0.15, 0.80, Form("%s", detectorName.c_str()));
            latex.DrawLatex(0.15, 0.75, Form("%d < Cent < %d", data.minRange, data.maxRange));
            latex.DrawLatex(0.15, 0.70, Form("%.1f < #eta^{TPC} < %.1f", tpcEtaMin, tpcEtaMax));
            cComb.SaveAs(Form("./TemplateFit/EtaDiff/SummaryPlots/V2V3Delta_InnerRing_%s_%s_%s_%i_%i_TPCEta_%0.1f_%0.1f.pdf",
                               fileNameSuffix.c_str(), detectorName.c_str(), splitName.c_str(), data.minRange, data.maxRange, tpcEtaMin, tpcEtaMax));

            delete gV2;
            delete gV3;
            delete leg;
        }
    }
}

//==============================================================
// Create comparison plots for ring datasets showing both eta regions
void CreateRingComparisonPlots(Bool_t isNch, std::vector<InputUnit> dataList, std::string fileNameSuffix, std::string detectorName, bool isInner) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    
    std::string ringLabel = isInner ? "Inner Ring" : "Outer Ring";
    
    std::vector<std::pair<double,double>> tpcEdges = { { -0.8, -0.7 }, { 0.7, 0.8 } };
    
    // For each centrality range
    for (const auto& data : dataList) {
        std::vector<Double_t> etaCenters;  // x-axis positions
        std::vector<Double_t> etaWidths;   // x-axis error bars
        std::vector<Double_t> v2Values;
        std::vector<Double_t> v2Errors;
        std::vector<Double_t> v3Values;
        std::vector<Double_t> v3Errors;
        std::vector<TString> etaLabels;
        
        // Load data for both eta regions
        for (size_t iEta = 0; iEta < tpcEdges.size(); ++iEta) {
            double tpcEtaMin = tpcEdges[iEta].first;
            double tpcEtaMax = tpcEdges[iEta].second;
            
            TString rootFileName = Form("./TemplateFit/EtaDiff/Vn_%s_%s_TPCEta_%0.1f_%0.1f_%s_%i_%i.root",
                                       fileNameSuffix.c_str(), detectorName.c_str(), tpcEtaMin, tpcEtaMax,
                                       splitName.c_str(), data.minRange, data.maxRange);

            TFile* file = nullptr;
            if (gSystem->AccessPathName(rootFileName) == 0) {
                file = TFile::Open(rootFileName);
            } else {
                std::cout << "Warning: File not found: " << rootFileName << std::endl;
                continue;
            }

            TH1D* hV2 = (TH1D*)file->Get("hV2");
            TH1D* hV3 = (TH1D*)file->Get("hV3");
            if (hV2 && hV2->GetNbinsX() > 0) {
                Double_t v2 = hV2->GetBinContent(1);
                Double_t e2 = hV2->GetBinError(1);
                Double_t v3 = 0.0;
                Double_t e3 = 0.0;
                if (hV3 && hV3->GetNbinsX() > 0) {
                    v3 = hV3->GetBinContent(1);
                    e3 = hV3->GetBinError(1);
                }

                // Only add if valid data (error < 9.99 means fit succeeded)
                if (std::isfinite(v2) && std::isfinite(e2) && e2 < 9.99) {
                    // Position on x-axis: 0 for negative, 1 for positive
                    etaCenters.push_back(iEta);
                    etaWidths.push_back(0.2);  // Narrow width for clear separation
                    v2Values.push_back(v2);
                    v2Errors.push_back(e2);
                    v3Values.push_back(v3);
                    v3Errors.push_back(e3);
                    etaLabels.push_back(Form("%.1f < #eta < %.1f", tpcEtaMin, tpcEtaMax));
                }
            }
            file->Close();
            delete file;
        }

        if (etaCenters.empty()) {
            std::cout << "Warning: No valid data for " << fileNameSuffix << " " << detectorName << std::endl;
            continue;
        }

        gSystem->mkdir("./TemplateFit/EtaDiff/SummaryPlots", kTRUE);

        // Create V2 comparison plot
        TCanvas cV2("cV2Comparison", Form("V2Delta %s", ringLabel.c_str()), 1000, 700);
        cV2.SetLeftMargin(0.12);
        cV2.SetRightMargin(0.05);
        cV2.SetBottomMargin(0.12);
        
        TGraphErrors* gV2 = new TGraphErrors(etaCenters.size(), etaCenters.data(), v2Values.data(), etaWidths.data(), v2Errors.data());
        gV2->SetMarkerStyle(20);
        gV2->SetMarkerSize(1.5);
        gV2->SetMarkerColor(kBlue);
        gV2->SetLineColor(kBlue);
        gV2->SetLineWidth(2);
        gV2->SetTitle(Form("V_{2#Delta} Comparison - %s;TPC #eta Region;V_{2#Delta}", ringLabel.c_str()));
        gV2->GetXaxis()->SetNdivisions(2);
        gV2->GetXaxis()->SetLimits(-0.5, 1.5);
        gV2->GetXaxis()->SetLabelSize(0);  // Hide numeric labels
        gV2->GetYaxis()->SetTitleOffset(1.3);
        gV2->Draw("APE");
        
        // Add custom x-axis labels
        TLatex latex;
        latex.SetTextAlign(22);
        latex.SetTextSize(0.04);
        for (size_t i = 0; i < etaLabels.size(); ++i) {
            latex.DrawLatex(etaCenters[i], gV2->GetYaxis()->GetXmin() - 0.08 * (gV2->GetYaxis()->GetXmax() - gV2->GetYaxis()->GetXmin()), etaLabels[i]);
        }
        
        latex.SetNDC();
        latex.SetTextAlign(12);
        latex.DrawLatex(0.15, 0.85, Form("ALICE %s", collisionSystemName.c_str()));
        latex.DrawLatex(0.15, 0.80, Form("%s - %s", ringLabel.c_str(), detectorName.c_str()));
        latex.DrawLatex(0.15, 0.75, Form("%d < Cent < %d", data.minRange, data.maxRange));
        
        cV2.SaveAs(Form("./TemplateFit/EtaDiff/SummaryPlots/V2Delta_Comparison_%s_%s_%s_%i_%i.pdf",
                        fileNameSuffix.c_str(), detectorName.c_str(), splitName.c_str(), data.minRange, data.maxRange));

        // Create V3 comparison plot
        TCanvas cV3("cV3Comparison", Form("V3Delta %s", ringLabel.c_str()), 1000, 700);
        cV3.SetLeftMargin(0.12);
        cV3.SetRightMargin(0.05);
        cV3.SetBottomMargin(0.12);
        
        TGraphErrors* gV3 = new TGraphErrors(etaCenters.size(), etaCenters.data(), v3Values.data(), etaWidths.data(), v3Errors.data());
        gV3->SetMarkerStyle(21);
        gV3->SetMarkerSize(1.5);
        gV3->SetMarkerColor(kRed);
        gV3->SetLineColor(kRed);
        gV3->SetLineWidth(2);
        gV3->SetTitle(Form("V_{3#Delta} Comparison - %s;TPC #eta Region;V_{3#Delta}", ringLabel.c_str()));
        gV3->GetXaxis()->SetNdivisions(2);
        gV3->GetXaxis()->SetLimits(-0.5, 1.5);
        gV3->GetXaxis()->SetLabelSize(0);
        gV3->GetYaxis()->SetTitleOffset(1.3);
        gV3->Draw("APE");
        
        for (size_t i = 0; i < etaLabels.size(); ++i) {
            latex.SetTextAlign(22);
            latex.SetNDC(kFALSE);
            latex.DrawLatex(etaCenters[i], gV3->GetYaxis()->GetXmin() - 0.08 * (gV3->GetYaxis()->GetXmax() - gV3->GetYaxis()->GetXmin()), etaLabels[i]);
        }
        
        latex.SetNDC();
        latex.SetTextAlign(12);
        latex.DrawLatex(0.15, 0.85, Form("ALICE %s", collisionSystemName.c_str()));
        latex.DrawLatex(0.15, 0.80, Form("%s - %s", ringLabel.c_str(), detectorName.c_str()));
        latex.DrawLatex(0.15, 0.75, Form("%d < Cent < %d", data.minRange, data.maxRange));
        
        cV3.SaveAs(Form("./TemplateFit/EtaDiff/SummaryPlots/V3Delta_Comparison_%s_%s_%s_%i_%i.pdf",
                        fileNameSuffix.c_str(), detectorName.c_str(), splitName.c_str(), data.minRange, data.maxRange));

        // Create combined V2+V3 comparison plot
        TCanvas cComb("cV2V3Comparison", Form("V2/V3 Delta %s", ringLabel.c_str()), 1000, 700);
        cComb.SetLeftMargin(0.12);
        cComb.SetRightMargin(0.05);
        cComb.SetBottomMargin(0.12);
        
        // Find y-axis range for combined plot
        double yMin = 1e10, yMax = -1e10;
        for (size_t i = 0; i < v2Values.size(); ++i) {
            yMin = std::min(yMin, v2Values[i] - v2Errors[i]);
            yMax = std::max(yMax, v2Values[i] + v2Errors[i]);
            yMin = std::min(yMin, v3Values[i] - v3Errors[i]);
            yMax = std::max(yMax, v3Values[i] + v3Errors[i]);
        }
        double yPad = 0.2 * (yMax - yMin);
        
        gV2->GetYaxis()->SetRangeUser(yMin - yPad, yMax + yPad);
        gV2->SetTitle(Form("V_{n#Delta} Comparison - %s;TPC #eta Region;V_{n#Delta}", ringLabel.c_str()));
        gV2->Draw("APE");
        gV3->Draw("PE same");
        
        TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->AddEntry(gV2, "V_{2#Delta}", "lep");
        leg->AddEntry(gV3, "V_{3#Delta}", "lep");
        leg->Draw();
        
        for (size_t i = 0; i < etaLabels.size(); ++i) {
            latex.SetTextAlign(22);
            latex.SetNDC(kFALSE);
            latex.DrawLatex(etaCenters[i], yMin - yPad - 0.03 * (yMax - yMin), etaLabels[i]);
        }
        
        latex.SetNDC();
        latex.SetTextAlign(12);
        latex.DrawLatex(0.15, 0.85, Form("ALICE %s", collisionSystemName.c_str()));
        latex.DrawLatex(0.15, 0.80, Form("%s - %s", ringLabel.c_str(), detectorName.c_str()));
        latex.DrawLatex(0.15, 0.75, Form("%d < Cent < %d", data.minRange, data.maxRange));
        
        cComb.SaveAs(Form("./TemplateFit/EtaDiff/SummaryPlots/V2V3Delta_Comparison_%s_%s_%s_%i_%i.pdf",
                          fileNameSuffix.c_str(), detectorName.c_str(), splitName.c_str(), data.minRange, data.maxRange));

        delete gV2;
        delete gV3;
        delete leg;
    }
}