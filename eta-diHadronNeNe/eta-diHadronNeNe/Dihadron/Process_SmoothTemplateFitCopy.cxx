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
//#include <TRandom3.h>
#include "TMath.h"
#include "THnSparse.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
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
    Double_t v1;
    Double_t v1_err;
    Double_t chi2ndf;
    Double_t chi2ndf_err;
    VnUnit(Double_t _v2, Double_t _v2_err, Double_t _v3, Double_t _v3_err, Double_t _v4, Double_t _v4_err, Double_t _F, Double_t _F_err, Double_t _G, Double_t _G_err, Double_t _v1, Double_t _v1_err, Double_t _chi2ndf, Double_t _chi2ndf_err) :
        v2(_v2), v2_err(_v2_err), v3(_v3), v3_err(_v3_err), v4(_v4), v4_err(_v4_err), F(_F), F_err(_F_err), G(_G), G_err(_G_err), v1(_v1), v1_err(_v1_err), chi2ndf(_chi2ndf), chi2ndf_err(_chi2ndf_err) {}
};

// declare functions
void ProcessConfig(Bool_t isNch, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName);
void ProcessConfig_EtaDiff(Bool_t isNch, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName);
VnUnit* TemplateFit(Bool_t isNch, InputUnit templ, InputUnit data, Bool_t cn2Tovn2, Double_t etaMin=-101., Double_t etaMax=-101.);
std::vector<Int_t> CheckAndMergeRanges(const std::vector<InputUnit>& inputUnits);
VnUnit* fitSample(Bool_t isNch, TFile* templatefile, InputUnit templ, TFile* datafile, InputUnit data, int sample = -1, Double_t etaMin=-101., Double_t etaMax=-101.);
void RooTempFitter(TH1 *lm, TH1 *hm, std::vector<Double_t>& fParamVal, std::vector<Double_t>& fParamErr, Bool_t kRefit);
void PlotFitting(TH1 *lm, TH1 *hm, Bool_t isNch, std::string fileSuffix, Int_t minRange, Int_t maxRange, const std::vector<Double_t>& par, const std::vector<Double_t>& parerr, Double_t etaMin=-101., Double_t etaMax=-101.);
void writeToCSV(const std::string &filename, const std::vector<std::vector<std::string>> &data, const std::string &folder, bool createHeader=false); //mine
double findChi2ndf(TH1 *lm, TH1 *hm, double& F_val, double& G_val, double& v11_val, double& v21_val, double& v31_val, double& v41_val);
std::vector<std::string> toStringV2(std::vector<std::vector<std::vector<double>>> vec, int Nelements, double etaMin, double etaMax);
std::vector<std::string> toStringV3(std::vector<std::vector<std::vector<double>>> vec, int Nelements, double etaMin, double etaMax);
std::vector<std::string> toStringV4(std::vector<std::vector<std::vector<double>>> vec, int Nelements, double etaMin, double etaMax);
void PlotTemplate(TH1 *hm, TH1 *fittedTemp, Bool_t isNch, std::string fileSuffix, Int_t minRange, Int_t maxRange, const std::vector<Double_t>& par, const std::vector<Double_t>& parerr, Double_t etaMin=-101., Double_t etaMax=-101.);
TFitResultPtr RooFourierFit(TH1 *hist, std::vector<Double_t>& fParamVal, std::vector<Double_t>& fParamErr);
TH1D* getFittedTemplate(TFitResultPtr fitResult, std::vector<Double_t> par, int Nbins);
void DrawText(double xmin, double ymin, double textSize, TString text);

// global variables
std::string collisionSystemName = "peripheral Ne-Ne"; //skal byttes for OO
std::string FITused = "TPC-FT0C";
collisionSystemName = "Ne-Ne";

//==============================================================
void Process_SmoothTemplateFit() {
    // 不显示窗口
    gROOT->SetBatch(kTRUE);
    std::vector<ConfigUnit> configList;


    configList.push_back(ConfigUnit(kCent, kTPCFT0C,  kEtaDiffOff, InputUnit("LHC25ae_pass2_598682", 80, 100), 
    {InputUnit("LHC25ae_pass2_598682", 0, 20)},
    "LHC25ae_pass2_598682"));

    configList.push_back(ConfigUnit(kCent, kTPCFT0A,  kEtaDiffOff, InputUnit("LHC25ae_pass2_598682", 80, 100), 
    {InputUnit("LHC25ae_pass2_598682", 0, 20)},
    "LHC25ae_pass2_598682"));

    configList.push_back(ConfigUnit(kCent, kTPCFT0C,  kEtaDiffOn, InputUnit("LHC25ae_pass2_598682", 80, 100), 
    {InputUnit("LHC25ae_pass2_598682", 0, 20)},
    "LHC25ae_pass2_598682"));

    configList.push_back(ConfigUnit(kCent, kTPCFT0A,  kEtaDiffOn, InputUnit("LHC25ae_pass2_598682", 80, 100), 
    {InputUnit("LHC25ae_pass2_598682", 0, 20)},
    "LHC25ae_pass2_598682"));

    //her itererer vi over konfigurasjonene
    for (auto config : configList) {
        if (config.isEtaDiff) { //eta-diff case
            ProcessConfig_EtaDiff(config.isNch, config.templ, config.dataList, config.outputFileName);
        } else { //ikke eta-diff case
            ProcessConfig(config.isNch, config.templ, config.dataList, config.outputFileName);
        }
    }

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
        vnResults.push_back(TemplateFit(isNch, templ, data, kFALSE)); //her finner vi hva v2 er, kTRUE betyr at cn2 to vn2 er aktivert?
    } //her står vi igjen med V_Delta n (ref, ref)

    // 创建输出文件
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    TFile outputFile(Form("./TemplateFit/Vn_%s_%s.root", outputFileName.c_str(), splitName.c_str()), "RECREATE"); //har lages .root for V_Delta n (ref, ref)

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
void ProcessConfig_EtaDiff(Bool_t isNch, InputUnit templ, std::vector<InputUnit> dataList, std::string outputFileName) {
    // Just looping the list
    for (const auto& data : dataList) {
        // 执行模板拟合获取所有结果
        std::vector<VnUnit*> vnResults;
        for (Int_t iEta = 0; iEta < etaTPC.size() - 1; iEta++) {
            double etaMin = etaTPC[iEta];
            double etaMax = etaTPC[iEta + 1];
            vnResults.push_back(TemplateFit(isNch, templ, data, kFALSE, etaMin, etaMax)); //men ikke her cn2 to vn2
        }

        // 创建输出文件
        std::string splitName = "Mult";
        if (!isNch) splitName = "Cent";
        TFile outputFile(Form("./TemplateFit/EtaDiff/Vn_%s_%s_%i_%i.root", outputFileName.c_str(), splitName.c_str(), data.minRange, data.maxRange), "RECREATE"); //har skapes .root for diff vns

        // 初始化直方图
        TH1D* hV2 = new TH1D("hV2", "V_{ #Delta 2};#eta^{TPC};V_{ #Delta 2}", 
                            etaTPC.size()-1, etaTPC.data());
        TH1D* hV3 = new TH1D("hV3", "V_{ #Delta 3};#eta^{TPC};V_{ #Delta 3}", 
                            etaTPC.size()-1, etaTPC.data());
        TH1D* hV4 = new TH1D("hV4", "V_{ #Delta 4};#eta^{TPC};V_{ #Delta 4}", 
                            etaTPC.size()-1, etaTPC.data());

        // 填充数据
        for (size_t i = 0; i < vnResults.size(); ++i) {
            hV2->SetBinContent(i+1, vnResults[i]->v2); //her hentes verdi, igjen uten noe normalisering !!!sjekk hvor den hentes fra
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

        std::cout << "Output file: " << Form("./TemplateFit/EtaDiff/Vn_%s_%s_%i_%i.root", outputFileName.c_str(), splitName.c_str(), data.minRange, data.maxRange) << std::endl;
        outputFile.Close();

    }

}



//==============================================================
VnUnit* TemplateFit(Bool_t isNch, InputUnit templ, InputUnit data, Bool_t cn2Tovn2, Double_t etaMin=-101., Double_t etaMax=-101.) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    TFile* templatefile = new TFile(Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d.root", templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange), "READ");
    if (!templatefile || !templatefile->IsOpen()) { //error
        std::cerr << "Cannot open template file: " << Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d.root", templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange) << std::endl;
        exit(1);
    }
    TFile* templatefile_etaDiff = nullptr;
    if (etaMin > -100. && etaMax > -100.) { //dersom pt-diff er aktivert
        templatefile_etaDiff = new TFile(Form("./ProcessOutput/etaDiff/BootstrapSample_%s_%s_%d_%d_Eta_%0.1f_%0.1f.root", templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange, etaMin, etaMax), "READ");
        if (!templatefile_etaDiff || !templatefile_etaDiff->IsOpen()) { //error
            std::cerr << "Cannot open template file: " << Form("./ProcessOutput/etaDiff/BootstrapSample_%s_%s_%d_%d_Eta_%0.1f_%0.1f.root", templ.fileNameSuffix.c_str(), splitName.c_str(), templ.minRange, templ.maxRange, etaMin, etaMax) << std::endl;
            exit(1);
        }
    }

    TFile* datafile = new TFile(Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange), "READ");
    if (!datafile || !datafile->IsOpen()) { //error
        std::cerr << "Cannot open input file: " << Form("./ProcessOutput/BootstrapSample_%s_%s_%d_%d.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange) << std::endl;
        exit(1);
    }
    TFile* datafile_etaDiff = nullptr;
    if (etaMin > -100. && etaMax > -100.) { //dersom eta-diff er aktivert
        datafile_etaDiff = new TFile(Form("./ProcessOutput/etaDiff/BootstrapSample_%s_%s_%i_%i_Eta_%0.1f_%0.1f.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, etaMin, etaMax), "READ");
        if (!datafile_etaDiff || !datafile_etaDiff->IsOpen()) { //error
            std::cerr << "Cannot open input file: " << Form("./ProcessOutput/etaDiff/BootstrapSample_%s_%s_%i_%i_Eta_%0.1f_%0.1f.root", data.fileNameSuffix.c_str(), splitName.c_str(), data.minRange, data.maxRange, etaMin, etaMax) << std::endl;
            exit(1);
        }
    }

    VnUnit* vnResult = fitSample(isNch, templatefile, templ, datafile, data, -1); //her finner vi fit for ref
    if (!vnResult) { //error
        std::cerr << "Cannot fit sample: " << data.fileNameSuffix << std::endl;
        exit(1);
    }
    VnUnit* vnResult_etaDiff = nullptr; //dersom etaDiff ikke er aktivt forblir denne nullptr
    if (etaMin > -100. && etaMax > -100.) { //dersom eta-diff er aktivert
        vnResult_etaDiff = fitSample(isNch, templatefile_etaDiff, templ, datafile_etaDiff, data, -1, etaMin, etaMax); //her finner vi fit for etaDiff
        if (!vnResult_etaDiff) {
            std::cerr << "Cannot fit eta-diff sample: " << data.fileNameSuffix << std::endl;
            exit(1);
        }
    }
    //EDIT START
    std::vector<std::vector<std::vector<double>>> ValueArray;
    std::vector<std::vector<std::vector<double>>> ValueErrorArray;
    std::vector<std::vector<double>> ErrorArray;
    int Nobs=7;//v22,v32,v42...v1,F,G,chi2ndf
    int NofSample = maxSample*maxSample; //100 DENNE VAR 1000
    int Nbin = 1;
    ResizeValueArray(ValueArray,ValueErrorArray,ErrorArray,Nobs,NofSample,Nbin);

    //her finner vi bootstrap samples (100 times)
    for(int sample=0;sample<NofSample;sample++) { //edited
        VnUnit* vnTemp = fitSample(isNch, templatefile, templ, datafile, data, sample);
        if (!vnTemp) { //error
            std::cerr << "Cannot fit sample: " << data.fileNameSuffix << " sample: " << sample << std::endl;
            exit(1);
        }
        //Here we might estimate the uncertainties for v11, F, G too
        ValueArray[0][sample][0] = vnTemp->v2; //3x100x1
        ValueErrorArray[0][sample][0] = vnTemp->v2_err;
        ValueArray[1][sample][0] = vnTemp->v3;
        ValueErrorArray[1][sample][0] = vnTemp->v3_err;
        ValueArray[2][sample][0] = vnTemp->v4;
        ValueErrorArray[2][sample][0] = vnTemp->v4_err;
        ValueArray[3][sample][0] = vnTemp->v1; 
        ValueErrorArray[3][sample][0] = vnTemp->v1_err;
        ValueArray[4][sample][0] = vnTemp->F;
        ValueErrorArray[4][sample][0] = vnTemp->F_err;
        ValueArray[5][sample][0] = vnTemp->G;
        ValueErrorArray[5][sample][0] = vnTemp->G_err;
        ValueArray[6][sample][0] = vnTemp->chi2ndf;
        ValueErrorArray[6][sample][0] = vnTemp->chi2ndf_err;
        if (etaMin > -100. && etaMax > -100.) { //dersom eta-diff er aktivert
            VnUnit* vnTemp_etaDiff = fitSample(isNch, templatefile_etaDiff, templ, datafile_etaDiff, data, sample, etaMin, etaMax);
            if (!vnTemp_etaDiff) { //error
                std::cerr << "Cannot fit eta-diff sample: " << data.fileNameSuffix << std::endl;
                exit(1);
            }
            ValueArray[0][sample][0] = vnTemp_etaDiff->v2;
            ValueErrorArray[0][sample][0] = vnTemp_etaDiff->v2_err;
            ValueArray[1][sample][0] = vnTemp_etaDiff->v3;
            ValueErrorArray[1][sample][0] = vnTemp_etaDiff->v3_err;
            ValueArray[2][sample][0] = vnTemp_etaDiff->v4;
            ValueErrorArray[2][sample][0] = vnTemp_etaDiff->v4_err;
            ValueArray[3][sample][0] = vnTemp_etaDiff->v1; 
            ValueErrorArray[3][sample][0] = vnTemp_etaDiff->v1_err;
            ValueArray[4][sample][0] = vnTemp_etaDiff->F;
            ValueErrorArray[4][sample][0] = vnTemp_etaDiff->F_err;
            ValueArray[5][sample][0] = vnTemp_etaDiff->G;
            ValueErrorArray[5][sample][0] = vnTemp_etaDiff->G_err;
            ValueArray[6][sample][0] = vnTemp_etaDiff->chi2ndf;
            ValueErrorArray[6][sample][0] = vnTemp_etaDiff->chi2ndf_err;
            delete vnTemp_etaDiff;
        }
        delete vnTemp;
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

    //her regnes bootstrap errors
    for(int iobs = 0;iobs < Nobs;iobs++){
        CalculateBootstrapError(ValueArray[iobs],ValueErrorArray[iobs],ErrorArray[iobs],1.); //regn ut errors
    }
    // we assign the boot-strap error to vnResult (the one shoed in the root files)
    vnResult->v2_err = ErrorArray[0][0];
    vnResult->v3_err = ErrorArray[1][0];
    vnResult->v4_err = ErrorArray[2][0];
    vnResult->v1_err = ErrorArray[3][0];
    vnResult->F_err = ErrorArray[4][0];
    vnResult->G_err = ErrorArray[5][0];
    vnResult->chi2ndf_err = ErrorArray[6][0];
    if (etaMin > -100. && etaMax > -100.) { //dersom eta-diff er aktivert
        vnResult_etaDiff->v2_err = ErrorArray[0][0];
        vnResult_etaDiff->v3_err = ErrorArray[1][0];
        vnResult_etaDiff->v4_err = ErrorArray[2][0];
        vnResult_etaDiff->v1_err = ErrorArray[3][0];
        vnResult_etaDiff->F_err = ErrorArray[4][0];
        vnResult_etaDiff->G_err = ErrorArray[5][0];
        vnResult_etaDiff->chi2ndf_err = ErrorArray[6][0];
    }

    if (cn2Tovn2) { //desrom kTURE da tar vi roten av Vn2. Bemerk at V er ikke cn per say siden V kommer fra FT, mens cn er observabel, denne er satt til false i alle kall
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
        double v11 = vnResult_etaDiff->v1;
        double v11e = vnResult_etaDiff->v1_err;
        double v21 = vnResult_etaDiff->v2;
        double v21e = vnResult_etaDiff->v2_err;
        double v31 = vnResult_etaDiff->v3;
        double v31e = vnResult_etaDiff->v3_err;
        double v41 = vnResult_etaDiff->v4;
        double v41e = vnResult_etaDiff->v4_err;
        double chi2ndf = vnResult_etaDiff->chi2ndf;
        double chi2ndfe = vnResult_etaDiff->chi2ndf_err;
        //fill the vector with data
        csvData.push_back({std::to_string(etaMin), std::to_string(etaMax), std::to_string(F), std::to_string(Fe), std::to_string(G), std::to_string(Ge), std::to_string(v11), std::to_string(v11e), std::to_string(v21), std::to_string(v21e), std::to_string(v31), std::to_string(v31e), std::to_string(v41), std::to_string(v41e), std::to_string(chi2ndf), std::to_string(chi2ndfe)});
        writeToCSV("BootStrapValues.csv", csvData, "./TemplateFit/etaDiff");
    } 
    else {
        double F = vnResult->F;
        double Fe = vnResult->F_err;
        double G = vnResult->G;
        double Ge = vnResult->G_err;
        double v11 = vnResult->v1;
        double v11e = vnResult->v1_err;
        double v21 = vnResult->v2;
        double v21e = vnResult->v2_err;
        double v31 = vnResult->v3;
        double v31e = vnResult->v3_err;
        double v41 = vnResult->v4;
        double v41e = vnResult->v4_err;
        double chi2ndf = vnResult->chi2ndf;
        double chi2ndfe = vnResult->chi2ndf_err;
        //fill the vector with data
        csvData.push_back({std::to_string(etaMin), std::to_string(etaMax), std::to_string(F), std::to_string(Fe), std::to_string(G), std::to_string(Ge), std::to_string(v11), std::to_string(v11e), std::to_string(v21), std::to_string(v21e), std::to_string(v31), std::to_string(v31e), std::to_string(v41), std::to_string(v41e), std::to_string(chi2ndf), std::to_string(chi2ndfe)});
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
VnUnit* fitSample(Bool_t isNch, TFile* templatefile, InputUnit templ, TFile* datafile, InputUnit data, int sample = -1, Double_t etaMin=-101., Double_t etaMax=-101.) {
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
     //Fitting template:
    std::vector<Double_t> fParamTemplate;
    std::vector<Double_t> fParamErrTemplate;
    TFitResultPtr fitResult = RooFourierFit(lm, fParamTemplate, fParamErrTemplate);
    TH1D* hFitTemp = getFittedTemplate(fitResult, fParamTemplate, lm->GetNbinsX());
    if (sample == -1){
      PlotTemplate(lm, hFitTemp, isNch, templ.fileNameSuffix, templ.minRange, templ.maxRange, fParamTemplate, fParamErrTemplate, etaMin, etaMax);}

    RooTempFitter(hFitTemp, hm, fParamVal, fParamErr, kFALSE); //FT magic happens here
    if (sample == -1) { //if -1 we load all of the samples, and only then will .pdf be created
        PlotFitting(hFitTemp, hm, isNch, data.fileNameSuffix, data.minRange, data.maxRange, fParamVal, fParamErr, etaMin, etaMax);
    }
    VnUnit* vnResult = new VnUnit(fParamVal[0], fParamErr[0], fParamVal[1], fParamErr[1], fParamVal[2], fParamErr[2], fParamVal[3], fParamErr[3], fParamVal[4], fParamErr[4], fParamVal[5], fParamErr[5], fParamVal[6], fParamErr[6]);
    return vnResult; //her sendes resultatet, altså verdier av v2,v3,v4 og errorer videre, altså fit slutter her
}


TH1D* getFittedTemplate(TFitResultPtr fitResult, std::vector<Double_t> par, int Nbins) {
    
    double xFitMin = -TMath::Pi() / 2.;
    double xFitMax = TMath::Pi() * 3. / 2.;

    TH1D* hFitTemp = new TH1D(
        "hFitTemp",
        "Fitted Template",
        Nbins,
        xFitMin,
        xFitMax
    );

    // Extract parameters
    int nParams = par.size();

    // Extract covariance matrix
    TMatrixDSym cov = fitResult->GetCovarianceMatrix();

    // Fill histogram and propagate errors
    for (int j = 1; j <= Nbins; j++) {
        double x = hFitTemp->GetBinCenter(j);

        // Compute Fourier series value: f(x) = p0 + 2*sum_{i=1}^{N-1} p_i * cos(i*x)
        double funcVal = par[0];
        for (int i = 1; i < nParams; i++) {
            funcVal += 2 * par[i] * TMath::Cos(i * x);
        }
        hFitTemp->SetBinContent(j, funcVal);

        // Compute propagated uncertainty using covariance matrix
        double sigma2 = 0.0;
        for (int i = 0; i < nParams; i++) {
            double dfdp_i = (i == 0) ? 1.0 : 2 * TMath::Cos(i * x);
            for (int k = 0; k < nParams; k++) {
                double dfdp_k = (k == 0) ? 1.0 : 2 * TMath::Cos(k * x);
                sigma2 += dfdp_i * dfdp_k * cov(i, k);
            }
        }

        hFitTemp->SetBinError(j, TMath::Sqrt(sigma2));
    }

    return hFitTemp;
}

TFitResultPtr RooFourierFit(TH1 *hist, std::vector<Double_t>& fParamVal, std::vector<Double_t>& fParamErr) {
    // 检查输入直方图是否有效
    if (!hist) {
        Error("FourierFit", "Invalid histogram pointer!");
    }

    int Nbins = hist->GetNbinsX();
    //int Npara = Nbins / 12.;
    int Npara = 5; //including a0 !!!
    
    fParamVal.clear();
    fParamErr.clear();
    fParamVal.resize(Npara);
    fParamErr.resize(Npara);

    // 定义傅里叶级数函数 - 5参数形式
    auto fourierFunc = [Npara](double *x, double *p) {
        auto funcVal = p[0];
        for (int i = 1; i<Npara; i++) {
          funcVal += 2*p[i]*TMath::Cos(i*x[0]);
        }
        return funcVal;
    };

    // 创建TF1对象
    TF1 *fitFunc = new TF1("fourierFit", fourierFunc, 
                          hist->GetXaxis()->GetXmin(),  // 使用直方图X范围
                          hist->GetXaxis()->GetXmax(), 
                          Npara);  // 参数数量

    std::vector<double> Coeff;
    HistFFT(hist, Coeff);

    fitFunc->SetParameter(0, Coeff[0]);  // 常数项初始化为平均值
    for (int i = 1; i < Npara; ++i) {
        fitFunc->SetParameter(i, Coeff[i]);  // 谐波项初始值设为RMS的10%
    }

    // 设置参数名称（可选）
    for (int i=0; i<Npara; i++) {
      fitFunc->SetParName(i, Form("a%i", i));
    }

    // 执行拟合（使用最小二乘法并抑制输出）
    // hist->Fit(fitFunc, "QN");
    TFitResultPtr fitResult = hist->Fit(fitFunc, "S");
    // 检查拟合状态
    if (static_cast<int>(fitResult)) {
        Error("FourierFit", "Fit failed with status %d", static_cast<int>(fitResult));
        for (int i = 0; i < Npara; ++i) {
            fParamVal[i] = -1;
            fParamErr[i] = 10;
        }
        return fitResult;
    }

    // 计算拟合质量指标
    double chi2 = fitFunc->GetChisquare();
    double ndf = fitFunc->GetNDF();
    double chi2_ndf = (ndf > 0) ? chi2 / ndf : 0;
    double pvalue = TMath::Prob(chi2, static_cast<int>(ndf));
    if (chi2_ndf > 2000) {
        std::cout << "WARNING: Chi2/NDF > 20, fit may be inaccurate!" << std::endl;
        for (int i = 0; i < Npara; ++i) {
            fParamVal[i] = -1;
            fParamErr[i] = 10;
        }
        return fitResult;
    }

    double a0 = fitFunc->GetParameter(0);
    double a0e = fitFunc->GetParError(0);

    fParamVal[0] = a0;
    fParamErr[0] = a0e;

    for (int i=1; i<Npara; i++) {
        fParamVal[i] = fitFunc->GetParameter(i);
        fParamErr[i] = fitFunc->GetParError(i);
    }
    
    return fitResult;
}


void PlotTemplate(TH1 *hm, TH1 *fittedTemp, Bool_t isNch, std::string fileSuffix, Int_t minRange, Int_t maxRange, const std::vector<Double_t>& par, const std::vector<Double_t>& parerr, Double_t etaMin=-101., Double_t etaMax=-101.) {
    gStyle->SetOptStat(0); 
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    
    double v12 = par[1] / par[0];
    double v21 = par[2] / par[0];
    double v31 = par[3] / par[0];
    double v41 = par[4] / par[0];
    double a0 = par[0];

    double v21e = Error_Ratio(a2, a2e, a0, a0e, 0);
    double v31e = Error_Ratio(a3, a3e, a0, a0e, 0);
    double v41e = Error_Ratio(a4, a4e, a0, a0e, 0);


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

    double a1 = v12 * a0;
    double a2 = v21 * a0;
    double a3 = v31 * a0;
    double a4 = v41 * a0;
    
    // =============== 新增：计算chi2/ndf ===============
    double chi2 = 0.0;
    int nBins = hm->GetNbinsX();
    int nParams = par.size(); // 参数个数: a0, a1, a2, a3, a4
    int ndf = nBins - nParams;
    
    for (int i = 1; i <= nBins; i++) {
        double data = hm->GetBinContent(i);
        double error = hm->GetBinError(i);
        double x = hm->GetBinCenter(i);
        double fit = fittedTemp->GetBinContent(i);
        
        if (error > 0) { // 忽略误差为0的bin
            double residual = data - fit;
            chi2 += (residual * residual) / (error * error);
        }
    }
    double chi2ndf = (ndf > 0) ? chi2 / ndf : 0;

    // 创建拟合曲线
    const int pointBin = (int)hm->GetNbinsX();
    Double_t CopyPointX[pointBin];
    Double_t CopyPointY[pointBin];
    Double_t CopyPointErrorY[pointBin];
    Double_t CopyPointErrorX[pointBin];
    for (int i=0; i<pointBin; ++i){
      CopyPointX[i] = hm->GetBinCenter(i+1);
      double x = hm->GetBinCenter(i+1);
      CopyPointY[i] = fittedTemp->GetBinContent(i+1);
      CopyPointErrorY[i] = fittedTemp->GetBinError(i+1);
      CopyPointErrorX[i] = 0.;
    };

    TGraphErrors* gCopy = new TGraphErrors(pointBin,CopyPointX,CopyPointY,CopyPointErrorX,CopyPointErrorY);

    Double_t PeriPointY[pointBin];
    for (int i=0; i<pointBin; ++i){
      PeriPointY[i] = a0;
    };
    TGraph* gPeri = new TGraph(pointBin,CopyPointX,PeriPointY);

    gCopy->SetLineColor(TColor::GetColor("#b30000"));
    gCopy->SetLineWidth(2);
    gCopy->Draw("same");

    gPeri->SetLineColor(colors[1]);
    gPeri->SetLineWidth(2);
    gPeri->SetLineStyle(5);
    gPeri->Draw("same");


    TF1* fit_p1 = new TF1("fit_p1","[0] + 2*[1]*cos(x)", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    TF1* fit_p2 = new TF1("fit_p2","[0] + 2*[1]*cos(2*x)", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    TF1* fit_p3 = new TF1("fit_p3","[0] + 2*[1]*cos(3*x)", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    TF1* fit_p4 = new TF1("fit_p4","[0] + 2*[1]*cos(4*x)", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    fit_p1->SetParameters(a0, a1);
    fit_p2->SetParameters(a0, a2);
    fit_p3->SetParameters(a0, a3);
    fit_p4->SetParameters(a0, a4);
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
    //fit_p2->Draw("same");
    //fit_p3->Draw("same");
    //fit_p4->Draw("same");

    // 添加图例
    TLegend* leg = new TLegend(0.5, 0.75, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->AddEntry(hm, "Data", "lep");
    leg->AddEntry(gCopy, Form("a_{0} + #Sigma_{n=1}^{%i}2a_{n}cos(n#Delta#phi)", par.size()-1), "l");
    //leg->AddEntry(gPeri, "Baseline a0", "l");
    //leg->AddEntry(fit_p1, Form("a_{0} + 2a_{1}cos(#Delta#phi), v1^{2}#times10^{3} = %0.2f", v12*1e3), "l");
    //leg->AddEntry(fit_p2, Form("a_{0} + 2a_{2}cos(2#Delta#phi), v2^{2}#times10^{3} = %0.2f", v21*1e3), "l");
    //leg->AddEntry(fit_p3, Form("a_{0} + 2a_{3}cos(3#Delta#phi), v3^{2}#times10^{3} = %0.2f", v31*1e3), "l");
    //leg->AddEntry(fit_p4, Form("a_{0} + 2a_{4}cos(4#Delta#phi), v4^{2}#times10^{3} = %0.2f", v41*1e3), "l");
    leg->Draw();

    // =============== 新增：添加chi2/ndf标签 ===============
    TLatex* chi2Label = new TLatex();
    chi2Label->SetNDC();
    chi2Label->SetTextFont(43);
    chi2Label->SetTextSize(20);
    chi2Label->DrawLatex(0.50, 0.70, Form("#chi^{2}/ndf = %.1f/%d = %.2f", chi2, ndf, chi2ndf));

    // 添加文本标签
    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(43);
    tex->SetTextSize(20);
    DrawText(0.2, 0.85, 20, Form("ALICE %s", collisionSystemName.c_str()));
    DrawText(0.2, 0.80, 20, Form("%s", FITused.c_str()));
    if (isNch) {
        DrawText(0.2, 0.75, 20, Form("%d < N_{ch} < %d", minRange, maxRange));
    }
    else {
        DrawText(0.2, 0.75, 20, Form("%d < Cent < %d", minRange, maxRange));
    }

    if (etaMin > -100. && etaMax > -100.)
        DrawText(0.2, 0.70, 20, Form("%.1f < #eta^{TPC} < %.1f", etaMin, etaMax));
    DrawText(0.2, 0.65, 20, Form("V_{2#Delta} = %.5f #pm %.5f", v21, v21e));
    DrawText(0.2, 0.60, 20, Form("V_{3#Delta} = %.5f #pm %.5f", v31, v31e));

    // 绘制底部残差图
    pad2->cd();


    // 创建残差直方图
    TH1D* hsubtract = (TH1D*)hm->Clone(Form("subtract"));
    // hsubtract->Add(gPeri, -1);

    TH1D* hResidual = (TH1D*)hsubtract->Clone(Form("pull"));
    double ymax_pull = 1, ymin_pull = 1;
    for (int i = 1; i <= hResidual->GetXaxis()->GetNbins(); i++) {
      double dat = hResidual->GetBinContent(i);
      double err = hResidual->GetBinError(i);
      double lin = fittedTemp->GetBinContent(i);
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
    if (etaMin > -100. && etaMax > -100.) {
        canvas->SaveAs(Form("./TemplateFit/EtaDiff/PDFs/TemplateFourierFit_%s_%s_%d_%d_Eta_%0.1f_%0.1f.pdf", fileSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax));
    } else {
        canvas->SaveAs(Form("./TemplateFit/PDFs/TemplateFourierFit_%s_%s_%d_%d.pdf", fileSuffix.c_str(), splitName.c_str(), minRange, maxRange));
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
    //Initialize fitter with given projections
    TemplateFitter *ft = new TemplateFitter(hm);
    //Setting up variable ( = delta phi, or just "x"):
    ft->AddVariable("x", "x", -TMath::Pi()/2.0, 1.5*TMath::Pi());

    if (!kRefit){ //dersom fit failer, skal man forandre disse: (if (!kRefit){)
        // Pb-Pb initial value
        ft->AddParameter("Fa","Fa",5., 0., 300.); //ft->AddParameter("Fa","Fa",4.5,0,30);
        ft->AddParameter("Ga","Ga",30000.0,0.,3000000000.); //ft->AddParameter("Ga","Ga",15,0,1000);
        ft->AddParameter("v2","v2",5e-3,-1.0,1.0); //ft->AddParameter("v2","v2",4e-3,-1.0,1.0); //kanskje større initial her?
        ft->AddParameter("v3","v3",1e-3,-1.0,1.0); //ft->AddParameter("v3","v3",6e-4,-1.0,1.0);
        ft->AddParameter("v4","v4",1e-4,-1.0,1.0);
        ft->AddParameter("v1","v1", 0., 0., 0.); // (5e-2, -1.0, 1.0)


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
    Double_t v11 = ft->getVal(5);
    Double_t v11e = ft->getErr(5);
    Double_t chi2ndf = findChi2ndf(lm, hm, F, G, v11, v21, v31, v41);

    printf("Values from fit:\n");
    printf("F  = %f +- %f\n",F,Fe);
    printf("G  = %f +- %f\n",G,Ge); //skrives ut i terminalen
    printf("V1  = %f +- %f\n",v11,v11e);
    printf("V2 = %f +- %f\n",v21,v21e);
    printf("V3 = %f +- %f\n",v31,v31e);
    printf("V4 = %f +- %f\n",v41,v41e);

    fParamVal[0] = v21; fParamErr[0] = v21e;
    fParamVal[1] = v31; fParamErr[1] = v31e;
    fParamVal[2] = v41; fParamErr[2] = v41e;
    fParamVal[3] = F;   fParamErr[3] = Fe;
    fParamVal[4] = G;   fParamErr[4] = Ge; //4
    fParamVal[5] = v11; fParamErr[5] = v11e;
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

void PlotFitting(TH1 *lm, TH1 *hm, Bool_t isNch, std::string fileSuffix, Int_t minRange, Int_t maxRange, const std::vector<Double_t>& par, const std::vector<Double_t>& parerr, Double_t etaMin=-101., Double_t etaMax=-101.) {
    gStyle->SetOptStat(0); 
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    
    double v21 = par[1];
    double v31 = par[2];
    double v41 = par[3];
    double F =   par[3];
    double G =   par[4];
    double v11 = par[5];

    double v21e = parerr[0];
    double v31e = parerr[1];
    double v41e = parerr[2];
    double Fe = parerr[3];
    double Ge = parerr[4];
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

    // Generate fit curves
    const int pointBin = (int)hm->GetNbinsX();
    Double_t CopyPointX[pointBin];
    Double_t CopyPointY[pointBin];
    for (int i=0; i<pointBin; ++i){
      CopyPointX[i] = hm->GetBinCenter(i+1);
      double x = hm->GetBinCenter(i+1);
      // CopyPointY[i] = F*lm->GetBinContent(i+1)+G*(1+2*v21*cos(2*x)+2*v31*cos(3*x));
      CopyPointY[i] = F*lm->GetBinContent(i+1) + G*(1+2*v11*cos(x)+2*v21*cos(2*x)+2*v31*cos(3*x)+2*v41*cos(4*x));
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

    TF1* fit_p5 = new TF1("fit_p5","[0]*[1] + [2]*(1 + 2*[3]*cos(x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    TF1* fit_p2 = new TF1("fit_p2","[0]*[1] + [2]*(1 + 2*[3]*cos(2*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    TF1* fit_p3 = new TF1("fit_p3","[0]*[1] + [2]*(1 + 2*[3]*cos(3*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    TF1* fit_p4 = new TF1("fit_p4","[0]*[1] + [2]*(1 + 2*[3]*cos(4*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    fit_p5->SetParameters(F,Y0,G,v11);
    fit_p2->SetParameters(F,Y0,G,v21);
    fit_p3->SetParameters(F,Y0,G,v31);
    fit_p4->SetParameters(F,Y0,G,v41);
    fit_p5->SetLineColor(colors[5]);
    fit_p5->SetLineWidth(2);
    fit_p5->SetLineStyle(2);
    fit_p2->SetLineColor(colors[2]);
    fit_p2->SetLineWidth(2);
    fit_p2->SetLineStyle(2);
    fit_p3->SetLineColor(colors[3]);
    fit_p3->SetLineWidth(2);
    fit_p3->SetLineStyle(3);
    fit_p4->SetLineColor(colors[4]);
    fit_p4->SetLineWidth(2);
    fit_p4->SetLineStyle(4);
    fit_p5->Draw("same");
    fit_p2->Draw("same");
    fit_p3->Draw("same");
    fit_p4->Draw("same");

    // 添加图例
    TLegend* leg = new TLegend(0.5, 0.65, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->AddEntry(hm, "Data", "lep");
    leg->AddEntry(gCopy, "FY(#Delta#phi)^{peri} + G(1+#Sigma_{n=1}^{4}2V_{n#Delta}cos(n#Delta#phi))", "l");
    leg->AddEntry(gPeri, "FY(#Delta#phi)^{peri} + G", "l");
    leg->AddEntry(fit_p5, "FY(0)^{peri} + G(1+2V_{1#Delta}cos(#Delta#phi))", "l");
    leg->AddEntry(fit_p2, "FY(0)^{peri} + G(1+2V_{2#Delta}cos(2#Delta#phi))", "l");
    leg->AddEntry(fit_p3, "FY(0)^{peri} + G(1+2V_{3#Delta}cos(3#Delta#phi))", "l");
    leg->AddEntry(fit_p4, "FY(0)^{peri} + G(1+2V_{4#Delta}cos(4#Delta#phi))", "l");
    leg->Draw();

    // 添加文本标签
    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(43);
    tex->SetTextSize(20);
    DrawText(0.2, 0.85, 20, Form("ALICE %s", collisionSystemName.c_str()));
    DrawText(0.2, 0.80, 20, Form("%s", FITused.c_str()));
    if (isNch) {
        DrawText(0.2, 0.75, 20, Form("%d < N_{ch} < %d", minRange, maxRange));
    }
    else {
        DrawText(0.2, 0.75, 20, Form("%d < Cent < %d", minRange, maxRange));
    }

    if (etaMin > -100. && etaMax > -100.)
        DrawText(0.2, 0.70, 20, Form("%.1f < #eta^{TPC} < %.1f", etaMin, etaMax));
    DrawText(0.2, 0.65, 20, Form("F = %.5f #pm %.5f", F, Fe));
    DrawText(0.2, 0.60, 20, Form("V_{1#Delta} = %.5f #pm %.5f", v11, v11e));
    DrawText(0.2, 0.55, 20, Form("V_{2#Delta} = %.5f #pm %.5f", v21, v21e));
    //DrawText(0.2, 0.60, 20, Form("V_{3#Delta} = %.5f #pm %.5f", v31, v31e));

    // 绘制底部残差图
    pad2->cd();

    TF1* fit_p1m = new TF1("fit_p1m","[0]*(1 + 2*[4]*cos(x) + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    fit_p1m->SetParameters(G,v21,v31,v41,v11);

    
    // here we subtract lm from hm
    TH1D* hsubtract = (TH1D*)hm->Clone(Form("subtract"));
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

    //create the vector that we will write into csv
    std::vector<std::vector<std::string>> csvData;
    //fill the vector with data
    csvData.push_back({std::to_string(etaMin), std::to_string(etaMax), std::to_string(F), std::to_string(Fe), std::to_string(G), std::to_string(Ge), std::to_string(v11), std::to_string(v11e), std::to_string(v21), std::to_string(v21e), std::to_string(v31), std::to_string(v31e), std::to_string(v41), std::to_string(v41e), std::to_string(chi2ndf)});

    // 保存结果
    if (etaMin > -100. && etaMax > -100.) {
        canvas->SaveAs(Form("./TemplateFit/etaDiff/PDFs/TemplateFit_%s_%s_%d_%d_Eta_%0.1f_%0.1f.pdf", fileSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax)); //her lageres diff .pdf
        writeToCSV("TemplateFit.csv", csvData, "./TemplateFit/etaDiff/PDFs");
    } else {
        canvas->SaveAs(Form("./TemplateFit/PDFs/TemplateFit_%s_%s_%d_%d.pdf", fileSuffix.c_str(), splitName.c_str(), minRange, maxRange)); //det samme, bare uten ikke diff
        writeToCSV("TemplateFit.csv", csvData, "./TemplateFit/PDFs");
    }
    
    // 清理内存
    delete hbkg1;
    delete gCopy;
    delete gPeri;
    delete fit_p2;
    delete fit_p3;
    delete fit_p4;
    delete fit_p5;
    delete leg;
    delete tex;
    delete hResidual;
    delete line;
  
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
        header.push_back({"etaMin", "etaMax", "F", "Fe", "G", "Ge", "v1", "v1e", "v2", "v2e", "v3", "v3e", "v4", "v4e", "chi2ndf", "chi2ndfe"});
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




double findChi2ndf(TH1 *lm, TH1 *hm, double& F_val, double& G_val, double& v11_val, double& v21_val, double& v31_val, double& v41_val) {
    double chi2;
    int nBins;
    int nParams;
    int ndf;
    double data;
    double error;
    double x;
    double fit;
    double residual;
    double chi2ndf;

    TF1* fit_p1m = new TF1("fit_p1m","[0]*(1 + 2*[4]*cos(x) + 2*[1]*cos(2*x) + 2*[2]*cos(3*x) + 2*[3]*cos(4*x))", -TMath::Pi()/2.0, 1.5*TMath::Pi());
    fit_p1m->SetParameters(G_val,v21_val,v31_val,v41_val,v11_val);

    TH1D* hsubtract = (TH1D*)hm->Clone(Form("subtract"));
    hsubtract->Add(lm, -F_val);

    chi2 = 0.0;
    nBins = hsubtract->GetNbinsX();
    nParams = 6; // 参数个数: F,G,v11,v21,v31,v41 !!!Needs to be changed for FT
    ndf = nBins - nParams; //degrees of freedom (the number of data points minus number of fit parameters)
    
    for (int i = 1; i <= nBins; i++) {
        data = hsubtract->GetBinContent(i);
        error = hsubtract->GetBinError(i);
        x = hsubtract->GetBinCenter(i);
        fit = fit_p1m->Eval(x);
        
        if (error > 0) { // 忽略误差为0的bin
            residual = data - fit;
            chi2 += (residual * residual) / (error * error); //Yes chi is waighted by the errors, i.e. large error low chi
        }
    }
    chi2ndf = (ndf > 0) ? chi2 / ndf : 0;

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