#ifndef BASICFORDIHADRON_H
#define BASICFORDIHADRON_H

#include "TH1D.h"
#include <complex>
#include <vector>
#include <cmath>


enum {corrAxis_kSample, corrAxis_kVz, corrAxis_kPt_TPC_trig, corrAxis_kPt_TPC_asso, corrAxis_kdPhiTPCTPC, corrAxis_kdEtaTPCTPC};
enum {trigAxis_sample, trigAxis_Vz, trigAxis_pT};

//Associated particle pT range
Int_t maxSample = 10; //Number of histograms we create, i.e. the number of times we loop over the processing functions (ref/diff) in dptdphi

//trigger particle pT bins
std::vector<float> etaTPC = {-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}; //trigger of TPC etas

// FT0 detector eta slices (4 slices per detector for original 598682 dataset)
std::vector<float> etaFT0C = {-2.2, -2.4, -2.6, -2.8, -3}; // FT0C: 4 slices
std::vector<float> etaFT0A = {3.5, 3.8, 4.1, 4.5, 4.9};      // FT0A: 4 slices

// FT0 detector full range (no slicing for ring datasets 616549/616551)
std::vector<float> etaFT0C_NoSlice = {-3.0, -2.1};  // FT0C full range: one "slice"
std::vector<float> etaFT0A_NoSlice = {3.5, 4.9};    // FT0A full range: one "slice"

double minPt = 0.2;
double maxPt = 10.0;

Int_t numberOfPlots = 6;


enum kDihadronCorrType {
    kTPCFT0A,
    kTPCFT0C,
    kFT0AFT0C,
    kNCorrType
};

std::map<int, std::string> DihadronCorrTypeName = {
    {kTPCFT0A, "TPC_FT0A"},
    {kTPCFT0C, "TPC_FT0C"},
};

std::map<int, std::vector<float>> DihadrondEtaRange = {
    {kTPCFT0A,{-5.8, -2.6}},
    {kTPCFT0C, {1.2, 4.2}},
};

std::map<int, std::vector<float>> FitEtaCoverage = {
    {kTPCFT0A,{3.5,4.9}},
    {kTPCFT0C, {-3.3,-2.1}},
};


enum {
    kCent = 0,
    kNch = 1
};

enum {
    kData = 0,
    kMC = 1
};

enum {
    kEtaDiffOff = 0,
    kEtaDiffOn = 1
};

enum kDihadronMethod{
    kFourierFit,
    kTemplateFit,
    kImprovedTemplateFit,
    kPeripheralSubFit,
    kNMethod
};

std::map<int, std::string> DihadronMethodName = {
    {kFourierFit, "FourierFit"},
    {kTemplateFit, "TemplateFit"},
    {kImprovedTemplateFit, "ImprovedTemplateFit"},
    {kPeripheralSubFit, "PeripheralSubFit"}
};

double BarlowRequirement = 1.;
double BarlowNbinsFraction = 0.01;
double MaxpT = 10.;

enum kObservable{
    kV22,
    kV32,
    kV42,
    kpTDiffv22,
    kpTDiffv32,
    kpTDiffv42,
    kNObservable
};

std::map<int, std::string> ObservableFilesMap = {
    {kV22, "Vn"},
    {kV32, "Vn"},
    {kV42, "Vn"},
    {kpTDiffv22, "PtDiff/Vn"},
    {kpTDiffv32, "PtDiff/Vn"},
    {kpTDiffv42, "PtDiff/Vn"}
};

std::map<int, std::vector<std::string>> ObservableNamesMap = {
    {kV22, {"hV2"}},
    {kV32, {"hV3"}},
    {kV42, {"hV4"}},
    {kpTDiffv22, {"hV2"}},
    {kpTDiffv32, {"hV3"}},
    {kpTDiffv42, {"hV4"}}
};

std::map<int, std::vector<std::string>> ObservableOutputNamesMap = {
    {kV22, {"v22"}},
    {kV32, {"v32"}},
    {kV42, {"v42"}},
    {kpTDiffv22, {"pTDiffv22"}},
    {kpTDiffv32, {"pTDiffv32"}},
    {kpTDiffv42, {"pTDiffv42"}}
};

std::map<int, std::vector<std::string>> ObservablePrintNamesMap = {
    {kV22, {"v_2\\{2\\}"}},
    {kV32, {"v_3\\{2\\}"}},
    {kV42, {"v_4\\{2\\}"}},
    {kpTDiffv22, {"v_2\\{2\\}(p_{T})"}},
    {kpTDiffv32, {"v_3\\{2\\}(p_{T})"}},
    {kpTDiffv42, {"v_4\\{2\\}(p_{T})"}}
};

void HistFFT(TH1* hist, std::vector<double>& Coeff) {
    if (!hist) return;
    
    Coeff.resize(6, 0.0); // 存储常数项(a0)和n=1到5的系数
    
    const double pi = TMath::Pi();
    const int nbins = hist->GetNbinsX();
    const double xmin = -0.5*pi;
    const double xmax = 1.5*pi;
    const double bin_width = (xmax - xmin) / nbins;
    
    // 计算归一化因子
    double norm = 0.0;
    for (int i = 1; i <= nbins; ++i) {
        norm += hist->GetBinContent(i) * bin_width;
    }
    
    // 计算常数项(a0)
    double a0 = 0.0;
    for (int i = 1; i <= nbins; ++i) {
        double x = xmin + (i - 0.5) * bin_width;
        a0 += hist->GetBinContent(i) * bin_width;
    }
    Coeff[0] = a0 / (2*pi);
    
    // 计算n=1到5的傅里叶系数
    for (int n = 1; n <= 5; ++n) {
        std::complex<double> sum(0, 0);
        for (int i = 1; i <= nbins; ++i) {
            double x = xmin + (i - 0.5) * bin_width;
            double y = hist->GetBinContent(i);
            std::complex<double> exp_term(cos(n*x), -sin(n*x));
            sum += y * bin_width * exp_term;
        }
        Coeff[n] = std::abs(sum) / pi;
    }
}    


#endif // BASICFORDIHADRON_H
