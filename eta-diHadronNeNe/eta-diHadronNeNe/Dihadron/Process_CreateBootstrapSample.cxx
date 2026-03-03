/*
 * @Author: Zhiyong Lu (zhiyong.lu@cern.ch)  
 * @Date: 2025-05-18 17:23:23 
 * @Last Modified by: Zhiyong Lu
 * @Last Modified time: 2025-05-18 18:12:39
 */
//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "TList.h"
#include "TH2D.h"
//#include "TRandom3.h"
#include <iostream>
#include <string>
#include <vector>
#include "./include/BasicForDihadron.h"
#include <fstream>
#include <sstream>
#include <algorithm>



struct InputUnit {
    std::string fileNameSuffix;
    Int_t corrType;
    Bool_t isNch;
    Bool_t isEtadiff;
    Int_t minRange;
    Int_t maxRange;

    InputUnit(std::string _fileNameSuffix, Int_t _corrType, Bool_t _isNch, Bool_t _isEtadiff, Int_t _minRange, Int_t _maxRange) :
        fileNameSuffix(_fileNameSuffix), corrType(_corrType), isNch(_isNch), isEtadiff(_isEtadiff), minRange(_minRange), maxRange(_maxRange) {}
};

void CreateBootstrapSample(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.);
void CreateBootstrapSample_EtaDiff(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax, Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.);
void readFakeFull(std::vector<int> &histoIndices, std::vector<double> &histoWaights, int maxSample, const std::string &fileName);

void Process_CreateBootstrapSample() {

    std::vector<InputUnit> inputList;

    // Simple switch: run one detector at a time to avoid overwriting
    // Change to kTPCFT0A for FT0A, or kTPCFT0C for FT0C
    const Int_t corrTypeToRun = kTPCFT0C;

    //Ne-Ne

    // inputList.push_back(InputUnit("LHC25af_pass2_598673", kTPCFT0C, kCent, kEtaDiffOff, 0, 20)); 
    // inputList.push_back(InputUnit("LHC25af_pass2_598673", kTPCFT0C, kCent, kEtaDiffOff, 80, 100)); 

    // inputList.push_back(InputUnit("LHC25af_pass2_598673", kTPCFT0A, kCent, kEtaDiffOff, 0, 20)); 
    // inputList.push_back(InputUnit("LHC25af_pass2_598673", kTPCFT0A, kCent, kEtaDiffOff, 80, 100)); 

    // inputList.push_back(InputUnit("LHC25af_pass2_598673", kTPCFT0C, kCent, kEtaDiffOn, 0, 20)); 
    // inputList.push_back(InputUnit("LHC25af_pass2_598673", kTPCFT0C, kCent, kEtaDiffOn, 80, 100)); 

    // inputList.push_back(InputUnit("LHC25af_pass2_598673", kTPCFT0A, kCent, kEtaDiffOn, 0, 20)); 
    // inputList.push_back(InputUnit("LHC25af_pass2_598673", kTPCFT0A, kCent, kEtaDiffOn, 80, 100)); 

    //O-O

    // inputList.push_back(InputUnit("LHC25ae_pass2_598682", kTPCFT0C, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_598682", kTPCFT0C, kCent, kEtaDiffOff, 80, 100));

    // inputList.push_back(InputUnit("LHC25ae_pass2_598682", kTPCFT0A, kCent, kEtaDiffOff, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_598682", kTPCFT0A, kCent, kEtaDiffOff, 80, 100));

    // inputList.push_back(InputUnit("LHC25ae_pass2_598682", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_598682", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));

    // inputList.push_back(InputUnit("LHC25ae_pass2_598682", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    // inputList.push_back(InputUnit("LHC25ae_pass2_598682", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));

    // Ne-Ne outer ring (615818) - full outer ring, correlate with TPC edges
    inputList.push_back(InputUnit("LHC25af_pass2_615818", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_615818", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));

    inputList.push_back(InputUnit("LHC25af_pass2_615818", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_615818", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));

    // Ne-Ne inner ring (615817) - full inner ring, correlate with TPC edges
    inputList.push_back(InputUnit("LHC25af_pass2_615817", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_615817", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));

    inputList.push_back(InputUnit("LHC25af_pass2_615817", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_615817", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));

    // O-O outer ring (616549) - full outer ring, correlate with TPC edges
    inputList.push_back(InputUnit("LHC25ae_pass2_616549", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25ae_pass2_616549", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));

    inputList.push_back(InputUnit("LHC25ae_pass2_616549", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25ae_pass2_616549", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));

    // O-O inner ring (618685) - full inner ring, correlate with TPC edges
    inputList.push_back(InputUnit("LHC25ae_pass2_618685", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25ae_pass2_618685", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));

    inputList.push_back(InputUnit("LHC25ae_pass2_618685", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25ae_pass2_618685", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));

    // Ne-Ne inner ring (617826) - new dataset
    inputList.push_back(InputUnit("LHC25af_pass2_617826", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_617826", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));

    inputList.push_back(InputUnit("LHC25af_pass2_617826", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_617826", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));

    // Ne-Ne outer ring (617910) - new dataset
    inputList.push_back(InputUnit("LHC25af_pass2_617910", kTPCFT0C, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_617910", kTPCFT0C, kCent, kEtaDiffOn, 80, 100));

    inputList.push_back(InputUnit("LHC25af_pass2_617910", kTPCFT0A, kCent, kEtaDiffOn, 0, 20));
    inputList.push_back(InputUnit("LHC25af_pass2_617910", kTPCFT0A, kCent, kEtaDiffOn, 80, 100));

    for (auto input : inputList) {
        if (input.corrType != corrTypeToRun) continue;

        // Check if this is a ring dataset
        bool isRingDataset = (input.fileNameSuffix.find("615817") != std::string::npos || 
                              input.fileNameSuffix.find("615818") != std::string::npos ||
                              input.fileNameSuffix.find("616549") != std::string::npos ||
                              input.fileNameSuffix.find("618685") != std::string::npos ||
                              input.fileNameSuffix.find("617826") != std::string::npos ||
                              input.fileNameSuffix.find("617910") != std::string::npos);

        // For ring datasets: process with TPC edges, no FT0 slicing
        if (isRingDataset && input.isEtadiff) {
            std::vector<std::pair<double,double>> tpcEdges = { { -0.8, -0.7 }, { 0.7, 0.8 } };
            bool isInner = (input.fileNameSuffix.find("615817") != std::string::npos || 
                            input.fileNameSuffix.find("618685") != std::string::npos ||
                            input.fileNameSuffix.find("617826") != std::string::npos);
            std::string ringType = isInner ? "inner" : "outer";
            // For ring datasets, use whole FT0 detector without subdivision
            for (auto rng : tpcEdges) {
                double etaMin = rng.first;
                double etaMax = rng.second;
                std::cout << "Processing Bootstrap Sample ring " << ringType << " [" << etaMin << "," << etaMax << "]: " << input.fileNameSuffix << std::endl;
                CreateBootstrapSample_EtaDiff(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange, etaMin, etaMax);
            }
        }
        // For original FT0-sliced datasets (598682): use FT0 slices
        else {
            // Get the appropriate FT0 eta slices for this detector type
            const std::vector<float>* ft0Slices = nullptr;
            if (input.corrType == kTPCFT0C) {
                ft0Slices = &etaFT0C;
            } else if (input.corrType == kTPCFT0A) {
                ft0Slices = &etaFT0A;
            }

            // Loop over FT0 slices
            if (ft0Slices) {
                for (int iFT0 = 0; iFT0 < ft0Slices->size() - 1; iFT0++) {
                    double ft0EtaMin = (*ft0Slices)[iFT0];
                    double ft0EtaMax = (*ft0Slices)[iFT0 + 1];

                    if (input.isEtadiff) {
                        std::cout << "Processing Bootstrap Sample eta diff FT0[" << ft0EtaMin << "," << ft0EtaMax << "]: " << input.fileNameSuffix << std::endl;
                        // For FT0 slicing, process the wider TPC ranges: [-0.8, -0.6] and [0.6, 0.8]
                        std::vector<std::pair<double,double>> tpcRanges = { { -0.8, -0.6 }, { 0.6, 0.8 } };
                        for (auto rng : tpcRanges) {
                            double etaMin = rng.first;
                            double etaMax = rng.second;
                            CreateBootstrapSample_EtaDiff(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange, etaMin, etaMax, ft0EtaMin, ft0EtaMax);
                        }
                    } else {
                        std::cout << "Processing Bootstrap Sample FT0[" << ft0EtaMin << "," << ft0EtaMax << "]: " << input.fileNameSuffix << std::endl;
                        CreateBootstrapSample(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange, ft0EtaMin, ft0EtaMax);
                    }
                }
            } else {
                // No FT0 slicing for other correlation types
                if (input.isEtadiff) {
                    std::cout << "Processing Bootstrap Sample eta diff: " << input.fileNameSuffix << std::endl;
                    for (int iEta = 0; iEta < etaTPC.size() - 1; iEta++) {
                        double etaMin = etaTPC[iEta];
                        double etaMax = etaTPC[iEta + 1];
                        CreateBootstrapSample_EtaDiff(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange, etaMin, etaMax);
                    }
                } else {
                    std::cout << "Processing Bootstrap Sample: " << input.fileNameSuffix << std::endl;
                    CreateBootstrapSample(input.fileNameSuffix, input.corrType, input.isNch, input.minRange, input.maxRange);
                }
            }
        }
    }

}

void CreateBootstrapSample(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    
    //Load values from mixed
    const char* inPath1;
    if (ft0EtaMin > -900. && ft0EtaMax > -900.) {
        inPath1 = Form("./ProcessOutput/Mixed_%s_%s_%i_%i_%s_FT0Eta_%0.1f_%0.1f.root",
                       fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, 
                       DihadronCorrTypeName[corrType].c_str(), ft0EtaMin, ft0EtaMax);
    } else {
        inPath1 = Form("./ProcessOutput/Mixed_%s_%s_%i_%i_%s.root",
                       fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str());
    }
    const char* inPath2 = Form("./Dihadron/ProcessOutput/Mixed_%s_%s_%i_%i_%s.root",
                               fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, DihadronCorrTypeName[corrType].c_str());
    TFile* file = TFile::Open(inPath1, "READ");
    if (!file || file->IsZombie()) {
        if (file) { file->Close(); delete file; }
        file = TFile::Open(inPath2, "READ");
    }

    if (!file || file->IsZombie()) { //error
        std::cerr << "Error opening input file! Tried: " << inPath1 << " and " << inPath2 << std::endl;
        return;
    }

    // 读取所有样本直方图
    std::vector<TH1D*> hists;
    for (Int_t sample = 0; sample < maxSample; ++sample) { //10 samples
        TH1D* h = dynamic_cast<TH1D*>(
            file->Get(Form("hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, sample)) //samples collected
        );
        if (!h) { //error
            std::cerr << "Error loading histogram for sample " << sample << std::endl;
            continue;
        }
        hists.push_back(h);
    }

    if (hists.empty()) { //error
        std::cerr << "No histograms loaded!" << std::endl;
        file->Close();
        return;
    }

    //Make bootstrap output
    TString outFileName;
    if (ft0EtaMin > -900. && ft0EtaMax > -900.) {
        outFileName = Form("./ProcessOutput/BootstrapSample_%s_%s_%i_%i_FT0Eta_%0.1f_%0.1f.root",
                          fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange,
                          ft0EtaMin, ft0EtaMax);
    } else {
        outFileName = Form("./ProcessOutput/BootstrapSample_%s_%s_%i_%i.root",
                          fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange);
    }
    TFile* outFile = TFile::Open(outFileName, "RECREATE");

    if (!outFile || outFile->IsZombie()) { //error
        std::cerr << "Error creating output file!" << std::endl;
        file->Close();
        return;
    }

    TH1D* hAll = dynamic_cast<TH1D*>(file->Get(Form("hPhiSameOverMixed_%d_%d", minRange, maxRange)));
    if (!hAll) {
        std::cerr << "Error loading all-sample histogram!" << std::endl;
        return;
    }
    hAll->SetName(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange));
    hAll->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange));
    hAll->GetXaxis()->SetTitle("#Delta#varphi");
    hAll->Write();

    // 初始化随机数生成器

    //Mine:
    std::vector<int> histoIndices;
    std::vector<double> histoWaights;
    std::string fileName = "FakeFullSamples/fake_full_sample.csv"; //form suffix for EtaDiff

    readFakeFull(histoIndices, histoWaights, maxSample, fileName);
    //End mine


    // 生成 maxSample^2 个bootstrap样本
    for (Int_t bs = 0; bs < maxSample * maxSample; ++bs) { //The 10 sub-samples 10 times... EDIT: now 10^2 = 100 samples
        // 随机选择样本索引（允许重复）
        std::vector<Int_t> selectedIndices;
        selectedIndices.clear();
        for (Int_t i = 0; i < numberOfPlots; ++i) { // randomly select samples with replacement
            int histoIndex = bs * 10 + i; //mine
            selectedIndices.push_back(histoIndices[histoIndex]); //pick 10 random integers EDIT MINE!!!
            //selectedIndices.push_back(i);
        }

        //Merge histograms
        TH1D* hmerge = nullptr;

        double totalWeight = 0.0; //?

        for (Int_t i=0; i < selectedIndices.size(); ++i) { //loop over the 10 selected indices

            TH1D* current = hists[selectedIndices[i]]; //load the histogram (sample) corresponding to the random index
            Double_t weight = current->GetEntries(); //So we treat each histogram same (commented out later)
            int histoIndex = bs * 10 + i; //mine]

            if (!hmerge) { //First time initializing
                hmerge = dynamic_cast<TH1D*>(current->Clone( //clone the current
                    Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs) 
                ));
                //hmerge->Scale(histoWaights[histoIndex]); //her brukte jeg vekter i tilleg til å velge tilfeldige histograms
                hmerge->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs));
                // hmerge->Scale(weight);
                // std::cout << "bs: " << bs << " idx: " << idx << " weight: " << weight << std::endl;
                totalWeight = histoWaights[histoIndex]; // * weight;
            } else {       // after the first initialization, we add
                // current->Scale(weight);
                // hmerge->Add(current, weight);
                hmerge->Add(current); //MED VEKTER
                //hmerge->Add(current);
                totalWeight += histoWaights[histoIndex]; // * weight;
            }
        }

        // 归一化并保存
        if (hmerge && totalWeight > 0) {
            //hmerge->Scale(1.0 / totalWeight);
            hmerge->Scale(1.0 / selectedIndices.size()); //divide by numberOfPlots...
            // Printf("size: %d", selectedIndices.size());
            // hmerge->SetDirectory(outFile);
            hmerge->Write();
            // delete hmerge; // 释放内存
        }
    }

    // 清理资源
    outFile->Close();
    file->Close();
    delete outFile;
    delete file;
}

/*
diHadron.cxx spits 10 samples (i.e. 10 correlation with smaller number of collisions)
dPhidEta creates 10 samples 1D dEta projection
bootstrap picks 10 random indices
creates a new histogram by mergining 10 samples based on the random indices
does this 100 times, i.e. we use 100 elemenets in total
*/

void CreateBootstrapSample_EtaDiff(std::string fileNameSuffix, Int_t corrType, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax, Double_t ft0EtaMin=-999., Double_t ft0EtaMax=-999.) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";

    const char* inPath1;
    const char* inPath2;
    if (ft0EtaMin > -900. && ft0EtaMax > -900.) {
        inPath1 = Form("./ProcessOutput/EtaDiff/Mixed_%s_%s_%i_%i_Eta_%0.1f_%0.1f_%s_FT0Eta_%0.1f_%0.1f.root",
                       fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax, 
                       DihadronCorrTypeName[corrType].c_str(), ft0EtaMin, ft0EtaMax);
        inPath2 = Form("./Dihadron/ProcessOutput/EtaDiff/Mixed_%s_%s_%i_%i_Eta_%0.1f_%0.1f_%s_FT0Eta_%0.1f_%0.1f.root",
                       fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax, 
                       DihadronCorrTypeName[corrType].c_str(), ft0EtaMin, ft0EtaMax);
    } else {
        inPath1 = Form("./ProcessOutput/EtaDiff/Mixed_%s_%s_%i_%i_Eta_%0.1f_%0.1f_%s.root",
                       fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax, DihadronCorrTypeName[corrType].c_str());
        inPath2 = Form("./Dihadron/ProcessOutput/EtaDiff/Mixed_%s_%s_%i_%i_Eta_%0.1f_%0.1f_%s.root",
                       fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax, DihadronCorrTypeName[corrType].c_str());
    }
    TFile* file = TFile::Open(inPath1, "READ");
    if (!file || file->IsZombie()) {
        if (file) { file->Close(); delete file; }
        file = TFile::Open(inPath2, "READ");
    }
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening input file! Tried: " << inPath1 << " and " << inPath2 << std::endl;
        return;
    }

    // 读取所有样本直方图
    std::vector<TH1D*> hists;
    for (Int_t sample = 0; sample < maxSample; ++sample) {
        TH1D* h = dynamic_cast<TH1D*>(
            file->Get(Form("hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, sample))
        );
        if (!h) {
            std::cerr << "Error loading histogram for sample " << sample << std::endl;
            continue;
        }
        hists.push_back(h);
    }

    if (hists.empty()) {
        std::cerr << "No histograms loaded!" << std::endl;
        file->Close();
        return;
    }

    TString outFileName;
    if (ft0EtaMin > -900. && ft0EtaMax > -900.) {
        outFileName = Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%i_%i_Eta_%0.1f_%0.1f_FT0Eta_%0.1f_%0.1f.root",
                          fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, 
                          etaMin, etaMax, ft0EtaMin, ft0EtaMax);
    } else {
        outFileName = Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%i_%i_Eta_%0.1f_%0.1f.root",
                          fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax);
    }
    TFile* outFile = TFile::Open(outFileName, "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error creating output file!" << std::endl;
        file->Close();
        return;
    }

    TH1D* hAll = dynamic_cast<TH1D*>(file->Get(Form("hPhiSameOverMixed_%d_%d", minRange, maxRange)));
    if (!hAll) {
        std::cerr << "Error loading all-sample histogram!" << std::endl;
        return;
    }
    hAll->SetName(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange));
    hAll->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d", minRange, maxRange));
    hAll->GetXaxis()->SetTitle("#Delta#varphi");
    hAll->Write();

    // 初始化随机数生成器
    //TRandom3 randGen;

    //Mine:
    std::vector<int> histoIndices;
    std::vector<double> histoWaights;
    int etaBin = std::upper_bound(etaTPC.begin(), etaTPC.end(), etaMin) - etaTPC.begin() - 1;
    if (etaBin > 14) etaBin = 14; // avoid missing fake_full_sample_15.csv
    std::string fileName = Form("FakeFullSamples/EtaDiff/fake_full_sample_%i.csv", etaBin); //form suffix for EtaDiff

    readFakeFull(histoIndices, histoWaights, maxSample, fileName);
    if (histoIndices.empty() || histoWaights.empty()) {
        std::cerr << "Error: no bootstrap weights loaded from '" << fileName << "'. Skipping eta-diff bootstrap." << std::endl;
        outFile->Close();
        file->Close();
        delete outFile;
        delete file;
        return;
    }
    //End mine


    // 生成 maxSample^2 个bootstrap样本
    for (Int_t bs = 0; bs < maxSample * maxSample; ++bs) { //EDIT from 10^2 to 10^3
        // 随机选择样本索引（允许重复）
        std::vector<Int_t> selectedIndices;
        selectedIndices.clear();
        for (Int_t i = 0; i < numberOfPlots; ++i) {
            int histoIndex = bs * 10 + i; //mine
            selectedIndices.push_back(histoIndices[histoIndex]); //pick 10 random integers EDIT !!!MINE
            //selectedIndices.push_back(i);
        }

        // 合并选中的直方图
        TH1D* hmerge = nullptr;
        double totalWeight = 0.0;

        for (Int_t i=0; i < selectedIndices.size(); ++i) { //loop over the 10 selected indices

            TH1D* current = hists[selectedIndices[i]]; //load the histogram (sample) corresponding to the random index
            Double_t weight = current->GetEntries(); //So we treat each histogram same (commented out later)
            int histoIndex = bs * 10 + i; //mine]

            if (!hmerge) { //First time initializing
                hmerge = dynamic_cast<TH1D*>(current->Clone( //clone the current
                    Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs) 
                ));
                //hmerge->Scale(histoWaights[histoIndex]); //MED VEKTER
                hmerge->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs));
                // hmerge->Scale(weight);
                // std::cout << "bs: " << bs << " idx: " << idx << " weight: " << weight << std::endl;
                totalWeight = histoWaights[histoIndex]; // * weight;
            } else {       // after the first initialization, we add
                // current->Scale(weight);
                // hmerge->Add(current, weight);
                hmerge->Add(current); //MED VEKTER
                //hmerge->Add(current);
                totalWeight += histoWaights[histoIndex]; // * weight;
            }
        }

        // 归一化并保存
        if (hmerge && totalWeight > 0) {
            hmerge->Scale(1.0 / selectedIndices.size());
            //hmerge->Scale(1.0 / totalWeight);
            // Printf("size: %d", selectedIndices.size());
            // hmerge->SetDirectory(outFile);
            hmerge->Write();
            // delete hmerge; // 释放内存
        }
    }

    // 清理资源
    outFile->Close();
    file->Close();
    delete outFile;
    delete file;
}




void readFakeFull(std::vector<int> &histoIndices, std::vector<double> &histoWaights, int maxSample, const std::string &fileName) {

    std::ifstream file(fileName);
    if (!file.is_open()) { // error
        std::cerr << "Error: could not open CSV file '" << fileName << "'.\n";
        return;
    }

    // Prepare vectors: clear previous contents and reserve expected size
    histoIndices.clear();
    histoWaights.clear();
    const size_t expected = static_cast<size_t>(maxSample) * maxSample * maxSample;
    histoIndices.reserve(expected);
    histoWaights.reserve(expected);

    std::string line;
    // Read header line and ignore it (if present)
    std::getline(file, line);
    // Read each subsequent line: integer,float
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string intStr, floatStr;

        if (!std::getline(ss, intStr, ',')) continue;
        if (!std::getline(ss, floatStr, ',')) continue;

        try {
            histoIndices.push_back(std::stoi(intStr));
        } catch (...) {
            // skip malformed integer
            continue;
        }
        try {
            histoWaights.push_back(std::stof(floatStr));
        } catch (...) {
            // keep vectors consistent: remove last int if float is bad
            histoIndices.pop_back();
            continue;
        }
    }
    file.close();
    if (histoIndices.size() != histoWaights.size()) {
        std::cerr << "Error: read different number of ints and floats (" << histoIndices.size()
                  << " vs " << histoWaights.size() << ").\n";
                  return;
    }

}