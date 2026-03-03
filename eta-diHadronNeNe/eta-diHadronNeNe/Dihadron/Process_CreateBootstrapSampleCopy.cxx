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
#include "TRandom3.h"
#include <iostream>
#include <string>
#include <vector>
#include "./include/BasicForDihadron.h"

struct InputUnit {
    std::string fileNameSuffix;
    Bool_t isNch;
    Bool_t isEtadiff;
    Int_t minRange;
    Int_t maxRange;

    InputUnit(std::string _fileNameSuffix, Bool_t _isNch, Bool_t _isEtadiff, Int_t _minRange, Int_t _maxRange) :
        fileNameSuffix(_fileNameSuffix), isNch(_isNch), isEtadiff(_isEtadiff), minRange(_minRange), maxRange(_maxRange) {}
};

void CreateBootstrapSample(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange);
void CreateBootstrapSample_PtDiff(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t pTMin, Double_t pTMax);
void CreateBootstrapSample_EtaDiff(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax);
void Process_CreateBootstrapSample() {

    std::vector<InputUnit> inputList;
    
    inputList.push_back(InputUnit("LHC25ae_pass2_604826", kCent, kEtaDiffOff, 0, 20));
    inputList.push_back(InputUnit("LHC25ae_pass2_604826", kCent, kEtaDiffOn, 0, 20));

    inputList.push_back(InputUnit("LHC25ae_pass2_604826", kCent, kEtaDiffOff, 80, 100));    
    inputList.push_back(InputUnit("LHC25ae_pass2_604826", kCent, kEtaDiffOn, 80, 100));
    for (auto input : inputList) {
        if (input.isEtadiff) {
            std::cout << "Processing Bootstrap Sample eta diff: " << input.fileNameSuffix << std::endl;
            for (int iEta = 0; iEta < etaBins.size() - 1; iEta++) {
                double etaMin = etaBins[iEta];
                double etaMax = etaBins[iEta + 1];
                CreateBootstrapSample_EtaDiff(input.fileNameSuffix, input.isNch, input.minRange, input.maxRange, etaMin, etaMax);
            }
        } else {
            std::cout << "Processing Bootstrap Sample: " << input.fileNameSuffix << std::endl;
            CreateBootstrapSample(input.fileNameSuffix, input.isNch, input.minRange, input.maxRange);
        }
    }

}

void CreateBootstrapSample(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    
    // 打开输入文件
    TFile* file = TFile::Open(
        Form("./ProcessOutput/Mixed_%s_%s_%i_%i.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange), 
        "READ"
    );
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
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

    // 创建输出文件
    TFile* outFile = TFile::Open(
        Form("./ProcessOutput/BootstrapSample_%s_%s_%i_%i.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange), 
        "RECREATE"
    );
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
    TRandom3 randGen;

    // 生成 maxSample^2 个bootstrap样本
    for (Int_t bs = 0; bs < maxSample * maxSample; ++bs) {
        // 随机选择样本索引（允许重复）
        std::vector<Int_t> selectedIndices;
        selectedIndices.clear();
        for (Int_t i = 0; i < maxSample; ++i) {
            selectedIndices.push_back(randGen.Integer(maxSample));
        }

        // 合并选中的直方图
        TH1D* hmerge = nullptr;
        Double_t totalWeight = 0.0;

        for (Int_t idx : selectedIndices) {
            if (idx < 0 || idx >= hists.size()) {
                std::cerr << "Invalid sample index: " << idx << std::endl;
                continue;
            }

            TH1D* current = hists[idx];
            Double_t weight = current->GetEntries();

            if (!hmerge) { // 第一次初始化合并直方图
                hmerge = dynamic_cast<TH1D*>(current->Clone(
                    Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs)
                ));
                hmerge->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs));
                // hmerge->Scale(weight);
                // std::cout << "bs: " << bs << " idx: " << idx << " weight: " << weight << std::endl;
                totalWeight = weight;
            } else {       // 后续累加
                // current->Scale(weight);
                // hmerge->Add(current, weight);
                hmerge->Add(current);
                totalWeight += weight;
            }
        }

        // 归一化并保存
        if (hmerge && totalWeight > 0) {
            // hmerge->Scale(1.0 / totalWeight);
            hmerge->Scale(1.0 / selectedIndices.size());
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

void CreateBootstrapSample_PtDiff(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t pTMin, Double_t pTMax) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    
    // 打开输入文件
    TFile* file = TFile::Open(
        Form("./ProcessOutput/PtDiff/Mixed_%s_%s_%i_%i_Pt_%0.1f_%0.1f.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, pTMin, pTMax), 
        "READ"
    );
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
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

    // 创建输出文件
    TFile* outFile = TFile::Open(
        Form("./ProcessOutput/PtDiff/BootstrapSample_%s_%s_%i_%i_Pt_%0.1f_%0.1f.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, pTMin, pTMax), 
        "RECREATE"
    );
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
    TRandom3 randGen;

    // 生成 maxSample^2 个bootstrap样本
    for (Int_t bs = 0; bs < maxSample * maxSample; ++bs) {
        // 随机选择样本索引（允许重复）
        std::vector<Int_t> selectedIndices;
        selectedIndices.clear();
        for (Int_t i = 0; i < maxSample; ++i) {
            selectedIndices.push_back(randGen.Integer(maxSample));
        }

        // 合并选中的直方图
        TH1D* hmerge = nullptr;
        Double_t totalWeight = 0.0;

        for (Int_t idx : selectedIndices) {
            if (idx < 0 || idx >= hists.size()) {
                std::cerr << "Invalid sample index: " << idx << std::endl;
                continue;
            }

            TH1D* current = hists[idx];
            Double_t weight = current->GetEntries();

            if (!hmerge) { // 第一次初始化合并直方图
                hmerge = dynamic_cast<TH1D*>(current->Clone(
                    Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs)
                ));
                hmerge->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs));
                // hmerge->Scale(weight);
                // std::cout << "bs: " << bs << " idx: " << idx << " weight: " << weight << std::endl;
                totalWeight = weight;
            } else {       // 后续累加
                // current->Scale(weight);
                // hmerge->Add(current, weight);
                hmerge->Add(current);
                totalWeight += weight;
            }
        }

        // 归一化并保存
        if (hmerge && totalWeight > 0) {
            // hmerge->Scale(1.0 / totalWeight);
            hmerge->Scale(1.0 / selectedIndices.size());
            // Printf("size: %d", selectedIndices.size());
            // hmerge->SetDirectory(outFile);
            hmerge->Write();
            delete hmerge; // 释放内存
        }
    }

    // 清理资源
    outFile->Close();
    file->Close();
    delete outFile;
    delete file;
}

void CreateBootstrapSample_EtaDiff(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax) {
    std::string splitName = "Mult";
    if (!isNch) splitName = "Cent";
    
    // 打开输入文件
    TFile* file = TFile::Open(
        Form("./ProcessOutput/EtaDiff/Mixed_%s_%s_%i_%i_Eta_%0.1f_%0.1f.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax), 
        "READ"
    );
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
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

    // 创建输出文件
    TFile* outFile = TFile::Open(
        Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_%s_%i_%i_Eta_%0.1f_%0.1f.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax), 
        "RECREATE"
    );
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
    TRandom3 randGen;

    // 生成 maxSample^2 个bootstrap样本
    for (Int_t bs = 0; bs < maxSample * maxSample; ++bs) {
        // 随机选择样本索引（允许重复）
        std::vector<Int_t> selectedIndices;
        selectedIndices.clear();
        for (Int_t i = 0; i < maxSample; ++i) {
            selectedIndices.push_back(randGen.Integer(maxSample));
        }

        // 合并选中的直方图
        TH1D* hmerge = nullptr;
        Double_t totalWeight = 0.0;

        for (Int_t idx : selectedIndices) {
            if (idx < 0 || idx >= hists.size()) {
                std::cerr << "Invalid sample index: " << idx << std::endl;
                continue;
            }

            TH1D* current = hists[idx];
            Double_t weight = current->GetEntries();

            if (!hmerge) {
                hmerge = dynamic_cast<TH1D*>(current->Clone(
                    Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs)
                ));
                hmerge->SetTitle(Form("bsSample_hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, bs));
                totalWeight = weight;
            } else {
                hmerge->Add(current);
                totalWeight += weight;
            }
        }

        // 归一化并保存
        if (hmerge && totalWeight > 0) {
            hmerge->Scale(1.0 / selectedIndices.size());
            hmerge->Write();
        }
    }

    // 清理资源
    outFile->Close();
    file->Close();
    delete outFile;
    delete file;
}
