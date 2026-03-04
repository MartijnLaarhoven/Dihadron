/*
 * @Author: Martijn Laarhoven (martijn.laarhoven@cern.ch)
 * @Date: 2026-02-06
 * @Last Modified by: Martijn Laarhoven
 * @Last Modified time: 2026-02-06
 */
//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "TSystem.h"
#include "TH1D.h"
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
    Int_t periMin;
    Int_t periMax;

    InputUnit(std::string _fileNameSuffix, Bool_t _isNch, Bool_t _isEtadiff, Int_t _minRange, Int_t _maxRange, Int_t _periMin, Int_t _periMax) :
        fileNameSuffix(_fileNameSuffix), isNch(_isNch), isEtadiff(_isEtadiff), minRange(_minRange), maxRange(_maxRange), periMin(_periMin), periMax(_periMax) {}
};

void CreateBootstrapSample_PeripheralSub(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Int_t periMin, Int_t periMax);
void CreateBootstrapSample_PeripheralSub_EtaDiff(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax, Int_t periMin, Int_t periMax);

void Process_CreateBootstrapSample_PeripheralSubtracted() {
    // Available eta bins in current ProcessOutput/EtaDiff (FT0C)
    std::vector<float> etaBinsNeg = {-0.8,-0.7,-0.6};
    std::vector<float> etaBinsPos = {0.6,0.7,0.8};

    std::vector<InputUnit> inputList;
    const Int_t periMin = 80;
    const Int_t periMax = 100;

    std::vector<std::string> datasets = {
        "LHC25ae_pass2_616549",
        "LHC25ae_pass2_618685",
        "LHC25af_pass2_615818",
        "LHC25af_pass2_615817"
    };

    for (const auto& ds : datasets) {
        inputList.push_back(InputUnit(ds, kCent, kEtaDiffOn, 0, 20, periMin, periMax));
    }

    for (auto input : inputList) {
        std::cout << "Processing Bootstrap Sample (PeripheralSub) eta diff: " << input.fileNameSuffix << std::endl;
        for (int iEta = 0; iEta < (int)etaBinsNeg.size() - 1; iEta++) {
            double etaMin = etaBinsNeg[iEta];
            double etaMax = etaBinsNeg[iEta + 1];
            CreateBootstrapSample_PeripheralSub_EtaDiff(input.fileNameSuffix, input.isNch, input.minRange, input.maxRange, etaMin, etaMax, input.periMin, input.periMax);
        }
        for (int iEta = 0; iEta < (int)etaBinsPos.size() - 1; iEta++) {
            double etaMin = etaBinsPos[iEta];
            double etaMax = etaBinsPos[iEta + 1];
            CreateBootstrapSample_PeripheralSub_EtaDiff(input.fileNameSuffix, input.isNch, input.minRange, input.maxRange, etaMin, etaMax, input.periMin, input.periMax);
        }
    }
}

void CreateBootstrapSample_PeripheralSub(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Int_t periMin, Int_t periMax) {
    std::string splitName = isNch ? "Mult" : "Cent";
    
    // Read peripheral-subtracted file
    TFile* file = TFile::Open(
        Form("./ProcessOutput/PeripheralSubtraction/PeripheralSub_%s_%s_%i_%i_Peri_%i_%i.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, periMin, periMax), 
        "READ"
    );
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening peripheral subtraction file!" << std::endl;
        return;
    }

    // Get the subtracted histogram
    TH1D* hSub = dynamic_cast<TH1D*>(file->Get("hPhiSameOverMixed_Sub"));
    if (!hSub) {
        std::cerr << "Error loading subtracted histogram!" << std::endl;
        file->Close();
        delete file;
        return;
    }

    // Create output file for bootstrap samples
    gSystem->mkdir("./ProcessOutput/PeripheralSubtraction", kTRUE);
    TFile* outFile = TFile::Open(
        Form("./ProcessOutput/PeripheralSubtraction/BootstrapSample_%s_%s_%i_%i_PeripheralSub.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange), 
        "RECREATE"
    );
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error creating output file!" << std::endl;
        file->Close();
        delete file;
        return;
    }

    TRandom3 rnd(0);
    Int_t nbins = hSub->GetNbinsX();

    // Generate bootstrap samples
    for (Int_t sample = 0; sample < maxSample * maxSample; ++sample) {
        TH1D* hBootstrap = new TH1D(
            Form("hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, sample),
            Form("Bootstrap sample %d", sample),
            nbins, hSub->GetXaxis()->GetXmin(), hSub->GetXaxis()->GetXmax()
        );
        hBootstrap->Sumw2();

        for (Int_t bin = 1; bin <= nbins; ++bin) {
            Double_t content = hSub->GetBinContent(bin);
            Double_t error = hSub->GetBinError(bin);
            Double_t value = rnd.Gaus(content, error);
            hBootstrap->SetBinContent(bin, value);
            hBootstrap->SetBinError(bin, error);
        }

        hBootstrap->Write();
        delete hBootstrap;
    }

    outFile->Close();
    file->Close();
    delete outFile;
    delete file;
    
    std::cout << "Bootstrap samples (PeripheralSub) created: " 
              << Form("./ProcessOutput/PeripheralSubtraction/BootstrapSample_%s_%s_%i_%i_PeripheralSub.root", 
                     fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange) << std::endl;
}

void CreateBootstrapSample_PeripheralSub_EtaDiff(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax, Int_t periMin, Int_t periMax) {
    std::string splitName = isNch ? "Mult" : "Cent";
    
    // Read peripheral-subtracted file
    TFile* file = TFile::Open(
        Form("./ProcessOutput/PeripheralSubtraction/EtaDiff/PeripheralSub_%s_%s_%i_%i_Eta_%0.1f_%0.1f_Peri_%i_%i.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax, periMin, periMax), 
        "READ"
    );
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening peripheral subtraction file (eta diff)!" << std::endl;
        return;
    }

    // Get the subtracted histogram
    TH1D* hSub = dynamic_cast<TH1D*>(file->Get("hPhiSameOverMixed_Sub"));
    if (!hSub) {
        std::cerr << "Error loading subtracted histogram (eta diff)!" << std::endl;
        file->Close();
        delete file;
        return;
    }

    // Create output file for bootstrap samples
    gSystem->mkdir("./ProcessOutput/PeripheralSubtraction/EtaDiff", kTRUE);
    TFile* outFile = TFile::Open(
        Form("./ProcessOutput/PeripheralSubtraction/EtaDiff/BootstrapSample_%s_%s_%i_%i_Eta_%0.1f_%0.1f_PeripheralSub.root", 
             fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax), 
        "RECREATE"
    );
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error creating output file (eta diff)!" << std::endl;
        file->Close();
        delete file;
        return;
    }

    TRandom3 rnd(0);
    Int_t nbins = hSub->GetNbinsX();

    // Generate bootstrap samples
    for (Int_t sample = 0; sample < maxSample * maxSample; ++sample) {
        TH1D* hBootstrap = new TH1D(
            Form("hPhiSameOverMixed_%d_%d_%d", minRange, maxRange, sample),
            Form("Bootstrap sample %d", sample),
            nbins, hSub->GetXaxis()->GetXmin(), hSub->GetXaxis()->GetXmax()
        );
        hBootstrap->Sumw2();

        for (Int_t bin = 1; bin <= nbins; ++bin) {
            Double_t content = hSub->GetBinContent(bin);
            Double_t error = hSub->GetBinError(bin);
            Double_t value = rnd.Gaus(content, error);
            hBootstrap->SetBinContent(bin, value);
            hBootstrap->SetBinError(bin, error);
        }

        hBootstrap->Write();
        delete hBootstrap;
    }

    outFile->Close();
    file->Close();
    delete outFile;
    delete file;
    
    std::cout << "Bootstrap samples (PeripheralSub, eta diff) created: " 
              << Form("./ProcessOutput/PeripheralSubtraction/EtaDiff/BootstrapSample_%s_%s_%i_%i_Eta_%0.1f_%0.1f_PeripheralSub.root", 
                     fileNameSuffix.c_str(), splitName.c_str(), minRange, maxRange, etaMin, etaMax) << std::endl;
}
