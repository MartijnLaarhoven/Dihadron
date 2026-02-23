/*
 * Extract v2 from eta-differential bootstrap samples using direct method
 * Instead of template fitting (which fails with sparse data), extract v2 as:
 * v2 = 2 * (y[peak] - y[valley]) where peak is near pi/2 and valley is near pi
 */
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "./include/BasicForDihadron.h"

struct VnUnit {
    Double_t v2; Double_t v2Err;
    Double_t v3; Double_t v3Err;
    VnUnit(Double_t _v2 = 0, Double_t _v2Err = 0, Double_t _v3 = 0, Double_t _v3Err = 0) :
        v2(_v2), v2Err(_v2Err), v3(_v3), v3Err(_v3Err) {}
};

Double_t ExtractV2FromHistogram(TH1D* hist) {
    if (!hist || hist->Integral() <= 0) return 0;
    
    // Find peak near pi/2 and valley near pi
    Double_t peakValue = 0, valleyValue = 0;
    Int_t nBins = hist->GetNbinsX();
    
    // Find peak (bin closest to pi/2)
    Int_t peakBin = hist->FindBin(TMath::Pi() / 2.0);
    peakValue = hist->GetBinContent(peakBin);
    
    // Find valley (bin closest to pi)
    Int_t valleyBin = hist->FindBin(TMath::Pi());
    if (valleyBin > nBins) valleyBin = nBins;
    valleyValue = hist->GetBinContent(valleyBin);
    
    // v2 = 2 * (peak - valley) / peak (normalized)
    if (peakValue <= 0) return 0;
    Double_t v2 = 2.0 * (peakValue - valleyValue) / peakValue;
    return v2;
}

void ProcessEtaDiffExtraction() {
    std::vector<float> etaBinsNeg = {-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0};
    std::vector<float> etaBinsPos = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    
    std::vector<std::string> datasets = {
        "LHC25ae_pass2_604826",
        "LHC25ae_pass2_604830",
        "LHC25af_pass2_604820"
    };
    
    std::vector<std::vector<float>*> etaBinsList = {
        &etaBinsNeg,  // 604826 standard
        &etaBinsPos,  // 604830 reversed
        &etaBinsPos   // 604820 reversed
    };
    
    // Process each dataset and eta bin
    for (size_t dIdx = 0; dIdx < datasets.size(); ++dIdx) {
        std::string dataset = datasets[dIdx];
        std::vector<float>& etaBins = *etaBinsList[dIdx];
        
        std::cout << "\n========== Processing " << dataset << " ==========" << std::endl;
        
        // Process 0-20% and 80-100%
        for (Int_t cent : {0, 80}) {
            Int_t centMax = (cent == 0) ? 20 : 100;
            std::cout << "Centrality: " << cent << "-" << centMax << "%" << std::endl;
            
            // Create output file for v2 results
            TFile* outFile = TFile::Open(
                Form("./TemplateFit/EtaDiff/Vn_%s_Cent_%d_%d.root", 
                     dataset.c_str(), cent, centMax), 
                "RECREATE"
            );
            
            TGraphErrors* gV2 = new TGraphErrors();
            gV2->SetName("gV2Delta");
            gV2->SetTitle(Form("v2 vs eta - %s Cent %d-%d%%", dataset.c_str(), cent, centMax));
            gV2->SetMarkerStyle(20);
            gV2->SetMarkerSize(0.5);
            
            int pointIdx = 0;
            
            // Process each eta bin
            for (size_t iEta = 0; iEta < etaBins.size() - 1; ++iEta) {
                Double_t etaMin = etaBins[iEta];
                Double_t etaMax = etaBins[iEta + 1];
                Double_t etaMid = (etaMin + etaMax) / 2.0;
                
                // Read bootstrap samples
                TFile* bootFile = TFile::Open(
                    Form("./ProcessOutput/EtaDiff/BootstrapSample_%s_Cent_%d_%d_Eta_%0.1f_%0.1f.root",
                         dataset.c_str(), cent, centMax, etaMin, etaMax),
                    "READ"
                );
                
                if (!bootFile || bootFile->IsZombie()) {
                    std::cerr << "Cannot open bootstrap file for eta [" << etaMin << ", " << etaMax << "]" << std::endl;
                    continue;
                }
                
                std::vector<Double_t> v2Values;
                
                // Extract v2 from each bootstrap sample
                for (Int_t sample = 0; sample < 900; ++sample) {
                    TH1D* hist = dynamic_cast<TH1D*>(
                        bootFile->Get(Form("bsSample_hPhiSameOverMixed_%d_%d_%d", cent, centMax, sample))
                    );
                    
                    if (!hist) continue;
                    
                    Double_t v2 = ExtractV2FromHistogram(hist);
                    if (v2 > 0) {
                        v2Values.push_back(v2);
                    }
                }
                
                bootFile->Close();
                delete bootFile;
                
                if (v2Values.empty()) {
                    std::cerr << "No valid v2 extracted for eta [" << etaMin << ", " << etaMax << "]" << std::endl;
                    continue;
                }
                
                // Calculate mean and error
                Double_t v2Mean = 0, v2Err = 0;
                for (Double_t v : v2Values) v2Mean += v;
                v2Mean /= v2Values.size();
                
                for (Double_t v : v2Values) {
                    Double_t diff = v - v2Mean;
                    v2Err += diff * diff;
                }
                v2Err = std::sqrt(v2Err / v2Values.size());
                
                std::cout << "  Eta [" << etaMin << ", " << etaMax << "]: v2 = " 
                          << v2Mean << " +/- " << v2Err << std::endl;
                
                gV2->SetPoint(pointIdx, etaMid, v2Mean);
                gV2->SetPointError(pointIdx, 0, v2Err);
                pointIdx++;
            }
            
            // Write graph
            outFile->cd();
            gV2->Write();
            outFile->Close();
            delete outFile;
            
            std::cout << "Output written to TemplateFit/EtaDiff/Vn_" << dataset << "_Cent_" 
                      << cent << "_" << centMax << ".root" << std::endl;
        }
    }
}

void Process_ExtractV2_EtaDiff() {
    ProcessEtaDiffExtraction();
}
