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
#include "TH2D.h"
#include "TParameter.h"
#include "TROOT.h"
#include "TCanvas.h"
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
	Bool_t isMc;

	InputUnit(std::string _fileNameSuffix, Bool_t _isNch, Bool_t _isEtadiff, Int_t _minRange, Int_t _maxRange, Bool_t _isMc=false) :
		fileNameSuffix(_fileNameSuffix), isNch(_isNch), isEtadiff(_isEtadiff), minRange(_minRange), maxRange(_maxRange), isMc(_isMc) {}
};

// global variables
std::string collisionSystemName = "peripheral NeNe";
std::string additionalSuffix = "";

// Helper function to set collision system name based on filename
std::string GetCollisionSystemName(const std::string& filename) {
	if (filename.find("ae") != std::string::npos) {
		return "O-O";
	} else if (filename.find("af") != std::string::npos) {
		return "Ne-Ne";
	}
	return "Unknown";
}

double ComputeScaleByIntegral(TH1D* hCent, TH1D* hPeri) {
	if (!hCent || !hPeri) return 1.0;
	double intCent = hCent->Integral();
	double intPeri = hPeri->Integral();
	if (intPeri <= 0.0) return 1.0;
	return intCent / intPeri;
}

// ZYAM: Find the minimum (baseline) of the correlation function around ∆φ ∼ 0
double ComputeZYAMBaseline(TH1D* hPeri) {
	if (!hPeri) return 0.0;
	
	// Find minimum in range around ∆φ ∼ 0
	// Typically within ±π/6 (±30 degrees) from ∆φ = 0
	double minValue = hPeri->GetBinContent(1);
	Int_t nbins = hPeri->GetNbinsX();
	
	// Search in the first and last bins (around ∆φ = 0, ±2π)
	Int_t searchRange = TMath::Max(1, (Int_t)(nbins / 12)); // ±π/6
	
	for (Int_t i = 1; i <= searchRange; ++i) {
		minValue = TMath::Min(minValue, hPeri->GetBinContent(i));
		minValue = TMath::Min(minValue, hPeri->GetBinContent(nbins - i + 1));
	}
	
	return minValue;
}

void PeripheralSubtraction_givenRange(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Int_t periMin, Int_t periMax, Bool_t isMc);
void PeripheralSubtraction_givenRange_EtaDiff(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax, Int_t periMin, Int_t periMax, Bool_t isMc);

void Process_PeripheralSubtraction() {
	// 不显示窗口
	gROOT->SetBatch(kTRUE);

	// Define eta bins for both configurations
	std::vector<float> etaBinsNeg = {-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0};
	std::vector<float> etaBinsPos = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

	std::vector<InputUnit> inputList;
	additionalSuffix = "";

	// Set collision system name based on filename
	std::string inputFileName = "LHC25ae_pass2_604826";
	collisionSystemName = GetCollisionSystemName(inputFileName);

	// Centrality bins to subtract peripheral from
	inputList.push_back(InputUnit(inputFileName, kCent, kEtaDiffOff, 0, 20));
	inputList.push_back(InputUnit(inputFileName, kCent, kEtaDiffOn, 0, 20));

	// Also run for the new dataset (same configuration)
	std::string inputFileNameNew = "LHC25ae_pass2_604830";
	inputList.push_back(InputUnit(inputFileNameNew, kCent, kEtaDiffOff, 0, 20));
	inputList.push_back(InputUnit(inputFileNameNew, kCent, kEtaDiffOn, 0, 20));

	// Also run for Ne-Ne datasets
	std::string inputFileNameNeStd = "LHC25af_pass2_611697";
	inputList.push_back(InputUnit(inputFileNameNeStd, kCent, kEtaDiffOff, 0, 20));
	inputList.push_back(InputUnit(inputFileNameNeStd, kCent, kEtaDiffOn, 0, 20));

	std::string inputFileNameNeRev = "LHC25af_pass2_604820";
	inputList.push_back(InputUnit(inputFileNameNeRev, kCent, kEtaDiffOff, 0, 20));
	inputList.push_back(InputUnit(inputFileNameNeRev, kCent, kEtaDiffOn, 0, 20));

	const Int_t periMin = 80;
	const Int_t periMax = 100;

	for (auto input : inputList) {
		if (input.isEtadiff) {
			// Use positive bins for reversed datasets, negative bins for standard
			std::vector<float>& etaBinsToUse = (input.fileNameSuffix == inputFileNameNew || input.fileNameSuffix == inputFileNameNeRev) ? etaBinsPos : etaBinsNeg;
			for (int iEta = 0; iEta < etaBinsToUse.size() - 1; iEta++) {
				double etaMin = etaBinsToUse[iEta];
				double etaMax = etaBinsToUse[iEta + 1];
				PeripheralSubtraction_givenRange_EtaDiff(input.fileNameSuffix, input.isNch, input.minRange, input.maxRange, etaMin, etaMax, periMin, periMax, input.isMc);
			}
		}
		else {
			PeripheralSubtraction_givenRange(input.fileNameSuffix, input.isNch, input.minRange, input.maxRange, periMin, periMax, input.isMc);
		}
	}
}

void PeripheralSubtraction_givenRange(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Int_t periMin, Int_t periMax, Bool_t isMc) {
	std::string splitName = isNch ? "Mult" : "Cent";
	std::string inputCentral = Form("./ProcessOutput/Mixed_%s%s_%s_%i_%i.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange));
	std::string inputPeripheral = Form("./ProcessOutput/Mixed_%s%s_%s_%i_%i.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(periMin), int(periMax));

	TFile* fCent = TFile::Open(inputCentral.c_str(), "READ");
	if (!fCent || !fCent->IsOpen()) {
		std::cerr << "Cannot open central file: " << inputCentral << std::endl;
		return;
	}
	TFile* fPeri = TFile::Open(inputPeripheral.c_str(), "READ");
	if (!fPeri || !fPeri->IsOpen()) {
		std::cerr << "Cannot open peripheral file: " << inputPeripheral << std::endl;
		fCent->Close();
		delete fCent;
		return;
	}

	TH1D* hCent = (TH1D*)fCent->Get(Form("hPhiSameOverMixed_%d_%d", minRange, maxRange));
	TH1D* hPeri = (TH1D*)fPeri->Get(Form("hPhiSameOverMixed_%d_%d", periMin, periMax));
	if (!hCent || !hPeri) {
		std::cerr << "Missing hPhiSameOverMixed histograms for peripheral subtraction." << std::endl;
		fPeri->Close();
		fCent->Close();
		delete fPeri;
		delete fCent;
		return;
	}

	double scale = ComputeScaleByIntegral(hCent, hPeri);
	double y0 = ComputeZYAMBaseline(hPeri);
	
	// ZYAM method: subtract (F_peri - Y0) from central
	// This gives: F_central - F_peri + Y0
	TH1D* hSub = (TH1D*)hCent->Clone(Form("hPhiSameOverMixed_Sub_%d_%d", minRange, maxRange));
	hSub->SetTitle(Form("Central - (Peripheral - Y0) (ZYAM, Y0=%.4f); #Delta#phi [rad]; S/M", y0));
	
	// Subtract peripheral with ZYAM baseline subtraction
	for (Int_t bin = 1; bin <= hPeri->GetNbinsX(); ++bin) {
		Double_t centContent = hCent->GetBinContent(bin);
		Double_t periContent = hPeri->GetBinContent(bin);
		Double_t periError = hPeri->GetBinError(bin);
		Double_t centError = hCent->GetBinError(bin);
		
		// Apply ZYAM: subtract (peri - y0)
		Double_t subContent = centContent - (periContent - y0);
		Double_t subError = TMath::Sqrt(centError*centError + periError*periError);
		
		hSub->SetBinContent(bin, subContent);
		hSub->SetBinError(bin, subError);
	}

	gSystem->mkdir("./ProcessOutput/PeripheralSubtraction", kTRUE);
	std::string outFile = Form("./ProcessOutput/PeripheralSubtraction/PeripheralSub_%s%s_%s_%i_%i_Peri_%i_%i.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), int(periMin), int(periMax));
	TFile* fout = TFile::Open(outFile.c_str(), "RECREATE");
	if (!fout || !fout->IsOpen()) {
		std::cerr << "Cannot create output file: " << outFile << std::endl;
	} else {
		TParameter<double> scaleParam("PeripheralScale", scale);
		TParameter<double> zyamBaseline("ZYAMBaseline", y0);
		hCent->Write("hPhiSameOverMixed_Central");
		hPeri->Write("hPhiSameOverMixed_Peripheral");
		hSub->Write("hPhiSameOverMixed_Sub");
		scaleParam.Write();
		zyamBaseline.Write();
		fout->Close();
		delete fout;
		std::cout << "Output file: " << outFile << " (ZYAM Y0=" << y0 << ")" << std::endl;
	}

	delete hSub;
	fPeri->Close();
	fCent->Close();
	delete fPeri;
	delete fCent;
}

void PeripheralSubtraction_givenRange_EtaDiff(std::string fileNameSuffix, Bool_t isNch, Int_t minRange, Int_t maxRange, Double_t etaMin, Double_t etaMax, Int_t periMin, Int_t periMax, Bool_t isMc) {
	std::string splitName = isNch ? "Mult" : "Cent";
	std::string inputCentral = Form("./ProcessOutput/EtaDiff/Mixed_%s%s_%s_%i_%i_Eta_%0.1f_%0.1f.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), etaMin, etaMax);
	std::string inputPeripheral = Form("./ProcessOutput/EtaDiff/Mixed_%s%s_%s_%i_%i_Eta_%0.1f_%0.1f.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(periMin), int(periMax), etaMin, etaMax);

	TFile* fCent = TFile::Open(inputCentral.c_str(), "READ");
	if (!fCent || !fCent->IsOpen()) {
		std::cerr << "Cannot open central file: " << inputCentral << std::endl;
		return;
	}
	TFile* fPeri = TFile::Open(inputPeripheral.c_str(), "READ");
	if (!fPeri || !fPeri->IsOpen()) {
		std::cerr << "Cannot open peripheral file: " << inputPeripheral << std::endl;
		fCent->Close();
		delete fCent;
		return;
	}

	TH1D* hCent = (TH1D*)fCent->Get(Form("hPhiSameOverMixed_%d_%d", minRange, maxRange));
	TH1D* hPeri = (TH1D*)fPeri->Get(Form("hPhiSameOverMixed_%d_%d", periMin, periMax));
	TH2D* hCent2D = (TH2D*)fCent->Get("hPhiEtaSMsum");
	TH2D* hPeri2D = (TH2D*)fPeri->Get("hPhiEtaSMsum");
	if (!hCent || !hPeri) {
		std::cerr << "Missing hPhiSameOverMixed histograms for peripheral subtraction (eta diff)." << std::endl;
		fPeri->Close();
		fCent->Close();
		delete fPeri;
		delete fCent;
		return;
	}

	double scale = ComputeScaleByIntegral(hCent, hPeri);
	double y0 = ComputeZYAMBaseline(hPeri);
	
	// ZYAM method: subtract (F_peri - Y0) from central
	TH1D* hSub = (TH1D*)hCent->Clone(Form("hPhiSameOverMixed_Sub_%d_%d", minRange, maxRange));
	hSub->SetTitle(Form("Central - (Peripheral - Y0) (ZYAM, Y0=%.4f); #Delta#phi [rad]; S/M", y0));
	
	// Subtract peripheral with ZYAM baseline subtraction
	for (Int_t bin = 1; bin <= hPeri->GetNbinsX(); ++bin) {
		Double_t centContent = hCent->GetBinContent(bin);
		Double_t periContent = hPeri->GetBinContent(bin);
		Double_t periError = hPeri->GetBinError(bin);
		Double_t centError = hCent->GetBinError(bin);
		
		// Apply ZYAM: subtract (peri - y0)
		Double_t subContent = centContent - (periContent - y0);
		Double_t subError = TMath::Sqrt(centError*centError + periError*periError);
		
		hSub->SetBinContent(bin, subContent);
		hSub->SetBinError(bin, subError);
	}

	TH2D* hSub2D = nullptr;
	if (hCent2D && hPeri2D) {
		hSub2D = (TH2D*)hCent2D->Clone("hPhiEtaSMsum_Sub");
		hSub2D->SetTitle(Form("Central - (Peripheral - Y0) (ZYAM, Y0=%.4f); #Delta#phi [rad]; #Delta#eta", y0));
		
		// Apply ZYAM to 2D histogram
		for (Int_t binx = 1; binx <= hPeri2D->GetNbinsX(); ++binx) {
			for (Int_t biny = 1; biny <= hPeri2D->GetNbinsY(); ++biny) {
				Double_t centContent = hCent2D->GetBinContent(binx, biny);
				Double_t periContent = hPeri2D->GetBinContent(binx, biny);
				Double_t periError = hPeri2D->GetBinError(binx, biny);
				Double_t centError = hCent2D->GetBinError(binx, biny);
				
				Double_t subContent = centContent - (periContent - y0);
				Double_t subError = TMath::Sqrt(centError*centError + periError*periError);
				
				hSub2D->SetBinContent(binx, biny, subContent);
				hSub2D->SetBinError(binx, biny, subError);
			}
		}
	}

	if (hSub2D) {
		TCanvas* c2D = new TCanvas(Form("c2D_Sub_%d_%d_%0.1f_%0.1f", int(minRange), int(maxRange), etaMin, etaMax),
			"Peripheral Subtraction 2D", 1200, 400);
		c2D->Divide(3, 1);
		c2D->cd(1);
		TH2D* hCent2D_draw = (TH2D*)hCent2D->Clone("hCent2D_draw");
		hCent2D_draw->GetYaxis()->SetRangeUser(-1.5, 1.5);
		hCent2D_draw->Draw("COLZ");
		c2D->cd(2);
		TH2D* hPeri2D_draw = (TH2D*)hPeri2D->Clone("hPeri2D_draw");
		hPeri2D_draw->GetYaxis()->SetRangeUser(-1.5, 1.5);
		hPeri2D_draw->Draw("COLZ");
		c2D->cd(3);
		TH2D* hSub2D_draw = (TH2D*)hSub2D->Clone("hSub2D_draw");
		hSub2D_draw->GetYaxis()->SetRangeUser(-1.5, 1.5);
		hSub2D_draw->Draw("COLZ");

		// Write canvas and cleanup before opening output file
		gSystem->mkdir("./ProcessOutput/PeripheralSubtraction/EtaDiff", kTRUE);
		std::string outFile = Form("./ProcessOutput/PeripheralSubtraction/EtaDiff/PeripheralSub_%s%s_%s_%i_%i_Eta_%0.1f_%0.1f_Peri_%i_%i.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), etaMin, etaMax, int(periMin), int(periMax));
		TFile* fout = TFile::Open(outFile.c_str(), "RECREATE");
		if (!fout || !fout->IsOpen()) {
			std::cerr << "Cannot create output file: " << outFile << std::endl;
			delete c2D;
			delete hCent2D_draw;
			delete hPeri2D_draw;
			delete hSub2D_draw;
		} else {
			c2D->Write();
			TParameter<double> scaleParam("PeripheralScale", scale);
			TParameter<double> zyamBaseline("ZYAMBaseline", y0);
			hCent->Write("hPhiSameOverMixed_Central");
			hPeri->Write("hPhiSameOverMixed_Peripheral");
			hSub->Write("hPhiSameOverMixed_Sub");
			hCent2D->Write("hPhiEtaSMsum_Central");
			hPeri2D->Write("hPhiEtaSMsum_Peripheral");
			hSub2D->Write("hPhiEtaSMsum_Sub");
			scaleParam.Write();
			zyamBaseline.Write();
			fout->Close();
			delete fout;
			std::cout << "Output file: " << outFile << " (ZYAM Y0=" << y0 << ")" << std::endl;
		}
		delete c2D;
		delete hCent2D_draw;
		delete hPeri2D_draw;
		delete hSub2D_draw;
	} else {
		// No 2D histograms, just write 1D results
		gSystem->mkdir("./ProcessOutput/PeripheralSubtraction/EtaDiff", kTRUE);
		std::string outFile = Form("./ProcessOutput/PeripheralSubtraction/EtaDiff/PeripheralSub_%s%s_%s_%i_%i_Eta_%0.1f_%0.1f_Peri_%i_%i.root", fileNameSuffix.c_str(), additionalSuffix.c_str(), splitName.c_str(), int(minRange), int(maxRange), etaMin, etaMax, int(periMin), int(periMax));
		TFile* fout = TFile::Open(outFile.c_str(), "RECREATE");
		if (!fout || !fout->IsOpen()) {
			std::cerr << "Cannot create output file: " << outFile << std::endl;
		} else {
			TParameter<double> scaleParam("PeripheralScale", scale);
			TParameter<double> zyamBaseline("ZYAMBaseline", y0);
			hCent->Write("hPhiSameOverMixed_Central");
			hPeri->Write("hPhiSameOverMixed_Peripheral");
			hSub->Write("hPhiSameOverMixed_Sub");
			scaleParam.Write();
			zyamBaseline.Write();
			fout->Close();
			delete fout;
			std::cout << "Output file: " << outFile << " (ZYAM Y0=" << y0 << ")" << std::endl;
		}
	}

	delete hSub2D;
	delete hSub;
	fPeri->Close();
	fCent->Close();
	delete fPeri;
	delete fCent;
}
