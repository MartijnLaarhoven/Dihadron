/*
 * @Author: Martijn Laarhoven (martijn.laarhoven@cern.ch)
 * @Date: 2026-03-03
 * Compare Uncorrected FourierFit vs TemplateFit methods for V2Delta, V3Delta, V4Delta
 */

#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

void Process_CompareFourierMethods() {
    // First run the Uncorrected FourierFit if not already done
    std::cout << "Running Uncorrected FourierFit..." << std::endl;
    gROOT->ProcessLine(".L Process_Uncorrected_FourierFit.cxx");
    gROOT->ProcessLine("Process_Uncorrected_FourierFit()");
    
    // Then create comparison plots
    std::cout << "\nCreating comparison plots..." << std::endl;
    gROOT->ProcessLine("CompareUncorrectedVsTemplateFit()");
}
