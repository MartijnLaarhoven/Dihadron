void run_zyam() {
    gROOT->LoadMacro("Process_ZYAM_Refit.cxx+");
    Process_LMHM_from_ZYAM("LHC25ae_pass2_604826", 0, 20, 80, 100);
}
