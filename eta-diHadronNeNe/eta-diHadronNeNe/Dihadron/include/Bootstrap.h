#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

#include <vector>
#include <TMath.h>
#include <iostream>

double epsilon = 1e-10;

bool isZero(double x) {
  return (x < epsilon && x > -epsilon);
}

//Passes values and error from the fit. Error array is empty and reshaped by the second function in this file. TotalOverSubSample is 1.
//the arrays are passed individualy for the three observables
void CalculateBootstrapError(const std::vector<std::vector<double>>& ValueArray, const std::vector<std::vector<double>>& ValueErrorArray, std::vector<double>& ErrorArray, Double_t TotalOverSubSample){
    int Nsample = ValueArray.size(); //100
    int Nbin = ValueArray[0].size(); //1
    std::vector<int> Count;
    std::vector<double> Mean;
    std::vector<double> SumWeight;
    std::vector<std::vector<bool>> MaskArray;
    Count.resize(Nbin);
    Mean.resize(Nbin);
    SumWeight.resize(Nbin);
    MaskArray.resize(Nsample); //Do we wish to compute the error here?
    // first calculate the mean and standard deviation for each bin
    for(int isample=0;isample<Nsample;isample++){
        MaskArray[isample].resize(Nbin);
        for(int ibin=0;ibin<Nbin;ibin++){ //fit failed scneario
            // Skip samples marked as failed (-999)
            if (ValueArray[isample][ibin] < -900) {
                MaskArray[isample][ibin] = false;
                continue;
            }
            if (ValueArray[isample][ibin]==-1 && ValueErrorArray[isample][ibin]==10) {
                MaskArray[isample][ibin] = false;
                continue;
            }
            if (isZero(ValueArray[isample][ibin]) || isZero(ValueErrorArray[isample][ibin])) { //Value is 0, i.e. turned off
                MaskArray[isample][ibin] = false;
                continue;
            }
            if (ValueArray[isample][ibin]!=ValueArray[isample][ibin] || ValueErrorArray[isample][ibin]!=ValueErrorArray[isample][ibin]) { //not sure what this is?
                MaskArray[isample][ibin] = false;
                continue;
            }
            MaskArray[isample][ibin] = true;
        }
    }
    ///If only some fit samples failed, here is the reason
    for(int ibin=0;ibin<Nbin;ibin++){ //just once
        SumWeight[ibin]=0;
        Mean[ibin]=0;
        for(int isample=0;isample<Nsample;isample++){ //100 times
            if (!MaskArray[isample][ibin]) continue; //do we wish to compute the error here?
            Count[ibin]++;
            Mean[ibin] += ValueArray[isample][ibin] * (1./(ValueErrorArray[isample][ibin]*ValueErrorArray[isample][ibin])); //this average is waighted by inverse error squared (avg(x)/var(x))
            // Printf("ValueArray[%d][%d]=%f, ValueErrorArray[%d][%d]=%f",j,i,ValueArray[j][i],j,i,ValueErrorArray[j][i]);
            SumWeight[ibin] += (1./(ValueErrorArray[isample][ibin]*ValueErrorArray[isample][ibin])); //sum of weights 1/var(x)
        }
        if(Count[ibin]>0) //this will always be true except if all sample fit failed
            Mean[ibin] = Mean[ibin]/SumWeight[ibin]; //here we multiply by the total variance: sum(1/var(xi))*sum(avg(xi)/sum(var(xi)))
    }
    for (Int_t isample = 0; isample < Nsample; isample++) //100 times
    {
        for (Int_t ibin = 0; ibin < Nbin; ibin++) //once
        {
            if (!MaskArray[isample][ibin]) continue; //if we don't wish to fit
            ErrorArray[ibin] += (ValueArray[isample][ibin] - Mean[ibin]) * (ValueArray[isample][ibin] - Mean[ibin]); //compute (xi - avg(x))^2 = std^2*N
        }
    }
    for (Int_t ibin = 0; ibin < Nbin; ibin++) //once
    {
        if (Count[ibin] > 2){ //if at least 3 samples succeeded
            ErrorArray[ibin] = TMath::Sqrt(ErrorArray[ibin] / (Count[ibin] - 1)); //here we compute the std
        }
        else {
            ErrorArray[ibin] = 10; //set the total error to 10 if two or less samples succeeded
            Printf("Bin %d has only %d samples, set error to 10",ibin,Count[ibin]);
        }
        // Printf("ErrorArray[%d]=%f",i,ErrorArray[i]);
    }

    // remove samples beyond 3 sigma
    for(int ibin=0;ibin<Nbin;ibin++) {
        for(int isample=0;isample<Nsample;isample++) {
            if (!MaskArray[isample][ibin]) continue;
            if (TMath::Abs(ValueArray[isample][ibin] - Mean[ibin]) > 3 * ErrorArray[ibin]) {
                MaskArray[isample][ibin] = false;
            }
        }
    }

    // re-calculate the mean and standard deviation for each bin
    for(int ibin=0;ibin<Nbin;ibin++){
        SumWeight[ibin]=0;
        Mean[ibin]=0;
        Count[ibin]=0;
        ErrorArray[ibin]=0;
        for(int isample=0;isample<Nsample;isample++){
            if (!MaskArray[isample][ibin]) continue;
            Count[ibin]++;
            double err = ValueErrorArray[isample][ibin];
            Mean[ibin] += ValueArray[isample][ibin] * (1./(ValueErrorArray[isample][ibin]*ValueErrorArray[isample][ibin]));
            SumWeight[ibin] += (1./(ValueErrorArray[isample][ibin]*ValueErrorArray[isample][ibin]));
        }
        if(Count[ibin]>0)
            Mean[ibin] = Mean[ibin]/SumWeight[ibin];
    }
    for (Int_t isample = 0; isample < Nsample; isample++)
    {
        for (Int_t ibin = 0; ibin < Nbin; ibin++)
        {
            if (!MaskArray[isample][ibin]) continue;
            ErrorArray[ibin] += (ValueArray[isample][ibin] - Mean[ibin]) * (ValueArray[isample][ibin] - Mean[ibin]);
        }
    }
    for (Int_t ibin = 0; ibin < Nbin; ibin++)
    {
        if (Count[ibin] > 2){
            ErrorArray[ibin] = TMath::Sqrt(ErrorArray[ibin] / (Count[ibin] - 1));
        }
        else {
            ErrorArray[ibin] = 10;
            Printf("Bin %d has only %d samples, set error to 10",ibin,Count[ibin]);
        }
        // Printf("ErrorArray[%d]=%f",i,ErrorArray[i]);
    }
    // SEM = SD/sqrt(N)
    // because ErrorArray[ibin] above is the error of subsample, which is smaller than the total sample,
    for (Int_t ibin = 0; ibin < Nbin; ibin++)
    {
        ErrorArray[ibin] = ErrorArray[ibin] / TMath::Sqrt(TotalOverSubSample); //TotalOverSubSample = 1
    }
    
}

/*
Loop over the samples, and find the means (waighted by the inverse of their respective error)
Use them to measure the std deviation of all samples, rejected of less then 2 samples succeeded
Then remove samples beyond 3 sigma and recalculate the std deviation
Find the standard diviation again
*/


void ResizeValueArray(std::vector<std::vector<std::vector<double>>>& ValueArray,
    std::vector<std::vector<std::vector<double>>>& ValueErrorArray,
    std::vector<std::vector<double>>& ErrorArray,
    int Nobs=1, int NofSample=10, int Nbin=9){
    ValueArray.clear();
    ValueErrorArray.clear();
    ErrorArray.clear();
    ValueArray.resize(Nobs);
    ValueErrorArray.resize(Nobs);
    ErrorArray.resize(Nobs);    
    for(int i=0;i<Nobs;i++){
        ErrorArray[i].resize(Nbin); //3x1
        ValueArray[i].resize(NofSample); //3x100
        ValueErrorArray[i].resize(NofSample); //3x100
        for(int j=0;j<NofSample;j++)ValueArray[i][j].resize(Nbin); //3x100x1
        for(int j=0;j<NofSample;j++)ValueErrorArray[i][j].resize(Nbin); //3x100x1
    }

}


#endif // BOOTSTRAP_H
