#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <iomanip>
#include <iostream>

void rel_diff() {

  TFile* f_eff = new TFile("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/eff_2016.root");
  TFile* f_wei = new TFile("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/wei_2016.root");

  TEfficiency* eff_2016 = (TEfficiency*)f_eff->Get("eff_den_clone");
  TEfficiency* wei_2016 = (TEfficiency*)f_wei->Get("eff_wei_den_clone");

  double efficiency_2016[8];
  double wei_efficiency_2016[8];

  for(int i=0; i<8; i++) {

    efficiency_2016[i] = eff_2016->GetEfficiency(i+1);
    //cout << "eff - " << efficiency_2016[i] << endl;
    wei_efficiency_2016[i]= wei_2016->GetEfficiency(i+1);
    //cout << "wei - " << wei_efficiency_2016[i] << endl;

  }

}
