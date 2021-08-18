#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TGraph.h>
#include <iomanip>
#include <iostream>

void rel_diff() {

  TFile* f_eff = new TFile("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/eff_2016.root");
  TFile* f_wei = new TFile("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/wei_2016.root");

  TEfficiency* eff_2016 = (TEfficiency*)f_eff->Get("eff_den_clone");
  TEfficiency* wei_2016 = (TEfficiency*)f_wei->Get("eff_wei_den_clone");

  Double_t efficiency_2016[8];
  Double_t wei_efficiency_2016[8];
  
  TCanvas *c1 = new TCanvas("c1","Efficiency",200,10,500,300);
  TCanvas *c2 = new TCanvas("c2","Weighted Efficiency",200,10,500,300);
  
  for(int i=0; i<8; i++) {

    efficiency_2016[i] = eff_2016->GetEfficiency(i+1);
    cout << "eff - " << efficiency_2016[i] << endl;
    wei_efficiency_2016[i]= wei_2016->GetEfficiency(i+1);
    cout << "wei - " << wei_efficiency_2016[i] << endl;
    cout << " " << endl;

  }
  
  cout << "graph" << endl;
  
  TGraph* gr = new TGraph(8,efficiency_2016,wei_efficiency_2016);
  gr->Draw("AC*");
  gr->SaveAs("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/comparison.gif");

}
