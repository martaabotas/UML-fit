#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TGraph.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <TAttMarker.h>
#include <TChain.h>
#include <TGraphErrors.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <TH1F.h>
#include <cstdlib>

void rel_diff() {

  /*TString input_file_eff = Form("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/eff_%i.root",year,year);
  TString input_file_wei = Form("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/wei_%i.root",year,year);

  TFile* f_eff = new TFile(input_file_eff);
  TFile* f_wei = new TFile(input_file_wei);*/

  TFile* f_eff = new TFile("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/eff_2016.root");
  TFile* f_wei = new TFile("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/wei_2016.root");

  TEfficiency* eff_2016 = (TEfficiency*)f_eff->Get("eff_den_clone");
  TEfficiency* wei_2016 = (TEfficiency*)f_wei->Get("eff_wei_den_clone");

  Double_t efficiency_2016[8];
  Double_t wei_efficiency_2016[8];
  Double_t diff[8];
  Double_t relative_diff[8];
  Double_t bins[] = {0, 1, 2, 3, 4, 5, 6, 7};
  Int_t n = 8;

  TCanvas *c1 = new TCanvas("c1","Relative Difference",200,10,500,300);
  //TCanvas *c2 = new TCanvas("c2","Weighted Efficiency",200,10,500,300);
  
  for(int i=0; i<n; i++) {

    efficiency_2016[i] = eff_2016->GetEfficiency(i+1);
    cout << "eff - " << efficiency_2016[i] << endl;
    wei_efficiency_2016[i]= wei_2016->GetEfficiency(i+1);
    cout << "wei - " << wei_efficiency_2016[i] << endl;
    diff[i] = abs(efficiency_2016[i]-wei_efficiency_2016[i]);
    cout << "diff - " << diff[i] << endl;
    relative_diff[i] = diff[i]/efficiency_2016[i];
    cout << "rel dif - " << relative_diff[i]*100 << " %" << endl;
    cout << " " << endl;
    
  }
  
  TGraph* gr1 = new TGraph(n,bins,relative_diff);
  gr1->Draw("AP");
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(kFullDotLarge);
  gr1->SetTitle("Efficiency relative difference"); //eff and wei eff diff
  gr1->GetXaxis()->SetTitle("q^{2} bins");
  gr1->GetYaxis()->SetTitle("Relative difference");
  c1->SaveAs("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/relative_diff_2016.gif");

}
