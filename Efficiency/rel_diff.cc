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

void rel_diff(int year) {

  if(year < 2016 || year > 2018) {return;}

  TString input_file_eff = Form("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/eff_%i.root",year);
  TString input_file_wei = Form("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/wei_%i.root",year);

  TFile* f_eff = new TFile(input_file_eff);
  TFile* f_wei = new TFile(input_file_wei);

  /*TFile* f_eff = new TFile("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/eff_2016.root");
    TFile* f_wei = new TFile("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/wei_2016.root");*/

  TEfficiency* eff = (TEfficiency*)f_eff->Get("eff_den_clone"); //these are for 2016
  TEfficiency* wei = (TEfficiency*)f_wei->Get("eff_wei_den_clone");

  Double_t efficiency[8];
  Double_t wei_efficiency[8];
  Double_t diff[8];
  Double_t relative_diff[8];
  Double_t bins[] = {0, 1, 2, 3, 4, 5, 6, 7};
  Int_t n = 8;

  TCanvas *c1 = new TCanvas("c1","Relative Difference",200,10,500,300);
  
  for(int i=0; i<n; i++) {

    efficiency[i] = eff->GetEfficiency(i+1);
    cout << "eff - " << efficiency[i] << endl;
    wei_efficiency[i]= wei->GetEfficiency(i+1);
    cout << "wei - " << wei_efficiency[i] << endl;
    diff[i] = abs(efficiency[i]-wei_efficiency[i]);
    cout << "diff - " << diff[i] << endl;
    relative_diff[i] = diff[i]/efficiency[i];
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
  c1->SaveAs(Form("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/relative_diff_%i.gif",year));

}
