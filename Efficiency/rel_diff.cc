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
  
  TEfficiency* eff = (TEfficiency*)f_eff->Get("eff_den_clone");
  TEfficiency* wei = (TEfficiency*)f_wei->Get("eff_wei_den_clone");
  
  Double_t efficiency[8];
  Double_t wei_efficiency[8];
  Double_t diff[8];
  Double_t relative_diff[8];
  
  Double_t bins[] = {0, 1, 2, 3, 4, 5, 6, 7};
  Int_t n = 8;

  TCanvas *c1 = new TCanvas("c1","All",200,10,500,300);
  
  TMultiGraph *mg = new TMultiGraph("mg","mg");
  
  for(int i=0; i<n; i++) {

    efficiency[i] = eff->GetEfficiency(i+1);
    cout << "efficiency - " << efficiency[i] << endl;
    wei_efficiency[i] = wei->GetEfficiency(i+1);
    cout << "weighted efficiency - " << wei_efficiency[i] << endl;
    diff[i] = abs(efficiency[i]-wei_efficiency[i]);
    cout << "diff - " << diff[i] << endl;
    relative_diff[i] = diff[i]/efficiency[i];
    cout << "rel dif - " << relative_diff[i]*100 << " %" << endl;
    cout << " " << endl;
    
  }
  
  TGraph* gr1 = new TGraph(n,bins,relative_diff);
  gr1->SetName("gr1");
  gr1->Draw("AP");
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(kFullDotLarge);
  gr1->SetTitle("Efficiency relative difference");
  
  TGraph* gr2 = new TGraph(n,bins,efficiency);
  gr2->SetName("gr2");
  gr2->Draw("AP");
  gr2->SetMarkerColor(kRed);
  gr2->SetMarkerStyle(kFullDotLarge);
  gr2->SetTitle("Efficiency");
  
  TGraph* gr3 = new TGraph(n,bins,wei_efficiency);
  gr3->SetName("gr3");
  gr3->Draw("AP");
  gr3->SetMarkerColor(kGreen);
  gr3->SetMarkerStyle(kFullDotLarge);
  gr3->SetTitle("Weighted Efficiency");
  
  mg->Add(gr1);
  mg->Add(gr2);
  mg->Add(gr3);
  mg->SetTitle("Efficiencies; q^{2} bins; Efficiency");
  mg->Draw("AP");
  c1->BuildLegend();
  c1->SaveAs(Form("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/rel_diff_%i.gif",year));

}
