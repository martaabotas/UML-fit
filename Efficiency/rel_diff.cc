#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TGraphPainter.h>
#include <TLine.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TMultiGraph.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <TAttMarker.h>
#include <TChain.h>
#include <TGraphErrors.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
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

  TCanvas c1;
  
  for(int i=0; i<n; i++) {

    efficiency[i] = eff->GetEfficiency(i+1);
    cout << "efficiency - " << efficiency[i] << endl;
    wei_efficiency[i] = wei->GetEfficiency(i+1);
    cout << "weighted efficiency - " << wei_efficiency[i] << endl;
    diff[i] = abs(efficiency[i]-wei_efficiency[i]);
    cout << "diff - " << diff[i] << endl;
    relative_diff[i] = diff[i]/efficiency[i];
    cout << "rel diff - " << relative_diff[i]*100 << " %" << endl;
    cout << " " << endl;
    
  }
  
  TGraph* gr1 = new TGraph(n,bins,relative_diff);
  gr1->SetName("gr1");
  gr1->Draw("AP");
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(kFullDotLarge);
  gr1->SetTitle("");
  gr1->GetYaxis()->SetTitle("Relative Difference");
  gr1->GetXaxis()->SetTitle("q^{2} bins");
  c1.SaveAs(Form("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/rel_diff_%i.gif",year));
  
  TMultiGraph *mg = new TMultiGraph();
  
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
  
  mg->Add(gr2);
  mg->Add(gr3);
  
  auto legend = new TLegend(0.65,0.7,0.9,0.88);
  legend->Draw();
  legend->AddEntry("gr2","Efficiency","p");
  legend->AddEntry("gr3","Weighted Efficiency","p");
  legend->SetBorderSize(0);
  
  TCanvas c2;
  c2.cd();

  mg->Draw("AP");
  legend->Draw();
  mg->GetYaxis()->SetTitle("Efficiency");
  mg->GetXaxis()->SetTitle("q^{2} bins");

  c2.SaveAs(Form("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/all_%i.gif",year));
  
}
