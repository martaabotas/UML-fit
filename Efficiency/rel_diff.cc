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

  TString wei_variable("bPt");
  
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
  Double_t q2Bins[] = {1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16};
  Int_t n_bins = 8;

  TCanvas c1;
  
  for(int i=0; i<n_bins; i++) {

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

  //PLOT RELATIVE DIFERENCE
  
  TGraph* gr1 = new TGraph(n_bins,q2Bins,relative_diff);
  gr1->SetName("gr1");
  gr1->Draw("AP");
  gr1->SetMarkerColor(kGreen); //CHANGE
  gr1->SetMarkerStyle(kFullDotLarge);
  gr1->SetTitle("");
  gr1->GetYaxis()->SetTitle("Relative Difference");
  gr1->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");
  
  auto leg1 = new TLegend(0.65,0.8,0.9,0.88);
  leg1->Draw();
  leg1->AddEntry("gr1","weights_"+wei_variable,"p");
  leg1->SetBorderSize(0);

  c1.SaveAs(Form("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/rel_diff_%i_bPt.gif",year));

  TFile* f_rel = new TFile(Form("./rootfiles/rel_diff_%i.root",year),"RECREATE");
  f_rel->cd();
  gr1->Write();
  f_rel->Close();
  
  //PLOT EFFICIENCIES IN MULTIGRAPH

  TMultiGraph *mg = new TMultiGraph();
  
  TGraph* gr2 = new TGraph(n_bins,q2Bins,efficiency);
  gr2->SetName("gr2");
  gr2->Draw("AP");
  gr2->SetMarkerColor(kRed);
  gr2->SetMarkerStyle(kFullDotLarge);
  gr2->SetTitle("");
  
  TGraph* gr3 = new TGraph(n_bins,q2Bins,wei_efficiency);
  gr3->SetName("gr3");
  gr3->Draw("AP");
  gr3->SetMarkerColor(kGreen);
  gr3->SetMarkerStyle(kFullDotLarge);
  gr3->SetTitle("Weighted Efficiency");
  
  mg->Add(gr2);
  mg->Add(gr3);
  
  auto leg2 = new TLegend(0.65,0.7,0.9,0.88);
  leg2->Draw();
  leg2->AddEntry("gr2","Efficiency","p");
  leg2->AddEntry("gr3","Weighted Efficiency","p");
  leg2->SetBorderSize(0);
  
  TCanvas c2;
  c2.cd();

  mg->Draw("AP");
  leg2->Draw();
  mg->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");

  c2.SaveAs(Form("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/all_%i_bPt.gif",year));
  
  TFile* f_cp = new TFile(Form("./rootfiles/cp_%i.root",year),"RECREATE");
  f_cp->cd();
  mg->Write();
  f_cp->Close();
  
}
