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
  
  TCanvas c1;
  c1.cd();
  gr1->SetName("gr1");
  gr1->Draw("AP");
  gr1->SetMarkerColor(kGreen); //CHANGE
  gr1->SetMarkerStyle(kFullDotLarge);
  gr1->SetTitle("");
  c1.SaveAs(Form("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/rel_diff_%i_%s.gif",year,wei_variable.Data()));

  TFile* f_rel = new TFile(Form("./rootfiles/rel_diff_%i_%s.root",year,wei_variable.Data()),"RECREATE");
  f_rel->cd();
  gr1->Write();
  f_rel->Close();
  
  //PLOT NOMINAL EFFICIENCY

  TGraph* gr2 = new TGraph(n_bins,q2Bins,efficiency);
  
  TCanvas c2;
  c2.cd();
  gr2->SetName("gr2");
  gr2->Draw("AP");
  gr2->SetMarkerColor(kBlue);
  gr2->SetMarkerStyle(kFullDotLarge);
  gr2->SetTitle("");
  gr2->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");
  gr2->GetYaxis()->SetTitle("Nominal Efficiency");
  c2.SaveAs(Form("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/rootfiles/nomeff/nom_eff_%i.gif",year));

  TFile* f_egr = new TFile(Form("./rootfiles/nom_eff_%i.root",year),"RECREATE");
  f_egr->cd();
  gr2->Write();
  f_egr->Close();
    
  //PLOT WEIGHTED EFFICIENCY
  
  TGraph* gr3 = new TGraph(n_bins,q2Bins,wei_efficiency);
  
  TCanvas c3;
  c3.cd();
  gr3->SetName("gr3");
  gr3->Draw("AP");
  gr3->SetMarkerColor(kOrange); //CHANGE
  gr3->SetMarkerStyle(kFullDotLarge);
  gr3->SetTitle("");
  c3.SaveAs(Form("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/wei_eff_%i_%s.gif",year,wei_variable.Data()));

  TFile* f_wgr = new TFile(Form("./rootfiles/wei_eff_%i_%s.root",year,wei_variable.Data()),"RECREATE");
  f_wgr->cd();
  gr3->Write();
  f_wgr->Close();
    
}
