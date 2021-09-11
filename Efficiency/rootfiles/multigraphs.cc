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

void multigraphs(int year) {
  
  if(year < 2016 || year > 2018) {return;}

  //get efficiencies

  TFile* f_nom_eff = new TFile(Form("./nom_eff_%i.root",year));
  TFile* f_wei_bEta = new TFile(Form("./wei_eff_%i_bEta.root",year));
  TFile* f_wei_k = new TFile(Form("./wei_eff_%i_kstTrkpEta.root",year));
  TFile* f_wei_bPt = new TFile(Form("./wei_eff_%i_bPt.root",year));

  TGraph* nom_eff = (TGraph*)f_nom_eff->Get("gr2");
  TGraph* wei_bEta = (TGraph*)f_wei_bEta->Get("gr3");
  TGraph* wei_k = (TGraph*)f_wei_k->Get("gr3");
  TGraph* wei_bPt = (TGraph*)f_wei_bPt->Get("gr3");

  //plot multigraph

  TMultiGraph* mg1 = new TMultiGraph();

  mg1->Add(nom_eff);
  mg1->Add(wei_bEta);
  mg1->Add(wei_k);
  mg1->Add(wei_bPt);
  
  auto leg1 = new TLegend(0.65,0.7,0.9,0.88);
  leg1->Draw();
  leg1->AddEntry(nom_eff,"Unweighted", "p");
  leg1->AddEntry(wei_bEta, "Weight: bEta", "p");
  leg1->AddEntry(wei_k, "Weight: kstTrkpEta", "p");
  leg1->AddEntry(wei_bPt, "Weight: bPt", "p");
  leg1->SetBorderSize(0);

  TCanvas c1;
  c1.cd();

  mg1->Draw("AP");
  leg1->Draw();
  mg1->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");

  c1.SaveAs(Form("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/rootfiles/effs_%i.gif",year));

  //get relative differences

  TFile* f_rel_bEta = new TFile(Form("./rel_diff_%i_bEta.root",year));
  TFile* f_rel_k = new TFile(Form("./rel_diff_%i_kstTrkpEta.root",year));
  TFile* f_rel_bPt = new TFile(Form("./rel_diff_%i_bPt.root",year));

  TGraph* rel_bEta = (TGraph*)f_rel_bEta->Get("gr1");
  TGraph* rel_k = (TGraph*)f_rel_k->Get("gr1");
  TGraph* rel_bPt = (TGraph*)f_rel_bPt->Get("gr1");

  //plot multigraph

  TMultiGraph* mg2 = new TMultiGraph();

  mg2->Add(rel_bEta);
  mg2->Add(rel_k);
  mg2->Add(rel_bPt);
  
  auto leg2 = new TLegend(0.65,0.7,0.9,0.88);
  leg2->Draw();
  leg2->AddEntry(rel_bEta, "Weight: bEta", "p");
  leg2->AddEntry(rel_k, "Weight: kstTrkpEta", "p");
  leg2->AddEntry(rel_bPt, "Weight: bPt", "p");
  leg2->SetBorderSize(0);

  TCanvas c2;
  c2.cd();

  mg2->Draw("AP");
  leg2->Draw();
  mg2->GetYaxis()->SetTitle("Relative Difference");
  mg2->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");

  c2.SaveAs(Form("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/rootfiles/rel_diffs_%i.gif",year));
  
}
  



  
  
