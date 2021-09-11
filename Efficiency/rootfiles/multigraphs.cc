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
  
  //break segmentation

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
  leg1->AddEntry("nom_eff","Nominal Efficiency", "p");
  leg1->AddEntry("wei_bETa", "weight: bPt", "p");
  leg1->AddEntry("wei_k", "weight: kstTrkpEta", "p");
  leg1->AddEntry("wei_bPt", "weight: bPt", "p");
  leg1->SetBorderSize(0);

  TCanvas c1;
  c1.cd();

  mg1->Draw("AP");
  leg1->Draw();
  mg1->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");

  c1.SaveAs(Form("/home/t3cms/u21mbotas/efficiency/UML-fit/Efficiency/rootfiles/effs_%i.gif",year));

}
  



  
  
