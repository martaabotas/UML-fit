#include <TH1F.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <TAttMarker.h>
#include <TChain.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <TH1D.h>
#include <TEfficiency.h>

double read_weights(TH1F* histo_variable, double var_value);
double getWeight(double var_value, TH1F* h_weight);

void efficiency2(int year){

  double q2Bins[] = {1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16};
  const int n_q2Bins = 8;

  if(year < 2016 || year > 2018){return;}

  // RECO 
  TString input_file_reco_jpsi = Form("/home/t3cms/mcarolina/samples/Anomalies/%iMC_JPSI.root",year);
  TString input_file_reco_psi = Form("/home/t3cms/mcarolina/samples/Anomalies/%iMC_PSI.root",year);
  TString input_file_reco_lmnr = Form("/home/t3cms/mcarolina/samples/Anomalies/%iMC_LMNR.root",year);
  TFile* f_reco_jpsi = new TFile(input_file_reco_jpsi);
  TFile* f_reco_psi = new TFile(input_file_reco_psi);
  TFile* f_reco_lmnr = new TFile(input_file_reco_lmnr);

  TTree* t_reco_jpsi = (TTree*)f_reco_jpsi->Get("ntuple");
  TTree* t_reco_psi = (TTree*)f_reco_psi->Get("ntuple");
  TTree* t_reco_lmnr = (TTree*)f_reco_lmnr->Get("ntuple");

  // GEN PU
  TString input_file_genpu_jpsi = Form("/home/t3cms/mcarolina/samples/Anomalies/%iGEN_MC_JPSI.root",year);
  TString input_file_genpu_psi = Form("/home/t3cms/mcarolina/samples/Anomalies/%iGEN_MC_PSI.root",year);
  TString input_file_genpu_lmnr = Form("/home/t3cms/mcarolina/samples/Anomalies/%iGEN_MC_LMNR.root",year);
  TFile* f_genpu_jpsi = new TFile(input_file_genpu_jpsi);
  TFile* f_genpu_psi = new TFile(input_file_genpu_psi);
  TFile* f_genpu_lmnr = new TFile(input_file_genpu_lmnr);

  TTree* t_genpu_jpsi = (TTree*)f_genpu_jpsi->Get("ntuple");
  TTree* t_genpu_psi = (TTree*)f_genpu_psi->Get("ntuple");
  TTree* t_genpu_lmnr = (TTree*)f_genpu_lmnr->Get("ntuple");

  // GEN
  TString input_file_mc_gen_jpsi = "/home/t3cms/mcarolina/samples/Anomalies/GEN_BFilter_B0JpsiKstar.root";
  TString input_file_mc_gen_psi = "/home/t3cms/mcarolina/samples/Anomalies/GEN_BFilter_B0PsiKstar.root";
  TString input_file_mc_gen_lmnr = "/home/t3cms/mcarolina/samples/Anomalies/GEN_BFilter_B0MuMuKstar_p*.root/ntuple";
  TFile* f_mc_gen_jpsi = new TFile(input_file_mc_gen_jpsi);
  TFile* f_mc_gen_psi = new TFile(input_file_mc_gen_psi);
  TChain* t_gen_lmnr = new TChain();

  TTree* t_gen_jpsi = (TTree*)f_mc_gen_jpsi->Get("ntuple");
  TTree* t_gen_psi = (TTree*)f_mc_gen_psi->Get("ntuple");
  t_gen_lmnr->Add(input_file_mc_gen_lmnr);

  TFile* weight_b[8];

  // MC vs SP weights (bEta)
  for(int w=0; w<n_q2Bins; w++) {
    weight_b[w] = new TFile(Form("/home/t3cms/mcarolina/samples/Anomalies/weights_%i_b%i.root",year,w));
  }

  TString weight_var("kstTrkpEta");
  TH1F* histo_wei_b[8];

  for(int wei=0; wei<n_q2Bins; wei++) {
    histo_wei_b[wei] = (TH1F*)weight_b[wei]->Get("weights_"+weight_var);
  }

  // EFFICIENCY
  TH1D* eff_num = new TH1D("eff_num", "eff_num", n_q2Bins, q2Bins);
  TH1D* eff_den = new TH1D("eff_den", "eff_den", n_q2Bins, q2Bins);

  TH1D* eff_wei_num = new TH1D("eff_wei_num", "eff_wei_num", n_q2Bins, q2Bins);
  TH1D* eff_wei_den = new TH1D("eff_wei_den", "eff_wei_den", n_q2Bins, q2Bins);

  double reco_mumuMass_jpsi;
  double reco_ctK_jpsi;
  double reco_ctL_jpsi;
  double reco_phi_jpsi;
  double reco_bEta_jpsi;
  double reco_runN_jpsi;

  double wei_variable;

  t_reco_jpsi->SetBranchAddress("mumuMass",&reco_mumuMass_jpsi);
  t_reco_jpsi->SetBranchAddress("cos_theta_k",&reco_ctK_jpsi);
  t_reco_jpsi->SetBranchAddress("cos_theta_l",&reco_ctL_jpsi);
  t_reco_jpsi->SetBranchAddress("phi_kst_mumu",&reco_phi_jpsi);
  t_reco_jpsi->SetBranchAddress("bEta",&reco_bEta_jpsi);
  t_reco_jpsi->SetBranchAddress("runN",&reco_runN_jpsi);

  t_reco_jpsi->SetBranchAddress(weight_var,&wei_variable);

  cout << "RECO JPSI" << endl;
  for(int evt = 0; evt < t_reco_jpsi->GetEntries()/100; evt++){
    t_reco_jpsi->GetEntry(evt);

    if( (year == 2016) && (reco_runN_jpsi >= 272007) && (reco_runN_jpsi <=278801) ){continue;}

    if( (pow(reco_mumuMass_jpsi,2) > q2Bins[4]) && (pow(reco_mumuMass_jpsi,2) < q2Bins[5]) ){

      eff_num->Fill(pow(reco_mumuMass_jpsi,2));
      eff_wei_num->Fill(pow(reco_mumuMass_jpsi,2),read_weights(histo_wei_b[4],wei_variable));
    } // 
  }

  double genpu_mumuMass_jpsi;
  double genpu_ctK_jpsi;
  double genpu_ctL_jpsi;
  double genpu_phi_jpsi;
  double genpu_runN_jpsi;

  t_genpu_jpsi->SetBranchAddress("genQ",&genpu_mumuMass_jpsi);
  t_genpu_jpsi->SetBranchAddress("gen_cos_theta_k",&genpu_ctK_jpsi);
  t_genpu_jpsi->SetBranchAddress("gen_cos_theta_l",&genpu_ctL_jpsi);
  t_genpu_jpsi->SetBranchAddress("gen_phi_kst_mumu",&genpu_phi_jpsi);
  t_genpu_jpsi->SetBranchAddress("runN",&genpu_runN_jpsi);

  t_genpu_jpsi->SetBranchAddress("gen"+weight_var,&wei_variable);

  double genpubEta_jpsi;
  double genpumupEta_jpsi;
  double genpumumEta_jpsi;
  double genpukstTrkpEta_jpsi;
  double genpukstTrkmEta_jpsi;
  double genpumupPt_jpsi;
  double genpumumPt_jpsi;
  double genpukstTrkpPt_jpsi;
  double genpukstTrkmPt_jpsi;

  t_genpu_jpsi->SetBranchAddress("genbEta",&genpubEta_jpsi);
  t_genpu_jpsi->SetBranchAddress("genmupEta",&genpumupEta_jpsi);
  t_genpu_jpsi->SetBranchAddress("genmumEta",&genpumumEta_jpsi);
  t_genpu_jpsi->SetBranchAddress("genkstTrkpEta",&genpukstTrkpEta_jpsi);
  t_genpu_jpsi->SetBranchAddress("genkstTrkmEta",&genpukstTrkmEta_jpsi);
  t_genpu_jpsi->SetBranchAddress("genmupPt",&genpumupPt_jpsi);
  t_genpu_jpsi->SetBranchAddress("genmumPt",&genpumumPt_jpsi);
  t_genpu_jpsi->SetBranchAddress("genkstTrkpPt",&genpukstTrkpPt_jpsi);
  t_genpu_jpsi->SetBranchAddress("genkstTrkmPt",&genpukstTrkmPt_jpsi);

  cout << "GEN PU JPSI" << endl;
  for(int evt = 0; evt < t_genpu_jpsi->GetEntries()/100; evt++){
    t_genpu_jpsi->GetEntry(evt);

    if( (year == 2016) && (genpu_runN_jpsi >= 272007) && (genpu_runN_jpsi <= 278801) ){continue;}

    if( (pow(genpu_mumuMass_jpsi,2) > q2Bins[4]) && (pow(genpu_mumuMass_jpsi,2) < q2Bins[5]) ){
      if(genpubEta_jpsi < 3){
        if(fabs(genpumupEta_jpsi)<2.5 && fabs(genpumumEta_jpsi)<2.5 &&
          fabs(genpukstTrkpEta_jpsi)<2.5 && fabs(genpukstTrkmEta_jpsi)<2.5 &&
          genpumupPt_jpsi>2.5 && genpumumPt_jpsi>2.5 &&
          genpukstTrkpPt_jpsi>0.4 && genpukstTrkmPt_jpsi>0.4){

	  eff_den->Fill(pow(genpu_mumuMass_jpsi,2));
	  eff_wei_den->Fill(pow(genpu_mumuMass_jpsi,2),read_weights(histo_wei_b[4],wei_variable));
        }
      }
    }
  }

  double reco_mumuMass_psi;
  double reco_ctK_psi;
  double reco_ctL_psi;
  double reco_phi_psi;
  double reco_bEta_psi;
  double reco_runN_psi;

  t_reco_psi->SetBranchAddress("mumuMass",&reco_mumuMass_psi);
  t_reco_psi->SetBranchAddress("cos_theta_k",&reco_ctK_psi);
  t_reco_psi->SetBranchAddress("cos_theta_l",&reco_ctL_psi);
  t_reco_psi->SetBranchAddress("phi_kst_mumu",&reco_phi_psi);
  t_reco_psi->SetBranchAddress("bEta",&reco_bEta_psi);
  t_reco_psi->SetBranchAddress("runN",&reco_runN_psi);

  t_reco_psi->SetBranchAddress(weight_var,&wei_variable);

  cout << "RECO PSI" << endl;
  for(int evt = 0; evt < t_reco_psi->GetEntries()/100; evt++){
    t_reco_psi->GetEntry(evt);

    if( (year == 2016) && (reco_runN_psi >= 272007) && (reco_runN_psi <= 278801) ){continue;}

    if( (pow(reco_mumuMass_psi,2) > q2Bins[6]) && (pow(reco_mumuMass_psi,2) < q2Bins[7]) ){

      eff_num->Fill(pow(reco_mumuMass_psi,2));    
      eff_wei_num->Fill(pow(reco_mumuMass_psi,2),read_weights(histo_wei_b[6],wei_variable));
    }
  }

  double genpu_mumuMass_psi;
  double genpu_ctK_psi;
  double genpu_ctL_psi;
  double genpu_phi_psi;
  double genpu_runN_psi;

  t_genpu_psi->SetBranchAddress("genQ",&genpu_mumuMass_psi);
  t_genpu_psi->SetBranchAddress("gen_cos_theta_k",&genpu_ctK_psi);
  t_genpu_psi->SetBranchAddress("gen_cos_theta_l",&genpu_ctL_psi);
  t_genpu_psi->SetBranchAddress("gen_phi_kst_mumu",&genpu_phi_psi);
  t_genpu_psi->SetBranchAddress("runN",&genpu_runN_psi);

  double genpubEta_psi;
  double genpumupEta_psi;
  double genpumumEta_psi;
  double genpukstTrkpEta_psi;
  double genpukstTrkmEta_psi;
  double genpumupPt_psi;
  double genpumumPt_psi;
  double genpukstTrkpPt_psi;
  double genpukstTrkmPt_psi;

  t_genpu_psi->SetBranchAddress("genbEta",&genpubEta_psi);
  t_genpu_psi->SetBranchAddress("genmupEta",&genpumupEta_psi);
  t_genpu_psi->SetBranchAddress("genmumEta",&genpumumEta_psi);
  t_genpu_psi->SetBranchAddress("genkstTrkpEta",&genpukstTrkpEta_psi);
  t_genpu_psi->SetBranchAddress("genkstTrkmEta",&genpukstTrkmEta_psi);
  t_genpu_psi->SetBranchAddress("genmupPt",&genpumupPt_psi);
  t_genpu_psi->SetBranchAddress("genmumPt",&genpumumPt_psi);
  t_genpu_psi->SetBranchAddress("genkstTrkpPt",&genpukstTrkpPt_psi);
  t_genpu_psi->SetBranchAddress("genkstTrkmPt",&genpukstTrkmPt_psi);

  t_genpu_psi->SetBranchAddress("gen"+weight_var,&wei_variable);

  cout << "GEN PU PSI" << endl;
  for(int evt = 0; evt < t_genpu_psi->GetEntries()/100; evt++){
    t_genpu_psi->GetEntry(evt);

    if( (year == 2016) && (genpu_runN_psi >= 272007) && (genpu_runN_psi <= 278801) ){continue;}

    if( (pow(genpu_mumuMass_psi,2) > q2Bins[6]) && (pow(genpu_mumuMass_psi,2) < q2Bins[7]) ){
      if(genpubEta_psi < 3){
        if(fabs(genpumupEta_psi)<2.5 && fabs(genpumumEta_psi)<2.5 &&
          fabs(genpukstTrkpEta_psi)<2.5 && fabs(genpukstTrkmEta_psi)<2.5 &&
          genpumupPt_psi>2.5 && genpumumPt_psi>2.5 &&
          genpukstTrkpPt_psi>0.4 && genpukstTrkmPt_psi>0.4){

          eff_den->Fill(pow(genpu_mumuMass_psi,2));
	  eff_wei_den->Fill(pow(genpu_mumuMass_psi,2),read_weights(histo_wei_b[6],wei_variable));
        }
      }
    }
  }

  double reco_mumuMass_lmnr;
  double reco_ctK_lmnr;
  double reco_ctL_lmnr;
  double reco_phi_lmnr;
  double reco_bEta_lmnr;
  double reco_runN_lmnr;

  t_reco_lmnr->SetBranchAddress("mumuMass",&reco_mumuMass_lmnr);
  t_reco_lmnr->SetBranchAddress("cos_theta_k",&reco_ctK_lmnr);
  t_reco_lmnr->SetBranchAddress("cos_theta_l",&reco_ctL_lmnr);
  t_reco_lmnr->SetBranchAddress("phi_kst_mumu",&reco_phi_lmnr);
  t_reco_lmnr->SetBranchAddress("bEta",&reco_bEta_lmnr);
  t_reco_lmnr->SetBranchAddress("runN",&reco_runN_lmnr);

  t_reco_lmnr->SetBranchAddress(weight_var,&wei_variable);

  cout << "RECO LMNR" << endl;
  for(int evt = 0; evt < t_reco_lmnr->GetEntries()/100; evt++){
    t_reco_lmnr->GetEntry(evt);

    if( (year == 2016) && (reco_runN_lmnr >= 272007) && (reco_runN_lmnr <= 278801) ){continue;}

    if( (pow(reco_mumuMass_lmnr,2) > q2Bins[0]) && (pow(reco_mumuMass_lmnr,2) < q2Bins[1]) ){

      eff_num->Fill(pow(reco_mumuMass_lmnr,2));
      eff_wei_num->Fill(pow(reco_mumuMass_lmnr,2),read_weights(histo_wei_b[0],wei_variable));
    }

    else if( (pow(reco_mumuMass_lmnr,2) > q2Bins[1]) && (pow(reco_mumuMass_lmnr,2) < q2Bins[2]) ){

      eff_num->Fill(pow(reco_mumuMass_lmnr,2));
      eff_wei_num->Fill(pow(reco_mumuMass_lmnr,2),read_weights(histo_wei_b[1],wei_variable));
    }

    else if( (pow(reco_mumuMass_lmnr,2) > q2Bins[2]) && (pow(reco_mumuMass_lmnr,2) < q2Bins[3]) ){

      eff_num->Fill(pow(reco_mumuMass_lmnr,2));
      eff_wei_num->Fill(pow(reco_mumuMass_lmnr,2),read_weights(histo_wei_b[2],wei_variable));
    }

    else if( (pow(reco_mumuMass_lmnr,2) > q2Bins[3]) && (pow(reco_mumuMass_lmnr,2) < q2Bins[4]) ){

      eff_num->Fill(pow(reco_mumuMass_lmnr,2));
      eff_wei_num->Fill(pow(reco_mumuMass_lmnr,2),read_weights(histo_wei_b[3],wei_variable));
    }

    else if( (pow(reco_mumuMass_lmnr,2) > q2Bins[5]) && (pow(reco_mumuMass_lmnr,2) < q2Bins[6]) ){

      eff_num->Fill(pow(reco_mumuMass_lmnr,2));
      eff_wei_num->Fill(pow(reco_mumuMass_lmnr,2),read_weights(histo_wei_b[5],wei_variable));
    }

    else if( (pow(reco_mumuMass_lmnr,2) > q2Bins[7]) && (pow(reco_mumuMass_lmnr,2) < q2Bins[8]) ){

      eff_num->Fill(pow(reco_mumuMass_lmnr,2));
      eff_wei_num->Fill(pow(reco_mumuMass_lmnr,2),read_weights(histo_wei_b[7],wei_variable));
    }

  }

  double genpu_mumuMass_lmnr;
  double genpu_ctK_lmnr;
  double genpu_ctL_lmnr;
  double genpu_phi_lmnr;
  double genpu_runN_lmnr;

  t_genpu_lmnr->SetBranchAddress("genQ",&genpu_mumuMass_lmnr);
  t_genpu_lmnr->SetBranchAddress("gen_cos_theta_k",&genpu_ctK_lmnr);
  t_genpu_lmnr->SetBranchAddress("gen_cos_theta_l",&genpu_ctL_lmnr);
  t_genpu_lmnr->SetBranchAddress("gen_phi_kst_mumu",&genpu_phi_lmnr);
  t_genpu_lmnr->SetBranchAddress("runN",&genpu_runN_lmnr);

  double genpubEta_lmnr;
  double genpumupEta_lmnr;
  double genpumumEta_lmnr;
  double genpukstTrkpEta_lmnr;
  double genpukstTrkmEta_lmnr;
  double genpumupPt_lmnr;
  double genpumumPt_lmnr;
  double genpukstTrkpPt_lmnr;
  double genpukstTrkmPt_lmnr;

  t_genpu_lmnr->SetBranchAddress("genbEta",&genpubEta_lmnr);
  t_genpu_lmnr->SetBranchAddress("genmupEta",&genpumupEta_lmnr);
  t_genpu_lmnr->SetBranchAddress("genmumEta",&genpumumEta_lmnr);
  t_genpu_lmnr->SetBranchAddress("genkstTrkpEta",&genpukstTrkpEta_lmnr);
  t_genpu_lmnr->SetBranchAddress("genkstTrkmEta",&genpukstTrkmEta_lmnr);
  t_genpu_lmnr->SetBranchAddress("genmupPt",&genpumupPt_lmnr);
  t_genpu_lmnr->SetBranchAddress("genmumPt",&genpumumPt_lmnr);
  t_genpu_lmnr->SetBranchAddress("genkstTrkpPt",&genpukstTrkpPt_lmnr);
  t_genpu_lmnr->SetBranchAddress("genkstTrkmPt",&genpukstTrkmPt_lmnr);

  t_genpu_lmnr->SetBranchAddress("gen"+weight_var,&wei_variable);

  cout << "GEN PU LMNR" << endl;
  for(int evt = 0; evt < t_genpu_lmnr->GetEntries()/100; evt++){
    t_genpu_lmnr->GetEntry(evt);

    if( (year == 2016) && (genpu_runN_lmnr >= 272007) && (genpu_runN_lmnr <= 278801) ){continue;}

    if( (pow(genpu_mumuMass_lmnr,2) > q2Bins[0]) && (pow(genpu_mumuMass_lmnr,2) < q2Bins[1]) ){
      if(genpubEta_lmnr < 3){
        if(fabs(genpumupEta_lmnr)<2.5 && fabs(genpumumEta_lmnr)<2.5 &&
          fabs(genpukstTrkpEta_lmnr)<2.5 && fabs(genpukstTrkmEta_lmnr)<2.5 &&
          genpumupPt_lmnr>2.5 && genpumumPt_lmnr>2.5 &&
          genpukstTrkpPt_lmnr>0.4 && genpukstTrkmPt_lmnr>0.4){

          eff_den->Fill(pow(genpu_mumuMass_lmnr,2));
	  eff_wei_den->Fill(pow(genpu_mumuMass_lmnr,2),read_weights(histo_wei_b[0],wei_variable));
        }
      }
    }

    else if( (pow(genpu_mumuMass_lmnr,2) > q2Bins[1]) && (pow(genpu_mumuMass_lmnr,2) < q2Bins[2]) ){
      if(genpubEta_lmnr < 3){
        if(fabs(genpumupEta_lmnr)<2.5 && fabs(genpumumEta_lmnr)<2.5 &&
          fabs(genpukstTrkpEta_lmnr)<2.5 && fabs(genpukstTrkmEta_lmnr)<2.5 &&
          genpumupPt_lmnr>2.5 && genpumumPt_lmnr>2.5 &&
          genpukstTrkpPt_lmnr>0.4 && genpukstTrkmPt_lmnr>0.4){

          eff_den->Fill(pow(genpu_mumuMass_lmnr,2));
          eff_wei_den->Fill(pow(genpu_mumuMass_lmnr,2),read_weights(histo_wei_b[1],wei_variable));
        }
      }
    }

    else if( (pow(genpu_mumuMass_lmnr,2) > q2Bins[2]) && (pow(genpu_mumuMass_lmnr,2) < q2Bins[3]) ){
      if(genpubEta_lmnr < 3){
        if(fabs(genpumupEta_lmnr)<2.5 && fabs(genpumumEta_lmnr)<2.5 &&
          fabs(genpukstTrkpEta_lmnr)<2.5 && fabs(genpukstTrkmEta_lmnr)<2.5 &&
          genpumupPt_lmnr>2.5 && genpumumPt_lmnr>2.5 &&
          genpukstTrkpPt_lmnr>0.4 && genpukstTrkmPt_lmnr>0.4){

          eff_den->Fill(pow(genpu_mumuMass_lmnr,2));
          eff_wei_den->Fill(pow(genpu_mumuMass_lmnr,2),read_weights(histo_wei_b[2],wei_variable));
        }
      }
    }

    else if( (pow(genpu_mumuMass_lmnr,2) > q2Bins[3]) && (pow(genpu_mumuMass_lmnr,2) < q2Bins[4]) ){
      if(genpubEta_lmnr < 3){
        if(fabs(genpumupEta_lmnr)<2.5 && fabs(genpumumEta_lmnr)<2.5 &&
          fabs(genpukstTrkpEta_lmnr)<2.5 && fabs(genpukstTrkmEta_lmnr)<2.5 &&
          genpumupPt_lmnr>2.5 && genpumumPt_lmnr>2.5 &&
          genpukstTrkpPt_lmnr>0.4 && genpukstTrkmPt_lmnr>0.4){

          eff_den->Fill(pow(genpu_mumuMass_lmnr,2));
          eff_wei_den->Fill(pow(genpu_mumuMass_lmnr,2),read_weights(histo_wei_b[3],wei_variable));
        }
      }
    }

    else if( (pow(genpu_mumuMass_lmnr,2) > q2Bins[5]) && (pow(genpu_mumuMass_lmnr,2) < q2Bins[6]) ){
      if(genpubEta_lmnr < 3){
        if(fabs(genpumupEta_lmnr)<2.5 && fabs(genpumumEta_lmnr)<2.5 &&
          fabs(genpukstTrkpEta_lmnr)<2.5 && fabs(genpukstTrkmEta_lmnr)<2.5 &&
          genpumupPt_lmnr>2.5 && genpumumPt_lmnr>2.5 &&
          genpukstTrkpPt_lmnr>0.4 && genpukstTrkmPt_lmnr>0.4){

          eff_den->Fill(pow(genpu_mumuMass_lmnr,2));
          eff_wei_den->Fill(pow(genpu_mumuMass_lmnr,2),read_weights(histo_wei_b[5],wei_variable));
        }
      }
    }

    else if( (pow(genpu_mumuMass_lmnr,2) > q2Bins[7]) && (pow(genpu_mumuMass_lmnr,2) < q2Bins[8]) ){
      if(genpubEta_lmnr < 3){
        if(fabs(genpumupEta_lmnr)<2.5 && fabs(genpumumEta_lmnr)<2.5 &&
          fabs(genpukstTrkpEta_lmnr)<2.5 && fabs(genpukstTrkmEta_lmnr)<2.5 &&
          genpumupPt_lmnr>2.5 && genpumumPt_lmnr>2.5 &&
          genpukstTrkpPt_lmnr>0.4 && genpukstTrkmPt_lmnr>0.4){

          eff_den->Fill(pow(genpu_mumuMass_lmnr,2));
          eff_wei_den->Fill(pow(genpu_mumuMass_lmnr,2),read_weights(histo_wei_b[7],wei_variable));
        }
      }
    }

  }

  // ACCEPTANCE
  TH1D* acc_num = new TH1D("acc_num", "acc_num", n_q2Bins, q2Bins);
  TH1D* acc_den = new TH1D("acc_den", "acc_den", n_q2Bins, q2Bins);

  double gen_mumuMass_jpsi;
  double gen_ctK_jpsi;
  double gen_ctL_jpsi;
  double gen_phi_jpsi;
  double gen_runN_jpsi;

  t_gen_jpsi->SetBranchAddress("genQ",&gen_mumuMass_jpsi);
  t_gen_jpsi->SetBranchAddress("cos_theta_k",&gen_ctK_jpsi);
  t_gen_jpsi->SetBranchAddress("cos_theta_l",&gen_ctL_jpsi);
  t_gen_jpsi->SetBranchAddress("phi_kst_mumu",&gen_phi_jpsi);
  t_gen_jpsi->SetBranchAddress("runN",&gen_runN_jpsi);

  double genbEta_jpsi;
  double genmupEta_jpsi;
  double genmumEta_jpsi;
  double genkstTrkpEta_jpsi;
  double genkstTrkmEta_jpsi;
  double genmupPt_jpsi;
  double genmumPt_jpsi;
  double genkstTrkpPt_jpsi;
  double genkstTrkmPt_jpsi;

  t_gen_jpsi->SetBranchAddress("genbEta",&genbEta_jpsi);
  t_gen_jpsi->SetBranchAddress("genmupEta",&genmupEta_jpsi);
  t_gen_jpsi->SetBranchAddress("genmumEta",&genmumEta_jpsi);
  t_gen_jpsi->SetBranchAddress("genkstTrkpEta",&genkstTrkpEta_jpsi);
  t_gen_jpsi->SetBranchAddress("genkstTrkmEta",&genkstTrkmEta_jpsi);
  t_gen_jpsi->SetBranchAddress("genmupPt",&genmupPt_jpsi);
  t_gen_jpsi->SetBranchAddress("genmumPt",&genmumPt_jpsi);
  t_gen_jpsi->SetBranchAddress("genkstTrkpPt",&genkstTrkpPt_jpsi);
  t_gen_jpsi->SetBranchAddress("genkstTrkmPt",&genkstTrkmPt_jpsi);

  cout << "GEN JPSI" << endl;
  for(int evt = 0; evt < t_gen_jpsi->GetEntries()/100; evt++){
    t_gen_jpsi->GetEntry(evt);

    if( (year == 2016) && (gen_runN_jpsi >= 272007) && (gen_runN_jpsi <= 278801) ){continue;};

    if( (pow(gen_mumuMass_jpsi,2) > q2Bins[4]) && (pow(gen_mumuMass_jpsi,2) < q2Bins[5]) ){
      if(genbEta_jpsi < 3){
        acc_den->Fill(pow(gen_mumuMass_jpsi,2));

        if(fabs(genmupEta_jpsi)<2.5 && fabs(genmumEta_jpsi)<2.5 &&
          fabs(genkstTrkpEta_jpsi)<2.5 && fabs(genkstTrkmEta_jpsi)<2.5 &&
          genmupPt_jpsi>2.5 && genmumPt_jpsi>2.5 &&
          genkstTrkpPt_jpsi>0.4 && genkstTrkmPt_jpsi>0.4){

          acc_num->Fill(pow(gen_mumuMass_jpsi,2));
        }
      }
    }
  }

  double gen_mumuMass_psi;
  double gen_ctK_psi;
  double gen_ctL_psi;
  double gen_phi_psi;
  double gen_runN_psi;

  t_gen_psi->SetBranchAddress("genQ",&gen_mumuMass_psi);
  t_gen_psi->SetBranchAddress("cos_theta_k",&gen_ctK_psi);
  t_gen_psi->SetBranchAddress("cos_theta_l",&gen_ctL_psi);
  t_gen_psi->SetBranchAddress("phi_kst_mumu",&gen_phi_psi);
  t_gen_psi->SetBranchAddress("runN",&gen_runN_psi);

  double genbEta_psi;
  double genmupEta_psi;
  double genmumEta_psi;
  double genkstTrkpEta_psi;
  double genkstTrkmEta_psi;
  double genmupPt_psi;
  double genmumPt_psi;
  double genkstTrkpPt_psi;
  double genkstTrkmPt_psi;

  t_gen_psi->SetBranchAddress("genbEta",&genbEta_psi);
  t_gen_psi->SetBranchAddress("genmupEta",&genmupEta_psi);
  t_gen_psi->SetBranchAddress("genmumEta",&genmumEta_psi);
  t_gen_psi->SetBranchAddress("genkstTrkpEta",&genkstTrkpEta_psi);
  t_gen_psi->SetBranchAddress("genkstTrkmEta",&genkstTrkmEta_psi);
  t_gen_psi->SetBranchAddress("genmupPt",&genmupPt_psi);
  t_gen_psi->SetBranchAddress("genmumPt",&genmumPt_psi);
  t_gen_psi->SetBranchAddress("genkstTrkpPt",&genkstTrkpPt_psi);
  t_gen_psi->SetBranchAddress("genkstTrkmPt",&genkstTrkmPt_psi);

  cout << "GEN PSI" << endl;
  for(int evt = 0; evt < t_gen_psi->GetEntries()/100; evt++){
    t_gen_psi->GetEntry(evt);

    if( (year == 2016) && (gen_runN_psi >= 272007) && (gen_runN_psi <= 278801) ){continue;}

    if( (pow(gen_mumuMass_psi,2) > q2Bins[6]) && (pow(gen_mumuMass_psi,2) < q2Bins[7]) ){
      if(genbEta_psi < 3){
        acc_den->Fill(pow(gen_mumuMass_psi,2));

        if(fabs(genmupEta_psi)<2.5 && fabs(genmumEta_psi)<2.5 &&
          fabs(genkstTrkpEta_psi)<2.5 && fabs(genkstTrkmEta_psi)<2.5 &&
          genmupPt_psi>2.5 && genmumPt_psi>2.5 &&
          genkstTrkpPt_psi>0.4 && genkstTrkmPt_psi>0.4){

          acc_num->Fill(pow(gen_mumuMass_psi,2));
        }
      }
    }
  }

  double gen_mumuMass_lmnr;
  double gen_ctK_lmnr;
  double gen_ctL_lmnr;
  double gen_phi_lmnr;
  double gen_runN_lmnr;

  t_gen_lmnr->SetBranchAddress("genQ",&gen_mumuMass_lmnr);
  t_gen_lmnr->SetBranchAddress("cos_theta_k",&gen_ctK_lmnr);
  t_gen_lmnr->SetBranchAddress("cos_theta_l",&gen_ctL_lmnr);
  t_gen_lmnr->SetBranchAddress("phi_kst_mumu",&gen_phi_lmnr);
  t_gen_lmnr->SetBranchAddress("runN",&gen_runN_lmnr);

  double genbEta_lmnr;
  double genmupEta_lmnr;
  double genmumEta_lmnr;
  double genkstTrkpEta_lmnr;
  double genkstTrkmEta_lmnr;
  double genmupPt_lmnr;
  double genmumPt_lmnr;
  double genkstTrkpPt_lmnr;
  double genkstTrkmPt_lmnr;

  t_gen_lmnr->SetBranchAddress("genbEta",&genbEta_lmnr);
  t_gen_lmnr->SetBranchAddress("genmupEta",&genmupEta_lmnr);
  t_gen_lmnr->SetBranchAddress("genmumEta",&genmumEta_lmnr);
  t_gen_lmnr->SetBranchAddress("genkstTrkpEta",&genkstTrkpEta_lmnr);
  t_gen_lmnr->SetBranchAddress("genkstTrkmEta",&genkstTrkmEta_lmnr);
  t_gen_lmnr->SetBranchAddress("genmupPt",&genmupPt_lmnr);
  t_gen_lmnr->SetBranchAddress("genmumPt",&genmumPt_lmnr);
  t_gen_lmnr->SetBranchAddress("genkstTrkpPt",&genkstTrkpPt_lmnr);
  t_gen_lmnr->SetBranchAddress("genkstTrkmPt",&genkstTrkmPt_lmnr);

  cout << "GEN LMNR" << endl;
  for(int evt = 0; evt < t_gen_lmnr->GetEntries()/100; evt++){
    t_gen_lmnr->GetEntry(evt);

    if( (year == 2016) && (gen_runN_lmnr >= 272007) && (gen_runN_lmnr <= 278801) ){continue;}

    if(genbEta_lmnr < 3){
      acc_den->Fill(pow(gen_mumuMass_lmnr,2));

      if(fabs(genmupEta_lmnr)<2.5 && fabs(genmumEta_lmnr)<2.5 &&
        fabs(genkstTrkpEta_lmnr)<2.5 && fabs(genkstTrkmEta_lmnr)<2.5 &&
        genmupPt_lmnr>2.5 && genmumPt_lmnr>2.5 &&
        genkstTrkpPt_lmnr>0.4 && genkstTrkmPt_lmnr>0.4){
   
        acc_num->Fill(pow(gen_mumuMass_lmnr,2));
      }
    }
  }

  TCanvas c_eff_num;
  c_eff_num.cd();
  eff_num->Draw();
  c_eff_num.SaveAs(Form("./Efficiency/histo_eff_num_%i.gif",year));

  TCanvas c_eff_den;
  c_eff_den.cd();
  eff_den->Draw();
  c_eff_den.SaveAs(Form("./Efficiency/histo_eff_den_%i.gif",year));

  TCanvas c_eff_wei_num;
  c_eff_wei_num.cd();
  eff_wei_num->Draw();
  c_eff_wei_num.SaveAs(Form("./Efficiency/histo_eff_wei_num_%i.gif",year));

  TCanvas c_eff_wei_den;
  c_eff_wei_den.cd();
  eff_wei_den->Draw();
  c_eff_wei_den.SaveAs(Form("./Efficiency/histo_eff_wei_den_%i.gif",year));

  TCanvas c_acc_num;
  c_acc_num.cd();
  acc_num->Draw();
  c_acc_num.SaveAs(Form("./Efficiency/histo_acc_num_%i.gif",year));

  TCanvas c_acc_den;
  c_acc_den.cd();
  acc_den->Draw();
  c_acc_den.SaveAs(Form("./Efficiency/histo_acc_den_%i.gif",year));

  TEfficiency* efficiency = new TEfficiency(*eff_num, *eff_den);
  TEfficiency* weighted = new TEfficiency(*eff_wei_num, *eff_wei_den);
  TEfficiency* acceptance = new TEfficiency(*acc_num, *acc_den);

  TH1D* eff_x_acc_num = new TH1D("eff_x_acc_num", "eff_x_acc_num", n_q2Bins, q2Bins);
  TH1D* eff_x_acc_den = new TH1D("eff_x_acc_den", "eff_x_acc_den", n_q2Bins, q2Bins);

  eff_x_acc_num->Multiply(eff_num,acc_num);
  eff_x_acc_den->Multiply(eff_den,acc_den);

  TEfficiency* eff_x_acc = new TEfficiency(*eff_x_acc_num, *eff_x_acc_den);
  
  TCanvas c_eff;
  c_eff.cd();
  efficiency->SetTitle(Form("Efficiency - %i",year));
  efficiency->Draw("AP");  
  c_eff.SaveAs(Form("./Efficiency/efficiency_%i.gif",year));

  TCanvas c_wei;
  c_wei.cd();
  weighted->SetTitle(Form("Weighted efficiency - %i",year));
  weighted->Draw("AP");   
  c_wei.SaveAs(Form("./Efficiency/weighted_%i.gif",year));

  TCanvas c_acc;
  c_acc.cd();
  acceptance->SetTitle(Form("Acceptance - %i",year));
  acceptance->Draw("AP");  
  c_acc.SaveAs(Form("./Efficiency/acceptance_%i.gif",year));

  TCanvas c_eff_x_acc;
  c_eff_x_acc.cd();
  eff_x_acc->SetTitle(Form("Efficiency x Acceptance - %i",year));
  eff_x_acc->Draw("AP");
  c_eff_x_acc.SaveAs(Form("./Efficiency/eff_x_acc_%i.gif",year));

  TFile* f_eff = new TFile(Form("./Efficiency/eff_%i.root",year),"RECREATE");
  f_eff->cd();
  efficiency->Write();
  f_eff->Close();

  TFile* f_wei = new TFile(Form("./Efficiency/wei_%i.root",year),"RECREATE");
  f_wei->cd();
  weighted->Write();
  f_wei->Close();

  TFile* f_acc = new TFile(Form("./Efficiency/acc_%i.root",year),"RECREATE");
  f_acc->cd();
  acceptance->Write();
  f_acc->Close();

  TFile* f_eff_x_acc = new TFile(Form("./Efficiency/eff_x_acc_%i.root",year),"RECREATE");
  f_eff_x_acc->cd();
  eff_x_acc->Write();
  f_eff_x_acc->Close();

}

double read_weights(TH1F* histo_variable, double var_value){

  double weight;
  double variable_min;
  double variable_max;

  variable_min = histo_variable->GetXaxis()->GetXmin();
  variable_max = histo_variable->GetXaxis()->GetXmax();

  //if the event is not in the range its weight is 1.
  if(var_value>=variable_min && var_value<=variable_max){weight = getWeight(var_value,histo_variable);}
  else{weight = 1;}

  return weight;
}

double getWeight(double var_value, TH1F* h_weight){
  int bin = h_weight->FindBin(var_value);
  return h_weight->GetBinContent(bin);
}


