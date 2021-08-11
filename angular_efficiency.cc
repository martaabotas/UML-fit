#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <RooRealVar.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>

double DecayRate(double ctK, double ctL, double phi, double Fl, double P1, double P2, double P3, double P4p, double P5p, double P6p, double P8p);
double read_weights(TH1F* histo_variable, double var_value);
double getWeight(double var_value, TH1F* h_weight);

void angular_efficiency(int year){

  double q2Bins[] = {1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16};
  const int n_q2Bins = 8;

  // RECO 
  TString input_file_reco_jpsi = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iMC_JPSI.root",year,year);
  TString input_file_reco_psi = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iMC_PSI.root",year,year);
  TString input_file_reco_lmnr = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iMC_LMNR.root",year,year);
  TFile* f_reco_jpsi = new TFile(input_file_reco_jpsi);
  TFile* f_reco_psi = new TFile(input_file_reco_psi);
  TFile* f_reco_lmnr = new TFile(input_file_reco_lmnr);  

  TTree* t_reco_jpsi = (TTree*)f_reco_jpsi->Get("ntuple");
  TTree* t_reco_psi = (TTree*)f_reco_psi->Get("ntuple");
  TTree* t_reco_lmnr = (TTree*)f_reco_lmnr->Get("ntuple");

  // GEN PU
  TString input_file_genpu_jpsi = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iGEN_MC_JPSI.root",year,year);
  TString input_file_genpu_psi = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iGEN_MC_PSI.root",year,year);
  TString input_file_genpu_lmnr = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iGEN_MC_LMNR.root",year,year);
  TFile* f_genpu_jpsi = new TFile(input_file_genpu_jpsi);
  TFile* f_genpu_psi = new TFile(input_file_genpu_psi);
  TFile* f_genpu_lmnr = new TFile(input_file_genpu_lmnr);

  TTree* t_genpu_jpsi = (TTree*)f_genpu_jpsi->Get("ntuple");
  TTree* t_genpu_psi = (TTree*)f_genpu_psi->Get("ntuple");
  TTree* t_genpu_lmnr = (TTree*)f_genpu_lmnr->Get("ntuple");

  // GEN
  TString input_file_gen_jpsi = "/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN_NoFilter/newphi/GEN_BFilter_B0JpsiKstar.root";
  TString input_file_gen_psi = "/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN_NoFilter/newphi/GEN_BFilter_B0PsiKstar.root";
  TString input_file_gen_lmnr = "/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN_NoFilter/newphi/GEN_BFilter_B0MuMuKstar_p*.root/ntuple";
  TFile* f_gen_jpsi = new TFile(input_file_gen_jpsi);
  TFile* f_gen_psi = new TFile(input_file_gen_psi);
  TChain* t_gen_lmnr = new TChain();

  TTree* t_gen_jpsi = (TTree*)f_gen_jpsi->Get("ntuple");
  TTree* t_gen_psi = (TTree*)f_gen_psi->Get("ntuple");
  t_gen_lmnr->Add(input_file_gen_lmnr);

  // MC vs SP weights (bEta)
  TFile* weight_b0 = new TFile(Form("~/public/UML-fit/results/mc_validation_plots/weights/weights_%i_b0.root",year));
  TFile* weight_b1 = new TFile(Form("~/public/UML-fit/results/mc_validation_plots/weights/weights_%i_b1.root",year));
  TFile* weight_b2 = new TFile(Form("~/public/UML-fit/results/mc_validation_plots/weights/weights_%i_b2.root",year));
  TFile* weight_b3 = new TFile(Form("~/public/UML-fit/results/mc_validation_plots/weights/weights_%i_b3.root",year));
  TFile* weight_b4 = new TFile(Form("~/public/UML-fit/results/mc_validation_plots/weights/weights_%i_b4.root",year));
  TFile* weight_b5 = new TFile(Form("~/public/UML-fit/results/mc_validation_plots/weights/weights_%i_b5.root",year));
  TFile* weight_b6 = new TFile(Form("~/public/UML-fit/results/mc_validation_plots/weights/weights_%i_b6.root",year));
  TFile* weight_b7 = new TFile(Form("~/public/UML-fit/results/mc_validation_plots/weights/weights_%i_b7.root",year));

  TH1F* histo_wei_b0 = (TH1F*)weight_b0->Get("weights_bEta");
  TH1F* histo_wei_b1 = (TH1F*)weight_b1->Get("weights_bEta");
  TH1F* histo_wei_b2 = (TH1F*)weight_b2->Get("weights_bEta");
  TH1F* histo_wei_b3 = (TH1F*)weight_b3->Get("weights_bEta");
  TH1F* histo_wei_b4 = (TH1F*)weight_b4->Get("weights_bEta");
  TH1F* histo_wei_b5 = (TH1F*)weight_b5->Get("weights_bEta");
  TH1F* histo_wei_b6 = (TH1F*)weight_b6->Get("weights_bEta");
  TH1F* histo_wei_b7 = (TH1F*)weight_b7->Get("weights_bEta");

  // angular parameters and respective statistical errors
  // non-resonant values are from LHCb (need to wait for the unblinding of the angular analysis)
  double FL_data[] = {0.655, 0.756, 0.684, 0.645, 0.55501, 0.461, 0.52136, 0.352};
  double P1_data[] = {-0.617, 0.168, 0.088, -0.071, -0.015484, -0.460, -0.039087, -0.511};
  double P2_data[] = {-0.443, -0.191, 0.105, 0.207, -0.0012001, 0.411, 0.0059624, 0.396};
  double P3_data[] = {0.324, 0.049, -0.090, -0.068, 0.23983, -0.078, 0.46670, -0.000};
  double P4p_data[] = {-0.080, -0.435, -0.312, -0.574, -0.95271, -0.491, -0.77781, -0.626};  
  double P5p_data[] = {0.365, -0.150, -0.439, -0.583, -0.0063346, -0.622, 0.0022117, -0.714};
  double P6p_data[] = {-0.226, -0.155, -0.293, -0.155, 0.0018248, -0.193, 0.010996, 0.061};
  double P8p_data[] = {-0.366, 0.037, 0.166, -0.129, -0.21960, 0.018, -0.35748, 0.007};

  double FL_MC[] = {0.714, 0.813, 0.740, 0.618, 0.6000, 0.456, 0.6000, 0.370};
  double P1_MC[] = {0.004, -0.052, -0.113, -0.147, -0.198, -0.236, -0.199, -0.402};
  double P2_MC[] = {-0.388, -0.258, 0.205, 0.418, -0.0003, 0.461, -0.0004, 0.445};
  double P3_MC[] = {0.00242, 0.00088, 0.00325, 0.00198, 0.0004, 0.00115, 0.0007, 0.00080};
  double P4p_MC[] = {0.163, -0.522, -0.917, -1.029, -0.870, -1.098, -0.816, -1.178};
  double P5p_MC[] = {0.345, -0.310, -0.723, -0.861, 0.0012, -0.852, 0.0005, -0.761};
  double P6p_MC[] = {-0.00237, -0.00095, -0.00074, 0.00002, -0.0003, -0.00091, 0.0000, -0.00114};
  double P8p_MC[] = {0.00226, 0.00051, -0.00032, -0.00143, 0.0004, 0.00037, 0.0002, 0.00324};

  // EFFICIENCY
  double eff_num[n_q2Bins];
  double eff_den[n_q2Bins];

  double eff_num_wei[n_q2Bins];
  double eff_den_wei[n_q2Bins];

  double den_data_eff = 0;
  double den_mc_eff = 0;

  double reco_mumuMass_jpsi;
  double reco_ctK_jpsi;
  double reco_ctL_jpsi;
  double reco_phi_jpsi;
  double reco_bEta_jpsi;
  double reco_runN_jpsi;

  t_reco_jpsi->SetBranchAddress("mumuMass",&reco_mumuMass_jpsi);
  t_reco_jpsi->SetBranchAddress("cos_theta_k",&reco_ctK_jpsi);
  t_reco_jpsi->SetBranchAddress("cos_theta_l",&reco_ctL_jpsi);
  t_reco_jpsi->SetBranchAddress("phi_kst_mumu",&reco_phi_jpsi);
  t_reco_jpsi->SetBranchAddress("bEta",&reco_bEta_jpsi);
  t_reco_jpsi->SetBranchAddress("runN",&reco_runN_jpsi);

  cout << "RECO JPSI" << endl;
  for(int evt = 0; evt < t_reco_jpsi->GetEntries(); evt++){
    t_reco_jpsi->GetEntry(evt);

    if( (year == 2016) && (reco_runN_jpsi >= 272007) && (reco_runN_jpsi <=278801) ){continue;}

    if( (pow(reco_mumuMass_jpsi,2) > q2Bins[4]) && (pow(reco_mumuMass_jpsi,2) < q2Bins[5]) ){

      den_data_eff = DecayRate(reco_ctK_jpsi,reco_ctL_jpsi,reco_phi_jpsi,FL_data[4],P1_data[4],P2_data[4],P3_data[4],P4p_data[4],P5p_data[4],P6p_data[4],P8p_data[4]);
      den_mc_eff = DecayRate(reco_ctK_jpsi,reco_ctL_jpsi,reco_phi_jpsi,FL_MC[4],P1_MC[4],P2_MC[4],P3_MC[4],P4p_MC[4],P5p_MC[4],P6p_MC[4],P8p_MC[4]);

      eff_num[4] += den_data_eff/den_mc_eff;
      eff_num_wei[4] += (den_data_eff/den_mc_eff)*read_weights(histo_wei_b4,reco_bEta_jpsi);
    }
  }

  double num_data_eff = 0;
  double num_mc_eff = 0;

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
  for(int evt = 0; evt < t_genpu_jpsi->GetEntries(); evt++){
    t_genpu_jpsi->GetEntry(evt);

    if( (year == 2016) && (genpu_runN_jpsi >= 272007) && (genpu_runN_jpsi <= 278801) ){continue;}

    if( (pow(genpu_mumuMass_jpsi,2) > q2Bins[4]) && (pow(genpu_mumuMass_jpsi,2) < q2Bins[5]) ){
      if(genpubEta_jpsi < 3){
        if(fabs(genpumupEta_jpsi)<2.5 && fabs(genpumumEta_jpsi)<2.5 &&
          fabs(genpukstTrkpEta_jpsi)<2.5 && fabs(genpukstTrkmEta_jpsi)<2.5 &&
          genpumupPt_jpsi>2.5 && genpumumPt_jpsi>2.5 &&
          genpukstTrkpPt_jpsi>0.4 && genpukstTrkmPt_jpsi>0.4){

            num_data_eff = DecayRate(genpu_ctK_jpsi,genpu_ctL_jpsi,genpu_phi_jpsi,FL_data[4],P1_data[4],P2_data[4],P3_data[4],P4p_data[4],P5p_data[4],P6p_data[4],P8p_data[4]);
            num_mc_eff = DecayRate(genpu_ctK_jpsi,genpu_ctL_jpsi,genpu_phi_jpsi,FL_MC[4],P1_MC[4],P2_MC[4],P3_MC[4],P4p_MC[4],P5p_MC[4],P6p_MC[4],P8p_MC[4]);

            eff_den[4] += num_data_eff/num_mc_eff;
            eff_den_wei[4] += (num_data_eff/num_mc_eff)*read_weights(histo_wei_b4,genpubEta_jpsi);
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

  cout << "RECO PSI" << endl;
  for(int evt = 0; evt < t_reco_psi->GetEntries(); evt++){
    t_reco_psi->GetEntry(evt);

    if( (year == 2016) && (reco_runN_psi >= 272007) && (reco_runN_psi <= 278801) ){continue;}

    if( (pow(reco_mumuMass_psi,2) > q2Bins[6]) && (pow(reco_mumuMass_psi,2) < q2Bins[7]) ){
    
      den_data_eff = DecayRate(reco_ctK_psi,reco_ctL_psi,reco_phi_psi,FL_data[6],P1_data[6],P2_data[6],P3_data[6],P4p_data[6],P5p_data[6],P6p_data[6],P8p_data[6]);
      den_mc_eff = DecayRate(reco_ctK_psi,reco_ctL_psi,reco_phi_psi,FL_MC[6],P1_MC[6],P2_MC[6],P3_MC[6],P4p_MC[6],P5p_MC[6],P6p_MC[6],P8p_MC[6]);

      eff_num[6] += den_data_eff/den_mc_eff;
      eff_num_wei[6] += (den_data_eff/den_mc_eff)*read_weights(histo_wei_b6,reco_bEta_psi);
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

  cout << "GEN PU PSI" << endl;
  for(int evt = 0; evt < t_genpu_psi->GetEntries(); evt++){
    t_genpu_psi->GetEntry(evt);

    if( (year == 2016) && (genpu_runN_psi >= 272007) && (genpu_runN_psi <= 278801) ){continue;}

    if( (pow(genpu_mumuMass_psi,2) > q2Bins[6]) && (pow(genpu_mumuMass_psi,2) < q2Bins[7]) ){
      if(genpubEta_psi < 3){
        if(fabs(genpumupEta_psi)<2.5 && fabs(genpumumEta_psi)<2.5 &&
          fabs(genpukstTrkpEta_psi)<2.5 && fabs(genpukstTrkmEta_psi)<2.5 &&
          genpumupPt_psi>2.5 && genpumumPt_psi>2.5 &&
          genpukstTrkpPt_psi>0.4 && genpukstTrkmPt_psi>0.4){

            num_data_eff = DecayRate(genpu_ctK_psi,genpu_ctL_psi,genpu_phi_psi,FL_data[6],P1_data[6],P2_data[6],P3_data[6],P4p_data[6],P5p_data[6],P6p_data[6],P8p_data[6]);
            num_mc_eff = DecayRate(genpu_ctK_psi,genpu_ctL_psi,genpu_phi_psi,FL_MC[6],P1_MC[6],P2_MC[6],P3_MC[6],P4p_MC[6],P5p_MC[6],P6p_MC[6],P8p_MC[6]);

            eff_den[6] += num_data_eff/num_mc_eff;
            eff_den_wei[6] += (num_data_eff/num_mc_eff)*read_weights(histo_wei_b6,genpubEta_psi);
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

  cout << "RECO LMNR" << endl;
  for(int evt = 0; evt < t_reco_lmnr->GetEntries(); evt++){
    t_reco_lmnr->GetEntry(evt);

    if( (year == 2016) && (reco_runN_lmnr >= 272007) && (reco_runN_lmnr <= 278801) ){continue;}

    if( (pow(reco_mumuMass_lmnr,2) > q2Bins[0]) && (pow(reco_mumuMass_lmnr,2) < q2Bins[1]) ){

      den_data_eff = DecayRate(reco_ctK_lmnr,reco_ctL_lmnr,reco_phi_lmnr,FL_data[0],P1_data[0],P2_data[0],P3_data[0],P4p_data[0],P5p_data[0],P6p_data[0],P8p_data[0]);
      den_mc_eff = DecayRate(reco_ctK_lmnr,reco_ctL_lmnr,reco_phi_lmnr,FL_MC[0],P1_MC[0],P2_MC[0],P3_MC[0],P4p_MC[0],P5p_MC[0],P6p_MC[0],P8p_MC[0]);

      eff_num[0] += den_data_eff/den_mc_eff;
      eff_num_wei[0] += (den_data_eff/den_mc_eff)*read_weights(histo_wei_b0,reco_bEta_lmnr);
    }

    else if( (pow(reco_mumuMass_lmnr,2) > q2Bins[1]) && (pow(reco_mumuMass_lmnr,2) < q2Bins[2]) ){

      den_data_eff = DecayRate(reco_ctK_lmnr,reco_ctL_lmnr,reco_phi_lmnr,FL_data[1],P1_data[1],P2_data[1],P3_data[1],P4p_data[1],P5p_data[1],P6p_data[1],P8p_data[1]);
      den_mc_eff = DecayRate(reco_ctK_lmnr,reco_ctL_lmnr,reco_phi_lmnr,FL_MC[1],P1_MC[1],P2_MC[1],P3_MC[1],P4p_MC[1],P5p_MC[1],P6p_MC[1],P8p_MC[1]);

      eff_num[1] += den_data_eff/den_mc_eff;
      eff_num_wei[1] += (den_data_eff/den_mc_eff)*read_weights(histo_wei_b1,reco_bEta_lmnr);
    }    

    else if( (pow(reco_mumuMass_lmnr,2) > q2Bins[2]) && (pow(reco_mumuMass_lmnr,2) < q2Bins[3]) ){

      den_data_eff = DecayRate(reco_ctK_lmnr,reco_ctL_lmnr,reco_phi_lmnr,FL_data[2],P1_data[2],P2_data[2],P3_data[2],P4p_data[2],P5p_data[2],P6p_data[2],P8p_data[2]);
      den_mc_eff = DecayRate(reco_ctK_lmnr,reco_ctL_lmnr,reco_phi_lmnr,FL_MC[2],P1_MC[2],P2_MC[2],P3_MC[2],P4p_MC[2],P5p_MC[2],P6p_MC[2],P8p_MC[2]);

      eff_num[2] += den_data_eff/den_mc_eff;
      eff_num_wei[2] += (den_data_eff/den_mc_eff)*read_weights(histo_wei_b2,reco_bEta_lmnr);
    }   

    else if( (pow(reco_mumuMass_lmnr,2) > q2Bins[3]) && (pow(reco_mumuMass_lmnr,2) < q2Bins[4]) ){

      den_data_eff = DecayRate(reco_ctK_lmnr,reco_ctL_lmnr,reco_phi_lmnr,FL_data[3],P1_data[3],P2_data[3],P3_data[3],P4p_data[3],P5p_data[3],P6p_data[3],P8p_data[3]);
      den_mc_eff = DecayRate(reco_ctK_lmnr,reco_ctL_lmnr,reco_phi_lmnr,FL_MC[3],P1_MC[3],P2_MC[3],P3_MC[3],P4p_MC[3],P5p_MC[3],P6p_MC[3],P8p_MC[3]);

      eff_num[3] += den_data_eff/den_mc_eff;
      eff_num_wei[3] += (den_data_eff/den_mc_eff)*read_weights(histo_wei_b3,reco_bEta_lmnr);
    }

    else if( (pow(reco_mumuMass_lmnr,2) > q2Bins[5]) && (pow(reco_mumuMass_lmnr,2) < q2Bins[6]) ){

      den_data_eff = DecayRate(reco_ctK_lmnr,reco_ctL_lmnr,reco_phi_lmnr,FL_data[5],P1_data[5],P2_data[5],P3_data[5],P4p_data[5],P5p_data[5],P6p_data[5],P8p_data[5]);
      den_mc_eff = DecayRate(reco_ctK_lmnr,reco_ctL_lmnr,reco_phi_lmnr,FL_MC[5],P1_MC[5],P2_MC[5],P3_MC[5],P4p_MC[5],P5p_MC[5],P6p_MC[5],P8p_MC[5]);

      eff_num[5] += den_data_eff/den_mc_eff;
      eff_num_wei[5] += (den_data_eff/den_mc_eff)*read_weights(histo_wei_b5,reco_bEta_lmnr);
    }

    else if( (pow(reco_mumuMass_lmnr,2) > q2Bins[7]) && (pow(reco_mumuMass_lmnr,2) < q2Bins[8]) ){

      den_data_eff = DecayRate(reco_ctK_lmnr,reco_ctL_lmnr,reco_phi_lmnr,FL_data[7],P1_data[7],P2_data[7],P3_data[7],P4p_data[7],P5p_data[7],P6p_data[7],P8p_data[7]);
      den_mc_eff = DecayRate(reco_ctK_lmnr,reco_ctL_lmnr,reco_phi_lmnr,FL_MC[7],P1_MC[7],P2_MC[7],P3_MC[7],P4p_MC[7],P5p_MC[7],P6p_MC[7],P8p_MC[7]);

      eff_num[7] += den_data_eff/den_mc_eff;
      eff_num_wei[7] += (den_data_eff/den_mc_eff)*read_weights(histo_wei_b7,reco_bEta_lmnr);
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

  cout << "GEN PU LMNR" << endl;
  for(int evt = 0; evt < t_genpu_lmnr->GetEntries(); evt++){
    t_genpu_lmnr->GetEntry(evt);

    if( (year == 2016) && (genpu_runN_lmnr >= 272007) && (genpu_runN_lmnr <= 278801) ){continue;}

    if( (pow(genpu_mumuMass_lmnr,2) > q2Bins[0]) && (pow(genpu_mumuMass_lmnr,2) < q2Bins[1]) ){
      if(genpubEta_lmnr < 3){
        if(fabs(genpumupEta_lmnr)<2.5 && fabs(genpumumEta_lmnr)<2.5 &&
          fabs(genpukstTrkpEta_lmnr)<2.5 && fabs(genpukstTrkmEta_lmnr)<2.5 &&
          genpumupPt_lmnr>2.5 && genpumumPt_lmnr>2.5 &&
          genpukstTrkpPt_lmnr>0.4 && genpukstTrkmPt_lmnr>0.4){

            num_data_eff = DecayRate(genpu_ctK_lmnr,genpu_ctL_lmnr,genpu_phi_lmnr,FL_data[0],P1_data[0],P2_data[0],P3_data[0],P4p_data[0],P5p_data[0],P6p_data[0],P8p_data[0]);
            num_mc_eff = DecayRate(genpu_ctK_lmnr,genpu_ctL_lmnr,genpu_phi_lmnr,FL_MC[0],P1_MC[0],P2_MC[0],P3_MC[0],P4p_MC[0],P5p_MC[0],P6p_MC[0],P8p_MC[0]);

            eff_den[0] += num_data_eff/num_mc_eff;
            eff_den_wei[0] += (num_data_eff/num_mc_eff)*read_weights(histo_wei_b0,genpubEta_lmnr);
        }
      }
    }

    else if( (pow(genpu_mumuMass_lmnr,2) > q2Bins[1]) && (pow(genpu_mumuMass_lmnr,2) < q2Bins[2]) ){
      if(genpubEta_lmnr < 3){
        if(fabs(genpumupEta_lmnr)<2.5 && fabs(genpumumEta_lmnr)<2.5 &&
          fabs(genpukstTrkpEta_lmnr)<2.5 && fabs(genpukstTrkmEta_lmnr)<2.5 &&
          genpumupPt_lmnr>2.5 && genpumumPt_lmnr>2.5 &&
          genpukstTrkpPt_lmnr>0.4 && genpukstTrkmPt_lmnr>0.4){

            num_data_eff = DecayRate(genpu_ctK_lmnr,genpu_ctL_lmnr,genpu_phi_lmnr,FL_data[1],P1_data[1],P2_data[1],P3_data[1],P4p_data[1],P5p_data[1],P6p_data[1],P8p_data[1]);
            num_mc_eff = DecayRate(genpu_ctK_lmnr,genpu_ctL_lmnr,genpu_phi_lmnr,FL_MC[1],P1_MC[1],P2_MC[1],P3_MC[1],P4p_MC[1],P5p_MC[1],P6p_MC[1],P8p_MC[1]);

            eff_den[1] += num_data_eff/num_mc_eff;
            eff_den_wei[1] += (num_data_eff/num_mc_eff)*read_weights(histo_wei_b1,genpubEta_lmnr);
        }
      }
    }

    else if( (pow(genpu_mumuMass_lmnr,2) > q2Bins[2]) && (pow(genpu_mumuMass_lmnr,2) < q2Bins[3]) ){
      if(genpubEta_lmnr < 3){
        if(fabs(genpumupEta_lmnr)<2.5 && fabs(genpumumEta_lmnr)<2.5 &&
          fabs(genpukstTrkpEta_lmnr)<2.5 && fabs(genpukstTrkmEta_lmnr)<2.5 &&
          genpumupPt_lmnr>2.5 && genpumumPt_lmnr>2.5 &&
          genpukstTrkpPt_lmnr>0.4 && genpukstTrkmPt_lmnr>0.4){

            num_data_eff = DecayRate(genpu_ctK_lmnr,genpu_ctL_lmnr,genpu_phi_lmnr,FL_data[2],P1_data[2],P2_data[2],P3_data[2],P4p_data[2],P5p_data[2],P6p_data[2],P8p_data[2]);
            num_mc_eff = DecayRate(genpu_ctK_lmnr,genpu_ctL_lmnr,genpu_phi_lmnr,FL_MC[2],P1_MC[2],P2_MC[2],P3_MC[2],P4p_MC[2],P5p_MC[2],P6p_MC[2],P8p_MC[2]);

            eff_den[2] += num_data_eff/num_mc_eff;
            eff_den_wei[2] += (num_data_eff/num_mc_eff)*read_weights(histo_wei_b2,genpubEta_lmnr);
        }
      }
    }

    else if( (pow(genpu_mumuMass_lmnr,2) > q2Bins[3]) && (pow(genpu_mumuMass_lmnr,2) < q2Bins[4]) ){
      if(genpubEta_lmnr < 3){
        if(fabs(genpumupEta_lmnr)<2.5 && fabs(genpumumEta_lmnr)<2.5 &&
          fabs(genpukstTrkpEta_lmnr)<2.5 && fabs(genpukstTrkmEta_lmnr)<2.5 &&
          genpumupPt_lmnr>2.5 && genpumumPt_lmnr>2.5 &&
          genpukstTrkpPt_lmnr>0.4 && genpukstTrkmPt_lmnr>0.4){

            num_data_eff = DecayRate(genpu_ctK_lmnr,genpu_ctL_lmnr,genpu_phi_lmnr,FL_data[3],P1_data[3],P2_data[3],P3_data[3],P4p_data[3],P5p_data[3],P6p_data[3],P8p_data[3]);
            num_mc_eff = DecayRate(genpu_ctK_lmnr,genpu_ctL_lmnr,genpu_phi_lmnr,FL_MC[3],P1_MC[3],P2_MC[3],P3_MC[3],P4p_MC[3],P5p_MC[3],P6p_MC[3],P8p_MC[3]);

            eff_den[3] += num_data_eff/num_mc_eff;
            eff_den_wei[3] += (num_data_eff/num_mc_eff)*read_weights(histo_wei_b3,genpubEta_lmnr);
        }
      }
    }

    else if( (pow(genpu_mumuMass_lmnr,2) > q2Bins[5]) && (pow(genpu_mumuMass_lmnr,2) < q2Bins[6]) ){
      if(genpubEta_lmnr < 3){
        if(fabs(genpumupEta_lmnr)<2.5 && fabs(genpumumEta_lmnr)<2.5 &&
          fabs(genpukstTrkpEta_lmnr)<2.5 && fabs(genpukstTrkmEta_lmnr)<2.5 &&
          genpumupPt_lmnr>2.5 && genpumumPt_lmnr>2.5 &&
          genpukstTrkpPt_lmnr>0.4 && genpukstTrkmPt_lmnr>0.4){

            num_data_eff = DecayRate(genpu_ctK_lmnr,genpu_ctL_lmnr,genpu_phi_lmnr,FL_data[5],P1_data[5],P2_data[5],P3_data[5],P4p_data[5],P5p_data[5],P6p_data[5],P8p_data[5]);
            num_mc_eff = DecayRate(genpu_ctK_lmnr,genpu_ctL_lmnr,genpu_phi_lmnr,FL_MC[5],P1_MC[5],P2_MC[5],P3_MC[5],P4p_MC[5],P5p_MC[5],P6p_MC[5],P8p_MC[5]);

            eff_den[5] += num_data_eff/num_mc_eff;
            eff_den_wei[5] += (num_data_eff/num_mc_eff)*read_weights(histo_wei_b5,genpubEta_lmnr);
        }
      }
    }

    else if( (pow(genpu_mumuMass_lmnr,2) > q2Bins[7]) && (pow(genpu_mumuMass_lmnr,2) < q2Bins[8]) ){
      if(genpubEta_lmnr < 3){
        if(fabs(genpumupEta_lmnr)<2.5 && fabs(genpumumEta_lmnr)<2.5 &&
          fabs(genpukstTrkpEta_lmnr)<2.5 && fabs(genpukstTrkmEta_lmnr)<2.5 &&
          genpumupPt_lmnr>2.5 && genpumumPt_lmnr>2.5 &&
          genpukstTrkpPt_lmnr>0.4 && genpukstTrkmPt_lmnr>0.4){
  
            num_data_eff = DecayRate(genpu_ctK_lmnr,genpu_ctL_lmnr,genpu_phi_lmnr,FL_data[7],P1_data[7],P2_data[7],P3_data[7],P4p_data[7],P5p_data[7],P6p_data[7],P8p_data[7]);
            num_mc_eff = DecayRate(genpu_ctK_lmnr,genpu_ctL_lmnr,genpu_phi_lmnr,FL_MC[7],P1_MC[7],P2_MC[7],P3_MC[7],P4p_MC[7],P5p_MC[7],P6p_MC[7],P8p_MC[7]);

            eff_den[7] += num_data_eff/num_mc_eff;
            eff_den_wei[7] += (num_data_eff/num_mc_eff)*read_weights(histo_wei_b7,genpubEta_lmnr);
        }
      }
    }

  }

  // ACCEPTANCE
  double acc_num[n_q2Bins];
  double acc_den[n_q2Bins];

  double den_data_acc = 0;
  double den_mc_acc = 0;

  double num_data_acc = 0;
  double num_mc_acc = 0;

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
  for(int evt = 0; evt < t_gen_jpsi->GetEntries(); evt++){
    t_gen_jpsi->GetEntry(evt);

    if( (year == 2016) && (gen_runN_jpsi >= 272007) && (gen_runN_jpsi <= 278801) ){continue;};  

    if( (pow(gen_mumuMass_jpsi,2) > q2Bins[4]) && (pow(gen_mumuMass_jpsi,2) < q2Bins[5]) ){
      if(genbEta_jpsi < 3){
        num_data_acc = DecayRate(gen_ctK_jpsi,gen_ctL_jpsi,gen_phi_jpsi,FL_data[4],P1_data[4],P2_data[4],P3_data[4],P4p_data[4],P5p_data[4],P6p_data[4],P8p_data[4]);     
        num_mc_acc = DecayRate(gen_ctK_jpsi,gen_ctL_jpsi,gen_phi_jpsi,FL_MC[4],P1_MC[4],P2_MC[4],P3_MC[4],P4p_MC[4],P5p_MC[4],P6p_MC[4],P8p_MC[4]);  

        acc_den[4] += num_data_acc/num_mc_acc;

        if(fabs(genmupEta_jpsi)<2.5 && fabs(genmumEta_jpsi)<2.5 &&
          fabs(genkstTrkpEta_jpsi)<2.5 && fabs(genkstTrkmEta_jpsi)<2.5 &&
          genmupPt_jpsi>2.5 && genmumPt_jpsi>2.5 &&
          genkstTrkpPt_jpsi>0.4 && genkstTrkmPt_jpsi>0.4){

            den_data_acc = DecayRate(gen_ctK_jpsi,gen_ctL_jpsi,gen_phi_jpsi,FL_data[4],P1_data[4],P2_data[4],P3_data[4],P4p_data[4],P5p_data[4],P6p_data[4],P8p_data[4]);
            den_mc_acc = DecayRate(gen_ctK_jpsi,gen_ctL_jpsi,gen_phi_jpsi,FL_MC[4],P1_MC[4],P2_MC[4],P3_MC[4],P4p_MC[4],P5p_MC[4],P6p_MC[4],P8p_MC[4]);           
            acc_num[4] += den_data_acc/den_mc_acc;
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
  for(int evt = 0; evt < t_gen_psi->GetEntries(); evt++){
    t_gen_psi->GetEntry(evt);

    if( (year == 2016) && (gen_runN_psi >= 272007) && (gen_runN_psi <= 278801) ){continue;}

    if( (pow(gen_mumuMass_psi,2) > q2Bins[6]) && (pow(gen_mumuMass_psi,2) < q2Bins[7]) ){
      if(genbEta_psi < 3){
        num_data_acc = DecayRate(gen_ctK_psi,gen_ctL_psi,gen_phi_psi,FL_data[6],P1_data[6],P2_data[6],P3_data[6],P4p_data[6],P5p_data[6],P6p_data[6],P8p_data[6]);      
        num_mc_acc = DecayRate(gen_ctK_psi,gen_ctL_psi,gen_phi_psi,FL_MC[6],P1_MC[6],P2_MC[6],P3_MC[6],P4p_MC[6],P5p_MC[6],P6p_MC[6],P8p_MC[6]);

        acc_den[6] += num_data_acc/num_mc_acc;

        if(fabs(genmupEta_psi)<2.5 && fabs(genmumEta_psi)<2.5 &&
          fabs(genkstTrkpEta_psi)<2.5 && fabs(genkstTrkmEta_psi)<2.5 &&
          genmupPt_psi>2.5 && genmumPt_psi>2.5 &&
          genkstTrkpPt_psi>0.4 && genkstTrkmPt_psi>0.4){

            den_data_acc = DecayRate(gen_ctK_psi,gen_ctL_psi,gen_phi_psi,FL_data[6],P1_data[6],P2_data[6],P3_data[6],P4p_data[6],P5p_data[6],P6p_data[6],P8p_data[6]);
            den_mc_acc = DecayRate(gen_ctK_psi,gen_ctL_psi,gen_phi_psi,FL_MC[6],P1_MC[6],P2_MC[6],P3_MC[6],P4p_MC[6],P5p_MC[6],P6p_MC[6],P8p_MC[6]);

            acc_num[6] += den_data_acc/den_mc_acc;     
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
  for(int evt = 0; evt < t_gen_lmnr->GetEntries(); evt++){
    t_gen_lmnr->GetEntry(evt);

    if( (year == 2016) && (gen_runN_lmnr >= 272007) && (gen_runN_lmnr <= 278801) ){continue;}

    if( (pow(gen_mumuMass_lmnr,2) > q2Bins[0]) && (pow(gen_mumuMass_lmnr,2) < q2Bins[1]) ){
      if(genbEta_lmnr < 3){
        num_data_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_data[0],P1_data[0],P2_data[0],P3_data[0],P4p_data[0],P5p_data[0],P6p_data[0],P8p_data[0]);
        num_mc_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_MC[0],P1_MC[0],P2_MC[0],P3_MC[0],P4p_MC[0],P5p_MC[0],P6p_MC[0],P8p_MC[0]);

        acc_den[0] += num_data_acc/num_mc_acc;

        if(fabs(genmupEta_lmnr)<2.5 && fabs(genmumEta_lmnr)<2.5 &&
          fabs(genkstTrkpEta_lmnr)<2.5 && fabs(genkstTrkmEta_lmnr)<2.5 &&
          genmupPt_lmnr>2.5 && genmumPt_lmnr>2.5 &&
          genkstTrkpPt_lmnr>0.4 && genkstTrkmPt_lmnr>0.4){

            den_data_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_data[0],P1_data[0],P2_data[0],P3_data[0],P4p_data[0],P5p_data[0],P6p_data[0],P8p_data[0]);
            den_mc_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_MC[0],P1_MC[0],P2_MC[0],P3_MC[0],P4p_MC[0],P5p_MC[0],P6p_MC[0],P8p_MC[0]);

            acc_num[0] += den_data_acc/den_mc_acc;
        }
      }
    }

    else if( (pow(gen_mumuMass_lmnr,2) > q2Bins[1]) && (pow(gen_mumuMass_lmnr,2) < q2Bins[2]) ){
      if(genbEta_lmnr < 3){
        num_data_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_data[1],P1_data[1],P2_data[1],P3_data[1],P4p_data[1],P5p_data[1],P6p_data[1],P8p_data[1]);
        num_mc_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_MC[1],P1_MC[1],P2_MC[1],P3_MC[1],P4p_MC[1],P5p_MC[1],P6p_MC[1],P8p_MC[1]);

        acc_den[1] += num_data_acc/num_mc_acc;

        if(fabs(genmupEta_lmnr)<2.5 && fabs(genmumEta_lmnr)<2.5 &&
          fabs(genkstTrkpEta_lmnr)<2.5 && fabs(genkstTrkmEta_lmnr)<2.5 &&
          genmupPt_lmnr>2.5 && genmumPt_lmnr>2.5 &&
          genkstTrkpPt_lmnr>0.4 && genkstTrkmPt_lmnr>0.4){

            den_data_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_data[1],P1_data[1],P2_data[1],P3_data[1],P4p_data[1],P5p_data[1],P6p_data[1],P8p_data[1]);
            den_mc_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_MC[1],P1_MC[1],P2_MC[1],P3_MC[1],P4p_MC[1],P5p_MC[1],P6p_MC[1],P8p_MC[1]);

            acc_num[1] += den_data_acc/den_mc_acc;
        }
      }
    }

    else if( (pow(gen_mumuMass_lmnr,2) > q2Bins[2]) && (pow(gen_mumuMass_lmnr,2) < q2Bins[3]) ){
      if(genbEta_lmnr < 3){
        num_data_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_data[2],P1_data[2],P2_data[2],P3_data[2],P4p_data[2],P5p_data[2],P6p_data[2],P8p_data[2]);
        num_mc_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_MC[2],P1_MC[2],P2_MC[2],P3_MC[2],P4p_MC[2],P5p_MC[2],P6p_MC[2],P8p_MC[2]);

        acc_den[2] += num_data_acc/num_mc_acc;

        if(fabs(genmupEta_lmnr)<2.5 && fabs(genmumEta_lmnr)<2.5 &&
          fabs(genkstTrkpEta_lmnr)<2.5 && fabs(genkstTrkmEta_lmnr)<2.5 &&
          genmupPt_lmnr>2.5 && genmumPt_lmnr>2.5 &&
          genkstTrkpPt_lmnr>0.4 && genkstTrkmPt_lmnr>0.4){

            den_data_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_data[2],P1_data[2],P2_data[2],P3_data[2],P4p_data[2],P5p_data[2],P6p_data[2],P8p_data[2]);
            den_mc_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_MC[2],P1_MC[2],P2_MC[2],P3_MC[2],P4p_MC[2],P5p_MC[2],P6p_MC[2],P8p_MC[2]);

            acc_num[2] += den_data_acc/den_mc_acc;
        }
      }
    }

    else if( (pow(gen_mumuMass_lmnr,2) > q2Bins[3]) && (pow(gen_mumuMass_lmnr,2) < q2Bins[4]) ){
      if(genbEta_lmnr < 3){
        num_data_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_data[3],P1_data[3],P2_data[3],P3_data[3],P4p_data[3],P5p_data[3],P6p_data[3],P8p_data[3]);
        num_mc_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_MC[3],P1_MC[3],P2_MC[3],P3_MC[3],P4p_MC[3],P5p_MC[3],P6p_MC[3],P8p_MC[3]);

        acc_den[3] += num_data_acc/num_mc_acc;
  
        if(fabs(genmupEta_lmnr)<2.5 && fabs(genmumEta_lmnr)<2.5 &&
          fabs(genkstTrkpEta_lmnr)<2.5 && fabs(genkstTrkmEta_lmnr)<2.5 &&
          genmupPt_lmnr>2.5 && genmumPt_lmnr>2.5 &&
          genkstTrkpPt_lmnr>0.4 && genkstTrkmPt_lmnr>0.4){

            den_data_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_data[3],P1_data[3],P2_data[3],P3_data[3],P4p_data[3],P5p_data[3],P6p_data[3],P8p_data[3]);
            den_mc_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_MC[3],P1_MC[3],P2_MC[3],P3_MC[3],P4p_MC[3],P5p_MC[3],P6p_MC[3],P8p_MC[3]);

            acc_num[3] += den_data_acc/den_mc_acc;
        }
      }
    }

    else if( (pow(gen_mumuMass_lmnr,2) > q2Bins[5]) && (pow(gen_mumuMass_lmnr,2) < q2Bins[6]) ){
      if(genbEta_lmnr < 3){
        num_data_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_data[5],P1_data[5],P2_data[5],P3_data[5],P4p_data[5],P5p_data[5],P6p_data[5],P8p_data[5]);
        num_mc_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_MC[5],P1_MC[5],P2_MC[5],P3_MC[5],P4p_MC[5],P5p_MC[5],P6p_MC[5],P8p_MC[5]);

        acc_den[5] += num_data_acc/num_mc_acc;

        if(fabs(genmupEta_lmnr)<2.5 && fabs(genmumEta_lmnr)<2.5 &&
          fabs(genkstTrkpEta_lmnr)<2.5 && fabs(genkstTrkmEta_lmnr)<2.5 &&
          genmupPt_lmnr>2.5 && genmumPt_lmnr>2.5 &&
          genkstTrkpPt_lmnr>0.4 && genkstTrkmPt_lmnr>0.4){

            den_data_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_data[5],P1_data[5],P2_data[5],P3_data[5],P4p_data[5],P5p_data[5],P6p_data[5],P8p_data[5]);
            den_mc_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_MC[5],P1_MC[5],P2_MC[5],P3_MC[5],P4p_MC[5],P5p_MC[5],P6p_MC[5],P8p_MC[5]);

            acc_num[5] += den_data_acc/den_mc_acc;
        }
      }
    }

    else if( (pow(gen_mumuMass_lmnr,2) > q2Bins[7]) && (pow(gen_mumuMass_lmnr,2) < q2Bins[8]) ){
      if(genbEta_lmnr < 3){
        num_data_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_data[7],P1_data[7],P2_data[7],P3_data[7],P4p_data[7],P5p_data[7],P6p_data[7],P8p_data[7]);
        num_mc_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_MC[7],P1_MC[7],P2_MC[7],P3_MC[7],P4p_MC[7],P5p_MC[7],P6p_MC[7],P8p_MC[7]);

        acc_den[7] += num_data_acc/num_mc_acc;
  
        if(fabs(genmupEta_lmnr)<2.5 && fabs(genmumEta_lmnr)<2.5 &&
          fabs(genkstTrkpEta_lmnr)<2.5 && fabs(genkstTrkmEta_lmnr)<2.5 &&
          genmupPt_lmnr>2.5 && genmumPt_lmnr>2.5 &&
          genkstTrkpPt_lmnr>0.4 && genkstTrkmPt_lmnr>0.4){

            den_data_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_data[7],P1_data[7],P2_data[7],P3_data[7],P4p_data[7],P5p_data[7],P6p_data[7],P8p_data[7]);
            den_mc_acc = DecayRate(gen_ctK_lmnr,gen_ctL_lmnr,gen_phi_lmnr,FL_MC[7],P1_MC[7],P2_MC[7],P3_MC[7],P4p_MC[7],P5p_MC[7],P6p_MC[7],P8p_MC[7]);

            acc_num[7] += den_data_acc/den_mc_acc;
        }
      }
    }

  }

  double efficiency[n_q2Bins];
  double weighted[n_q2Bins];
  double acceptance[n_q2Bins];
  double eff_x_acc[n_q2Bins];

  double efficiency_error[n_q2Bins];
  double weighted_error[n_q2Bins];
  double acceptance_error[n_q2Bins];
  double eff_x_acc_error[n_q2Bins];
  
  double q2Bins_half[n_q2Bins];
  double q2Bins_err[n_q2Bins];

  cout << '|' << setw(15) << "q2Bin" << '|' << setw(15) << "Eff error" << '|' << setw(15) << "Acc error" << '|' << setw(15) << "EffxAcc error" << '|' << endl;

  for(int i = 0; i < n_q2Bins; i++){
    q2Bins_half[i] = q2Bins[i] + 0.5* (q2Bins[i+1]-q2Bins[i]);
    q2Bins_err[i] = 0.5* (q2Bins[i+1]-q2Bins[i]);   

    efficiency[i] = eff_num[i]/eff_den[i];
    weighted[i] = eff_num_wei[i]/eff_den_wei[i];
    acceptance[i] = acc_num[i]/acc_den[i];    
    eff_x_acc[i] = efficiency[i]*acceptance[i];

    // binomial error
    efficiency_error[i] = sqrt( ( efficiency[i]*(1-efficiency[i]) )/eff_den[i] );
    weighted_error[i] = sqrt( ( weighted[i]*(1-weighted[i]) )/eff_den_wei[i] );
    acceptance_error[i] = sqrt( ( acceptance[i]*(1-acceptance[i]) )/acc_den[i] );
    eff_x_acc_error[i] = sqrt( pow(efficiency[i],2)*pow(acceptance_error[i],2) + pow(acceptance[i],2)*pow(efficiency_error[i],2) );

    cout << '|' << setw(15) << i << '|' << setw(15) << (efficiency_error[i]/efficiency[i])*100 << '|' << setw(15) << (acceptance_error[i]/acceptance[i])*100 << '|' << setw(15) << (eff_x_acc_error[i]/eff_x_acc[i])*100 << '|' << endl; 

  }

  TGraphErrors* gr_eff = new TGraphErrors(n_q2Bins,q2Bins_half,efficiency,q2Bins_err,efficiency_error);
  TGraphErrors* gr_wei = new TGraphErrors(n_q2Bins,q2Bins_half,weighted,q2Bins_err,weighted_error);
  TGraphErrors* gr_acc = new TGraphErrors(n_q2Bins,q2Bins_half,acceptance,q2Bins_err,acceptance_error);
  TGraphErrors* gr_eff_x_acc = new TGraphErrors(n_q2Bins,q2Bins_half,eff_x_acc,q2Bins_err,eff_x_acc_error);

  TCanvas c_eff;
  c_eff.cd();
  gr_eff->SetTitle(Form("Efficiency - %i",year));
  gr_eff->GetXaxis()->SetTitle("q^{2} [GeV]");
  gr_eff->SetMarkerStyle(8);
  gr_eff->SetMarkerColor(1);
  gr_eff->Draw("AP");
  c_eff.SaveAs(Form("~/public/UML-fit/Angular_efficiency/eff_%i.gif",year));
  c_eff.SaveAs(Form("~/public/UML-fit/Angular_efficiency/eff_%i.pdf",year));

  TCanvas c_wei;
  c_wei.cd();
  gr_wei->SetTitle(Form("Weighted efficiency - %i",year));
  gr_wei->GetXaxis()->SetTitle("q^{2} [GeV]");
  gr_wei->SetMarkerStyle(8);
  gr_wei->SetMarkerColor(1);
  gr_wei->Draw("AP");
  c_wei.SaveAs(Form("~/public/UML-fit/Angular_efficiency/wei_%i.gif",year));
  c_wei.SaveAs(Form("~/public/UML-fit/Angular_efficiency/wei_%i.pdf",year));

  TCanvas c_acc;
  c_acc.cd();
  gr_acc->SetTitle(Form("Acceptance - %i",year));
  gr_acc->GetXaxis()->SetTitle("q^{2} [GeV]");
  gr_acc->SetMarkerStyle(8);
  gr_acc->SetMarkerColor(1);
  gr_acc->Draw("AP");
  c_acc.SaveAs(Form("~/public/UML-fit/Angular_efficiency/acc_%i.gif",year));
  c_acc.SaveAs(Form("~/public/UML-fit/Angular_efficiency/acc_%i.pdf",year));

  TCanvas c_eff_x_acc;
  c_eff_x_acc.cd();
  gr_eff_x_acc->SetTitle(Form("Efficiency x Acceptance - %i",year));
  gr_eff_x_acc->GetXaxis()->SetTitle("q^{2} [GeV]");
  gr_eff_x_acc->GetYaxis()->SetTitle("#epsilon");
  gr_eff_x_acc->SetMarkerStyle(8);
  gr_eff_x_acc->SetMarkerColor(1);
  gr_eff_x_acc->Draw("AP");
  c_eff_x_acc.SaveAs(Form("~/public/UML-fit/Angular_efficiency/eff_x_acc_%i.gif",year));
  c_eff_x_acc.SaveAs(Form("~/public/UML-fit/Angular_efficiency/eff_x_acc_%i.pdf",year));

  TFile* f;
  f = new TFile(Form("~/public/UML-fit/Angular_efficiency/eff_x_acc_%i.root",year),"RECREATE");
  f->cd();
  gr_eff->Write();
  gr_wei->Write();
  gr_acc->Write();
  gr_eff_x_acc->Write();
  f->Close();

  return;
}

double DecayRate(double ctK, double ctL, double phi, double Fl, double P1, double P2, double P3, double P4p, double P5p, double P6p, double P8p){

  double dec = ( 0.75 * (1-Fl) * (1-ctK*ctK)
                 + Fl * ctK*ctK
                 + ( 0.25 * (1-Fl) * (1-ctK*ctK) - Fl * ctK*ctK ) * ( 2 * ctL*ctL -1 )
                 + 0.5 * P1 * (1-Fl) * (1-ctK*ctK) * (1-ctL*ctL) * cos(2*phi)
                 + 2 * cos(phi) * ctK * sqrt(Fl * (1-Fl) * (1-ctK*ctK)) * ( P4p * ctL * sqrt(1-ctL*ctL) + P5p * sqrt(1-ctL*ctL) )
                 + 2 * sin(phi) * ctK * sqrt(Fl * (1-Fl) * (1-ctK*ctK)) * ( P8p * ctL * sqrt(1-ctL*ctL) - P6p * sqrt(1-ctL*ctL) )
                 + 2 * P2 * (1-Fl) * (1-ctK*ctK) * ctL
                 - P3 * (1-Fl) * (1-ctK*ctK) * (1-ctL*ctL) * sin(2*phi) );

  double ret = (9./(32 * 3.14159265) * dec); // no penalty term

  return ret;
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





