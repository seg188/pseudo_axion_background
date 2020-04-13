#include <iostream>

#ifndef GJetsReader_C
#define GJetsReader_C


#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TSelector.h"

const Int_t kMaxfParticles = 1293;

class GJetsReader : public TSelector {

public :

   TTreeReader fReader;
   TTreeReaderArray<double> TwoProngLoose_pt;
   TTreeReaderArray<double> TwoProngLoose_eta;
   TTreeReaderArray<double> TwoProngLoose_phi;
   TTreeReaderArray<double> TwoProngLoose_MassPi0;
   TTreeReaderArray<double> TwoProngLoose_MassEta;
   TTreeReaderArray<double> TwoProng_MassEta;
   TTreeReaderArray<double> TwoProngLoose_mass;
   TTreeReaderArray<double> TwoProng_mass;
   TTreeReaderArray<double> TwoProng_MassPi0;
   TTreeReaderArray<double> TwoProng_pt;
   TTreeReaderArray<double> TwoProng_eta;
   TTreeReaderArray<double> TwoProng_phi;
   TTreeReaderArray<double> IDPhoton_mass;
   TTreeReaderArray<double> IDPhoton_pt;
   TTreeReaderArray<double> IDPhoton_eta;
   TTreeReaderArray<double> IDPhoton_phi;
   TTreeReaderArray<double> Loose1IDPhoton_mass;
   TTreeReaderArray<double> Loose1IDPhoton_pt;
   TTreeReaderArray<double> Loose1IDPhoton_eta;
   TTreeReaderArray<double> Loose1IDPhoton_phi;
   TTreeReaderValue<double> HT;
   TTreeReaderArray<double> MET;
   TTreeReaderValue<double> HT_gen;
   TTreeReaderValue<double> HT_bare;
   TTreeReaderValue<double> pthat;

   TTreeReaderArray<double> Jet_mass;
   TTreeReaderArray<double> Jet_pt;
   TTreeReaderArray<double> Jet_eta;
   TTreeReaderArray<double> Jet_phi;

   TTreeReaderArray<double> TwoProng_CHpos_mass;
   TTreeReaderArray<double> TwoProng_CHpos_pt;
   TTreeReaderArray<double> TwoProng_CHpos_eta;
   TTreeReaderArray<double> TwoProng_CHpos_phi;

   TTreeReaderArray<double> TwoProngLoose_CHpos_mass;
   TTreeReaderArray<double> TwoProngLoose_CHpos_pt;
   TTreeReaderArray<double> TwoProngLoose_CHpos_eta;
   TTreeReaderArray<double> TwoProngLoose_CHpos_phi;

   TTreeReaderArray<double> TwoProng_CHneg_mass;
   TTreeReaderArray<double> TwoProng_CHneg_pt;
   TTreeReaderArray<double> TwoProng_CHneg_eta;
   TTreeReaderArray<double> TwoProng_CHneg_phi;
   TTreeReaderArray<double> TwoProng_trackAsym;
   TTreeReaderArray<double> TwoProngLoose_trackAsym;
   TTreeReaderArray<double> TwoProng_photon_pt;
   TTreeReaderArray<double> TwoProngLoose_photon_pt;
   TTreeReaderArray<double> TwoProng_photon_eta;
   TTreeReaderArray<double> TwoProngLoose_photon_eta;
   TTreeReaderArray<double> TwoProng_photonAsym;
   TTreeReaderArray<double> TwoProngLoose_photonAsym;

   TTreeReaderArray<double> TwoProngLoose_CHneg_mass;
   TTreeReaderArray<double> TwoProngLoose_CHneg_pt;
   TTreeReaderArray<double> TwoProngLoose_CHneg_eta;
   TTreeReaderArray<double> TwoProngLoose_CHneg_phi;
   
   TTreeReaderArray<double> Photon_mass;
   TTreeReaderArray<double> Photon_pt;
   TTreeReaderArray<double> Photon_eta;
   TTreeReaderArray<double> Photon_phi;
   TTreeReaderArray<double> Photon_HE;
   TTreeReaderArray<double> Photon_isoCh;
   TTreeReaderArray<double> Photon_isoGamma;
   TTreeReaderArray<double> Photon_sigmaIetaIeta;
   TTreeReaderArray<int> Photon_passVeto;
   
   TTreeReaderArray<double> TwoProngLoose_chargedIso;
   TTreeReaderArray<double> TwoProngLoose_neutralIso;
   TTreeReaderArray<double> TwoProngLoose_egammaIso;
   TTreeReaderArray<double> TwoProng_chargedIso;
   TTreeReaderArray<double> TwoProng_neutralIso;
   TTreeReaderArray<double> TwoProng_egammaIso;
   TTreeReaderValue<double> mcN;
   TTreeReaderValue<int>    nPV;
   TTreeReaderValue<int>    nJets;
   TTreeReaderValue<double> mcXS;
   TTreeReaderValue<int>    eventNum;
   TTreeReaderValue<int>    nTwoProngs;
   TTreeReaderValue<int>    nTwoProngCands;





GJetsReader(TTree * = 0) : 
      TwoProngLoose_MassEta(fReader, "TwoProngLoose_MassEta"),
      TwoProng_MassEta(fReader, "TwoProng_MassEta"),
      TwoProngLoose_MassPi0(fReader, "TwoProngLoose_MassPi0"),
      TwoProngLoose_pt(fReader, "TwoProngLoose_pt"), 
      TwoProngLoose_eta(fReader, "TwoProngLoose_eta"), 
      TwoProngLoose_phi(fReader, "TwoProngLoose_phi"),
      TwoProngLoose_mass(fReader, "TwoProngLoose_mass"),
      TwoProng_mass(fReader, "TwoProng_mass"),  
      TwoProng_pt(fReader, "TwoProng_pt"), 
      TwoProng_eta(fReader, "TwoProng_eta"), 
      TwoProng_phi(fReader, "TwoProng_phi"),
      IDPhoton_mass(fReader, "IDPhoton_mass"),
      TwoProng_MassPi0(fReader, "TwoProng_MassPi0"),
      IDPhoton_pt(fReader, "IDPhoton_pt"),
      IDPhoton_eta(fReader, "IDPhoton_eta"),
      IDPhoton_phi(fReader, "IDPhoton_phi"),
      HT(fReader, "HT"),
      HT_bare(fReader, "HT_bare"),
      MET(fReader, "MET"),
      HT_gen(fReader, "HT_gen"),
      Jet_mass(fReader, "Jet_mass"),  
      Jet_pt(fReader, "Jet_pt"), 
      Jet_eta(fReader, "Jet_eta"), 
      Jet_phi(fReader, "Jet_phi"),

      TwoProng_CHpos_mass(fReader, "TwoProng_CHpos_mass"),
      TwoProng_CHpos_pt(fReader, "TwoProng_CHpos_pt"),
      TwoProng_CHpos_eta(fReader, "TwoProng_CHpos_eta"),
      TwoProng_CHpos_phi(fReader, "TwoProng_CHpos_phi"),
      TwoProng_trackAsym(fReader, "TwoProng_trackAsym"),

      TwoProngLoose_CHpos_mass(fReader, "TwoProngLoose_CHpos_mass"),
      TwoProngLoose_CHpos_pt(fReader, "TwoProngLoose_CHpos_pt"),
      TwoProngLoose_CHpos_eta(fReader, "TwoProngLoose_CHpos_eta"),
      TwoProngLoose_CHpos_phi(fReader, "TwoProngLoose_CHpos_phi"),
      TwoProngLoose_trackAsym(fReader, "TwoProngLoose_trackAsym"),

      TwoProng_CHneg_mass(fReader, "TwoProng_CHneg_mass"),
      TwoProng_CHneg_pt(fReader, "TwoProng_CHneg_pt"),
      TwoProng_CHneg_eta(fReader, "TwoProng_CHneg_eta"),
      TwoProng_CHneg_phi(fReader, "TwoProng_CHneg_phi"),
      TwoProng_photon_pt(fReader, "TwoProng_photon_pt"),
      TwoProngLoose_photon_pt(fReader, "TwoProngLoose_photon_pt"),
      TwoProng_photon_eta(fReader, "TwoProng_photon_eta"),
      TwoProngLoose_photon_eta(fReader, "TwoProngLoose_photon_eta"),
      TwoProng_photonAsym(fReader, "TwoProng_photonAsym"),
      TwoProngLoose_photonAsym(fReader, "TwoProngLoose_photonAsym"),

      TwoProngLoose_CHneg_mass(fReader, "TwoProngLoose_CHneg_mass"),
      TwoProngLoose_CHneg_pt(fReader, "TwoProngLoose_CHneg_pt"),
      TwoProngLoose_CHneg_eta(fReader, "TwoProngLoose_CHneg_eta"),
      TwoProngLoose_CHneg_phi(fReader, "TwoProngLoose_CHneg_phi"),
       
      Loose1IDPhoton_mass(fReader, "Loose1IDPhoton_mass"),
      Loose1IDPhoton_pt(fReader, "Loose1IDPhoton_pt"),
      Loose1IDPhoton_eta(fReader, "Loose1IDPhoton_eta"),
      Loose1IDPhoton_phi(fReader, "Loose1IDPhoton_phi"),
     
      pthat(fReader, "pthat"),
      Photon_mass(fReader, "Photon_mass"),
      Photon_pt(fReader, "Photon_pt"),
      Photon_eta(fReader, "Photon_eta"),
      Photon_phi(fReader, "Photon_phi"),
      Photon_HE(fReader, "Photon_HE"),
      Photon_isoCh(fReader, "Photon_isoCh"),
      Photon_isoGamma(fReader, "Photon_isoGamma"),
      Photon_sigmaIetaIeta(fReader, "Photon_sigmaIetaIeta"),
      Photon_passVeto(fReader, "Photon_passVeto"),

      
      nTwoProngs(fReader, "nTwoProngs"),
      TwoProngLoose_chargedIso(fReader, "TwoProngLoose_chargedIso"),
      TwoProngLoose_neutralIso(fReader, "TwoProngLoose_neutralIso"),
      TwoProngLoose_egammaIso(fReader, "TwoProngLoose_egammaIso"),
      TwoProng_chargedIso(fReader, "TwoProng_chargedIso"),
      TwoProng_neutralIso(fReader, "TwoProng_neutralIso"),
      TwoProng_egammaIso(fReader, "TwoProng_egammaIso"),
      mcN(fReader, "mcN"),
      nPV(fReader, "nPV"),
      nJets(fReader, "nJets"),
      mcXS(fReader, "mcXS"),
      nTwoProngCands(fReader, "nTwoProngCands"),
      eventNum(fReader, "eventNum")
   
   {
   }
   //USER HISTOGRAMS

   static TH1* fTwoProng_MassPi0_pt40_60_test;
   static TH1* fTwoProng_MassPi0_pt40_60;
   static TH1* fTwoProng_pt_signal;
   static TH1* fTwoProng_pt_sideband;
   static TH1* fTwoProng_pt_distr;

   static TH1* fTwoProng_eta_signal;
   static TH1* fTwoProng_eta_sideband;
   static TH1* fTwoProng_eta_distr;

   static TH1* fTwoProng_phi_signal;
   static TH1* fTwoProng_phi_sideband;
   static TH1* fTwoProng_phi_distr;

   static TH1* fTwoProng_MassPi0_signal;
   static TH1* fTwoProng_MassPi0_sideband;
   static TH1* fTwoProng_MassPi0_distr;

   static TH1* fPhi_pt_signal;
   static TH1* fPhi_pt_sideband;
   static TH1* fPhi_pt_distr;

   static TH1* fPhi_eta_signal;
   static TH1* fPhi_eta_sideband;
   static TH1* fPhi_eta_distr;

   static TH1* fPhi_phi_signal;
   static TH1* fPhi_phi_sideband;
   static TH1* fPhi_phi_distr;

   static TH1* fPhi_mass_signal;
   static TH1* fPhi_mass_sideband;
   static TH1* fPhi_mass_distr;

   static TH1* fIDPhoton_pt_signal;
   static TH1* fIDPhoton_pt_sideband;
   static TH1* fIDPhoton_pt_distr;

   static TH1* fIDPhoton_eta_signal;
   static TH1* fIDPhoton_eta_sideband;
   static TH1* fIDPhoton_eta_distr;

   static TH1* fIDPhoton_phi_signal;
   static TH1* fIDPhoton_phi_sideband;
   static TH1* fIDPhoton_phi_distr;

   static TH1* fTightIDPhoton_pt_signal;
   static TH1* fTightIDPhoton_pt_sideband;
   static TH1* fTightIDPhoton_pt_distr;

   static TH1* fTightIDPhoton_eta_signal;
   static TH1* fTightIDPhoton_eta_sideband;
   static TH1* fIDPhoton_eta_distr;

   static TH1* fIDPhoton_phi_signal;
   static TH1* fIDPhoton_phi_sideband;
   static TH1* fIDPhoton_phi_distr;

   static TH1* fTightIDPhoton_pt_signal;
   static TH1* fTightIDPhoton_pt_sideband;
   static TH1* fTightIDPhoton_pt_distr;

   static TH1* fTightIDPhoton_eta_signal;
   static TH1* fTightIDPhoton_eta_sideband;
   static TH1* fTightIDPhoton_eta_distr;

   static TH1* fTightIDPhoton_phi_signal;
   static TH1* fTightIDPhoton_phi_sideband;
   static TH1* fTightIDPhoton_phi_distr;

   static TH1* fHT_distr;
   static TH1* fHT_signal;
   static TH1* fHT_sideband;

   static TH1* fHTbare_distr;
   static TH1* fHTbare_signal;
   static TH1* fHTbare_sideband;

   static TH1* MET_pt_distr;
   static TH1* MET_phi_distr;

  // static TH1* mT_distr;
  // static TH1* mT_signal;

  // static TH1* mu_w_deltaphi;
  // static TH1* mu_met_deltaphi;

   static TH2* fTwoProng_MassPi0_PT;
   static TH2* fTwoProng_MassPi0_PT_sideband;

   static TH2* fTwoProng_MassPi0_PT_lowPTIDPhoton;
   static TH2* fTwoProng_MassPi0_PT_sideband_lowPTIDPhoton;

   static TH2* fTwoProng_MassPi0_PT_highPTIDPhoton;
   static TH2* fTwoProng_MassPi0_PT_sideband_highPTIDPhoton;

   static bool is_data;

   static TH2* PT_Mass_Pi0_01;
   static TH2* PT_Mass_Pi0_02;
   static TH2* PT_Mass_Pi0_03;
   static TH2* PT_Mass_Pi0_04;
   static TH2* PT_Mass_Pi0_05;
   static TH2* PT_Mass_Pi0_06;
   static TH2* PT_Mass_Pi0_07;
   static TH2* PT_Mass_Pi0_08;
   static TH2* PT_Mass_Pi0_09;
   static TH2* PT_Mass_Pi0_10;

   static TH2* PT_Mass_Pi0_01_side;
   static TH2* PT_Mass_Pi0_02_side;
   static TH2* PT_Mass_Pi0_03_side;
   static TH2* PT_Mass_Pi0_04_side;
   static TH2* PT_Mass_Pi0_05_side;
   static TH2* PT_Mass_Pi0_06_side;
   static TH2* PT_Mass_Pi0_07_side;
   static TH2* PT_Mass_Pi0_08_side;
   static TH2* PT_Mass_Pi0_09_side;
   static TH2* PT_Mass_Pi0_10_side;

   static TH2* Get_Mass_PT_plot_side(double pt_val);
   static TH2* Get_Mass_PT_plot_tight(double pt_val);
   static TH2* get_plot(int num);
   static TH2* get_plot_side(int num);
   static int get_sample_num(double weight);
   static TH1* sample_weights;
   static int N_SAMPLES_PASSED;


   virtual ~GJetsReader() { }

   virtual void    Init(TTree *tree);

   virtual void    SlaveBegin(TTree *tree);

   virtual Bool_t  Process(Long64_t entry);

   virtual void    Terminate();

   virtual Int_t   Version() const { return 2; }

   ClassDef(GJetsReader,0);

};




#endif