#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TSelector.h"
#include "UserIn.cpp"

using namespace user_input;

#ifndef WJETSREADER_C
#define WJETSREADER_C

const Int_t kMaxfParticles = 1293;

class WJetsReader : public TSelector {
public :

   virtual ~WJetsReader() { }

   virtual void    Init(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual Bool_t  Process(Long64_t entry);
   virtual void    Terminate();
   virtual Int_t   Version() const { return 2; }

   //TREE VARIABLES
   ////////////////////////////////////////////////////////////////////////////////
   TTreeReader fReader;
   TTreeReaderValue<double> HT;
   TTreeReaderArray<int>    Imag_W;
   TTreeReaderValue<double> HT_bare;
   TTreeReaderValue<double> mT;
   TTreeReaderValue<double> MET_pt;
   TTreeReaderValue<double> MET_phi;
   //TTreeReaderArray<double> MET;
   TTreeReaderArray<double> Jet_mass;
   TTreeReaderArray<double> Jet_pt;
   TTreeReaderArray<double> Jet_eta;
   TTreeReaderArray<double> Jet_phi;
   TTreeReaderArray<double> Muon_mass;
   TTreeReaderArray<double> Muon_pt;
   TTreeReaderArray<double> Muon_eta;
   TTreeReaderArray<double> Muon_phi;
   TTreeReaderArray<double> TightMuon_mass;
   TTreeReaderArray<double> TightMuon_pt;
   TTreeReaderArray<double> TightMuon_eta;
   TTreeReaderArray<double> TightMuon_phi;
   TTreeReaderArray<double> W_mass;
   TTreeReaderArray<double> W_pt;
   TTreeReaderArray<double> W_eta;
   TTreeReaderArray<double> W_phi;
   TTreeReaderArray<double> TwoProng_MassEta;
   TTreeReaderArray<double> TwoProng_MassPi0;
   TTreeReaderArray<double> TwoProng_pt;
   TTreeReaderArray<double> TwoProng_eta;
   TTreeReaderArray<double> TwoProng_phi;
   TTreeReaderArray<double> TwoProngLoose_MassEta;
   TTreeReaderArray<double> TwoProngLoose_MassPi0;
   TTreeReaderArray<double> TwoProngLoose_pt;
   TTreeReaderArray<double> TwoProngLoose_eta;
   TTreeReaderArray<double> TwoProngLoose_phi;
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
   TTreeReaderArray<double> TwoProngLoose_chargedIso;
   TTreeReaderArray<double> TwoProngLoose_neutralIso;
   TTreeReaderArray<double> TwoProngLoose_egammaIso;
   TTreeReaderArray<double> TwoProng_chargedIso;
   TTreeReaderArray<double> TwoProng_neutralIso;
   TTreeReaderArray<double> TwoProng_egammaIso;
   TTreeReaderValue<double> mcN;
   TTreeReaderValue<int>    nPV;
   //TTreeReaderValue<int>    nJets;
   TTreeReaderValue<double> mcXS;
   TTreeReaderValue<int>    nTwoProngs;

   WJetsReader(TTree * = 0) :

      HT(fReader, "HT"),
      Imag_W(fReader, "Imag_W"),
      HT_bare(fReader, "HT"), //HT BARE NOT PROPERLY INITLIZED IN TREE, CHANGE THIS LATER
      MET_phi(fReader, "MET_phi"),
      MET_pt(fReader, "MET_pt"),
      mT(fReader, "mT"),

      Jet_mass(fReader, "Jet_mass"),
      Jet_pt(fReader, "Jet_pt"),
      Jet_eta(fReader, "Jet_eta"),
      Jet_phi(fReader, "Jet_phi"),

      Muon_mass(fReader, "Muon_mass"),  
      Muon_pt(fReader, "Muon_pt"), 
      Muon_eta(fReader, "Muon_eta"), 
      Muon_phi(fReader, "Muon_phi"),

      TightMuon_mass(fReader, "TightMuon_mass"),  
      TightMuon_pt(fReader, "TightMuon_pt"), 
      TightMuon_eta(fReader, "TightMuon_eta"), 
      TightMuon_phi(fReader, "TightMuon_phi"),

      W_mass(fReader, "W_mass"),  
      W_pt(fReader, "W_pt"), 
      W_eta(fReader, "W_eta"), 
      W_phi(fReader, "W_phi"),

      TwoProng_MassEta(fReader, "TwoProng_MassEta"),  
      TwoProng_pt(fReader, "TwoProng_pt"), 
      TwoProng_eta(fReader, "TwoProng_eta"), 
      TwoProng_phi(fReader, "TwoProng_phi"),
      TwoProng_MassPi0(fReader, "TwoProng_MassPi0"),

      TwoProngLoose_MassEta(fReader, "TwoProngLoose_MassEta"),  
      TwoProngLoose_pt(fReader, "TwoProngLoose_pt"), 
      TwoProngLoose_eta(fReader, "TwoProngLoose_eta"), 
      TwoProngLoose_phi(fReader, "TwoProngLoose_phi"),
      TwoProngLoose_MassPi0(fReader, "TwoProngLoose_MassPi0"),

      TwoProng_CHpos_mass(fReader, "TwoProng_CHpos_mass"),
      TwoProng_CHpos_pt(fReader, "TwoProng_CHpos_pt"),
      TwoProng_CHpos_eta(fReader, "TwoProng_CHpos_eta"),
      TwoProng_CHpos_phi(fReader, "TwoProng_CHpos_phi"),

      TwoProngLoose_CHpos_mass(fReader, "TwoProngLoose_CHpos_mass"),
      TwoProngLoose_CHpos_pt(fReader, "TwoProngLoose_CHpos_pt"),
      TwoProngLoose_CHpos_eta(fReader, "TwoProngLoose_CHpos_eta"),
      TwoProngLoose_CHpos_phi(fReader, "TwoProngLoose_CHpos_phi"),

      TwoProng_CHneg_mass(fReader, "TwoProng_CHneg_mass"),
      TwoProng_CHneg_pt(fReader, "TwoProng_CHneg_pt"),
      TwoProng_CHneg_eta(fReader, "TwoProng_CHneg_eta"),
      TwoProng_CHneg_phi(fReader, "TwoProng_CHneg_phi"),
      TwoProng_trackAsym(fReader, "TwoProng_trackAsym"),
      TwoProngLoose_trackAsym(fReader, "TwoProngLoose_trackAsym"),
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
   
      TwoProngLoose_chargedIso(fReader, "TwoProngLoose_chargedIso"),
      TwoProngLoose_neutralIso(fReader, "TwoProngLoose_neutralIso"),
      TwoProngLoose_egammaIso(fReader, "TwoProngLoose_egammaIso"),

      TwoProng_chargedIso(fReader, "TwoProng_chargedIso"),
      TwoProng_neutralIso(fReader, "TwoProng_neutralIso"),
      TwoProng_egammaIso(fReader, "TwoProng_egammaIso"),

      mcN(fReader, "mcN"),
      nPV(fReader, "nPV"),
      //nJets(fReader, "nJets"),
      mcXS(fReader, "mcXS"),

      nTwoProngs(fReader, "nTwoProngs")

   { }

   ///////////////////////////////////////////////////////////////////////////////////////
   //USER HISTOGRAMS TO BE WRITTEN TO FILE, RUN FILE

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

   static TH1* fW_pt_signal;
   static TH1* fW_pt_signal_noW;
   static TH1* fW_pt_signal_highptu;

   static TH1* fW_pt_sideband;
   static TH1* fW_pt_distr;

   static TH1* fW_eta_signal;
   static TH1* fW_eta_sideband;
   static TH1* fW_eta_distr;

   static TH1* fW_phi_signal;
   static TH1* fW_phi_sideband;
   static TH1* fW_phi_distr;

   static TH1* fW_mass_signal;
   static TH1* fW_mass_sideband;
   static TH1* fW_mass_distr;
   static TH1* fW_mass_distr_Imag_W;

   static TH1* fMuon_pt_signal;
   static TH1* fMuon_pt_sideband;
   static TH1* fMuon_pt_distr;

   static TH1* fMuon_eta_signal;
   static TH1* fMuon_eta_sideband;
   static TH1* fMuon_eta_distr;

   static TH1* fMuon_phi_signal;
   static TH1* fMuon_phi_sideband;
   static TH1* fMuon_phi_distr;

   static TH1* fTightMuon_pt_signal;
   static TH1* fTightMuon_pt_sideband;
   static TH1* fTightMuon_pt_distr;

   static TH1* fTightMuon_eta_signal;
   static TH1* fTightMuon_eta_sideband;
   static TH1* fTightMuon_eta_distr;

   static TH1* fTightMuon_phi_signal;
   static TH1* fTightMuon_phi_sideband;
   static TH1* fTightMuon_phi_distr;

   static TH1* fHT_distr;
   static TH1* fHT_signal;
   static TH1* fHT_sideband;

   static TH1* fHTbare_distr;
   static TH1* fHTbare_signal;
   static TH1* fHTbare_sideband;

   static TH1* MET_pt_distr;
   static TH1* MET_phi_distr;

   static TH1* mT_distr;
   static TH1* mT_signal;
   
   static TH1* mu_w_deltaphi;
   static TH1* mu_met_deltaphi;

   static TH2* fTwoProng_MassPi0_PT;
   static TH2* fTwoProng_MassPi0_PT_sideband;

   static TH2* fTwoProng_MassPi0_PT_lowPTmuon;
   static TH2* fTwoProng_MassPi0_PT_sideband_lowPTmuon;

   static TH2* fTwoProng_MassPi0_PT_highPTmuon;
   static TH2* fTwoProng_MassPi0_PT_sideband_highPTmuon;

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

   ClassDef(WJetsReader,0);
};





#endif