   #include "../../include/WJetsReader.hh"
   #include "../../include/UserIn.cpp"

   #ifndef WJETS_READER
   #define WJETS_READER

   TH1* WJetsReader::fTwoProng_pt_signal = nullptr;
   TH1* WJetsReader::fTwoProng_pt_sideband = nullptr;
   TH1* WJetsReader::fTwoProng_pt_distr = nullptr;

   TH1* WJetsReader::fTwoProng_eta_signal = nullptr;
   TH1* WJetsReader::fTwoProng_eta_sideband = nullptr;
   TH1* WJetsReader::fTwoProng_eta_distr = nullptr;

   TH1* WJetsReader::fTwoProng_phi_signal = nullptr;
   TH1* WJetsReader::fTwoProng_phi_sideband = nullptr;
   TH1* WJetsReader::fTwoProng_phi_distr = nullptr;

   TH1* WJetsReader::fTwoProng_MassPi0_signal = nullptr;
   TH1* WJetsReader::fTwoProng_MassPi0_sideband = nullptr;
   TH1* WJetsReader::fTwoProng_MassPi0_distr = nullptr;

   TH1* WJetsReader::fW_pt_signal = nullptr;
   TH1* WJetsReader::fW_pt_signal_noW = nullptr;
   TH1* WJetsReader::fW_pt_signal_highptu = nullptr;

   TH1* WJetsReader::fW_pt_sideband = nullptr;
   TH1* WJetsReader::fW_pt_distr = nullptr;

   TH1* WJetsReader::fW_eta_signal = nullptr;
   TH1* WJetsReader::fW_eta_sideband = nullptr;
   TH1* WJetsReader::fW_eta_distr = nullptr;

   TH1* WJetsReader::fW_phi_signal = nullptr;
   TH1* WJetsReader::fW_phi_sideband = nullptr;
   TH1* WJetsReader::fW_phi_distr = nullptr;

   TH1* WJetsReader::fW_mass_signal = nullptr;
   TH1* WJetsReader::fW_mass_sideband = nullptr;
   TH1* WJetsReader::fW_mass_distr = nullptr;
   TH1* WJetsReader::fW_mass_distr_Imag_W = nullptr;

   TH1* WJetsReader::fMuon_pt_signal = nullptr;
   TH1* WJetsReader::fMuon_pt_sideband = nullptr;
   TH1* WJetsReader::fMuon_pt_distr = nullptr;

   TH1* WJetsReader::fMuon_eta_signal = nullptr;
   TH1* WJetsReader::fMuon_eta_sideband = nullptr;
   TH1* WJetsReader::fMuon_eta_distr = nullptr;

   TH1* WJetsReader::fMuon_phi_signal = nullptr;
   TH1* WJetsReader::fMuon_phi_sideband = nullptr;
   TH1* WJetsReader::fMuon_phi_distr = nullptr;

   TH1* WJetsReader::fTightMuon_pt_signal = nullptr;
   TH1* WJetsReader::fTightMuon_pt_sideband = nullptr;
   TH1* WJetsReader::fTightMuon_pt_distr = nullptr;

   TH1* WJetsReader::fTightMuon_eta_signal = nullptr;
   TH1* WJetsReader::fTightMuon_eta_sideband = nullptr;
   TH1* WJetsReader::fTightMuon_eta_distr = nullptr;

   TH1* WJetsReader::fTightMuon_phi_signal = nullptr;
   TH1* WJetsReader::fTightMuon_phi_sideband = nullptr;
   TH1* WJetsReader::fTightMuon_phi_distr = nullptr;

   TH1* WJetsReader::fHT_distr = nullptr;
   TH1* WJetsReader::fHT_signal = nullptr;
   TH1* WJetsReader::fHT_sideband = nullptr;

   TH1* WJetsReader::fHTbare_distr = nullptr;
   TH1* WJetsReader::fHTbare_signal = nullptr;
   TH1* WJetsReader::fHTbare_sideband = nullptr;

   TH1* WJetsReader::MET_pt_distr = nullptr;
   TH1* WJetsReader::MET_phi_distr = nullptr;

   TH1* WJetsReader::mT_distr = nullptr;
   TH1* WJetsReader::mT_signal = nullptr;
   
   TH1* WJetsReader::mu_w_deltaphi = nullptr;
   TH1* WJetsReader::mu_met_deltaphi = nullptr;

   TH2* WJetsReader::fTwoProng_MassPi0_PT = nullptr;
   TH2* WJetsReader::fTwoProng_MassPi0_PT_sideband = nullptr;

   TH2* WJetsReader::fTwoProng_MassPi0_PT_lowPTmuon = nullptr;
   TH2* WJetsReader::fTwoProng_MassPi0_PT_sideband_lowPTmuon = nullptr;

   TH2* WJetsReader::fTwoProng_MassPi0_PT_highPTmuon = nullptr;
   TH2* WJetsReader::fTwoProng_MassPi0_PT_sideband_highPTmuon = nullptr;

   bool WJetsReader::is_data = false;

   TH2* WJetsReader::PT_Mass_Pi0_01 = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_02 = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_03 = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_04 = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_05 = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_06 = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_07 = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_08 = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_09 = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_10 = nullptr;

   TH2* WJetsReader::PT_Mass_Pi0_01_side = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_02_side = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_03_side = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_04_side = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_05_side = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_06_side = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_07_side = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_08_side = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_09_side = nullptr;
   TH2* WJetsReader::PT_Mass_Pi0_10_side = nullptr;

   TH1* WJetsReader::sample_weights = nullptr;
   int WJetsReader::N_SAMPLES_PASSED =0;

   void WJetsReader::Init(TTree *tree)
{
   fReader.SetTree(tree); 
}

void WJetsReader::SlaveBegin(TTree *tree)
{

}

Bool_t WJetsReader::Process(Long64_t entry)
{
   fReader.SetLocalEntry(entry);

   if (*mT < MT_CUT) return kTRUE;
   
   double eventWeight(1);
   if (!is_data) eventWeight = (*mcXS) / (*mcN) * constants::LUMI2016;

   //EVENTS ARE ALREADY FILTERD FOR MUON ID, NO NEED TO CHECK IN ANALYSIS

   TLorentzVector fMuon, fTwoProng, fTwoProngLoose, fW;

   fMuon.SetPtEtaPhiM(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0]);
   fW.SetPtEtaPhiM(W_pt[0], W_eta[0], W_phi[0], W_mass[0]);
   mu_w_deltaphi->Fill(fMuon.DeltaPhi(fW), eventWeight);
   mu_met_deltaphi->Fill(Muon_phi[0] - *MET_phi, eventWeight);

   //TWO PRONG LOOPS, USED FOR CHECKING THAT THERE IS A TWO PRONG IN THE EVENT

   bool is_tight_twoprong(false);
   bool is_loose_twoprong(false);


   for (double j : TwoProng_pt){
      fTwoProng.SetPtEtaPhiM(TwoProng_pt[0], TwoProng_eta[0], TwoProng_phi[0], TwoProng_MassPi0[0]);
      is_tight_twoprong = true;
      break;
   }

   if (!is_tight_twoprong){
      for (double j : TwoProngLoose_pt){
         fTwoProngLoose.SetPtEtaPhiM(TwoProngLoose_pt[0], TwoProngLoose_eta[0], TwoProngLoose_phi[0], TwoProngLoose_MassPi0[0]);
         is_loose_twoprong = true;
         break;
      }
   }

   TLorentzVector W_vector, mu_vector;
   W_vector.SetPtEtaPhiM(W_pt[0], W_eta[0], W_phi[0], W_mass[0]);
   mu_vector.SetPtEtaPhiM(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0]);

   if (is_tight_twoprong and Muon_pt[0] > MUON_PT_CUT and mu_vector.DeltaR(fTwoProng) > GAMMA_TWOPRONG_DPHI_CUT and W_pt[0] > W_PT_CUT){

      //if (W_vector.DeltaR(fTwoProng) < 0.3) {
         fTwoProng_pt_signal->Fill(fTwoProng.Pt(), eventWeight);
         Get_Mass_PT_plot_tight(fTwoProng.Pt())->Fill(TwoProng_MassPi0[0], get_sample_num(eventWeight));
         fW_pt_signal->Fill(fW.Pt(), eventWeight);
         fW_mass_signal->Fill(fW.M(), eventWeight);
         fHT_signal->Fill(*HT, eventWeight);
         fHTbare_signal->Fill(*HT_bare, eventWeight);
         fTwoProng_eta_signal->Fill(TwoProng_eta[0], eventWeight);
         fTwoProng_phi_signal->Fill(TwoProng_phi[0], eventWeight);
         fTwoProng_MassPi0_signal->Fill(TwoProng_MassPi0[0], eventWeight);
         fMuon_pt_signal->Fill(Muon_pt[0], eventWeight);
         fMuon_phi_signal->Fill(Muon_phi[0], eventWeight);
         fMuon_eta_signal->Fill(Muon_eta[0], eventWeight);
         fTightMuon_pt_signal->Fill(TightMuon_pt[0], eventWeight);
         fTightMuon_phi_signal->Fill(TightMuon_phi[0], eventWeight);
         fTightMuon_eta_signal->Fill(TightMuon_eta[0], eventWeight);
         fW_eta_signal->Fill(W_eta[0], eventWeight);
         fW_phi_signal->Fill(W_phi[0], eventWeight);
         fTwoProng_MassPi0_PT->Fill(TwoProng_MassPi0[0], TwoProng_pt[0], eventWeight);
         mT_signal->Fill(*mT, eventWeight);
         if (TwoProng_eta[0] < 1.43) fTwoProng_MassPi0_PT_lowPTmuon->Fill(TwoProng_MassPi0[0], TwoProng_pt[0], eventWeight);
         else fTwoProng_MassPi0_PT_highPTmuon->Fill(TwoProng_MassPi0[0], TwoProng_pt[0], eventWeight);

      //}
      
   }

   if (is_loose_twoprong and Muon_pt[0] > MUON_PT_CUT and fTwoProngLoose.DeltaR(mu_vector) > GAMMA_TWOPRONG_DPHI_CUT and W_pt[0] > W_PT_CUT){
     // if (W_vector.DeltaR(fTwoProngLoose) > 0.3){
         Get_Mass_PT_plot_side(fTwoProngLoose.Pt())->Fill(TwoProngLoose_MassPi0[0], get_sample_num(eventWeight)); ///////CHANGE THIS SO THAT IT REFLECTS AN INTEGER CORRESPONDING TO EACH SAMPLE WEIGHT
         fTwoProng_pt_sideband->Fill(fTwoProngLoose.Pt(), eventWeight);
         fW_pt_sideband->Fill(fW.Pt(), eventWeight);
         fW_mass_sideband->Fill(fW.M(), eventWeight);
         fMuon_pt_sideband->Fill(fMuon.Pt(), eventWeight);
         fHT_sideband->Fill(*HT, eventWeight);
         fHTbare_sideband->Fill(*HT_bare, eventWeight);
         fTwoProng_eta_sideband->Fill(TwoProngLoose_eta[0], eventWeight);
         fTwoProng_phi_sideband->Fill(TwoProngLoose_phi[0], eventWeight);
         fTwoProng_MassPi0_sideband->Fill(TwoProngLoose_MassPi0[0], eventWeight);
         fMuon_pt_sideband->Fill(Muon_pt[0], eventWeight);
         fMuon_phi_sideband->Fill(Muon_phi[0], eventWeight);
         fMuon_eta_sideband->Fill(Muon_eta[0], eventWeight);
         fTightMuon_pt_sideband->Fill(TightMuon_pt[0], eventWeight);
         fTightMuon_phi_sideband->Fill(TightMuon_phi[0], eventWeight);
         fTightMuon_eta_sideband->Fill(TightMuon_eta[0], eventWeight);
         fW_eta_sideband->Fill(W_eta[0], eventWeight);
         fW_phi_sideband->Fill(W_phi[0], eventWeight);
         fTwoProng_MassPi0_PT_sideband->Fill(TwoProngLoose_MassPi0[0], TwoProngLoose_pt[0], eventWeight);
         if (TwoProngLoose_eta[0] < 1.43) fTwoProng_MassPi0_PT_sideband_lowPTmuon->Fill(TwoProngLoose_MassPi0[0], TwoProngLoose_pt[0], eventWeight);
         else fTwoProng_MassPi0_PT_sideband_highPTmuon->Fill(TwoProngLoose_MassPi0[0], TwoProngLoose_pt[0], eventWeight);

   //   }
   }

   fW_pt_distr->Fill(fW.Pt(), eventWeight);
   fW_mass_distr->Fill(fW.M(), eventWeight);
   fHT_distr->Fill(*HT, eventWeight);
   fHTbare_distr->Fill(*HT_bare, eventWeight);
   fMuon_pt_distr->Fill(Muon_pt[0], eventWeight);
   fMuon_phi_distr->Fill(Muon_phi[0], eventWeight);
   fMuon_eta_distr->Fill(Muon_eta[0], eventWeight);
   fTightMuon_pt_distr->Fill(TightMuon_pt[0], eventWeight);
   fTightMuon_phi_distr->Fill(TightMuon_phi[0], eventWeight);
   fTightMuon_eta_distr->Fill(TightMuon_eta[0], eventWeight);
   fW_eta_distr->Fill(W_eta[0], eventWeight);
   fW_phi_distr->Fill(W_phi[0], eventWeight);

   if (is_tight_twoprong){
      //fill 2p codes here
      fTwoProng_pt_distr->Fill(TwoProng_pt[0], eventWeight);
      fTwoProng_eta_distr->Fill(TwoProng_eta[0], eventWeight);
      fTwoProng_phi_distr->Fill(TwoProng_phi[0], eventWeight);
      fTwoProng_MassPi0_distr->Fill(TwoProng_MassPi0[0], eventWeight);

   } else if (is_loose_twoprong) {

      fTwoProng_pt_distr->Fill(TwoProngLoose_pt[0], eventWeight);
      fTwoProng_eta_distr->Fill(TwoProngLoose_eta[0], eventWeight);
      fTwoProng_phi_distr->Fill(TwoProngLoose_phi[0], eventWeight);
      fTwoProng_MassPi0_distr->Fill(TwoProngLoose_MassPi0[0], eventWeight);

   }

   MET_pt_distr->Fill(*MET_pt, eventWeight);
   MET_phi_distr->Fill(*MET_phi, eventWeight);

   mT_distr->Fill(*mT, eventWeight);


   return kTRUE;
}

void WJetsReader::Terminate()
{
   

}

TH2* WJetsReader::Get_Mass_PT_plot_tight(double pt_val){
   for (int j = 0; j < PT_BOUNDS1.size() - 1; j++){
      if (pt_val > PT_BOUNDS1[j] and pt_val < PT_BOUNDS1[j + 1]) return get_plot(j);
   }
   if (pt_val > PT_BOUNDS1[PT_BOUNDS1.size()-1]) return PT_Mass_Pi0_10;
   else return nullptr;
}

TH2* WJetsReader::Get_Mass_PT_plot_side(double pt_val){
   for (int j = 0; j < PT_BOUNDS1.size() - 1; j++){
      if (pt_val > PT_BOUNDS1[j] and pt_val < PT_BOUNDS1[j + 1]) return get_plot_side(j);
   }
   if (pt_val > PT_BOUNDS1[PT_BOUNDS1.size()-1]) return PT_Mass_Pi0_10_side;
   else return nullptr;
}

TH2* WJetsReader::get_plot(int num){

   if (num == 0) return WJetsReader::PT_Mass_Pi0_01;
   if (num == 1) return WJetsReader::PT_Mass_Pi0_02;
   if (num == 2) return WJetsReader::PT_Mass_Pi0_03;
   if (num == 3) return WJetsReader::PT_Mass_Pi0_04;
   if (num == 4) return WJetsReader::PT_Mass_Pi0_05;
   if (num == 5) return WJetsReader::PT_Mass_Pi0_06;
   if (num == 6) return WJetsReader::PT_Mass_Pi0_07;
   if (num == 7) return WJetsReader::PT_Mass_Pi0_08;
   if (num == 8) return WJetsReader::PT_Mass_Pi0_09;
   if (num == 9) return WJetsReader::PT_Mass_Pi0_10;

   return nullptr;
}

TH2* WJetsReader::get_plot_side(int num){

   if (num == 0) return WJetsReader::PT_Mass_Pi0_01_side;
   if (num == 1) return WJetsReader::PT_Mass_Pi0_02_side;
   if (num == 2) return WJetsReader::PT_Mass_Pi0_03_side;
   if (num == 3) return WJetsReader::PT_Mass_Pi0_04_side;
   if (num == 4) return WJetsReader::PT_Mass_Pi0_05_side;
   if (num == 5) return WJetsReader::PT_Mass_Pi0_06_side;
   if (num == 6) return WJetsReader::PT_Mass_Pi0_07_side;
   if (num == 7) return WJetsReader::PT_Mass_Pi0_08_side;
   if (num == 8) return WJetsReader::PT_Mass_Pi0_09_side;
   if (num == 9) return WJetsReader::PT_Mass_Pi0_10_side;

   return nullptr;
}

int WJetsReader::get_sample_num(double weight){
   for (int j = 1; j <= WJetsReader::N_SAMPLES_PASSED; j++){
      if ( WJetsReader::sample_weights->GetBinContent(j) == weight ) return j;
   }
   WJetsReader::N_SAMPLES_PASSED++;
   WJetsReader::sample_weights->SetBinContent(N_SAMPLES_PASSED, weight);
   return WJetsReader::N_SAMPLES_PASSED;
}


#endif