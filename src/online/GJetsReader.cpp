#include "../../include/GJetsReader.hh"
#include "../../include/UserIn.cpp"

#ifndef GJETS_READER
#define GJETS_READER

   TH1* GJetsReader::fTwoProng_pt_signal = nullptr;
   TH1* GJetsReader::fTwoProng_pt_sideband = nullptr;
   TH1* GJetsReader::fTwoProng_pt_distr = nullptr;

   TH1* GJetsReader::fTwoProng_eta_signal = nullptr;
   TH1* GJetsReader::fTwoProng_eta_sideband = nullptr;
   TH1* GJetsReader::fTwoProng_eta_distr = nullptr;

   TH1* GJetsReader::fTwoProng_phi_signal = nullptr;
   TH1* GJetsReader::fTwoProng_phi_sideband = nullptr;
   TH1* GJetsReader::fTwoProng_phi_distr = nullptr;

   TH1* GJetsReader::fTwoProng_MassPi0_signal = nullptr;
   TH1* GJetsReader::fTwoProng_MassPi0_sideband = nullptr;
   TH1* GJetsReader::fTwoProng_MassPi0_distr = nullptr;

   TH1* GJetsReader::fPhi_pt_signal = nullptr;
   TH1* GJetsReader::fPhi_pt_sideband = nullptr;
   TH1* GJetsReader::fPhi_pt_distr = nullptr;

   TH1* GJetsReader::fPhi_eta_signal = nullptr;
   TH1* GJetsReader::fPhi_eta_sideband = nullptr;
   TH1* GJetsReader::fPhi_eta_distr = nullptr;

   TH1* GJetsReader::fPhi_phi_signal = nullptr;
   TH1* GJetsReader::fPhi_phi_sideband = nullptr;
   TH1* GJetsReader::fPhi_phi_distr = nullptr;

   TH1* GJetsReader::fPhi_mass_signal = nullptr;
   TH1* GJetsReader::fPhi_mass_sideband = nullptr;
   TH1* GJetsReader::fPhi_mass_distr = nullptr;


   TH1* GJetsReader::fIDPhoton_pt_signal = nullptr;
   TH1* GJetsReader::fIDPhoton_pt_sideband = nullptr;
   TH1* GJetsReader::fIDPhoton_pt_distr = nullptr;

   TH1* GJetsReader::fIDPhoton_eta_signal = nullptr;
   TH1* GJetsReader::fIDPhoton_eta_sideband = nullptr;
   TH1* GJetsReader::fIDPhoton_eta_distr = nullptr;

   TH1* GJetsReader::fIDPhoton_phi_signal = nullptr;
   TH1* GJetsReader::fIDPhoton_phi_sideband = nullptr;
   TH1* GJetsReader::fIDPhoton_phi_distr = nullptr;

   TH1* GJetsReader::fHT_distr = nullptr;
   TH1* GJetsReader::fHT_signal = nullptr;
   TH1* GJetsReader::fHT_sideband = nullptr;

   TH1* GJetsReader::fHTbare_distr = nullptr;
   TH1* GJetsReader::fHTbare_signal = nullptr;
   TH1* GJetsReader::fHTbare_sideband = nullptr;

   TH1* GJetsReader::MET_pt_distr = nullptr;
   TH1* GJetsReader::MET_phi_distr = nullptr;

   TH2* GJetsReader::fTwoProng_MassPi0_PT = nullptr;
   TH2* GJetsReader::fTwoProng_MassPi0_PT_sideband = nullptr;

   bool GJetsReader::is_data = false;

   TH2* GJetsReader::PT_Mass_Pi0_01 = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_02 = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_03 = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_04 = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_05 = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_06 = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_07 = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_08 = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_09 = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_10 = nullptr;

   TH2* GJetsReader::PT_Mass_Pi0_01_side = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_02_side = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_03_side = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_04_side = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_05_side = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_06_side = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_07_side = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_08_side = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_09_side = nullptr;
   TH2* GJetsReader::PT_Mass_Pi0_10_side = nullptr;

   TH2* GJetsReader::fPhiMass_PT_signal = nullptr;
   TH2* GJetsReader::fPhiMass_PT_sideband = nullptr;

   TH1* GJetsReader::sample_weights = nullptr;
   int GJetsReader::N_SAMPLES_PASSED =0;

bool getID(double fgHE, double fgG_Iso, double fgC_Iso, double fgSIE)
{
  bool process = false;

  if (fgHE < 0.05)
  {
    if (fgC_Iso < 5)
    {
      if (fgG_Iso < 2.75)
      {
        if (fgSIE < 0.0105) 
        {
          process = true;
          //
        }
      }
    }

  }

  return process;
}


void GJetsReader::Init(TTree *tree)
{
   fReader.SetTree(tree);
}

void GJetsReader::SlaveBegin(TTree *tree)
{

}

Bool_t GJetsReader::Process(Long64_t entry)
{
   using namespace user_input;
   fReader.SetLocalEntry(entry);

   //if (*mT < 50.0) return kTRUE;
  // std::cout << "In Process Method"  << std::endl; 
   double eventWeight(1);
   if (!is_data) eventWeight = *mcXS / *mcN * constants::LUMI2016;

   //EVENTS ARE ALREADY FILTERD FOR IDPhoton ID, NO NEED TO CHECK IN ANALYSIS
   bool ID_photon = false;
   TLorentzVector fIDPhoton, fTwoProng, fTwoProngLoose;



   double gamma_pt = 0;
   for (double k : Photon_pt){
        fIDPhoton.SetPtEtaPhiM(Photon_pt[0], Photon_eta[0], Photon_phi[0], Photon_mass[0]);
        gamma_pt  = Photon_pt[0];
        ID_photon = true;
        break;
   }

   if (!ID_photon) return kTRUE;


   //TWO PRONG LOOPS, USED FOR CHECKING THAT THERE IS A TWO PRONG IN THE EVENT

   bool is_tight_twoprong(false);
   bool is_loose_twoprong(false);
//**************************************************************************************8
//RIGHT HERE WE NEED TO FIX THE DEFINITION OF THE LOOSE VS TIGHT TWO PRONGS....
//THE CURRENT SELECTION IS INCLUSIVE OF N_ISO > 0.1 FOR LOOSE... THESE NEED TO GO

   for (double j : TwoProng_pt){
      fTwoProng.SetPtEtaPhiM(TwoProng_pt[0], TwoProng_eta[0], TwoProng_phi[0], TwoProng_MassPi0[0]);
//      std::cout << TwoProng_pt[0] << std::endl;
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

   if (!is_tight_twoprong and !is_loose_twoprong) return kTRUE;


   TLorentzVector Phi_vector, gamma_vector;

   gamma_vector = fIDPhoton;

   if (is_tight_twoprong) {
      Phi_vector = gamma_vector + fTwoProng;
    } else {
      Phi_vector = gamma_vector + fTwoProngLoose;
    }


  if (is_tight_twoprong){
   if (gamma_pt > PHOTON_PT_CUT and gamma_vector.DeltaR(fTwoProng) > GAMMA_TWOPRONG_DPHI_CUT){

      //if (W_vector.DeltaR(fTwoProng) < 0.3) {
         fTwoProng_pt_signal->Fill(fTwoProng.Pt(), eventWeight);
         Get_Mass_PT_plot_tight(fTwoProng.Pt())->Fill(TwoProng_MassPi0[0], get_sample_num(eventWeight));
        //W_pt_signal->Fill(fW.Pt(), eventWeight);
        // fW_mass_signal->Fill(fW.M(), eventWeight);
         fHT_signal->Fill(*HT, eventWeight);
         fHTbare_signal->Fill(*HT_bare, eventWeight);
         fTwoProng_eta_signal->Fill(TwoProng_eta[0], eventWeight);
         fTwoProng_phi_signal->Fill(TwoProng_phi[0], eventWeight);
         fTwoProng_MassPi0_signal->Fill(TwoProng_MassPi0[0], eventWeight);
         fIDPhoton_pt_signal->Fill(Photon_pt[0], eventWeight);
         fIDPhoton_phi_signal->Fill(Photon_phi[0], eventWeight);
         fIDPhoton_eta_signal->Fill(Photon_eta[0], eventWeight);

         double Phi_Mass = (gamma_vector + fTwoProng).M();
   
        // fW_eta_signal->Fill(W_eta[0], eventWeight);
        // fW_phi_signal->Fill(W_phi[0], eventWeight);
         fTwoProng_MassPi0_PT->Fill(TwoProng_MassPi0[0], TwoProng_pt[0], eventWeight);

         fPhiMass_PT_signal->Fill(Phi_Mass,TwoProng_pt[0], eventWeight );
        // mT_signal->Fill(*mT, eventWeight);

   }
 }

if (is_loose_twoprong){
  if (gamma_pt > PHOTON_PT_CUT and fTwoProngLoose.DeltaR(gamma_vector) > GAMMA_TWOPRONG_DPHI_CUT){
     // if (W_vector.DeltaR(fTwoProngLoose) > 0.3){

         fTwoProng_pt_sideband->Fill(fTwoProngLoose.Pt(), eventWeight);
         Get_Mass_PT_plot_side(fTwoProngLoose.Pt())->Fill(TwoProngLoose_MassPi0[0], get_sample_num(eventWeight));
         //fW_pt_sideband->Fill(fW.Pt(), eventWeight);
         //fW_mass_sideband->Fill(fW.M(), eventWeight);
         fIDPhoton_pt_sideband->Fill(fIDPhoton.Pt(), eventWeight);
         fHT_sideband->Fill(*HT, eventWeight);
         fHTbare_sideband->Fill(*HT_bare, eventWeight);
         fTwoProng_eta_sideband->Fill(TwoProngLoose_eta[0], eventWeight);
         fTwoProng_phi_sideband->Fill(TwoProngLoose_phi[0], eventWeight);
         fTwoProng_MassPi0_sideband->Fill(TwoProngLoose_MassPi0[0], eventWeight);
         fIDPhoton_pt_sideband->Fill(Photon_pt[0], eventWeight);
         fIDPhoton_phi_sideband->Fill(Photon_phi[0], eventWeight);
         fIDPhoton_eta_sideband->Fill(Photon_eta[0], eventWeight);
         double Phi_Mass = (gamma_vector + fTwoProngLoose).M();
 
        // fW_eta_sideband->Fill(W_eta[0], eventWeight);
        // fW_phi_sideband->Fill(W_phi[0], eventWeight);
         fTwoProng_MassPi0_PT_sideband->Fill(TwoProngLoose_MassPi0[0], TwoProngLoose_pt[0], eventWeight);
         fPhiMass_PT_sideband->Fill(Phi_Mass,TwoProngLoose_pt[0], eventWeight );
   //   }
     }
   }


   //fW_pt_distr->Fill(fW.Pt(), eventWeight);
   //fW_mass_distr->Fill(fW.M(), eventWeight);
   fHT_distr->Fill(*HT, eventWeight);
   fHTbare_distr->Fill(*HT_bare, eventWeight);
   fIDPhoton_pt_distr->Fill(Photon_pt[0], eventWeight);
   fIDPhoton_phi_distr->Fill(Photon_phi[0], eventWeight);
   fIDPhoton_eta_distr->Fill(Photon_eta[0], eventWeight);

   //fW_eta_distr->Fill(W_eta[0], eventWeight);
   //fW_phi_distr->Fill(W_phi[0], eventWeight);

//   if (Imag_W[0] == 1.0) fW_mass_distr_Imag_W->Fill(fW.M(), eventWeight); 

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

  // MET_pt_distr->Fill(MET_pt[0], eventWeight);
  // MET_phi_distr->Fill(MET_phi[0], eventWeight);

 //  mT_distr->Fill(*mT, eventWeight);


   return kTRUE;
}

void GJetsReader::Terminate()
{
 // std::cout << "IN END LOOP" << std::endl; 

}

TH2* GJetsReader::Get_Mass_PT_plot_tight(double pt_val){
   for (int j = 0; j < PT_BOUNDS1.size() - 1; j++){
      if (pt_val > PT_BOUNDS1[j] and pt_val < PT_BOUNDS1[j + 1]) return get_plot(j);
   }
   if (pt_val > PT_BOUNDS1[PT_BOUNDS1.size()-1]) return PT_Mass_Pi0_10;
   else return nullptr;
}

TH2* GJetsReader::Get_Mass_PT_plot_side(double pt_val){
   for (int j = 0; j < PT_BOUNDS1.size() - 1; j++){
      if (pt_val > PT_BOUNDS1[j] and pt_val < PT_BOUNDS1[j + 1]) return get_plot_side(j);
   }
   if (pt_val > PT_BOUNDS1[PT_BOUNDS1.size()-1]) return PT_Mass_Pi0_10_side;
   else return nullptr;
}

TH2* GJetsReader::get_plot(int num){

   if (num == 0) return GJetsReader::PT_Mass_Pi0_01;
   if (num == 1) return GJetsReader::PT_Mass_Pi0_02;
   if (num == 2) return GJetsReader::PT_Mass_Pi0_03;
   if (num == 3) return GJetsReader::PT_Mass_Pi0_04;
   if (num == 4) return GJetsReader::PT_Mass_Pi0_05;
   if (num == 5) return GJetsReader::PT_Mass_Pi0_06;
   if (num == 6) return GJetsReader::PT_Mass_Pi0_07;
   if (num == 7) return GJetsReader::PT_Mass_Pi0_08;
   if (num == 8) return GJetsReader::PT_Mass_Pi0_09;
   if (num == 9) return GJetsReader::PT_Mass_Pi0_10;

   return nullptr;
}

TH2* GJetsReader::get_plot_side(int num){

   if (num == 0) return GJetsReader::PT_Mass_Pi0_01_side;
   if (num == 1) return GJetsReader::PT_Mass_Pi0_02_side;
   if (num == 2) return GJetsReader::PT_Mass_Pi0_03_side;
   if (num == 3) return GJetsReader::PT_Mass_Pi0_04_side;
   if (num == 4) return GJetsReader::PT_Mass_Pi0_05_side;
   if (num == 5) return GJetsReader::PT_Mass_Pi0_06_side;
   if (num == 6) return GJetsReader::PT_Mass_Pi0_07_side;
   if (num == 7) return GJetsReader::PT_Mass_Pi0_08_side;
   if (num == 8) return GJetsReader::PT_Mass_Pi0_09_side;
   if (num == 9) return GJetsReader::PT_Mass_Pi0_10_side;

   return nullptr;
}

int GJetsReader::get_sample_num(double weight){
   for (int j = 1; j <= GJetsReader::N_SAMPLES_PASSED; j++){
      if ( GJetsReader::sample_weights->GetBinContent(j) == weight ) return j;
   }
   GJetsReader::N_SAMPLES_PASSED++;
   GJetsReader::sample_weights->SetBinContent(N_SAMPLES_PASSED, weight);
   return GJetsReader::N_SAMPLES_PASSED;
}


#endif