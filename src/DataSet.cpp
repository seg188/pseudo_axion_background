#include "DataSet.hh"
#include "../util/Traverse.cpp"
#include "online/WJetsReader.cpp"
#include "online/GJetsReader.cpp"

using namespace online;
using namespace histogram;

void WJetsDataSet::GetHistograms(){
	std::cout << "allocating histograms!!!! " << std::endl;
	global_tags::ReaderTag = "online/WJetsReader.cpp";
	global_tags::tree_tag = "Analyzer/fTree";

	_fTwoProng_pt_signal = new TH1D("fTwoProng_pt_signal", "fTwoProng_pt_signal", PT_NUMBINS, PT_MIN, PT_MAX);
	_fTwoProng_pt_sideband = new TH1D("fTwoProng_pt_sideband", "fTwoProng_pt_sideband", PT_NUMBINS, PT_MIN, PT_MAX);
	_fTwoProng_pt_distr = new TH1D("fTwoProng_pt_distr", "fTwoProng_pt_distr", PT_NUMBINS, PT_MIN, PT_MAX);

	_fW_pt_signal = new TH1D("fW_pt_signal", "fW_pt_signal", PT_NUMBINS, PT_MIN, PT_MAX);
	_fW_pt_sideband = new TH1D("fW_pt_sideband", "fW_pt_sideband", PT_NUMBINS, PT_MIN, PT_MAX);
	_fW_pt_distr = new TH1D("fW_pt_distr", "fW_pt_distr", PT_NUMBINS, PT_MIN, PT_MAX);

	_fW_mass_signal = new TH1D("fW_mass_signal", "fW_mass_signal", W_MASS_NUMBINS, W_MASS_MIN, W_MASS_MAX);
	_fW_mass_sideband = new TH1D("fW_mass_sideband", "fW_mass_sideband", W_MASS_NUMBINS, W_MASS_MIN, W_MASS_MAX);
	_fW_mass_distr = new TH1D("fW_mass_distr", "fW_mass_distr", W_MASS_NUMBINS, W_MASS_MIN, W_MASS_MAX);
	_fW_mass_distr_Imag_W = new TH1D("fW_mass_distr_Imag_W", "fW_mass_distr", W_MASS_NUMBINS, W_MASS_MIN, W_MASS_MAX);

	_fMuon_pt_signal = new TH1D("fMuon_pt_signal", "fMuon_pt_signal", PT_NUMBINS, PT_MIN, PT_MAX);
	_fMuon_pt_sideband = new TH1D("fMuon_pt_sideband", "fMuon_pt_sideband", PT_NUMBINS, PT_MIN, PT_MAX);
	_fMuon_pt_distr = new TH1D("fMuon_pt_distr", "fMuon_pt_distr", PT_NUMBINS, PT_MIN, PT_MAX);

	_fMuon_eta_signal = new TH1D("fMuon_eta_signal", "fMuon_eta_signal", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);
	_fMuon_eta_sideband = new TH1D("fMuon_eta_sideband", "fMuon_eta_sideband", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);
	_fMuon_eta_distr = new TH1D("fMuon_eta_distr", "fMuon_eta_distr", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);

	_fMuon_phi_signal = new TH1D("fMuon_phi_signal", "fMuon_phi_signal", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);
	_fMuon_phi_sideband = new TH1D("fMuon_phi_sideband", "fMuon_phi_sideband", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);
	_fMuon_phi_distr = new TH1D("fMuon_phi_distr", "fMuon_phi_distr", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);

	_fTightMuon_pt_signal = new TH1D("fTightMuon_pt_signal", "fMuon_pt_signal", PT_NUMBINS, PT_MIN, PT_MAX);
	_fTightMuon_pt_sideband = new TH1D("fTightMuon_pt_sideband", "fMuon_pt_sideband", PT_NUMBINS, PT_MIN, PT_MAX);
	_fTightMuon_pt_distr = new TH1D("fTightMuon_pt_distr", "fMuon_pt_distr", PT_NUMBINS, PT_MIN, PT_MAX);

	_fTightMuon_eta_signal = new TH1D("fTightMuon_eta_signal", "fTightMuon_eta_signal", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);
	_fTightMuon_eta_sideband = new TH1D("fTightMuon_eta_sideband", "fTightMuon_eta_sideband", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);
	_fTightMuon_eta_distr = new TH1D("fTightMuon_eta_distr", "fTightMuon_eta_distr", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);

	_fTightMuon_phi_signal = new TH1D("fTightMuon_phi_signal", "fTightMuon_phi_signal", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);
	_fTightMuon_phi_sideband = new TH1D("fTightMuon_phi_sideband", "fTightMuon_phi_sideband", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);
	_fTightMuon_phi_distr = new TH1D("fTightMuon_phi_distr", "fTightMuon_phi_distr", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);

 	_fHT_distr = new TH1D("fHT_distr", "fHT_distr", HT_NUMBINS, HT_MIN, HT_MAX);
 	_fHT_signal = new TH1D("fHT_signal", "fHT_distr", HT_NUMBINS, HT_MIN, HT_MAX);
 	_fHT_sideband = new TH1D("fHT_sideband", "fHT_distr", HT_NUMBINS, HT_MIN, HT_MAX);

 	_fTwoProng_eta_signal = new TH1D("fTwoProng_eta_signal", "fTwoProng_eta_signal", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);
 	_fTwoProng_eta_sideband = new TH1D("fTwoProng_eta_sideband", "fTwoProng_eta_sideband", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);
 	_fTwoProng_eta_distr = new TH1D("fTwoProng_eta_distr", "fTwoProng_eta_distr", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);

 	_fTwoProng_phi_signal = new TH1D("fTwoProng_phi_signal", "fTwoProng_phi_signal", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);
 	_fTwoProng_phi_sideband  = new TH1D("fTwoProng_phi_sideband", "fTwoProng_phi_sideband", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);
 	_fTwoProng_phi_distr = new TH1D("fTwoProng_phi_distr", "fTwoProng_phi_distr", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);

 	_fTwoProng_MassPi0_signal = new TH1D("fTwoProng_MassPi0_signal", "fTwoProng_MassPi0_signal", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN);
 	_fTwoProng_MassPi0_sideband  = new TH1D("fTwoProng_MassPi0_sideband", "fTwoProng_MassPi0_sideband", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN);
 	_fTwoProng_MassPi0_distr  = new TH1D("fTwoProng_MassPi0_distr", "fTwoProng_MassPi0_distr", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN);

 	_fW_eta_signal = new TH1D("fW_eta_signal", "fW_eta_signal", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);
 	_fW_eta_sideband = new TH1D("fW_eta_sideband", "fW_eta_sideband", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);
 	_fW_eta_distr = new TH1D("fW_eta_distr", "fW_eta_distr", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);

 	_fW_phi_signal = new TH1D("fW_phi_signal", "fW_phi_signal", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);
 	_fW_phi_sideband = new TH1D("fW_phi_sideband", "fW_phi_sideband", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);
 	_fW_phi_distr = new TH1D("fW_phi_distr", "fW_phi_distr", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);

	_fHTbare_distr = new TH1D("fHTBare_distr", "fHTBare_distr", HT_NUMBINS, HT_MIN, HT_MAX);
	_fHTbare_signal = new TH1D("fHTBare_signal", "fHTBare_signal", HT_NUMBINS, HT_MIN, HT_MAX);
	_fHTbare_sideband = new TH1D("fHTBare_sideband", "fHTBare_sideband", HT_NUMBINS, HT_MIN, HT_MAX);

	_MET_phi_distr= new TH1D("MET_phi_distr", "MET_phi_distrl", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);
	_MET_pt_distr = new TH1D("MET_pt_distr", "MET_pt_distr", PT_NUMBINS, PT_MIN, PT_MAX);

	_mT_distr = new TH1D("mT_distr", "mT_distr", MT_NUMBINS, MT_MIN, MT_MAX);
	_mT_signal = new TH1D("mT_signal", "mT_signal", MT_NUMBINS, MT_MIN, MT_MAX);

	_mu_w_deltaphi = new TH1D("Muon_W_dphi_distr", "Muon_W_dphi_distr", DPHI_NUMBINS, DPHI_MINBIN, DPHI_MAXBIN);
	_mu_met_deltaphi = new TH1D("Muon_met_dphi_distr", "Muon_met_dphi_distr", DPHI_NUMBINS, DPHI_MINBIN, DPHI_MAXBIN);

	_fTwoProng_MassPi0_PT = new TH2D("Mass_pi0_pt", "TwoProng_mass_pt", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, PT_NUMBINS, PT_MIN, PT_MAX);
	_fTwoProng_MassPi0_PT_sideband = new TH2D("Mass_pi0_pt_sideband", "TwoProng_mass_pt_sideband", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, PT_NUMBINS, PT_MIN, PT_MAX);

	_fTwoProng_MassPi0_PT_lowPTmuon = new TH2D("Mass_pi0_pt_lowPTmuon", "TwoProng_mass_pt", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, PT_NUMBINS, PT_MIN, PT_MAX);
	_fTwoProng_MassPi0_PT_sideband_lowPTmuon = new TH2D("Mass_pi0_pt_sideband_lowPTmuon", "TwoProng_mass_pt_sideband", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, PT_NUMBINS, PT_MIN, PT_MAX);

	_fTwoProng_MassPi0_PT_highPTmuon = new TH2D("Mass_pi0_pt_highPTmuon", "TwoProng_mass_pt", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, PT_NUMBINS, PT_MIN, PT_MAX);
	_fTwoProng_MassPi0_PT_sideband_highPTmuon = new TH2D("Mass_pi0_pt_sideband_highPTmuon", "TwoProng_mass_pt_sideband", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, PT_NUMBINS, PT_MIN, PT_MAX);

	_PT_Mass_Pi0_01 = new TH2D("PT_Mass_pi0_01", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_02 = new TH2D("PT_Mass_pi0_02", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_03 = new TH2D("PT_Mass_pi0_03", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_04 = new TH2D("PT_Mass_pi0_04", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_05 = new TH2D("PT_Mass_pi0_05", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_06 = new TH2D("PT_Mass_pi0_06", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_07 = new TH2D("PT_Mass_pi0_07", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_08 = new TH2D("PT_Mass_pi0_08", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_09 = new TH2D("PT_Mass_pi0_09", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_10 = new TH2D("PT_Mass_pi0_10", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);

	_PT_Mass_Pi0_01_side = new TH2D("PT_Mass_pi0_01_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_02_side = new TH2D("PT_Mass_pi0_02_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_03_side = new TH2D("PT_Mass_pi0_03_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_04_side = new TH2D("PT_Mass_pi0_04_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_05_side = new TH2D("PT_Mass_pi0_05_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_06_side = new TH2D("PT_Mass_pi0_06_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_07_side = new TH2D("PT_Mass_pi0_07_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_08_side = new TH2D("PT_Mass_pi0_08_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_09_side = new TH2D("PT_Mass_pi0_09_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
	_PT_Mass_Pi0_10_side = new TH2D("PT_Mass_pi0_10_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);

	_fW_pt_signal_noW = new TH1D("fW_pt_signal_noW", "fW_pt_signal_noW_cut", PT_NUMBINS, PT_MIN, PT_MAX);
	_fW_pt_signal_highptu = new TH1D("fW_pt_highptu", "fW_pt_signal", PT_NUMBINS, PT_MIN, PT_MAX);

	_sample_weights = new TH1D("sample_weights", "weights_by_sample_int", 10, 0.5, 10.5);	

}

void WJetsDataSet::InitializeHistograms(){
	std::cout << "setting static histograms in reader" << std::endl;
    WJetsReader::fTwoProng_pt_signal = _fTwoProng_pt_signal;
    WJetsReader::fTwoProng_pt_sideband = _fTwoProng_pt_sideband;
    WJetsReader::fTwoProng_pt_distr = _fTwoProng_pt_distr;

    WJetsReader::fTwoProng_eta_signal = _fTwoProng_eta_signal;
    WJetsReader::fTwoProng_eta_sideband = _fTwoProng_eta_sideband;
    WJetsReader::fTwoProng_eta_distr = _fTwoProng_eta_distr;

    WJetsReader::fTwoProng_phi_signal = _fTwoProng_phi_signal;
    WJetsReader::fTwoProng_phi_sideband = _fTwoProng_phi_sideband;
    WJetsReader::fTwoProng_phi_distr = _fTwoProng_phi_distr;

    WJetsReader::fTwoProng_MassPi0_signal = _fTwoProng_MassPi0_signal;
    WJetsReader::fTwoProng_MassPi0_sideband = _fTwoProng_MassPi0_sideband;
    WJetsReader::fTwoProng_MassPi0_distr = _fTwoProng_MassPi0_distr;

    WJetsReader::fW_pt_signal = _fW_pt_signal;
    WJetsReader::fW_pt_signal_noW = _fW_pt_signal_noW;
    WJetsReader::fW_pt_signal_highptu = _fW_pt_signal_highptu;

    WJetsReader::fW_pt_sideband = _fW_pt_sideband;
    WJetsReader::fW_pt_distr = _fW_pt_distr;

    WJetsReader::fW_eta_signal = _fW_eta_signal;
    WJetsReader::fW_eta_sideband = _fW_eta_sideband;
    WJetsReader::fW_eta_distr = _fW_eta_distr;

    WJetsReader::fW_phi_signal = _fW_phi_signal;
    WJetsReader::fW_phi_sideband = _fW_phi_sideband;
    WJetsReader::fW_phi_distr = _fW_phi_distr;

    WJetsReader::fW_mass_signal = _fW_mass_signal;
    WJetsReader::fW_mass_sideband = _fW_mass_sideband;
    WJetsReader::fW_mass_distr = _fW_mass_distr;

    WJetsReader::fMuon_pt_signal = _fMuon_pt_signal;
    WJetsReader::fMuon_pt_sideband = _fMuon_pt_sideband;
    WJetsReader::fMuon_pt_distr = _fMuon_pt_distr;

    WJetsReader::fMuon_eta_signal = _fMuon_eta_signal;
    WJetsReader::fMuon_eta_sideband = _fMuon_eta_sideband;
    WJetsReader::fMuon_eta_distr = _fMuon_eta_distr;

    WJetsReader::fMuon_phi_signal = _fMuon_phi_signal;
    WJetsReader::fMuon_phi_sideband = _fMuon_phi_sideband;
    WJetsReader::fMuon_phi_distr = _fMuon_phi_distr;

    WJetsReader::fTightMuon_pt_signal = _fTightMuon_pt_signal;
    WJetsReader::fTightMuon_pt_sideband = _fTightMuon_pt_sideband;
    WJetsReader::fTightMuon_pt_distr = _fTightMuon_pt_distr;

    WJetsReader::fTightMuon_eta_signal = _fTightMuon_eta_signal;
    WJetsReader::fTightMuon_eta_sideband = _fTightMuon_eta_sideband;
    WJetsReader::fTightMuon_eta_distr = _fTightMuon_eta_distr;

    WJetsReader::fTightMuon_phi_signal = _fTightMuon_phi_signal;
    WJetsReader::fTightMuon_phi_sideband = _fTightMuon_phi_sideband;
    WJetsReader::fTightMuon_phi_distr = _fTightMuon_phi_distr;

    WJetsReader::fHT_distr = _fHT_distr;
    WJetsReader::fHT_signal = _fHT_signal;
    WJetsReader::fHT_sideband = _fHT_sideband;

    WJetsReader::fHTbare_distr = _fHTbare_distr;
    WJetsReader::fHTbare_signal = _fHTbare_signal;
    WJetsReader::fHTbare_sideband = _fHTbare_sideband;

    WJetsReader::MET_pt_distr = _MET_pt_distr;
    WJetsReader::MET_phi_distr = _MET_phi_distr; 

    WJetsReader::mT_distr = _mT_distr;
    WJetsReader::mT_signal = _mT_signal;
   
    WJetsReader::mu_w_deltaphi = _mu_w_deltaphi;
    WJetsReader::mu_met_deltaphi = _mu_met_deltaphi;

    WJetsReader::fTwoProng_MassPi0_PT = _fTwoProng_MassPi0_PT;
    WJetsReader::fTwoProng_MassPi0_PT_sideband = _fTwoProng_MassPi0_PT_sideband;

    WJetsReader::fTwoProng_MassPi0_PT_lowPTmuon = _fTwoProng_MassPi0_PT_lowPTmuon;
    WJetsReader::fTwoProng_MassPi0_PT_sideband_lowPTmuon = _fTwoProng_MassPi0_PT_sideband_lowPTmuon;

    WJetsReader::fTwoProng_MassPi0_PT_highPTmuon = _fTwoProng_MassPi0_PT_highPTmuon;
    WJetsReader::fTwoProng_MassPi0_PT_sideband_highPTmuon = _fTwoProng_MassPi0_PT_sideband_highPTmuon;

    WJetsReader::PT_Mass_Pi0_01 = _PT_Mass_Pi0_01;
    WJetsReader::PT_Mass_Pi0_02 = _PT_Mass_Pi0_02;
    WJetsReader::PT_Mass_Pi0_03 = _PT_Mass_Pi0_03;
    WJetsReader::PT_Mass_Pi0_04 = _PT_Mass_Pi0_04;
    WJetsReader::PT_Mass_Pi0_05 = _PT_Mass_Pi0_05;
    WJetsReader::PT_Mass_Pi0_06 = _PT_Mass_Pi0_06;
    WJetsReader::PT_Mass_Pi0_07 = _PT_Mass_Pi0_07;
    WJetsReader::PT_Mass_Pi0_08 = _PT_Mass_Pi0_08;
    WJetsReader::PT_Mass_Pi0_09 = _PT_Mass_Pi0_09;
    WJetsReader::PT_Mass_Pi0_10 = _PT_Mass_Pi0_10;

    WJetsReader::PT_Mass_Pi0_01_side = _PT_Mass_Pi0_01_side;
    WJetsReader::PT_Mass_Pi0_02_side = _PT_Mass_Pi0_02_side;
    WJetsReader::PT_Mass_Pi0_03_side = _PT_Mass_Pi0_03_side;
    WJetsReader::PT_Mass_Pi0_04_side = _PT_Mass_Pi0_04_side;
    WJetsReader::PT_Mass_Pi0_05_side = _PT_Mass_Pi0_05_side;
    WJetsReader::PT_Mass_Pi0_06_side = _PT_Mass_Pi0_06_side;
    WJetsReader::PT_Mass_Pi0_07_side = _PT_Mass_Pi0_07_side;
    WJetsReader::PT_Mass_Pi0_08_side = _PT_Mass_Pi0_08_side;
    WJetsReader::PT_Mass_Pi0_09_side = _PT_Mass_Pi0_09_side;
    WJetsReader::PT_Mass_Pi0_10_side = _PT_Mass_Pi0_10_side;

    WJetsReader::sample_weights = _sample_weights;
    WJetsReader::N_SAMPLES_PASSED = 0;

}


void GJetsDataSet::GetHistograms(){
   std::cout << "allocating histograms!!!! " << std::endl;
   global_tags::ReaderTag = "online/GJetsReader.cpp";
   global_tags::tree_tag = "twoprongNtuplizer/fTree";

   _fTwoProng_pt_signal = new TH1D("fTwoProng_pt_signal", "fTwoProng_pt_signal", PT_NUMBINS, PT_MIN, PT_MAX);
   _fTwoProng_pt_sideband = new TH1D("fTwoProng_pt_sideband", "fTwoProng_pt_sideband", PT_NUMBINS, PT_MIN, PT_MAX);
   _fTwoProng_pt_distr = new TH1D("fTwoProng_pt_distr", "fTwoProng_pt_distr", PT_NUMBINS, PT_MIN, PT_MAX);

   _fPhi_pt_signal = new TH1D("fPhi_pt_signal", "fPhi_pt_signal", PT_NUMBINS, PT_MIN, PT_MAX);
   _fPhi_pt_sideband = new TH1D("fPhi_pt_sideband", "fPhi_pt_sideband", PT_NUMBINS, PT_MIN, PT_MAX);
   _fPhi_pt_distr = new TH1D("fPhi_pt_distr", "fPhi_pt_distr", PT_NUMBINS, PT_MIN, PT_MAX);

   _fPhi_mass_signal = new TH1D("fPhi_mass_signal", "fPhi_mass_signal", PHI_MASS_NUMBINS, PHI_MASS_MIN, PHI_MASS_MAX);
   _fPhi_mass_sideband = new TH1D("fPhi_mass_sideband", "fPhi_mass_sideband", PHI_MASS_NUMBINS, PHI_MASS_MIN, PHI_MASS_MAX);
   _fPhi_mass_distr = new TH1D("fPhi_mass_distr", "fPhi_mass_distr", PHI_MASS_NUMBINS, PHI_MASS_MIN, PHI_MASS_MAX);

   _fIDPhoton_pt_signal = new TH1D("fIDPhoton_pt_signal", "fIDPhoton_pt_signal", PT_NUMBINS, PT_MIN, PT_MAX);
   _fIDPhoton_pt_sideband = new TH1D("fIDPhoton_pt_sideband", "fIDPhoton_pt_sideband", PT_NUMBINS, PT_MIN, PT_MAX);
   _fIDPhoton_pt_distr = new TH1D("fIDPhoton_pt_distr", "fIDPhoton_pt_distr", PT_NUMBINS, PT_MIN, PT_MAX);

   _fIDPhoton_eta_signal = new TH1D("fIDPhoton_eta_signal", "fIDPhoton_eta_signal", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);
   _fIDPhoton_eta_sideband = new TH1D("fIDPhoton_eta_sideband", "fIDPhoton_eta_sideband", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);
   _fIDPhoton_eta_distr = new TH1D("fIDPhoton_eta_distr", "fIDPhoton_eta_distr", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);

   _fIDPhoton_phi_signal = new TH1D("fIDPhoton_phi_signal", "fIDPhoton_phi_signal", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);
   _fIDPhoton_phi_sideband = new TH1D("fIDPhoton_phi_sideband", "fIDPhoton_phi_sideband", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);
   _fIDPhoton_phi_distr = new TH1D("fIDPhoton_phi_distr", "fIDPhoton_phi_distr", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);

   _fHT_distr = new TH1D("fHT_distr", "fHT_distr", HT_NUMBINS, HT_MIN, HT_MAX);
   _fHT_signal = new TH1D("fHT_signal", "fHT_distr", HT_NUMBINS, HT_MIN, HT_MAX);
   _fHT_sideband = new TH1D("fHT_sideband", "fHT_distr", HT_NUMBINS, HT_MIN, HT_MAX);

   _fTwoProng_eta_signal = new TH1D("fTwoProng_eta_signal", "fTwoProng_eta_signal", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);
   _fTwoProng_eta_sideband = new TH1D("fTwoProng_eta_sideband", "fTwoProng_eta_sideband", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);
   _fTwoProng_eta_distr = new TH1D("fTwoProng_eta_distr", "fTwoProng_eta_distr", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);

   _fTwoProng_phi_signal = new TH1D("fTwoProng_phi_signal", "fTwoProng_phi_signal", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);
   _fTwoProng_phi_sideband  = new TH1D("fTwoProng_phi_sideband", "fTwoProng_phi_sideband", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);
   _fTwoProng_phi_distr = new TH1D("fTwoProng_phi_distr", "fTwoProng_phi_distr", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);

   _fTwoProng_MassPi0_signal = new TH1D("fTwoProng_MassPi0_signal", "fTwoProng_MassPi0_signal", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN);
   _fTwoProng_MassPi0_sideband  = new TH1D("fTwoProng_MassPi0_sideband", "fTwoProng_MassPi0_sideband", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN);
   _fTwoProng_MassPi0_distr  = new TH1D("fTwoProng_MassPi0_distr", "fTwoProng_MassPi0_distr", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN);

   _fPhi_eta_signal = new TH1D("fPhi_eta_signal", "fPhi_eta_signal", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);
   _fPhi_eta_sideband = new TH1D("fPhi_eta_sideband", "fPhi_eta_sideband", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);
   _fPhi_eta_distr = new TH1D("fPhi_eta_distr", "fPhi_eta_distr", ETA_NUMBINS, ETA_MINBIN, ETA_MAXBIN);

   _fPhi_phi_signal = new TH1D("fPhi_phi_signal", "fPhi_phi_signal", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);
   _fPhi_phi_sideband = new TH1D("fPhi_phi_sideband", "fPhi_phi_sideband", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);
   _fPhi_phi_distr = new TH1D("fPhi_phi_distr", "fPhi_phi_distr", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);

   _fHTbare_distr = new TH1D("fHTBare_distr", "fHTBare_distr", HT_NUMBINS, HT_MIN, HT_MAX);
   _fHTbare_signal = new TH1D("fHTBare_signal", "fHTBare_signal", HT_NUMBINS, HT_MIN, HT_MAX);
   _fHTbare_sideband = new TH1D("fHTBare_sideband", "fHTBare_sideband", HT_NUMBINS, HT_MIN, HT_MAX);

   _MET_phi_distr= new TH1D("MET_phi_distr", "MET_phi_distrl", PHI_NUMBINS, PHI_MINBIN, PHI_MAXBIN);
   _MET_pt_distr = new TH1D("MET_pt_distr", "MET_pt_distr", PT_NUMBINS, PT_MIN, PT_MAX);

   _fTwoProng_MassPi0_PT = new TH2D("Mass_pi0_pt", "TwoProng_mass_pt", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, PT_NUMBINS, PT_MIN, PT_MAX);
   _fTwoProng_MassPi0_PT_sideband = new TH2D("Mass_pi0_pt_sideband", "TwoProng_mass_pt_sideband", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, PT_NUMBINS, PT_MIN, PT_MAX);

   _fPhiMass_PT_sideband = new TH2D("fPhiMass_PT_sideband", "fPhiMass_PT_sideband", PHI_MASS_NUMBINS, PHI_MASS_MIN, PHI_MASS_MAX, PT_BOUNDS1.size(), PT_BOUNDS);
   _fPhiMass_PT_signal = new TH2D("fPhiMass_PT_signal", "fPhiMass_PT_signal", PHI_MASS_NUMBINS, PHI_MASS_MIN, PHI_MASS_MAX, PT_BOUNDS1.size(), PT_BOUNDS);

   _PT_Mass_Pi0_01 = new TH2D("PT_Mass_pi0_01", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_02 = new TH2D("PT_Mass_pi0_02", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_03 = new TH2D("PT_Mass_pi0_03", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_04 = new TH2D("PT_Mass_pi0_04", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_05 = new TH2D("PT_Mass_pi0_05", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_06 = new TH2D("PT_Mass_pi0_06", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_07 = new TH2D("PT_Mass_pi0_07", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_08 = new TH2D("PT_Mass_pi0_08", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_09 = new TH2D("PT_Mass_pi0_09", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_10 = new TH2D("PT_Mass_pi0_10", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);

   _PT_Mass_Pi0_01_side = new TH2D("PT_Mass_pi0_01_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_02_side = new TH2D("PT_Mass_pi0_02_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_03_side = new TH2D("PT_Mass_pi0_03_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_04_side = new TH2D("PT_Mass_pi0_04_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_05_side = new TH2D("PT_Mass_pi0_05_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_06_side = new TH2D("PT_Mass_pi0_06_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_07_side = new TH2D("PT_Mass_pi0_07_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_08_side = new TH2D("PT_Mass_pi0_08_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_09_side = new TH2D("PT_Mass_pi0_09_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);
   _PT_Mass_Pi0_10_side = new TH2D("PT_Mass_pi0_10_side", "plots_w_weight", TP_MASS_NUMBINS, TP_MASS_MINBIN, TP_MASS_MAXBIN, 10, 0.5000, 10.500);


   _sample_weights = new TH1D("sample_weights", "weights_by_sample_int", 10, 0.5, 10.5);  

}

void GJetsDataSet::InitializeHistograms(){
   std::cout << "setting static histograms in reader" << std::endl;
    GJetsReader::fTwoProng_pt_signal = _fTwoProng_pt_signal;
    GJetsReader::fTwoProng_pt_sideband = _fTwoProng_pt_sideband;
    GJetsReader::fTwoProng_pt_distr = _fTwoProng_pt_distr;

    GJetsReader::fTwoProng_eta_signal = _fTwoProng_eta_signal;
    GJetsReader::fTwoProng_eta_sideband = _fTwoProng_eta_sideband;
    GJetsReader::fTwoProng_eta_distr = _fTwoProng_eta_distr;

    GJetsReader::fTwoProng_phi_signal = _fTwoProng_phi_signal;
    GJetsReader::fTwoProng_phi_sideband = _fTwoProng_phi_sideband;
    GJetsReader::fTwoProng_phi_distr = _fTwoProng_phi_distr;

    GJetsReader::fTwoProng_MassPi0_signal = _fTwoProng_MassPi0_signal;
    GJetsReader::fTwoProng_MassPi0_sideband = _fTwoProng_MassPi0_sideband;
    GJetsReader::fTwoProng_MassPi0_distr = _fTwoProng_MassPi0_distr;

    GJetsReader::fPhi_pt_signal = _fPhi_pt_signal;

    GJetsReader::fPhi_pt_sideband = _fPhi_pt_sideband;
    GJetsReader::fPhi_pt_distr = _fPhi_pt_distr;

    GJetsReader::fPhi_eta_signal = _fPhi_eta_signal;
    GJetsReader::fPhi_eta_sideband = _fPhi_eta_sideband;
    GJetsReader::fPhi_eta_distr = _fPhi_eta_distr;

    GJetsReader::fPhi_phi_signal = _fPhi_phi_signal;
    GJetsReader::fPhi_phi_sideband = _fPhi_phi_sideband;
    GJetsReader::fPhi_phi_distr = _fPhi_phi_distr;

    GJetsReader::fPhi_mass_signal = _fPhi_mass_signal;
    GJetsReader::fPhi_mass_sideband = _fPhi_mass_sideband;
    GJetsReader::fPhi_mass_distr = _fPhi_mass_distr;

    GJetsReader::fIDPhoton_pt_signal = _fIDPhoton_pt_signal;
    GJetsReader::fIDPhoton_pt_sideband = _fIDPhoton_pt_sideband;
    GJetsReader::fIDPhoton_pt_distr = _fIDPhoton_pt_distr;

    GJetsReader::fIDPhoton_eta_signal = _fIDPhoton_eta_signal;
    GJetsReader::fIDPhoton_eta_sideband = _fIDPhoton_eta_sideband;
    GJetsReader::fIDPhoton_eta_distr = _fIDPhoton_eta_distr;

    GJetsReader::fIDPhoton_phi_signal = _fIDPhoton_phi_signal;
    GJetsReader::fIDPhoton_phi_sideband = _fIDPhoton_phi_sideband;
    GJetsReader::fIDPhoton_phi_distr = _fIDPhoton_phi_distr;


    GJetsReader::fHT_distr = _fHT_distr;
    GJetsReader::fHT_signal = _fHT_signal;
    GJetsReader::fHT_sideband = _fHT_sideband;

    GJetsReader::fHTbare_distr = _fHTbare_distr;
    GJetsReader::fHTbare_signal = _fHTbare_signal;
    GJetsReader::fHTbare_sideband = _fHTbare_sideband;

    GJetsReader::MET_pt_distr = _MET_pt_distr;
    GJetsReader::MET_phi_distr = _MET_phi_distr; 
  

    GJetsReader::fTwoProng_MassPi0_PT = _fTwoProng_MassPi0_PT;
    GJetsReader::fTwoProng_MassPi0_PT_sideband = _fTwoProng_MassPi0_PT_sideband;

    GJetsReader::fPhiMass_PT_sideband = _fPhiMass_PT_sideband;
    GJetsReader::fPhiMass_PT_signal = _fPhiMass_PT_signal;

    GJetsReader::PT_Mass_Pi0_01 = _PT_Mass_Pi0_01;
    GJetsReader::PT_Mass_Pi0_02 = _PT_Mass_Pi0_02;
    GJetsReader::PT_Mass_Pi0_03 = _PT_Mass_Pi0_03;
    GJetsReader::PT_Mass_Pi0_04 = _PT_Mass_Pi0_04;
    GJetsReader::PT_Mass_Pi0_05 = _PT_Mass_Pi0_05;
    GJetsReader::PT_Mass_Pi0_06 = _PT_Mass_Pi0_06;
    GJetsReader::PT_Mass_Pi0_07 = _PT_Mass_Pi0_07;
    GJetsReader::PT_Mass_Pi0_08 = _PT_Mass_Pi0_08;
    GJetsReader::PT_Mass_Pi0_09 = _PT_Mass_Pi0_09;
    GJetsReader::PT_Mass_Pi0_10 = _PT_Mass_Pi0_10;

    GJetsReader::PT_Mass_Pi0_01_side = _PT_Mass_Pi0_01_side;
    GJetsReader::PT_Mass_Pi0_02_side = _PT_Mass_Pi0_02_side;
    GJetsReader::PT_Mass_Pi0_03_side = _PT_Mass_Pi0_03_side;
    GJetsReader::PT_Mass_Pi0_04_side = _PT_Mass_Pi0_04_side;
    GJetsReader::PT_Mass_Pi0_05_side = _PT_Mass_Pi0_05_side;
    GJetsReader::PT_Mass_Pi0_06_side = _PT_Mass_Pi0_06_side;
    GJetsReader::PT_Mass_Pi0_07_side = _PT_Mass_Pi0_07_side;
    GJetsReader::PT_Mass_Pi0_08_side = _PT_Mass_Pi0_08_side;
    GJetsReader::PT_Mass_Pi0_09_side = _PT_Mass_Pi0_09_side;
    GJetsReader::PT_Mass_Pi0_10_side = _PT_Mass_Pi0_10_side;

    GJetsReader::sample_weights = _sample_weights;
    GJetsReader::N_SAMPLES_PASSED = 0;

}