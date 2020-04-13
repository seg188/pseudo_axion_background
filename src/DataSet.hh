#include <iostream>
#include "../include/UserIn.cpp"
#include "../util/Traverse.cpp"

#ifndef DATA_SET_CLASS
#define DATA_SET_CLASS

namespace online {


class DataSet{

public:

	TString _tag;
	const char* _data_folder_name;
	const char* _output_file_name;
	//TFile* _newFile;
	bool _is_data;
	TSelector* _reader;

	DataSet(){}
	DataSet(const char* data_folder_name, TString tag, bool is_data){
		_tag = tag;
		_data_folder_name = data_folder_name;
		_is_data = is_data;
	}

	void InitializeHistograms();

	void SetReader(TSelector* reader){
		_reader = reader;
	}


}; // CLASS DATA SET

class WJetsDataSet : public DataSet{
public:

	void GetHistograms();
	void InitializeHistograms();


    WJetsDataSet(const char* data_folder_name, TString tag, bool is_data){
		_tag = tag;
		_data_folder_name = data_folder_name;
		_is_data = is_data;
		
		GetHistograms();

		InitializeHistograms();
	}

	void RunReader(){
		TH1::SetDefaultSumw2(true);
    	std::vector<TString> fFiles;
    	std::cout << "Starting process..." << std::endl;
    	ProcessDirectory("", _data_folder_name, &fFiles, 1);

    //	_newFile->Write();
   // 	_newFile->Close();

	}

///////////////////////////////////////////////////////////////////////////////////////
//USER HISTOGRAMS TO BE WRITTEN TO FILE, RUN FILE
public:
	TH1* _fTwoProng_pt_signal;
    TH1* _fTwoProng_pt_sideband;
    TH1* _fTwoProng_pt_distr;

    TH1* _fTwoProng_eta_signal;
    TH1* _fTwoProng_eta_sideband;
    TH1* _fTwoProng_eta_distr;

    TH1* _fTwoProng_phi_signal;
    TH1* _fTwoProng_phi_sideband;
    TH1* _fTwoProng_phi_distr;

    TH1* _fTwoProng_MassPi0_signal;
    TH1* _fTwoProng_MassPi0_sideband;
    TH1* _fTwoProng_MassPi0_distr;

    TH1* _fW_pt_signal;
    TH1* _fW_pt_signal_noW;
    TH1* _fW_pt_signal_highptu;

    TH1* _fW_pt_sideband;
    TH1* _fW_pt_distr;

    TH1* _fW_eta_signal;
    TH1* _fW_eta_sideband;
    TH1* _fW_eta_distr;

    TH1* _fW_phi_signal;
    TH1* _fW_phi_sideband;
    TH1* _fW_phi_distr;

    TH1* _fW_mass_signal;
    TH1* _fW_mass_sideband;
    TH1* _fW_mass_distr;
    TH1* _fW_mass_distr_Imag_W;

    TH1* _fMuon_pt_signal;
    TH1* _fMuon_pt_sideband;
    TH1* _fMuon_pt_distr;

    TH1* _fMuon_eta_signal;
    TH1* _fMuon_eta_sideband;
    TH1* _fMuon_eta_distr;

    TH1* _fMuon_phi_signal;
    TH1* _fMuon_phi_sideband;
    TH1* _fMuon_phi_distr;

    TH1* _fTightMuon_pt_signal;
    TH1* _fTightMuon_pt_sideband;
    TH1* _fTightMuon_pt_distr;

    TH1* _fTightMuon_eta_signal;
    TH1* _fTightMuon_eta_sideband;
    TH1* _fTightMuon_eta_distr;

    TH1* _fTightMuon_phi_signal;
    TH1* _fTightMuon_phi_sideband;
    TH1* _fTightMuon_phi_distr;

    TH1* _fHT_distr;
    TH1* _fHT_signal;
    TH1* _fHT_sideband;

    TH1* _fHTbare_distr;
    TH1* _fHTbare_signal;
    TH1* _fHTbare_sideband;

    TH1* _MET_pt_distr;
    TH1* _MET_phi_distr; 

    TH1* _mT_distr;
    TH1* _mT_signal;
   
    TH1* _mu_w_deltaphi;
    TH1* _mu_met_deltaphi;

    TH2* _fTwoProng_MassPi0_PT;
    TH2* _fTwoProng_MassPi0_PT_sideband;

    TH2* _fTwoProng_MassPi0_PT_lowPTmuon;
    TH2* _fTwoProng_MassPi0_PT_sideband_lowPTmuon;

    TH2* _fTwoProng_MassPi0_PT_highPTmuon;
    TH2* _fTwoProng_MassPi0_PT_sideband_highPTmuon;

    TH2* _PT_Mass_Pi0_01;
    TH2* _PT_Mass_Pi0_02;
    TH2* _PT_Mass_Pi0_03;
    TH2* _PT_Mass_Pi0_04;
    TH2* _PT_Mass_Pi0_05;
    TH2* _PT_Mass_Pi0_06;
    TH2* _PT_Mass_Pi0_07;
    TH2* _PT_Mass_Pi0_08;
    TH2* _PT_Mass_Pi0_09;
    TH2* _PT_Mass_Pi0_10;

    TH2* _PT_Mass_Pi0_01_side;
    TH2* _PT_Mass_Pi0_02_side;
    TH2* _PT_Mass_Pi0_03_side;
    TH2* _PT_Mass_Pi0_04_side;
    TH2* _PT_Mass_Pi0_05_side;
    TH2* _PT_Mass_Pi0_06_side;
    TH2* _PT_Mass_Pi0_07_side;
    TH2* _PT_Mass_Pi0_08_side;
    TH2* _PT_Mass_Pi0_09_side;
    TH2* _PT_Mass_Pi0_10_side;

    TH1* _sample_weights;

    
}; //CLASS WJETSDATASET


class GJetsDataSet : public DataSet{
	public:

	void GetHistograms();
	void InitializeHistograms();


    GJetsDataSet(const char* data_folder_name, TString tag, bool is_data){
		_tag = tag;
		_data_folder_name = data_folder_name;
		_is_data = is_data;
		
		GetHistograms();

		InitializeHistograms();
	}

	void RunReader(){
		TH1::SetDefaultSumw2(true);
    	std::vector<TString> fFiles;
    	std::cout << "Starting process..." << std::endl;
    	ProcessDirectory("", _data_folder_name, &fFiles, 1);

    //	_newFile->Write();
   // 	_newFile->Close();

	}

public:
   TH1* _fTwoProng_MassPi0_pt40_60_test;
   TH1* _fTwoProng_MassPi0_pt40_60;
   TH1* _fTwoProng_pt_signal;
   TH1* _fTwoProng_pt_sideband;
   TH1* _fTwoProng_pt_distr;

   TH1* _fTwoProng_eta_signal;
   TH1* _fTwoProng_eta_sideband;
   TH1* _fTwoProng_eta_distr;

   TH1* _fTwoProng_phi_signal;
   TH1* _fTwoProng_phi_sideband;
   TH1* _fTwoProng_phi_distr;

   TH1* _fTwoProng_MassPi0_signal;
   TH1* _fTwoProng_MassPi0_sideband;
   TH1* _fTwoProng_MassPi0_distr;

   TH1* _fPhi_pt_signal;
   TH1* _fPhi_pt_sideband;
   TH1* _fPhi_pt_distr;

   TH1* _fPhi_eta_signal;
   TH1* _fPhi_eta_sideband;
   TH1* _fPhi_eta_distr;

   TH1* _fPhi_phi_signal;
   TH1* _fPhi_phi_sideband;
   TH1* _fPhi_phi_distr;

   TH1* _fPhi_mass_signal;
   TH1* _fPhi_mass_sideband;
   TH1* _fPhi_mass_distr;

   TH1* _fIDPhoton_pt_signal;
   TH1* _fIDPhoton_pt_sideband;
   TH1* _fIDPhoton_pt_distr;

   TH1* _fIDPhoton_eta_signal;
   TH1* _fIDPhoton_eta_sideband;
   TH1* _fIDPhoton_eta_distr;

   TH1* _fIDPhoton_phi_signal;
   TH1* _fIDPhoton_phi_sideband;
   TH1* _fIDPhoton_phi_distr;

   TH1* _fHT_distr;
   TH1* _fHT_signal;
   TH1* _fHT_sideband;

   TH1* _fHTbare_distr;
   TH1* _fHTbare_signal;
   TH1* _fHTbare_sideband;

   TH1* _MET_pt_distr;
   TH1* _MET_phi_distr;

   TH2* _fTwoProng_MassPi0_PT;
   TH2* _fTwoProng_MassPi0_PT_sideband;

   TH2* _fPhiMass_PT_signal;
   TH2* _fPhiMass_PT_sideband;

   bool _is_data;

   TH2* _PT_Mass_Pi0_01;
   TH2* _PT_Mass_Pi0_02;
   TH2* _PT_Mass_Pi0_03;
   TH2* _PT_Mass_Pi0_04;
   TH2* _PT_Mass_Pi0_05;
   TH2* _PT_Mass_Pi0_06;
   TH2* _PT_Mass_Pi0_07;
   TH2* _PT_Mass_Pi0_08;
   TH2* _PT_Mass_Pi0_09;
   TH2* _PT_Mass_Pi0_10;

   TH2* _PT_Mass_Pi0_01_side;
   TH2* _PT_Mass_Pi0_02_side;
   TH2* _PT_Mass_Pi0_03_side;
   TH2* _PT_Mass_Pi0_04_side;
   TH2* _PT_Mass_Pi0_05_side;
   TH2* _PT_Mass_Pi0_06_side;
   TH2* _PT_Mass_Pi0_07_side;
   TH2* _PT_Mass_Pi0_08_side;
   TH2* _PT_Mass_Pi0_09_side;
   TH2* _PT_Mass_Pi0_10_side;

   TH1* _sample_weights;
   int _N_SAMPLES_PASSED;


}; // CLASS GJetsDataSet


}; //NAMESPACE ONLINE
#endif