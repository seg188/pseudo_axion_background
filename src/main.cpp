#include "DataSet.hh"
#include "DataSet.cpp"
#include "analysis/GenerateAnalysis.cpp"
//#include "fit/fit.cpp"
//#include "fit/GPR/GP1D_Fitter.hh"
//#include "fit/GPR/GP1D_Fitter.cpp"


using namespace online;

const char* wj_folder_name = "/cms/seg188/eos/twoprong/trees/wjets_nov18";
TString tag = "wjets2016_mc";
bool is_data = false;

const char* gj_folder_name = "/cms/chiarito/eos/twoprong/trees/no_filter/gjets2016/Mar8";
TString gj_tag = "gjets2016_mc";
bool gj_is_data = false;

const char* sm_folder_name = "/cms/seg188/eos/twoprong/trees/SingleMuon_oct04";
TString sm_tag = "wjets_data";
bool sm_is_data = true; 

const char* sp_folder_name = "/cms/chiarito/eos/twoprong/trees/no_filter/photon2016/Mar31_fullAsym/SinglePhoton";
TString sp_tag = "sp_data";
bool sp_is_data = true;


void ReadGJetsData(const char* folder_name, TString tag,bool is_data ){
	TFile* f = new TFile(tag + ".root", "RECREATE");
	GJetsDataSet gj_data(folder_name, tag, is_data);
	gj_data.RunReader();

	f->cd();

	auto Omega_Mass_Distributions = analysis::GetMassDistributions(gj_data);

	for (auto vec : Omega_Mass_Distributions){
		for (auto plot : vec) plot->Write();
	}

	auto Phi_Mass_PT = gj_data._fPhiMass_PT_signal;
	Phi_Mass_PT->Write();

	auto GJets_Ratios = analysis::GenerateAnalysis(gj_data);

	for (auto plot : GJets_Ratios) plot->Write();
}

void ReadSPData(const char* folder_name, TString tag,bool is_data ){
	TFile* f = new TFile(tag + ".root", "RECREATE");
	GJetsDataSet gj_data(folder_name, tag, is_data);
	gj_data.RunReader();

	f->cd();

	auto Omega_Mass_Distributions = analysis::GetMassDistributions(gj_data);

	for (auto plot : Omega_Mass_Distributions[1]){
		plot->Write();
	}
}

void ReadWJetsData(const char* folder_name, TString tag,bool is_data ){
	TFile* f = new TFile(tag + ".root", "RECREATE");
	WJetsDataSet wj_data(folder_name, tag, is_data);
	wj_data.RunReader();

	auto WJets_Ratios = analysis::GenerateAnalysis(wj_data);
	f->cd();

	for (auto plot : WJets_Ratios) plot->Write();
}

void ReadSMData(const char* folder_name, TString tag,bool is_data ){
	TFile* f = new TFile(tag + ".root", "RECREATE");
	WJetsDataSet wj_data(folder_name, tag, is_data);
	wj_data.RunReader();

	auto WJets_Ratios = analysis::GenerateAnalysis(wj_data);
	f->cd();

	for (auto plot : WJets_Ratios) plot->Write();
}

int main(){

	ReadGJetsData(gj_folder_name, gj_tag, gj_is_data);

//	ReadWJetsData(wj_folder_name, tag, is_data);

//	ReadWJetsData(sm_folder_name, sm_tag, sm_is_data);

//	ReadSPData(sp_folder_name, sp_tag, sp_is_data);


	


	//fit::FitRatios(gj_data);
}