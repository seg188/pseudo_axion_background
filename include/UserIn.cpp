
#ifndef USER_IN
#define USER_IN


namespace user_input{

	const double MT_CUT = 50.0;
	const double MUON_PT_CUT = 50.0;
	const double W_PT_CUT = 50.0;
	const double GAMMA_TWOPRONG_DPHI_CUT = 0.10;
	const std::vector<int> PT_BOUNDS1{20, 40, 60, 80, 120, 160, 240, 300, 460, 700}; //IF YOU CHANGE THIS, MAKE SURE TO CHANGE THE CORRESPONDING MEM BER IN ::histogram


	//GJETS

	const double PHOTON_PT_CUT = 150.0;
	
}; //NAMESPACE user_inputs

namespace global_tags{
	TString ReaderTag = "";
	TString tree_tag = "";
};

namespace histogram{

	const double PT_BOUNDS[11] = {20.0, 40.0, 60.0, 80.0, 120.0, 160.0, 240.0, 300.0, 460.0, 700.0, 9999999.9};

	int PT_NUMBINS(400);
	int W_MASS_NUMBINS(200);
	int HT_NUMBINS(100);

	int MT_NUMBINS(400);
	int MT_MIN(0);
	int MT_MAX(200);

	double PT_MAX(2000);
	double PT_MIN(0);

	double W_MASS_MIN(50);
	double W_MASS_MAX(450);

	double HT_MIN(0);
	double HT_MAX(4000);

	int ETA_NUMBINS(200);
	double ETA_MINBIN(-5);
	double ETA_MAXBIN(5);

	int PHI_NUMBINS(200);
	double PHI_MINBIN(-4);
	double PHI_MAXBIN(4);

	int DPHI_NUMBINS(250);
	double DPHI_MINBIN(6.3);
	double DPHI_MAXBIN(6.3);

	int TP_MASS_NUMBINS(30);
	double TP_MASS_MINBIN(0);
	double TP_MASS_MAXBIN(5.0);

	const int PHI_MASS_NUMBINS = 10;
	const double PHI_MASS_MIN = 0.0; 
	const double PHI_MASS_MAX = 2000.0;

}; //NAMESPACE histogram

namespace constants{
	const double LUMI2016 = 36814;
}

#endif