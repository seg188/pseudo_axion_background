#include "GP1D_Fitter.hh"
#include "GP1D_Fitter.cpp"
#include "../../analysis/GenerateAnalysis.cpp"


double ker(double x1, double x2, std::vector<double> parameters){
	double val = TMath::Exp(-1.000 * parameters[0] * TMath::Power(x1-x2, 2.000));
	return val;
}

int NPARS = 1;

double cin_param(int num){
	 while (true) // Loop until user enters a valid input
    {
        std::cout << "param " << num << ": ";
        double x;
        std::cin >> x;
 
        if (std::cin.fail()) // has a previous extraction failed?
        {
            // yep, so let's handle the failure
            std::cin.clear(); // put us back in 'normal' operation mode
            std::cin.ignore(32767,'\n'); // and remove the bad input
        }
        else // else our extraction succeeded
        {
            std::cin.ignore(32767, '\n'); // clear (up to 32767) characters out of the buffer until a '\n' character is removed
            return x; // so return the value we extracted
        }
    }
}

bool approve_fit(){
	while (true) // Loop until user enters a valid input
    {
        std::cout << "good?: (1/0): " ;
        int x;
        std::cin >> x;
 
        if (std::cin.fail()) // has a previous extraction failed?
        {
            // yep, so let's handle the failure
            std::cin.clear(); // put us back in 'normal' operation mode
            std::cin.ignore(32767,'\n'); // and remove the bad input
        }
        else // else our extraction succeeded
        {
            std::cin.ignore(32767, '\n'); // clear (up to 32767) characters out of the buffer until a '\n' character is removed
            if (x == 1 ) return true;
            if (x == 0) return false;
            
            std::cout << "type 1 (yes) or 0 (no)" << std::endl;
        }
    }
}

std::vector<fit::GP1D_Fitter*> fit_vector(std::vector<TH1*> plot_vector, bool wait_approve = true){
	std::vector<fit::GP1D_Fitter*> fits;
	TCanvas* fit_canvas = new TCanvas("fit_canvas", "fits");
	for (int plt = 0; plt < plot_vector.size(); plt++){
		TString name_tag = TString("../../../out/plots/") + plot_vector[plt]->GetName();

		bool approved = false;
		fits.push_back(new fit::GP1D_Fitter());
		std::vector<double> params1;
		for (int n = 0; n < NPARS; n++) params1.push_back(1.0);

		while (!approved){
			delete fits[plt];
			fits[plt] = new fit::GP1D_Fitter();
			fits[plt]->SetData(plot_vector[plt]);
			fits[plt]->SetKernel(ker, params1);
			int try_n = 0;
			std::vector<double> params;

			if (wait_approve){
				for (int n = 0; n < NPARS; n++) params.push_back(cin_param(n));
					
			} else {

				for (int n = 0; n < NPARS; n++) params.push_back(1.00);	
			}
		
			fits[plt]->_hyperparams = params;
			fits[plt]->LearnParams();

			if (try_n == 0)fits[plt]->Fit(true);
			else {
				fits[plt]->Update();
				fits[plt]->UpdateGraphs();
			}

			try_n++;
			
			if (wait_approve){
				fit_canvas->cd();
				plot_vector[plt]->Draw();
				fits[plt]->_sigma_1_graph->Draw("SAME f");
				fits[plt]->_mean_graph->Draw("SAME");
				fit_canvas->Print(name_tag + TString(".pdf"), ".pdf");
				approved = approve_fit();
			} else {
				approved = true;
			}
			
		}//while

	} //for plt

	return fits;

}




void test_gpr(){
	TRandom* ranGen = new TRandom(time(NULL));
/*
	int n_bins = 30;
	double x_min = 0.0;
	double x_max = 5.0;
	
	TH1* data1 = new TH1D("test", "test", n_bins, x_min, x_max);
	TH1* errs = new TH1D("errs", "test", n_bins, x_min, x_max);

	for (int k = 1; k <= data1->GetNbinsX() - 3; k++){
		double x = data1->GetBinCenter(k);
		data1->SetBinContent(k, TMath::Gaus(x, 3)*TMath::Exp(3-x) + ranGen->Gaus(0.00, 0.05));
		data1->SetBinError( k, TMath::Abs(ranGen->Gaus(0.0, 0.10)) );
		errs->SetBinContent(k, data1->GetBinError(k));
	
	} 

*/

	TFile* file = TFile::Open("~/hex/sp_data.root");

	std::vector<TH1*> plots;
	for (int k = 0; k < 9; k++){
		std::ostringstream strs;
		strs << k;
		plots.push_back((TH1*)file->Get(TString("sp_data_MassPi0_PT_") + TString(strs.str()) + TString("_side") ));
	}
	
   auto fits = fit_vector(plots);

	auto data = plots[1];
	fit::GP1D_Fitter fitter;
	fitter.SetData(data);
	fitter.SetKernel(ker, {1.0});
	fitter.LearnParams();
	fitter.Fit(true);
	TCanvas* canvas = new TCanvas("canv");
	data->Draw();
	fitter._sigma_1_graph->Draw("SAME f");
	fitter._mean_graph->Draw("SAME");

	//errs->Draw();
	//fitter._error_fitter->_mean_graph->Draw("SAME");


/*


	
	
	//1.0});//1.0, 1.0}); //FIRST PARAMETER IS ALWAYS BETA
	
	");
	

*/
	//
//	std::cout << "min: " << min._current_params[0] << std::endl;

}