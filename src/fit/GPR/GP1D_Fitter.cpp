#include "GP1D_Fitter.hh"
#include "minimizer.cpp"
#ifndef GP_FITTER_C
#define GP_FITTER_C

namespace fit{
template <typename vector, typename matrix>
vector product(matrix first, vector second){
	vector result;
	result.resize(first.GetNcols());
	for (int i = 0; i < second.size(); i++){
		double sum = 0;
		for (int j = 0; j < first.GetNrows(); j++){
			sum += first[i][j] * second[j];
		}
		result[i] = sum;	
	}

	return result;
}

void checkErrors(GP1D_Fitter* fit){
	//TCanvas* canv = new TCanvas("error_check");
	//TH1* error_plot = new TH1D("errors", "errors", fit->_data.size(), fit->_X_AXIS[0], fit->_X_AXIS[fit->_X_AXIS.size() - 1]);
	//for (int x = 0; x < fit->_error.size(); x++){
	//	error_plot->SetBinContent(x+1, fit->_error[x]);
	//}

	//error_plot->Draw();
	//fit->_error_fitter->_mean_graph->Draw("SAME");

}


GP1D_Fitter::GP1D_Fitter(TH1* data){
	std::cout << "setting data..." << std::endl;
	int num_zeros = 0;
	for (int x =1; x <= data->GetNbinsX(); x++){
		if (data->GetBinError(x) == 0) num_zeros++;
	}
	_data.resize(data->GetNbinsX() - num_zeros );
	_error.resize(data->GetNbinsX() - num_zeros);
	int db = 0;
	
	for (int x = 1; x <= data->GetNbinsX(); x++){
		if (data->GetBinError(x) == 0) continue;
		std::vector<double> point = {data->GetBinCenter(x), data->GetBinContent(x)};
		_data[db] = point;
		_error[db] = TMath::Power(data->GetBinError(x), 2);
		db++;
	}
	std::cout << "...done" << std::endl;
	_target.resize(_data.size());
	for (int j = 0; j < _target.size(); j++) _target[j] = _data[j][1];
}

GP1D_Fitter::GP1D_Fitter(std::vector<double> datax, std::vector<double> datay, std::vector<double> err = {}){
		_data.resize(datax.size());

		for (int x = 0; x < _data.size(); x++){
			_data[x] = {datax[x], datay[x]};
		}

		if (err.size() == datax.size()){
			_error = err;
		} else {
			_error.resize(datax.size());
			for (int x = 0; x < _data.size(); x++){
				_error[x] = 1.00/BETA;
			}
		}
		_target.resize(_data.size());
		for (int j = 0; j < _target.size(); j++) _target[j] = _data[j][1];
}

void GP1D_Fitter::SetData(std::vector<double> datax, std::vector<double> datay, std::vector<double> err = {}){
		_data.resize(datax.size());

		for (int x = 0; x < _data.size(); x++){
			_data[x] = {datax[x], datay[x]};
		}

		if (err.size() == datax.size()){
			_error = err;
		} else {
			_error.resize(datax.size());
			for (int x = 0; x < _data.size(); x++){
				_error[x] = 1.00/BETA;
			}
		}
	_target.resize(_data.size());
	for (int j = 0; j < _target.size(); j++) _target[j] = _data[j][1];
}
void GP1D_Fitter::SetData(std::vector<std::vector<double>> data){
		_data = data;
		_target.resize(_data.size());
		for (int j = 0; j < _target.size(); j++) _target[j] = _data[j][1];
}


void GP1D_Fitter::SetParameters(std::vector<double> params){
	_hyperparams = params;
	Update();
}

void GP1D_Fitter::Update(){
	if (_error_fitter)_error_fitter->Update();
	Kernel::_params = _hyperparams;
	SetCMatrix();
	SetMean();
	GetVariance();
}

void GP1D_Fitter::SetData(TH1* data){
	std::cout << "setting data..." << std::endl;
	int num_zeros = 0;
	for (int x =1; x <= data->GetNbinsX(); x++){
		if (data->GetBinError(x) == 0) num_zeros++;
	}
	_data.resize(data->GetNbinsX() - num_zeros );
	_error.resize(data->GetNbinsX() - num_zeros);
	int db = 0;
	
	for (int x = 1; x <= data->GetNbinsX(); x++){
		if (data->GetBinError(x) == 0) continue;
		std::vector<double> point = {data->GetBinCenter(x), data->GetBinContent(x)};
		_data[db] = point;
		_error[db] = TMath::Power(data->GetBinError(x), 2);
		db++;
	}
	_target.resize(_data.size());
	for (int j = 0; j < _target.size(); j++) _target[j] = _data[j][1];
	std::cout << "...done" << std::endl;
}
void GP1D_Fitter::GetBeta(bool with_error){

	double variance = 0;
	for (int x = 0; x < _data.size(); x++){
		variance += TMath::Power(mean(_data[x][0]) - _data[x][1], 2.0);
	}

	BETA = 1.0/variance;

	std::cout << "BETA: " << BETA << std::endl;


}

void GP1D_Fitter::SetKernel(double (*ker)(double, double, std::vector<double>), std::vector<double> params){
	_ker = ker;
	_hyperparams = params;
	NPARS = _hyperparams.size();
	Kernel::_params = _hyperparams;
	SetKerVector();
}


void GP1D_Fitter::SetKerVector(){
	std::cout << "setting kernel vectors..." << std::endl;
	for (int x = 0; x < _data.size(); x++){
		double point = _data[x][0];
		Kernel* ker = new Kernel(_ker, point);
		_ker_vector.push_back(ker);
	}

	std::cout << "...done" << std::endl;
}

void GP1D_Fitter::SetCMatrix(){
	Kernel::_params = _hyperparams;
	TMatrixD C_N(_data.size(), _data.size());
	_C_N.ResizeTo(C_N);
	_C_N_inv.ResizeTo(C_N);
	for (int i = 0; i<_data.size(); i++){
		for (int j = 0; j < _data.size(); j++){
			C_N[i][j] = (*(_ker_vector[j]))(_data[i][0]);
		}
	}

	for (int j = 0; j < _data.size(); j++){
			C_N[j][j] +=  1.0/BETA;
		

	}

	_C_N = C_N;
	_C_N_inv = C_N;
	_C_N_inv.Invert();

}
void GP1D_Fitter::clear(){
	_data = {{}};
	_error = {};
	_C_N = TMatrixD();
	_C_N_inv = TMatrixD();
	//for (auto ptr : _ker_vector) delete ptr;


}

void GP1D_Fitter::SetMean(){
	_C_inv_t_prod = product(_C_N_inv, _target );
}

double GP1D_Fitter::mean(double x){
	Kernel::_params = _hyperparams;
	double val = 0;
	for (int k = 0; k < _ker_vector.size(); k++){
		val += (*(_ker_vector[k]))(x)*_C_inv_t_prod[k];
	}
	return val;
}

void GP1D_Fitter::UpdateGraphs(){
	delete _mean_graph;
	delete _sigma_1_graph;
	InitializeGraphs();
}

void GP1D_Fitter::InitializeGraphs(){
	_mean_graph = new TGraph(NUM_AXIS_POINTS);
	for (int x = 0; x < NUM_AXIS_POINTS; x++){
		_mean_graph->SetPoint(x, _X_AXIS[x], mean(_X_AXIS[x]) );
	}
	_mean_graph->SetLineColor(kRed);

	_sigma_1_graph = new TGraph(NUM_AXIS_POINTS);

    for (int i = 0; i < NUM_AXIS_POINTS; i++) {
    	double x = _X_AXIS[i];
   		double x2 = _X_AXIS[NUM_AXIS_POINTS-i-1];
     	_sigma_1_graph->SetPoint(i, x, mean(x) + error(x));
     
     	_sigma_1_graph->SetPoint(NUM_AXIS_POINTS + i, x2, mean(x2) - error(x2));
  	}
   	_sigma_1_graph->SetFillStyle(3002);
    _sigma_1_graph->SetFillColor(2);
 


}
void GP1D_Fitter::SetXAxis(){
	_X_AXIS.resize(NUM_AXIS_POINTS);
	double d_npoints = static_cast<double>(NUM_AXIS_POINTS);
	for (int x = 0; x < NUM_AXIS_POINTS; x++){
		_X_AXIS[x] = _data[0][0] + (_data[_data.size()-1][0] - _data[0][0]) * static_cast<double>(x) / d_npoints;
	}
}
void GP1D_Fitter::Fit(bool with_error){ 
	////BASIC FIT WITH NO HYERPARAM MINIMIZATION-----------------------
	SetXAxis();
	SetCMatrix();
	SetMean();
	//GetBeta(true);
	InitializeErrorFit();
	Update();
	InitializeGraphs();
}

void GP1D_Fitter::Fit(){ 
	////BASIC FIT WITH NO HYERPARAM MINIMIZATION-----------------------
	SetXAxis();
	SetCMatrix();
	SetMean();
	//GetBeta(true);
	Update();
	InitializeGraphs();
}

void GP1D_Fitter::CheckKer(int num){
	

    TGraph* graph = new TGraph(NUM_AXIS_POINTS);
	for (int x = 0; x < NUM_AXIS_POINTS; x++){
		graph->SetPoint(x, _X_AXIS[x], (*(_ker_vector[num]))(_X_AXIS[x]) );
		std::cout  << (*(_ker_vector[num]))(_X_AXIS[x]) << std::endl;
	}


	TGraph* graph2 = new TGraph(NUM_AXIS_POINTS);
	for (int x = 0; x < NUM_AXIS_POINTS; x++){
		graph2->SetPoint(x, _X_AXIS[x], _ker(_data[num][0], _X_AXIS[x], _hyperparams ) );
	}
	graph2->SetLineColor(kBlue);
	TCanvas* canvas2 = new TCanvas("kercheck2", "kercheck");
	graph2->Draw();

	TCanvas* canvas = new TCanvas("kercheck", "kercheck");
	graph->SetLineColor(kGreen);
	graph->Draw();
}

double GP1D_Fitter::KT_CNinv_K(double x){
	std::vector<double> CN_inv_K;
	Kernel::_params = _hyperparams;
	CN_inv_K.resize(_ker_vector.size());
	for (int j = 0; j < _ker_vector.size(); j++){
		CN_inv_K[j] = 0;
		for (int k = 0; k < _ker_vector.size(); k++){
			CN_inv_K[j] += (*(_ker_vector[k]))(x) * _C_N_inv[j][k];
		}
	}

	double val = 0;

	for (int j = 0; j < _ker_vector.size(); j++){
		val += (*(_ker_vector[j]))(x) * CN_inv_K[j];
	}

	return val;
}

double GP1D_Fitter::_error_fit(double x){
	if (_error_fitter) {
		//std::cout << _error_fitter->mean(x) << std::endl;
		return TMath::Abs(_error_fitter->mean(x));
		//return 0.0;
	}

	else return 0.0;
}

/*
void GP1D_Fitter::HijackErrorFunction(){
	_data_store = _data;
	_error_store = _error;

	std::vector<double> xdata;
	xdata.resize(_data.size());

	for (int k = 0; k < xdata.size(); k++){
		xdata[k] = _data[k][0];
	}

	SetData(xdata, _error);
	_param_learner.Minimize();

	_data = _data_store;
	_error = _error_store;

}
*/

void GP1D_Fitter::InitializeErrorFit(){
	std::vector<double> x_data;
	x_data.resize(_data.size());

	for (int j = 0; j < x_data.size(); j++) x_data[j] = _data[j][0];

	_error_fitter = new GP1D_Fitter(x_data, _error);
	_error_fitter->SetKernel(_ker, _hyperparams);
	_error_fitter->BETA = BETA;
	//HijackErrorFunction();
	/*
	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "");
	min->SetMaxFunctionCalls(1000); // for Minuit/Minuit2
  	min->SetMaxIterations(1000);  // for GSL
  	min->SetTolerance(0.00001);
  	min->SetPrintLevel(0);
    
  	std::cout << "learning beta...." << std::endl;
  	fit::Minimizer::fitter = _error_fitter;
  	ROOT::Math::Functor log_lbeta(&fit::Minimizer::FitBeta, 1);
  	double beta_step = 1.0;
  	min->SetFunction(log_lbeta);
  	min->SetVariable(0, "beta", BETA, beta_step);
  	min->Minimize();

  	BETA = min->X()[0];

  	std::cout << "LEARNED ERROR BETA: " << BETA << std::endl;
  	*/
	
	_error_fitter->Fit();
	checkErrors(this);

}
void GP1D_Fitter::GetVariance(){
	double variance = 0;
	for (int x = 0; x < _data.size(); x++){
		variance += TMath::Power(_data[x][1] - mean(_data[x][0]), 2.0);
	}
	_variance = variance;

}

double GP1D_Fitter::error(double x){
	double val = _error_fit(x) + _variance ;// + _ker(x, x, _hyperparams);// 1.00/(BETA) ;
//	val = val - KT_CNinv_K(x);


	return TMath::Sqrt(val);
}
/*

void GP1D_Fitter::LearnErrors(double (*function)(std::vector<double>)){
	_param_learner.Set(function, _hyperparams);
	InitializeErrorFit();

}
*/

void GP1D_Fitter::LearnParams(double (*function)(std::vector<double>) = nullptr){
//	InitializeErrorFit();
	std::cout << "Learning hyperparameters..." << std::endl;
	fit::Minimizer::fitter = this;

	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
	min->SetMaxFunctionCalls(10000); // for Minuit/Minuit2
  	min->SetMaxIterations(10000);  // for GSL
  	min->SetTolerance(0.0001);
  	min->SetPrintLevel(1);
  	
	std::cout << "learning beta...." << std::endl;
  	ROOT::Math::Functor log_lbeta(&fit::Minimizer::FitBeta, 1);
  	double beta_step = 1.0;
  	min->SetFunction(log_lbeta);
  	if (isnan(BETA)) BETA = 100.000;
  	min->SetVariable(0, "beta", BETA, beta_step);
  	min->Minimize();


  	BETA = min->X()[0]; 

  	std::cout << "LEARNED BETA: " << BETA << std::endl;
  	


  	ROOT::Math::Functor log_liklihood(&fit::Minimizer::FitChi1, _hyperparams.size());

  	 std::vector<double> step;
  	 for (auto x : _hyperparams) step.push_back(0.01);

  	min->SetFunction(log_liklihood);
  	int k =0;
  	std::ostringstream strs;
  	for (auto x : _hyperparams){
  		strs<<"a";
  		min->SetVariable(k,strs.str(),_hyperparams[k], step[k]);
  		k++;
  	}

  	min->Minimize();
  	const double *xs = min->X();
  	k=0;
  	std::cout << "learned: " << std::endl << "{";
  	for (auto x : _hyperparams){
  		std::cout << xs[k] << ", ";
  		_hyperparams[k] = xs[k];
  		k++;
  	}
  	std::cout << "}" << std::endl;

/*
  	std::cout << "learning beta...." << std::endl;
  	
  	min->SetFunction(log_lbeta);
  	min->SetVariable(0, "beta", BETA, beta_step);
  	min->Minimize();

  	BETA = min->X()[0];
	*/
  	

  	std::cout << "BETA: " << BETA << std::endl;
 
 	std::cout << "learned: " << std::endl << "{";
  	for (auto x : _hyperparams){
  		std::cout << x << ", ";
  	}
  	std::cout << "}" << std::endl;



  	Update();
  
/*
	_param_learner.Set(_hyperparams, function);
	_param_learner.Minimize();
	_hyperparams = _param_learner._min_params;
	std::cout << "Found {";
	for (double val : _hyperparams) std::cout << val << ", ";
	std::cout << "}" << std::endl;
	Update();
	std::cout << "...done" << std::endl;
	*/
}


}; //NAMESPACE FIT

#endif