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


GP1D_Fitter::GP1D_Fitter(std::vector<double> datax, std::vector<double> datay, std::vector<double> err = {}){
		_data.resize(datax.size());

		for (int x = 0; x < _data.size(); x++){
			_data[x] = {datax[x], datay[x]};
		}

		if (err.size() == datax.size()){
			for (auto x : err){
				_error.push_back(TMath::Log(x*x));
			}
			
		} else {
			_error.resize(datax.size());
			for (int x = 0; x < _data.size(); x++){
				_error[x] = TMath::Log(BETA);
			}
		}
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
		_error[db] = TMath::Log(TMath::Power(data->GetBinError(x), 2));
		db++;
	}
	std::cout << "...done" << std::endl;
}


void GP1D_Fitter::SetData(std::vector<std::vector<double>> data){
		_data = data;
}


void GP1D_Fitter::SetParameters(std::vector<double> params){
	_hyperparams = params;
	Update();
}

void GP1D_Fitter::Update(){

	SetCMatrix();
	Kernel::_params = _hyperparams;
	SetMean();
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
		_error[db] = TMath::Log(TMath::Power(data->GetBinError(x), 2));
		db++;
	}
	
	std::cout << "...done" << std::endl;
}
void GP1D_Fitter::SetKernel(double (*ker)(double, double, std::vector<double>), std::vector<double> params){
	_ker = ker;
	_hyperparams = params;
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

double GP1D_Fitter::KT_CNinv_K(double x){
	std::vector<double> CN_inv_K;
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

void GP1D_Fitter::SetCMatrix(){
	TMatrixD C_N(_data.size(), _data.size());
	_C_N.ResizeTo(C_N);
	_C_N_inv.ResizeTo(C_N);
	for (int i = 0; i<_data.size(); i++){
		for (int j = 0; j < _data.size(); j++){
			C_N[i][j] = (*(_ker_vector[i]))(_data[j][0]);
		}
	}

	for (int j = 0; j < _data.size(); j++){
			C_N[j][j] += _error[j] + 1.0/BETA;

	}

	_C_N = C_N;
	_C_N_inv = C_N;
	_C_N_inv.Invert();

}


void GP1D_Fitter::SetMean(){
	_target.resize(_data.size());
	for (int j = 0; j < _target.size(); j++) _target[j] = _data[j][1];
	_C_inv_t_prod = product(_C_N_inv, _target );
}

double GP1D_Fitter::mean(double x){
	double val = 0;
	for (int k = 0; k < _ker_vector.size(); k++){
		val += (*(_ker_vector[k]))(x)*_C_inv_t_prod[k];
	}
	return val;
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
     	_sigma_1_graph->SetPoint(i, x, mean(x) + TMath::Sqrt(error(x)));
     
     	_sigma_1_graph->SetPoint(NUM_AXIS_POINTS + i, x2, mean(x2) - TMath::Sqrt(error(x2)));
  	}
   	_sigma_1_graph->SetFillStyle(3005);
    _sigma_1_graph->SetFillColor(2);
 


}
void GP1D_Fitter::SetXAxis(){
	_X_AXIS.resize(NUM_AXIS_POINTS);
	double d_npoints = static_cast<double>(NUM_AXIS_POINTS);
	for (int x = 0; x < NUM_AXIS_POINTS; x++){
		_X_AXIS[x] = _data[0][0] + (_data[_data.size()-1][0] - _data[0][0]) * static_cast<double>(x) / d_npoints;
	}
}

void GP1D_Fitter::InitializeErrorFit(){
	std::vector<double> x_data;
	x_data.resize(_data.size());

	for (int j = 0; j < x_data.size(); j++) x_data[j] = _data[j][0];

	_error_fitter = new GP1D_Fitter(x_data, _error);
	//HijackErrorFunction();
	_error_fitter->SetKernel(_ker, _hyperparams);
	_error_fitter->BETA = BETA;
	_error_fitter->LearnParams();
	_error_fitter->Fit();

}

double GP1D_Fitter::error(double x){
	double val = _ker(x, x, _hyperparams);
	if (_error_fitter) val += TMath::Exp(_error_fitter->mean(x));
	else val += 1.0/BETA;
	val = val - KT_CNinv_K(x);
	return val;
}

void GP1D_Fitter::Fit(bool with_error){ 
	////BASIC FIT WITH NO HYERPARAM MINIMIZATION-----------------------
	SetXAxis();
	InitializeErrorFit();
	SetCMatrix();
	SetMean();
	Update();
	InitializeGraphs();
}

void GP1D_Fitter::Fit(){ 
	////BASIC FIT WITH NO HYERPARAM MINIMIZATION-----------------------
	SetXAxis();
	//InitializeErrorFit();
	SetCMatrix();
	SetMean();
	Update();
	InitializeGraphs();
}

void GP1D_Fitter::LearnParams(double (*function)(std::vector<double>) = nullptr){
//	InitializeErrorFit();
	std::cout << "Learning hyperparameters..." << std::endl;
	fit::Minimizer::fitter = this;
	_param_learner.Set(_hyperparams, function);
	_param_learner.Minimize();
	_hyperparams = _param_learner._min_params;
	std::cout << "Found {";
	for (double val : _hyperparams) std::cout << val << ", ";
	std::cout << "}" << std::endl;
	Update();
	std::cout << "...done" << std::endl;
}


}; //namespace fit



#endif