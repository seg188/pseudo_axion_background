#ifndef GP_FITTER_H
#define GP_FITTER_H

namespace fit{

class GP1D_Fitter;

class Kernel
{
 	public:
 	double (*_function)(double, double, std::vector<double>);
 	double _fixed_point;
 	
 	static std::vector<double> _params;
 	
  	Kernel( double (*function)(double, double, std::vector<double>) , double point){
  		_fixed_point = point;
  		_function = function;
  	}


  	double operator ()(double x){
  		return _function(_fixed_point, x, _params);
  	}
//int (*ptr)( int ) = la; // E R R O R. See the workaround:  
	
};

std::vector<double> Kernel::_params = {};    

class Minimizer{
public:
	int NUM_GRID_POINTS = 200;
	double FULL_RANGE = 3.0;
	double BETA_RANGE = 5.0;
	double GR_STEP_SIZE = 0.1;
	double LEARNING_RATE = 0.001;
	double THRESHOLD = 0.01;
	int NUM_GRID_CHECKS = 20;
	double _global_min = 99999.0;
	void SetThreshold(double t){THRESHOLD = t;}
	void SetNumGridChecks(int n){NUM_GRID_CHECKS = n;}
	double (*_function)(std::vector<double>);
	std::vector<double> _current_params;
	std::vector<double> _min_params;
	static GP1D_Fitter* fitter;
	static double FitChi2(std::vector<double>);
	static double FitChi1(const double *);
	static double FitBeta(const double *);


	Minimizer(std::vector<double> initial, double (*function)(std::vector<double>)=nullptr){
		if (function) _function = function;
		else _function = FitChi2;
		_current_params = initial;
		_min_params = _current_params;
	}

	void Set(std::vector<double> initial, double (*function)(std::vector<double>)=nullptr){
		if (function) _function = function;
		else _function = FitChi2;
		_current_params = initial;
		_min_params = _current_params;
	}


	Minimizer(){}

	double MinimizeParameter(int par){
		_current_params = _min_params;
		double par_val = _current_params[par];
		double center = _current_params[par];
		double step_size = FULL_RANGE / static_cast<double>(NUM_GRID_POINTS);
		if (par == 0) step_size = BETA_RANGE / static_cast<double>(NUM_GRID_POINTS);

		for (int k = 0; k < NUM_GRID_POINTS; k++){
			_current_params[par] = center - 0.50*FULL_RANGE + step_size*static_cast<double>(k);
			double val = _function(_current_params);
			if (val < _global_min){
	
				_global_min = val;
				_min_params = _current_params;
			}
		}

		return par_val;
	}

	bool CheckGrad(std::vector<double> vector){
		bool all_minimized = true;
		int j = 0;
		for (double val : vector){ 
			if (TMath::Abs(val) > THRESHOLD) all_minimized = false;
			//if (TMath::Abs(val) > 999888899.) MinimizeParameter(j);
			j++;
		
		}

		return all_minimized;
	}

	void GradientDecent(){
		std::vector<double> grad = gradient();
		std::cout << "{";
		//for (int j = 0; j < _current_params.size(); j++){
		//		std::cout <<  grad[j] << ", ";
		//	}
		//	std::cout << "}" << std::endl;

		while (!CheckGrad(grad)){
			std::cout << "{";
			double momentum = 0;
			for (double val : grad) momentum += val*val;
			for (int j = 0; j < _current_params.size(); j++){
				std::cout <<  grad[j] << ", ";
				_current_params[j] += LEARNING_RATE * grad[j];// + LEARNING_RATE * grad[j] * momentum;
			}
			std::cout << "}" << std::endl;
			grad = gradient();
		}

	//	_min_params = _current_params;
	}

	void GridSearchMin(int j){
		//if (j > _current_params.size()) return;
		for (int k = j; k < _current_params.size(); k++){
			MinimizeParameter(k);
			//GridSearchMin(j+1);
		}
		

	}

	void Minimize(){

		for (int j = 0; j < NUM_GRID_CHECKS; j++) GridSearchMin(0);
		//GradientDecent();

		_current_params = _min_params;

	}

	std::vector<double> gradient(){
		std::vector<double> parameters;
		std::vector<double> grad;
		grad.resize(_current_params.size());

		for (int j = 0; j < _current_params.size(); j++){
			parameters = _current_params;
			parameters[j] += GR_STEP_SIZE;
			double diff = _function(parameters) - _function(_current_params);
			grad[j] = diff / GR_STEP_SIZE;
		}

		return grad;
	}


}; //class minimizer

class GP1D_Fitter{

	public:
	int NUM_AXIS_POINTS = 200;
	static double BETA;
	TGraph* _mean_graph;
	TGraph* _sigma_1_graph;
	TGraph* _sigma_2_graph;
	void SetKernel(double (*ker)(double, double, std::vector<double> parameters), std::vector<double> hyperparams);
	void Fit();
	void Fit(bool);
	GP1D_Fitter(TH1*);
	GP1D_Fitter(std::vector<double>, std::vector<double>, std::vector<double>);
	double mean(double);
	double error(double);
	void CheckKer(int);
	std::vector<double> _hyperparams;
	void SetData(TH1* data);
	void SetData(std::vector<double>, std::vector<double>, std::vector<double>);
	void SetData(std::vector<std::vector<double>>);
	GP1D_Fitter(){}
	void SetParameters(std::vector<double> params);
	void Update();
	void LearnParams(double (*function)(std::vector<double>)=nullptr);
	fit::Minimizer _param_learner;
	//bool _mini_set = false;
	std::vector<std::vector<double>> _data_store;
	std::vector<double> _error_store;
	void GetBeta(bool with_error);
	void clear();
	void UpdateGraphs();
	static int NPARS;
	void GetVariance();
	double _variance = 0;
//	void HijackErrorFunction();
//	void LearnErrors(double (*function)(std::vector<double>));

	public:
	std::vector<std::vector<double>> _data;
	std::vector<double> _error;
	TMatrixD _C_N;
	TMatrixD _C_N_inv;
	std::vector<double> _C_inv_t_prod;
	double (*_ker)(double, double, std::vector<double>);
	double _x_min;
	double _x_max;
	void SetKerVector();
	void SetCMatrix();
	void SetMean();
	std::vector<double> _target;
	std::vector<Kernel*> _ker_vector;
	void InitializeGraphs();
	std::vector<double> _X_AXIS;
	void SetXAxis();
	double _error_fit(double);
	double KT_CNinv_K(double x);
	void InitializeErrorFit();
	GP1D_Fitter* _error_fitter = nullptr;
	
};

fit::GP1D_Fitter* fit::Minimizer::fitter = nullptr;

double fit::Minimizer::FitChi2(std::vector<double> params = {}){
	if (params.size()==0) params = fitter->_hyperparams;
	fitter->SetParameters(params);
	double val = fitter->_C_N.Determinant();
	//val = 1.000*TMath::Power(val, 0.50);

	for (int i = 0; i < fitter->_target.size(); i++){
		for (int j = 0; j < fitter->_target.size(); j++){
			val = val * TMath::Exp(-0.500*(fitter->_target[i] * fitter->_C_N_inv[i][j] * fitter->_target[j]));
		}
	}


	return 1.000000000*val;
	
}

double GP1D_Fitter::BETA = 10.000;
double fit::Minimizer::FitBeta(const double* xx){
	int x = GP1D_Fitter::NPARS;
	//std::cout << "npasrs" << GP1D_Fitter::NPARS << std::endl;
	std::vector<double> pars;
	for (int x = 0; x < GP1D_Fitter::NPARS; x++){
		pars.push_back(xx[x]);
	}

	GP1D_Fitter::BETA = xx[0];
	
	return FitChi2();
}

double fit::Minimizer::FitChi1(const double* xx){
	int x = GP1D_Fitter::NPARS;
	//std::cout << "npasrs" << GP1D_Fitter::NPARS << std::endl;
	std::vector<double> pars;
	for (int x = 0; x < GP1D_Fitter::NPARS; x++){
		pars.push_back(xx[x]);
	}


	return FitChi2(pars);
}

int GP1D_Fitter::NPARS = 1;

}; //NAMESPACE FIT

#endif
