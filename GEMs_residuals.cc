#include <iostream>
#include <unordered_map>
#include <filesystem>
#include <cmath>
#include <vector>
#include <fstream>
#include <format>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"

// uncomment below to show residual plots
//#define PLOT_RESIDUALS
//#define FILTER_BY_ANGLE /* ONLY FOR AFTER RUNNING THIS SCRIPT ON DATA AND RUNNING THE OUTPUT THROUGH ANGLE3D.CC*/

/* Basic Definitions */
gErrorIgnoreLevel = 6001;
const string GEMS[] = {"GEM-I", "GEM-II", "GEM-III"};
const string COORDS[] = {"X", "Y", "Z"};

enum coordinate_t { X, Y, Z };
enum GEM_t 	  { GEM1, GEM2, GEM3 };
enum DATA_t 	  { ERR, VAL};

/* Definition of Constants */
constexpr double Z3 = 0.0;
constexpr double Z2 = 230.0;
constexpr double Z1 = Z2 + 200;   
GEM_t target_gem = GEM2; // GEM to shift during gradient descent

/* SETTINGS */
const int run_number = 4000; // set run number for input files
const string corrected_filename_y =  "corrected_y1y2y3.txt";// output file of corrected y1y2y3
const string corrected_filename_x =  "corrected_x1x2x3.txt";// output file of corrected x1x2x3
const string filter_file = "3DAngles.txt"; // for filtering normal events from infile - for after running Angle3D.C

/* for output files */
double offset_x;
double offset_y;

void read_filter (const string &filename, vector<uint32_t>& filter)
{
	std::ifstream infile;
	infile.open(filename);
	uint32_t a1,a2;
	while (infile >> a1 >> a2)
		filter.push_back(a1);
	
}

tuple<unordered_map<uint32_t, vector<vector<double>>>, unordered_map<uint32_t, vector<uint16_t>>> gem_residuals_sorting() {
	const TString infiles[] = {
		Form("output_file_run_%d_x1y1.txt", run_number),
		Form("output_file_run_%d_x2y2.txt", run_number),
		Form("output_file_run_%d_x3y3.txt", run_number)
	};
	int counter = 0;
	std::vector<double> run_number_sorted;
	std::unordered_map<uint32_t, std::vector<uint16_t>> adc_;
	auto fill_hadc = true;
	
	double x1, x2, x3;
	double run_number_y1, run_number_y2, run_number_y3;
	double y1, y2, y3;
	std::unordered_map<uint32_t, std::vector<double>> gems[3];
	
	uint32_t run_num;
	uint16_t hadc, ladc; 
	double  a2, a3, a4, a5;
	for (int i = 0; i < 3; i++){
		std::ifstream infile1;
		infile1.open(infiles[i]);
		while (infile1 >> run_num >> a2 >> a3 >> a4 >> a5 >> hadc >> ladc) {
		    gems[i][run_num] = {a2 *0.39062500, a3 * 0.39062500};
		    if (fill_hadc)
			    adc_[run_num]  = {hadc, ladc};
		}
		if (fill_hadc)
			fill_hadc = false;
	}
	int foo = 0;
	
	uint32_t largest_run = 0;
	for (auto &pair: gems[2])
	{
		if (foo++ < 4)
			cout << pair.first << " " << pair.second[0] / .391 << " " << pair.second[1]/.391 << endl;
	    	if (pair.first > largest_run)
	    		largest_run = pair.first;
	}

	unordered_map<uint32_t, vector<vector<double>>> events;
	unordered_map<uint32_t, vector<uint16_t>> adc;
	
	std::ofstream infile1("sorted_data_file_run_number_x1_x2_x3.txt");
	std::ofstream infile2("sorted_data_file_run_number_y1_y2_y3.txt");
	
	for (uint32_t i = 0; i < largest_run; i++)
	{
	        if ((gems[0].find(i) != gems[0].end()) && (gems[1].find(i) != gems[1].end()) && (gems[2].find(i) != gems[2].end())){
			events[i] = vector<vector<double>>(3);
			events[i][X] = {gems[GEM1][i][X], gems[GEM2][i][X], gems[GEM3][i][X]};
			events[i][Y] = {gems[GEM1][i][Y], gems[GEM2][i][Y] - 8, gems[GEM3][i][Y]};
			events[i][Z] = {Z1, Z2, Z3};
			adc[i] 	     = adc_[i];

	}
	}

	for (auto &elem: events){
		infile1 << elem.first << "\t" << elem.second[X][GEM1] << "\t" << elem.second[Z][GEM1] << "\t" << elem.second[X][GEM2] << "\t" << elem.second[Z][GEM2] <<"\t" << elem.second[X][GEM3] << "\t" << elem.second[Z][GEM3] << endl;
		infile2 << elem.first << "\t" << elem.second[Y][GEM1] << "\t" << elem.second[Z][GEM1] << "\t" << elem.second[Y][GEM2] << "\t" << elem.second[Z][GEM2] <<"\t" << elem.second[Y][GEM3] << "\t" << elem.second[Z][GEM3] << endl;
	}

	//int i = 0;
	//for (auto &elem : events)
	//{
	//	if (i++ < 4){
	//	Double_t x[3] = {elem.second[X][GEM1], elem.second[X][GEM2], elem.second[X][GEM3]};
	//	Double_t y[3] = {elem.second[Y][GEM1], elem.second[Y][GEM2], elem.second[Y][GEM3]};
	//	Double_t z[3] = {elem.second[Z][GEM1], elem.second[Z][GEM2], elem.second[Z][GEM3]};
	//	TGraph *xz = new TGraph(3,z, x);
	//	TGraph *yz = new TGraph(3,z, y);

	//	TCanvas *c1 = new TCanvas();
	//	TString Titlex = Form("XZ plot event %d", elem.first);
	//	TString Titley = Form("YZ plot event %d", elem.first);
	//	xz->SetTitle(Titlex);
	//	yz->SetTitle(Titley);
	//	xz->GetXaxis()->SetTitle("Z");
	//	yz->GetXaxis()->SetTitle("Z");
	//	xz->GetYaxis()->SetTitle("X");
	//	yz->GetYaxis()->SetTitle("Y");
	//	xz->SetMarkerStyle(20);
	//	xz->SetMarkerSize(1);
	//	xz->SetMarkerColor(kBlue);
	//	yz->SetMarkerStyle(20);
	//	yz->SetMarkerSize(1);
	//	yz->SetMarkerColor(kBlue);
	//	xz->Draw("AP");
	//	TCanvas *c2 = new TCanvas();
	//	yz->Draw("AP");
	//	}
	//}

	return {events,adc};
}

void write_sorted_file(const string &filename, unordered_map<uint32_t, vector<vector<double>>>&data, unordered_map<uint32_t, vector<uint16_t>>* adc = nullptr )
{

	std::ofstream ofile;
	ofile.open(filename);
	if (adc){
			for (auto &elem: data)
				ofile << elem.first << "\t"
				      << elem.second[X][0] << "\t" << Z1 << "\t"
				      << elem.second[X][1] << "\t" << Z2 << "\t"
				      << elem.second[X][2] << "\t" << Z3 << "\t"
				      << (*adc)[elem.first][0] << "\t" << (*adc)[elem.first][1] << endl;
		ofile << "event_num" << "\t"
		      << "x1" << "\t" << "Z1" << "\t"
		      << "x2" << "\t" << "Z2" << "\t"
		      << "x3" << "\t" << "Z3" << "\t"
		      << "hadc"  << "\t" << "ladc" << endl;
		ofile << GEMS[target_gem] << " offset by " << offset_x << " (mm)" << endl;
	}

	else{
		for (auto &elem:data)
			ofile << elem.first << "\t"
			      << elem.second[Y][0] << "\t" << Z1 << "\t"
			      << elem.second[Y][1] << "\t" << Z2 << "\t"
			      << elem.second[Y][2] << "\t" << Z3 << "\t"
			      << endl;
		ofile << "event_num" << "\t"
		      << "y1" << "\t" << "Z1" << "\t"
		      << "y2" << "\t" << "Z2" << "\t"
		      << "y3" << "\t" << "Z3" << "\t"
		      << endl;
		ofile << GEMS[target_gem] << " offset by " << offset_y << " (mm)" << endl;
	}
}

// takes two known positions (v1,v2) and gives error of the third position (v3)
void collect_errors(const unordered_map<uint32_t, vector<vector<double>>> &data, vector<vector<vector<double>*>> &error, coordinate_t C, GEM_t G)
{
	auto len = data.size();

	vector<double> *target_error = error[C][G];

	double slope, intercept, expected;
	int i = 0;
	for (auto &elem:data)
	{
		auto v1     = elem.second[C][(G + 1) %3];
		auto v2     = elem.second[C][(G + 2) %3];
		auto target = elem.second[C][G];

		auto z1     = elem.second[Z][(G + 1) %3];
		auto z2     = elem.second[Z][(G + 2) %3];
		auto z3     = elem.second[Z][G];
		slope = ( v1- v2)/(z1 - z2);
		intercept = v1 - z1 * slope;
		expected = slope * z3 + intercept;
		(*target_error)[i++] = expected - target;
	}
}


// fit data to gaussian and return mean. has option to make plot of the distribution
double extract_gauss_mean(const string var_name, vector<double> &v1, const double plt_lower,  const double plt_upper, bool show_plot = false, const uint32_t num_bins = 250)
{
	// Plot distribution of errors.
	string hist_title = "Distribution of errors of " + var_name;
	string hist_name = "error_dist_" + var_name;
	TH1D *error_distribution = new TH1D(hist_name.c_str(), hist_title.c_str(), num_bins, plt_lower,plt_upper);


	// fill distribution
	for (auto &elem: v1){
		error_distribution->Fill(elem);
	}

	gStyle->SetOptFit(1); // show mean on legend

	//fit to gaussian
	TF1 *g1 = new TF1 ("g1","gaus",plt_lower,plt_upper);
	if (show_plot)
		error_distribution->Fit(g1, "rQ");
	else
		error_distribution->Fit(g1, "RQ0");



	auto mean = g1->GetParameter(1);


	if (show_plot){
		TCanvas *c1 = new TCanvas();
		error_distribution->GetXaxis()->SetTitle("residual (mm)");
		error_distribution->GetYaxis()->SetTitle("counts");
		error_distribution->Draw();
	}

	if (!show_plot)
		delete error_distribution;
	return mean;
}

inline double total_error(vector<vector<vector<double>*>> &error, coordinate_t C)
{
	
	return abs(extract_gauss_mean("dummy", *error[C][0], -20,20)) +  // correct returns the mean of the error dist
               abs(extract_gauss_mean("dummy", *error[C][1], -20,20)) + 
               abs(extract_gauss_mean("dummy", *error[C][2], -20,20)) ; 		
}

/* v1 is perturbed vector */
double compute_gradient_1d(
	unordered_map<uint32_t, vector<vector<double>>> &data, coordinate_t C,
	vector<vector<vector<double>*>> &error,
	double (*error_function)(unordered_map<uint32_t, vector<vector<double>>>&,coordinate_t, vector<vector<vector<double>*>>&),
	double epsilon = 1e-6)
{
	// Perturb v1
	auto  pv1 = data, mv1 = data; 
	for (auto& val : pv1)
		val.second[C][target_gem] += epsilon;
	for (auto& val : mv1) 
		val.second[C][target_gem] -= epsilon;


	/* Compute finite difference approximation to the gradient with midpoint rule */
	double dfdv = (error_function(pv1, C,error)
		     - error_function(mv1, C,error))/(2 * epsilon);

	return dfdv;
}

/* error vectors are just dummies so I can pass the function pointer to gradient descent */
/* NOTE probably faster if I calculate chi2 by hand instead of making a plot and fitting it */
int plot = 0;
inline double chi2_fit(unordered_map<uint32_t, vector<vector<double>>> &data, coordinate_t C,vector<vector<vector<double>*>> &error)
{
	const Int_t nPts = 3;
	double total_chi_sqr = 0;

	//TF1 fitFunc = TF1("fitFunc", "pol1", 0, 300);
	//for (auto &elem:data){
	//	double z_vals[nPts]     = { elem.second[Z][GEM1], elem.second[Z][GEM2], elem.second[Z][GEM3]};


	//	double coord_vals[nPts] = { elem.second[C][GEM1], elem.second[C][GEM2], elem.second[C][GEM3]};
	//	TGraph graph = TGraph(nPts, coord_vals, z_vals);


	//	graph.Fit(&fitFunc, "Q0"); 

	//	double chi2 = fitFunc.GetChisquare();

	//	total_chi_sqr += chi2;

	//}
	double test_total_chi_sqr = 0;
	for (auto &elem:data){
		double z_vals[nPts]     = { elem.second[Z][GEM1], elem.second[Z][GEM2], elem.second[Z][GEM3]};


		double coord_vals[nPts] = { elem.second[C][GEM1], elem.second[C][GEM2], elem.second[C][GEM3]};

		TVectorD v(nPts);
		TMatrixD A(nPts,2);

		for (auto i = 0; i < nPts; i++){
			A(i, 1) = 1.0;
			A(i, 0) = coord_vals[i];
			v(i)    = z_vals[i];
		}

		TDecompSVD svd(A);
		Bool_t ok;
		TVectorD coeffs = svd.Solve(v, ok);

		double slope = coeffs[0];
		double intercept = coeffs[1];
		double sigma = 1;
		
		double chi2 = 0;
		for (auto i = 0; i < nPts; i++)
		{
			chi2 += (z_vals[i] - slope * coord_vals[i] - intercept) * (z_vals[i] - slope * coord_vals[i] - intercept)
				/ (sigma * sigma);
		}

		test_total_chi_sqr += chi2;

	}
	//if (test_total_chi_sqr != total_chi_sqr)
	//	cout << test_total_chi_sqr - total_chi_sqr << endl;
	return test_total_chi_sqr / data.size();
}

double gradient_descent(
	unordered_map<uint32_t, vector<vector<double>>> &data,
	vector<vector<vector<double>*>> &error, coordinate_t C,
	double (*error_function)(unordered_map<uint32_t, vector<vector<double>>>&, coordinate_t,vector<vector<vector<double>*>>& ),
	double learning_rate = 0.00005, uint32_t num_itr = 1000, double TOL = 0.0001
)
{
	auto n = data.size();
	double err = error_function(data, C, error);

	double gradient = 9999;

	double shift = 0;

	uint32_t counter = 1;

	while ( fabs(gradient) > TOL && num_itr > counter)
	{
		gradient = compute_gradient_1d(data, C, error, error_function);

		for (auto &elem: data) 
		{
			//data[X][i] -= learning_rate * gradient; 
			elem.second[C][target_gem] -= learning_rate * gradient;
			//data[Z][i] -= learning_rate * gradient;
		}

		shift -= learning_rate * gradient;

		/* collect errors */
		//collect_errors(data,error, C,GEM3);
		//collect_errors(data,error, C,GEM2);
		//collect_errors(data,error, C,GEM1);
		//err = error_function(data, C, error);
		counter++;
		cout << counter << endl;
	}
		cout << endl << "@@@@@@@ GRADIENT DESCENT COMPLETED @@@@@@@" << endl
		     << "GEM2 shifted " << shift << " in the " << COORDS[C] << " direction after " << counter << " iterations!" << endl;

		return shift;
}

void GEMs_residuals() {
	// *** POPULATE DATA VECTORS *** //
	auto output = gem_residuals_sorting();
	unordered_map<uint32_t, vector<vector<double>>> all_event = get<0>(output);// vector is in form V[X][GEMII]
	auto adc = get<1>(output);

#ifdef FILTER_BY_ANGLE 
	cout << "Running for all filtered events" << endl;
	unordered_map<uint32_t, vector<vector<double>>> filtered_events; // vector is in form V[X][GEMII]
	std::ifstream infile;
	infile.open(filter_file);
	uint32_t event_num;
	double dumy;
	while(infile >> event_num >> dumy){
	        filtered_events[event_num] = all_event[event_num];
	}

	
	auto event_map = filtered_events;
#else 
	auto event_map = all_event;

#endif // #ifndef FILTER_BY_ANGLE

	// *** CORRECT VALUES *** //
	auto n = event_map.size();
	vector<double> x1_err(n), x2_err(n), x3_err(n);
	vector<double> y1_err(n), y2_err(n), y3_err(n);
	vector<vector<double>*> y_err = {&y1_err, &y2_err, &y3_err};
	vector<vector<double>*> x_err = {&x1_err, &x2_err, &x3_err};
	vector<vector<vector<double>*>> err = {x_err, y_err};
	vector<double> x1_corrected(n),x2_corrected(n),x3_corrected(n);

#ifdef PLOT_RESIDUALS	
	collect_errors(event_map, err, X, GEM1);
	collect_errors(event_map, err, X, GEM2);
	collect_errors(event_map, err, X, GEM3);

	collect_errors(event_map, err, Y, GEM1);
	collect_errors(event_map, err, Y, GEM2);
	collect_errors(event_map, err, Y, GEM3);

	extract_gauss_mean("GEM-I X before"  , x1_err, -200, 200, true);
	extract_gauss_mean("GEM-II X before" , x2_err, -200, 200, true);
	extract_gauss_mean("GEM-III X before", x3_err, -200, 200, true);

	extract_gauss_mean("GEM-I Y before"  , y1_err, -200, 200, true);
	extract_gauss_mean("GEM-II Y before" , y2_err, -200, 200, true);
	extract_gauss_mean("GEM-III Y before", y3_err, -200, 200, true);
#endif // #ifdef PLOT_RESIDUALS

	offset_x = gradient_descent(event_map, err, X, &chi2_fit);
	offset_y = gradient_descent(event_map, err, Y, &chi2_fit);

	collect_errors(event_map, err, X, GEM1);
	collect_errors(event_map, err, X, GEM2);
	collect_errors(event_map, err, X, GEM3);

	collect_errors(event_map, err, Y, GEM1);
	collect_errors(event_map, err, Y, GEM2);
	collect_errors(event_map, err, Y, GEM3);

	TH2D *histox = new TH2D("delta x vs y plot", "delta x vs y plot", 250, 0,256, 250, -100, 100);
	TH2D *histoy = new TH2D("delta y vs x plot", "delta y vs x plot", 250, 0,128, 250, -100, 100);

	auto i_ = 0;
	for (auto &elem: event_map){
		for (int i = 2; i < 3; i++){
			histox->Fill(elem.second[Y][i], (*x_err[i])[i_++]);
			histoy->Fill(elem.second[X][i], (*y_err[i])[i_++]);
		}
	}

	TCanvas *c99 = new TCanvas();
	TGraph* gr =  new TGraph();
	int point = 0;

	for (int ix = 1; ix <= histox->GetNbinsX(); ++ix) {
	    for (int iy = 1; iy <= histox->GetNbinsY(); ++iy) {
		double content = histox->GetBinContent(ix, iy);
		if (content > 0) {
		    double x = histox->GetXaxis()->GetBinCenter(ix);
		    double y = histox->GetYaxis()->GetBinCenter(iy);
		    for (int k = 0; k < (int)content; ++k) {  // add content-weighted points
			gr->SetPoint(point++, x, y);
		    }
		}
	    }
	}

	gStyle->SetOptFit(1); // show mean on legend
	// Fit the graph to a straight line
	TF1* line = new TF1("line", "pol1", histox->GetXaxis()->GetXmin(), histox->GetXaxis()->GetXmax());
	gr->Fit(line, "Q");  // "Q" to suppress fit output
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1);   
	gr->SetMarkerColor(kBlue);
	gr->SetTitle("X residual as a function of Y for GEM-I");
	gr->GetXaxis()->SetTitle("Y (mm)");
	gr->GetYaxis()->SetTitle("X residual (mm)");
	gr->Draw("AP");


	TCanvas *c9r = new TCanvas();
	TGraph* gr2 = new TGraph();
	point = 0;

	for (int ix = 1; ix <= histoy->GetNbinsX(); ++ix) {
	    for (int iy = 1; iy <= histoy->GetNbinsY(); ++iy) {
		double content = histoy->GetBinContent(ix, iy);
		if (content > 0) {
		    double x = histoy->GetXaxis()->GetBinCenter(ix);
		    double y = histoy->GetYaxis()->GetBinCenter(iy);
		    for (int k = 0; k < (int)content; ++k) {  // add content-weighted points
			gr2->SetPoint(point++, x, y);
		    }
		}
	    }
	}
	// Fit the gr2aph to a straight line
	TF1* line2 = new TF1("line", "pol1", histoy->GetXaxis()->GetXmin(), histoy->GetXaxis()->GetXmax());
	gr2->SetMarkerStyle(20);
	gr2->SetMarkerSize(1);   
	gr2->SetMarkerColor(kBlue);
	gr2->Fit(line2, "Q");  // "Q" to suppress fit output
	gr2->SetTitle("Y residual as a function of X for GEM-I");
	gr2->GetXaxis()->SetTitle("X (mm)");
	gr2->GetYaxis()->SetTitle("Y residual (mm)");
	gr2->Draw("AP");

#ifndef FILTER_BY_ANGLE

        // *** WRITE TO OUTPUT FILES *** //
        write_sorted_file(corrected_filename_y, event_map);
        write_sorted_file(corrected_filename_x, event_map, &adc);

#endif // #ifdef FILTER_BY_ANGLE

#ifdef PLOT_RESIDUALS
	extract_gauss_mean("GEM-I X after"  , x1_err, -100, 100, true);
	extract_gauss_mean("GEM-II X after" , x2_err, -100, 100, true);
	extract_gauss_mean("GEM-III X after", x3_err, -100, 100, true);

	extract_gauss_mean("GEM-I Y after"  , y1_err, -100, 100, true);
	extract_gauss_mean("GEM-II Y after" , y2_err, -100, 100, true);
	extract_gauss_mean("GEM-III Y after", y3_err, -100, 100, true);
#endif // #ifdef PLOT_RESIDUALS

}
