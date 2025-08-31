#include <iostream>
#include <tuple>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <fstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"
using namespace std;
using namespace TMath;


#define PLOT_ANGLE_XY

const string file_x = "corrected_x1x2x3.txt";
const string file_y = "corrected_y1y2y3.txt";
const string outfile = "3DAngles.txt";

/* adc values that are only important for the output file */
vector<uint16_t> hadc;
vector<uint16_t> ladc;


void read_file_x(string filename, vector<double> *x1,vector<double> *x2,vector<double> *x3,vector<double> *z1,vector<double> *z2,vector<double> *z3,vector<unsigned int> *event_number_x)
{
	std::ifstream infile;
	infile.open(filename);
        double b1, b2, b3, b4, b5, b6, b7;
	uint16_t b8, b9;
        while (infile >> b1 >> b2 >> b3 >> b4 >> b5 >> b6 >> b7 >> b8 >> b9) {
		event_number_x->push_back(b1);
		x1->push_back(b2);
		z1->push_back(b3);
		x2->push_back(b4);
		z2->push_back(b5);
		x3->push_back(b6);
		z3->push_back(b7);
		hadc.push_back(b8);
		ladc.push_back(b9);
	    }
}

void read_file_y(string filename, vector<double>* y1,vector<double>* y2,vector<double>* y3,vector<unsigned int>* event_number_y)
{
	std::ifstream infile;
	infile.open(filename);
        double b1, b2, b3, b4, b5, b6, b7;
        while (infile >> b1 >> b2 >> b3 >> b4 >> b5 >> b6 >> b7 ) {
		event_number_y->push_back(b1);
		y1->push_back(b2);
		y2->push_back(b4);
		y3->push_back(b6);
	    }
}





inline double calculate_mean(vector<double>* n)
{
	double total = 0;
	
	for (auto& elem : *n)
		total += elem;
	return total / n->size();
}

inline tuple<TVectorD,vector<double>> get_eigenvector(const vector<vector<double>*> &events)
{
	vector<double> *x, *y, *z;
	x = events[0];
	y = events[1];
	z = events[2];

	double mean_x = calculate_mean(x); 
	double mean_y = calculate_mean(y);
	double mean_z = calculate_mean(z);

	vector<double> means = {mean_x, mean_y, mean_z};
	
	
	for (auto& element : *x)
		element -= mean_x;
	for (auto& element : *y)
		element -= mean_y;
	for (auto& element : *z)
		element -= mean_z;

	TArrayD event(9);

	int j = 0, k = 0, l = 0;
	for (int i = 0; i < 9; i++){
		if (i % 3 == 2){
			event[i] = (*z)[l++];
		}
		else if (i % 3 == 1){
			event[i] = (*y)[k++];
		}
		else{
			event[i] = (*x)[j++];
		}
	}

	

	TMatrixD A(3,3); TMatrixD V(3,3);  TVectorD v(3);
	A.SetMatrixArray(event.GetArray());

	TDecompSVD svd(A);
	svd.Decompose();
	V = svd.GetV();
	v = TMatrixDColumn(V,0);
	return tuple<TVectorD, vector<double>> {v, means};

}

inline tuple<TVectorD, vector<double>> dummy_fit(const vector<vector<double>*> &events)
{
	vector<double> *x, *y, *z;
	x = events[0];
	y = events[1];
	z = events[2];

	TVectorD v(3);

	v(0) = (*x)[0] - (*x)[2];
	v(1) = (*y)[0] - (*y)[2];
	v(2) = (*z)[0] - (*z)[2];

	double norm = 0;
	for (int i = 0; i < v.GetNrows(); ++i) {
	    norm += v(i) * v(i);
	}
	norm = std::sqrt(norm);

	v *= 1/norm;

	return {v, {0,0,0}};

}

inline tuple<vector<double>, TVectorD> calculate_slope(const vector<vector<double>*> &event)
{
	TVectorD eigen_vector(3);
	vector<double> means(3);
	double polar_angle; tuple<TVectorD,vector<double>> _eigenvector = get_eigenvector(event);
	eigen_vector = get<0>(_eigenvector);
	means = get<1>(_eigenvector);
	polar_angle = eigen_vector[2]; // ACos

	return tuple<vector<double>, TVectorD> {means, eigen_vector};
}



void Angle3D()
{

	auto histo_lowr_bound = -40  ; // for setting bound on histo and fitting function
	auto histo_uppr_bound =  40  ;
	TH1D *histo_angles = new TH1D ("polar angle distribution", "Polar Angle Distribution(deg)", 100, histo_lowr_bound, histo_uppr_bound);

#ifdef PLOT_ANGLE_XY
	TH1D *histo_angles_x = new TH1D ("X Z angle distribution", "X Z angle distribution (rad)", 100, histo_lowr_bound, histo_uppr_bound);
	TH1D *histo_angles_y = new TH1D ("Y Z angle distribution", "Y Z angle distribution (rad)", 100, histo_lowr_bound, histo_uppr_bound);
#endif // #ifdef PLOT_ANGLE_XY

	vector<double> x1,x2,x3;
	vector<double> y1,y2,y3;
	vector<double> z1,z2,z3;
	vector<unsigned int> event_num_x, event_num_y;

	bool plot_line = true;


	read_file_x(file_x,&x1,&x2,&x3,&z1,&z2,&z3,&event_num_x);
	read_file_y(file_y,&y1,&y2,&y3,&event_num_y);


	if(y1.size() == x1.size()  && event_num_x == event_num_y)
		std::cout << "successfully read files" << std::endl;
	else{
		cout << "event numbers or x and y vector sizes differ: aborting" << endl;
		abort();
	}


	TVectorD v(3); 
	vector<double> meanz; 

	ofstream ofile;

	vector<vector<double>*> test_event(3);


	ofile.open(outfile);
	for (unsigned int i = 0; i < x1.size(); i++){
		vector<double> xvec = {x1[i], x2[i], x3[i]}; 
                vector<double> yvec = {y1[i], y2[i], y3[i]};
                vector<double> zvec = {z1[i], z2[i], z3[i]};

		vector<vector<double>*> event = {&xvec,&yvec,&zvec};
		tuple<vector<double>,TVectorD> output = calculate_slope(event);
		auto eigen_vector = get<1>(output);

		/* calculate polar angle */
		double z = eigen_vector[2]; 

		double angle;
		angle = atan2(eigen_vector[0], z) * 180/TMath::Pi(); 

		ofile << event_num_x[i] << "\t" << double(angle)<< "\t" << hadc[i]<< "\t" << ladc[i] << endl;
		histo_angles->Fill(angle);
	#ifdef PLOT_ANGLE_XY
		double anglex = ACos(eigen_vector[0]) * 180/ TMath::Pi();
		double angley = ACos(eigen_vector[1]) * 180/ TMath::Pi();

		if  ((anglex < 88 || anglex > 92) && (angley < 88 || angley > 92) && (fabs(ACos(sqrt( abs(1 - Cos(anglex) * Cos(anglex) - Cos(angley) * Cos(angley)))) - angle) > 2))
			cout << anglex << " " << angley << " " << angle << endl;

		histo_angles_x->Fill(anglex);
		histo_angles_y->Fill(angley);
	#endif // #ifdef PLOT_ANGLE_XY
		if (plot_line)
		{
			meanz = get<0>(output);
			v = get<1>(output);
			for (int i = 0; i < 3; i++)
				test_event[i] = event[i];
			plot_line = false;
		}
	}


	
	TGraph2D * gr1 = new TGraph2D();
	for (int i = 0; i < 10000; i++)
	{
		 double x_,y_,z_;
		 x_ = meanz[0] + v[0] * (i - 5000) * 10;  
		 y_ = meanz[1] + v[1] * (i - 5000) * 10;
		 z_ = meanz[2] + v[2] * (i - 5000) * 10;

		 gr1->SetPoint(i, x_, y_, z_);
	}


	// Plot points with mean 0
	TGraph2D * gr = new TGraph2D();

	gr->SetPoint(0,x1[0],y1[0],z1[0]);
	gr->SetPoint(1,x2[0],y2[0],z2[0]);
	gr->SetPoint(2,x3[0],y3[0],z3[0]);


	ofile << "event_number" << "\t" << "polar_angle(deg)" << "\t" <<  "hadc" << "\t" << "ladc" << endl;
	ofile.close();


	
	// Plot polar angle distribution
	TF1 *g1 = new TF1("g1", "gaus", histo_lowr_bound,histo_uppr_bound);
	TF1 *p1 = new TF1("p1", "[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1.)", 0, 3.14);
	g1->SetParameters(1,1);
	p1->SetParameters(1,1,1);
	histo_angles->Fit(g1,"r");

	histo_angles-> GetXaxis()->SetTitle("Polar Angle (deg)");
	histo_angles-> GetYaxis()->SetTitle("Counts");
	gStyle->SetOptFit(1);
	histo_angles->Draw();

#ifdef PLOT_ANGLE_XY

	TCanvas *c9 = new TCanvas();
	histo_angles_x->Fit(g1,"r");
	histo_angles_x-> GetXaxis()->SetTitle("Polar Angle (rads)");
	histo_angles_x-> GetYaxis()->SetTitle("Counts");
	histo_angles_x->Draw();

	TCanvas *c7 = new TCanvas();
	histo_angles_y->Fit(g1,"r");
	histo_angles_y-> GetXaxis()->SetTitle("Polar Angle (rads)");
	histo_angles_y-> GetYaxis()->SetTitle("Counts");
	histo_angles_y->Draw();
#endif //  #ifdef PLOT_ANGLE_XY

	TView *view = TView::CreateView(1);
	TCanvas *c1 = new TCanvas();
	//gr1->SetMargin(0.2);
	gr->SetMarkerStyle(2);
	gr->SetMarkerColor(kBlue);
	gr->GetXaxis()->SetTitle("X (mm)");
	gr->GetYaxis()->SetTitle("Y (mm)");
	gr->GetZaxis()->SetTitle("Z (mm)");
	gr->Draw("P0");
	gr1->SetMarkerStyle(2);
	gr1->SetMarkerColor(kRed);
	gr1->SetMarkerSize(0.5);
	gr1->Draw("P0 SAME");
}
