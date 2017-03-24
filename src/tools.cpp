#include <iostream>
#include "tools.h"

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    vector<VectorXd> estimations_4d;
	rmse << 0,0,0,0;

    // check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size	
    if ((estimations.size()==0) |(estimations.size()!=ground_truth.size())){
    	cout<<"estimation vector size should not be zero OR estimation vector size should be equal to ground truth vector size"<<"\n";
    }

    for (int i=0; i < estimations.size(); ++i) {
      VectorXd converted(4);
      converted << estimations[i][0],                           // px
                   estimations[i][1],                           // py
                   cos(float(estimations[i][3]))*estimations[i][2],    // vx
                   sin(float(estimations[i][3]))*estimations[i][2];    // vy
      estimations_4d.push_back(converted);
    }




	//accumulate squared residuals
	for(int i=0; i < estimations_4d.size(); ++i){

        VectorXd residual=estimations_4d[i]-ground_truth[i];

        residual = residual.array()*residual.array();
		rmse += residual;
		//cout<<"RMSE: "<<rmse<<"\n";
        
	}

	//calculate the mean
	int n = estimations_4d.size();
	rmse = rmse/n;

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

