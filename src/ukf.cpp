#include <iostream>
#include "ukf.h"
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <string>

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.55;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // State dimension
  n_x_= 5;

  // Augmented state dimension
  n_aug_= 7;

  // Weights of sigma points

  weights_= VectorXd(2*n_aug_+1);

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // the current NIS for radar
  NIS_radar_ = 0;

  // the current NIS for laser
  NIS_laser_ = 0 ;

  // time when the state is true, in us
  time_us_ = 0;

  Xsig_pred_ = MatrixXd(15, 5);

  //Initialize P
  P_ <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
           -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
            0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
            -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
            -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
	if (!is_initialized_){



		if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ ){
			float theta = meas_package.raw_measurements_[1];
			float rho = meas_package.raw_measurements_[0];
			float px = rho*cos(theta);
			float py = rho*sin(theta);


			x_(0) = px;
			x_(1) = py;
			x_(2) = 0.1;
			x_(3) = 0.1;
			x_(4) = 0.1;
			time_us_ = meas_package.timestamp_;
			//std::cout<<"X_initial: \n"<<x_<<"\n";


		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ ) {

			x_(0) = meas_package.raw_measurements_[0];
			x_(1) = meas_package.raw_measurements_[1];
			x_(2) = 0.1;
			x_(3) = 0.1;
			x_(4) = 0.1;
			time_us_ = meas_package.timestamp_;


		}

		else{

		  return;
		}

		is_initialized_ = true;
	}

	else{

	  double delta_t = (meas_package.timestamp_ - time_us_) /1000000.0;

	  time_us_= meas_package.timestamp_;


	  if(meas_package.sensor_type_==MeasurementPackage::LASER && use_laser_ ){
	        //std::cout<<"Meas:\n"<<x_<<"\n";
		   	Prediction(delta_t);
	        //Prediction(delta_t);
	        UpdateLidar(meas_package);
	      }
	  else if (meas_package.sensor_type_==MeasurementPackage::RADAR && use_radar_ )

	  {

	         Prediction(delta_t);
	         UpdateRadar(meas_package);
	      }

	  else return;

	}
}



MatrixXd UKF::GenerateSigmaPoints() {

  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);


  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  //set first column of sigma point matrix
  Xsig.col(0) = x_;

  //set remaining sigma points

  for(int i=0; i<n_x_; i++)
  {
      Xsig.col(i+1) = x_+ sqrt(lambda_+n_x_)*A.col(i);
      Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_)*A.col(i);

  }

  //write result
  return Xsig;
}


MatrixXd UKF::AugmentedSigmaPoints() {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);


  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();


  double lambda = 3 - n_aug_;
  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda+n_aug_) * L.col(i);
  }

  return Xsig_aug;

}


MatrixXd UKF::SigmaPointPrediction(MatrixXd Xsig_aug, double delta_t) {

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //std::cout<<"Delta_t:\n"<<delta_t<<"\n";

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  //print result

  //write result
  return Xsig_pred;

}



/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */

void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //Obtain sigma points
  MatrixXd Xsig = MatrixXd(11, 5);
  MatrixXd Xsig_aug = MatrixXd(15, 7);

  //Generate sigma points in the original state
  Xsig = GenerateSigmaPoints();
  //std::cout<<"Sigma points: \n"<<Xsig<<"\n";

  //Add v noise
  Xsig_aug = AugmentedSigmaPoints();
  //std::cout<<"Sigma augmented:\n"<<Xsig_aug<<"\n";

  //Predict sigma points.
  Xsig_pred_= SigmaPointPrediction(Xsig_aug, delta_t);
  //std::cout<<"X_sig_pred:\n"<<Xsig_pred_<<"\n";

  //Calculate the mean and the covariance of the predicted state
  //create vector for predicted state
  VectorXd x_pred = VectorXd(5);

  //create vector for predicted covariance
  MatrixXd P_pred = MatrixXd(5, 5);

  //define spreading parameter
  double lambda = 3 - n_aug_;

  // set weights
  double weight_0 = lambda/(lambda+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda);
    weights_(i) = weight;
  }



  //predicted state mean
  x_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_pred = x_pred+ weights_(i) * Xsig_pred_.col(i);
   }

  //predicted state covariance matrix
  P_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_pred;
    //angle normalization

    x_diff(3) = x_diff(3)- std::ceil(double((x_diff(3)-M_PI)/(2.*M_PI)))*2.*M_PI;
    //std::cout<<x_diff(3)<<"\n";

    P_pred = P_pred + weights_(i) * x_diff * x_diff.transpose() ;



  }
  P_ = P_pred;
  x_ = x_pred;

  //std::cout<<"x_predicted:\n"<<x_<<"\n";
  //std::cout<<"P_predicted:\n"<<P_<<"\n";

  }


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  //set measurement dimension, radar can measure px, py
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);


  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);


    // measurement model
    Zsig(0,i) = p_x ;  //px
    Zsig(1,i) = p_y ;  //py

  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga point
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<    std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;

    S = S + R;

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0);

    for(int i=0; i<2 * n_aug_ + 1; i++) {
    	VectorXd z_diff = Zsig.col(i) - z_pred;
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
        }



    MatrixXd K = Tc*(S.inverse());

    x_ = x_ + K*(meas_package.raw_measurements_-z_pred);
    P_ = P_ - K*S*K.transpose();


    NIS_laser_ = ((meas_package.raw_measurements_-z_pred).transpose())*S.inverse()*(meas_package.raw_measurements_-z_pred);






}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  //std::cout<<"Radar";
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //define spreading parameter
  double lambda = 3 - n_aug_;

  //set vector for weights
  double weight_0 = lambda/(lambda+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_+lambda);
    weights_(i) = weight;
  }

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r



    if(p_x>0.01){
    	 Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    	 }
    else {Zsig(1,i) = M_PI_2;}

    if(sqrt(p_x*p_x + p_y*p_y)>0.01){
    	Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }
    else {Zsig(2,i) =0;}



  }
  	//std::cout<<"Z_sig:\n"<<Zsig <<"\n";



    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);

    for (int i=0; i < 2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
      VectorXd z_diff = Zsig.col(i) - z_pred;
      S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<    std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_,0,
            0, 0,std_radrd_*std_radrd_;

    S = S + R;
    //std::cout<<"z_pred:\n"<<z_pred<<"\n";
    //std::cout<<"S:\n"<<S<<"\n";



    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0);
    for(int i=0; i<2 * n_aug_ + 1; i++) {

      VectorXd z_diff = Zsig.col(i) - z_pred;
      //angle normalization
      z_diff(1) = z_diff(1)- std::ceil(float((z_diff(1)-M_PI)/(2.*M_PI)))*2.*M_PI;

      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      //angle normalization
      x_diff(3) = x_diff(3)- std::ceil(float((x_diff(3)-M_PI)/(2.*M_PI)))*2.*M_PI;

      Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
    }


    MatrixXd K = Tc*(S.inverse());

    x_ = x_ + K*(meas_package.raw_measurements_-z_pred);
    P_ = P_ - K*S*K.transpose();
    NIS_radar_ = ((meas_package.raw_measurements_-z_pred).transpose())*S.inverse()*(meas_package.raw_measurements_-z_pred);
    //std::cout<<"NIS:\n"<<NIS_radar_<<"\n";
    //std::cout<<"x:\n"<<x_<<"\n";
    //std::cout<<"P:\n"<<P_ <<"\n";



}
