#include "kalman_filter.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_=F_*x_;
  MatrixXd Ft=F_.transpose();
  P_=F_*P_*Ft+Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y=z-H_*x_;
  MatrixXd Ht=H_.transpose();
  MatrixXd S=H_*P_*Ht+R_;
  MatrixXd Si=S.inverse();
  MatrixXd K=P_*Ht*Si;
  //New State
  x_=x_+(K*y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_=(I-K*H_)*P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  
  // If rho == 0, skip the update step to avoid dividing by zero
  if( px == 0. && py == 0. ) {
    return;
  }
  
  float rho=sqrt(px*px+py*py);
  float theta=atan2(py,px);
  float rho_dot=(px*vx+py*vy)/rho;
  
  VectorXd h_func(3);
  h_func<<rho,theta,rho_dot;
  VectorXd y=z-h_func;
  
  // Normalize the angle
  while (y(1)>M_PI) {
    y(1) -= 2 * M_PI;
  }
  while (y(1)<-M_PI) {
    y(1) += 2 * M_PI;
  }
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_*P_*Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_*Ht*Si;

  // Compute new state
  x_ = x_ + ( K*y );
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = ( I - K*H_ )*P_;
  

}


