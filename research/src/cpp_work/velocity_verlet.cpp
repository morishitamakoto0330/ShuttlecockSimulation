/*
 * Oblique projection by using Velocity Verlet algorithm
 */

#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/opencv.hpp"

#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>

#include "./matrix.hpp"

int X_INIT, Y_INIT, x_max, y_max;
int DATA = 1;
double VX_INIT, VY_INIT;

std::string img_str = "../../res/image_simulate/0226/data" + std::to_string(DATA) + "-";
std::string tmp_str;

const double m = 0.005;      // mass
const double g = 9.80665;    // gravity
const double dt = 1.0/600;   // delta time
const double mpp = 4.17/1000;   // [m/pixel]
const double v_terminal = 6.80;  // terminal velocity

const double I = 0.00005;    // inertia
const double v_mag_x = 1.05;    // velocity magnification(x)
const double v_mag_y = 1.07;    // velocity magnification(y)
const double angle = 15.0;     // shuttle angle


// distance between center of gravity and working point of resistance force
const double l = 0.005;
// set ratio considering air anisotropy
//const double ratio = 0.1;
double ratio;
// adjustable parameter
double GAMMMA_ver, GAMMMA_hor, GAMMMA_ine;



void fix_xy(double* x, double* y, int* x_pixel, int* y_pixel);
void calc_T(double theta, Matrix GAMMMA_matrix, Matrix* m);
void update_velocity(Matrix* velocity, Matrix gravity, Matrix gammma);
void file_read(std::vector<int> *v, int data_num);
void fit_line(std::vector<int> *v, int n);
void scale(std::vector<int> *v, double theta);


int main(void)
{
	std::vector<int> v;  // reality shuttle position
	// output image
	cv::Mat img = cv::imread("../../res/white_image.png");
	// draw net line
	cv::line(img, cv::Point(1090,0), cv::Point(1090,img.rows), cv::Scalar(255,0,0), 2);

	x_max = img.cols;
	y_max = img.rows;
	
	// get shuttle position
	file_read(&v, DATA);
	scale(&v, angle/180*M_PI);
	
	// adjust VX_INIT, VY_INIT
	fit_line(&v, 5);
	
	std::vector<std::vector<double>> theta_vector, velocity_theta_vector;
	std::vector<double> error, _theta_vector, _velocity_theta_vector;

	std::vector<std::vector<std::pair<double, double>>> velocity_vector, velocity_dash_vector, g_vector;
	std::vector<std::pair<double, double>> _velocity_vector, _velocity_dash_vector, _g_vector;

	int step;
	int x1, y1, x2, y2;
	double x, y, vx, vy, x_prev, y_prev, vx_prev, vy_prev;
	double fx, fy;
	double theta, omega, torque, theta_prev, omega_prev, torque_prev;
	double sum_error, optimal_ratio, min_error;

	Matrix v_matrix(2, 1);
	Matrix g_matrix(2, 1);
	Matrix E_matrix(2, 2);
	Matrix T_matrix(2, 2);
	Matrix tmp_matrix(2, 2);
	Matrix f_matrix(2, 2);
	Matrix GAMMMA_matrix(2, 2);
	
	Matrix v_dash_matrix(2, 1);
	Matrix g_dash_matrix(2, 1);

	min_error = 0.0;


	// set parameter ---------------------------------------------------
	GAMMMA_ine = 0.0016;
	GAMMMA_ver = 0.022;
	//GAMMMA_hor = 0.0023;
	//ratio = 0.1;
	// -----------------------------------------------------------------


	// simulate until the object gets out of (x,y) range
	for(int i = 1; i <= 1; i++) {
		//GAMMMA_ine = 0.0010 + 0.0001*i;
		//GAMMMA_ver = 0.010 + 0.001*i;
		ratio = i*0.1;
		GAMMMA_hor = GAMMMA_ver*ratio;
		//GAMMMA_ver = GAMMMA_hor*ratio;
		//GAMMMA_ine = (m*g - GAMMMA_hor*v_terminal)/std::pow(v_terminal, 2);
		//if(GAMMMA_ine < 0) break;

		
		GAMMMA_matrix.setMatrix({{-1*GAMMMA_ver, 0}, {0, -1*GAMMMA_hor}});
		g_matrix.setMatrix({{0.0},{-1*dt*g}});
		
		step = 0;
		sum_error = 0.0;
		
		_theta_vector.clear();
		_velocity_theta_vector.clear();
		_velocity_vector.clear();
		_velocity_dash_vector.clear();
		
		x = 0.0;
		y = 0.0;
		x1 = 0;
		y1 = 0;
		x2 = 0;
		y2 = 0;

		vx = VX_INIT;
		vy = VY_INIT;
		
		fx = 0.0;
		fy = 0.0;

		theta = atan(vy/vx);
		omega = 0.0;

		x_prev = x;
		y_prev = y;
		vx_prev = vx;
		vy_prev = vy;

		v_matrix.setMatrix({{vx}, {vy}});
		
		theta_prev = theta;
		omega_prev = omega;
		torque = (-1)*GAMMMA_ver*l*(vx*cos(M_PI/2 - theta) - vy*sin(M_PI/2 - theta));
		torque_prev = torque;

		while((0 <= x2) && (x2 < x_max) && (0 <= y2) && (y2 < y_max)) {
		//while((0 <= x2) && (x2 < x_max) && (0 <= y2) && (y2 < 20000)) {
			
			// set Matrix
			f_matrix.unitMatrix(2, 2);
			E_matrix.unitMatrix(2, 2);
			T_matrix.unitMatrix(2, 2);
			
			g_dash_matrix.zeroMatrix(2, 1);
			
			step++;
			if(step/10.0 - v.size()/2.0 + 1.0 > 0) break;
			
			// update theta(t)--------------------------------------------------------------
			theta = theta_prev + omega_prev*dt + torque_prev*dt*dt/(2*I);
			
			// update vx(t), vy(t)---------------------------------------------------------
			
			// rotate velocity (v_x, v_y) -> (v_ver, v_hor)
			T_matrix.rotateMatrix(M_PI/2 - theta);
			T_matrix.productMatrix(v_matrix.getMatrix());
			v_dash_matrix.setMatrix(T_matrix.getMatrix());

			// update (v_ver, v_hor)
			g_dash_matrix.setMatrix({{m*g*cos(theta)}, {-1*m*g*sin(theta)}});
			update_velocity(&v_dash_matrix, g_dash_matrix, GAMMMA_matrix);

			// undo rotation (v_ver, v_hor) -> (v_x, v_y)
			E_matrix.rotateMatrix(-1*(M_PI/2 - theta));
			E_matrix.productMatrix(v_dash_matrix.getMatrix());

			v_matrix.setMatrix(E_matrix.getMatrix());
			vx = v_matrix.getElement(0, 0);
			vy = v_matrix.getElement(1, 0);


			// update torque(t)--------------------------------------------------------------
			torque = (-1)*GAMMMA_ver*l*v_dash_matrix.getElement(0, 0);
			
			
			// update omega(t)---------------------------------------------------------------
			omega = omega_prev + dt*(torque + torque_prev)/(2*I);
			
			
			// set f(t)--------------------------------------------------------------------
			calc_T(theta_prev, GAMMMA_matrix, &f_matrix);
			f_matrix.productMatrix({{vx_prev}, {vy_prev}});
			fx = f_matrix.getElement(0, 0);
			fy = f_matrix.getElement(1, 0);
			
			// update x(t), y(t) ([m] -> [pixel])---------------------------------------------
			// resistance considering rotation
			x = x_prev + vx_prev*dt + fx*dt*dt/(2*m);
			y = y_prev + vy_prev*dt + fy*dt*dt/(2*m);

			// set pixel position (x[m], y[m]) -> (x[pixel], y[pixel])
			x1 = X_INIT;
			x2 = X_INIT;
			y1 = Y_INIT;
			y2 = Y_INIT;
			
			fix_xy(&x_prev, &y_prev, &x1, &y1);
			fix_xy(&x, &y, &x2, &y2);

			std::cout << "step:" << step << "-------------------------------" << std::endl;
			std::cout << "(x,y)=(" << x << "," << y << ")";
			std::cout << "(vx,vy)=(" << vx << "," << vy << ")" << std::endl;
			std::cout << "(fx,fy)=(" << fx << "," << fy << ")" << std::endl;
			std::cout << "θ=" << theta/M_PI*180 << ",_θ=" << atan(vy/vx)/M_PI*180 << ",ω=" << omega << ",T=" << torque << std::endl;

			// calculate error
			if(step%10 == 0) {
				if((step/10.0 - v.size()/2.0 + 1.0) <= 0) {
					int _step = step/10;
					int _x = v[_step*2];
					int _y = v[_step*2+1];
					sum_error += sqrt(pow(x1-_x, 2) + pow(y1-_y, 2));
			
					// draw shuttle trajectory
					// simulate
					cv::circle(img, cv::Point(x1, y1), 4, cv::Scalar(i*25,0,0), 2);
					// reality
					cv::circle(img, cv::Point(_x, _y), 4, cv::Scalar(0,0,255), 2);
				}
				/*
				if((0 <= x1) && (x1 < x_max) && (0 <= y1) && (y1 < y_max)) {
					cv::circle(img, cv::Point(x1, y1), 4, cv::Scalar(i*25,0,0), 2);
				}
				*/


				// save value
				_theta_vector.push_back(theta);
				if(vx < 0) _velocity_theta_vector.push_back(atan(vy/vx) - M_PI);
				else _velocity_theta_vector.push_back(atan(vy/vx));
				_velocity_vector.push_back(std::make_pair(vx, vy));
				_velocity_dash_vector.push_back(std::make_pair(v_dash_matrix.getElement(0, 0), v_dash_matrix.getElement(1, 0)));
				_g_vector.push_back(std::make_pair(g_dash_matrix.getElement(0, 0), g_dash_matrix.getElement(1, 0)));
			}


			// next step
			x_prev = x;
			y_prev = y;
			vx_prev = vx;
			vy_prev = vy;
			theta_prev = theta;
			omega_prev = omega;
			torque_prev = torque;
		}
		
		// save vector
		theta_vector.push_back(_theta_vector);
		velocity_theta_vector.push_back(_velocity_theta_vector);
		velocity_dash_vector.push_back(_velocity_dash_vector);
		velocity_vector.push_back(_velocity_vector);
		g_vector.push_back(_g_vector);


		// distance error in actual trajectory and simulation trajectory
		sum_error = sum_error*mpp;
		error.push_back(sum_error);
		if(min_error == 0.0) min_error = sum_error;
		if(min_error > sum_error) {
			min_error = sum_error;
		}
		std::cout << "distance error:" << sum_error << std::endl;
	}

	// disp vector value
	double _v_x, _v_y, _v_ver, _v_hor, _g_ver, _g_hor, _f_ver, _f_hor;

	// output stream
	std::ofstream ofs("../../res/graph_data/0226/tmp.dat", std::ios::out);

	std::cout << "vector value" << std::endl;
	for(int i = 0; i < velocity_dash_vector.size(); i++) {
		std::cout << i << " ------------" <<  std::endl;
		for(int j = 0; j < velocity_dash_vector[i].size(); j++) {
			_v_x = velocity_vector[i][j].first;
			_v_y = velocity_vector[i][j].second;
			_v_ver = velocity_dash_vector[i][j].first;
			_v_hor = velocity_dash_vector[i][j].second;
			_g_ver = g_vector[i][j].first;
			_g_hor = g_vector[i][j].second;
			
			
			std::cout << j << " ";
			std::cout << _v_ver << " ";
			std::cout << _v_hor << " ";
			std::cout << theta_vector[i][j]/M_PI*180 << " ";
			std::cout << velocity_theta_vector[i][j]/M_PI*180 << " ";
			//std::cout << (theta_vector[i][j] - velocity_theta_vector[i][j])/M_PI*180 << " ";
			//std::cout << _v_x << " ";
			//std::cout << _v_y << " ";
			//std::cout << sqrt(pow(_v_ver, 2) + pow(_v_hor, 2)) << " "; 
			//std::cout << sqrt(pow(_v_x, 2) + pow(_v_y, 2)) << " "; 
			//std::cout << _g_ver << " ";
			//std::cout << _g_hor << " ";
			
			_f_ver = -1*GAMMMA_ver*_v_ver;
			_f_hor = -1*GAMMMA_hor*_v_hor - GAMMMA_ine*abs(_v_hor)*_v_hor;
			//std::cout << sqrt(_f_ver*_f_ver + _f_hor*_f_hor) << " ";
			std::cout << std::endl;
			

			ofs << j << " ";
			ofs << _v_ver << " ";
			ofs << _v_hor << " ";
			ofs << theta_vector[i][j]/M_PI*180 << " ";
			ofs << velocity_theta_vector[i][j]/M_PI*180 << " ";
			ofs << std::endl;
		}
		std::cout << std::endl;
		ofs << std::endl;
	}

	std::cout << "index : distance error, average" << std::endl;
	
	for(int i = 0; i < error.size(); i++) {
		std::cout << i << ": " << error[i] << ", " << error[i]/(v.size()/2) << std::endl;
	}
	
	cv::circle(img, cv::Point(X_INIT,Y_INIT), 4, cv::Scalar(0,0,0), 3);
	
	tmp_str = std::to_string(v_mag_x);
	tmp_str.erase(tmp_str.begin() + 1);
	img_str += "v_mag_x=" + tmp_str;
	
	tmp_str = std::to_string(abs(v_mag_y));
	tmp_str.erase(tmp_str.begin() + 1);
	img_str += "v_mag_y=" + tmp_str;
	
	tmp_str = std::to_string(ratio);
	if(ratio >= 10.0) tmp_str.erase(tmp_str.begin() + 2);
	else tmp_str.erase(tmp_str.begin() + 1);
	img_str += "ratio=" + tmp_str;
	
	tmp_str = std::to_string(GAMMMA_ver);
	tmp_str.erase(tmp_str.begin() + 1);
	img_str += "GAMMMA_ver=" + tmp_str;

	tmp_str = std::to_string(GAMMMA_ine);
	tmp_str.erase(tmp_str.begin() + 1);
	img_str += "GAMMMA_ine=" + tmp_str;

	tmp_str = std::to_string(l);
	tmp_str.erase(tmp_str.begin() + 1);
	img_str += "l=" + tmp_str;
	
	tmp_str = std::to_string(angle);
	if(angle >= 10.0) tmp_str.erase(tmp_str.begin() + 2);
	else tmp_str.erase(tmp_str.begin() + 1);
	img_str += "angle=" + tmp_str;

	tmp_str = std::to_string(I);
	tmp_str.erase(tmp_str.begin() + 1);
	img_str += "I=" + tmp_str;
	
	img_str += ".png";
	cv::imwrite(img_str, img);
	

	
	std::cout << std::endl;
	std::cout << "X_INIT=" << X_INIT << std::endl;
	std::cout << "Y_INIT=" << Y_INIT << std::endl;
	std::cout << "VX_INIT=" << VX_INIT << std::endl;
	std::cout << "VY_INIT=" << VY_INIT << std::endl;
	
	return 0;
}



void fix_xy(double* x, double* y, int* x_pixel, int* y_pixel)
{
	*x_pixel += (int)((*x)/mpp);
	*y_pixel += (int)((*y)/mpp);

	*y_pixel = Y_INIT*2 - (*y_pixel);
}

void calc_T(double theta, Matrix GAMMMA_matrix, Matrix* m)
{
	(*m).rotateMatrix(-1*(M_PI/2 - theta));
	(*m).productMatrix(GAMMMA_matrix.getMatrix());
	(*m).rotateMatrix(M_PI/2 - theta);
}

void update_velocity(Matrix *velocity, Matrix gravity, Matrix gammma)
{
	double v_ver = (*velocity).getElement(0, 0);
	double v_hor = (*velocity).getElement(1, 0);
	double grav_ver = gravity.getElement(0, 0);
	double grav_hor = gravity.getElement(1, 0);
	double gam_ver = abs(gammma.getElement(0, 0));
	double gam_hor = abs(gammma.getElement(1, 1));
	
	// vertical----------------------------
	v_ver = v_ver + (-1*gam_ver*v_ver + grav_ver)*dt/m;
	
	
	// horizontal--------------------------
	v_hor = v_hor + (-1*GAMMMA_ine*abs(v_hor)*v_hor - gam_hor*v_hor + grav_hor)*dt/m;
	
	// update velocity
	(*velocity).setMatrix({{v_ver}, {v_hor}});
}


void file_read(std::vector<int> *v, int data_num)
{
	std::ifstream ifs("../python_work/_data_0121.txt");
	std::string str, buf;
	std::vector<int> _v;
	int data_count = 1;
	
	if(ifs.fail()) {
		std::cerr << "failed to read file." << std::endl;
		exit(1);
	}
	
	while(getline(ifs, str)) {
		std::stringstream _str;
		if(str == "")  data_count++;
		_str << str;
		
		if(data_count == data_num) {
			while(getline(_str, buf, ',')) {
				_v.push_back(std::stoi(buf));
			}
		} else if(data_count > data_num) break;
	}

	// set parameter
	X_INIT = _v[0];
	Y_INIT = _v[1];
	VX_INIT = (_v[2] - _v[0])*60*mpp;
	VY_INIT = -1*(_v[3] - _v[1])*60*mpp;
	
	*v = _v;
}


void fit_line(std::vector<int> *v, int n)
{
	int sum_x = 0, sum_y = 0, sum_xx = 0, sum_xy = 0, x, y;
	int sign_vx = signbit(VX_INIT) ? -1 : 1;
	
	for(int i = 0; i < n*2; i += 2) {
		x = (*v)[i];
		y = (*v)[i + 1];
		
		sum_x += x;
		sum_y += y;
		sum_xy += x*y;
		sum_xx += x*x;
	}
	
	double a = (double)(sum_x*sum_y - n*sum_xy)/(sum_x*sum_x - n*sum_xx);
	
	VX_INIT = VX_INIT/mpp/60;
	VY_INIT = VY_INIT/mpp/60;
	double v_size = sqrt(VX_INIT*VX_INIT + VY_INIT*VY_INIT);
	
	a *= -1;
	
	// set (vx, vy) [pixel/(1/60)s]
	VX_INIT = v_size/sqrt(1 + a*a)*sign_vx;
	VY_INIT = a*VX_INIT;

	// convert -> [m/s]
	VX_INIT = VX_INIT*60*mpp;
	VY_INIT = VY_INIT*60*mpp;

	VX_INIT *= v_mag_x;
	VY_INIT *= v_mag_y;


	// invert shuttle direction if vx < 0
	if(VX_INIT < 0) {
		for(int i = 0; i < (*v).size(); i += 2) {
			(*v)[i] = x_max - (*v)[i];
		}
		
		VX_INIT *= -1;
		X_INIT = (*v)[0];
	}
}



void scale(std::vector<int> *v, double theta)
{
	int x, y, prev_x = (*v)[0], prev_y = (*v)[1];
	double dx;
	std::vector<int> delta_v;

	for(int i = 2; i < (*v).size(); i = i + 2) {
		x = (*v)[i];
		y = (*v)[i + 1];

		dx = (double)(x - prev_x);
		dx /= cos(theta);

		delta_v.push_back(std::round(dx - (x - prev_x)));

		prev_x = x;
		prev_y = y;
	}
	
	
	int delta_sum = 0;
	
	for(int i = 2; i < (*v).size(); i = i + 2) {
		delta_sum += delta_v[i/2 - 1];
		(*v)[i] += delta_sum;
	}
}










