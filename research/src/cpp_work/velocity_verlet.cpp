/*
 * Oblique projection by using Velocity Verlet algorithm
 */

#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/opencv.hpp"

#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>
//#include <sys/time.h>
//#include <sys/resource.h>

#include "./matrix.hpp"

// data1
//static const int X_INIT = 534;
//static const int Y_INIT = 463;
// data3
//static const int X_INIT = 956;
//static const int Y_INIT = 582;
int X_INIT, Y_INIT;
int DATA = 3;
double VX_INIT, VY_INIT;

const double m = 0.005;      // mass
const double g = 9.80665;    // gravity
const double dt = 1.0/600;   // delta time
// distance between center of gravity and working point of resistance force
const double l = 0.005;
const double I = 0.00001;    // inertia
// set ratio considering air anisotropy
const double ratio = 1.5;
const double _ratio = 0.1;


// adjustable parameter
double GAMMMA_ine, GAMMMA_ver, GAMMMA_hor;



void fix_xy(double* x, double* y, int* x_pixel, int* y_pixel);
void calc_T(double theta, Matrix GAMMMA_matrix, Matrix* m);
void update_velocity(Matrix* velocity, Matrix gravity, Matrix gammma);
void calc_equation(double* v, double gammma, double gravity, double sign);
void file_read(std::vector<int> *v, int data_num);
void fit_line(std::vector<int> v, int n);



int main(void)
{
	std::vector<int> v;
	file_read(&v, DATA);
	
	
	std::vector<double> error;

	int step, x_max, y_max;
	int x1, y1, x2, y2;
	double x, y, vx, vy, x_prev, y_prev, vx_prev, vy_prev;
	double fx, fy;
	double theta,  omega, torque, theta_prev, omega_prev, torque_prev;
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

	// output image
	cv::Mat img = cv::imread("../../res/white_image.png");

	// set parameter
	x_max = img.cols;
	y_max = img.rows;

	// draw net line
	cv::line(img, cv::Point(1090,0), cv::Point(1090,img.rows), cv::Scalar(255,0,0), 2);

	min_error = 0.0;

	// simulate until the object gets out of (x,y) range
	for(int i = 1; i <= 10; i++) {
		// initialize parameter
		GAMMMA_ver = i*0.001;
		//GAMMMA_ver = i*0.0001 + 0.008;
		GAMMMA_hor = GAMMMA_ver*ratio;
		GAMMMA_ine = GAMMMA_ver*_ratio;
		
		GAMMMA_matrix.setMatrix({{-1*GAMMMA_ver, 0}, {0, -1*GAMMMA_hor}});
		g_matrix.setMatrix({{0.0},{-1*dt*g}});
		
		step = 0;
		sum_error = 0.0;
		
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
		torque = (-1)*GAMMMA_ver*l*(vx*sin(theta) - vy*cos(theta));
		torque_prev = torque;

		while((0 <= x2) && (x2 < x_max) && (0 <= y2) && (y2 < y_max)) {
			
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
			torque = (-1)*GAMMMA_ver*l*(vx*sin(theta) - vy*cos(theta));
			
			
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

			// output theta value
			//writing_file << step << " " << theta/M_PI*180 << " " << atan(vy/vx)/M_PI*180 << std::endl;

			// calculate error
			if(step%10 == 0) {
				int _step = step/10;
				int _x = v[_step*2];
				int _y = v[_step*2+1];
				sum_error += sqrt(pow(x1-_x, 2) + pow(y1-_y, 2));
			
				// draw shuttle trajectory
				// simulate
				//cv::circle(img, cv::Point(x1, y1), 4, cv::Scalar(i*25,0,0), 2);
				//cv::line(img, cv::Point(x1,y1), cv::Point(x2,y2), cv::Scalar(i*25,0,255));
				// reality
				cv::circle(img, cv::Point(_x, _y), 4, cv::Scalar(0,0,255), 2);
			}
			cv::line(img, cv::Point(x1,y1), cv::Point(x2,y2), cv::Scalar(i*25,0,255));

			// next step
			x_prev = x;
			y_prev = y;
			vx_prev = vx;
			vy_prev = vy;
			theta_prev = theta;
			omega_prev = omega;
			torque_prev = torque;
		}
		
		// distance error in actual trajectory and simulation trajectory
		sum_error = sum_error*4.31/1000;
		error.push_back(sum_error);
		if(min_error == 0.0) min_error = sum_error;
		if(min_error > sum_error) {
			min_error = sum_error;
		}
		std::cout << "distance error:" << sum_error << std::endl;
	}
	
	
	std::cout << "index : distance error, average" << std::endl;
	
	for(int i = 0; i < error.size(); i++) {
		std::cout << i << ": " << error[i] << ", " << error[i]/(step/10) << std::endl;
	}
	
	cv::circle(img, cv::Point(X_INIT,Y_INIT), 4, cv::Scalar(0,0,0), 3);
	cv::imwrite("../../res/image_simulate/1215/shuttle_point_1.png", img);

	
	std::cout << std::endl;
	std::cout << "X_INIT=" << X_INIT << std::endl;
	std::cout << "Y_INIT=" << Y_INIT << std::endl;
	std::cout << "VX_INIT=" << VX_INIT << std::endl;
	std::cout << "VY_INIT=" << VY_INIT << std::endl;
	
	return 0;
}



void fix_xy(double* x, double* y, int* x_pixel, int* y_pixel)
{
	*x_pixel += (int)((*x)*1000/4.31);
	*y_pixel += (int)((*y)*1000/4.31);

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
	
	double sign_ver = 1;
	double sign_hor = 1;
	if(v_ver > 0) sign_ver = -1;
	if(v_hor > 0) sign_hor = -1;

	// vertical----------------------------
	//calc_equation(&v_ver, gam_ver, grav_ver, sign_ver);
	
	v_ver = v_ver + (-1*gam_ver*v_ver + grav_ver)*dt/m;
	
	
	// horizontal--------------------------
	//calc_equation(&v_hor, gam_hor, grav_hor, sign_hor);
	
	v_hor = v_hor + (-1*GAMMMA_ine*abs(v_hor)*v_hor - gam_hor*v_hor + grav_hor)*dt/m;
	
	
	// update velocity
	(*velocity).setMatrix({{v_ver}, {v_hor}});
}

void calc_equation(double* v, double gammma, double gravity, double sign)
{
	// viscous resistance--------------------------------
	/*
	(*v) = gravity/gammma + (_v - gravity/gammma)*exp(-1*gammma*dt/m);
	*/
	
	
	// inertial resistance------------------------------
	/*
	double _v = (*v);
	double C1, C2;
	if( ((_v > 0)&&(gravity < 0)) || ((_v < 0)&&(gravity > 0)) ) {
		gravity = abs(gravity);
		C1 = sqrt(gravity/gammma);
		
		(*v) = C1*tan(sign*sqrt(gammma*gravity)*dt/m + atan(_v/C1));
	} else {
		gravity = abs(gravity);
		C1 = sqrt(gravity/gammma);
		
		C2 = ((_v - C1)/(_v + C1))*exp(sign*sqrt(gammma*gravity)*dt/m);
		(*v) = C1*(1 + C2)/(1 - C2);
	}
	*/

	// both viscous and inertial resistance----------------
	/*
	double _v = (*v);
	double gammma_v = gammma;
	double gammma_i = gammma*_GAMMMA_ratio;
	double C = gravity/(sign*gammma_i) - gammma_v*gammma_v/(4*gammma_i*gammma_i);
	
	_v -= gammma_v/(2*sign*gammma_i);
	
	if(C > 0) {
		C = sqrt(C);
		
		(*v) = C*tan(atan(_v/C) + sign*gammma_i*C*dt/m);
	} else {
		C = sqrt(abs(C));
		
		double _C = ((_v - C)/(_v + C))*exp(2*sign*gammma_i*C*dt/m);
		(*v) = C*(1 + _C)/(1 - _C);
	}
	
	(*v) += gammma_v/(2*sign*gammma_i);
	*/

}

void file_read(std::vector<int> *v, int data_num)
{
	std::ifstream ifs("../python_work/_data_1215.txt");
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

	for(int i = 0; i < _v.size(); i++) {
		std::cout << _v[i] << std::endl;
	}
	
	// set parameter
	X_INIT = _v[0];
	Y_INIT = _v[1];
	VX_INIT = (_v[2] - _v[0])*60*4.31/1000;
	VY_INIT = -1*(_v[3] - _v[1])*60*4.31/1000;
	
	*v = _v;
}


void fit_line(std::vector<int> v, int n)
{
	int sum_x = 0, sum_y = 0, sum_xx = 0, sum_xy = 0, x, y;
	
	for(int i = 0; i < n*2; i += 2) {
		x = v[i];
		y = v[i + 1];
		
		sum_x += x;
		sum_y += y;
		sum_xy += x*y;
		sum_xx += x*x;
	}
	
	double a = (double)(sum_x*sum_y - n*sum_xy)/(sum_x*sum_x - n*sum_xx);
	double v_size = sqrt(VX_INIT*VX_INIT + VY_INIT*VY_INIT);
	
	a *= -1;
	
	// set (vx, vy)
	VX_INIT = v_size/sqrt(1 + a*a);
	VY_INIT = a*VX_INIT;
}

















