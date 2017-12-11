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
static const int X_INIT = 956;
static const int Y_INIT = 582;

const double m = 0.005;      // mass
const double g = 9.80665;    // gravity
const double dt = 1.0/600;   // delta time
// distance between center of gravity and working point of resistance force
const double l = 0.005;
const double I = 0.00001;    // inertia
// set ratio considering air anisotropy
const double GAMMMA_ratio = 1.0;
const double _GAMMMA_ratio = 1.0;

void fix_xy(double* x, double* y, int* x_pixel, int* y_pixel);
void calc_T(double theta, Matrix GAMMMA_matrix, Matrix* m);
void update_velocity(Matrix* velocity, Matrix gravity, Matrix gammma);
void calc_equation(double* v, double gammma, double gravity, double sign);



int main(void)
{
	// data1
	//std::vector<int> v{534, 463, 666, 442, 754, 437, 831, 441, 901, 444, 967, 444, 1027, 446, 1083, 449, 1136, 453, 1183, 457, 1232, 461, 1274, 464, 1316, 468, 1355, 473, 1394, 479, 1431, 484, 1464, 490, 1497, 495, 1529, 501, 1559, 509, 1589, 516, 1618, 523};
	// data2
	//std::vector<int> v{568, 547, 710, 494, 796, 477, 875, 471, 949, 461, 1012, 451, 1071, 444, 1126, 438, 1177, 433, 1224, 429, 1269, 426, 1311, 423, 1351, 421, 1388, 419, 1423, 418, 1456, 417, 1488, 418, 1519, 418, 1547, 419, 1576, 421, 1601, 423};
	// data3
	std::vector<int> v{956, 582, 987, 556, 1017, 529, 1041, 505, 1066, 480, 1092, 456, 1114, 436, 1137, 418, 1160, 399, 1179, 385, 1200, 370, 1221, 356, 1239, 343, 1257, 331, 1274, 320, 1289, 311, 1305, 301, 1324, 292, 1336, 285, 1351, 277, 1366, 270, 1380, 264, 1395, 258, 1409, 254, 1426, 250, 1439, 247, 1452, 244, 1464, 241, 1476, 242, 1485, 240, 1502, 239, 1512, 239, 1523, 239, 1537, 242, 1546, 242, 1561, 243, 1569, 245, 1581, 249, 1591, 252, 1606, 258, 1617, 262, 1627, 268, 1637, 272, 1647, 278, 1658, 285, 1668, 290, 1678, 297, 1687, 304, 1694, 310, 1707, 321, 1717, 330, 1726, 339, 1735, 349, 1744, 359, 1753, 370, 1761, 380, 1770, 392, 1779, 403};
	// data4
	//std::vector<int> v{838, 601, 888, 595, 938, 594, 985, 591, 1029, 586, 1066, 582, 1109, 580, 1141, 579, 1177, 578, 1208, 578, 1240, 577, 1267, 577, 1298, 577, 1324, 577, 1349, 578, 1373, 580, 1396, 582, 1417, 583, 1440, 586, 1460, 590, 1485, 593, 1502, 596, 1519, 597};
	
	std::vector<double> error;

	int step, x_max, y_max;
	int x1, y1, x2, y2;
	double x, y, vx, vy, x_prev, y_prev, vx_prev, vy_prev;
	double fx, fy;
	double theta,  omega, torque, theta_prev, omega_prev, torque_prev;
	double GAMMMA_ver, GAMMMA_hor;
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
	//cv::Mat img = cv::imread("../../res/image_simulate/shuttle_trajectory/shuttle_point_1.png");
	cv::Mat img = cv::imread("../../res/image_simulate/shuttle_trajectory/shuttle_point_3.png");

	/*
	// output file name
	std::string file_name = "data.txt";

	std::ofstream writing_file;
	writing_file.open(file_name, std::ios::out);
	*/

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
		GAMMMA_hor = GAMMMA_ver*GAMMMA_ratio;
		
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

		// data1 (n=2)
		//vx = 132*60*4.31/1000;
		//vy = 21*60*4.31/1000;
		// data3 (n=7)
		vx = 7.63772;
		vy = 7.15113;
		
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
			if(step/10*2 >= v.size()) break;
			
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
				cv::circle(img, cv::Point(x1, y1), 4, cv::Scalar(i*25,0,0), 2);
				//cv::line(img, cv::Point(x1,y1), cv::Point(x2,y2), cv::Scalar(i*25,0,255));
				// reality
				//cv::circle(img, cv::Point(_x, _y), 4, cv::Scalar(0,0,255), 2);
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
	cv::imwrite("../../res/image_simulate/1130/shuttle_point_6.png", img);

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
	calc_equation(&v_ver, gam_ver, grav_ver, sign_ver);
	
	// horizontal--------------------------
	calc_equation(&v_hor, gam_hor, grav_hor, sign_hor);
	
	
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
}










