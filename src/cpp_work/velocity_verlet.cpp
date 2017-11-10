/*
 * Oblique projection by using Velocity Verlet algorithm
 */

#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/opencv.hpp"

#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>

// data1
//static const int X_INIT = 534;
//static const int Y_INIT = 463;
// data3
static const int X_INIT = 956;
static const int Y_INIT = 582;


void fix_xy(double* x, double* y, int* x_pixel, int* y_pixel);

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
	
	int step, x_max, y_max;
	int x1, y1, x2, y2;
	double x, y, vx, vy, x_prev, y_prev, vx_prev, vy_prev;
	double fx, fy;
	double theta,  omega, torque, theta_prev, omega_prev, torque_prev;
	double m, g, dt;
	double GAMMMA, C, D, density, area;
	double GAMMMA_ver, GAMMMA_hor, l, I, G_diff, GAMMMA_ratio;
	double c1, c2, c3, c4, constant_A, constant_B, constant_C, constant_D;
	double sum_error;

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
	sum_error = 0;
	x_max = img.cols;
	y_max = img.rows;
	m = 0.005;    // mass
	g = 9.80665;  // gravity
	//dt = 0.01;    // delta time
	dt = 1.0/60;
	l = 0.005;     // distance between center of gravity and working point of resistance force
	I = 0.00001;  // inertia
	
	//density = 1.293;
	//area = M_PI*0.0365*0.0365;

	// draw net line
	cv::line(img, cv::Point(1090,0), cv::Point(1090,img.rows), cv::Scalar(255,0,0), 2);

	// set ratio considering air anisotropy
	GAMMMA_ratio = 1.5;


	// simulate until the object gets out of (x,y) range
	for(int i = 1; i <= 1; i++) {
		// initialize parameter

		//GAMMMA = i*0.001;
		//C = i*0.1;
		//D = C*area*density;
		//GAMMMA_hor = i*0.001;
		//GAMMMA_ver = i*0.001;
		GAMMMA_ver = 0.0057;
		//GAMMMA_ver = GAMMMA_hor*GAMMMA_ratio;
		GAMMMA_hor = GAMMMA_ver*GAMMMA_ratio;
		G_diff = GAMMMA_ver - GAMMMA_hor;
		
		
		step = 0;
		
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

		theta_prev = theta;
		omega_prev = omega;
		torque = (-1)*GAMMMA_ver*l*(vx*sin(theta) - vy*cos(theta));
		torque_prev = torque;

		while((0 <= x2) && (x2 < x_max) && (0 <= y2) && (y2 < y_max)) {
			
			step++;
			if(step*2 > v.size()) break;
			
			// update theta(t)--------------------------------------------------------------
			theta = theta_prev + omega_prev*dt + torque_prev*dt*dt/(2*I);
			
			
			// update vx(t), vy(t)---------------------------------------------------------
			
			// viscous resistance
			//vx = (2*m - GAMMMA*dt)/(2*m + GAMMMA*dt)*_vx;
			//vy = (2*m - GAMMMA*dt)/(2*m + GAMMMA*dt)*_vy - (2*m*g*dt/(2*m + GAMMMA*dt));
			
			// inertial resistance
			//vx = (sqrt(m*m + 2*m*D*dt*vx - (D*D*dt*dt*vx*vx)) - m)/(D*dt);
			//vy = (sqrt(m*m + 2*m*D*dt*vy - (D*D*dt*dt*vy*vy) - 2*m*g*D*dt*dt) - m)/(D*dt);
			
			// resistance considering rotation
			/*
			constant_E = 2*m + dt*(GAMMMA_hor*sin(theta)*sin(theta) + GAMMMA_ver*cos(theta)*cos(theta));
			constant_D = constant_E*(2*m + dt*(GAMMMA_hor*cos(theta)*cos(theta) + GAMMMA_ver*sin(theta)*sin(theta)));
			constant_C = dt*(G_diff)*sin(theta_prev)*cos(theta_prev)/constant_E + dt*G_diff*(2*m - dt*(GAMMMA_hor*cos(theta_prev)*cos(theta_prev) + GAMMMA_ver*sin(theta_prev)*sin(theta_prev)))*sin(theta)*cos(theta)/constant_D;
			constant_B = (2*m - dt*(GAMMMA_hor*sin(theta_prev)*sin(theta_prev) + GAMMMA_ver*cos(theta_prev)*cos(theta_prev)))/constant_E + G_diff*G_diff*dt*dt*sin(theta_prev)*cos(theta_prev)*sin(theta)*cos(theta)/constant_D;
			constant_A = 1 - G_diff*G_diff*dt*dt*sin(theta)*sin(theta)*cos(theta)*cos(theta)/constant_D;
			
			vx = constant_B*vx_prev/constant_A + constant_C*vy_prev/constant_A - G_diff*g*dt*dt/constant_A;
			
			constant_A = 2*m + dt*(GAMMMA_hor*cos(theta)*cos(theta) + GAMMMA_ver*sin(theta)*sin(theta));
			vy = (2*m - dt*(GAMMMA_hor*cos(theta_prev)*cos(theta_prev) + GAMMMA_ver*sin(theta_prev)*sin(theta_prev)))/constant_A*vy_prev + dt*G_diff*(vx_prev*sin(theta_prev)*cos(theta_prev) + vx*sin(theta)*cos(theta))/constant_A - 2*m*g*dt/constant_A;
			
			*/
			c1 = (2*m - dt*(GAMMMA_ver*sin(theta_prev)*sin(theta_prev) + GAMMMA_hor*cos(theta_prev)*cos(theta_prev)))/(2*m + dt*(GAMMMA_ver*sin(theta)*sin(theta) + GAMMMA_hor*cos(theta)*cos(theta)));
			c2 = dt*G_diff/(2*m + dt*(GAMMMA_ver*sin(theta)*sin(theta) + GAMMMA_hor*cos(theta)*cos(theta)));
			
			c3 = (2*m - dt*(GAMMMA_ver*cos(theta_prev)*cos(theta_prev) + GAMMMA_hor*sin(theta_prev)*sin(theta_prev)))/(2*m + dt*(GAMMMA_ver*cos(theta)*cos(theta) + GAMMMA_hor*sin(theta)*sin(theta)));
			c4 = dt*G_diff/(2*m + dt*(GAMMMA_ver*cos(theta)*cos(theta) + GAMMMA_hor*sin(theta)*sin(theta)));
			
			constant_A = 1 - c2*c4*sin(theta)*sin(theta)*cos(theta)*cos(theta);
			constant_B = c1 + c2*c4*sin(theta_prev)*sin(theta)*cos(theta_prev)*cos(theta);
			constant_C = c2*(sin(theta_prev)*cos(theta_prev) + c3*sin(theta)*cos(theta));
			constant_D = -c2*g*dt*sin(theta)*cos(theta);
			
			vx = (constant_B*vx_prev + constant_C*vy_prev + D)/constant_A;
			vy = c3*vy_prev + c4*(vx_prev*sin(theta_prev)*cos(theta_prev) + vx*sin(theta)*cos(theta)) - g*dt;
			
			// update torque(t)--------------------------------------------------------------
			torque = (-1)*GAMMMA_ver*l*(vx*sin(theta) - vy*cos(theta));
			
			
			// update omega(t)---------------------------------------------------------------
			omega = omega_prev + dt*(torque + torque_prev)/(2*I);
			
			
			// set f(t)--------------------------------------------------------------------
			// viscous resistance
			//fx = -1*GAMMMA*_vx;
			//fy = -1*GAMMMA*_vy - m*g;
			// resistance considering rotation
			fx = G_diff*vy_prev*sin(theta_prev)*cos(theta_prev) - (GAMMMA_ver*sin(theta_prev)*sin(theta_prev) + GAMMMA_hor*cos(theta_prev)*cos(theta_prev))*vx_prev;
			fy = G_diff*vx_prev*sin(theta_prev)*cos(theta_prev) - (GAMMMA_ver*cos(theta_prev)*cos(theta_prev) + GAMMMA_hor*sin(theta_prev)*sin(theta_prev))*vy_prev - m*g;
			
			// update x(t), y(t) ([m] -> [pixel])---------------------------------------------
			
			// viscous resistance
			//x = x + ((dt - GAMMMA/(2*m)*dt*dt)*_vx)*1000/4.31;
			//y = y + ((dt - (GAMMMA + m*g)/(2*m)*dt*dt)*_vy)*1000/4.31;
			
			// inertial resistance
			//x = x + ((1 - D*vx*dt/(2*m))*vx*dt)*1000/4.31;
			//y = y + ((1 - D*vy*dt/(2*m))*vy*dt - g*dt*dt/2)*1000/4.31;

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

			// draw shuttle trajectory
			// evaluate error
			//cv::line(img, cv::Point(x1,y1), cv::Point(x2,y2), cv::Scalar(i*25,0,255));
			int _x = v[step*2];
			int _y = v[step*2+1];
			sum_error += sqrt(pow(x1-_x, 2) + pow(y1-_y, 2));
			cv::circle(img, cv::Point(x1, y1), 4, cv::Scalar(255,0,0), 2);
			cv::circle(img, cv::Point(_x, _y), 4, cv::Scalar(0,0,255), 2);

			// next step
			x_prev = x;
			y_prev = y;
			vx_prev = vx;
			vy_prev = vy;
			theta_prev = theta;
			omega_prev = omega;
			torque_prev = torque;

		}
	}

	sum_error = sum_error*4.31/1000;
	std::cout << "error:" << sum_error << std::endl;
	
	cv::circle(img, cv::Point(X_INIT,Y_INIT), 4, cv::Scalar(0,0,0), 3);
	cv::imwrite("../../res/image_simulate/shuttle_points/shuttle_point_11_09.png", img);

	return 0;
}



void fix_xy(double* x, double* y, int* x_pixel, int* y_pixel)
{
	*x_pixel += (int)((*x)*1000/4.31);
	*y_pixel += (int)((*y)*1000/4.31);

	*y_pixel = Y_INIT*2 - (*y_pixel);
}










