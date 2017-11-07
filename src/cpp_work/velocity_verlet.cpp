/*
 * Oblique projection by using Velocity Verlet algorithm
 */

#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/opencv.hpp"

#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>

static const int X_INIT = 956;
static const int Y_INIT = 582;


void fix_xy(double* x, double* y, int* x_pixel, int* y_pixel);

int main(void)
{
	int step, x_max, y_max;
	int x1, y1, x2, y2;
	double x, y, vx, vy, x_prev, y_prev, vx_prev, vy_prev;
	double fx, fy;
	double theta,  omega, torque, theta_prev, omega_prev, torque_prev;
	double m, g, dt;
	double GAMMMA, C, D, density, area;
	double GAMMMA_ver, GAMMMA_hor, l, I, G_diff;
	double c1, c2, c3, c4, constant_A, constant_B, constant_C, constant_D;

	// output image
	cv::Mat img = cv::imread("../../res/image_simulate/shuttle_trajectory/shuttle_point_3.png");

	/*
	// output file name
	std::string file_name = "data.txt";

	std::ofstream writing_file;
	writing_file.open(file_name, std::ios::out);
	*/

	// initialize
	x_max = img.cols;
	y_max = img.rows;
	m = 0.005;    // mass
	g = 9.80665;  // gravity
	dt = 0.01;    // delta time
	l = 0.005;     // distance between center of gravity and working point of resistance force
	I = 0.00005;  // inertia
	
	GAMMMA_ver = 0.007;
	GAMMMA_hor = 0.007;
	G_diff = GAMMMA_ver - GAMMMA_hor;

	//density = 1.293;
	//area = M_PI*0.0365*0.0365;

	cv::line(img, cv::Point(1090,0), cv::Point(1090,img.rows), cv::Scalar(255,0,0), 2);

	// simulate until the object gets out of (x,y) range
	for(int i = 1; i <= 10; i++) {
		// initialize parameter

		//GAMMMA = i*0.001;
		//C = i*0.1;
		//D = C*area*density;
		GAMMMA_ver = i*0.002;
		GAMMMA_hor = i*0.001;
		G_diff = GAMMMA_ver - GAMMMA_hor;
		
		
		step = 0;
		
		x = 0.0;
		y = 0.0;
		x1 = 0;
		y1 = 0;
		x2 = 0;
		y2 = 0;

		// n=7
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
		torque = GAMMMA_ver*l*(vx*sin(theta) - vy*cos(theta));
		torque_prev = torque;

		while((0 <= x2) && (x2 < x_max) && (0 <= y2) && (y2 < y_max)) {
			
			step++;
			
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
			torque = GAMMMA_ver*l*(vx*sin(theta) - vy*cos(theta));
			
			
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
	}

	cv::circle(img, cv::Point(X_INIT,Y_INIT), 4, cv::Scalar(255,0,0), 3);
	cv::imwrite("../../res/image_simulate/shuttle_points/shuttle_point_11_06_4.png", img);

	return 0;
}



void fix_xy(double* x, double* y, int* x_pixel, int* y_pixel)
{
	*x_pixel += (int)((*x)*1000/4.31);
	*y_pixel += (int)((*y)*1000/4.31);

	*y_pixel = Y_INIT*2 - (*y_pixel);
}










