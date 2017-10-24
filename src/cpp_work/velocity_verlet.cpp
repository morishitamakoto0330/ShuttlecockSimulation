/*
 * Oblique projection by using Velocity Verlet algorithm
 */

#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/opencv.hpp"

#include <iostream>
#include <math.h>
#include <string.h>

int main(void)
{
	int step, x_max, y_max, y_init;
	double x, y, vx, vy, x_prev, y_prev, vx_prev, vy_prev;
	double fx, fy;
	double _x, _y, _vx, _vy, _fx, _fy, _x_prev, _y_prev;
	double theta,  omega, torque, theta_prev, omega_prev, torque_prev;
	double m, g, dt;
	//double Fx, Fy;
	double GAMMMA, C, D, density, area;
	double GAMMMA_ver, GAMMMA_hor, l, I, G_diff;
	double constant_A, constant_B, constant_C, constant_D, constant_E;

	// result image
	cv::Mat img = cv::imread("../../res/image_simulate/shuttle_trajectory/shuttle_point_3.png");

	// initialize
	x_max = img.cols;
	y_max = img.rows;
	m = 0.005;    // mass
	g = 9.80665;  // gravity
	dt = 0.01;    // delta time
	l = 0.01;     // distance between center of gravity and working point of resistance force
	I = 0.00005;  // inertia
	y_init = 582; // initial y value

	//density = 1.293;
	//area = M_PI*0.0365*0.0365;

	cv::line(img, cv::Point(1090,0), cv::Point(1090,img.rows), cv::Scalar(255,0,0), 2);

	// simulate until the object gets out of (x,y) range
	for(int i = 1; i <= 10; i++) {
		GAMMMA = i*0.001;
		//C = i*0.1;
		//D = C*area*density;
		GAMMMA_ver = i*0.001;
		GAMMMA_hor = i*0.001;
		G_diff = GAMMMA_hor - GAMMMA_ver;
		
		step = 0;
		
		x = 956.0;
		y = 582.0;
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

		_x = x;
		_y = y;
		_vx = vx;
		_vy = vy;
		_fx = fx;
		_fy = fy;
		_x_prev = _x;
		_y_prev = _y;

		theta_prev = theta;
		omega_prev = omega;
		torque = GAMMMA_hor*l*(vx*sin(theta) - vy*cos(theta));
		torque_prev = torque;

		while((0 < x) && (x < x_max) && (0 < y) && (y < y_init*2)) {
			
			step++;
			
			// update theta(t)--------------------------------------------------------------
			theta = theta_prev + omega_prev*dt + torque_prev*dt*dt/(2*I);
			
			
			// update vx(t), vy(t)---------------------------------------------------------
			
			// viscous resistance
			_vx = (2*m - GAMMMA*dt)/(2*m + GAMMMA*dt)*_vx;
			_vy = (2*m - GAMMMA*dt)/(2*m + GAMMMA*dt)*_vy - (2*m*g*dt/(2*m + GAMMMA*dt));
			
			// inertial resistance
			//vx = (sqrt(m*m + 2*m*D*dt*vx - (D*D*dt*dt*vx*vx)) - m)/(D*dt);
			//vy = (sqrt(m*m + 2*m*D*dt*vy - (D*D*dt*dt*vy*vy) - 2*m*g*D*dt*dt) - m)/(D*dt);
			
			// resistance considering rotation
			constant_E = 2*m + dt*(GAMMMA_hor*sin(theta)*sin(theta) + GAMMMA_ver*cos(theta)*cos(theta));
			constant_D = constant_E*(2*m + dt*(GAMMMA_hor*cos(theta)*cos(theta) + GAMMMA_ver*sin(theta)*sin(theta)));
			constant_C = dt*(G_diff)*sin(theta_prev)*cos(theta_prev)/constant_E + dt*G_diff*(2*m - dt*(GAMMMA_hor*cos(theta_prev)*cos(theta_prev) + GAMMMA_ver*sin(theta_prev)*sin(theta_prev)))*sin(theta)*cos(theta)/constant_D;
			constant_B = (2*m - dt*(GAMMMA_hor*sin(theta_prev)*sin(theta_prev) + GAMMMA_ver*cos(theta_prev)*cos(theta_prev)))/constant_E + G_diff*G_diff*dt*dt*sin(theta_prev)*cos(theta_prev)*sin(theta)*cos(theta)/constant_D;
			constant_A = 1 - G_diff*G_diff*dt*dt*sin(theta)*sin(theta)*cos(theta)*cos(theta)/constant_D;
			
			vx = constant_B*vx_prev/constant_A + constant_C*vy_prev/constant_A - G_diff*g*dt*dt/constant_A;
			
			constant_A = 2*m + dt*(GAMMMA_hor*cos(theta)*cos(theta) + GAMMMA_ver*sin(theta)*sin(theta));
			vy = (2*m - dt*(GAMMMA_hor*cos(theta_prev)*cos(theta_prev) + GAMMMA_ver*sin(theta_prev)*sin(theta_prev)))/constant_A*vy_prev + dt*G_diff*(vx_prev*sin(theta_prev)*cos(theta_prev) + vx*sin(theta)*cos(theta))/constant_A - 2*m*g*dt/constant_A;
			
			
			// update torque(t)--------------------------------------------------------------
			torque = GAMMMA_hor*l*(vx*sin(theta) - vy*cos(theta));
			
			
			// update omega(t)---------------------------------------------------------------
			omega = omega_prev + dt*(torque + torque_prev)/(2*I);
			
			
			// set f(t)--------------------------------------------------------------------
			// viscous resistance
			_fx = -1*GAMMMA*_vx;
			_fy = -1*GAMMMA*_vy - m*g;
			// resistance considering rotation
			fx = G_diff*vy_prev*sin(theta_prev)*cos(theta_prev) - (GAMMMA_hor*sin(theta_prev)*sin(theta_prev) + GAMMMA_ver*cos(theta_prev)*cos(theta_prev))*vx_prev;
			fy = G_diff*vx_prev*sin(theta_prev)*cos(theta_prev) - (GAMMMA_hor*cos(theta_prev)*cos(theta_prev) + GAMMMA_ver*sin(theta_prev)*sin(theta_prev))*vy_prev;
			fy -= m*g;
			
			// update x(t), y(t) ([m] -> [pixel])---------------------------------------------
			
			// viscous resistance
			_x = _x + ((dt - GAMMMA/(2*m)*dt*dt)*_vx)*1000/4.31;
			_y = _y + ((dt - (GAMMMA + m*g)/(2*m)*dt*dt)*_vy)*1000/4.31;
			
			// inertial resistance
			//x = x + ((1 - D*vx*dt/(2*m))*vx*dt)*1000/4.31;
			//y = y + ((1 - D*vy*dt/(2*m))*vy*dt - g*dt*dt/2)*1000/4.31;

			// resistance considering rotation
			x = x_prev + (vx_prev*dt + fx*dt*dt/(2*m))*1000/4.31;
			y = y_prev + (vy_prev*dt + fy*dt*dt/(2*m))*1000/4.31;
			
			
			std::cout << "step:" << step << "-------------------------------" << std::endl;
			std::cout << "(x,y)=(" << x << "," << y_init*2-y << ")";
			std::cout << "(vx,vy)=(" << vx << "," << vy << ")" << std::endl;
			std::cout << "(fx,fy)=(" << fx << "," << fy << ")" << std::endl;
			std::cout << "(_fx,_fy)=(" << _fx << "," << _fy << ")" << std::endl;
			std::cout << "θ=" << theta << ",ω=" << omega << ",T=" << torque << std::endl;

			cv::line(img, cv::Point(x_prev,y_init*2 - y_prev), cv::Point(x,y_init*2 - y), cv::Scalar(i*25,0,255));
			//cv::line(img, cv::Point(_x_prev,y_init*2 - _y_prev), cv::Point(_x,y_init*2 - _y), cv::Scalar(0,i*25,255));

			// next step
			x_prev = x;
			y_prev = y;
			vx_prev = vx;
			vy_prev = vy;
			theta_prev = theta;
			omega_prev = omega;
			torque_prev = torque;

			_x_prev = _x;
			_y_prev = _y;
		}
	}

	cv::circle(img, cv::Point(956,y_init), 4, cv::Scalar(255,0,0), 3);
	cv::imwrite("../../res/image_simulate/shuttle_points/shuttle_point_10_24.png", img);

	return 0;
}



