#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/core.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/highgui.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/imgproc.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/opencv.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/imgproc/types_c.h"

#include <iostream>
#include <string>
#include <vector>

#include "./debug.hpp"


void dispString(std::string str)
{
	std::cout << str << std::endl;
}

void dispInt(int i)
{
	std::cout << i << std::endl;
}

// print BGR and HSV value
void dispPixelValue(cv::Mat frame, int x, int y)
{
	int B = frame.at<cv::Vec3b>(y, x)[0];
	int G = frame.at<cv::Vec3b>(y, x)[1];
	int R = frame.at<cv::Vec3b>(y, x)[2];

	cv::Mat hsv_frame;
	cv::cvtColor(frame, hsv_frame, cv::COLOR_BGR2HSV);

	int H = hsv_frame.at<cv::Vec3b>(y, x)[0];
	int S = hsv_frame.at<cv::Vec3b>(y, x)[1];
	int V = hsv_frame.at<cv::Vec3b>(y, x)[2];
	/*
	int H = hsv_frame.data[y*hsv_frame.step + x*hsv_frame.elemSize() + 0];
	int S = hsv_frame.data[y*hsv_frame.step + x*hsv_frame.elemSize() + 1];
	int V = hsv_frame.data[y*hsv_frame.step + x*hsv_frame.elemSize() + 2];
	*/
	std::cout << "(x,y)=" << "(" << x << "," << y << ")  ";
	std::cout << "(B,G,R)=" << "(" << B << "," << G << "," << R << ")  ";
	std::cout << "(H,S,V)=" << "(" << H << "," << S << "," << V << ")  ";
	std::cout << std::endl;
}




