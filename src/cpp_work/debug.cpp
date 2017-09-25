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

// write point on (x,y) in certain image file
void writePoint(int x, int y)
{
	std::string img_name = "../../res/image/result/shuttle_point.png";
	
	cv::Mat img = cv::imread(img_name);
	cv::circle(img, cv::Point(x,y), 5, cv::Scalar(0,0,255), 5, 4);
	cv::line(img, cv::Point(670,0), cv::Point(670,img.rows-1), cv::Scalar(0,255,0), 1, cv::LINE_AA, 0);
	cv::line(img, cv::Point(1350,0), cv::Point(1350,img.rows-1), cv::Scalar(0,255,0), 1, cv::LINE_AA, 0);
	cv::imwrite(img_name, img);
}





































