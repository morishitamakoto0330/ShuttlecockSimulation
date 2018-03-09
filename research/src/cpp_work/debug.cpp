#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/core.hpp"
#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/highgui.hpp"
#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/imgproc.hpp"
#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/opencv.hpp"
#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/imgproc/types_c.h"

#include <iostream>
#include <string>
#include <vector>

#include "./debug.hpp"



void disp_grav_pos(std::vector<std::pair<int, int>> v)
{
	std::cout << "grav_pos" << std::endl;
	for(int i = 0; i < v.size(); i++) {
		std::cout << "index: " << i << ", (x,y)=(" << v[i].first << "," << v[i].second << ")" << std::endl;
	}
}

void disp_area_value(std::vector<int> v)
{
	std::cout << "area_value" << std::endl;
	for(int i = 0; i < v.size(); i++) {
		std::cout << "index: " << i << ", area=" << v[i] << std::endl;
	}
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

// write point on (x,y) in white color image(shuttle_point.png)
void writePoint(int x, int y, cv::Mat* img)
{
	

	cv::circle(*img, cv::Point(x,y), 4, cv::Scalar(0,0,255), 3, 4);

	/*
	cv::Mat img = cv::imread(img_name);
	cv::circle(img, cv::Point(x,y), 5, cv::Scalar(0,0,255), 5, 4);
	cv::line(img, cv::Point(670,0), cv::Point(670,img.rows-1), cv::Scalar(0,255,0), 1, cv::LINE_AA, 0);
	cv::line(img, cv::Point(1350,0), cv::Point(1350,img.rows-1), cv::Scalar(0,255,0), 1, cv::LINE_AA, 0);
	cv::imwrite(img_name, img);
	*/
}


// write string
void write_string_to_file(std::string str, std::string file_name)
{
	// output stream
	std::ofstream ofs(file_name, std::ios::out);
	
	// output
	ofs << str << std::endl;
}










