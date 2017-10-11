#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/core.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/highgui.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/imgproc.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/opencv.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/imgproc/types_c.h"

#include <iostream>
#include <string>
#include <vector>
#include <math.h>

#include "./extraction.hpp"
#include "./debug.hpp"

int main(void)
{

	int key;
	int x1 = 1090;
	int x2 = 850;

	std::vector<std::pair<int, int>> v;

	for(int i = 0; i < 10000; i++) {
		v.push_back(std::make_pair(0, 0));
	}

	/*
	std::string img_str = "../../res/image_progressive/shuttle_point/shuttle_point.png";

	cv::Mat frame, img;
	//cv::Mat img = cv::imread(img_str);
	
	cv::VideoCapture cap("../../res/movie/progressive_1.MTS");

	if(!cap.isOpened()) return -1;
	std::cout << "MTS file open." << std::endl;

	
	while(1) {
		key = cv::waitKey(100);
		
		cap >> frame;
		img = frame.clone();

		// show image
		cv::line(img, cv::Point(x1, 0), cv::Point(x1, img.rows), cv::Scalar(0, 0, 255), 3);
		//cv::line(img, cv::Point(x2, 0), cv::Point(x2, img.rows-1), cv::Scalar(0, 0, 255), 3);
		//cv::imwrite(img_str, img);
		cv::resize(img, img, cv::Size(), 0.6, 0.6);
		cv::imshow("img", img);
		
		if(key == 27) break;
	}
	*/
	return 0;
}




