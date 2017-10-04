#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/core.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/highgui.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/imgproc.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/opencv.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/imgproc/types_c.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "./extraction.hpp"
#include "./debug.hpp"

int main(void)
{

	std::map<int, int> m;

	m.insert(std::make_pair(0, 0));
	/*
	std::vector<std::vector<int>> v;
	std::vector<int> _v(10, 5);

	
	for(int i = 0; i != _v.size(); i++)
	{
		v.push_back(_v);
	}

	for(int i = 0; i != v.size(); i++)
	{
		for(int j = 0; j != v[i].size(); j++)
		{
			std::cout << i << "," << j << ":" << v[i][j] << " ";
		}
		std::cout << std::endl;
	}
	*/


	/*
	int key;

	cv::Mat frame, img;
	cv::Mat img_back = cv::imread("../../res/image_interlace/capture_36400.png");
	
	cv::VideoCapture cap("../../res/movie/progressive_1.MTS");

	if(!cap.isOpened()) return -1;
	std::cout << "MTS file open." << std::endl;

	
	while(1)
	{
		key = cv::waitKey(100);
		
		cap >> frame;
		img = frame.clone();

		// show image
		cv::resize(img, img, cv::Size(), 0.6, 0.6);
		cv::imshow(output, img);
		
		if(key == 27) break;
	}
	*/
	
	return 0;
}




