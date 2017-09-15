#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/core.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/highgui.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/imgproc.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/opencv.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/imgproc/types_c.h"

#include <iostream>
#include <string>
#include <vector>

#include "./extraction.hpp"
#include "./debug.hpp"

int main(void)
{
	int x, y, _x, _y;
	int old_x, old_y;
	int key;
	int H, S, V;
	int range;
	cv::Mat img;
	std::string win_name = "test";

	mouseParam mouseEvent;

	cv::namedWindow(win_name);
	
	cv::setMouseCallback(win_name, CallBackFunc, &mouseEvent);

	cv::createTrackbar("H", win_name, &H, 360, 0, 0);
	cv::createTrackbar("S", win_name, &S, 255, 0, 0);
	cv::createTrackbar("V", win_name, &V, 255, 0, 0);
	cv::createTrackbar("range", win_name, &range, 255, 0, 0);

	while(1)
	{
		key = cv::waitKey(100);
		
		// get mouse position
		x = mouseEvent.x;
		y = mouseEvent.y;

		// reset image
		//img = cv::imread("../../res/image/capture_1.png");
		img = cv::imread("../../res/image/capture_5.png");

		// extract color
		//colorExtraction(&img, &img, cv::COLOR_BGR2HSV, H-range, H+range, S-range, S+range, V-range, V+range);
		//colorExtraction(&img, &img, cv::COLOR_BGR2HSV, H-range, H+range, 0, 255, 0, 255);
		colorExtraction(&img, &img, cv::COLOR_BGR2HSV, 0, 360, S-range, S+range, 0, 255);
		//colorExtraction(&img, &img, cv::COLOR_BGR2HSV, 0, 0, 0, 0, V-range, V+range);

		_x = img.cols;
		_y = img.rows;

		// draw mousePos lines
		if((0 <= x) && (x < _x) && (0 <= y) && (y < _y))
		{
			cv::line(img, cv::Point(0,y), cv::Point(_x,y), cv::Scalar(0,255,0), 1, cv::LINE_AA, 0);
			cv::line(img, cv::Point(x,0), cv::Point(x,_y), cv::Scalar(0,255,0), 1, cv::LINE_AA, 0);
		}


		// show image
		cv::imshow(win_name, img);

		// print BGR and HSV value of click position
		if(mouseEvent.event == cv::EVENT_LBUTTONDOWN)
		{
			if((old_x != x) &&(old_y != y))
			{
				dispPixelValue(img, x, y);
			}
			old_x = x;
			old_y = y;
		}

		// print present current HSV and range value if put s key
		if(key == 115)
		{
			std::cout << "Trackbar:(H,S,V)=(" << H << "," << S << "," << V << ")";
			std::cout << "   range=" << range << std::endl;
		}

		// finish if put esc key
		if(key == 27) break;
	}
	return 0;
}




