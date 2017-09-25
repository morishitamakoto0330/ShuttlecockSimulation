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

static const int threshold = 3;


//void extractFrame(cv::Mat* frame, cv::Mat* im);
void myMask(cv::Mat* src, cv::Mat* dst);

int main(int argc, char* argv[])
{
	// open file
	cv::VideoCapture cap("../../res/movie/00000.MTS");
	
	// check file open
	if(!cap.isOpened()) {
		std::cout << "Failed to open movie file." << std::endl;
		return 1;
	}

	// valiable
	int index = 1;
	int test_index = 1;
	int x_lower = 670;
	int x_upper = 1350;
	
	std::string input_win = "input";
	std::string output_win = "output";
	std::string label_win = "label";

	cv::Mat im1, im2, im3, frame;
	cv::Mat im_mask;
	cv::Mat labeledImage;

	mouseParam mouseEvent;
	
	// set mouse event
	cv::setMouseCallback(input_win, CallBackFunc, &mouseEvent);

	// create window
	cv::namedWindow(input_win);
	cv::namedWindow(output_win);
	cv::namedWindow(label_win);

	// get 3 flame
	cap >> frame;
	myMask(&frame, &frame);
	cvtColor(frame, im1, cv::COLOR_BGR2GRAY);
	cap >> frame;
	myMask(&frame, &frame);
	cvtColor(frame, im2, cv::COLOR_BGR2GRAY);
	cap >> frame;
	myMask(&frame, &frame);
	cvtColor(frame, im3, cv::COLOR_BGR2GRAY);

	while(1) {

		// get input key
		int key = cv::waitKey(20);

		// finish processing
		if(key == 27) {
			cv::destroyAllWindows();
			break;
		}
		
		// move object detection from 3 frames
		moveObjDetection(im1, im2, im3, &im_mask);

		// label image
		labeling(&im_mask, &labeledImage);

		// show input and output
		cv::resize(frame, frame, cv::Size(), 0.6, 0.6);
		cv::imshow(input_win, frame);
		cv::resize(im_mask, im_mask, cv::Size(), 0.6, 0.6);
		cv::imshow(output_win, im_mask);
		cv::resize(labeledImage, labeledImage, cv::Size(), 0.6, 0.6);
		cv::imshow(label_win, labeledImage);

		
		// save
		std::string img_name = "../../res/image/result/capture_masked_deinterlace_" + std::to_string(index) + ".png";
		index++;
		cv::imwrite(img_name, labeledImage);





		// shift 3 frames
		im2.copyTo(im1, im2);
		im3.copyTo(im2, im3);
		cap >> frame;
		cv::cvtColor(frame, im3, cv::COLOR_BGR2GRAY);
	}
	return 0;
}

// background image mask
void myMask(cv::Mat* src, cv::Mat* dst)
{
	cv::Mat diff;
	cv::Mat img_back = cv::imread("../../res/image/capture_36400.png");
	
	deinterlace(&img_back, &img_back);
	deinterlace(src, src);

	cv::absdiff(*src, img_back, diff);
	cv::bitwise_and(*src, diff, *dst);
}




