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


int main(int argc, char* argv[])
{
	// open file
	//cv::VideoCapture cap("../../res/movie/interlace_1.MTS");
	cv::VideoCapture cap("../../res/movie/progressive_1.MTS");
	
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
	cv::Mat im_mask1, im_mask2, im_mask3;
	cv::Mat im_mask, im_mask_old, im_mask_init;
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
	cvtColor(frame, im1, cv::COLOR_BGR2GRAY);
	cap >> frame;
	cvtColor(frame, im2, cv::COLOR_BGR2GRAY);
	cap >> frame;
	cvtColor(frame, im3, cv::COLOR_BGR2GRAY);

	// init mask image
	im_mask1 = im1.clone();
	im_mask2 = im2.clone();
	im_mask3 = im3.clone();

	while(1) {
		// get input key
		int key = cv::waitKey(20);

		// finish processing if put ESC
		if(key == 27) {
			cv::destroyAllWindows();
			break;
		}
		
		// move object detection from 3 frames
		moveObjDetection(im1, im2, im3, &im_mask);

		// remove noise by erode and dilate
		erode_dilate(im_mask, &im_mask, 1);

		// label image
		labeling(im_mask, &labeledImage);

		/*
		// show input and output
		cv::resize(frame, frame, cv::Size(), 0.6, 0.6);
		cv::imshow(input_win, frame);
		cv::resize(im_mask, im_mask, cv::Size(), 0.6, 0.6);
		cv::imshow(output_win, im_mask);
		cv::resize(labeledImage, labeledImage, cv::Size(), 0.6, 0.6);
		cv::imshow(label_win, labeledImage);
		*/
		im_mask_init = cv::Mat::zeros(cv::Size(im1.cols, im1.rows), CV_8UC1);
		
		cv::bitwise_or(im_mask1, im_mask_init, im_mask_init);
		cv::bitwise_or(im_mask2, im_mask_init, im_mask_init);
		cv::bitwise_or(im_mask3, im_mask_init, im_mask_init);
		
		cv::resize(im_mask_init, im_mask_init, cv::Size(), 0.5, 0.5);
		cv::imshow("mask", im_mask_init);
		/*
		cv::Mat a,b;

		create_image(im_mask, &a, im_mask);
		combine_image(a, labeledImage, &b);
		cv::resize(b, b, cv::Size(), 0.33, 0.33);
		cv::imshow("moving object <-     -> labeling", b);
		*/

		/*
		// save capture
		//std::string img_name = "../../res/image_progressive/result/MOD_label_" + std::to_string(index) + ".png";
		std::string img_name = "../../res/image_progressive/result/MOD_label_erodeDilate_" + std::to_string(index) + ".png";
		index++;
		cv::imwrite(img_name, labeledImage);
		*/

		im_mask2.copyTo(im_mask1);
		im_mask3.copyTo(im_mask2);
		im_mask3 = im_mask.clone();

		// shift 3 frames
		im2.copyTo(im1, im2);
		im3.copyTo(im2, im3);
		cap >> frame;
		cv::cvtColor(frame, im3, cv::COLOR_BGR2GRAY);
	}
	return 0;
}






