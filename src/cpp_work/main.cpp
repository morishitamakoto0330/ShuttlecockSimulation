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


void extractFrame(cv::Mat* frame, cv::Mat* im);

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
	int x, y;
	int x_lower = 670;
	int x_upper = 1350;
	
	std::string input_win = "input";
	std::string output_win = "output";
	std::string label_win = "label";

	cv::Mat im1, im2, im3, frame;
	cv::Mat im_mask;
	cv::Mat labeledImage;

	mouseParam mouseEvent;
	
	cv::setMouseCallback(input_win, CallBackFunc, &mouseEvent);

	// create window
	cv::namedWindow(input_win);
	cv::namedWindow(output_win);
	cv::namedWindow(label_win);

	while(1) {

		// get input key
		int key = cv::waitKey(20);

		// finish processing
		if(key == 27) {
			cv::destroyAllWindows();
			break;
		}
		
		// flame difference
		// get 3 flame
		cap >> frame;
		extractFrame(&frame, &im1);
		cap >> frame;
		extractFrame(&frame, &im2);
		cap >> frame;
		extractFrame(&frame, &im3);

		moveObjDetection(im1, im2, im3, &im_mask);

		// label image
		labeling(&im_mask, &labeledImage);

		// test
		cv::Mat test_img = labeledImage.clone();
		cv::line(test_img, cv::Point(x_lower,0), cv::Point(x_lower,test_img.rows-1), cv::Scalar(0,255,0), 1, cv::LINE_AA, 0);
		cv::line(test_img, cv::Point(x_upper,0), cv::Point(x_upper,test_img.rows-1), cv::Scalar(0,255,0), 1, cv::LINE_AA, 0);
		std::string test_str = "../../res/image/test_";
		test_str += std::to_string(test_index) +  ".png";
		cv::imwrite(test_str, test_img);

		test_index++;

		// show input and output
		cv::resize(frame, frame, cv::Size(), 0.6, 0.6);
		cv::imshow(input_win, frame);
		cv::resize(im_mask, im_mask, cv::Size(), 0.6, 0.6);
		cv::imshow(output_win, im_mask);
		cv::resize(labeledImage, labeledImage, cv::Size(), 0.6, 0.6);
		cv::imshow(label_win, labeledImage);

		// save image if put "s" key
		if(key == 115) {
			std::string img_name = "../../res/image/capture_" + std::to_string(index) + ".png";
			cv::imwrite(img_name, frame);

			std::cout << "succeed in saving capture." << std::endl;
			
			index++;
		}

		/*
		// get mouse position
		x = mouseEvent.x;
		y = mouseEvent.y;

		// print (x,y) (B,G,R) (H,S,V)
		dispPixelValue(frame, x, y);
		*/

		// shift 3 frames
		im2.copyTo(im1, im2);
		im3.copyTo(im2, im3);
		cap >> frame;
		cv::cvtColor(frame, im3, cv::COLOR_BGR2GRAY);
	}
	return 0;
}

void extractFrame(cv::Mat* frame, cv::Mat* im)
{
	cv::Mat extracted;
	colorExtraction(frame, &extracted, cv::COLOR_BGR2HSV, 0, 360, 0, 255, 0, 255);
	cvtColor(extracted, *im, cv::COLOR_BGR2GRAY);
}






