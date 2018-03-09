#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/opencv.hpp"

#include <iostream>
#include <string>
#include <vector>

#include "./extraction.hpp"
#include "./debug.hpp"


// global variable
extern int area_threshold;
extern int x_threshold;

extern std::string img_white_name;
extern std::vector<std::pair<int, int>> shuttle_trajectory;
extern std::vector<int> shuttle_area_value;
extern std::vector<std::vector<std::pair<int, int>>> prev_gravity_pos;
extern std::vector<std::vector<int>> prev_area_value;
extern std::vector<int> track_index;


int main(int argc, char* argv[])
{
	// open file
	//cv::VideoCapture cap("../../res/movie/progressive_1.mp4");
	//cv::VideoCapture cap("../../res/movie/2017_01_18.mp4");
	cv::VideoCapture cap("../../res/movie/00007.mp4");
	
	// check file open
	if(!cap.isOpened()) {
		std::cout << "Failed to open movie file." << std::endl;
		return 1;
	}

	// valiable
	int x, y, area;
	int index = 1;
	int test_index = 1;
	double frame_count = 0, prop_frame_count;
	
	// get frame num
	prop_frame_count = cap.get(cv::CAP_PROP_FRAME_COUNT);

	
	std::string input_win = "input";
	std::string output_win = "output";
	std::string label_win = "label";

	cv::Mat im1, im2, im3, frame;
	cv::Mat im_mask1, im_mask2, im_mask3;
	cv::Mat im_mask, im_mask_old, im_mask_init;
	cv::Mat labeledImage;


	mouseParam mouseEvent;



	// skip frames
	while(1) {
		cap >> frame;
		frame_count++;
		
		if(frame_count > 120) break;
	}
	std::cout << "skip frames : " << frame_count << std::endl;




	// set mouse event
	cv::setMouseCallback(input_win, CallBackFunc, &mouseEvent);

	// create window
	cv::namedWindow(input_win);
	cv::namedWindow(output_win);
	cv::namedWindow(label_win);

	// get initial 3 flames
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
		if(key == 27) break;

		/*
		// finish processing if put ESC
		if(key == 27) {
			int n = shuttle_trajectory.size();
			// erase last data
			shuttle_trajectory.erase(shuttle_trajectory.begin() + n - 1);
			shuttle_area_value.erase(shuttle_area_value.begin() + n - 1);
			// insert head data
			shuttle_trajectory.insert(shuttle_trajectory.begin(), std::make_pair(-1, -1));
			shuttle_area_value.insert(shuttle_area_value.begin(), -1);
			
			int target_index = n - 1;
			int erase_num_front = 0, erase_num_back = 0;
			bool isFront = true;
			int b, g, r;
			double angle;

			b = 0;
			g = 0;
			r = 255;
			int count = 0;
			cv::Mat img_white = cv::imread(img_white_name);
			cv::line(img_white, cv::Point(x_threshold, 0), cv::Point(x_threshold, img_white.rows), cv::Scalar(255,0,0), 2);
			
			
			for(int i = n - 1; i >= 0; i--) {
				x = shuttle_trajectory[i].first;
				y = shuttle_trajectory[i].second;
				area = shuttle_area_value[i];
				
				
				if(x == -1) {
					b = rand()&255;
					g = rand()&255;
					r = rand()&255;
					
					// erase not shuttle_trajectory point
					// front
					for(int j = 0; j < erase_num_front; j++) {
						shuttle_trajectory.erase(shuttle_trajectory.begin() + target_index);
						shuttle_area_value.erase(shuttle_area_value.begin() + target_index);
						target_index--;
					}
					// back
					for(int j = 0; j < erase_num_back; j++) {
						shuttle_trajectory.erase(shuttle_trajectory.begin() + i + 1);
						shuttle_area_value.erase(shuttle_area_value.begin() + i + 1);
						target_index--;
					}

					// draw shuttle_trajectory
					for(int j = target_index; j > i; j--) {
						cv::circle(img_white, cv::Point(shuttle_trajectory[j].first, shuttle_trajectory[j].second), 4, cv::Scalar(b,g,r), 2);
					}
					
					// update valiabl to erase noise
					target_index = i - 1;
					erase_num_front = 0;
					erase_num_back = 0;
					isFront = true;
					
					continue;
				}

				
				// calc angle 
				if(((i-2) >= 0) && (x != -1) && (shuttle_trajectory[i-1].first != -1) && (shuttle_trajectory[i-2].first != -1)) {
					calc_angle(shuttle_trajectory[i], shuttle_trajectory[i-1], shuttle_trajectory[i-2], &angle);
					x = shuttle_trajectory[i-1].first;
					y = shuttle_trajectory[i-1].second;
					cv::putText(img_white, std::to_string((int)angle), cv::Point(x+10, y+10), 0, 0.5, cv::Scalar(0,255,0), 2);
					// check erase index
					if(isFront && (angle >= 14)) {
						erase_num_front++;
					} else if(isFront && (angle < 14)) {
						isFront = false;
					} else if(!isFront && (angle >= 14)) {
						erase_num_back++;
					}
				}
				
			}

			// reverse
			shuttle_trajectory.erase(shuttle_trajectory.begin());
			shuttle_area_value.erase(shuttle_area_value.begin());
			n = shuttle_trajectory.size();
			shuttle_trajectory.insert(shuttle_trajectory.begin() + n, std::make_pair(-1, -1));
			shuttle_area_value.insert(shuttle_area_value.begin() + n, -1);


			// erase noise (the number of points <= 20)-----------------------------
			count = -1;
			for(int i = shuttle_trajectory.size() - 1; i >= 0; i--) {
				x = shuttle_trajectory[i].first;

				if(x != -1) count++;
				else {
					if(count <= 20) {
						for(int j = 0; j < count + 1; j++) {
							shuttle_trajectory.erase(shuttle_trajectory.begin() + i + 1);
							shuttle_area_value.erase(shuttle_area_value.begin() + i + 1);
						}
					}
					count = 0;
				}
			}
			if(count <= 20) {
				for(int j = 0; j < count + 1; j++) {
					shuttle_trajectory.erase(shuttle_trajectory.begin());
					shuttle_area_value.erase(shuttle_area_value.begin());
				}
			}
			// -------------------------------------------------------------------


			// disp shuttle_trajectory
			count = 1;
			std::cout << "shuttle [No." << count << "]------" << std::endl;
			
			for(int i = 0; i < shuttle_trajectory.size(); i++) {
				x = shuttle_trajectory[i].first;
				y = shuttle_trajectory[i].second;
				area = shuttle_area_value[i];
				
				std::cout << "(x,y)=(" << x << "," << y << "), area=" << area;
				std::cout << std::endl;
				
				if(x == -1) {
					count++;
					std::cout << "shuttle [No." << count << "]------" << std::endl;
				}
			}
			cv::imwrite("../../trash/image_hoge/hoge0117.png", img_white);
			
			cv::destroyAllWindows();
			break;
		}
	*/
		
		// move object detection from 3 frames
		moveObjDetection(im1, im2, im3, &im_mask);

		// remove noise by erode and dilate
		//erode_dilate(im_mask, &im_mask, 1);

		// label image
		labeling(im_mask, &labeledImage);

		cv::Mat a,b;

		create_image(im_mask, &a, im_mask);
		combine_image(a, labeledImage, &b);
		cv::resize(b, b, cv::Size(), 0.33, 0.33);
		cv::imshow("moving object <-     -> labeling", b);


		im_mask2.copyTo(im_mask1);
		im_mask3.copyTo(im_mask2);
		im_mask3 = im_mask.clone();

		// shift 3 frames
		im2.copyTo(im1, im2);
		im3.copyTo(im2, im3);
		cap >> frame;
		cv::cvtColor(frame, im3, cv::COLOR_BGR2GRAY);

		// check frame
		frame_count++;
		if((frame_count-5) == prop_frame_count) {
			std::cout << "finish movie file." << std::endl;
			break;
		}
	}






	std::cout << "prop_frame_count=" << prop_frame_count << std::endl;
	std::cout << "frame_count=" << frame_count << std::endl;







			int n = shuttle_trajectory.size();
			// erase last data
			shuttle_trajectory.erase(shuttle_trajectory.begin() + n - 1);
			shuttle_area_value.erase(shuttle_area_value.begin() + n - 1);
			// insert head data
			shuttle_trajectory.insert(shuttle_trajectory.begin(), std::make_pair(-1, -1));
			shuttle_area_value.insert(shuttle_area_value.begin(), -1);
			
			int target_index = n - 1;
			int erase_num_front = 0, erase_num_back = 0;
			bool isFront = true;
			int b, g, r;
			double angle;

			b = 0;
			g = 0;
			r = 255;
			int count = 0;
			cv::Mat img_white = cv::imread(img_white_name);
			cv::line(img_white, cv::Point(x_threshold, 0), cv::Point(x_threshold, img_white.rows), cv::Scalar(255,0,0), 2);
			
			
			for(int i = n - 1; i >= 0; i--) {
				x = shuttle_trajectory[i].first;
				y = shuttle_trajectory[i].second;
				area = shuttle_area_value[i];
				
				
				if(x == -1) {
					b = rand()&255;
					g = rand()&255;
					r = rand()&255;
					
					// erase not shuttle_trajectory point
					// front
					for(int j = 0; j < erase_num_front; j++) {
						shuttle_trajectory.erase(shuttle_trajectory.begin() + target_index);
						shuttle_area_value.erase(shuttle_area_value.begin() + target_index);
						target_index--;
					}
					// back
					for(int j = 0; j < erase_num_back; j++) {
						shuttle_trajectory.erase(shuttle_trajectory.begin() + i + 1);
						shuttle_area_value.erase(shuttle_area_value.begin() + i + 1);
						target_index--;
					}

					// draw shuttle_trajectory
					for(int j = target_index; j > i; j--) {
						cv::circle(img_white, cv::Point(shuttle_trajectory[j].first, shuttle_trajectory[j].second), 4, cv::Scalar(b,g,r), 2);
					}
					
					// update valiabl to erase noise
					target_index = i - 1;
					erase_num_front = 0;
					erase_num_back = 0;
					isFront = true;
					
					continue;
				}

				
				// calc angle 
				if(((i-2) >= 0) && (x != -1) && (shuttle_trajectory[i-1].first != -1) && (shuttle_trajectory[i-2].first != -1)) {
					calc_angle(shuttle_trajectory[i], shuttle_trajectory[i-1], shuttle_trajectory[i-2], &angle);
					x = shuttle_trajectory[i-1].first;
					y = shuttle_trajectory[i-1].second;
					cv::putText(img_white, std::to_string((int)angle), cv::Point(x+10, y+10), 0, 0.5, cv::Scalar(0,255,0), 2);
					// check erase index
					if(isFront && (angle >= 14)) {
						erase_num_front++;
					} else if(isFront && (angle < 14)) {
						isFront = false;
					} else if(!isFront && (angle >= 14)) {
						erase_num_back++;
					}
				}
				
			}

			// reverse
			shuttle_trajectory.erase(shuttle_trajectory.begin());
			shuttle_area_value.erase(shuttle_area_value.begin());
			n = shuttle_trajectory.size();
			shuttle_trajectory.insert(shuttle_trajectory.begin() + n, std::make_pair(-1, -1));
			shuttle_area_value.insert(shuttle_area_value.begin() + n, -1);


			// erase noise (the number of points <= 20)-----------------------------
			count = -1;
			for(int i = shuttle_trajectory.size() - 1; i >= 0; i--) {
				x = shuttle_trajectory[i].first;

				if(x != -1) count++;
				else {
					if(count <= 10) {
						for(int j = 0; j < count + 1; j++) {
							shuttle_trajectory.erase(shuttle_trajectory.begin() + i + 1);
							shuttle_area_value.erase(shuttle_area_value.begin() + i + 1);
						}
					}
					count = 0;
				}
			}
			if(count <= 10) {
				for(int j = 0; j < count + 1; j++) {
					shuttle_trajectory.erase(shuttle_trajectory.begin());
					shuttle_area_value.erase(shuttle_area_value.begin());
				}
			}
			// -------------------------------------------------------------------


			// disp shuttle_trajectory
			count = 1;
			std::cout << "shuttle [No." << count << "]------" << std::endl;
			
			for(int i = 0; i < shuttle_trajectory.size(); i++) {
				x = shuttle_trajectory[i].first;
				y = shuttle_trajectory[i].second;
				area = shuttle_area_value[i];
				
				std::cout << "(x,y)=(" << x << "," << y << "), area=" << area;
				std::cout << std::endl;
				
				if(x == -1) {
					count++;
					std::cout << "shuttle [No." << count << "]------" << std::endl;
				}
			}
			cv::imwrite("../../trash/image_hoge/hoge0216.png", img_white);
			
			cv::destroyAllWindows();



	return 0;
}






