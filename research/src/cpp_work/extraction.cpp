#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/core.hpp"
#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/highgui.hpp"
#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/imgproc.hpp"
#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/opencv.hpp"
#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/imgproc/types_c.h"

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <math.h>

#include "./extraction.hpp"


// global variable--------------------------------------------
int area_threshold = 50;
int x_threshold = 1090;
int y_threshold = 100;

// white image name to debug
std::string img_white_name = "../../res/white_image.png";
// shuttle detection
std::vector<std::pair<int, int>> shuttle_trajectory;
std::vector<int> shuttle_area_value;
// previous condition
std::vector<std::vector<std::pair<int, int>>> prev_gravity_pos;
std::vector<std::vector<int>> prev_area_value;
// track
std::vector<int> track_index;

/*
 * extract color------------------------------------------------
 */
void colorExtraction(cv::Mat* src, cv::Mat* dst,
		int code,
		int ch1Lower, int ch1Upper,
		int ch2Lower, int ch2Upper,
		int ch3Lower, int ch3Upper
)
{
	cv::Mat colorImage;
	int lower[3];
	int upper[3];

	cv::Mat lut = cv::Mat(256, 1, CV_8UC3);

	cv::cvtColor(*src, colorImage, code);

	lower[0] = ch1Lower;
	lower[1] = ch2Lower;
	lower[2] = ch3Lower;

	upper[0] = ch1Upper;
	upper[1] = ch2Upper;
	upper[2] = ch3Upper;

	// fix HSV value range
	if(lower[0] < 0) lower[0] = 0;
	if(lower[1] < 0) lower[1] = 0;
	if(lower[2] < 0) lower[2] = 0;

	// H: 0~360 -> 0~180
	upper[0] = upper[0] / 2;
	if(upper[0] > 180) upper[0] = 180;
	if(upper[1] > 255) upper[1] = 255;
	if(upper[2] > 255) upper[2] = 255;

	for(int i = 0; i < 256; i++) {
		for(int k = 0; k < 3; k++) {
			if(lower[k] <= upper[k]) {
				if((lower[k] <= i) && (i <= upper[k])) {
					lut.data[i*lut.step+k] = 255;
				} else {
					lut.data[i*lut.step+k] = 0;
				}
			} else {
				if((i <= upper[k]) || (lower[k] <= i)) {
					lut.data[i*lut.step+k] = 255;
				} else {
					lut.data[i*lut.step+k] = 0;
				}
			}
		}
	}

	// binarization by using LUT
	cv::LUT(colorImage, lut, colorImage);

	// disassemble per channel
	std::vector<cv::Mat> planes;
	cv::split(colorImage, planes);

	// crate mask
	cv::Mat maskImage;
	cv::bitwise_and(planes[0], planes[1], maskImage);
	cv::bitwise_and(maskImage, planes[2], maskImage);
	
	// output
	cv::Mat maskedImage;
	src->copyTo(maskedImage, maskImage);
	*dst = maskedImage.clone();
}






/* 
 * labeling--------------------------------------------------------
 */

void labeling(cv::Mat src, cv::Mat* dst)
{
	cv::Mat bin;
	cv::Mat stats;
	cv::Mat centroids;

	// binarization
	cv::threshold(src, bin, 0, 255, cv::THRESH_BINARY | cv::THRESH_OTSU);
	// create image to label
	cv::Mat labelImage(src.size(), CV_32S);

	// labeling(simple)
	//int labelNum = cv::connectedComponents(bin, labelImage, 8);
	// labeling(detail)
	int labelNum = cv::connectedComponentsWithStats(bin, labelImage, stats, centroids);
	
	// label coloring
	std::vector<cv::Vec3b> colors(labelNum);
	colors[0] = cv::Vec3b(0, 0, 0);
	for(int label = 1; label < labelNum; label++) {
		colors[label] = cv::Vec3b((rand()&255), (rand()&255), (rand()&255));
	}

	// create result image
	cv::Mat _dst(src.size(), CV_8UC3);
	for(int y = 0; y < _dst.rows; y++) {
		for(int x = 0; x < _dst.cols; x++) {
			int label = labelImage.at<int>(y, x);
			cv::Vec3b &pixel = _dst.at<cv::Vec3b>(y, x);
			pixel = colors[label];
		}
	}

	// get parameter(center of gravity, area value)
	//
	// gravity position (x,y) -> vector<pair<int, int>>
	// area value             -> vector<int>
	
	std::pair<int, int> tmp = std::make_pair(0, 0);
	std::vector<std::pair<int, int>> gravity_pos(labelNum, tmp);

	std::vector<int> area_value(labelNum, 0);

	
	int x, y, max_x, max_y;
	int area, max_area;

	for(int i = 1; i < labelNum; i++) {
		double *param = centroids.ptr<double>(i);

		x = static_cast<int>(param[0]);
		y = static_cast<int>(param[1]);

		gravity_pos[i-1] = std::make_pair(x, y);
	}

	for(int i = 1; i < labelNum; i++) {
		int *param = stats.ptr<int>(i);
		area_value[i-1] = param[cv::ConnectedComponentsTypes::CC_STAT_AREA];
	}


	max_area = 0;

	for(int i = area_value.size()-1; i >= 0; i--) {
		
		x = gravity_pos[i].first;
		y = gravity_pos[i].second;
		area = area_value[i];
		
		if(area >= area_threshold) {
			//cv::circle(img_white, cv::Point(x,y), 4, cv::Scalar(0,0,255), 2);
			//cv::imwrite(img_white_name, img_white);
			// get max area and its position(x,y)
			if(max_area < area) {
				max_area = area;
				max_x = x;
				max_y = y;
			}
		} else {
			area_value.erase(area_value.begin()+i);
			gravity_pos.erase(gravity_pos.begin()+i);
		}
	}
	
	
	
	
	// -------------------------------------------------------
	// check shuttle
	// no processing for the first time
	// After second time, check position and area change
	// -------------------------------------------------------

	x = 0;
	y = 0;
	area = 0;
	int _x, _y, _area, direction, s;
	_x = 0;
	_y = 0;
	_area = 0;

	// debug
	cv::Mat img_shuttle;

	// following tracking--------------------
	for(int i = 0; i < track_index.size(); i++) {
		// tracking point parameter
		s = shuttle_trajectory.size();
		x = shuttle_trajectory[track_index[i]-1].first;
		y = shuttle_trajectory[track_index[i]-1].second;
		area = shuttle_area_value[track_index[i]-1];
		//direction = signbit(x - shuttle_trajectory[track_index[i]-2].first);
		direction = signbit(y - shuttle_trajectory[track_index[i]-2].second);
		
		std::cout << "following track" << std::endl;
		std::cout << "(x,y)=(" << x << "," << y << "), area=" << area;
		std::cout << " direction=" << direction << std::endl;
		for(int j = 0; j < area_value.size(); j++) {
			// point candidate
			_x = gravity_pos[j].first;
			_y = gravity_pos[j].second;
			_area = area_value[j];
			std::cout << "(_x,_y)=(" << _x << "," << _y << "), _area=" << _area;
			std::cout << " direction=" << signbit(_x-x) << std::endl;

			//if((std::abs(x - _x) <= 150) && (std::abs(y - _y) <= 70) && ((std::abs(area - _area) <= 200) || (area>1000 && _area>1000))  && (signbit(_x - x) == direction)) {
			if((std::abs(x - _x) <= 10) && (std::abs(y - _y) <= 70) && ((std::abs(area - _area) <= 100) || (area>1000 && _area>1000))  && (signbit(_y - y) == direction)) {
				shuttle_trajectory.insert(shuttle_trajectory.begin() + track_index[i], std::make_pair(_x, _y));
				shuttle_area_value.insert(shuttle_area_value.begin() + track_index[i], _area);
				track_index[i]++;
				break;
			}
		}

		if(s == shuttle_trajectory.size()) {
			track_index.erase(track_index.begin() + i);
		}
	}
	
	// previous tracking----------------------------------

	if(prev_area_value.size() == 0) {
		//no processing...
	} else {
		int prev_index = prev_area_value.size() - 1;
		for(int i = 0; i < prev_area_value[prev_index].size(); i++) {
			// previous frame parameters
			x = prev_gravity_pos[prev_index][i].first;
			y = prev_gravity_pos[prev_index][i].second;
			area = prev_area_value[prev_index][i];
			
			for(int j = 0; j < area_value.size(); j++) {
				// prevent frame parameters
				_x = gravity_pos[j].first;
				_y = gravity_pos[j].second;
				_area = area_value[j];

				// check position and area change
				// and check cross net (x -> _x) or (_x -> x)
				//if((std::abs(x - _x) <= 150) && (std::abs(y - _y) <= 50) && (std::abs(area - _area) <= 50)) {
				if((std::abs(x - _x) <= 50) && (std::abs(y - _y) <= 70) && (std::abs(area - _area) <= 100)) {
					//if(((x<=x_threshold)&&(_x>x_threshold)) || ((x>x_threshold)&&(_x<=x_threshold))) {
					if(((y<=y_threshold)&&(_y>y_threshold)) || ((y>y_threshold)&&(_y<=y_threshold))) {
					//if(((y<=100)&&(_y>100)) || ((y>100)&&(_y<=100))) {
						std::cout << "cross net at ";
						std::cout << "(x,y)=(" << x << "," << y << "), ";
						std::cout << "(_x,_y)=(" << _x << "," << _y << ")" << std::endl;
						// (x,y)->(_x,_y) : cross net!
						// (x,y)   : start point (previous)
						// (_x,_y) : start point (following)
						
						// previous tracking (start at (x,y))
						// track shuttle prev_gravity_pos and prev_area_value
						// index: prev_index-1, prev_index-2, ..., 0
						// finish if no shuttle point
						int now_index = shuttle_trajectory.size();
						//direction = signbit(x - _x); // 1:L->R 0:R->L
						direction = signbit(y - _y);
						shuttle_trajectory.push_back(std::make_pair(x, y));
						shuttle_trajectory.push_back(std::make_pair(_x, _y));
						shuttle_area_value.push_back(area);
						shuttle_area_value.push_back(_area);

						// hoge1:target parameters , hoge2:check parameters
						int x1, x2, y1, y2, area1, area2, shuttle_size;
						x1 = x;
						y1 = y;
						area1 = area;
						shuttle_size = shuttle_trajectory.size();
						
						std::cout << "previous track" << std::endl;
						std::cout << "(x1,y1)=(" << x1 << "," << y1 << "), area1=" << area1 << std::endl;
						std::cout << "direction1=" << direction << std::endl;
						
						for(int k = prev_index-1; k >= 0; k--) {
							// check previous frame
							for(int l = 0; l < prev_area_value[k].size(); l++) {
								x2 = prev_gravity_pos[k][l].first;
								y2 = prev_gravity_pos[k][l].second;
								area2 = prev_area_value[k][l];
								std::cout << "(x2,y2)=(" << x2 << "," << y2 << "), area2=" << area2 << std::endl;
								std::cout << "direction2=" << direction << std::endl;
								// check around (x,y)
								//if((std::abs(x1 - x2) <= 150) && (std::abs(y1 - y2) <= 70) && ((std::abs(area1 - area2) <= 200) || ((area1>1000)&&(area2>1000))) && (signbit(x2 - x1) == direction)) {
								if((std::abs(x1 - x2) <= 50) && (std::abs(y1 - y2) <= 70) && ((std::abs(area1 - area2) <= 200) || ((area1>1000)&&(area2>1000))) && (signbit(y2 - y1) == direction)) {
									std::cout << "うーん、gooooooooooooooooooooooooood！" << std::endl;
									// insert (x2, y2) to head in shuttle_trajectory
									// insert area2 to head in shuttle_area_value
									shuttle_trajectory.insert(shuttle_trajectory.begin()+now_index, std::make_pair(x2, y2));
									shuttle_area_value.insert(shuttle_area_value.begin()+now_index, area2);
									break;
								}
							}
							// finish previous position tracking if no change of size
							if(shuttle_size == shuttle_trajectory.size()) break;
							
							// next
							x1 = x2;
							y1 = y2;
							area1 = area2;
							shuttle_size = shuttle_trajectory.size();
						}
						// following tracking (start at (_x, _y))
						track_index.push_back(shuttle_trajectory.size());
						
						// punctuation    pos:(-1,-1) area:-1
						shuttle_trajectory.push_back(std::make_pair(-1, -1));
						shuttle_area_value.push_back(-1);
					}
				}
			}

		}

	}
	
	prev_area_value.push_back(area_value);
	prev_gravity_pos.push_back(gravity_pos);

	// set previous frame limit (30 frames)
	if(prev_area_value.size() > 30) {
		prev_area_value.erase(prev_area_value.begin());
		prev_gravity_pos.erase(prev_gravity_pos.begin());
	}
	
	*dst = _dst.clone();
}





/*
 * mouse function----------------------------------------
 */
void CallBackFunc(int eventType, int x, int y, int flags, void* userdata)
{
	mouseParam *ptr = static_cast<mouseParam*> (userdata);

	ptr->x = x;
	ptr->y = y;
	ptr->event = eventType;
	ptr->flags = flags;
}





/*
 * diff 3 frames
 */
void moveObjDetection(cv::Mat im1, cv::Mat im2, cv::Mat im3, cv::Mat* dst)
{
	cv::Mat d1, d2, diff;
	cv::Mat im_mask, mask;

	cv::absdiff(im1, im2, d1);
	cv::absdiff(im2, im3, d2);
	cv::bitwise_and(d1, d2, diff);

	cv::threshold(diff, mask, 5, 1, cv::THRESH_BINARY);
	cv::threshold(mask, im_mask, 0, 255, cv::THRESH_BINARY);
	cv::medianBlur(im_mask, im_mask, 5);

	*dst = im_mask.clone();
}


/*
 * deinterlace
 */
void deinterlace(cv::Mat src, cv::Mat* dst)
{
	cv::Mat img = src.clone();

	char* line1;
	char* line2;
	char* line3;
	char* data = (char*)img.data;

	int index;
	int step = img.step;
	int height = img.rows;
	int width = img.cols;

	for(int y = 1; y < height-1; y+=2)
	{
		// set focus 3 lines
		line1 = data + ((y-1) * step);
		line2 = data + ((y) * step);
		line3 = data + ((y+1) * step);

		for(int x = 0; x < width; x++)
		{
			for(int c = 0; c < img.channels(); c++)
			{
				// center line value = upper/2 + lower/2;
				index = x * img.elemSize() + c;
				line2[index] = 0.5f*(unsigned char)line1[index] + 0.5f*(unsigned char)line2[index];
				//int a = img.data[y*(img.step) + x*(img.elemSize()) + c];
				//img.data[y*(img.step) + x*(img.elemSize()) + c] = 0;
			}
		}
	}

	if(height > 1 && height % 2 == 0)
	{
		line1 = data + ((height-2) * step);
		line2 = data + ((height-1) * step);
		memcpy(line2, line1, width);
	}

	// output
	dst->data = (unsigned char*)data;
}


/*
 * remove noise by erode and dilate 
 */
void erode_dilate(cv::Mat src, cv::Mat* dst, int num)
{
	cv::Mat img = src.clone();
	
	cv::erode(img, img, cv::Mat(), cv::Point(-1, -1), num);
	cv::dilate(img, img, cv::Mat(), cv::Point(-1, -1), num);

	*dst = img.clone();
}


/*
 * combine 2 images
 */
void combine_image(cv::Mat im1, cv::Mat im2, cv::Mat* dst)
{
	int width, height;

	width = im1.cols;
	height = im1.rows;

	// set output image
	cv::Mat output = cv::Mat::zeros(cv::Size(width*2, height), CV_8UC3);

	// set area
	cv::Rect roi1(0, 0, width, height);
	cv::Rect roi2(width, 0, width, height);

	// create roi image
	// copy to certain area
	cv::Mat output_roi = output(roi1);
	cv::Mat im1_roi = im1(roi1);

	im1_roi.copyTo(output_roi);
	
	cv::Mat output_roi2 = output(roi2);
	cv::Mat im2_roi = im2(roi1);
	
	im2_roi.copyTo(output_roi2);

	// set rect
	cv::rectangle(output, roi2, cv::Scalar(0, 0, 255), 2);
	
	*dst = output.clone();
}


/*
 * create mask image(8UC3)
 */
void create_image(cv::Mat src, cv::Mat* dst, cv::Mat mask)
{
	int width = src.cols;
	int height = src.rows;

	cv::Mat img = cv::Mat::zeros(cv::Size(width, height), CV_8UC3);

	for(int j = 0; j < height; j++) {
		for(int i = 0; i < width; i++) {
			for(int c = 0; c < 3; c++) {
				img.at<cv::Vec3b>(j, i)[c] = 255;
			}
		}
	}

	img.copyTo(*dst, mask);
}


/*
 * calculate angle
 */
void calc_angle(std::pair<int, int> p1, std::pair<int, int> p2, std::pair<int, int> p3, double* dst)
{
	std::pair<int, int> a = std::make_pair(p2.first - p1.first, p2.second - p1.second);
	std::pair<int, int> b = std::make_pair(p3.first - p2.first, p3.second - p2.second);

	int inner_product = a.first * b.first + a.second * b.second;
	double norm_a = sqrt(a.first * a.first + a.second * a.second);
	double norm_b = sqrt(b.first * b.first + b.second * b.second);

	*dst = acos(inner_product / (norm_a * norm_b)) / M_PI * 180;
}























