#include "/usr/local/Cellar/opencv3/3.3.1_1/include/opencv2/opencv.hpp"
#include <iostream>
#include <stdio.h>

void file_read(std::vector<int> *v, int data_num);

int main(void)
{
	std::string win_str = "debug_window";
	cv::namedWindow(win_str);

	int h_cross_count = 10;
	int v_cross_count = 7;
	int num_of_images = 7;

	cv::Mat capture_image = cv::imread("../../res/image_reality/shuttle_point_4.png");

	// calibrate points
	std::vector<std::vector<cv::Point3f>> object_points;
	std::vector<std::vector<cv::Point2f>> image_points;
	std::vector<cv::Point3f> obj;

	for(int i = 0; i < h_cross_count*v_cross_count; i++) {
		obj.push_back(cv::Point3f(i/h_cross_count, i%h_cross_count, 0.0f));
	}

	for(int i = 0; i < num_of_images; i++) {
		char filename[32];
		sprintf(filename, "../../res/calibrate%d.JPG", i);

		cv::Mat frame = cv::imread(filename);
		cv::Mat gray;

		cv::flip(frame, frame, -1);
		cv::cvtColor(frame, gray, CV_BGR2GRAY);

		// find chessbord(10*7)
		cv::Size chessbordPatterns(h_cross_count, v_cross_count);
		std::vector<cv::Point2f> centers;
		bool found = cv::findChessboardCorners(gray, chessbordPatterns, centers, cv::CALIB_CB_ADAPTIVE_THRESH | cv::CALIB_CB_NORMALIZE_IMAGE);

		if(found) {
			cv::cornerSubPix(gray, centers, cv::Size(11, 11), cv::Size(-1, -1), cv::TermCriteria(cv::TermCriteria::EPS + cv::TermCriteria::COUNT, 30, 0.1));
			object_points.push_back(obj);
			image_points.push_back(centers);

			// draw
			cv::drawChessboardCorners(gray, chessbordPatterns, cv::Mat(centers), true);
			cv::imshow(win_str, gray);
		} else {
			std::cout << "not found" << std::endl;
		}

		// exit if put 'q'
		int key = cv::waitKey(1);
		if(key == 'q') break;
	}

	// calibrate
	std::vector<cv::Mat> rvecs, tvecs;
	cv::Mat mtx(3, 3, CV_64F);
	cv::Mat dist(8, 1, CV_64F);
	cv::calibrateCamera(object_points, image_points, cv::Size(capture_image.cols, capture_image.rows), mtx, dist, rvecs, tvecs);


	cv::Mat undistorted;
	cv::undistort(capture_image, undistorted, mtx, dist);
	
	/*
	while(1) {
		cv::undistort(capture_image, undistorted, mtx, dist);

		cv::imshow(win_str, undistorted);
		
		
		// exit if put 'q'
		int key = cv::waitKey(1);
		if(key == 'q') break;
	}
	*/


	// real points
	std::vector<int> v;
	file_read(&v, 4);

	for(int i = 0; i < v.size(); i=i+2) {
		cv::circle(undistorted, cv::Point(v[i], v[i+1]), 3, cv::Scalar(0, 0, 0), 2);
	}

	cv::imwrite("../../res/hoge4.png", undistorted);



	return 0;
}








void file_read(std::vector<int> *v, int data_num)
{
	std::ifstream ifs("../python_work/_data_1215.txt");
	std::string str, buf;
	std::vector<int> _v;
	int data_count = 1;
	
	if(ifs.fail()) {
		std::cerr << "failed to read file." << std::endl;
		exit(1);
	}
	
	while(getline(ifs, str)) {
		std::stringstream _str;
		if(str == "")  data_count++;
		_str << str;
		
		if(data_count == data_num) {
			while(getline(_str, buf, ',')) {
				_v.push_back(std::stoi(buf));
			}
		} else if(data_count > data_num) break;
	}
	
	*v = _v;
}















