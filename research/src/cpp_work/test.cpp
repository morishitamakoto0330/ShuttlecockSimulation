#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/core.hpp"
#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/highgui.hpp"
#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/imgproc.hpp"
#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/opencv.hpp"
#include "/usr/local/Cellar/opencv/3.3.1_1/include/opencv2/imgproc/types_c.h"

#include <iostream>
#include <string>
#include <vector>
#include <math.h>

#include "./extraction.hpp"
#include "./debug.hpp"
#include "./matrix.hpp"

int main(void)
{
	Matrix m1(2, 2);
	Matrix m2(2, 1);
	m1.randomCreateMatrix();
	m2.randomCreateMatrix();

	m1.disp();
	m2.disp();
	
	m1.productMatrix(m2.getMatrix());
	m1.disp();
	/*
	std::vector<int> v{956, 582, 987, 556, 1017, 529, 1041, 505, 1066, 480, 1092, 456, 1114, 436, 1137, 418, 1160, 399, 1179, 385, 1200, 370, 1221, 356, 1239, 343, 1257, 331, 1274, 320, 1289, 311, 1305, 301, 1324, 292, 1336, 285, 1351, 277, 1366, 270, 1380, 264, 1395, 258, 1409, 254, 1426, 250, 1439, 247, 1452, 244, 1464, 241, 1476, 242, 1485, 240, 1502, 239, 1512, 239, 1523, 239, 1537, 242, 1546, 242, 1561, 243, 1569, 245, 1581, 249, 1591, 252, 1606, 258, 1617, 262, 1627, 268, 1637, 272, 1647, 278, 1658, 285, 1668, 290, 1678, 297, 1687, 304, 1694, 310, 1707, 321, 1717, 330, 1726, 339, 1735, 349, 1744, 359, 1753, 370, 1761, 380, 1770, 392, 1779, 403};

	int sum_x = 0;
	int sum_y = 0;
	int sum_xx = 0;
	int sum_xy = 0;
	int n = 2;

	for(int i = 0; i < n*2; i+=2) {
		sum_x += v[i];
		sum_y += v[i+1];
		sum_xy += v[i] * v[i+1];
		sum_xx += v[i] * v[i];
	}

	double a = (double)(sum_x*sum_y - n*sum_xy)/(sum_x*sum_x - n*sum_xx);
	std::cout << sum_x << "," << sum_y << "," << sum_xx << "," << sum_xy << std::endl;
	std::cout << "傾きa=" << a << std::endl;

	a *= -1;
	double v_size = 40.460;
	double vx = v_size/(sqrt(1 + a*a));
	double vy = vx*a;

	double theta = atan(a);

	std::cout << "(vx, vy)=(" << vx << "," << vy << ")" << std::endl;
	std::cout << "[pixel/((1/60)s)] -> [m/s]" << std::endl;
	std::cout << "(vx, vy)=(" << vx*60*4.31/1000 << "," << vy*60*4.31/1000 << ")" << std::endl;
	*/

	/*
	vx = v_size*cos(theta);
	vy = v_size*sin(theta);
	std::cout << "(vx, vy)=(" << vx << "," << vy << ")" << std::endl;
	*/




	/*
	* ---------------------------------------------------------------
	* extract data
	*/

	
	//std::string img_str = "../../res/image_simulate/shuttle_trajectory/shuttle_point.png";
	//std::string img_save_str = "../../res/image_simulate/shuttle_trajectory/shuttle_point_1_1.png";

	//cv::Mat img = cv::imread(img_str);
	
	// (x,y) data 1, 2, 3, 4
	
	// data1
	//std::vector<int> v{391, 537, 394, 510, 419, 476, 468, 451, 534, 463, 666, 442, 754, 437, 831, 441, 901, 444, 967, 444, 1027, 446, 1083, 449, 1136, 453, 1183, 457, 1232, 461, 1274, 464, 1316, 468, 1355, 473, 1394, 479, 1431, 484, 1464, 490, 1497, 495, 1529, 501, 1559, 509, 1589, 516, 1618, 523};
	//std::vector<int> v{534, 463, 666, 442, 754, 437, 831, 441, 901, 444, 967, 444, 1027, 446, 1083, 449, 1136, 453, 1183, 457, 1232, 461, 1274, 464, 1316, 468, 1355, 473, 1394, 479, 1431, 484, 1464, 490, 1497, 495, 1529, 501, 1559, 509, 1589, 516, 1618, 523};
	
	// data2
	//std::vector<int> v{405, 578, 432, 549, 489, 533, 568, 547, 710, 494, 796, 477, 875, 471, 949, 461, 1012, 451, 1071, 444, 1126, 438, 1177, 433, 1224, 429, 1269, 426, 1311, 423, 1351, 421, 1388, 419, 1423, 418, 1456, 417, 1488, 418, 1519, 418, 1547, 419, 1576, 421, 1601, 423};
	//std::vector<int> v{568, 547, 710, 494, 796, 477, 875, 471, 949, 461, 1012, 451, 1071, 444, 1126, 438, 1177, 433, 1224, 429, 1269, 426, 1311, 423, 1351, 421, 1388, 419, 1423, 418, 1456, 417, 1488, 418, 1519, 418, 1547, 419, 1576, 421, 1601, 423};
	
	// data3
	//std::vector<int> v{956, 582, 987, 556, 1017, 529, 1041, 505, 1066, 480, 1092, 456, 1114, 436, 1137, 418, 1160, 399, 1179, 385, 1200, 370, 1221, 356, 1239, 343, 1257, 331, 1274, 320, 1289, 311, 1305, 301, 1324, 292, 1336, 285, 1351, 277, 1366, 270, 1380, 264, 1395, 258, 1409, 254, 1426, 250, 1439, 247, 1452, 244, 1464, 241, 1476, 242, 1485, 240, 1502, 239, 1512, 239, 1523, 239, 1537, 242, 1546, 242, 1561, 243, 1569, 245, 1581, 249, 1591, 252, 1606, 258, 1617, 262, 1627, 268, 1637, 272, 1647, 278, 1658, 285, 1668, 290, 1678, 297, 1687, 304, 1694, 310, 1707, 321, 1717, 330, 1726, 339, 1735, 349, 1744, 359, 1753, 370, 1761, 380, 1770, 392, 1779, 403};
	
	// data4
	//std::vector<int> v{838, 601, 888, 595, 938, 594, 985, 591, 1029, 586, 1066, 582, 1109, 580, 1141, 579, 1177, 578, 1208, 578, 1240, 577, 1267, 577, 1298, 577, 1324, 577, 1349, 578, 1373, 580, 1396, 582, 1417, 583, 1440, 586, 1460, 590, 1485, 593, 1502, 596, 1519, 597};
	//std::vector<std::pair<int, int>> xy_data;

	// vector<int> -> pair<int, int>
	//for(int i = 0; i < v.size(); i+=2) {
		//xy_data.push_back(std::make_pair(v[i], v[i+1]));
		//cv::circle(img, cv::Point(v[i],v[i+1]), 4, cv::Scalar(0,0,255), 2);
	//}

	//cv::imwrite(img_save_str, img);
	
	
	/*
	* ------------------------------------------
	*/
	
	
	//cv::Mat img = cv::imread("../../res/image_simulate/shuttle_points/shuttle_point.png");
	//std::cout << "cols,rows=" << img.cols << "," << img.rows << std::endl;
	/*
	std::string img_str = "../../res/bayassan/image.JPG";
	std::string img_copy_str = "../../res/bayassan/image_copy.JPG";
	std::string tmp_str = "../../res/bayassan/template6.JPG";

	cv::Mat img = cv::imread(img_str);
	cv::Mat tmp = cv::imread(tmp_str);
	cv::Mat result;
	cv::Mat img_copy = img.clone();

	// template matching
	//cv::matchTemplate(img, tmp, result, cv::TM_SQDIFF);
	//cv::matchTemplate(img, tmp, result, cv::TM_SQDIFF_NORMED);
	//cv::matchTemplate(img, tmp, result, cv::TM_CCORR);
	cv::matchTemplate(img, tmp, result, cv::TM_CCORR_NORMED);
	//cv::matchTemplate(img, tmp, result, cv::TM_CCOEFF);
	//cv::matchTemplate(img, tmp, result, cv::TM_CCOEFF_NORMED);
	// get matching point
	double maxVal;
	cv::Point maxPt;
	cv::minMaxLoc(result, 0, &maxVal, 0, &maxPt);

	// save result image
	std::string s_result_str = "../../res/bayassan/result_save.JPG";
	cv::imwrite(s_result_str, result);

	// show matching area
	cv::rectangle(img_copy, maxPt, cv::Point(maxPt.x + tmp.cols, maxPt.y + tmp.rows), cv::Scalar(0, 255, 255), 2, 8, 0);
	cv::imwrite(img_copy_str, img_copy);
	*/
	
	
	
	
	
	/*
	int x, y, key;

	//std::string img_str = "../../res/image_progressive/shuttle_point/shuttle_point.png";

	
	cv::Mat frame, img;
	mouseParam mouseEvent;
	
	cv::VideoCapture cap("../../res/movie/progressive_1.MTS");

	if(!cap.isOpened()) return -1;
	std::cout << "MTS file open." << std::endl;

	cap >> frame;
	img = frame.clone();
	
	std::string win_name = "mouse_test";
	// create window
	cv::namedWindow(win_name);
	// set mouse event
	cv::setMouseCallback(win_name, CallBackFunc, &mouseEvent);
	
	cv::resize(img, img, cv::Size(), 0.5, 0.5);
	cv::imshow(win_name, img);
	
	
	while(1) {
		x = mouseEvent.x;
		y = mouseEvent.y;
		key = cv::waitKey(20);
		
		
		//cv::line(img, cv::Point(x,0), cv::Point(x,img.rows), cv::Scalar(0,0,255), 2);
		//cv::line(img, cv::Point(0,y), cv::Point(img.cols,y), cv::Scalar(0,0,255), 2);
		
		if(mouseEvent.event == cv::EVENT_LBUTTONDOWN) {
			std::cout << "(x,y)=(" << x << "," << y << ")" << std::endl;
		}


		// show image
		//cv::imwrite(img_str, img);
		//cv::resize(img, img, cv::Size(), 0.6, 0.6);
		//cv::imshow(win_name, img);
		
		if(key == 27) break;
	}
	*/
	return 0;
}




