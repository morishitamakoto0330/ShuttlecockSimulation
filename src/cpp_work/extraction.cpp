#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/core.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/highgui.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/imgproc.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/opencv.hpp"
#include "/usr/local/Cellar/opencv3/3.3.0_3/include/opencv2/imgproc/types_c.h"

#include <iostream>
#include <string>
#include <vector>

#include "./extraction.hpp"

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

	for(int i = 0; i < 256; i++)
	{
		for(int k = 0; k < 3; k++)
		{
			if(lower[k] <= upper[k])
			{
				if((lower[k] <= i) && (i <= upper[k]))
				{
					lut.data[i*lut.step+k] = 255;
				} else {
					lut.data[i*lut.step+k] = 0;
				}
			} else {
				if((i <= upper[k]) || (lower[k] <= i))
				{
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

void labeling(cv::Mat* input, cv::Mat* output)
{
	cv::Mat bin;
	cv::Mat stats;
	cv::Mat centroids;

	// binarization
	cv::threshold(*input, bin, 0, 255, cv::THRESH_BINARY | cv::THRESH_OTSU);
	// create image to label
	cv::Mat labelImage(input->size(), CV_32S);

	// labeling(simple)
	//int labelNum = cv::connectedComponents(bin, labelImage, 8);
	// labeling(detail)
	int labelNum = cv::connectedComponentsWithStats(bin, labelImage, stats, centroids);
	
	// label coloring
	std::vector<cv::Vec3b> colors(labelNum);
	colors[0] = cv::Vec3b(0, 0, 0);
	for(int label = 1; label < labelNum; label++)
	{
		colors[label] = cv::Vec3b((rand()&255), (rand()&255), (rand()&255));
	}

	// create result image
	cv::Mat dst(input->size(), CV_8UC3);
	for(int y = 0; y < dst.rows; y++)
	{
		for(int x = 0; x < dst.cols; x++)
		{
			int label = labelImage.at<int>(y, x);
			cv::Vec3b &pixel = dst.at<cv::Vec3b>(y, x);
			pixel = colors[label];
		}
	}

	// get parameter(center of gravity, area value)
	
	int x_gravity[labelNum];
	int y_gravity[labelNum];

	int area[labelNum];
	
	// center of gravity
	for(int i = 1; i < labelNum; i++)
	{
		double *param = centroids.ptr<double>(i);

		x_gravity[i-1] = static_cast<int>(param[0]);
		y_gravity[i-1] = static_cast<int>(param[0]);
	}

	// area
	for(int i = 1; i < labelNum; i++)
	{
		int *param = stats.ptr<int>(i);
		area[i-1] = param[cv::ConnectedComponentsTypes::CC_STAT_AREA];
	}
	
	/*
	// debug---------------
	int _x = 0;
	int count = 0;
	for(int i = 0; i < labelNum; i++)
	{
		_x = x_gravity[i];
		if((670 <= _x) && (_x <= 1350) && (10 <= area[i]))
		{
			count++;
			std::cout << i << ":" << area[i] << ", ";
		}
	}
	std::cout << std::endl;
	std::cout << "num_good=" << count << ", ";
	std::cout << "num=" << labelNum << std::endl;
	*/

	*output = dst.clone();
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
 * union lablel and delete noise
 */
void unionLabel(cv::Mat* src, cv::Mat* dst)
{
	// access pixel
	for(int y = 0; y < src->rows; y++)
	{
		for(int x = 0; x < src->cols; x++)
		{
			for(int c = 0; c < src->channels(); c++)
			{
				int a = src->data[y*(src->step) + x*(src->elemSize()) + c];
				std::cout << a << " ";
				dst->data[y*(src->step) + x*(src->elemSize()) + c] = 0;
			}
		}
		std::cout << std::endl;
		break;
	}
}



















