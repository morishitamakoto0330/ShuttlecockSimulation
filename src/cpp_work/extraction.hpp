// global variable
extern int point_number;
extern int area_threshold;
extern int x_threshold;

extern std::string img_white_name;
extern std::vector<std::pair<int, int>> shuttle_trajectory;
extern std::vector<int> shuttle_area_value;
extern std::vector<std::vector<std::pair<int, int>>> prev_gravity_pos;
extern std::vector<std::vector<int>> prev_area_value;
extern std::vector<int> track_index;



struct mouseParam {
	int x;
	int y;
	int event;
	int flags;
};

// mouse function
void CallBackFunc(int eventType, int x, int y, int flags, void* userdata);

void colorExtraction(cv::Mat* src, cv::Mat* dst,
		int code,
		int ch1Lower, int ch1Upper,
		int ch2Lower, int ch2Upper,
		int ch3Lower, int ch3Upper
);

void labeling(cv::Mat src, cv::Mat* dst);

void deinterlace(cv::Mat src, cv::Mat* dst);

void moveObjDetection(cv::Mat im1, cv::Mat im2, cv::Mat im3, cv::Mat* dst);

void erode_dilate(cv::Mat src, cv::Mat* dst, int num);

void combine_image(cv::Mat im1, cv::Mat im2, cv::Mat* dst);

void create_image(cv::Mat src, cv::Mat* dst, cv::Mat mask);

void calc_angle(std::pair<int, int> p1, std::pair<int, int> p2, std::pair<int, int> p3, double* dst);



