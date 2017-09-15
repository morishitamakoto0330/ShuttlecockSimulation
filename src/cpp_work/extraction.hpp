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

void labeling(cv::Mat* src, cv::Mat* dst);

void moveObjDetection(cv::Mat im1, cv::Mat im2, cv::Mat im3, cv::Mat* dst);


