#include "filter.h"

cFilter::cFilter(int device, int width, int height, int kernel_size, double std_dev, int block_size, int c) : filter_error(false),
													     balance_flag(false),
													     kernel_size(kernel_size),
													     std_dev(std_dev),
													     block_size(block_size),
													     c(c),
													     cCapture(device, width, height) {
	filter_error = captureError();
}

cFilter::~cFilter() {
}

// capture frame, convert to grayscale, apply gaussian blur, apply balance (if applicable), and apply adaptive threshold method
void cFilter::filter() {
	if (filter_error) return;
	capture();
	cvtColor(capture_frame, filter_frame, CV_BGR2GRAY);
	GaussianBlur(filter_frame, gaussian_frame, cv::Size(kernel_size, kernel_size), std_dev, std_dev);
	if (balance_flag) absdiff(gaussian_frame, balance_frame, gaussian_frame);
	adaptiveThreshold(gaussian_frame, threshold_frame, 255, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY, block_size, c);
}

// capture frame, convert to grayscale, apply gaussian blur --> this frame will be used in the filter method to balance the image before thresholding
void cFilter::balance(bool flag) {
	if (filter_error) return;
	balance_flag = flag;
	if (balance_flag) {
		capture();
		cvtColor(capture_frame, filter_frame, CV_BGR2GRAY);
		GaussianBlur(filter_frame, balance_frame, cv::Size(kernel_size, kernel_size), std_dev, std_dev);
	}
}

cv::Mat& cFilter::filterFrame() {
	return filter_frame;
}

cv::Mat& cFilter::gaussianFrame() {
	return gaussian_frame;
}

cv::Mat& cFilter::thresholdFrame() {
	return threshold_frame;
}

bool cFilter::filterError() {
	return filter_error;
}