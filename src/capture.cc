#include "capture.h"

cCapture::cCapture(int device, int width, int height) : capture_error(false), device(device), width(width), height(height) {
	cap = new cv::VideoCapture(device);
	capture_error = cap->isOpened() ? false : true;
	if (!capture_error) {
		cap->set(CV_CAP_PROP_FRAME_WIDTH, width);
		cap->set(CV_CAP_PROP_FRAME_HEIGHT, height);
	}
}

cCapture::~cCapture() {
	if (capture_error) return;
	delete cap;
}

void cCapture::capture() {
	if (capture_error) return;
	cap->operator>>(capture_frame);
}

cv::Mat& cCapture::captureFrame() {
	return capture_frame;
}

bool cCapture::captureError() {
	return capture_error;
}