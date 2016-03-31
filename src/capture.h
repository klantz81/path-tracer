#ifndef CAPTURE_H
#define CAPTURE_H

#include <opencv/cv.h>
#include <opencv/highgui.h>

class cCapture {
  private:
	bool capture_error;
	int device, width, height;
	cv::VideoCapture *cap;

  protected:
	cv::Mat capture_frame;

  public:
	cCapture(int device, int width, int height);
	~cCapture();

	void capture();

	cv::Mat& captureFrame();

	bool captureError();
};

#endif