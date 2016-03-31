#ifndef FILTER_H
#define FILTER_H

#include "capture.h"

class cFilter : public cCapture {
  private:
	bool filter_error, balance_flag;
	int kernel_size, block_size, c;
	double std_dev;

  protected:
	cv::Mat filter_frame, gaussian_frame, balance_frame, threshold_frame;

  public:
	cFilter(int device, int width, int height, int kernel_size, double std_dev, int block_size, int c);
	~cFilter();

	void filter();
	void balance(bool flag);

	cv::Mat& filterFrame();
	cv::Mat& gaussianFrame();
	cv::Mat& thresholdFrame();

	bool filterError();
};

#endif