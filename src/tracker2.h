#ifndef TRACKER2_H
#define TRACKER2_H

#include <opencv/cv.h>

#include "blob.h"
#include "disjointset.h"

class cTracker2 {
  private:
	double min_area, max_radius;
	
	node **labels; unsigned int width, height;

	// storage of the current blobs and the blobs from the previous frame
	vector<cBlob> blobs, blobs_previous;
	
	cDisjointSet ds;

  protected:

  public:

	cTracker2(double min_area, double max_radius);
	~cTracker2();

	void extractBlobs(cv::Mat &mat);
	void trackBlobs(cv::Mat &mat, bool history);
	void scaleBlobs();
	vector<cBlob>& getBlobs();
};


#endif