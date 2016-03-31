#include "tracker2.h"

cTracker2::cTracker2(double min_area, double max_radius) :
	min_area(min_area), max_radius(max_radius),
	labels(NULL) {
}

cTracker2::~cTracker2() {
	if (labels) delete [] labels;
}

void cTracker2::extractBlobs(cv::Mat &mat) {
	// mat.cols, mat.rows -- allocate vectors
	if (mat.cols != width || mat.rows != height) {
		width = mat.cols; height = mat.rows;
		if (labels) delete [] labels;
		labels = new node*[width*height];
	}

	// reset our data structure for reuse
	ds.Reset();
	int index;

	// generate equivalence sets -- connected component labeling (4-connected)
	labels[0] = ds.MakeSet(0);
	for (int j = 1; j < mat.cols; j++)
		labels[j] = mat.data[j] != mat.data[j-1] ? ds.MakeSet(0) : labels[j-1];
	for (int j = mat.cols; j < mat.rows*mat.cols; j++) {
		if (mat.data[j] == mat.data[j-1]) {
			labels[j] = labels[j-1];
			if (mat.data[j-1] == mat.data[j-mat.cols]) ds.Union(labels[j-1], labels[j-mat.cols]);
		}
		else if (mat.data[j] == mat.data[j-mat.cols]) labels[j] = labels[j-mat.cols];
		else labels[j] = ds.MakeSet(0);
	}

	// the representative elements in our disjoint set data struct are associated with indices
	// we reduce those indices to 0,1,...,n and allocate our blobs
	cBlob temp;
	temp.event = BLOB_NULL;
	blobs.clear();
	for (int i = 0; i < ds.Reduce(); i++)
		blobs.push_back(temp);

	// populate our blob vector
	for (int j = 0; j < mat.rows; j++) {
		for (int i = 0; i < mat.cols; i++) {
			index = ds.Find(labels[j*mat.cols+i])->i;
			if (blobs[index].event == BLOB_NULL) {
				blobs[index].min.x = blobs[index].max.x = i;
				blobs[index].min.y = blobs[index].max.y = j;				 
				blobs[index].event  = BLOB_DOWN;
				blobs[index].height = 0;
			} else {
				if      (blobs[index].min.x > i) blobs[index].min.x = i;
				else if (blobs[index].max.x < i) blobs[index].max.x = i;
				blobs[index].max.y = j;
			}
		}
	}

	// apply blob filter
	for (int i = 0; i < blobs.size(); i++) {
		if ((blobs[i].max.x-blobs[i].min.x)*(blobs[i].max.y-blobs[i].min.y) < min_area) { blobs.erase(blobs.begin()+i); i--; }
	}

	// find blob centers
	for (int i = 0; i < blobs.size(); i++) {
		  blobs[i].location.x = blobs[i].origin.x = (blobs[i].max.x + blobs[i].min.x) / 2.0;
		  blobs[i].location.y = blobs[i].origin.y = (blobs[i].max.y + blobs[i].min.y) / 2.0;
	}

}

void cTracker2::trackBlobs(cv::Mat &mat, bool history) {
  	// clear the blobs from two frames ago
	blobs_previous.clear();
	
	// before we populate the blobs vector with the current frame, we need to store the live blobs in blobs_previous
	for (int i = 0; i < blobs.size(); i++)
		if (blobs[i].event != BLOB_UP)
			blobs_previous.push_back(blobs[i]);

	extractBlobs(mat);

	// initialize previous blobs to untracked
	for (int i = 0; i < blobs_previous.size(); i++) blobs_previous[i].tracked = false;

	// main tracking loop -- O(n^2) -- simply looks for a blob in the previous frame within a specified radius
	for (int i = 0; i < blobs.size(); i++) {
		for (int j = 0; j < blobs_previous.size(); j++) {
			if (blobs_previous[j].tracked) continue;

			if (sqrt(pow(blobs[i].location.x - blobs_previous[j].location.x, 2.0) + pow(blobs[i].location.y - blobs_previous[j].location.y, 2.0)) < max_radius) {
				blobs_previous[j].tracked = true;
				blobs[i].event = BLOB_MOVE;
				blobs[i].origin.x = history ? blobs_previous[j].origin.x : blobs_previous[j].location.x;
				blobs[i].origin.y = history ? blobs_previous[j].origin.y : blobs_previous[j].location.y;
			}
		}
	}

	// add any blobs from the previous frame that weren't tracked as having been removed
	for (int i = 0; i < blobs_previous.size(); i++) {
		if (!blobs_previous[i].tracked) {
			blobs_previous[i].event = BLOB_UP;
			blobs.push_back(blobs_previous[i]);
		}
	}
}

vector<cBlob>& cTracker2::getBlobs() {
	return blobs;
}
