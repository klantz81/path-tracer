#include "fiducial.h"

cFiducialIdentifier::cFiducialIdentifier(double orientation, bool value) : orientation(orientation), value(value) {
}

cFiducialIdentifier::~cFiducialIdentifier() {
}

bool cFiducialIdentifier::operator<(const cFiducialIdentifier& f) const {
	return (orientation < f.orientation);
}


cFiducial::cFiducial() {
}

cFiducial::~cFiducial() {
}

#include <iostream>

std::vector<fiducial>& cFiducial::extractFiducials(std::vector<cBlob>& blobs) {
	fiducials.clear();

	if (blobs.size() < 5) return fiducials;

	sort(blobs.begin(), blobs.end());

	for (int i = 0; i < blobs.size() - 1; i++) {
		if (blobs[i].event == BLOB_UP) continue;
		for (int j = i + 1; j < blobs.size(); j++) {
			if (blobs[j].event == BLOB_UP) continue;
			if (blobs[j].contains(blobs[i], i)) break;
		}
	}

	fiducial temp; // center, dimensions, orientation, id
	double theta;
	double x, y; int child1, child2, child3, child4;

	for (int i = 0; i < blobs.size(); i++) {
		identifier.clear();
		if (blobs[i].height == 4) {										// possible fiducial
			child1 = blobs[i].children[0];									// first child (ones)
			if (blobs[child1].height != 3) break;
 
			child2 = blobs[child1].children[blobs[child1].children.size()-1];				// second child (zeroes)
			if (blobs[child2].height != 2) break;

			child3 = blobs[child2].children[blobs[child2].children.size()-1];				// third child
			if (blobs[child3].height != 1) break;

			temp.center.x = blobs[child3].location.x;							// evaluate center
			temp.center.y = blobs[child3].location.y;

			child4 = blobs[child3].children[0];
			x = blobs[child4].location.x - temp.center.x;
			y = blobs[child4].location.y - temp.center.y;
			temp.orientation = atan2(y, x);									// evaluate orientation
			if (temp.orientation < 0) temp.orientation += 2 * M_PI;
			
			//process id
			for (int j = 0; j < blobs[child1].children.size() - 1; j++) {
				x = blobs[blobs[child1].children[j]].location.x - temp.center.x;
				y = blobs[blobs[child1].children[j]].location.y - temp.center.y;
				theta = atan2(y, x); if (theta < 0) theta += 2 * M_PI; theta -= temp.orientation; if (theta < 0) theta += 2 * M_PI;
				identifier.push_back(cFiducialIdentifier(theta, true));  // ones
			}
			for (int j = 0; j < blobs[child2].children.size() - 1; j++) {
				x = blobs[blobs[child2].children[j]].location.x - temp.center.x;
				y = blobs[blobs[child2].children[j]].location.y - temp.center.y;
				theta = atan2(y, x); if (theta < 0) theta += 2 * M_PI; theta -= temp.orientation; if (theta < 0) theta += 2 * M_PI;
				identifier.push_back(cFiducialIdentifier(theta, false));  // zeroes
			}

			sort(identifier.begin(), identifier.end());
			temp.id = generateID();
			fiducials.push_back(temp);
		}
	}
	return fiducials;
}

int cFiducial::generateID() {
	int id = 0;
	for (int i = 0; i < identifier.size(); i++) id += (identifier[i].value ? 1 : 0)*pow(2, i);
	return id;
}