#ifndef FIDUCIAL_H
#define FIDUCIAL_H

#include <algorithm>
#include <vector>
#include <math.h>
#include "blob.h"

struct fiducial {
	point center;
	point dimensions;
	double orientation;
	int id;
};

class cFiducialIdentifier {
  private:
	double orientation;
	bool   value;

  protected:
    
  public:
	cFiducialIdentifier(double orientation, bool value);
	~cFiducialIdentifier();
	bool operator<(const cFiducialIdentifier& f) const;
	
	friend class cFiducial;
};

class cFiducial {
  private:
	std::vector<fiducial> fiducials;
	std::vector<cFiducialIdentifier> identifier;
    
  protected:
    
  public:
	cFiducial();
	~cFiducial();
	
	std::vector<fiducial>& extractFiducials(std::vector<cBlob>& blobs);
	int generateID();
};

#endif