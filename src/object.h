#ifndef OBJECT_H
#define OBJECT_H

#include <vector>
#include "vector.h"
#include "matrix.h"
#include "vec.h"
#include "common.h"
#include "enum.h"


// abstract base class: the plane, sphere, and triangle objects are derived from this
class cObject {
    public:
	int type;
	int material;
	vector3 color;
	vector3 emissive;
	int texture;
	vector3 e0, e0_, e1, e1_, e2, e2_;

	cObject();
	cObject(vector3 color, vector3 emissive, int material, int texture, vector3 e0, vector3 e1, vector3 e2);
	~cObject();
	virtual void applyCamera(cMatrix& camera) = 0;
	virtual void save() = 0;
};

class cPlane : public cObject {
    public:
	vector3 normal, normal_;
	vector3 point, point_;

	cPlane();
	cPlane(vector3 normal, vector3 point, vector3 color, vector3 emissive, int material, int texture, vector3 e0, vector3 e1);
	~cPlane();
	void applyCamera(cMatrix& camera);
	void save();
};

class cSphere : public cObject {
    public:
	vector3 center, center_;
	double radius, radius_;

	cSphere();
	cSphere(vector3 center, double radius, vector3 color, vector3 emissive, int material, int texture, vector3 e0, vector3 e1);
	~cSphere();
	void applyCamera(cMatrix& camera);
	void save();
};

class cTriangle : public cObject {
    public:
	vector3 p0, p1, p2, p0_, p1_, p2_;
	vector3 n0, n1, n2, n0_, n1_, n2_;
	vector3 n, n_;

	cTriangle();
	cTriangle(vector3 p0, vector3 p1, vector3 p2, vector3 n0, vector3 n1, vector3 n2, vector3 color, vector3 emissive, int material, int texture, vector3 e0, vector3 e1, vector3 e2);
	~cTriangle();
	void applyCamera(cMatrix& camera);
	void save();
};

void applyCamera(vector3 origin, vector3 forward, vector3 up, std::vector<cObject*> objects);
void copyObjects(std::vector<cObject*> objects, __object _objects[], __bounds& bounds);

#endif