#include "object.h"



cObject::cObject() : color(vector3(1.0,1.0,1.0)), emissive(vector3(0.0,0.0,0.0)), material(DIFFUSE), texture(-1), e0(0.0,1.0,0.0), e1(0.0,0.0,1.0), e2(0.0,0.0,0.0), e0_(0.0,1.0,0.0), e1_(0.0,0.0,1.0), e2_(0.0,0.0,0.0) { }
cObject::cObject(vector3 color, vector3 emissive, int material, int texture, vector3 e0, vector3 e1, vector3 e2) : color(color), emissive(emissive), material(material), texture(texture), e0(e0), e1(e1), e2(e2), e0_(e0), e1_(e1), e2_(e2) { }
cObject::~cObject() { }



cPlane::cPlane() : normal(vector3(0.0,0.0,1.0)), point(vector3()), normal_(vector3(0.0,0.0,1.0)), point_(vector3()), cObject(vector3(1.0,1.0,1.0), vector3(0.0,0.0,0.0), DIFFUSE, -1, vector3(0.0,1.0,0.0), vector3(0.0,0.0,1.0), vector3()) { this->type = PLANE; }
cPlane::cPlane(vector3 normal, vector3 point, vector3 color, vector3 emissive, int material, int texture, vector3 e0, vector3 e1) : normal(normal), point(point), normal_(normal), point_(point), cObject(color, emissive, material, texture, e0, e1, vector3()) { this->type = PLANE; }
cPlane::~cPlane() { }
void cPlane::applyCamera(cMatrix& camera) {
	vector3 p = this->point_ + this->normal_;
	this->point = camera.mult(this->point_);
	p = camera.mult(p);
	this->normal = (p - this->point).unit();

	p = this->point_ + e0_;
	p = camera.mult(p);
	this->e0 = (p - this->point);
	p = this->point_ + e1_;
	p = camera.mult(p);
	this->e1 = (p - this->point);
}
void cPlane::save() {
	this->normal_ = this->normal;
	this->point_ = this->point;
	this->e0_ = this->e0;
	this->e1_ = this->e1;
}



cSphere::cSphere() : center(vector3(0.0,0.0,0.0)), center_(vector3(0.0,0.0,0.0)), radius(1.0), radius_(1.0), cObject(vector3(1.0,1.0,1.0), vector3(0.0,0.0,0.0), DIFFUSE, -1, vector3(0.0,1.0,0.0), vector3(0.0,0.0,1.0), vector3()) { this->type = SPHERE; }
cSphere::cSphere(vector3 center, double radius, vector3 color, vector3 emissive, int material, int texture, vector3 e0, vector3 e1) : center(center), center_(center), radius(radius), radius_(radius), cObject(color, emissive, material, texture, e0, e1, vector3()) { this->type = SPHERE; }
cSphere::~cSphere() { }
void cSphere::applyCamera(cMatrix& camera) {
	this->center = camera.mult(this->center_);

	vector3 p = this->center_ + e0_;
	p = camera.mult(p);
	this->e0 = (p - this->center);
	p = this->center_ + e1_;
	p = camera.mult(p);
	this->e1 = (p - this->center);
}
void cSphere::save() {
	this->center_ = this->center;
	this->radius_ = this->radius;
	this->e0_ = this->e0;
	this->e1_ = this->e1;
}



cTriangle::cTriangle() :
    p0(vector3(0.0,1.0,0.0)),  p1(vector3(-1.0,0.0,0.0)),  p2(vector3(1.0,0.0,0.0)),
    p0_(vector3(0.0,1.0,0.0)), p1_(vector3(-1.0,0.0,0.0)), p2_(vector3(1.0,0.0,0.0)),
    n0(vector3(0.0,0.0,1.0)),  n1(vector3(0.0,0.0,1.0)),  n2(vector3(0.0,0.0,1.0)),
    n0_(vector3(0.0,0.0,1.0)), n1_(vector3(0.0,0.0,1.0)), n2_(vector3(0.0,0.0,1.0)),
    cObject(vector3(1.0,1.0,1.0),vector3(0.0,0.0,0.0), DIFFUSE, -1, vector3(0.0,0.0,0.0), vector3(1.0,0.0,0.0), vector3(0.0,1.0,0.0)) {
	this->type = TRIANGLE;
	this->n = this->n_ = ((p1-p0).cross(p2-p0)).unit();
}
cTriangle::cTriangle(vector3 p0, vector3 p1, vector3 p2, vector3 n0, vector3 n1, vector3 n2, vector3 color, vector3 emissive, int material, int texture, vector3 e0, vector3 e1, vector3 e2) :
    p0(p0),  p1(p1),  p2(p2),
    p0_(p0), p1_(p1), p2_(p2),
    n0(n0),  n1(n1),  n2(n2),
    n0_(n0), n1_(n1), n2_(n2),
    cObject(color, emissive, material, texture, e0, e1, e2) {
	this->type = TRIANGLE;
	this->n = this->n_ = ((p1-p0).cross(p2-p0)).unit();
}
cTriangle::~cTriangle() { }
void cTriangle::applyCamera(cMatrix& camera) {
	vector3 _n0, _n1, _n2;
	vector3 _n;

	_n  = this->p0_ + this->n_;
	_n0 = this->p0_ + this->n0_;
	_n1 = this->p1_ + this->n1_;
	_n2 = this->p2_ + this->n2_;

	this->p0 = camera.mult(this->p0_);
	this->p1 = camera.mult(this->p1_);
	this->p2 = camera.mult(this->p2_);

	_n  = camera.mult(_n);
	_n0 = camera.mult(_n0);
	_n1 = camera.mult(_n1);
	_n2 = camera.mult(_n2);

	this->n  = (_n  - this->p0).unit();
	this->n0 = (_n0 - this->p0).unit();
	this->n1 = (_n1 - this->p1).unit();
	this->n2 = (_n2 - this->p2).unit();
}
void cTriangle::save() {
	this->p0_ = this->p0;
	this->p1_ = this->p1;
	this->p2_ = this->p2;
	this->n0_ = this->n0;
	this->n1_ = this->n1;
	this->n2_ = this->n2;
	this->n_  = this->n;
}

void applyCamera(vector3 origin, vector3 forward, vector3 up, std::vector<cObject*> objects) {
	vector3 right = (forward.cross(up)).unit();
		   up = (right.cross(forward)).unit();
//	forward = forward*-1;
//	vector3 right = (up.cross(forward)).unit();
//		   up = (forward.cross(right)).unit();
	double cam[] = {   right.x,    right.y,    right.z, -(right*origin),
			      up.x,       up.y,       up.z, -(up*origin),
			-forward.x, -forward.y, -forward.z,  (forward*origin),
				 0,         0,         0,    1 };

	cMatrix camera(4, 4, cam);

	for (int i = 0; i < objects.size(); i++) {
		objects[i]->applyCamera(camera);
		//objects[i]->save();
	}
}

void copyObjects(std::vector<cObject*> objects, __object _objects[], __bounds& bounds) {
	for (int i = 0; i < objects.size(); i++) {

		_objects[i].type     = objects[i]->type;
		_objects[i].material = objects[i]->material;
		_objects[i].texture  = objects[i]->texture;

		_objects[i].color.x = objects[i]->color.x;
		_objects[i].color.y = objects[i]->color.y;
		_objects[i].color.z = objects[i]->color.z;

		_objects[i].emission.x = objects[i]->emissive.x;
		_objects[i].emission.y = objects[i]->emissive.y;
		_objects[i].emission.z = objects[i]->emissive.z;

		_objects[i].e0.x = objects[i]->e0.x;
		_objects[i].e0.y = objects[i]->e0.y;
		_objects[i].e0.z = objects[i]->e0.z;
		_objects[i].e1.x = objects[i]->e1.x;
		_objects[i].e1.y = objects[i]->e1.y;
		_objects[i].e1.z = objects[i]->e1.z;
		_objects[i].e2.x = objects[i]->e2.x;
		_objects[i].e2.y = objects[i]->e2.y;
		_objects[i].e2.z = objects[i]->e2.z;

		if (objects[i]->type == PLANE) {
			_objects[i].normal.x = ((cPlane *)objects[i])->normal.x;
			_objects[i].normal.y = ((cPlane *)objects[i])->normal.y;
			_objects[i].normal.z = ((cPlane *)objects[i])->normal.z;
			_objects[i].point.x = ((cPlane *)objects[i])->point.x;
			_objects[i].point.y = ((cPlane *)objects[i])->point.y;
			_objects[i].point.z = ((cPlane *)objects[i])->point.z;
		} else if (objects[i]->type == SPHERE) {
			_objects[i].center.x = ((cSphere *)objects[i])->center.x;
			_objects[i].center.y = ((cSphere *)objects[i])->center.y;
			_objects[i].center.z = ((cSphere *)objects[i])->center.z;
			_objects[i].radius = ((cSphere *)objects[i])->radius;

			bounds.minx = MIN(bounds.minx,_objects[i].center.x-_objects[i].radius);
			bounds.miny = MIN(bounds.miny,_objects[i].center.y-_objects[i].radius);
			bounds.minz = MIN(bounds.minz,_objects[i].center.z-_objects[i].radius);
			bounds.maxx = MAX(bounds.maxx,_objects[i].center.x+_objects[i].radius);
			bounds.maxy = MAX(bounds.maxy,_objects[i].center.y+_objects[i].radius);
			bounds.maxz = MAX(bounds.maxz,_objects[i].center.z+_objects[i].radius);
		} else if (objects[i]->type == TRIANGLE) {
			_objects[i].p0.x = ((cTriangle *)objects[i])->p0.x;
			_objects[i].p0.y = ((cTriangle *)objects[i])->p0.y;
			_objects[i].p0.z = ((cTriangle *)objects[i])->p0.z;
			_objects[i].p1.x = ((cTriangle *)objects[i])->p1.x;
			_objects[i].p1.y = ((cTriangle *)objects[i])->p1.y;
			_objects[i].p1.z = ((cTriangle *)objects[i])->p1.z;
			_objects[i].p2.x = ((cTriangle *)objects[i])->p2.x;
			_objects[i].p2.y = ((cTriangle *)objects[i])->p2.y;
			_objects[i].p2.z = ((cTriangle *)objects[i])->p2.z;
			_objects[i].n0.x = ((cTriangle *)objects[i])->n0.x;
			_objects[i].n0.y = ((cTriangle *)objects[i])->n0.y;
			_objects[i].n0.z = ((cTriangle *)objects[i])->n0.z;
			_objects[i].n1.x = ((cTriangle *)objects[i])->n1.x;
			_objects[i].n1.y = ((cTriangle *)objects[i])->n1.y;
			_objects[i].n1.z = ((cTriangle *)objects[i])->n1.z;
			_objects[i].n2.x = ((cTriangle *)objects[i])->n2.x;
			_objects[i].n2.y = ((cTriangle *)objects[i])->n2.y;
			_objects[i].n2.z = ((cTriangle *)objects[i])->n2.z;
			_objects[i].n.x = ((cTriangle *)objects[i])->n.x;
			_objects[i].n.y = ((cTriangle *)objects[i])->n.y;
			_objects[i].n.z = ((cTriangle *)objects[i])->n.z;

			bounds.minx = MIN(bounds.minx,_objects[i].p0.x);
			bounds.miny = MIN(bounds.miny,_objects[i].p0.y);
			bounds.minz = MIN(bounds.minz,_objects[i].p0.z);
			bounds.maxx = MAX(bounds.maxx,_objects[i].p0.x);
			bounds.maxy = MAX(bounds.maxy,_objects[i].p0.y);
			bounds.maxz = MAX(bounds.maxz,_objects[i].p0.z);
			bounds.minx = MIN(bounds.minx,_objects[i].p1.x);
			bounds.miny = MIN(bounds.miny,_objects[i].p1.y);
			bounds.minz = MIN(bounds.minz,_objects[i].p1.z);
			bounds.maxx = MAX(bounds.maxx,_objects[i].p1.x);
			bounds.maxy = MAX(bounds.maxy,_objects[i].p1.y);
			bounds.maxz = MAX(bounds.maxz,_objects[i].p1.z);
			bounds.minx = MIN(bounds.minx,_objects[i].p2.x);
			bounds.miny = MIN(bounds.miny,_objects[i].p2.y);
			bounds.minz = MIN(bounds.minz,_objects[i].p2.z);
			bounds.maxx = MAX(bounds.maxx,_objects[i].p2.x);
			bounds.maxy = MAX(bounds.maxy,_objects[i].p2.y);
			bounds.maxz = MAX(bounds.maxz,_objects[i].p2.z);
		}
	}
}