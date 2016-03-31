#ifndef OBJ_H
#define OBJ_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <GL/gl.h>
#include <GL/glu.h>
#include <math.h>
#include "misc.h"

struct vertex {
	std::vector<float> v;
	void normalize() {
		float magnitude = 0.0f;
		for (int i = 0; i < v.size(); i++)
			magnitude += pow(v[i], 2.0f);
		magnitude = sqrt(magnitude);
		for (int i = 0; i < v.size(); i++)
			v[i] /= magnitude;
	}
	vertex operator-(vertex v2) {
		vertex v3;
		if (v.size() != v2.v.size()) {
			v3.v.push_back(0.0f);
			v3.v.push_back(0.0f);
			v3.v.push_back(0.0f);
		} else {
			for (int i = 0; i < v.size(); i++)
				v3.v.push_back(v[i] - v2.v[i]);
		}
		return v3;
	}
	vertex cross(vertex v2) {
		vertex v3;
		if (v.size() != 3 || v2.v.size() != 3) {
			v3.v.push_back(0.0f);
			v3.v.push_back(0.0f);
			v3.v.push_back(0.0f);
		} else {
			v3.v.push_back(v[1]*v2.v[2]-v[2]*v2.v[1]);
			v3.v.push_back(v[2]*v2.v[0]-v[0]*v2.v[2]);
			v3.v.push_back(v[0]*v2.v[1]-v[1]*v2.v[0]);
		}
		return v3;
	}
};

struct face {
	std::vector<int> vertex;
	std::vector<int> texture;
	std::vector<int> normal;
};

class cObj {
  private:
  protected:
  public:
	double minx, miny, minz, maxx, maxy, maxz;
	std::vector<vertex> vertices;
	std::vector<vertex> texcoords;
	std::vector<vertex> normals;
	std::vector<vertex> parameters;
	std::vector<face> faces;

	GLuint list;
	cObj(std::string filename);
	~cObj();
	void compileList();
	void render();
};

#endif