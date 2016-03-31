#include <sstream>
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <SDL/SDL.h>
#include <SDL/SDL_ttf.h>
#include <SDL/SDL_image.h>
#include <math.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include "glm-0.9.2.6/glm/glm.hpp"
#include "glm-0.9.2.6/glm/gtc/matrix_transform.hpp"
#include "glm-0.9.2.6/glm/gtc/type_ptr.hpp"

#include "src/vec.h"
#include "src/common.h"
#include "src/enum.h"
#include "src/vector.h"
#include "src/keyboard.h"
#include "src/timer.h"
#include "src/matrix.h"
#include "src/glhelper.h"
#include "src/obj.h"
#include "src/object.h"
#include "src/buffer.h"


struct vertex_ {
	double x,  y,   z;
	double tx, ty, tz;
	double t[6];
};

bool hyperplaneSeparation(__vector n,
			  __vector p0, __vector p1, __vector p2,
			  double halfWidth, double halfHeight, double halfDepth) {
	double _p0 = n * p0,
	       _p1 = n * p1,
	       _p2 = n * p2;
	double min = MIN(_p0, MIN(_p1, _p2)),
	       max = MAX(_p0, MAX(_p1, _p2));
	double r = halfWidth  * fabs(n.x) +
	           halfHeight * fabs(n.y) +
		   halfDepth  * fabs(n.z);
	return -r > max || r < min;
}

_box* deleteTree(_box* node) {
	if (node->child0) deleteTree(node->child0);
	if (node->child1) deleteTree(node->child1);
	delete node;
}

struct __splitting_plane {
	double value;
	char event;
};

bool __splitting_plane_sort(__splitting_plane a, __splitting_plane b) { return a.value < b.value; }

_box* buildTree(unsigned char kdType, __node& node, _box* parent, unsigned char max_depth, unsigned char max_objects, __object objects[], unsigned int object_count, const __bounds& bounds, unsigned char depth) {

	static unsigned int id = 0, leaf_id = 0;

	_box *b = new _box;
	b->id = id++;
	b->leaf_id = -1;
	b->parent = parent;
	b->child0 = b->child1 = 0;
	b->bounds = bounds;

	std::vector<double> list; // for median selection
	std::vector<__splitting_plane> listx, listy, listz; __splitting_plane sp; // for SAH split selection

	unsigned char axis = NOAXIS;

	if (kdType == KD_EVEN || kdType == KD_MEDIAN) { // split along the axis of the largest dimension
		double _x = bounds.maxx - bounds.minx;
		double _y = bounds.maxy - bounds.miny;
		double _z = bounds.maxz - bounds.minz;
		axis = _x > _y && _x > _z ? X : (_y > _z ? Y : Z);
	}
	
	if (depth == 0) { // for the root node, add all the primitives
		for (int i = 0; i < object_count; i++) {
			if (objects[i].type == TRIANGLE) {
				b->objects.push_back(i);

				if (kdType == KD_EVEN) {
				  
				} else if (kdType == KD_MEDIAN) {
					if (axis == Z) {
						list.push_back(objects[i].p0.z); list.push_back(objects[i].p1.z); list.push_back(objects[i].p2.z);
					} else if (axis == Y) {
						list.push_back(objects[i].p0.y); list.push_back(objects[i].p1.y); list.push_back(objects[i].p2.y);
					} else if (axis == X) {
						list.push_back(objects[i].p0.x); list.push_back(objects[i].p1.x); list.push_back(objects[i].p2.x);
					}
				} else if (kdType == KD_SAH) {
					sp.value = MIN(objects[i].p0.z,MIN(objects[i].p1.z,objects[i].p2.z)); sp.event = PRIMITIVE_START; listz.push_back(sp);
					sp.value = MAX(objects[i].p0.z,MAX(objects[i].p1.z,objects[i].p2.z)); sp.event = PRIMITIVE_END;   listz.push_back(sp);
					sp.value = MIN(objects[i].p0.y,MIN(objects[i].p1.y,objects[i].p2.y)); sp.event = PRIMITIVE_START; listy.push_back(sp);
					sp.value = MAX(objects[i].p0.y,MAX(objects[i].p1.y,objects[i].p2.y)); sp.event = PRIMITIVE_END;   listy.push_back(sp);
					sp.value = MIN(objects[i].p0.x,MIN(objects[i].p1.x,objects[i].p2.x)); sp.event = PRIMITIVE_START; listx.push_back(sp);
					sp.value = MAX(objects[i].p0.x,MAX(objects[i].p1.x,objects[i].p2.x)); sp.event = PRIMITIVE_END;   listx.push_back(sp);
				}
			} else if (objects[i].type == SPHERE) {
				b->objects.push_back(i);

				if (kdType == KD_EVEN) {
				  
				} else if (kdType == KD_MEDIAN) {
					if (axis == Z) {
						list.push_back(objects[i].center.z-objects[i].radius); list.push_back(objects[i].center.z+objects[i].radius);
					} else if (axis == Y) {
						list.push_back(objects[i].center.y-objects[i].radius); list.push_back(objects[i].center.y+objects[i].radius);
					} else if (axis == X) {
						list.push_back(objects[i].center.x-objects[i].radius); list.push_back(objects[i].center.x+objects[i].radius);
					}
				} else if (kdType == KD_SAH) {
					sp.value = objects[i].center.z-objects[i].radius; sp.event = PRIMITIVE_START; listz.push_back(sp);
					sp.value = objects[i].center.z+objects[i].radius; sp.event = PRIMITIVE_END;   listz.push_back(sp);
					sp.value = objects[i].center.y-objects[i].radius; sp.event = PRIMITIVE_START; listy.push_back(sp);
					sp.value = objects[i].center.y+objects[i].radius; sp.event = PRIMITIVE_END;   listy.push_back(sp);
					sp.value = objects[i].center.x-objects[i].radius; sp.event = PRIMITIVE_START; listx.push_back(sp);
					sp.value = objects[i].center.x+objects[i].radius; sp.event = PRIMITIVE_END;   listx.push_back(sp);
				}
			}
		}
	} else { // apply the hyperplane separation theorem to inner nodes
		int index;
		for (int i = 0; i < b->parent->objects.size(); i++) {
			index = b->parent->objects[i];
			if (objects[index].type == TRIANGLE) {
				// minimum bounding box
				double sx = MIN(objects[index].p0.x, MIN(objects[index].p1.x, objects[index].p2.x));
				double sy = MIN(objects[index].p0.y, MIN(objects[index].p1.y, objects[index].p2.y));
				double sz = MIN(objects[index].p0.z, MIN(objects[index].p1.z, objects[index].p2.z));
				double bx = MAX(objects[index].p0.x, MAX(objects[index].p1.x, objects[index].p2.x));
				double by = MAX(objects[index].p0.y, MAX(objects[index].p1.y, objects[index].p2.y));
				double bz = MAX(objects[index].p0.z, MAX(objects[index].p1.z, objects[index].p2.z));
				if (sx > b->bounds.maxx || bx < b->bounds.minx) continue;
				if (sy > b->bounds.maxy || by < b->bounds.miny) continue;
				if (sz > b->bounds.maxz || bz < b->bounds.minz) continue;

				double _width  = b->bounds.maxx - b->bounds.minx;
				double _height = b->bounds.maxy - b->bounds.miny;
				double _depth  = b->bounds.maxz - b->bounds.minz;
				__vector o(b->bounds.minx + _width / 2.0, b->bounds.miny + _height / 2.0, b->bounds.minz + _depth  / 2.0);
				__vector p0 = objects[index].p0 - o;
				__vector p1 = objects[index].p1 - o;
				__vector p2 = objects[index].p2 - o;
				__vector e0 = p1 - p0;
				__vector e1 = p2 - p0;
				__vector e2 = p2 - p1;
				double hw = _width/2.0, hh = _height/2.0, hd = _depth/2.0;

				// triangle normal
				if (hyperplaneSeparation(objects[index].n, p0, p1, p2, hw, hh, hd)) continue;

				// 9 cross products
				if (hyperplaneSeparation(e0.cross(__vector(1.0,0.0,0.0)), p0, p1, p2, hw, hh, hd)) continue;
				if (hyperplaneSeparation(e0.cross(__vector(0.0,1.0,0.0)), p0, p1, p2, hw, hh, hd)) continue;
				if (hyperplaneSeparation(e0.cross(__vector(0.0,0.0,1.0)), p0, p1, p2, hw, hh, hd)) continue;
				if (hyperplaneSeparation(e1.cross(__vector(1.0,0.0,0.0)), p0, p1, p2, hw, hh, hd)) continue;
				if (hyperplaneSeparation(e1.cross(__vector(0.0,1.0,0.0)), p0, p1, p2, hw, hh, hd)) continue;
				if (hyperplaneSeparation(e1.cross(__vector(0.0,0.0,1.0)), p0, p1, p2, hw, hh, hd)) continue;
				if (hyperplaneSeparation(e2.cross(__vector(1.0,0.0,0.0)), p0, p1, p2, hw, hh, hd)) continue;
				if (hyperplaneSeparation(e2.cross(__vector(0.0,1.0,0.0)), p0, p1, p2, hw, hh, hd)) continue;
				if (hyperplaneSeparation(e2.cross(__vector(0.0,0.0,1.0)), p0, p1, p2, hw, hh, hd)) continue;

				b->objects.push_back(index);

				if (kdType == KD_EVEN) {
				  
				} else if (kdType == KD_MEDIAN) {
					if (axis == Z) {
						list.push_back(objects[index].p0.z); list.push_back(objects[index].p1.z); list.push_back(objects[index].p2.z);
					} else if (axis == Y) {
						list.push_back(objects[index].p0.y); list.push_back(objects[index].p1.y); list.push_back(objects[index].p2.y);
					} else if (axis == X) {
						list.push_back(objects[index].p0.x); list.push_back(objects[index].p1.x); list.push_back(objects[index].p2.x);
					}
				} else if (kdType == KD_SAH) {
					sp.value = MIN(objects[index].p0.z,MIN(objects[index].p1.z,objects[index].p2.z)); sp.event = PRIMITIVE_START; listz.push_back(sp);
					sp.value = MAX(objects[index].p0.z,MAX(objects[index].p1.z,objects[index].p2.z)); sp.event = PRIMITIVE_END;   listz.push_back(sp);
					sp.value = MIN(objects[index].p0.y,MIN(objects[index].p1.y,objects[index].p2.y)); sp.event = PRIMITIVE_START; listy.push_back(sp);
					sp.value = MAX(objects[index].p0.y,MAX(objects[index].p1.y,objects[index].p2.y)); sp.event = PRIMITIVE_END;   listy.push_back(sp);
					sp.value = MIN(objects[index].p0.x,MIN(objects[index].p1.x,objects[index].p2.x)); sp.event = PRIMITIVE_START; listx.push_back(sp);
					sp.value = MAX(objects[index].p0.x,MAX(objects[index].p1.x,objects[index].p2.x)); sp.event = PRIMITIVE_END;   listx.push_back(sp);
				}
			} else if (objects[index].type == SPHERE) {
				// minimum bounding box
				double sx = objects[index].center.x - objects[index].radius;
				double sy = objects[index].center.y - objects[index].radius;
				double sz = objects[index].center.z - objects[index].radius;
				double bx = objects[index].center.x + objects[index].radius;
				double by = objects[index].center.y + objects[index].radius;
				double bz = objects[index].center.z + objects[index].radius;
				if (sx > b->bounds.maxx || bx < b->bounds.minx) continue;
				if (sy > b->bounds.maxy || by < b->bounds.miny) continue;
				if (sz > b->bounds.maxz || bz < b->bounds.minz) continue;
				
				b->objects.push_back(index);

				if (kdType == KD_EVEN) {
				  
				} else if (kdType == KD_MEDIAN) {
					if (axis == Z) {
						list.push_back(objects[index].center.z-objects[index].radius); list.push_back(objects[index].center.z+objects[index].radius);
					} else if (axis == Y) {
						list.push_back(objects[index].center.y-objects[index].radius); list.push_back(objects[index].center.y+objects[index].radius);
					} else if (axis == X) {
						list.push_back(objects[index].center.x-objects[index].radius); list.push_back(objects[index].center.x+objects[index].radius);
					}
				} else if (kdType == KD_SAH) {
					sp.value = objects[index].center.z-objects[index].radius; sp.event = PRIMITIVE_START; listz.push_back(sp);
					sp.value = objects[index].center.z+objects[index].radius; sp.event = PRIMITIVE_END;   listz.push_back(sp);
					sp.value = objects[index].center.y-objects[index].radius; sp.event = PRIMITIVE_START; listy.push_back(sp);
					sp.value = objects[index].center.y+objects[index].radius; sp.event = PRIMITIVE_END;   listy.push_back(sp);
					sp.value = objects[index].center.x-objects[index].radius; sp.event = PRIMITIVE_START; listx.push_back(sp);
					sp.value = objects[index].center.x+objects[index].radius; sp.event = PRIMITIVE_END;   listx.push_back(sp);
				}
			}
		}
	}

	if (depth < max_depth) { // termination critera
		if (b->objects.size() > max_objects) { // termination critera
			if (kdType == KD_EVEN) { // split the node in half

				__bounds bounds0 = bounds, bounds1 = bounds;

				b->axis = axis;
				if (axis == Z) {
					b->axis = Z;
					bounds0.maxz = bounds1.minz = (bounds.minz + bounds.maxz)/2.0;
					b->child0 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds0, depth + 1);
					b->child1 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds1, depth + 1);
				} else if (axis == Y) {
					b->axis = Y;
					bounds0.maxy = bounds1.miny = (bounds.miny + bounds.maxy)/2.0;
					b->child0 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds0, depth + 1);
					b->child1 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds1, depth + 1);
				} else if (axis == X) {
					b->axis = X;
					bounds0.maxx = bounds1.minx = (bounds.minx + bounds.maxx)/2.0;
					b->child0 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds0, depth + 1);
					b->child1 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds1, depth + 1);
				}

			} else if (kdType == KD_MEDIAN) { // split the node at the object median

				std::sort(list.begin(), list.end());

				double median = (list.size() % 2) ? list[(list.size() - 1) / 2] : ((list[list.size()/2] + list[list.size()/2-1])/2.0);

				__bounds bounds0 = bounds, bounds1 = bounds;

				b->axis = axis;
				if (axis == Z) {
					if (median > bounds0.minz && median < bounds1.maxz) { // ensure the median splits the node appropriately
						bounds0.maxz = bounds1.minz = median;
						b->child0 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds0, depth + 1);
						b->child1 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds1, depth + 1);
					}
				} else if (axis == Y) {
					if (median > bounds0.miny && median < bounds1.maxy) { // ensure the median splits the node appropriately
						bounds0.maxy = bounds1.miny = median;
						b->child0 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds0, depth + 1);
						b->child1 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds1, depth + 1);
					}
				} else if (axis == X) {
					if (median > bounds0.minx && median < bounds1.maxx) { // ensure the median splits the node appropriately
						bounds0.maxx = bounds1.minx = median;
						b->child0 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds0, depth + 1);
						b->child1 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds1, depth + 1);
					}
				}

			} else if (kdType == KD_SAH) {

				std::sort(listz.begin(), listz.end(), __splitting_plane_sort);
				std::sort(listy.begin(), listy.end(), __splitting_plane_sort);
				std::sort(listx.begin(), listx.end(), __splitting_plane_sort);

				// find the best split plane using surface area heuristics
				double traversal_cost = 1.0;
				double intersection_cost = 1.0;

				double best_cost = 1e8, cost, nosplit_cost = b->objects.size() * intersection_cost; // best cost, current cost, and the cost of not splitting the node
				double split_position;
				double _width  = bounds.maxx - bounds.minx;
				double _height = bounds.maxy - bounds.miny;
				double _depth  = bounds.maxz - bounds.minz;
				
				double surface_area = 2.0 * (_width * _height + _height * _depth + _depth * _width), // surface area of node
				       templ, tempr, surface_area_l, surface_area_r; // temp values for dimension split and surface area of child nodes
				int nl, nr; // number of primitives in child nodes
				char flag; // we are including primitives lying on a split plane in both children.. this should prevent a decrease in nr on the first pass

				nl = 0; nr = b->objects.size(); flag = 0;
				for (int i = 0; i < listz.size(); i++) {
					if (listz[i].event == PRIMITIVE_START) { if (flag == 1) { nr--; flag = 0; } nl++; }
					if (listz[i].event == PRIMITIVE_END)   { if (flag == 1) { nr--; } flag = 1; }

					if (listz[i].value <= bounds.minz || listz[i].value >= bounds.maxz) continue; // ensure the split position is a proper one

					templ = listz[i].value - bounds.minz;
					tempr = bounds.maxz - listz[i].value;
					surface_area_l = 2.0 * (_width * _height + _height * templ + templ * _width);
					surface_area_r = 2.0 * (_width * _height + _height * tempr + tempr * _width);

					cost = traversal_cost + surface_area_l / surface_area * nl * intersection_cost + surface_area_r / surface_area * nr * intersection_cost;
					if (cost < best_cost) {
						best_cost = cost;
						axis = Z;
						split_position = listz[i].value;
					}
				}
				nl = 0; nr = b->objects.size(); flag = 0;
				for (int i = 0; i < listy.size(); i++) {
					if (listy[i].event == PRIMITIVE_START) { if (flag == 1) { nr--; flag = 0; } nl++; }
					if (listy[i].event == PRIMITIVE_END)   { if (flag == 1) { nr--; } flag = 1; }

					if (listy[i].value <= bounds.miny || listy[i].value >= bounds.maxy) continue; // ensure the split position is a proper one

					templ = listy[i].value - bounds.miny;
					tempr = bounds.maxy - listy[i].value;
					surface_area_l = 2.0 * (_width * templ + templ * _depth + _depth * _width);
					surface_area_r = 2.0 * (_width * tempr + tempr * _depth + _depth * _width);

					cost = traversal_cost + surface_area_l / surface_area * nl * intersection_cost + surface_area_r / surface_area * nr * intersection_cost;
					if (cost < best_cost) {
						best_cost = cost;
						axis = Y;
						split_position = listy[i].value;
					}
				}
				nl = 0; nr = b->objects.size(); flag = 0;
				for (int i = 0; i < listx.size(); i++) {
					if (listx[i].event == PRIMITIVE_START) { if (flag == 1) { nr--; flag = 0; } nl++; }
					if (listx[i].event == PRIMITIVE_END)   { if (flag == 1) { nr--; } flag = 1; }

					if (listx[i].value <= bounds.minx || listx[i].value >= bounds.maxx) continue; // ensure the split position is a proper one

					templ = listx[i].value - bounds.minx;
					tempr = bounds.maxx - listx[i].value;
					surface_area_l = 2.0 * (templ * _height + _height * _depth + _depth * templ);
					surface_area_r = 2.0 * (tempr * _height + _height * _depth + _depth * tempr);

					cost = traversal_cost + surface_area_l / surface_area * nl * intersection_cost + surface_area_r / surface_area * nr * intersection_cost;
					if (cost < best_cost) {
						best_cost = cost;
						axis = X;
						split_position = listx[i].value;
					}
				}

				// split it
				__bounds bounds0 = bounds, bounds1 = bounds;

				if (best_cost < nosplit_cost) { // termination criteria
					b->axis = axis;
					if (axis == Z) {
						bounds0.maxz = bounds1.minz = split_position;
						b->child0 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds0, depth + 1);
						b->child1 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds1, depth + 1);
					} else if (axis == Y) {
						bounds0.maxy = bounds1.miny = split_position;
						b->child0 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds0, depth + 1);
						b->child1 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds1, depth + 1);
					} else if (axis == X) {
						bounds0.maxx = bounds1.minx = split_position;
						b->child0 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds0, depth + 1);
						b->child1 = buildTree(kdType, node, b, max_depth, max_objects, objects, object_count, bounds1, depth + 1);
					}
				}

			}
		}
	}

	if (b->child0 == 0) {
		b->leaf_id = leaf_id++;
		node.leaf_size++;
		node.max_leaf_objects = MAX(node.max_leaf_objects, b->objects.size());
	}
	
	node.size++;
	
	return b;
}

void copyKdTree(__node& node, _box* root) {
	int id = root->id;

	node.id[id] = id;
	node.leaf_id[id] = root->leaf_id;
	node.parent[id] = root->parent ? root->parent->id : -1;
	node.child0[id] = root->child0 ? root->child0->id : -1;
	node.child1[id] = root->child1 ? root->child1->id : -1;
	node.bounds[id] = root->bounds;
	node.axis[id] = root->axis;
	node.object_count[id] = root->objects.size();

	if (node.leaf_id[id] > -1) {
		int offset = node.leaf_id[id] * node.max_leaf_objects;
		for (int i = 0; i < root->objects.size(); i++)
			node.objects[offset + i] = root->objects[i];
	} else {
		copyKdTree(node, root->child0);
		copyKdTree(node, root->child1);
	}
	
}

__node buildKdTree(unsigned char kdType, unsigned char max_depth, unsigned char max_objects, __object objects[], unsigned int object_count, const __bounds& bounds) {

	// populate and return __node for passing to gpu
	__node node;
	node.size = 0;
	node.leaf_size = 0;
	node.max_leaf_objects = 0;

	_box* root = buildTree(kdType, node, 0, max_depth, max_objects, objects, object_count, bounds, 0);
	
	node.id = (int *)malloc(sizeof(int) * node.size);
	node.leaf_id = (int *)malloc(sizeof(int) * node.size);
	node.parent = (int *)malloc(sizeof(int) * node.size);
	node.child0 = (int *)malloc(sizeof(int) * node.size);
	node.child1 = (int *)malloc(sizeof(int) * node.size);
	node.bounds = (__bounds *)malloc(sizeof(__bounds) * node.size);
	node.axis = (char *)malloc(sizeof(char) * node.size);
	node.object_count = (int *)malloc(sizeof(int) * node.size);
	node.objects = (int *)malloc(sizeof(int) * node.leaf_size * node.max_leaf_objects);

	copyKdTree(node, root);
	
	deleteTree(root);

	return node;
}

void deleteKdTree(__node& node) {
	free(node.id);
	free(node.leaf_id);
	free(node.parent);
	free(node.child0);
	free(node.child1);
	free(node.bounds);
	free(node.axis);
	free(node.object_count);
	free(node.objects);
}

int main(int argc, char* argv[]) {
//	__dimensions dim = {1920,1080};
	__dimensions dim = {1280,720};
//	__dimensions dim = {640,360};
	const int BPP = 24;
	const int MAX_BOUNCES = 8;
	
	cBuffer buf(dim.width, dim.height);

	cKeyboard kb;
	cTimer t0, t1; double elapsed0, elapsed1;
  
	timespec current;
	clock_gettime(CLOCK_REALTIME, &current);
	srand(current.tv_sec + current.tv_nsec / 1000000000.0);
	
	bool active = true;
	unsigned int samples = 1;

	SDL_Init(SDL_INIT_EVERYTHING);
	SDL_Surface *screen = mySDLInit(dim.width, dim.height, BPP, false);
	SDL_Event event;

	std::vector<cObject*> objects;

// dragon scene ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	cObj cube("dragon_smooth.obj");
	cObj cube("cube_uv.obj");
	double scale = 8;
	vector3 translate(-7,0,3);

	//objects.push_back(new cPlane(vector3( 0, -1,  0), vector3( -16,  128,              -16), vector3(  1,  1,  1)*0.0,   vector3(1,1,1)*1.5, DIFFUSE,  2, vector3(0,0,80)*10, vector3(80,0,0)*10)); // top
	objects.push_back(new cPlane(vector3( 1,  0,  0), vector3(  -20,  cube.miny*scale,    -25), vector3(  1,  0.5,  0.5)*0.9,   vector3(),          DIFFUSE, 7, vector3(0,10, 0), vector3( 0,0,10))); // bottom
	objects.push_back(new cPlane(vector3(-1,  0,  0), vector3(   20,  cube.miny*scale,    -20), vector3(  0.5,  0.5,  1)*0.9,   vector3(),          DIFFUSE, 7, vector3(0,10, 0), vector3( 0,0,10))); // bottom
	objects.push_back(new cPlane(vector3( 0,  0,  1), vector3(  -20,  cube.miny*scale,    -20), vector3(  1,  1,  1)*0.9,   vector3(),          DIFFUSE, -1, vector3(1,0, 0), vector3( 0,1,0))); // bottom
	objects.push_back(new cPlane(vector3( 0,  0, -1), vector3(  -20,  cube.miny*scale,     45), vector3(  1,  1,  1)*0.9,   vector3(),          DIFFUSE, -1, vector3(1,0, 0), vector3( 0,1,0))); // bottom
	objects.push_back(new cPlane(vector3( 0,  1,  0), vector3(  -20,  cube.miny*scale,    -20), vector3(  1,  1,  1)*0.9,   vector3(),          DIFFUSE,  7, vector3(10,0, 0), vector3( 0,0,10))); // bottom
	objects.push_back(new cPlane(vector3( 0, -1,  0), vector3(  -20,  cube.miny*scale+20, -20), vector3(  1,  1,  1)*0.9,   vector3(1,1,1)*0.0,    DIFFUSE, -1, vector3(1,0, 0), vector3( 0,0,1))); // bottom

	double light = 4;
	objects.push_back(new cTriangle(vector3(-12,cube.miny*scale+19.9,-12), vector3(-12,cube.miny*scale+19.9, -2), vector3(-2,cube.miny*scale+19.9,-2),
					vector3(0,-1,0), vector3(0,-1,0), vector3(0,-1,0), vector3(1,1,1), vector3(1,1,1)*light, DIFFUSE, -1, vector3(), vector3(), vector3()));
	objects.push_back(new cTriangle(vector3(-12,cube.miny*scale+19.9,-12), vector3( -2,cube.miny*scale+19.9, -12), vector3(-2,cube.miny*scale+19.9, -2),
					vector3(0,-1,0), vector3(0,-1,0), vector3(0,-1,0), vector3(1,1,1), vector3(1,1,1)*light, DIFFUSE, -1, vector3(), vector3(), vector3()));

	objects.push_back(new cTriangle(vector3(2,cube.miny*scale+19.9,-12), vector3(2,cube.miny*scale+19.9, -2), vector3(12,cube.miny*scale+19.9,-2),
					vector3(0,-1,0), vector3(0,-1,0), vector3(0,-1,0), vector3(1,1,1), vector3(1,1,1)*light, DIFFUSE, -1, vector3(), vector3(), vector3()));
	objects.push_back(new cTriangle(vector3(2,cube.miny*scale+19.9,-12), vector3(12,cube.miny*scale+19.9,-12), vector3(12,cube.miny*scale+19.9,-2),
					vector3(0,-1,0), vector3(0,-1,0), vector3(0,-1,0), vector3(1,1,1), vector3(1,1,1)*light, DIFFUSE, -1, vector3(), vector3(), vector3()));

	objects.push_back(new cTriangle(vector3(-12,cube.miny*scale+19.9,2), vector3(-12,cube.miny*scale+19.9, 12), vector3(-2,cube.miny*scale+19.9,12),
					vector3(0,-1,0), vector3(0,-1,0), vector3(0,-1,0), vector3(1,1,1), vector3(1,1,1)*light, DIFFUSE, -1, vector3(), vector3(), vector3()));
	objects.push_back(new cTriangle(vector3(-12,cube.miny*scale+19.9,2), vector3(-2,cube.miny*scale+19.9,2), vector3(-2,cube.miny*scale+19.9,12),
					vector3(0,-1,0), vector3(0,-1,0), vector3(0,-1,0), vector3(1,1,1), vector3(1,1,1)*light, DIFFUSE, -1, vector3(), vector3(), vector3()));

	objects.push_back(new cTriangle(vector3(2,cube.miny*scale+19.9,2), vector3(2,cube.miny*scale+19.9, 12), vector3(12,cube.miny*scale+19.9,12),
					vector3(0,-1,0), vector3(0,-1,0), vector3(0,-1,0), vector3(1,1,1), vector3(1,1,1)*light, DIFFUSE, -1, vector3(), vector3(), vector3()));
	objects.push_back(new cTriangle(vector3(2,cube.miny*scale+19.9,2), vector3(12,cube.miny*scale+19.9,2), vector3(12,cube.miny*scale+19.9,12),
					vector3(0,-1,0), vector3(0,-1,0), vector3(0,-1,0), vector3(1,1,1), vector3(1,1,1)*light, DIFFUSE, -1, vector3(), vector3(), vector3()));

	for (int i = 0; i < MIN(1000000,cube.faces.size()); i++) {
		objects.push_back(new cTriangle(vector3(cube.vertices[cube.faces[i].vertex[0]].v[0], cube.vertices[cube.faces[i].vertex[0]].v[1], cube.vertices[cube.faces[i].vertex[0]].v[2]) * scale + translate,
						vector3(cube.vertices[cube.faces[i].vertex[1]].v[0], cube.vertices[cube.faces[i].vertex[1]].v[1], cube.vertices[cube.faces[i].vertex[1]].v[2]) * scale + translate,
						vector3(cube.vertices[cube.faces[i].vertex[2]].v[0], cube.vertices[cube.faces[i].vertex[2]].v[1], cube.vertices[cube.faces[i].vertex[2]].v[2]) * scale + translate,
						vector3(cube.normals[cube.faces[i].normal[0]].v[0], cube.normals[cube.faces[i].normal[0]].v[1], cube.normals[cube.faces[i].normal[0]].v[2]).unit(),
						vector3(cube.normals[cube.faces[i].normal[1]].v[0], cube.normals[cube.faces[i].normal[1]].v[1], cube.normals[cube.faces[i].normal[1]].v[2]).unit(),
						vector3(cube.normals[cube.faces[i].normal[2]].v[0], cube.normals[cube.faces[i].normal[2]].v[1], cube.normals[cube.faces[i].normal[2]].v[2]).unit(),
						vector3(1,1,1)*0.9, vector3(0,0,0), DIFFUSE, 6,
						vector3(cube.texcoords[cube.faces[i].texture[0]].v[0], cube.texcoords[cube.faces[i].texture[0]].v[1], cube.texcoords[cube.faces[i].texture[0]].v[2]),
						vector3(cube.texcoords[cube.faces[i].texture[1]].v[0], cube.texcoords[cube.faces[i].texture[1]].v[1], cube.texcoords[cube.faces[i].texture[1]].v[2]),
						vector3(cube.texcoords[cube.faces[i].texture[2]].v[0], cube.texcoords[cube.faces[i].texture[2]].v[1], cube.texcoords[cube.faces[i].texture[2]].v[2])));
	}
	double angle_x = 0.2, angle_y = 0.3, angle_z = -0.4;
	vector3 e0,e1;
/*
	double rotx[] = {       1,             0,            0, 0,
	                        0,  cos(angle_x), sin(angle_x), 0,
	                        0, -sin(angle_x), cos(angle_x), 0,
				0,             0,            0, 1 };
	double roty[] = {       cos(angle_y), 0, -sin(angle_y), 0,
				           0, 1,             0, 0,
	                        sin(angle_y), 0,  cos(angle_y), 0,
				           0, 0,             0, 1 };
	double rotz[] = {        cos(angle_z), sin(angle_z), 0, 0,
	                        -sin(angle_z), cos(angle_z), 0, 0,
				            0,            0, 1, 0,
	                                    0,            0, 0, 1 };
	cMatrix matx(4, 4, rotx);
	cMatrix maty(4, 4, roty);
	cMatrix matz(4, 4, rotz);
	e0 = matz.mult(maty.mult(matx.mult(vector3(0,1,0))));
	e1 = matz.mult(maty.mult(matx.mult(vector3(1,0,0))));

	objects.push_back(new cSphere(vector3(-5,scale*cube.miny+3,-3), 3, vector3(1,1,1), vector3(), DIFFUSE, 3, e0, e1));
*/
	angle_x = 0.4; angle_y = 0.3; angle_z = 0.2;
	double _rotx[] = {       1,             0,            0, 0,
	                        0,  cos(angle_x), sin(angle_x), 0,
	                        0, -sin(angle_x), cos(angle_x), 0,
				0,             0,            0, 1 };
	double _roty[] = {       cos(angle_y), 0, -sin(angle_y), 0,
				           0, 1,             0, 0,
	                        sin(angle_y), 0,  cos(angle_y), 0,
				           0, 0,             0, 1 };
	double _rotz[] = {        cos(angle_z), sin(angle_z), 0, 0,
	                        -sin(angle_z), cos(angle_z), 0, 0,
				            0,            0, 1, 0,
	                                    0,            0, 0, 1 };
	cMatrix _matx(4, 4, _rotx);
	cMatrix _maty(4, 4, _roty);
	cMatrix _matz(4, 4, _rotz);
	e0 = _matz.mult(_maty.mult(_matx.mult(vector3(0,1,0))));
	e1 = _matz.mult(_maty.mult(_matx.mult(vector3(1,0,0))));

	objects.push_back(new cSphere(vector3(5,scale*cube.miny+6,-1), 6, vector3(1,1,1), vector3(1,1,1)*0.4, DIFFUSE, 4, e0, e1));

	angle_x = M_PI/4.0; angle_y = M_PI+0.5; angle_z = 0;
	double __rotx[] = {       1,             0,            0, 0,
	                        0,  cos(angle_x), sin(angle_x), 0,
	                        0, -sin(angle_x), cos(angle_x), 0,
				0,             0,            0, 1 };
	double __roty[] = {       cos(angle_y), 0, -sin(angle_y), 0,
				           0, 1,             0, 0,
	                        sin(angle_y), 0,  cos(angle_y), 0,
				           0, 0,             0, 1 };
	double __rotz[] = {        cos(angle_z), sin(angle_z), 0, 0,
	                        -sin(angle_z), cos(angle_z), 0, 0,
				            0,            0, 1, 0,
	                                    0,            0, 0, 1 };
	cMatrix __matx(4, 4, __rotx);
	cMatrix __maty(4, 4, __roty);
	cMatrix __matz(4, 4, __rotz);
	e0 = __matz.mult(__maty.mult(__matx.mult(vector3(0,1,0))));
	e1 = __matz.mult(__maty.mult(__matx.mult(vector3(1,0,0))));

	objects.push_back(new cSphere(vector3(-7,scale*cube.maxy+4,3), 4, vector3(1,1,1), vector3(), SPECULAR, -1, e0, e1));
//	objects.push_back(new cSphere(vector3(-7,scale*cube.miny+6,2), 4, vector3(1,1,1), vector3(), DIFFUSE, 5, e0, e1));

	objects.push_back(new cSphere(vector3(12,scale*cube.miny+5,7), 5, vector3(1,1,1), vector3(1,1,1)*0.0, REFRACTIVE, -1, e0, e1));

//	objects.push_back(new cSphere(vector3(4,scale*cube.miny+6,-2),  2, vector3(1,1,1), vector3(1,1,1)*3, DIFFUSE, -1, vector3(0.0,1.0,0.0), vector3(1.0,0.0,0.0)));
//	objects.push_back(new cSphere(vector3(-10,scale*cube.miny+12,0),  5, vector3(1,1,1), vector3(1,1,1)*3, DIFFUSE, -1, vector3(0.0,1.0,0.0), vector3(1.0,0.0,0.0)));
//	objects.push_back(new cSphere(vector3(3,scale*cube.miny+1,3), 1, vector3(1,1,1), vector3(), REFRACTIVE, -1));
//	objects.push_back(new cSphere(vector3(0,scale*cube.miny+32,0), 8, vector3(1,1,1), vector3(1,1,1)*2, DIFFUSE, -1));


	// set up the camera matrix and apply it to the primitives
	int frame = 0;

//	vector3 origin(-8.5,0,16);
	vector3 origin(7,8,40);
	vector3 target(0,cube.miny*scale+9,0);
	vector3 forward = (target-origin).unit();
	vector3 up = vector3(0.05, 1, 0.0).unit();
	applyCamera(origin, forward, up, objects);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	__buffers buffers;
	buffers.object_count = objects.size();
	buffers.objects_host = (__object *)malloc(sizeof(__object) * objects.size());

	__bounds bounds = {1e8,1e8,1e8,-1e8,-1e8,-1e8};
	copyObjects(objects, buffers.objects_host, bounds);

	__node node = buildKdTree(KD_SAH, 18, 1, buffers.objects_host, objects.size(), bounds);
	std::cerr << "number of nodes: " << node.size << "   number of child nodes: " << node.leaf_size << "  max objects in a child: " << node.max_leaf_objects << std::endl;


	// set up the camera model
	__camera cmodel;
	cmodel.focal_length = 30.0 / 1000.0;
	cmodel.aperture = 0.25;//1.0;
	cmodel.focal_distance = sqrt(pow((bounds.minx+bounds.maxx)/2.0,2.0)+pow((bounds.miny+bounds.maxy)/2.0,2.0)+pow((bounds.minz+bounds.maxz)/2.0,2.0)); // focal point at center of objects
	cmodel.image_distance = 1.0 / (1.0 / cmodel.focal_length - 1.0 / cmodel.focal_distance); // thin lens
	cmodel.aperture_diameter = cmodel.focal_length / cmodel.aperture;

	__textures tex;
	tex.count = 0;
	tex.texture[tex.count] = loadImage("grid.png",    tex.width[tex.count], tex.height[tex.count], tex.bpp[tex.count]); tex.count++;
	tex.texture[tex.count] = loadImage("tile.jpg",    tex.width[tex.count], tex.height[tex.count], tex.bpp[tex.count]); tex.count++;
	tex.texture[tex.count] = loadImage("clouds.jpg",  tex.width[tex.count], tex.height[tex.count], tex.bpp[tex.count]); tex.count++;
	tex.texture[tex.count] = loadImage("9.png",       tex.width[tex.count], tex.height[tex.count], tex.bpp[tex.count]); tex.count++;
	tex.texture[tex.count] = loadImage("uv.png",      tex.width[tex.count], tex.height[tex.count], tex.bpp[tex.count]); tex.count++;
	tex.texture[tex.count] = loadImage("earth.png",   tex.width[tex.count], tex.height[tex.count], tex.bpp[tex.count]); tex.count++;
	tex.texture[tex.count] = loadImage("uv_save.jpg", tex.width[tex.count], tex.height[tex.count], tex.bpp[tex.count]); tex.count++;
	tex.texture[tex.count] = loadImage("test.jpg",    tex.width[tex.count], tex.height[tex.count], tex.bpp[tex.count]); tex.count++;



	initializePathTracer(node, buffers, MAX_BOUNCES, dim, bounds, tex);
        setupPathTracer(node, buffers, dim, bounds);


// this section to render progress
	// shaders
	GLuint glProgram, glShaderV, glShaderF;
	GLint vertex, texCoord, Projection, View, Model, textureSample;
	createProgram(glProgram, glShaderV, glShaderF, "src/vertex.sh", "src/fragment.sh");
	vertex        = glGetAttribLocation(glProgram, "vertex");
	texCoord      = glGetAttribLocation(glProgram, "texCoord");
	Projection    = glGetUniformLocation(glProgram, "Projection");
	View          = glGetUniformLocation(glProgram, "View");
	Model         = glGetUniformLocation(glProgram, "Model");
	textureSample = glGetUniformLocation(glProgram, "textureSample");

	// vertices
	const int vertex_count = 4;
	vertex_ vertices[vertex_count];
	vertices[0].x  = 0;         vertices[0].y  = 0;          vertices[0].z = 0;     vertices[0].tx = 0;     vertices[0].ty = 0;	
	vertices[1].x  = dim.width; vertices[1].y  = 0;          vertices[1].z = 0;     vertices[1].tx = 1;     vertices[1].ty = 0;	
	vertices[2].x  = dim.width; vertices[2].y  = dim.height; vertices[2].z = 0;     vertices[2].tx = 1;     vertices[2].ty = 1;
	vertices[3].x  = 0;         vertices[3].y  = dim.height; vertices[3].z = 0;     vertices[3].tx = 0;     vertices[3].ty = 1;
	GLuint vbo_vertices;
	glGenBuffers(1, &vbo_vertices);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertex_) * vertex_count, vertices, GL_DYNAMIC_DRAW);

	// indices
	const int indices_count = 6;
	unsigned int indices[indices_count] = {0, 1, 2, 2, 3, 0};
	GLuint vbo_indices;
	glGenBuffers(1, &vbo_indices);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_indices);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * indices_count, indices, GL_STATIC_DRAW);

	// create a texture for rendering the result
	GLuint texture;
	setupTexture(texture);

	glm::mat4 projection = glm::ortho(0.0f, (float)dim.width, (float)dim.height, 0.0f, -5.0f, 5.0f); 
	glm::mat4 view       = glm::mat4(1.0f);
	glm::mat4 model      = glm::mat4(1.0f);

	
// this section to render the kd tree
	// shaders
	GLuint glProgramc, glShaderVc, glShaderFc;
	GLint vertexc, colorc, Projectionc, Viewc, Modelc;
	createProgram(glProgramc, glShaderVc, glShaderFc, "src/vertexc.sh", "src/fragmentc.sh");
	vertexc        = glGetAttribLocation(glProgramc, "vertex");
	colorc         = glGetUniformLocation(glProgramc, "color");
	Projectionc    = glGetUniformLocation(glProgramc, "Projection");
	Viewc          = glGetUniformLocation(glProgramc, "View");
	Modelc         = glGetUniformLocation(glProgramc, "Model");

	// vertices
	int vertex_countc = 24 * node.size;
	vertex_ *verticesc = new vertex_[vertex_countc];
	unsigned int _c = 0;
	for (int i = 0; i < node.size; i++) {
		verticesc[_c].x = node.bounds[i].minx; verticesc[_c].y = node.bounds[i].miny; verticesc[_c].z = node.bounds[i].minz; _c++;
		verticesc[_c].x = node.bounds[i].maxx; verticesc[_c].y = node.bounds[i].miny; verticesc[_c].z = node.bounds[i].minz; _c++;
		verticesc[_c].x = node.bounds[i].minx; verticesc[_c].y = node.bounds[i].miny; verticesc[_c].z = node.bounds[i].minz; _c++;
		verticesc[_c].x = node.bounds[i].minx; verticesc[_c].y = node.bounds[i].maxy; verticesc[_c].z = node.bounds[i].minz; _c++;
		verticesc[_c].x = node.bounds[i].minx; verticesc[_c].y = node.bounds[i].miny; verticesc[_c].z = node.bounds[i].minz; _c++;
		verticesc[_c].x = node.bounds[i].minx; verticesc[_c].y = node.bounds[i].miny; verticesc[_c].z = node.bounds[i].maxz; _c++;

		verticesc[_c].x = node.bounds[i].minx; verticesc[_c].y = node.bounds[i].maxy; verticesc[_c].z = node.bounds[i].maxz; _c++;
		verticesc[_c].x = node.bounds[i].minx; verticesc[_c].y = node.bounds[i].miny; verticesc[_c].z = node.bounds[i].maxz; _c++;
		verticesc[_c].x = node.bounds[i].minx; verticesc[_c].y = node.bounds[i].maxy; verticesc[_c].z = node.bounds[i].maxz; _c++;
		verticesc[_c].x = node.bounds[i].maxx; verticesc[_c].y = node.bounds[i].maxy; verticesc[_c].z = node.bounds[i].maxz; _c++;
		verticesc[_c].x = node.bounds[i].minx; verticesc[_c].y = node.bounds[i].maxy; verticesc[_c].z = node.bounds[i].maxz; _c++;
		verticesc[_c].x = node.bounds[i].minx; verticesc[_c].y = node.bounds[i].maxy; verticesc[_c].z = node.bounds[i].minz; _c++;

		verticesc[_c].x = node.bounds[i].maxx; verticesc[_c].y = node.bounds[i].maxy; verticesc[_c].z = node.bounds[i].minz; _c++;
		verticesc[_c].x = node.bounds[i].minx; verticesc[_c].y = node.bounds[i].maxy; verticesc[_c].z = node.bounds[i].minz; _c++;
		verticesc[_c].x = node.bounds[i].maxx; verticesc[_c].y = node.bounds[i].maxy; verticesc[_c].z = node.bounds[i].minz; _c++;
		verticesc[_c].x = node.bounds[i].maxx; verticesc[_c].y = node.bounds[i].maxy; verticesc[_c].z = node.bounds[i].maxz; _c++;
		verticesc[_c].x = node.bounds[i].maxx; verticesc[_c].y = node.bounds[i].maxy; verticesc[_c].z = node.bounds[i].minz; _c++;
		verticesc[_c].x = node.bounds[i].maxx; verticesc[_c].y = node.bounds[i].miny; verticesc[_c].z = node.bounds[i].minz; _c++;

		verticesc[_c].x = node.bounds[i].maxx; verticesc[_c].y = node.bounds[i].miny; verticesc[_c].z = node.bounds[i].maxz; _c++;
		verticesc[_c].x = node.bounds[i].maxx; verticesc[_c].y = node.bounds[i].maxy; verticesc[_c].z = node.bounds[i].maxz; _c++;
		verticesc[_c].x = node.bounds[i].maxx; verticesc[_c].y = node.bounds[i].miny; verticesc[_c].z = node.bounds[i].maxz; _c++;
		verticesc[_c].x = node.bounds[i].minx; verticesc[_c].y = node.bounds[i].miny; verticesc[_c].z = node.bounds[i].maxz; _c++;
		verticesc[_c].x = node.bounds[i].maxx; verticesc[_c].y = node.bounds[i].miny; verticesc[_c].z = node.bounds[i].maxz; _c++;
		verticesc[_c].x = node.bounds[i].maxx; verticesc[_c].y = node.bounds[i].miny; verticesc[_c].z = node.bounds[i].minz; _c++;
	}
	GLuint vbo_verticesc;
	glGenBuffers(1, &vbo_verticesc);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_verticesc);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertex_) * vertex_countc, verticesc, GL_DYNAMIC_DRAW);
	delete [] verticesc;

	// indices
	int indices_countc = 24 * node.size;
	unsigned int *indicesc = new unsigned int[indices_countc];
	unsigned int _d = 0;
	for (int i = 0; i < _c; i++) indicesc[i] = i;
	GLuint vbo_indicesc;
	glGenBuffers(1, &vbo_indicesc);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_indicesc);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * indices_countc, indicesc, GL_STATIC_DRAW);
	delete [] indicesc;

	glm::mat4 projectionc = glm::perspective(2.0f*180.0f/3.14159265359f*(float)atan(0.036f*dim.height/dim.width/(2.0f*cmodel.image_distance)), dim.width/(float)dim.height, 1.0f, 1000.0f);
	glm::mat4 viewc       = glm::mat4(1.0f);
	glm::mat4 modelc      = glm::mat4(1.0f);



	__offsets offset = {0,0};
	bool sleep = true, grab = false;
	while (active) {

		while (SDL_PollEvent(&event)) {
			switch (event.type) {
			    case SDL_QUIT: active = false; break;
			}
		}

		elapsed0 = t0.elapsed(true);
		elapsed1 = t1.elapsed(false);

		if (kb.getKeyState(KEY_SPACE))    { sleep = true;  }
		if (kb.getKeyState(KEY_G))        { sleep = false; }
		if (kb.getKeyState(KEY_F))        { grab  = true;  }

		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

		if (sleep) {
			glEnable(GL_DEPTH_TEST);
			glUseProgram(glProgramc);

			// Project, View, and Model uniforms
			glUniform4f(colorc, 1.0,1.0,1.0,1.0);
			glUniformMatrix4fv(Projectionc, 1, GL_FALSE, glm::value_ptr(projectionc));
			viewc = glm::lookAt(glm::vec3(0.0f, 0.0f,  0.0f),
					    glm::vec3(0.0f, 0.0f, -1.0f),
					    glm::vec3(0.0f, 1.0f,  0.0f));
			glUniformMatrix4fv(Viewc, 1, GL_FALSE, glm::value_ptr(viewc));
			modelc = glm::mat4(1.0);
			glUniformMatrix4fv(Modelc, 1, GL_FALSE, glm::value_ptr(modelc));

			// vertex and texture attributes
			glBindBuffer(GL_ARRAY_BUFFER, vbo_verticesc);
			glDisableVertexAttribArray(vertex);
			glDisableVertexAttribArray(texture);
			glEnableVertexAttribArray(vertexc);
			glVertexAttribPointer(vertexc, 3, GL_DOUBLE, GL_FALSE, sizeof(vertex_), 0);
			
			// draw it
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_indicesc);
			glDrawElements(GL_LINES, indices_countc, GL_UNSIGNED_INT, 0);
			SDL_GL_SwapBuffers();

			if (grab) {
				buf.save(++frame);
				grab = false;
			}
		} else {
			runPathTracer(node, buffers, MAX_BOUNCES, dim, offset, samples, cmodel, tex);

			glDisable(GL_DEPTH_TEST);
			glUseProgram(glProgram);

			// texture
			glUniform1i(textureSample, 0);
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, texture);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, dim.width, dim.height, 0, GL_RGB, GL_UNSIGNED_BYTE, buffers.char_host);

			// Project, View, and Model uniforms
			glUniformMatrix4fv(Projection, 1, GL_FALSE, glm::value_ptr(projection));
			glUniformMatrix4fv(View, 1, GL_FALSE, glm::value_ptr(view));
			//model = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f,0.0f,-2.0f));
			glUniformMatrix4fv(Model, 1, GL_FALSE, glm::value_ptr(model));

			// vertex and texture attributes
			glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);
			glDisableVertexAttribArray(vertexc);
			glEnableVertexAttribArray(vertex);
			glVertexAttribPointer(vertex, 3, GL_DOUBLE, GL_FALSE, sizeof(vertex_), 0);
			glEnableVertexAttribArray(texture);
			glVertexAttribPointer(texCoord, 2, GL_DOUBLE, GL_FALSE, sizeof(vertex_), (char *)NULL + 24);

			// draw it
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_indices);
			glDrawElements(GL_TRIANGLES, indices_count, GL_UNSIGNED_INT, 0);
			SDL_GL_SwapBuffers();

			offset.x += 128;
			if (offset.x >= dim.width) { // (dim.width - 128)) {
				offset.x = 0;
				offset.y += 128;
				if (offset.y >= dim.height) { // (dim.height - 128)) {
					offset.x = offset.y = 0;
					/*if (samples % 50 == 0) {
						std::stringstream s;
						s << "screen" << (int)samples << ".ppm";
						savePPM(buffers.char_host, dim.width, dim.height, BPP, s.str().c_str());
					}*/
					std::cout << "samples: " << samples++ << std::endl;

					grabFrame(buffers, dim);
/* // this section needs to be updated to reflect the kd tree implmentation
					if (samples == 2000000) { // this is used to animate -- camera position changes every n frames

						samples = 1;
						frame += 1;
						std::stringstream s;
						s << "frame" << (frame < 10 ? "000" : (frame < 100 ? "00" : "0")) << frame << ".ppm";
						savePPM(buffers.char_host, dim.width, dim.height, BPP, s.str().c_str());

						vector3 origin(sin(frame/360.0*2.0*3.14159265359)*2.5,  2,  cos(frame/360.0*2.0*3.14159265359)*2.5);
						vector3 forward = (vector3(0,0.5,0)-origin).unit();
						vector3 up(0, 1,  0);
						applyCamera(origin, forward, up, objects);

						bounds = {1e8,1e8,1e8,-1e8,-1e8,-1e8};
						copyObjects(objects, buffers.objects_host, bounds);
// node -- tree needs to be rebuilt
					        setupPathTracer(node, buffers, dim, bounds);
					}
*/
				}
			}
		}

	}

	savePPM(buffers.char_host, dim.width, dim.height, BPP, "screen.ppm");

	for (int i = 0; i < tex.count; i++) delete tex.texture[i];
	
	releasePathTracer(node, buffers, tex);

	deleteKdTree(node);
	free(buffers.objects_host);

	// release the objects
	for (int i = 0; i < objects.size(); i++) delete objects[i];
	
	// release the texture, vertex buffer objects, and shader program.. also shut down SDL
	releaseProgram(glProgramc, glShaderVc, glShaderFc);
	
	deleteTexture(texture);
	glDeleteBuffers(1, &vbo_indices);
	glDeleteBuffers(1, &vbo_vertices);
	releaseProgram(glProgram, glShaderV, glShaderF);

	SDL_Quit();

	return 0;
}