#ifndef COMMON_H
#define COMMON_H

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define CLAMP(x, y, z) ((x)<(y)?(y):((x)>(z)?(z):(x)))

struct __textures {
	int count;
	int width[10], height[10], bpp[10];
	unsigned char *texture[10], *texture_device[10];
};

struct __intersection_t {
	bool intersects;
	double tmin, tmax;
};

struct __dimensions {
	int width, height;
};

struct __offsets {
	int x, y;
};

struct __bounds {
	double minx, miny, minz;
	double maxx, maxy, maxz;
};

struct _box {
	_box *parent, *child0, *child1;
	int id, leaf_id;
	__bounds bounds;
	char axis;
	std::vector<unsigned int> objects;
};

struct __node {
	int size, leaf_size;
	int max_leaf_objects; // for offset into objects property

	int *id, *id_device;
	int *leaf_id, *leaf_id_device;

	// applies to all nodes
	int *parent, *child0, *child1, *parent_device, *child0_device, *child1_device;
	
	__bounds *bounds, *bounds_device;
	char *axis, *axis_device;
	int *object_count, *object_count_device;

	// applies to leaf nodes
	int *objects, *objects_device;
};

struct __ray {
	__vector origin;
	__vector direction;
};

struct __intersection {
	bool intersects;
	double t;
	__ray ray;
	__vector normal;
	double temp0, temp1;
};

struct __object {
	int type; // geometry
	int material;
	int texture;
	__vector color;
	__vector emission;
	__vector e0;
	__vector e1;
	__vector e2;
	
	// for planes
	__vector normal;
	__vector point;

	// for spheres
	__vector center;
	double radius;
	
	// for triangles
	__vector p0, p1, p2;
	__vector n0, n1, n2;
	__vector n;
};

struct __buffers {
	unsigned char *char_host;
	unsigned char *char_device;
	double *doubles_host;
	double *doubles_device;
	__object *objects_host;
	__object *objects_device;
	float *rand_host;
	float *rand_device;
	unsigned int object_count;
};

struct __camera {
	double focal_length, aperture;
	double focal_distance, image_distance;
	double aperture_diameter;
};

bool initializePathTracer(__node& n, __buffers& b, int max_bounces, __dimensions dim, __bounds bounds, __textures& texture);
bool setupPathTracer     (__node& n, __buffers& b,                  __dimensions dim, __bounds bounds);
bool releasePathTracer   (__node& n, __buffers& b,                                                  __textures& texture);
bool runPathTracer       (__node& n, __buffers  b, int max_bounces, __dimensions dim,               __offsets offset, unsigned long long samples, __camera cmodel, __textures& texture);
bool grabFrame           (           __buffers  b, __dimensions dim);

#endif