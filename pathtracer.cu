#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <curand.h>
#include <vector>
#include <iostream>

#include "src/vec_device.h"
#include "src/stack_device.h"
#include "src/common.h"
#include "src/enum.h"


__device__ __intersection_t rayIntersectsAABB(const __node& n, unsigned short index, __ray& r) {

	__intersection_t isect;
	isect.intersects = false;
	isect.tmin = 0;
	isect.tmax = 1e8;

	char side = 0;

	if (r.origin.x >= n.bounds_device[index].minx && r.origin.x <= n.bounds_device[index].maxx &&
	    r.origin.y >= n.bounds_device[index].miny && r.origin.y <= n.bounds_device[index].maxy &&
	    r.origin.z >= n.bounds_device[index].minz && r.origin.z <= n.bounds_device[index].maxz) {
		isect.intersects = true;
		isect.tmin = 0;
	}
	
	if (!isect.intersects && r.origin.x < n.bounds_device[index].minx && r.direction.x > 0) {
		double t = (n.bounds_device[index].minx - r.origin.x) / r.direction.x;
		double y = r.origin.y + t * r.direction.y;
		double z = r.origin.z + t * r.direction.z;
		if (y >= n.bounds_device[index].miny && y <= n.bounds_device[index].maxy && z >= n.bounds_device[index].minz && z <= n.bounds_device[index].maxz) {
			isect.intersects = true;
			isect.tmin = t;
			side = 1;
		}
	}
	if (!isect.intersects && r.origin.x > n.bounds_device[index].maxx && r.direction.x < 0) {
		double t = (n.bounds_device[index].maxx - r.origin.x) / r.direction.x;
		double y = r.origin.y + t * r.direction.y;
		double z = r.origin.z + t * r.direction.z;
		if (y >= n.bounds_device[index].miny && y <= n.bounds_device[index].maxy && z >= n.bounds_device[index].minz && z <= n.bounds_device[index].maxz) {
			isect.intersects = true;
			isect.tmin = t;
			side = 2;
		}
	}
	if (!isect.intersects && r.origin.y < n.bounds_device[index].miny && r.direction.y > 0) {
		double t = (n.bounds_device[index].miny - r.origin.y) / r.direction.y;
		double x = r.origin.x + t * r.direction.x;
		double z = r.origin.z + t * r.direction.z;
		if (x >= n.bounds_device[index].minx && x <= n.bounds_device[index].maxx && z >= n.bounds_device[index].minz && z <= n.bounds_device[index].maxz) {
			isect.intersects = true;
			isect.tmin = t;
			side = 3;
		}
	}
	if (!isect.intersects && r.origin.y > n.bounds_device[index].maxy && r.direction.y < 0) {
		double t = (n.bounds_device[index].maxy - r.origin.y) / r.direction.y;
		double x = r.origin.x + t * r.direction.x;
		double z = r.origin.z + t * r.direction.z;
		if (x >= n.bounds_device[index].minx && x <= n.bounds_device[index].maxx && z >= n.bounds_device[index].minz && z <= n.bounds_device[index].maxz) {
			isect.intersects = true;
			isect.tmin = t;
			side = 4;
		}
	}
	if (!isect.intersects && r.origin.z < n.bounds_device[index].minz && r.direction.z > 0) {
		double t = (n.bounds_device[index].minz - r.origin.z) / r.direction.z;
		double x = r.origin.x + t * r.direction.x;
		double y = r.origin.y + t * r.direction.y;
		if (x >= n.bounds_device[index].minx && x <= n.bounds_device[index].maxx && y >= n.bounds_device[index].miny && y <= n.bounds_device[index].maxy) {
			isect.intersects = true;
			isect.tmin = t;
			side = 5;
		}
	}
	if (!isect.intersects && r.origin.z > n.bounds_device[index].maxz && r.direction.z < 0) {
		double t = (n.bounds_device[index].maxz - r.origin.z) / r.direction.z;
		double x = r.origin.x + t * r.direction.x;
		double y = r.origin.y + t * r.direction.y;
		if (x >= n.bounds_device[index].minx && x <= n.bounds_device[index].maxx && y >= n.bounds_device[index].miny && y <= n.bounds_device[index].maxy) {
			isect.intersects = true;
			isect.tmin = t;
			side = 6;
		}
	}

	if (isect.intersects) {
		double epsilon = 0.0000000000001;
		if (side != 1) {
			__vector n0(1.0,0.0,0.0);
			__vector p0(n.bounds_device[index].minx, n.bounds_device[index].miny, n.bounds_device[index].minz);
			double den = r.direction * n0;
			if (fabs(den) > epsilon) {
				double t = ((p0 - r.origin) * n0) / den;
				if (t > isect.tmin) {
					__vector _p0 = r.origin + r.direction * t;
					if (_p0.y >= n.bounds_device[index].miny && _p0.y <= n.bounds_device[index].maxy && _p0.z >= n.bounds_device[index].minz && _p0.z <= n.bounds_device[index].maxz) {
						isect.tmax = t;
					}
				}
			}
		}
		if (side != 2) {
			__vector n0(1.0,0.0,0.0);
			__vector p0(n.bounds_device[index].maxx, n.bounds_device[index].miny, n.bounds_device[index].minz);
			double den = r.direction * n0;
			if (fabs(den) > epsilon) {
				double t = ((p0 - r.origin) * n0) / den;
				if (t > isect.tmin) {
					__vector _p0 = r.origin + r.direction * t;
					if (_p0.y >= n.bounds_device[index].miny && _p0.y <= n.bounds_device[index].maxy && _p0.z >= n.bounds_device[index].minz && _p0.z <= n.bounds_device[index].maxz) {
						isect.tmax = t;
					}
				}
			}
		}
		if (side != 3) {
			__vector n0(0.0,1.0,0.0);
			__vector p0(n.bounds_device[index].minx, n.bounds_device[index].miny, n.bounds_device[index].minz);
			double den = r.direction * n0;
			if (fabs(den) > epsilon) {
				double t = ((p0 - r.origin) * n0) / den;
				if (t > isect.tmin) {
					__vector _p0 = r.origin + r.direction * t;
					if (_p0.x >= n.bounds_device[index].minx && _p0.x <= n.bounds_device[index].maxx && _p0.z >= n.bounds_device[index].minz && _p0.z <= n.bounds_device[index].maxz) {
						isect.tmax = t;
					}
				}
			}
		}
		if (side != 4) {
			__vector n0(0.0,1.0,0.0);
			__vector p0(n.bounds_device[index].minx, n.bounds_device[index].maxy, n.bounds_device[index].minz);
			double den = r.direction * n0;
			if (fabs(den) > epsilon) {
				double t = ((p0 - r.origin) * n0) / den;
				if (t > isect.tmin) {
					__vector _p0 = r.origin + r.direction * t;
					if (_p0.x >= n.bounds_device[index].minx && _p0.x <= n.bounds_device[index].maxx && _p0.z >= n.bounds_device[index].minz && _p0.z <= n.bounds_device[index].maxz) {
						isect.tmax = t;
					}
				}
			}
		}
		if (side != 5) {
			__vector n0(0.0,0.0,1.0);
			__vector p0(n.bounds_device[index].minx, n.bounds_device[index].miny, n.bounds_device[index].minz);
			double den = r.direction * n0;
			if (fabs(den) > epsilon) {
				double t = ((p0 - r.origin) * n0) / den;
				if (t > isect.tmin) {
					__vector _p0 = r.origin + r.direction * t;
					if (_p0.y >= n.bounds_device[index].miny && _p0.y <= n.bounds_device[index].maxy && _p0.x >= n.bounds_device[index].minx && _p0.x <= n.bounds_device[index].maxx) {
						isect.tmax = t;
					}
				}
			}
		}
		if (side != 6) {
			__vector n0(0.0,0.0,1.0);
			__vector p0(n.bounds_device[index].minx, n.bounds_device[index].miny, n.bounds_device[index].maxz);
			double den = r.direction * n0;
			if (fabs(den) > epsilon) {
				double t = ((p0 - r.origin) * n0) / den;
				if (t > isect.tmin) {
					__vector _p0 = r.origin + r.direction * t;
					if (_p0.y >= n.bounds_device[index].miny && _p0.y <= n.bounds_device[index].maxy && _p0.x >= n.bounds_device[index].minx && _p0.x <= n.bounds_device[index].maxx) {
						isect.tmax = t;
					}
				}
			}
		}
	}

	return isect;
}



__device__ __vector sampleRay(__node n, __ray ray, __object objects[], int object_count, int depth, float* rand_device, int rand_index, int rand_size, int max_bounces, __textures texture) {
	double epsilon = 0.0000000000001;

	__intersection intersect; // best intersection
	__intersection isect; // test intersection
	
	int which = -1;

	max_bounces = MIN(20, MAX(max_bounces, 1));
	__vector __a[20], __b[20], sample(0,0,0);

	for (int l = 0; l < max_bounces; l++) {
		__a[l].x = __a[l].y = __a[l].z = 0;
		__b[l].x = __b[l].y = __b[l].z = 0;
	}

	for (int l = 0; l < max_bounces; l++) {

		ray.direction = ray.direction.unit();
		intersect.intersects = false;

		__stack stack;
		__intersection_t root_intersection = rayIntersectsAABB(n, 0, ray);

		if (root_intersection.intersects) {
			stack.push(0, root_intersection.tmin, root_intersection.tmax);
			while (!stack.empty() && !intersect.intersects) {
				__stack_element se = stack.pop();
				while (n.leaf_id_device[se.id] < 0) {
					char axis = n.axis_device[se.id];
					double tsplit;
					int first, second;
					if (axis == Z) {
						tsplit = (n.bounds_device[n.child0_device[se.id]].maxz - ray.origin.z) / ray.direction.z; // what about rays parallel to the split plane?
						if (n.bounds_device[n.child1_device[se.id]].minz - ray.origin.z >= 0.0) {
							first = n.child0_device[se.id];
							second = n.child1_device[se.id];
						} else {
							first = n.child1_device[se.id];
							second = n.child0_device[se.id];
						}
					} else if (axis == Y) {
						tsplit = (n.bounds_device[n.child0_device[se.id]].maxy - ray.origin.y) / ray.direction.y; // what about rays parallel to the split plane?
						if (n.bounds_device[n.child1_device[se.id]].miny - ray.origin.y >= 0.0) {
							first = n.child0_device[se.id];
							second = n.child1_device[se.id];
						} else {
							first = n.child1_device[se.id];
							second = n.child0_device[se.id];
						}
					} else {
						tsplit = (n.bounds_device[n.child0_device[se.id]].maxx - ray.origin.x) / ray.direction.x; // what about rays parallel to the split plane?
						if (n.bounds_device[n.child1_device[se.id]].minx - ray.origin.x >= 0.0) {
							first = n.child0_device[se.id];
							second = n.child1_device[se.id];
						} else {
							first = n.child1_device[se.id];
							second = n.child0_device[se.id];
						}
					}

					if (tsplit >= se.tmax || tsplit < 0) {
						se.id = first;
					} else if (tsplit <= se.tmin) {
						se.id = second;
					} else {
						stack.push(second, tsplit, se.tmax);
						se.id = first;
						se.tmax = tsplit;
					}
				}

				// check intersections for se.id
				for (int j = 0; j < n.object_count_device[se.id]; j++) {
					int k = n.objects_device[n.leaf_id_device[se.id] * n.max_leaf_objects + j];

					isect.intersects = false;

					if (objects[k].type == SPHERE) {

						double a = ray.direction * ray.direction;
						double b = (ray.direction * ray.origin - ray.direction * objects[k].center) * 2.0;
						double c = ray.origin * ray.origin + objects[k].center * objects[k].center - ray.origin * objects[k].center * 2.0 - objects[k].radius * objects[k].radius;
						double det = b * b - 4 * a * c;
						if (det < epsilon) continue;

						double t0 = (-b + sqrt(det))/(2 * a);
						double t1 = (-b - sqrt(det))/(2 * a);
						if (t0 < epsilon && t1 < epsilon) continue;

						isect.intersects = true;
						isect.t = t0 < epsilon ? t1 : (t1 < epsilon ? t0 : (t0 < t1 ? t0 : t1));

						isect.ray.origin = ray.origin + ray.direction * isect.t;

						isect.normal = (isect.ray.origin - objects[k].center).unit();

					} else if (objects[k].type == TRIANGLE) {

						double den = ray.direction * objects[k].n;
						if (fabs(den) < epsilon) continue;

						double num = (objects[k].p0 - ray.origin) * objects[k].n;
						double num_den = num/den;
						if (num_den < epsilon) continue;

						__vector v0 = objects[k].p1 - objects[k].p0;
						__vector v1 = objects[k].p2 - objects[k].p0;
						__vector p = (ray.origin + ray.direction * num_den) - objects[k].p0;
							
						double pv0 = p*v0;
						double pv1 = p*v1;
						double v0v0 = v0*v0;
						double v0v1 = v0*v1;
						double v1v1 = v1*v1;
						den = v0v0*v1v1 - v0v1*v0v1;
						double s = (pv0*v1v1 - pv1*v0v1)/den;
						double t = (pv1*v0v0 - pv0*v0v1)/den;
						if (s >= 0 && t >= 0 && s+t<1.0) {

							isect.intersects = true;
							isect.t = num_den;
							isect.temp0 = s;
							isect.temp1 = t;

							isect.ray.origin = ray.origin + ray.direction * isect.t;

							isect.normal = (objects[k].n0 + (objects[k].n1 - objects[k].n0)*s + (objects[k].n2 - objects[k].n0)*t).unit();
						}

					}

					if (isect.intersects && isect.t < se.tmax) {
						if (!intersect.intersects || isect.t < intersect.t) {
							intersect = isect;
							which = k;
						}
					}
				}
			}
		}


		
		
		for (int j = 0; j < object_count; j++) {

			int k = j;

			isect.intersects = false;

			if (objects[k].type == PLANE) {

				double den = ray.direction * objects[k].normal;
				if (fabs(den) < epsilon) continue;

				double num = (objects[k].point - ray.origin) * objects[k].normal;
				double num_den = num/den;
				if (num_den < epsilon) continue;

				isect.intersects = true;
				isect.t = num_den;

				isect.ray.origin = ray.origin + ray.direction * isect.t;

				isect.normal = objects[k].normal.unit();
				
			} else break;

			if (isect.intersects) {
				if (!intersect.intersects || isect.t < intersect.t) {
					intersect = isect;
					which = k;
				}
			}
		}


		

		if (intersect.intersects) {

			__a[l] = objects[which].emission;
			__b[l] = objects[which].color;

			// texture mapping
			if (objects[which].type == PLANE && objects[which].texture > -1) {
				__vector v0 = objects[which].e0;
				__vector v1 = objects[which].e1;
				__vector p  = (intersect.ray.origin - objects[which].point);

				double pv0 = p*v0;
				double pv1 = p*v1;
				double v0v0 = v0*v0;
				double v0v1 = v0*v1;
				double v1v1 = v1*v1;
				double den = v0v0*v1v1 - v0v1*v0v1;
				double s = (pv0*v1v1 - pv1*v0v1)/den;
				double t = (pv1*v0v0 - pv0*v0v1)/den;

				s = fabs(s);
				t = fabs(t);
				s = s - int(s);
				t = t - int(t);
				int twidth  = texture.width[objects[which].texture];
				int theight = texture.height[objects[which].texture];
				int tb      = texture.bpp[objects[which].texture] / 8;
				s *= twidth - 1;
				t *= theight - 1;
				unsigned int  _t =  t,      _s =  s,
				             __t = _t + 1, __s = _s + 1;
				 _s =  _s > (twidth)  - 2 ? (twidth  - 2) :  _s;
				__s = __s > (twidth)  - 1 ? (twidth  - 1) : __s;
				 _t =  _t > (theight) - 2 ? (theight - 2) :  _t;
				__t = __t > (theight) - 1 ? (theight - 1) : __t;

				__vector _a = __vector(texture.texture_device[objects[which].texture][ _t*twidth*tb +  _s*tb + 0]/256.0,
						       texture.texture_device[objects[which].texture][ _t*twidth*tb +  _s*tb + 1]/256.0,
						       texture.texture_device[objects[which].texture][ _t*twidth*tb +  _s*tb + 2]/256.0) * (__t -  t) * (__s -  s);
				__vector _b = __vector(texture.texture_device[objects[which].texture][__t*twidth*tb +  _s*tb + 0]/256.0,
						       texture.texture_device[objects[which].texture][__t*twidth*tb +  _s*tb + 1]/256.0, 
						       texture.texture_device[objects[which].texture][__t*twidth*tb +  _s*tb + 2]/256.0) * (  t - _t) * (__s -  s);
				__vector _c = __vector(texture.texture_device[objects[which].texture][ _t*twidth*tb + __s*tb + 0]/256.0,
						       texture.texture_device[objects[which].texture][ _t*twidth*tb + __s*tb + 1]/256.0,
						       texture.texture_device[objects[which].texture][ _t*twidth*tb + __s*tb + 2]/256.0) * (__t -  t) * (  s - _s);
				__vector _d = __vector(texture.texture_device[objects[which].texture][__t*twidth*tb + __s*tb + 0]/256.0,
						       texture.texture_device[objects[which].texture][__t*twidth*tb + __s*tb + 1]/256.0,
						       texture.texture_device[objects[which].texture][__t*twidth*tb + __s*tb + 2]/256.0) * (  t - _t) * (  s - _s);

				__a[l] = __a[l].h(_a + _b + _c + _d);
				__b[l] = __b[l].h(_a + _b + _c + _d);
			} else if (objects[which].type == SPHERE && objects[which].texture > -1) {
				__vector v0 = objects[which].e0.unit();
				__vector v1 = objects[which].e1.unit();
				__vector p  = (intersect.ray.origin - objects[which].center).unit();

				// rework this
				double phi = acos(v0*p);
				double t = phi/3.14159265359;
				double theta = acos(v1*p/sin(phi))/(2*3.14159265359);
				double s = v0.cross(v1)*p > 0.0 ? theta : (1.0-theta);

				int twidth  = texture.width[objects[which].texture];
				int theight = texture.height[objects[which].texture];
				int tb      = texture.bpp[objects[which].texture] / 8;
				s *= twidth - 1;
				t *= theight - 1;
				unsigned int  _t =  t,      _s =  s,
				             __t = _t + 1, __s = _s + 1;
				 _s =  _s > (twidth)  - 2 ? (twidth  - 2) :  _s;
				__s = __s > (twidth)  - 1 ? (twidth  - 1) : __s;
				 _t =  _t > (theight) - 2 ? (theight - 2) :  _t;
				__t = __t > (theight) - 1 ? (theight - 1) : __t;

				__vector _a = __vector(texture.texture_device[objects[which].texture][ _t*twidth*tb +  _s*tb + 0]/256.0,
						       texture.texture_device[objects[which].texture][ _t*twidth*tb +  _s*tb + 1]/256.0,
						       texture.texture_device[objects[which].texture][ _t*twidth*tb +  _s*tb + 2]/256.0) * (__t -  t) * (__s -  s);
				__vector _b = __vector(texture.texture_device[objects[which].texture][__t*twidth*tb +  _s*tb + 0]/256.0,
						       texture.texture_device[objects[which].texture][__t*twidth*tb +  _s*tb + 1]/256.0, 
						       texture.texture_device[objects[which].texture][__t*twidth*tb +  _s*tb + 2]/256.0) * (  t - _t) * (__s -  s);
				__vector _c = __vector(texture.texture_device[objects[which].texture][ _t*twidth*tb + __s*tb + 0]/256.0,
						       texture.texture_device[objects[which].texture][ _t*twidth*tb + __s*tb + 1]/256.0,
						       texture.texture_device[objects[which].texture][ _t*twidth*tb + __s*tb + 2]/256.0) * (__t -  t) * (  s - _s);
				__vector _d = __vector(texture.texture_device[objects[which].texture][__t*twidth*tb + __s*tb + 0]/256.0,
						       texture.texture_device[objects[which].texture][__t*twidth*tb + __s*tb + 1]/256.0,
						       texture.texture_device[objects[which].texture][__t*twidth*tb + __s*tb + 2]/256.0) * (  t - _t) * (  s - _s);

				__a[l] = __a[l].h(_a + _b + _c + _d);
				__b[l] = __b[l].h(_a + _b + _c + _d);
			} else if (objects[which].type == TRIANGLE && objects[which].texture > -1) {
				double s = intersect.temp0;
				double t = intersect.temp1;

				__vector v0 = objects[which].e1 - objects[which].e0;
				__vector v1 = objects[which].e2 - objects[which].e0;

				__vector texcoord = objects[which].e0 + v0 * s + v1 * t;
				s = texcoord.x;
				t = 1.0 - texcoord.y;
				
				int twidth  = texture.width[objects[which].texture];
				int theight = texture.height[objects[which].texture];
				int tb      = texture.bpp[objects[which].texture] / 8;
				s *= twidth - 1;
				t *= theight - 1;
				unsigned int  _t =  t,      _s =  s,
				             __t = _t + 1, __s = _s + 1;
				 _s =  _s > (twidth)  - 2 ? (twidth  - 2) :  _s;
				__s = __s > (twidth)  - 1 ? (twidth  - 1) : __s;
				 _t =  _t > (theight) - 2 ? (theight - 2) :  _t;
				__t = __t > (theight) - 1 ? (theight - 1) : __t;

				__vector _a = __vector(texture.texture_device[objects[which].texture][ _t*twidth*tb +  _s*tb + 0]/256.0,
						       texture.texture_device[objects[which].texture][ _t*twidth*tb +  _s*tb + 1]/256.0,
						       texture.texture_device[objects[which].texture][ _t*twidth*tb +  _s*tb + 2]/256.0) * (__t -  t) * (__s -  s);
				__vector _b = __vector(texture.texture_device[objects[which].texture][__t*twidth*tb +  _s*tb + 0]/256.0,
						       texture.texture_device[objects[which].texture][__t*twidth*tb +  _s*tb + 1]/256.0, 
						       texture.texture_device[objects[which].texture][__t*twidth*tb +  _s*tb + 2]/256.0) * (  t - _t) * (__s -  s);
				__vector _c = __vector(texture.texture_device[objects[which].texture][ _t*twidth*tb + __s*tb + 0]/256.0,
						       texture.texture_device[objects[which].texture][ _t*twidth*tb + __s*tb + 1]/256.0,
						       texture.texture_device[objects[which].texture][ _t*twidth*tb + __s*tb + 2]/256.0) * (__t -  t) * (  s - _s);
				__vector _d = __vector(texture.texture_device[objects[which].texture][__t*twidth*tb + __s*tb + 0]/256.0,
						       texture.texture_device[objects[which].texture][__t*twidth*tb + __s*tb + 1]/256.0,
						       texture.texture_device[objects[which].texture][__t*twidth*tb + __s*tb + 2]/256.0) * (  t - _t) * (  s - _s);

				__a[l] = __a[l].h(_a + _b + _c + _d);
				__b[l] = __b[l].h(_a + _b + _c + _d);
			}


			if (objects[which].material == DIFFUSE) {

				__vector w = intersect.normal;

				// cosine weighted sampling
				double u1 = rand_device[rand_size*(l+1) + rand_index + 0];
				double u2 = rand_device[rand_size*(l+1) + rand_index + 1];
				double r1 = 2 * M_PI * u1;
				double r2 = sqrt(1 - u2);
				double r3 = sqrt(u2);

				__vector u(0,0,0);
				if      (fabs(w.x) < fabs(w.y) + epsilon && fabs(w.x) < fabs(w.z) + epsilon) u.x = 1;
				else if (fabs(w.y) < fabs(w.x) + epsilon && fabs(w.y) < fabs(w.z) + epsilon) u.y = 1;
				else u.z = 1;

				u = u.cross(w).unit();
				__vector v = w.cross(u).unit();
					 u = v.cross(w).unit();
				__vector d = ( u * (cos(r1) * r2) + v * (sin(r1) * r2) + w * r3 ) .unit();

				intersect.ray.direction = d;
				ray = intersect.ray;

			} else if (objects[which].material == SPECULAR) {

				intersect.ray.direction = (ray.direction - intersect.normal * (ray.direction * intersect.normal * 2.0)).unit();
				ray = intersect.ray;
				
			} else if (objects[which].material == REFRACTIVE) {

				bool into = ray.direction * intersect.normal < 0; // entering the medium

				double n1 = into ? 1.0 : 1.5;
				double n2 = into ? 1.5 : 1.0;
				double n1n2 = n1/n2;
				__vector n  = into ? intersect.normal : (intersect.normal * -1);
				__vector r = ray.direction;
				
				double n1n22 = n1n2 * n1n2;
				double rn   = r * n;
				double rn2  = rn * rn;
				
				double a = 1 - n1n22 * (1 - rn2);
				if (a >= 0) {
					// fresnel
					double a_sqrt = sqrt(a);
					double Rs = ((n1 * rn * -1) - (n2 * a_sqrt)) / ((n1 * rn * -1) + (n2 * a_sqrt));
					Rs *= Rs;
					double Rp = ((n1 * a_sqrt) - (n2 * rn * -1)) / ((n1 * a_sqrt) + (n2 * rn * -1));
					Rp *= Rp;
					double R = (Rs + Rp)/2.0;
					//double T = 1 - R;
					if (rand_device[rand_size*(l+1) + rand_index + 2] < R) {	// reflect
						intersect.ray.direction = (ray.direction - intersect.normal * (ray.direction * intersect.normal * 2.0)).unit();
						ray = intersect.ray;
					} else {							// transmit
						//ray.origin = intersect.ray.origin;
						//ray.direction = r * n1n2 - n * (n1n2 * rn + a_sqrt);
						intersect.ray.direction = r * n1n2 - n * (n1n2 * rn + a_sqrt);
						ray = intersect.ray;
					}
				} else {
					  intersect.ray.direction = (ray.direction - intersect.normal * (ray.direction * intersect.normal * 2.0)).unit();
					  ray = intersect.ray; // total internal reflection
				}
			}
			
		} else break;
	}

	sample = __a[max_bounces - 1];
	for (int l = max_bounces-2; l >= 0; l--) sample = __a[l] + __b[l].h(sample);

	return sample;
}

__global__ void kernel(__node n, __buffers b, int samples, int max_bounces, __dimensions dim, __offsets offset, __camera cmodel, __textures texture) {
  	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int y = threadIdx.y + blockIdx.y * blockDim.y;

	int rand_index = y * 8 * 16 * 3 + x * 3;

	int index = ((y + offset.y) * dim.width + (dim.width - x - offset.x - 1)) * 3;
	if (index < 0 || index >= dim.height * dim.width * 3) return;

	double u1 = b.rand_device[rand_index + 0];
	double u2 = b.rand_device[rand_index + 1];
	double r1 = 2 * M_PI * u1;
	double r2 = sqrt(1 - u2);
	__vector o = __vector( (x + offset.x) - dim.width  / 2,
			       (y + offset.y) - dim.height / 2,
			                                    0) + __vector(cos(r1)*r2, sin(r1)*r2, 0) * 0.5;
	
	o.x =   o.x/dim.width  * 36.0/1000.0;			// 36mm sensor
	o.y =   o.y/dim.height * 36.0/1000.0 * dim.height/dim.width;
	o.z =   cmodel.image_distance;

	__ray ray;
	ray.origin = o;
	ray.direction = o.unit() * -1;

	__vector p = ray.direction * (cmodel.focal_distance / fabs(ray.direction.z));

	u1 = b.rand_device[rand_index+1];
	u2 = b.rand_device[rand_index+2];
	r1 = 2 * M_PI * u1;
	r2 = u2;
	
	ray.origin = __vector(cos(r1)*r2, sin(r1)*r2, 0) * cmodel.aperture_diameter * 0.5;
	ray.direction = (p - ray.origin).unit();

	// sample the ray
	//__vector _sample = sampleRay(b, ray, objects, object_count, 0, rand_device, rand_index, dim.width*dim.height*3, max_bounces, texture);
	__vector _sample = sampleRay(n, ray, b.objects_device, b.object_count, 0, b.rand_device, rand_index, 16 * 8 * 16 * 8 * 3, max_bounces, texture);

	// add the sample to the  accumulation
	b.doubles_device[index+0] = (b.doubles_device[index + 0] * (samples - 1.0) + _sample.x) / samples;
	b.doubles_device[index+1] = (b.doubles_device[index + 1] * (samples - 1.0) + _sample.y) / samples;
	b.doubles_device[index+2] = (b.doubles_device[index + 2] * (samples - 1.0) + _sample.z) / samples;

	// save the current frame
	b.char_device[index+0] = CLAMP((int)(b.doubles_device[index + 0] * 255.0), 0, 255);
	b.char_device[index+1] = CLAMP((int)(b.doubles_device[index + 1] * 255.0), 0, 255);
	b.char_device[index+2] = CLAMP((int)(b.doubles_device[index + 2] * 255.0), 0, 255);
}

bool initializePathTracer(__node& n, __buffers& b, int max_bounces, __dimensions dim, __bounds bounds, __textures& texture) {

	int num_bytes;

// allocate tree -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	printf("%s\n", cudaGetErrorString(  cudaMalloc((void**)&n.id_device, (sizeof(int) * n.size))   ));
	printf("%s\n", cudaGetErrorString(  cudaMalloc((void**)&n.leaf_id_device, (sizeof(int) * n.size))   ));
	printf("%s\n", cudaGetErrorString(  cudaMalloc((void**)&n.parent_device, (sizeof(int) * n.size))   ));
	printf("%s\n", cudaGetErrorString(  cudaMalloc((void**)&n.child0_device, (sizeof(int) * n.size))   ));
	printf("%s\n", cudaGetErrorString(  cudaMalloc((void**)&n.child1_device, (sizeof(int) * n.size))   ));
	printf("%s\n", cudaGetErrorString(  cudaMalloc((void**)&n.bounds_device, (sizeof(__bounds) * n.size))   ));
	printf("%s\n", cudaGetErrorString(  cudaMalloc((void**)&n.axis_device, (sizeof(char) * n.size))   ));
	printf("%s\n", cudaGetErrorString(  cudaMalloc((void**)&n.object_count_device, (sizeof(int) * n.size))   ));
	printf("%s\n", cudaGetErrorString(  cudaMalloc((void**)&n.objects_device, (sizeof(int) * n.leaf_size * n.max_leaf_objects))   ));
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// for storing the image result
	num_bytes = sizeof(unsigned char) * dim.width * dim.height * 3;
	b.char_host = (unsigned char*)malloc(num_bytes);
	cudaMalloc((void**)&(b.char_device), num_bytes);

	// for accumulating samples
	num_bytes = sizeof(double) * dim.width * dim.height * 3;
	b.doubles_host = (double*)malloc(num_bytes);
	cudaMalloc((void**)&(b.doubles_device), num_bytes);

	// for storing the objects on the device
	num_bytes = sizeof(__object) * b.object_count;
	cudaMalloc((void**)&b.objects_device, num_bytes);

	// for generating random uniforms
	int m = 1 + max_bounces;
//	num_bytes = sizeof(float) * dim.width * dim.height * 3 * m;
	num_bytes = sizeof(float) * 16 * 8 * 16 * 8 * 3 * m;
	cudaMalloc((void**)&b.rand_device, num_bytes);
	
	// for storing the textures
	for (int i = 0; i < texture.count; i++) {
		num_bytes = sizeof(unsigned char) * texture.width[i] * texture.height[i] * texture.bpp[i] / 8;
		cudaMalloc((void**)&(texture.texture_device[i]), num_bytes);
		cudaMemcpy(texture.texture_device[i], texture.texture[i], num_bytes, cudaMemcpyHostToDevice);
	}

	return true;
}

bool setupPathTracer(__node& n, __buffers& b, __dimensions dim, __bounds bounds) {

	int num_bytes;

// copy tree -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	printf("%s\n", cudaGetErrorString(  cudaMemcpy(n.id_device, n.id, sizeof(int) * n.size, cudaMemcpyHostToDevice)  ));
	printf("%s\n", cudaGetErrorString(  cudaMemcpy(n.leaf_id_device, n.leaf_id, sizeof(int) * n.size, cudaMemcpyHostToDevice)  ));
	printf("%s\n", cudaGetErrorString(  cudaMemcpy(n.parent_device, n.parent, sizeof(int) * n.size, cudaMemcpyHostToDevice)  ));
	printf("%s\n", cudaGetErrorString(  cudaMemcpy(n.child0_device, n.child0, sizeof(int) * n.size, cudaMemcpyHostToDevice)  ));
	printf("%s\n", cudaGetErrorString(  cudaMemcpy(n.child1_device, n.child1, sizeof(int) * n.size, cudaMemcpyHostToDevice)  ));
	printf("%s\n", cudaGetErrorString(  cudaMemcpy(n.bounds_device, n.bounds, sizeof(__bounds) * n.size, cudaMemcpyHostToDevice)  ));
	printf("%s\n", cudaGetErrorString(  cudaMemcpy(n.axis_device, n.axis, sizeof(char) * n.size, cudaMemcpyHostToDevice)  ));
	printf("%s\n", cudaGetErrorString(  cudaMemcpy(n.object_count_device, n.object_count, sizeof(int) * n.size, cudaMemcpyHostToDevice)  ));
	printf("%s\n", cudaGetErrorString(  cudaMemcpy(n.objects_device, n.objects, sizeof(int) * n.leaf_size * n.max_leaf_objects, cudaMemcpyHostToDevice)  ));
// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// clear the image result
	num_bytes = sizeof(unsigned char) * dim.width * dim.height * 3;
	memset(b.char_host, 0, num_bytes);
	cudaMemset(b.char_device, 0, num_bytes);
	// zero out the sample accumulation
	num_bytes = sizeof(double) * dim.width * dim.height * 3;
	memset(b.doubles_host, 0, num_bytes);
	cudaMemset(b.doubles_device, 0, num_bytes);
	// store the objects on the device
	num_bytes = sizeof(__object) * b.object_count;
	cudaMemcpy(b.objects_device, b.objects_host, num_bytes, cudaMemcpyHostToDevice);

	return true;
}

bool releasePathTracer(__node& n, __buffers& b, __textures& texture) {

///////////////////
	cudaFree(n.id_device);
	cudaFree(n.leaf_id_device);
	cudaFree(n.parent_device);
	cudaFree(n.child0_device);
	cudaFree(n.child1_device);
	cudaFree(n.bounds_device);
	cudaFree(n.axis_device);
	cudaFree(n.object_count_device);
	cudaFree(n.objects_device);
///////////////////

	free(b.char_host);
	cudaFree(b.char_device);
	free(b.doubles_host);
	cudaFree(b.doubles_device);
	cudaFree(b.objects_device);
	cudaFree(b.rand_device);

	for (int i = 0; i < texture.count; i++) cudaFree(texture.texture_device[i]);

	return true;
}

bool runPathTracer(__node& n, __buffers b, int max_bounces, __dimensions dim, __offsets offset, unsigned long long samples, __camera cmodel, __textures& texture) {

	int dimx = 16;
	int dimy = 16;
	dim3 dimGrid(8,8);
//	dim3 dimGrid(dim.width/dimx, dim.height/dimy);
	dim3 dimBlock(dimx, dimy);

	// generate random uniforms
	if (offset.x == 0 && offset.y == 0) {
		int m = 1 + max_bounces;//4 + max_bounces * 4; // 1 plus 10 depth
		curandGenerator_t gen;
		curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_MTGP32);
		curandSetPseudoRandomGeneratorSeed(gen, samples); //samples=seed
//		curandGenerateUniform(gen, rand_device, dim.width * dim.height * 3 * m);
		curandGenerateUniform(gen, b.rand_device, 16 * 8 * 16 * 8 * 3 * m);
		curandDestroyGenerator(gen);
	}
	
	//cudaDeviceSetLimit (cudaLimitStackSize, );
	//cudaThreadSetLimit (cudaLimitStackSize, 8192*16)

	kernel<<<dimGrid, dimBlock>>>(n, b, samples, max_bounces, dim, offset, cmodel, texture);

	return true;
}

bool grabFrame(__buffers b, __dimensions dim) {
	int num_bytes;

	num_bytes = sizeof(unsigned char) * dim.width * dim.height * 3;
	printf("%s\n", cudaGetErrorString( cudaMemcpy(b.char_host, b.char_device, num_bytes, cudaMemcpyDeviceToHost) ));

//	num_bytes = sizeof(double) * dim.width * dim.height * 3;
//	printf("%s\n", cudaGetErrorString( cudaMemcpy(b.doubles_host, b.doubles_device, num_bytes, cudaMemcpyDeviceToHost) ));

	return true;
}