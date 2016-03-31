#ifndef VEC_DEVICE_H
#define VEC_DEVICE_H

class __vector {
  public:
	double x, y, z;

	__host__ __device__ __vector();
	__host__ __device__ __vector(double x, double y, double z);

	__host__ __device__ double dot(__vector w);
	__host__ __device__ double operator*(__vector w);
	__host__ __device__ __vector cross(__vector w);
	
	__host__ __device__ __vector add(__vector w);
	__host__ __device__ __vector operator+(__vector w);
	__host__ __device__ __vector sub(__vector w);
	__host__ __device__ __vector operator-(__vector w);
	__host__ __device__ __vector mult(double s);
	__host__ __device__ __vector operator*(double s);
	
	__host__ __device__ __vector h(__vector w);
	__host__ __device__ double length();
	__host__ __device__ __vector unit();
};

__host__ __device__ __vector::__vector() : x(0.0), y(0.0), z(0.0) { }
__host__ __device__ __vector::__vector(double x, double y, double z) : x(x), y(y), z(z) { }

__host__ __device__ double __vector::dot(__vector w) { return this->x*w.x + this->y*w.y + this->z*w.z; }
__host__ __device__ double __vector::operator*(__vector w) { return this->x*w.x + this->y*w.y + this->z*w.z; }
__host__ __device__ __vector __vector::cross(__vector w) { return __vector(this->y * w.z - this->z * w.y, this->z * w.x - this->x * w.z, this->x * w.y - this->y * w.x); }

__host__ __device__ __vector __vector::add(__vector w) { return __vector(this->x + w.x, this->y + w.y, this->z + w.z); }
__host__ __device__ __vector __vector::operator+(__vector w) { return __vector(this->x + w.x, this->y + w.y, this->z + w.z); }
__host__ __device__ __vector __vector::sub(__vector w) { return __vector(this->x - w.x, this->y - w.y, this->z - w.z); }
__host__ __device__ __vector __vector::operator-(__vector w) { return __vector(this->x - w.x, this->y - w.y, this->z - w.z); }
__host__ __device__ __vector __vector::mult(double s) { return __vector(s * this->x, s * this->y, s * this->z); }
__host__ __device__ __vector __vector::operator*(double s) { return __vector(s * this->x, s * this->y, s * this->z); }

__host__ __device__ __vector __vector::h(__vector w) { return __vector(this->x * w.x, this->y * w.y, this->z * w.z); }
__host__ __device__ double __vector::length() { return sqrt(this->x*this->x + this->y*this->y + this->z*this->z); }
__host__ __device__ __vector __vector::unit() {
	double mag = sqrt(this->x*this->x + this->y*this->y + this->z*this->z);
	return __vector(this->x/mag, this->y/mag, this->z/mag);
}

#endif