#ifndef VEC_H
#define VEC_H

class __vector {
  public:
	double x, y, z;

	__vector();
	__vector(double x, double y, double z);

	double dot(__vector w);
	double operator*(__vector w);
	__vector cross(__vector w);
	
	__vector add(__vector w);
	__vector operator+(__vector w);
	__vector sub(__vector w);
	__vector operator-(__vector w);
	__vector mult(double s);
	__vector operator*(double s);
	
	__vector h(__vector w);
	double length();
	__vector unit();
};

#endif