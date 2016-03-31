#ifndef ENUM_H
#define ENUM_H

enum event { PRIMITIVE_START, PRIMITIVE_END };
enum tree  { KD_EVEN, KD_MEDIAN, KD_SAH };
enum axis  { X, Y, Z, NOAXIS };

enum geometry  { NONE, PLANE, SPHERE, TRIANGLE };
enum materials { DIFFUSE, SPECULAR, REFRACTIVE };

#endif