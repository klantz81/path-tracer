#include "buffer.h"

cBuffer::cBuffer(const int WIDTH, const int HEIGHT) : width(WIDTH), height(HEIGHT) {
	buffer = 0;
	buffer = new unsigned char [WIDTH * HEIGHT * 4];
}

cBuffer::~cBuffer() {
	if (buffer) delete [] buffer;
}

void cBuffer::save(const int frame) {
	std::stringstream out;
	if      (frame < 10)
		out << "frame000" << frame << ".tga";
	else if (frame < 100)
		out << "frame00" << frame << ".tga";
	else if (frame < 1000)
		out << "frame0" << frame << ".tga";
	else if (frame < 10000)
		out << "frame" << frame << ".tga";
	std::string s = out.str();
	
	glReadPixels(0, 0, width, height, GL_BGRA, GL_UNSIGNED_BYTE, buffer);

	std::fstream of(s.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);

	char header[18] = { 0 };
	header[2]  = 2;
	header[12] = width & 0xff;
	header[13] = width >> 8;
	header[14] = height & 0xff;
	header[15] = height >> 8;
	header[16] = 32;

	of.write(header, 18);
	of.write((char *)buffer, width * height * 4);
}