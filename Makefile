all: main.cc lib/pathtracer.o lib/glhelper.o lib/vector.o lib/keyboard.o lib/timer.o lib/matrix.o lib/misc.o lib/complex.o lib/obj.o lib/object.o lib/buffer.o src/enum.h src/common.h src/vec.h
#	g++ main.cc /home/cuda-5.0/lib64/libcudart.so /home/cuda-5.0/lib64/libcurand.so lib/pathtracer.o lib/complex.o lib/glhelper.o lib/matrix.o lib/misc.o lib/vector.o lib/keyboard.o lib/timer.o -o bin/main -L/usr/lib `sdl-config --cflags --libs` -lGL -lGLEW
	g++ main.cc /home/keith/cuda5/lib64/libcudart.so /home/keith/cuda5/lib64/libcurand.so lib/obj.o lib/object.o lib/pathtracer.o lib/complex.o lib/glhelper.o lib/matrix.o lib/misc.o lib/vector.o lib/keyboard.o lib/timer.o lib/buffer.o -o bin/main -L/usr/lib `sdl-config --cflags --libs` -lGL -lGLEW -lSDL_image

lib/pathtracer.o: pathtracer.cu src/enum.h src/common.h src/vec_device.h
	nvcc pathtracer.cu -arch=sm_30 -c -o lib/pathtracer.o

lib/object.o: src/object.h src/object.cc src/vector.h src/matrix.h src/vec.h src/common.h src/enum.h
	g++ src/object.cc -c -o lib/object.o

lib/obj.o: src/obj.h src/obj.cc src/misc.h
	g++ src/obj.cc -c -o lib/obj.o

lib/perlin.o: src/perlin.h src/perlin.cc
	g++ src/perlin.cc -c -o lib/perlin.o

lib/buffer.o: src/buffer.h src/buffer.cc
	g++ src/buffer.cc -c -o lib/buffer.o

lib/glhelper.o: src/glhelper.h src/glhelper.cc src/misc.h
	g++ src/glhelper.cc -c -o lib/glhelper.o

lib/misc.o: src/misc.h src/misc.cc src/complex.h
	g++ src/misc.cc -c -o lib/misc.o

lib/complex.o: src/complex.h src/complex.cc
	g++ src/complex.cc -c -o lib/complex.o

lib/vector.o: src/vector.h src/vector.cc
	g++ src/vector.cc -c -o lib/vector.o

lib/timer.o: src/timer.h src/timer.cc
	g++ src/timer.cc -c -o lib/timer.o

lib/matrix.o: src/matrix.h src/matrix.cc src/vector.h
	g++ src/matrix.cc -c -o lib/matrix.o

lib/keyboard.o: src/keyboard.h src/keyboard.cc
	g++ src/keyboard.cc -c -o lib/keyboard.o

clean:
	@rm -f *~ src/*~ lib/* bin/*
