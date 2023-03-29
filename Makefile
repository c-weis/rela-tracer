SRC_PATH = src
FILES = test.cc tracer.cc scene.cc materials.cc math.cc
SOURCES = $(FILES:%.cc=$(SRC_PATH)/%.cc)
test:
	g++ -o test --std=c++17 -I. ${SOURCES}

clean:
	rm test test_image.ppm 
