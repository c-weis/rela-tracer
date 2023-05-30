SRC_PATH = src
FILES = main.cc test.cc colors.cc tracer.cc scene.cc materials.cc math.cc
SOURCES = $(FILES:%.cc=$(SRC_PATH)/%.cc)

build:
	g++ -O3 -o rela-tracer --std=c++17 -Isrc/include ${SOURCES}
debug:
	g++ -g3 -o debug --std=c++17 -Isrc/include ${SOURCES}
compiledebug:
	g++ -g3 --verbose -ftemplate-depth=20 -c --std=c++17 -Isrc/include ${SOURCES}
profile:
	g++ -O3 -pg -o profile --std=c++17 -Isrc/include ${SOURCES}
gputest:
	g++ -O3 -fopenacc -fcf-protection=none -no-pie -o test --std=c++17 -Isrc/include ${SOURCES}
