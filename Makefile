SRC_PATH = src
FILES = colors.cc tracer.cc scene.cc materials.cc math.cc
SOURCES = $(FILES:%.cc=$(SRC_PATH)/%.cc)

debug:
	g++ -g3 -o debug --std=c++17 -Isrc/include $(SRC_PATH)/test.cc ${SOURCES}
compiledebug:
	g++ -g3 --verbose -ftemplate-depth=20 -c --std=c++17 -Isrc/include $(SRC_PATH)/test.cc ${SOURCES}
profile:
	g++ -O3 -pg -o profile --std=c++17 -Isrc/include $(SRC_PATH)/test.cc ${SOURCES}
test:
	g++ -O3 -o test --std=c++17 -Isrc/include $(SRC_PATH)/test.cc ${SOURCES}
gputest:
	g++ -O3 -fopenacc -fcf-protection=none -no-pie -o test --std=c++17 -Isrc/include $(SRC_PATH)/test.cc ${SOURCES}
