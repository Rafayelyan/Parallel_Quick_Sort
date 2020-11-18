CXX:=mpicxx

CXX_OPTS:=

SRCS:=$(wildcard src/*.cpp)

HEADERS:=$(wildcard src/*h)

OBJECTS:=$(patsubst src/%.cpp,bin/%.o,$(SRCS))

bin/myprogram : $(OBJECTS)
	@$(CXX) $^ $(CXX_OPTS) -o $@

bin/%.o : src/%.cpp $(HEADERS)
	@$(CXX) $< $(CXX_OPTS) -c -o $@


