PROJECT_ROOT:= $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
include ./ext/sdsl-lite/Make.helper

GCC_PARANOID=-pedantic -Wcast-align -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Werror -Winline -Wno-error=unused-parameter -Wno-error=unused-variable
CLANG_PARANOID=-pedantic -Weverything -Wno-c++98-compat

CXX_FLAGS=-std=c++11 -Wall -Wextra -DPROJECT_ROOT="\"$(PROJECT_ROOT)\"" -O3 -DNDEBUG

INCLUDES=-isystem$(INC_DIR)
LIB=$(LIB_DIR)/libsdsl.a $(LIB_DIR)/libdivsufsort.a $(LIB_DIR)/libdivsufsort64.a

OBJECTS=src/index.o src/graph.o
BINS=concatenate wanda-build wanda-assemble

%.o: %.cpp
	# @$(CXX) $(CXX_FLAGS) $(GCC_PARANOID) $(INCLUDES) -c $< -o $@
	@$(CXX) $(CXX_FLAGS) $(CLANG_PARANOID) $(INCLUDES) -c $< -o $@

all: $(OBJECTS) $(BINS)

wanda-build: src/wanda-build.cpp $(OBJECTS)
	@$(CXX) $(CXX_FLAGS) $(INCLUDES) -o wanda-build src/wanda-build.cpp $(OBJECTS) $(LIB)

wanda-assemble: src/wanda-assemble.cpp $(OBJECTS)
	@$(CXX) $(CXX_FLAGS) $(INCLUDES) -o wanda-assemble src/wanda-assemble.cpp $(OBJECTS) $(LIB)

concatenate: src/concatenate.cpp
	@$(CXX) $(CXX_FLAGS) $(INCLUDES) -o concatenate src/concatenate.cpp $(LIB)

clean:
	rm -rf $(OBJECTS) $(BINS) *.dSYM
