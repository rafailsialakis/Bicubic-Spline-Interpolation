
# ============================================================
#  Makefile — Parallel Bicubic Spline Interpolation (C++)
#  Usage:
#    make              → build (release)
#    make debug        → build with debug symbols & sanitizers
#    make run          → build + run with default args
#    make clean        → remove build artefacts
#    make stb          → download stb headers automatically
# ============================================================
 
CXX      := g++
TARGET   := bin/bsi_parallel
SRC      := Cpp/bsi_parallel.cpp
 
CXXFLAGS_COMMON := -std=c++17 -Wall -Wextra -fopenmp
CXXFLAGS_REL    := -O2
CXXFLAGS_DBG    := -O0 -g -fsanitize=address,undefined
 
INCLUDES        := -I includes
CXXFLAGS_COMMON := -std=c++17 -Wall -Wextra -fopenmp $(INCLUDES) \
                   -Wno-missing-field-initializers   # silence stb warnings
CXXFLAGS_REL    := -O2
CXXFLAGS_DBG    := -O0 -g -fsanitize=address,undefined
 
LDFLAGS := -lm

 .PHONY: all debug run clean stb
 
all: $(TARGET)
 
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS_COMMON) $(CXXFLAGS_REL) -o $@ $^ $(LDFLAGS)
	@echo "Built $(TARGET) (release)"
 
run: all
	@mkdir -p Images
	./$(TARGET) $(INPUT) $(OUTPUT) $(SCALE)

clean:
	rm -f $(TARGET) $(TARGET)_debug