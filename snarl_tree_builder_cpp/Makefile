# Makefile for Snarl Tree Builder
# Simple direct compilation without CMake

CXX = g++
CXXFLAGS = -std=c++17 -O3 -pthread
LIBBDSG_DIR = ../libbdsg

INCLUDES = -I$(LIBBDSG_DIR)/bdsg/include \
           -I$(LIBBDSG_DIR)/bdsg/deps/libhandlegraph/src/include \
           -I$(LIBBDSG_DIR)/bdsg/deps/sdsl-lite/include \
           -I$(LIBBDSG_DIR)/bdsg/deps/sparsepp

LDFLAGS = -L$(LIBBDSG_DIR)/lib -L$(LIBBDSG_DIR)
LIBS = -lbdsg -lhandlegraph -lsdsl -ldivsufsort -ldivsufsort64
RPATH = -Wl,-rpath,$(LIBBDSG_DIR):$(LIBBDSG_DIR)/lib

TARGETS = construct_snarl_tree_map visualize_snarl_tree chunk_point_finding
SOURCE1 = construct_snarl_tree_map.cpp
SOURCE2 = visualize_snarl_tree.cpp
SOURCE3 = chunk_point_finding.cpp

.PHONY: all clean help

all: $(TARGETS)

construct_snarl_tree_map: $(SOURCE1)
	@echo "Compiling construct_snarl_tree_map..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o construct_snarl_tree_map $(SOURCE1) $(LDFLAGS) $(LIBS) $(RPATH)
	@echo "Build complete: construct_snarl_tree_map"

visualize_snarl_tree: $(SOURCE2)
	@echo "Compiling visualize_snarl_tree..."
	$(CXX) $(CXXFLAGS) -o visualize_snarl_tree $(SOURCE2)
	@echo "Build complete: visualize_snarl_tree"

chunk_point_finding: $(SOURCE3)
	@echo "Compiling chunk_point_finding..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o chunk_point_finding $(SOURCE3) $(LDFLAGS) $(LIBS) $(RPATH)
	@echo "Build complete: chunk_point_finding"

clean:
	@echo "Cleaning build artifacts..."
	rm -f $(TARGETS)
	rm -rf build
	@echo "Clean complete"

help:
	@echo "Snarl Tree Builder - Build Instructions"
	@echo ""
	@echo "Usage:"
	@echo "  make        - Build the executable"
	@echo "  make clean  - Remove build artifacts"
	@echo "  make help   - Show this help message"
	@echo ""
	@echo "Requirements:"
	@echo "  - g++ with C++17 support"
	@echo "  - libbdsg installed in ../libbdsg/"

