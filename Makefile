# simple makefile

MKDIR_P = mkdir -p

# main directories
BIN_DIR = ./bin
OBJ_DIR = ./build
SRC_DIR = ./src
INC_DIR = ./include
LIB_DIR = ./libs

# sub project directories
INC_DIRS = $(shell find $(INC_DIR) -type d)

# Default to using system's default version of python
PYTHON_BIN     ?= python3
PYTHON_CONFIG  := $(PYTHON_BIN)-config
PYTHON_INCLUDE ?= $(shell $(PYTHON_CONFIG) --includes)
EXTRA_FLAGS    := $(PYTHON_INCLUDE)
# NOTE: Since python3.8, the correct invocation is `python3-config --libs --embed`. 
# So of course the proper way to get python libs for embedding now is to
# invoke that, check if it crashes, and fall back to just `--libs` if it does.
LDFLAGS        += $(shell if $(PYTHON_CONFIG) --libs --embed >/dev/null; then $(PYTHON_CONFIG) --libs --embed; else $(PYTHON_CONFIG) --libs; fi)

# Either finds numpy or set -DWITHOUT_NUMPY
EXTRA_FLAGS     += $(shell $(PYTHON_BIN) $(CURDIR)/numpy_flags.py)
WITHOUT_NUMPY   := $(findstring $(EXTRA_FLAGS), WITHOUT_NUMPY)

# compiler and linker options
EXE_NAME = field_trans

CXX = g++

CXX_FLAGS = -W -Wall -Wextra -s -g -O3 -static -std=c++17 -DLIN64
INC_FLAGS = $(addprefix -I,$(INC_DIRS))
# LD_FLAGS = -static
LIB_FLAGS = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -lm -ldl -rdynamic -lpython3.6m

INC_FLAGS += $(PYTHON_INCLUDE)
INC_FLAGS += $(EXTRA_FLAGS)

# collect sources ...
SRCS = $(shell find $(SRC_DIR) -name "*.cpp")
OBJS = $(SRCS:%.cpp=$(OBJ_DIR)/%.o)
#DEPS = $(OBJS:.o=.d)


# rules for c++ files
$(OBJ_DIR)/%.o: %.cpp
		$(MKDIR_P) $(dir $@)
		$(CXX) $(INC_FLAGS) $(CXX_FLAGS) -c $< -o $@

#rules for target
$(BIN_DIR)/$(EXE_NAME): $(OBJS)
		$(MKDIR_P) $(BIN_DIR)
		$(CXX) $(LD_FLAGS) $(OBJS) $(LIBS) $(LIB_FLAGS) -o $@
# $(CXX) $(LD_FLAGS) $(OBJS) $(LIBS) $(LIB_FLAGS) -o $@

clean:
		rm -rf $(OBJ_DIR)/*
		rm -rf $(BIN_DIR)/*