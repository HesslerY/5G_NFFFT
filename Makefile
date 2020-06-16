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

# compiler and linker options
EXE_NAME = field_trans

CXX = g++

CXX_FLAGS = -W -Wall -Wextra -s -g -O3 -static -std=c++17 -DLIN64
INC_FLAGS = $(addprefix -I,$(INC_DIRS))
LD_FLAGS = -static
# LIB_FLAGS = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -lm -ldl -rdynamic -lpython3.6m
LIB_FLAGS = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -lm -ldl -rdynamic 

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