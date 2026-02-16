# Top-level directories
BIN_DIR = bin
DEP_DIR = dep
INC_DIR = include
OBJ_DIR = obj
SRC_DIR = src

# Subdirecotires
INC_SUBDIRS = $(shell find $(INC_DIR) -type d)
SRC_SUBDIRS = $(shell find $(SRC_DIR) -type d)

# SRC = $(shell find $(SRC_DIR) -name '*.cpp')
SRC = $(shell find $(SRC_DIR) -name '*.cpp')
OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
DEP = $(SRC:$(SRC_DIR)/%.cpp=$(DEP_DIR)/%.d)
TARGET = $(BIN_DIR)/main.exe

LIB_SRC := $(filter-out $(SRC_DIR)/main.cpp, $(SRC))
LIB_OBJ := $(filter-out $(OBJ_DIR)/main.o, $(OBJ))

TEST_DIR := tests
TEST_TARGET := $(BIN_DIR)/run_test.exe
TEST_SRC := $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJ := $(TEST_SRC:$(TEST_DIR)/%.cpp=$(OBJ_DIR)/$(TEST_DIR)/%.o)
TEST_DEP := $(TEST_SRC:$(TEST_DIR)/%.cpp=$(DEP_DIR)/$(TEST_DIR)/%.d)

INC = $(addprefix -I,$(INC_SUBDIRS))
CXX = g++
CXXFLAGS = -Wall -O2 -MMD -MP -std=c++17 $(INC) -Iexternal # directory to external codes, add more if needed

.PHONY: all clean test proj1

# make all
all: $(TARGET)

$(TARGET): $(OBJ)
	@mkdir -p $(BIN_DIR)
	@$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@) $(DEP_DIR)/$(dir $*)
	@$(CXX) $(CXXFLAGS) -c $< -o $@ -MF $(DEP_DIR)/$*.d

-include $(DEP)

# make test
test: $(TEST_TARGET)
	@$(TEST_TARGET)

$(TEST_TARGET): $(TEST_OBJ) $(LIB_OBJ)
	@mkdir -p $(BIN_DIR)
	@$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJ_DIR)/$(TEST_DIR)/%.o: $(TEST_DIR)/%.cpp
	@mkdir -p $(dir $@) $(DEP_DIR)/$(TEST_DIR)
	@$(CXX) $(CXXFLAGS) -c $< -o $@ -MF $(DEP_DIR)/$(TEST_DIR)/$*.d

# make clean
clean:
	@rm -rf $(OBJ_DIR) $(BIN_DIR) $(DEP_DIR)
