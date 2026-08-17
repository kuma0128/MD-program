CXX ?= c++
AR ?= ar
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -I.

BUILD_DIR := build
LIB_SOURCES := EWALD/ewald.cpp MC/monte_carlo.cpp RATTLE/rattle.cpp ROTATE/rotate.cpp SHAKE/shake.cpp
LIB_OBJECTS := $(LIB_SOURCES:%.cpp=$(BUILD_DIR)/%.o)

.PHONY: all test clean

all: $(BUILD_DIR)/libmd.a $(BUILD_DIR)/gaussian $(BUILD_DIR)/md_tests

test: $(BUILD_DIR)/md_tests
	$(BUILD_DIR)/md_tests

$(BUILD_DIR)/gaussian: GAUSS/main.cpp GAUSS/gaussian.hpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

$(BUILD_DIR)/libmd.a: $(LIB_OBJECTS)
	$(AR) rcs $@ $^

$(BUILD_DIR)/md_tests: tests/test_md.cpp $(BUILD_DIR)/libmd.a | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $@

clean:
	rm -rf $(BUILD_DIR)
