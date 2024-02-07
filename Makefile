CXX = g++
CXX_FLAGS = -std=c++17 -Wall -Wextra -Werror -pedantic -O3

SRC_DIR = src
INC_DIR = include
BUILD_DIR = build

TARGET = $(BUILD_DIR)/libooz.a

SRC_FILES = $(wildcard $(SRC_DIR)/*.c)
SRC_FILES += $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRC_FILES))
OBJ_FILES += $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRC_FILES))

$(TARGET): $(OBJ_FILES)
	ar rcs $@ $^

$(BUILD_DIR)/%.o: $(SRC_FILES)
	$(shell mkdir -p $(dir $@))
	$(CXX) $(CXX_FLAGS) -I$(INC_DIR) -c -o $@ $<

.PHONY: clean

clean:
	rm -rf $(BUILD_DIR)
