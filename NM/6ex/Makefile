.PHONY: all clean run

CC := g++  # Use the default g++ version on the system

CFLAGS = -std=c++11 -Wno-deprecated-enum-enum-conversion -O2 -Wall -Wextra
# Uncomment these if using OpenCV
# CFLAGS += `pkg-config --cflags opencv4`
# LDFLAGS = `pkg-config --libs opencv4`

SOURCES := $(wildcard *.cpp)
OBJECTS := $(patsubst %.cpp,%.o,$(SOURCES))
DEPENDS := $(patsubst %.cpp,%.d,$(SOURCES))

all: main

clean:
	$(RM) $(OBJECTS) $(DEPENDS) main

-include $(DEPENDS)

# Make object files
%.o: %.cpp Makefile
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

# Linking the executable from the object files
main: $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

# Optional run target to execute the program after building
run: main
	./main
