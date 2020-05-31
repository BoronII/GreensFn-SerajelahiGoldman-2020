# This is a makefile I wrote to play around with things 
CC=gcc

# compiler flags: 
# -g adds debugging information to the executable file 
# -Wall turns on most, but not all, compiler warnings 
CFLAGS = -g -Wall 

# the build target executable:
TARGET=gf

$(TARGET): GreensFn.cpp 
	$(CC) $(CFLAGS) -o $(TARGET) GreensFn.cpp -lm -g

clean:
	$(RM) $(TARGET)
