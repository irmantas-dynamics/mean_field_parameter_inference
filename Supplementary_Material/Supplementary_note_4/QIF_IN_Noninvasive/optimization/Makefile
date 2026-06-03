# ----------------------------
# Compiler and flags
# ----------------------------
CC = gcc
SRC_DIR := ../../liblsoda-master/src           # path to sourcecode of liblsoda
CFLAGS = -O3 -fopenmp -fPIC -Wall -I$(SRC_DIR)
LDFLAGS = -fopenmp -lm -lgsl -lgslcblas

# LSODA library location (do not leave white space after directory name)
LSODA_DIR := ../../liblsoda-master/src
LSODA_STATIC = $(LSODA_DIR)/liblsoda.a       # static library
#LSODA_SHARED = -L$(LSODA_DIR) -llsoda       # shared library alternative

# Source files
SRCS_LIB = losslib.c
OBJS_LIB = $(SRCS_LIB:.c=.o)

SRCS_TEST = test.c
OBJS_TEST = $(SRCS_TEST:.c=.o)

# ----------------------------
# Default target
# ----------------------------
all: liblosslib.so 
#test_program

# ----------------------------
# Build shared library for Python
# ----------------------------
liblosslib.so: $(OBJS_LIB) $(LSODA_STATIC)
	$(CC) -shared -o $@ $(OBJS_LIB) $(LSODA_STATIC) $(LDFLAGS)

# If LSODA is a shared library instead of static, use:
#liblosslib.so: $(OBJS_LIB)
#	$(CC) -shared -o $@ $(OBJS_LIB) $(LSODA_SHARED) $(LDFLAGS)

# ----------------------------
# Build test executable
# ----------------------------
#test_program: $(OBJS_TEST) $(OBJS_LIB) $(LSODA_STATIC)
#	$(CC) -o $@ $^ $(LDFLAGS)

# ----------------------------
# Compile object files
# ----------------------------
#%.o: %.c
#	$(CC) $(CFLAGS) -c $< -o $@

# ----------------------------
# Clean
# ----------------------------
clean:
	rm -f *.o liblosslib.so

.PHONY: all clean

