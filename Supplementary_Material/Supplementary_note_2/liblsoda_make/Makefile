include ../Make.inc

# ----------------------------
# Source files
# ----------------------------
SRCS =  \
cfode.c \
cfode_static.c \
common.c \
corfailure.c \
correction.c \
daxpy.c \
ddot.c \
dgefa.c \
dgesl.c \
dscal.c \
fnorm.c \
idamax.c \
intdy.c \
lsoda.c \
methodswitch.c \
orderswitch.c \
prja.c \
scaleh.c \
solsy.c \
stoda.c \
vmnorm.c \
strdup_printf.c

OBJS = $(SRCS:%.c=%.o)

# ----------------------------
# Compile object files
# ----------------------------
%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS-add) $(CFLAGS) -DNDEBUG -c $< -o $@

# ----------------------------
# Static library
# ----------------------------
STATIC_LIB = liblsoda.a

release debug: $(STATIC_LIB)

$(STATIC_LIB): $(OBJS)
	rm -f $@
	ar rcs $@ $^

# ----------------------------
# Clean
# ----------------------------
clean:
	rm -f *.o $(STATIC_LIB)

