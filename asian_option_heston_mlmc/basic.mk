VPATH = src include
SELF = basic.mk
CC = gcc
INC_FLAGS = -I./
INC_FLAGS += -I./include
INC_FLAGS += -I../../sdeWeakApproximation
INC_FLAGS += -I/opt/local/include
DEBUG = -02 -Wall
CFLAGS = $(DEBUG)
LIB_PATH = -L../
LIB_PATH += -L/opt/local/lib
LIBS = -lsde_wa
LIBS += -lgsl
LIBS += -lm
TARGET = asian_option_heston_mlmc
SRC = $(shell find . -name '*.c')
OBJS = $(SRC:.c=.o)

.SUFFIXES :
.SUFFIXES : .o .c
.c.o :%.c
	$(CC) $(CFLAGS) $(INC_FLAGS) -c $< -o $*.o

.PHONY: all
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(INC_FLAGS) $(LIB_PATH) $(LIBS) -o $@ $(OBJS) 

Makefile: $(SELF)
	rm -f $@
	cp $(SELF) $@
	chmod +w $@
	echo '# Automatically-generated dependencies list:' >>$@
	gcc ${CFLAGS} ${INC_FLAGS} -MM  \
	${SRC}	\
	>> $@

.PHONY: clean
clean :
	rm -f $(OBJS) $(TARGET)


