CC:= gcc
override CFLAGS += -O3 -Wall -g -I./include
TARGET = main

SRC_DIR = ./
SRC 	= $(wildcard $(SRC_DIR)*.c)
OBJ 	= $(SRC:%.c=%.o)



$(TARGET) : $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@ -lm

%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@ -lm


clean:
	rm *.o
