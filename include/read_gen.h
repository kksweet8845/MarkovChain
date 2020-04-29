#ifndef READ_GEN_H
#define READ_GEN_H

#include<regex.h>
#include <stdlib.h>
#include <stdio.h>

#define MAXLEN 100



typedef struct entry {
    char** content;
    unsigned int total_size;
    unsigned int size;
    char* entryName;

} entry_t;



typedef struct Island {
    int from, to;
    char* start_addr, *close_addr;
} island_t;

typedef island_t* island_t_ptr;


typedef entry_t* entry_ptr;

typedef int* island_array;


FILE* readFile(char*);
int freadline(char*, FILE*);
int extract(const regex_t*, char*, regmatch_t**, int);
void initEntry(entry_ptr);
void appendEntry(entry_ptr, const char*);
int extractEntry(const regex_t*, FILE*, entry_ptr);
int entryToTxt(const entry_ptr, char*);
char toLowercase(char);
int calculateGene(const char*, island_array, char*);
void initIsland(island_t_ptr);
island_t_ptr extractS(entry_ptr, int, int, char*, char*);
// void calculateIsland(entry_ptr, island_t_ptr);

#endif
