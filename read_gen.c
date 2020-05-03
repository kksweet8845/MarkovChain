#include "read_gen.h"
#include <string.h>


/**
 * freadline() - read a line unitil newline with a maxlen == 100
 * @chunk: A pointer to catch string
 * @input_stream: The input stream to read.
 */
int freadline(char *chunk, FILE *input_stream)
{
    memset(chunk, '\0', MAXLEN + 1);
    char *tmp = chunk;
    char single_char[2];

    int cur_index = 0;
    while (fread(single_char, sizeof(char), 1, input_stream) == 1 &&
           cur_index <= MAXLEN) {
        if (strcmp(single_char, "\n") != 0) {
            *tmp = single_char[0];
            cur_index++;
            tmp++;
        } else {
            break;
        }
    }
    if (!cur_index) {
        return -1;
    }
    return 0;
}



FILE *readFile(char *filename)
{
    FILE *ftp = fopen(filename, "r");
    return ftp;
}

/**
 * Extrac the pattern and return 0 for success. Otherwise, -1
 * @pattern: The regex to match
 * @string: The char* to be extract
 * @matchptr: The ptr to match the pointer
 * @nmatch: The int pointer to indicate the number of match
 * @cflags: The options
 */
int extract(const regex_t *compiled_pattern,
            char *string,
            regmatch_t **matchptr,
            int nmatch)
{
    int ok = 0;
    ok = regexec(compiled_pattern, string, nmatch, *matchptr, 0);
    return ok;
}

void initEntry(entry_ptr entry)
{
    (entry)->total_size = 10000;
    (entry)->size = 0;
    (entry)->content = malloc(sizeof(char *) * entry->total_size);
    (entry)->entryName = NULL;
}

void initIsland(island_t_ptr island)
{
    island->from = 0;
    island->to = 0;
    island->start_addr = island->close_addr = NULL;
}


void appendEntry(entry_ptr entry, const char *gene)
{
    if (entry->size >= entry->total_size) {
        entry->content = (char **) realloc(
            entry->content, sizeof(char *) * (entry->total_size) * 2);
        entry->total_size *= 2;
    }
    entry->content[(entry->size)] =
        (char *) malloc(sizeof(char) * (strlen(gene) + 1));
    strcpy(entry->content[(entry->size)++], gene);
}


// int calculateGene(const char* string, island_array arr, char* end){
//     char* tmp = string;
//     for(tmp=string;*tmp != '\0';tmp++){
//         arr[toLowercase(*tmp) - (int)'a']++;
//         if(end != NULL && tmp == end){
//             printf("Stop!\n");
//             return 1;
//         }
//     }
//     return 0;
// }



int extractEntry(const regex_t *entry_compiled_pattern,
                 FILE *ftp,
                 entry_ptr entry)
{
    char *buffer = malloc(sizeof(char) * (MAXLEN + 1));
    regmatch_t *matchptr = malloc(sizeof(regmatch_t));
    int ok = 0;
    int flag = 0;
    const char *gene_pattern = "^[atgcATGC]+";
    const char *another_entry_pattern = "^>[.]+";
    regex_t *gene_compiled_pattern = (regex_t *) malloc(sizeof(regex_t));
    regex_t *another_compiled_pattern = (regex_t *) malloc(sizeof(regex_t));
    // memeory leak problem- bug!
    ok = regcomp(gene_compiled_pattern, gene_pattern, REG_EXTENDED | REG_ICASE);
    ok = regcomp(another_compiled_pattern, another_entry_pattern,
                 REG_EXTENDED | REG_ICASE);
    initEntry(entry);
    while (freadline(buffer, ftp) != EOF) {
        switch (flag) {
        case 0:
            ok = extract(entry_compiled_pattern, buffer, &matchptr, 1);
            // zero for match
            if (!ok) {
                flag = 1;
                entry->entryName = malloc(sizeof(char) * (strlen(buffer) + 1));
                strcpy(entry->entryName, buffer);
                printf("%s\n", entry->entryName);
            }
            break;
        case 1:
            ok = extract(gene_compiled_pattern, buffer, &matchptr, 1);
            // printf("%s\n", buffer);
            if (!ok) {
                appendEntry(entry, buffer);
            }
            break;
        }
        if (flag == 1 &&
            !extract(another_compiled_pattern, buffer, &matchptr, 1)) {
            break;
        }
    }


    free(buffer);
    free(matchptr);
    free(gene_compiled_pattern);
    free(another_compiled_pattern);
    return 0;
}


int entryToTxt(const entry_ptr entry, char *filename)
{
    FILE *fpt = fopen(filename, "w");
    fprintf(fpt, "%s\n", entry->entryName);
    for (int i = 0; i < entry->size; i++) {
        fprintf(fpt, "%s\n", entry->content[i]);
    }
    return 0;
}


char toLowercase(char alpha)
{
    if (alpha >= 65 && alpha <= 90) {
        alpha += 32;
    }
    return alpha;
}


island_t_ptr extractS(entry_ptr entry,
                      int from,
                      int to,
                      char *startPattern,
                      char *closePattern)
{
    int count = 0;

    regex_t *start_compiled_pattern = (regex_t *) malloc(sizeof(regex_t));
    regex_t *close_compiled_pattern = (regex_t *) malloc(sizeof(regex_t));

    int ok = 0;
    ok =
        regcomp(start_compiled_pattern, startPattern, REG_EXTENDED | REG_ICASE);
    ok =
        regcomp(close_compiled_pattern, closePattern, REG_EXTENDED | REG_ICASE);

    regmatch_t *start_pos = (regmatch_t *) malloc(sizeof(regmatch_t));
    regmatch_t *close_pos = (regmatch_t *) malloc(sizeof(regmatch_t));

    island_t_ptr res = malloc(sizeof(island_t));
    initIsland(res);
    while (from <= to) {
        ok = extract(start_compiled_pattern, entry->content[from], &start_pos,
                     1);
        if (!ok) {
            printf("Find start\n");
            printf("%s\n", entry->content[from]);
            break;
        }
        from++;
    }
    res->from = from;
    res->start_addr = entry->content[from] + start_pos->rm_so;
    while (from <= to) {
        ok = extract(close_compiled_pattern, entry->content[from], &close_pos,
                     1);
        if (!ok) {
            printf("Find finish\n");
            printf("%s\n", entry->content[from]);
            break;
        }
        from++;
    }
    res->to = from;
    res->close_addr = entry->content[to] + close_pos->rm_eo;
    return res;
}


island_t_ptr extractS_from_file(char *name,
                                char *filename,
                                entry_ptr entry,
                                int from,
                                int to,
                                char *start_pattern,
                                char *end_pattern)
{
    regex_t *start_compiled_pattern = (regex_t *) malloc(sizeof(regex_t));
    regex_t *close_compiled_pattern = (regex_t *) malloc(sizeof(regex_t));

    int ok = 0;
    ok = regcomp(start_compiled_pattern, start_pattern,
                 REG_EXTENDED | REG_ICASE);
    ok = regcomp(close_compiled_pattern, end_pattern, REG_EXTENDED | REG_ICASE);
    regmatch_t *start_pos = (regmatch_t *) malloc(sizeof(regmatch_t));
    regmatch_t *close_pos = (regmatch_t *) malloc(sizeof(regmatch_t));
    printf("=================================================\n");
    FILE *obj_island = fopen(filename, "r");
    island_t_ptr res_island = NULL;
    if (!obj_island) {
        printf("Generating island...\n");
        res_island = extractS(entry, from, to, start_pattern, end_pattern);
        FILE *obj_island = fopen(filename, "w");
        fwrite(res_island, sizeof(island_t), 1, obj_island);
        printf("Save the obj as %s\n", filename);
        fclose(obj_island);
    } else {
        res_island = malloc(sizeof(island_t));
        while (fread(res_island, sizeof(island_t), 1, obj_island))
            ;
        fclose(obj_island);
        ok = extract(start_compiled_pattern, entry->content[res_island->from],
                     &start_pos, 1);
        ok = extract(close_compiled_pattern, entry->content[res_island->to],
                     &close_pos, 1);
        res_island->start_addr =
            entry->content[res_island->from] + start_pos->rm_so;
        res_island->close_addr = entry->content[to] + close_pos->rm_eo;
        printf("Finish loading from %s\n", filename);
    }
    printf("%s\n", name);
    printf("Start from line : %d\n", res_island->from);
    printf("End to line : %d\n", res_island->to);
    printf("Start position: %c\n", *res_island->start_addr);
    printf("End position: %c\n", *res_island->close_addr);

    free(start_compiled_pattern);
    free(close_compiled_pattern);
    free(start_pos);
    free(close_pos);

    return res_island;
}