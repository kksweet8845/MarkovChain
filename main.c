#include "read_gen.h"
#include "markov_chain.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>


int main(void){

    char* srcname = "./data/NC_000006.12_chromosome_6.txt";
    // char* srcname = "./data/custimized_.txt";
    FILE* ftp = readFile(srcname);

    char* buffer = malloc(sizeof(char)*(MAXLEN+1));
    const char* entry_pattern = ">NC_000006.12 Homo sapiens chromosome 6, GRCh38.p13 Primary Assembly";
    regex_t* compiled_pattern = (regex_t*) malloc(sizeof(regex_t));
    int ok = 0;
    ok = regcomp(compiled_pattern, entry_pattern, REG_EXTENDED);

    entry_ptr nc1_entry = malloc(sizeof(entry_t));

    extractEntry(compiled_pattern, ftp, nc1_entry);
    char* start_pattern = "ttggtaccat";
    char* end_pattern = "CTTTGCCTG";
    island_t_ptr res_island = extractS(nc1_entry,
             100000,
             199999,
             start_pattern,
             end_pattern);
    printf("Start from line : %d\n", res_island->from);
    printf("End to line : %d\n", res_island->to);
    printf("Start position: %c\n", *res_island->start_addr);
    printf("End position: %c\n", *res_island->close_addr);
    char* tmp;
    int i;
// #define SIMPLE
#define HIDDEN

#ifdef SIMPLE
/* ====================================================== */
/* ===================== Zero Order ===================== */
/* ====================================================== */
    char* states_name[4] = {"a","t","c","g"};
    // printf("%s %s %s %s\n", states_name[0], states_name[1], states_name[2], states_name[3]);
    smm_t_ptr zero_smm = SimpleMarkovChain(4, states_name, 0);
    simple_fit(zero_smm, res_island, nc1_entry);
    double prob = simple_genS(zero_smm, res_island, nc1_entry);
    for(int i=0;i<zero_smm->n_states;i++){
        printf("%lf ", zero_smm->trans_m[0][i]);
    }
    printf("\n");
    printf("probability : %lf\n", prob);

/* ====================================================== */
/* ===================== First Order ==================== */
/* ====================================================== */
    char* states_name_first[4] = { "a", "t", "c", "g" };
    smm_t_ptr first_smm = SimpleMarkovChain(4, states_name_first, 1);
    simple_fit(first_smm, res_island, nc1_entry);
    prob = simple_genS(first_smm, res_island, nc1_entry);

    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            printf("%lf ", first_smm->trans_m[i][j]);
        }
        printf("\n");
    }

    printf("first probability : %lf\n", prob);

    char* states_name_second[16] ={
        "aa", "at", "ac", "ag",
        "ta", "tt", "tc", "tg",
        "ca", "ct", "cc", "cg",
        "ga", "gt", "gc", "gg"
    };
    smm_t_ptr second_smm = SimpleMarkovChain(16, states_name_second, 2);
    simple_fit(second_smm, res_island, nc1_entry);
    prob = simple_genS(second_smm, res_island, nc1_entry);
    for(int i=0;i<16;i++){
        for(int j=0 ;j<16;j++){
            // printf("%s to %s : %lf\n", states_name_second[i], states_name_second[j], second_smm->trans_m[i][j]);
            if(i == -1 && j >=0){
                printf("%12s", states_name_second[j]);
            }else if(i >= 0){
                if(j == -1){
                    printf("%12s", states_name_second[i]);
                }
                printf("%12.4le", second_smm->trans_m[i][j]);
            }
        }
        printf("\n");
    }
    printf("Second probability : %lf\n", prob);

#endif

/* ====================================================== */
/* ================= Hidden Markov Chain ================ */
/* ====================================================== */
#ifdef HIDDEN

    char* e_name[4] = {"a","t","c","g"};
    double trans_m_22[2][2] = {
        {0.8, 0.2},
        {0.1, 0.9}
    };
    double **trans_m = malloc(sizeof(sizeof(double*) * 2));
    for(int i=0;i<2;i++){
        trans_m[i] = malloc(sizeof(double)*2);
    }
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            trans_m[i][j] = trans_m_22[i][j];
        }
    }
    double states_m_22[2][4] = {
        {0.2,0.2,0.3,0.3},
        {0.4,0.4,0.1,0.1}
    };

    double **states_m = malloc(sizeof(double*) * 2);
    for(int i=0;i<2;i++){
        states_m[i] = malloc(sizeof(double)*4);
    }
    for(int i=0;i<2;i++){
        for(int j=0;j<4;j++){
            states_m[i][j] = states_m_22[i][j];
        }
    }
    // printf("trans_m[0] : %p\n", &trans_m[0][1]);
    hmm_t_ptr hmm = HiddenMarkovChain(2, 4, e_name, trans_m, states_m);
    double for_prob = forward_algorithm(nc1_entry, res_island, hmm);
    double bac_prob = backward_algorithm(nc1_entry, res_island, hmm);
    printf("Forward prob: %le\nBackward prob: %le\n", for_prob, bac_prob);
#endif
    return 0;
}
