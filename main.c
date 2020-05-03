#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "markov_chain.h"
#include "read_gen.h"


int main(void)
{
    char *srcname = "./data/NC_000006.12_chromosome_6.txt";
    // char* srcname = "./data/custimized_.txt";
    FILE *ftp = readFile(srcname);

    char *buffer = malloc(sizeof(char) * (MAXLEN + 1));
    const char *entry_pattern =
        ">NC_000006.12 Homo sapiens chromosome 6, GRCh38.p13 Primary Assembly";
    regex_t *compiled_pattern = (regex_t *) malloc(sizeof(regex_t));
    int ok = 0;
    ok = regcomp(compiled_pattern, entry_pattern, REG_EXTENDED);

    entry_ptr nc1_entry = malloc(sizeof(entry_t));

    extractEntry(compiled_pattern, ftp, nc1_entry);
    // char* start_pattern = "tt";
    // char* end_pattern = "ca";

    int from = 100000;
    int to = 199999;
    // int from = 0;
    // int to = 0;
    // ? =========================================================

    char *start_pattern = "ttggtaccat";
    char *end_pattern = "CTTTGCCTG";

    char *filename = "obj_island_1.obj";
    island_t_ptr res_island =
        extractS_from_file("res 1 island", filename, nc1_entry, from, to,
                           start_pattern, end_pattern);

    // ? =========================================================
    char *start_pattern2 = "gagatccttc";
    char *end_pattern2 = "ctttaaaagaaaa";
    char *filename2 = "obj_island_2.obj";
    island_t_ptr res_island2 =
        extractS_from_file("res 2 island", filename2, nc1_entry, from, to,
                           start_pattern2, end_pattern2);

    printf("========================================================\n");
    printf("===================== Summary ==========================\n");
    printf("========================================================\n\n");
#define SIMPLE
#define HIDDEN

#ifdef SIMPLE
    /* ====================================================== */
    /* ===================== Zero Order ===================== */
    /* ====================================================== */
    printf("Simple Markov Chain\n");

    char *states_name[4] = {"a", "t", "c", "g"};
    // printf("%s %s %s %s\n", states_name[0], states_name[1], states_name[2],
    // states_name[3]);
    smm_t_ptr zero_smm = SimpleMarkovChain(4, states_name, 0);
    simple_fit(zero_smm, res_island, nc1_entry);
    double prob = simple_genS(zero_smm, res_island, nc1_entry);
    printf("Zero order Markov Chain\n");
    printf("A T C G sta\n");
    for (int i = 0; i < zero_smm->n_states; i++) {
        printf("%12s", states_name[i]);
    }
    printf("\n");
    for (int i = 0; i < zero_smm->n_states; i++) {
        printf("%12e ", zero_smm->trans_m[0][i]);
    }
    printf("\n");
    printf("zero order probability : %lf\n", prob);
    printf("\n\n\n\n");
    /* ====================================================== */
    /* ===================== First Order ==================== */
    /* ====================================================== */
    char *states_name_first[4] = {"a", "t", "c", "g"};
    smm_t_ptr first_smm = SimpleMarkovChain(4, states_name_first, 1);
    simple_fit(first_smm, res_island, nc1_entry);
    prob = simple_genS(first_smm, res_island, nc1_entry);

    printf("First order Markov Chain\n");
    printf("%8s", "");
    for (int i = 0; i < 4; i++) {
        printf("%8s", states_name_first[i]);
    }
    printf("\n");
    for (int i = 0; i < 4; i++) {
        printf("%-12s", states_name_first[i]);
        for (int j = 0; j < 4; j++) {
            printf("%12e ", first_smm->trans_m[i][j]);
        }
        printf("\n");
    }

    printf("first probability : %lf\n", prob);

    printf("\n\n\n\n");

    /* ====================================================== */
    /* ==================== Second Order ==================== */
    /* ====================================================== */

    char *states_name_second[16] = {"aa", "at", "ac", "ag", "ta", "tt",
                                    "tc", "tg", "ca", "ct", "cc", "cg",
                                    "ga", "gt", "gc", "gg"};
    smm_t_ptr second_smm = SimpleMarkovChain(16, states_name_second, 2);
    simple_fit(second_smm, res_island, nc1_entry);
    prob = simple_genS(second_smm, res_island, nc1_entry);

    printf("Second order Markov Chain\n");
    printf("%12s", "");
    for (int i = 0; i < 16; i++) {
        printf("%12s", states_name_second[i]);
    }
    printf("\n");
    for (int i = 0; i < 16; i++) {
        printf("%12s", states_name_second[i]);
        for (int j = 0; j < 16; j++) {
            // printf("%s to %s : %lf\n", states_name_second[i],
            // states_name_second[j], second_smm->trans_m[i][j]);
            if (i == -1 && j >= 0) {
                printf("%12s", states_name_second[j]);
            } else if (i >= 0) {
                if (j == -1) {
                    printf("%12s", states_name_second[i]);
                }
                printf("%12.4le", second_smm->trans_m[i][j]);
            }
        }
        printf("\n");
    }
    printf("Second probability : %lf\n", prob);

    printf("\n\n");
#endif

/* ====================================================== */
/* ================= Hidden Markov Chain ================ */
/* ====================================================== */
#ifdef HIDDEN

#define num_states 4


    printf("Hidden Markov Chain\n");

    char *e_name[4] = {"a", "t", "c", "g"};
    double trans_m_22[num_states][num_states] = {
        {0.8, 0.01, 0.07, 0.01, 0.01},  {0.01, 0.9, 0.02, 0.05, 0.02},
        {0.01, 0.06, 0.9, 0.01, 0.02},  {0.013, 0.024, 0.06, 0.9, 0.003},
    };
    double **trans_m = (double **) malloc(sizeof(double *) * num_states);
    for (int i = 0; i < num_states; i++) {
        trans_m[i] = malloc(sizeof(double) * 4);
    }
    for (int i = 0; i < num_states; i++) {
        for (int j = 0; j < num_states; j++) {
            trans_m[i][j] = trans_m_22[i][j];
        }
    }
    double states_m_22[num_states][4] = {
        {0.2, 0.2, 0.3, 0.3}, {0.4, 0.4, 0.1, 0.1},    {0.9, 0.01, 0.01, 0.08},
        {0.4, 0.4, 0.1, 0.1},
    };

    double **states_m = malloc(sizeof(double *) * num_states);
    for (int i = 0; i < num_states; i++) {
        states_m[i] = malloc(sizeof(double) * 4);
    }
    for (int i = 0; i < num_states; i++) {
        for (int j = 0; j < 4; j++) {
            states_m[i][j] = states_m_22[i][j];
        }
    }
    // printf("trans_m[0] : %p\n", &trans_m[0][1]);
    hmm_t_ptr hmm = HiddenMarkovChain(num_states, 4, e_name, trans_m, states_m);
    double for_prob = forward_algorithm(nc1_entry, res_island, hmm, NULL);
    double bac_prob = backward_algorithm(nc1_entry, res_island, hmm, NULL);
    printf("Forward prob: %e\nBackward prob: %e\n", for_prob, bac_prob);
    em_t_ptr em_1 = viterbi_algorithm(nc1_entry, res_island, hmm);

    FILE *output_fpt = fopen("path_output.txt", "w");


    for (int i = 0; i < em_1->len; i++) {
        fprintf(output_fpt, "%d", em_1->state_path[i]);
    }

    printf("\n");
    printf("After applying learning algorithm\n");
    learning_algorithm(nc1_entry, res_island, hmm);

    printf("State transisiton matrix\n");
    for (int i = 0; i < num_states; i++) {
        for (int j = 0; j < num_states; j++) {
            printf("%8e ", hmm->trans_m[i][j]);
        }
        printf("\n");
    }
    printf("State gene matrix\n");
    for (int i = 0; i < num_states; i++) {
        for (int j = 0; j < 4; j++) {
            printf("%8e ", hmm->states_m[i][j]);
        }
        printf("\n");
    }
    for_prob = forward_algorithm(nc1_entry, res_island, hmm, NULL);
    bac_prob = backward_algorithm(nc1_entry, res_island, hmm, NULL);
    printf("Forward prob: %e\nBackward prob: %e\n", for_prob, bac_prob);

    em_t_ptr em_2 = viterbi_algorithm(nc1_entry, res_island, hmm);

    FILE *output2_fpt = fopen("path_output2.txt", "w");


    for (int i = 0; i < em_1->len; i++) {
        fprintf(output2_fpt, "state : %d\n", em_2->state_path[i]);
    }

    printf("================================================\n");
    printf("Another 100k base section\n");
    for_prob = forward_algorithm(nc1_entry, res_island2, hmm, NULL);
    printf("Forward prob: %e\n", for_prob);
#endif
    return 0;
}
