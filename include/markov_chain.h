#ifndef SIMPLEMM_H
#define SIMPLEMM_H
#include "read_gen.h"


#define LZERO  (-1.0E10) // log(0)
#define LSMALL (-0.5E10) // log values < LSMALL are set to LZERO
#define minLogExp -log(-LZERO) // ~=-23



// ====================================================== //
// ================= Simple Markob Chain ================ //
// ====================================================== //
typedef struct Simple_MarkovChain {
    double** trans_m;
    int n_states, order;
    char** s_name;
} smm_t;

typedef smm_t* smm_t_ptr;

/**
 * Craete a simple markov chain
 * */
smm_t_ptr SimpleMarkovChain(int, char**, int);
void init_smm(smm_t_ptr, int, int);
/**
 * Fit the training sequence
 * */
void simple_fit(smm_t_ptr, island_t_ptr, entry_ptr);
void calculate_trans_gene(char*, smm_t_ptr, char*, char**);
void calculate_gene(char*, smm_t_ptr, char*);
void calculate_transGene_higher(char*, char*, smm_t_ptr, char*);


void zero_order_fit(smm_t_ptr, island_t_ptr, entry_ptr);
void first_order_fit(smm_t_ptr, island_t_ptr, entry_ptr);
void second_order_fit(smm_t_ptr, island_t_ptr, entry_ptr);

/**
 * Return the probability of generating the sequence
 * */
double simple_genS(smm_t_ptr, island_t_ptr, entry_ptr);

double cal_sum(double, double(*)(double));
double zero_order_gen(smm_t_ptr, island_t_ptr, entry_ptr);
double first_order_gen(smm_t_ptr, island_t_ptr, entry_ptr);
double second_order_gen(smm_t_ptr, island_t_ptr, entry_ptr);

int indexOf(smm_t_ptr, char*, int);




// ====================================================== //
// ================= Hidden Markov Chain ================ //
// ====================================================== //




typedef struct Hidden_MarkovChain {
    double **trans_m;
    double **states_m;
    char **e_name;
    int n_states;
    int e_num;
} hmm_t;


typedef hmm_t* hmm_t_ptr;


hmm_t_ptr HiddenMarkovChain(int, int, char**, double**, double**);
hmm_t_ptr init_hmm(int, int, char**, double**, double**);



double logsum(double*, int);
void forward_sum(hmm_t_ptr, double*, double*, char*, char*);
double forward_algorithm(entry_ptr, island_t_ptr, hmm_t_ptr, double**);
void backward_sum(hmm_t_ptr, double*, double* , double* , char*, int, double**, int*);
double backward_algorithm(entry_ptr, island_t_ptr, hmm_t_ptr, double**);
void learning_algorithm(entry_ptr, island_t_ptr, hmm_t_ptr);


void reverse_str(char*,const char*, char*);

void hmm_learning(entry_ptr, island_t_ptr);
int indexOf_h(hmm_t_ptr, char*, int);


typedef struct EM_path {
    int *state_path;
    double P_x;
    int len;
} em_t;

typedef em_t* em_t_ptr;

em_t_ptr viterbi_algorithm(entry_ptr, island_t_ptr, hmm_t_ptr);
double max(double*, int, int*);













#endif