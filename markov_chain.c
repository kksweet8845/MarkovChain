#include "markov_chain.h"
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* ====================================================== */
/* ================= Simple Markov Chain ================ */
/* ====================================================== */

smm_t_ptr SimpleMarkovChain(int states, char **states_name, int order)
{
    smm_t *sm = (smm_t *) malloc(sizeof(smm_t));
    init_smm(sm, states, order);
    sm->s_name = states_name;

    return sm;
}

void init_smm(smm_t_ptr smm, int states, int order)
{
    smm->order = order;
    // init first order
    smm->trans_m = (double **) malloc(sizeof(double *) * states);
    for (int i = 0; i < states; i++) {
        smm->trans_m[i] = (double *) malloc(sizeof(double) * states);
        for (int j = 0; j < states; j++) {
            smm->trans_m[i][j] = 0.0;
        }
    }

    smm->n_states = states;
}



void calculate_trans_gene(char *str, smm_t_ptr smm, char *end, char **prev)
{
    char *cur;

    if (*prev != NULL) {
        smm->trans_m[indexOf(smm, *prev, 1)][indexOf(smm, str, 1)]++;
    }
    for (cur = str; *cur != '\0'; cur++) {
        if (end != NULL && cur == end) {
            return;
        }
        smm->trans_m[indexOf(smm, cur, 1)][indexOf(smm, cur + 1, 1)]++;
    }
    *prev = cur - 1;
}

void calculate_gene(char *str, smm_t_ptr smm, char *end)
{
    char *tmp = str;
    for (tmp = str; *tmp != '\0'; tmp++) {
        if (end != NULL && tmp == end) {
            return;
        }
        smm->trans_m[0][indexOf(smm, tmp, 1)]++;
    }
}

void calculate_transGene_higher(char *str,
                                char *str_next,
                                smm_t_ptr smm,
                                char *end)
{
    char *cur, *next;
    for (cur = str, next = str + 2; *next != '\0' && *(next + 1) != '\0';
         cur++, next++) {
        if (end != NULL && next + 1 == end) {
            return;
        }
        smm->trans_m[indexOf(smm, cur, 2)][indexOf(smm, next, 2)]++;
    }
    char *tmp = next;
    next = (char *) malloc(sizeof(char) * 3);
    memset(next, '\0', 3);
    *(next) = *tmp;
    *(next + 1) = *(str_next);
    smm->trans_m[indexOf(smm, cur, 2)][indexOf(smm, next, 2)]++;
    tmp = str_next + 1;
    smm->trans_m[indexOf(smm, next, 2)][indexOf(smm, tmp, 2)]++;
    *(next) = *(str_next);
    *(next + 1) = *(str_next + 1);
    cur++;
    smm->trans_m[indexOf(smm, cur, 2)][indexOf(smm, next, 2)]++;
    free(next);
}

void simple_fit(smm_t_ptr smm, island_t_ptr island, entry_ptr entry)
{
    switch (smm->order) {
    case 0:
        zero_order_fit(smm, island, entry);
        break;
    case 1:
        first_order_fit(smm, island, entry);
        break;
    case 2:
        second_order_fit(smm, island, entry);
        break;
    }
}



void zero_order_fit(smm_t_ptr smm, island_t_ptr island, entry_ptr entry)
{
    char *str;
    int i;
    for (i = island->from, str = island->start_addr; i < island->to + 1;
         str = entry->content[++i]) {
        calculate_gene(str, smm, island->close_addr);
    }
    // str = entry->content[island->to];
    // calculate_gene(str, smm, island->close_addr);

    double sum = 0;
    for (int i = 0; i < smm->n_states; i++) {
        sum += smm->trans_m[0][i];
    }
    for (int i = 0; i < smm->n_states; i++) {
        smm->trans_m[0][i] /= sum;
    }
}

void first_order_fit(smm_t_ptr smm, island_t_ptr island, entry_ptr entry)
{
    char *str, *prev = NULL;
    int i;
    for (i = island->from, str = island->start_addr, prev = NULL;
         i < island->to + 1; str = entry->content[++i]) {
        calculate_trans_gene(str, smm, island->close_addr, &prev);
    }
    // str = entry->content[island->to];
    // calculate_trans_gene(str, smm, island->close_addr, &prev);

    double sum = 0;
    for (int i = 0; i < smm->n_states; i++) {
        for (int j = 0; j < smm->n_states; j++) {
            sum += smm->trans_m[i][j];
        }
    }
    for (int i = 0; i < smm->n_states; i++) {
        for (int j = 0; j < smm->n_states; j++) {
            smm->trans_m[i][j] /= sum;
        }
    }
}

void second_order_fit(smm_t_ptr smm, island_t_ptr island, entry_ptr entry)
{
    char *str, *prev = NULL, *str_next;
    int i;
    for (i = island->from, str = island->start_addr,
        str_next = entry->content[island->from + 1];
         i < island->to + 1; str = str_next, str_next = entry->content[++i]) {
        calculate_transGene_higher(str, str_next, smm, island->close_addr);
    }
    // str = entry->content[island->to];
    // calculate_transGene_higher(str, str_next, smm, island->close_addr);

    double sum = 0;
    for (int i = 0; i < smm->n_states; i++) {
        for (int j = 0; j < smm->n_states; j++) {
            sum += smm->trans_m[i][j];
        }
    }
    for (int i = 0; i < smm->n_states; i++) {
        for (int j = 0; j < smm->n_states; j++) {
            smm->trans_m[i][j] /= sum;
        }
    }
}



int indexOf(smm_t_ptr smm, char *key, int len)
{
    int value = -1;
    char *tmp = malloc(sizeof(char) * (len + 1));
    memset(tmp, '\0', len + 1);
    strncpy(tmp, key, len);
    for (int i = 0; i < len; i++) {
        tmp[i] = tolower(tmp[i]);
    }
    for (int i = 0; i < smm->n_states; i++) {
        if (strcmp(smm->s_name[i], tmp) == 0) {
            value = i;
            return value;
        }
    }
    free(tmp);
    return value;
}

double cal_sum(double num, double (*f)(double))
{
    double ans = (*f)(num);
    return ans;
}

double zero_order_gen(smm_t_ptr smm, island_t_ptr island, entry_ptr entry)
{
    char *cur = NULL, *tmp;
    int i = 0;
    double sum = 0.0;
    for (cur = island->start_addr, i = island->from; i < island->to + 1;
         cur = entry->content[++i]) {
        for (tmp = cur; *tmp != '\0'; tmp++) {
            if (tmp == island->close_addr) {
                break;
            }
            sum += cal_sum(smm->trans_m[0][indexOf(smm, tmp, 1)], log2);
        }
    }
    return sum;
}


double first_order_gen(smm_t_ptr smm, island_t_ptr island, entry_ptr entry)
{
    char *cur = NULL, *tmp;
    int i = 0;
    double sum = 0;
    for (cur = island->start_addr, i = island->from; i < island->to + 1;) {
        for (tmp = cur; *(tmp + 1) != '\0'; tmp++) {
            if (tmp + 1 == island->close_addr) {
                return sum;
            }
            sum += cal_sum(
                smm->trans_m[indexOf(smm, tmp, 1)][indexOf(smm, tmp + 1, 1)],
                log2);
        }
        cur = entry->content[++i];
        sum += cal_sum(smm->trans_m[indexOf(smm, tmp, 1)][indexOf(smm, cur, 1)],
                       log2);
    }
    return sum;
}

double second_order_gen(smm_t_ptr smm, island_t_ptr island, entry_ptr entry)
{
    char *cur = NULL, *next, *tmp, *row, *str_next;
    int i = 0;
    char *concat = (char *) malloc(sizeof(char) * 3);
    memset(concat, '\0', 3);
    double sum = 0;
    for (row = island->start_addr, i = island->from; i < island->to + 1;
         row = entry->content[++i]) {
        for (cur = row, next = row + 2; *next != '\0' && *(next + 1) != '\0';
             cur++, next++) {
            if (next + 1 == island->close_addr) {
                return sum;
            }
            sum += cal_sum(
                smm->trans_m[indexOf(smm, cur, 2)][indexOf(smm, next, 2)],
                log2);
        }
        *concat = *next;
        str_next = entry->content[i + 1];
        *(concat + 1) = *(str_next);
        tmp = str_next + 1;
        sum += cal_sum(
            smm->trans_m[indexOf(smm, cur, 2)][indexOf(smm, concat, 2)], log2);
        sum += cal_sum(
            smm->trans_m[indexOf(smm, concat, 2)][indexOf(smm, tmp, 2)], log2);
        cur++;
        tmp = str_next;
        sum += cal_sum(smm->trans_m[indexOf(smm, cur, 2)][indexOf(smm, tmp, 2)],
                       log2);
    }
    free(concat);
    return sum;
}

double simple_genS(smm_t_ptr smm, island_t_ptr island, entry_ptr entry)
{
    double prob = -1.0;
    switch (smm->order) {
    case 0:
        prob = zero_order_gen(smm, island, entry);
        break;
    case 1:
        prob = first_order_gen(smm, island, entry);
        break;
    case 2:
        prob = second_order_gen(smm, island, entry);
        break;
    }

    return prob;
}



/* ====================================================== */
/* ================= Hidden Markov Chain ================ */
/* ====================================================== */



hmm_t_ptr HiddenMarkovChain(int states,
                            int e_num,
                            char **e_name,
                            double **trans_m,
                            double **states_m)
{
    return init_hmm(states, e_num, e_name, trans_m, states_m);
}


hmm_t_ptr init_hmm(int states,
                   int e_num,
                   char **e_name,
                   double **trans_m,
                   double **states_m)
{
    hmm_t_ptr hmm = (hmm_t_ptr) malloc(sizeof(hmm_t));
    if (trans_m != NULL) {
        hmm->trans_m = trans_m;
    } else {
        hmm->trans_m = (double **) malloc(sizeof(double *) * states);
        for (int i = 0; i < states; i++)
            hmm->trans_m[i] = (double *) malloc(sizeof(double) * states);
        for (int i = 0; i < states; i++) {
            for (int j = 0; j < states; j++) {
                hmm->trans_m[i][j] = 0.0;
            }
        }
    }
    if (states_m != NULL) {
        hmm->states_m = states_m;
    } else {
        hmm->states_m = (double **) malloc(sizeof(double *) * states);
        for (int i = 0; i < states; i++) {
            hmm->states_m[i] = (double *) malloc(sizeof(double) * e_num);
        }
        for (int i = 0; i < states; i++) {
            for (int j = 0; j < e_num; j++) {
                hmm->states_m[i][j] = 0.0;
            }
        }
    }
    if (e_name == NULL) {
        printf("Not specify the e_name\n");
    } else {
        hmm->e_name = e_name;
    }
    hmm->n_states = states;
    hmm->e_num = e_num;
    return hmm;
}


int indexOf_h(hmm_t_ptr hmm, char *key, int len)
{
    int value = -1;
    char *tmp = malloc(sizeof(char) * (len + 1));
    memset(tmp, '\0', len + 1);
    strncpy(tmp, key, len);
    for (int i = 0; i < len; i++) {
        tmp[i] = tolower(tmp[i]);
    }
    for (int i = 0; i < hmm->e_num; i++) {
        if (strcmp(hmm->e_name[i], tmp) == 0) {
            value = i;
            return value;
        }
    }
    free(tmp);
    return value;
}

double logAdd(double x, double y)
{
    double temp, diff, z;
    if (x < y) {
        temp = x;
        x = y;
        y = temp;
    }
    diff = y - x;          // notice that diff <= 0
    if (diff < minLogExp)  // if y' is far smaller than x'
        return (x < LSMALL) ? LZERO : x;
    else {
        z = exp(diff);
        return x + log2(1.0 + z);
    }
}


double logsum(double *f, int states)
{
    // findmax
    double max = f[0];
    for (int i = 0; i < states; i++) {
        if (max < f[i])
            max = f[i];
    }
    double p = 0;
    for (int i = 0; i < states; i++) {
        f[i] -= max;
        f[i] = exp2(f[i]);
        p += f[i];
    }
    p = log2(p);

    return p + max;
}

void forward_sum(hmm_t_ptr hmm,
                 double *f_next,
                 double *f_cur,
                 char *cur,
                 char *start)
{
    double *f_tmp = (double *) malloc(sizeof(double) * hmm->n_states);
    double p = 0;
    for (int l = 0; l < hmm->n_states; l++) {
        p = 0;
        for (int k = 0; k < hmm->n_states; k++) {
            if (cur == start) {
                f_tmp[k] = log2((double) (1.0f / (double) hmm->n_states));
            } else {
                f_tmp[k] = f_cur[k] + log2(hmm->trans_m[k][l]);
                // p += f_cur[k] + log2(hmm->trans_m[k][l]);
            }
        }
        // * log version
        f_next[l] = logsum(f_tmp, hmm->n_states) +
                    log2(hmm->states_m[l][indexOf_h(hmm, cur, 1)]);
    }
    free(f_tmp);
}



double forward_algorithm(entry_ptr entry,
                         island_t_ptr island,
                         hmm_t_ptr hmm,
                         double **f)
{
    double *f_cur = (double *) malloc(sizeof(double) * hmm->n_states);
    double *f_next = (double *) malloc(sizeof(double) * hmm->n_states);
    double *tmp;
    int *f_iter = (int *) malloc(sizeof(int) * hmm->n_states);
    double s = 0.0f;
    for (int i = 0; i < hmm->n_states; i++) {
        f_cur[i] = 0;
        f_next[i] = 0;
        f_iter[i] = 0;
    }
    int i;
    double P_x = 0;
    char *row, *cur;
    int index = 0;
    for (row = island->start_addr, i = island->from; i < island->to + 1;
         row = entry->content[++i]) {
        for (cur = row; *cur != '\0' && cur != island->close_addr; cur++) {
            // * log version
            forward_sum(hmm, f_next, f_cur, cur, island->start_addr);
            tmp = f_cur;
            f_cur = f_next;
            f_next = tmp;
        }
        if (f != NULL) {
            for (int i = 0; i < hmm->n_states; i++) {
                f[i][index] = f_cur[i];
            }
            index++;
            for (int i = 0; i < hmm->n_states; i++) {
                f[i][index] = 0.0f;
            }
        }
    }
    P_x = logsum(f_cur, hmm->n_states);
    free(f_cur);
    free(f_iter);
    free(f_next);
    return P_x;
}


void reverse_str(char *dst, const char *src, char *end)
{
    char *tmp, *cur;
    memset(dst, '\0', strlen(dst));
    cur = dst;
    int len = 0;
    tmp = src;
    while (*tmp != '\0') {
        if (end != NULL && tmp == end) {
            break;
        }
        len++;
        tmp++;
    }
    tmp--;
    int i;
    i = len;
    while (i-- > 0) {
        *cur++ = *tmp--;
    }
    *cur = '\0';
}

void backward_sum(hmm_t_ptr hmm,
                  double *b_next,
                  double *b_cur,
                  double *b_tmp,
                  char *buffer,
                  int i,
                  double **b,
                  int *index)
{
    char *cur;
    double *tmp;
    for (cur = buffer; *(cur + i) != '\0'; cur++) {
        for (int k = 0; k < hmm->n_states; k++) {
            for (int l = 0; l < hmm->n_states; l++) {
                b_tmp[l] = log2(hmm->trans_m[k][l]) +
                           log2(hmm->states_m[l][indexOf_h(hmm, cur, 1)]) +
                           b_cur[l];
                // p += hmm->trans_m[k][l]*hmm->states_m[l][indexOf_h(hmm, cur,
                // 1)]*b_cur[l];
            }
            b_next[k] = logsum(b_tmp, hmm->n_states);
            if (b != NULL)
                b[k][*index] = b_next[k];
        }
        (*index)++;
        if (b != NULL)
            for (int k = 0; k < hmm->n_states; k++)
                b[k][*index] = 0.0f;
        tmp = b_cur;
        b_cur = b_next;
        b_next = tmp;
    }
}


double backward_algorithm(entry_ptr entry,
                          island_t_ptr island,
                          hmm_t_ptr hmm,
                          double **b)
{
    double *b_cur = (double *) malloc(sizeof(double) * hmm->n_states);
    double *b_next = (double *) malloc(sizeof(double) * hmm->n_states);
    double *b_tmp = (double *) malloc(sizeof(double) * hmm->n_states);
    double a_0l = 1.0f / hmm->n_states;
    for (int i = 0; i < hmm->n_states; i++) {
        b_cur[i] = 0.0f;
        b_next[i] = 0.0f;
    }
    char *row, *cur;
    char *reverse_buffer = (char *) malloc(sizeof(char) * (MAXLEN + 1));
    memset(reverse_buffer, '\0', MAXLEN + 1);
    int i, index = 0;
    for (row = entry->content[island->to], i = island->to; i > island->from - 1;
         row = entry->content[--i]) {
        if (i == island->from || island->to == island->from) {
            reverse_str(reverse_buffer, island->start_addr, island->close_addr);
            backward_sum(hmm, b_next, b_cur, b_tmp, reverse_buffer, 1, b,
                         &index);
        } else if (i == island->to && island->to != island->from) {
            reverse_str(reverse_buffer, row, island->close_addr);
            backward_sum(hmm, b_next, b_cur, b_tmp, reverse_buffer, 0, b,
                         &index);
        } else {
            reverse_str(reverse_buffer, row, NULL);
            backward_sum(hmm, b_next, b_cur, b_tmp, reverse_buffer, 0, b,
                         &index);
        }
    }
    double P_x = 0;
    for (int l = 0; l < hmm->n_states; l++) {
        // * P(X) = sum_l a_0l * e_l(x_1) * b_l(1)
        b_tmp[l] =
            log2((double) a_0l) +
            log2(hmm->states_m[l][indexOf_h(hmm, island->start_addr, 1)]) +
            b_cur[l];
    }
    P_x = logsum(b_tmp, hmm->n_states);
    free(b_cur);
    free(b_next);
    free(reverse_buffer);
    free(b_tmp);
    return P_x;
}

void learning_algorithm(entry_ptr entry, island_t_ptr island, hmm_t_ptr hmm)
{
    double **f_log, **b_log, **A, **E, ***tmp_a, ***tmp_e;

    int *e_index = (int *) malloc(sizeof(int) * hmm->e_num);
    for (int i = 0; i < hmm->e_num; i++)
        e_index[i] = 0;
    double P_x_log;
    int seq_length = (island->to - island->from + 1) * 80;
    f_log = (double **) malloc(sizeof(double *) * hmm->n_states);
    b_log = (double **) malloc(sizeof(double *) * hmm->n_states);
    A = (double **) malloc(sizeof(double *) * hmm->n_states);
    E = (double **) malloc(sizeof(double *) * hmm->n_states);
    tmp_a = (double ***) malloc(sizeof(double **) * hmm->n_states);
    tmp_e = (double ***) malloc(sizeof(double **) * hmm->n_states);
    for (int i = 0; i < hmm->n_states; i++) {
        f_log[i] = (double *) malloc(sizeof(double) * seq_length);
        b_log[i] = (double *) malloc(sizeof(double) * seq_length);
        tmp_a[i] = (double **) malloc(sizeof(double *) * hmm->n_states);
        tmp_e[i] = (double **) malloc(sizeof(double *) * hmm->e_num);
        for (int j = 0; j < hmm->n_states; j++) {
            tmp_a[i][j] = (double *) malloc(sizeof(double) * seq_length);
        }
        for (int j = 0; j < hmm->e_num; j++) {
            tmp_e[i][j] = (double *) malloc(sizeof(double) * seq_length);
        }
        A[i] = (double *) malloc(sizeof(double) * hmm->n_states);
        E[i] = (double *) malloc(sizeof(double) * hmm->e_num);
    }


    P_x_log = forward_algorithm(entry, island, hmm, f_log);
    backward_algorithm(entry, island, hmm, b_log);

    // * update a, e
    char *row, *cur;
    int i, index = 0;
    for (row = island->start_addr, i = island->from, index = 0;
         i < island->to + 1; row = entry->content[++i]) {
        for (cur = row; *cur != '\0' && cur != island->close_addr; cur++) {
            for (int k = 0; k < hmm->n_states; k++) {
                for (int l = 0; l < hmm->n_states; l++) {
                    tmp_a[k][l][index] =
                        f_log[k][index] + log2(hmm->trans_m[k][l]) +
                        log2(hmm->states_m[l][indexOf_h(hmm, cur, 1)]) +
                        b_log[l][index + 1];
                }
                tmp_e[k][indexOf_h(hmm, cur, 1)]
                     [e_index[indexOf_h(hmm, cur, 1)]++] =
                         f_log[k][index] + b_log[k][index];
            }
            index++;
        }
    }

    for (int k = 0; k < hmm->n_states; k++) {
        for (int l = 0; l < hmm->n_states; l++) {
            A[k][l] =
                (logsum(tmp_a[k][l], index) - P_x_log) / (double) seq_length;
            A[k][l] = exp2(A[k][l]);
        }
        for (int e_name = 0; e_name < hmm->e_num; e_name++) {
            E[k][e_name] =
                (logsum(tmp_e[k][e_name], e_index[e_name]) - P_x_log) /
                (double) seq_length;
            E[k][e_name] = exp2(E[k][e_name]);
        }
    }

    // * update a, e
    for (int k = 0; k < hmm->n_states; k++) {
        double A_sum = 0;
        for (int l = 0; l < hmm->n_states; l++) {
            A_sum += A[k][l];
        }
        for (int l = 0; l < hmm->n_states; l++) {
            hmm->trans_m[k][l] = A[k][l] / A_sum;
        }
    }

    for (int k = 0; k < hmm->n_states; k++) {
        double E_sum = 0;
        for (int e_name = 0; e_name < hmm->e_num; e_name++) {
            E_sum += E[k][e_name];
        }
        for (int e = 0; e < hmm->e_num; e++) {
            hmm->states_m[k][indexOf_h(hmm, hmm->e_name[e], 1)] =
                E[k][indexOf_h(hmm, hmm->e_name[e], 1)] / E_sum;
        }
    }

    free(f_log);
    free(b_log);
    free(A);
    free(E);
    free(tmp_a);
    free(tmp_e);
    free(e_index);
}

double max(double *array, int length, int *argmax)
{
    double max;
    max = array[0];
    *argmax = 0;
    for (int i = 0; i < length; i++) {
        if (max < array[i]) {
            max = array[i];
            *argmax = i;
        }
    }
    return max;
}

em_t_ptr viterbi_algorithm(entry_ptr entry, island_t_ptr island, hmm_t_ptr hmm)
{
    int seq_len = (island->to - island->from + 1) * 80;
    double *v_cur = (double *) malloc(sizeof(double) * hmm->n_states);
    double *v_next = (double *) malloc(sizeof(double) * hmm->n_states);
    double *v_tmp = (double *) malloc(sizeof(double) * hmm->n_states);
    int **ptr = (int **) malloc(sizeof(int *) * seq_len);

    for (int i = 0; i < hmm->n_states; i++) {
        v_cur[i] = v_next[i] = v_tmp[i] = 0;
    }
    for (int i = 0; i < seq_len; i++) {
        ptr[i] = (int *) malloc(sizeof(int) * hmm->n_states);
        for (int j = 0; j < hmm->n_states; j++) {
            ptr[i][j] = -1;
        }
    }

    char *row, *cur;
    double *tmp;
    int i, index;
    for (row = island->start_addr, i = island->from, index = 0;
         i < island->to + 1; row = entry->content[++i]) {
        for (cur = row; *cur != '\0' && cur != island->close_addr; cur++) {
            for (int l = 0; l < hmm->n_states; l++) {
                for (int k = 0; k < hmm->n_states; k++) {
                    v_tmp[k] = v_cur[k] + log2(hmm->trans_m[k][l]);
                }
                v_next[l] = log2(hmm->states_m[l][indexOf_h(hmm, cur, 1)]) +
                            max(v_tmp, hmm->n_states, &ptr[index][l]);
            }
            index++;
            tmp = v_cur;
            v_cur = v_next;
            v_next = tmp;
        }
    }

    double P_x = 0;
    for (int k = 0; k < hmm->n_states; k++) {
        v_cur[k] += log2((1.0f / hmm->n_states));
    }
    int *pi = malloc(sizeof(int) * seq_len);
    memset(pi, 0xf, sizeof(pi));
    em_t_ptr em = (em_t_ptr) malloc(sizeof(em_t));
    em->state_path = pi;
    em->len = index;
    P_x = max(v_cur, hmm->n_states, &ptr[index][0]);
    em->P_x = P_x;
    for (int i = index; i >= 0; i--) {
        if (i == index) {
            pi[i - 1] = ptr[index][0];
        } else {
            pi[i - 1] = ptr[i][pi[i]];
        }
    }

    free(v_cur);
    free(v_next);
    free(ptr);
    return em;
}
