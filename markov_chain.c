#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "markov_chain.h"
#include <math.h>
#include <ctype.h>


/* ====================================================== */
/* ================= Simple Markov Chain ================ */
/* ====================================================== */

smm_t_ptr SimpleMarkovChain(int states, char** states_name, int order){

    smm_t* sm = (smm_t*) malloc(sizeof(smm_t));
    init_smm(sm, states, order);
    sm->s_name = states_name;

    return sm;
}

void init_smm(smm_t_ptr smm, int states, int order){

    smm->order = order;
    // init first order
    smm->trans_m = (double**) malloc(sizeof(double*) * states);
    for(int i=0;i<states;i++){
        smm->trans_m[i] = (double*) malloc(sizeof(double) * states);
        for(int j=0;j<states;j++){
            smm->trans_m[i][j] = 0.0;
        }
    }

    smm->n_states = states;
}



void calculate_trans_gene(char* str, smm_t_ptr smm, char* end, char** prev){
    char* cur;

    if( *prev != NULL) {
        smm->trans_m[indexOf(smm, *prev, 1)][indexOf(smm, str, 1)]++;
    }
    for(cur=str;*cur != '\0'; cur++){
        if(end != NULL && cur == end){
            return;
        }
        smm->trans_m[indexOf(smm, cur, 1)][indexOf(smm, cur+1, 1)]++;
    }
    *prev = cur-1;
}

void calculate_gene(char* str, smm_t_ptr smm, char* end){
    char* tmp = str;
    for(tmp=str;*tmp != '\0';tmp++){
        if(end != NULL && tmp == end){
            return;
        }
        smm->trans_m[0][indexOf(smm, tmp, 1)]++;
    }
}

void calculate_transGene_higher(
    char* str,
    char*str_next,
    smm_t_ptr smm,
    char* end){
    char *cur, *next;
    for(cur=str, next=str+2;
        *next != '\0' && *(next+1) != '\0';
        cur++, next++){
        if(end != NULL && next+1 == end){
            return;
        }
        smm->trans_m[indexOf(smm, cur, 2)][indexOf(smm, next, 2)]++;
    }
    char *tmp = next;
    next = (char*) malloc(sizeof(char)*3);
    memset(next, '\0', 3);
    *(next) = *tmp;
    *(next+1) = *(str_next);
    smm->trans_m[indexOf(smm, cur, 2)][indexOf(smm, next, 2)]++;
    tmp = str_next+1;
    smm->trans_m[indexOf(smm, next, 2)][indexOf(smm, tmp, 2)]++;
    *(next) = *(str_next);
    *(next+1) = *(str_next+1);
    cur++;
    smm->trans_m[indexOf(smm, cur, 2)][indexOf(smm, next, 2)]++;
    free(next);
}

void simple_fit(smm_t_ptr smm, island_t_ptr island, entry_ptr entry){
    switch(smm->order){
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



void zero_order_fit(smm_t_ptr smm, island_t_ptr island, entry_ptr entry){
    char *str;
    int i;
    for(i=island->from, str=island->start_addr; i<island->to+1;str=entry->content[++i]){
        calculate_gene(str, smm, island->close_addr);
    }
    // str = entry->content[island->to];
    // calculate_gene(str, smm, island->close_addr);

    double sum = 0;
    for(int i=0;i<smm->n_states;i++){
        sum += smm->trans_m[0][i];
        // printf("%lf\n", smm->trans_m[0][i]);
    }
    for(int i=0;i<smm->n_states;i++){
        smm->trans_m[0][i] /= sum;
        // printf("%lf\n", smm->trans_m[0][i]);
    }
}

void first_order_fit(smm_t_ptr smm, island_t_ptr island, entry_ptr entry){
    char* str, *prev = NULL;
    int i;
    for(i=island->from,str=island->start_addr, prev=NULL;
        i<island->to+1;
        str=entry->content[++i]){
        calculate_trans_gene(str, smm, island->close_addr, &prev);
    }
    // str = entry->content[island->to];
    // calculate_trans_gene(str, smm, island->close_addr, &prev);

    double sum = 0;
    for(int i=0;i<smm->n_states;i++){
        for(int j=0;j<smm->n_states;j++){
            sum += smm->trans_m[i][j];
        }
    }
    for(int i=0;i<smm->n_states;i++){
        for(int j=0;j<smm->n_states;j++){
            smm->trans_m[i][j] /= sum;
        }
    }
}

void second_order_fit(smm_t_ptr smm, island_t_ptr island, entry_ptr entry){
    char *str, *prev=NULL, *str_next;
    int i;
    for(i=island->from, str=island->start_addr, str_next=entry->content[island->from+1];
        i<island->to+1;
        str=str_next,
        str_next=entry->content[++i]){
        calculate_transGene_higher(str, str_next, smm, island->close_addr);
    }
    // str = entry->content[island->to];
    // calculate_transGene_higher(str, str_next, smm, island->close_addr);

    double sum = 0;
    for(int i=0;i<smm->n_states;i++){
        for(int j=0;j<smm->n_states;j++){
            sum += smm->trans_m[i][j];
        }
    }
    for(int i=0;i<smm->n_states;i++){
        for(int j=0;j<smm->n_states;j++){
            smm->trans_m[i][j] /= sum;
        }
    }
}




int indexOf(smm_t_ptr smm, char* key, int len){

    int value = -1;
    char* tmp = malloc(sizeof(char) * (len+1));
    memset(tmp, '\0', len+1);
    strncpy(tmp, key, len);
    for(int i=0;i<len;i++){
        tmp[i] = tolower(tmp[i]);
    }
    for(int i=0;i<smm->n_states;i++){
        if(strcmp(smm->s_name[i], tmp) == 0){
            value = i;
            return value;
        }
    }
    free(tmp);
    return value;
}

double cal_sum(double num, double(*f)(double)){
    double ans = (*f)(num);
    // printf("log(%lf) = %lf\n", num, ans);
    return ans;
}

double zero_order_gen(smm_t_ptr smm, island_t_ptr island, entry_ptr entry){

    char* cur = NULL, *tmp;
    int i=0;
    double sum = 0.0;
    for(cur=island->start_addr, i=island->from; i<island->to+1 ; cur=entry->content[++i]){
        for(tmp = cur; *tmp != '\0'; tmp++){
            if(tmp == island->close_addr){
                break;
            }
            sum += cal_sum(smm->trans_m[0][indexOf(smm, tmp, 1)], log2);
        }
    }
    return sum;
}


double first_order_gen(smm_t_ptr smm, island_t_ptr island, entry_ptr entry){
    char* cur = NULL, *tmp;
    int i=0;
    double sum = 0;
    for(cur=island->start_addr, i=island->from; i<island->to+1;){
        for(tmp = cur; *(tmp+1) != '\0' ; tmp++){
            if(tmp+1 == island->close_addr){
                return sum;
            }
            sum += cal_sum(smm->trans_m[indexOf(smm, tmp, 1)][indexOf(smm, tmp+1, 1)], log2);
        }
        cur=entry->content[++i];
        sum += cal_sum(smm->trans_m[indexOf(smm, tmp, 1)][indexOf(smm, cur, 1)], log2);
    }
    return sum;
}

double second_order_gen(smm_t_ptr smm, island_t_ptr island, entry_ptr entry){

    char *cur = NULL, *next, *tmp, *row, *str_next;
    int i=0;
    char* concat = (char*) malloc(sizeof(char) * 3);
    memset(concat, '\0', 3);
    double sum = 0;
    for(row=island->start_addr, i=island->from; i<island->to+1;row=entry->content[++i]){
        for(cur=row, next=row+2; *next != '\0' && *(next+1) != '\0'; cur++, next++ ){
            if(next+1 == island->close_addr){
                return sum;
            }
            sum += cal_sum(smm->trans_m[indexOf(smm, cur, 2)][indexOf(smm, next, 2)], log2);
        }
        *concat = *next;
        str_next = entry->content[i+1];
        *(concat+1) = *(str_next);
        tmp = str_next+1;
        sum += cal_sum(smm->trans_m[indexOf(smm, cur, 2)][indexOf(smm, concat, 2)], log2);
        sum += cal_sum(smm->trans_m[indexOf(smm, concat, 2)][indexOf(smm, tmp, 2)], log2);
        cur++;
        tmp = str_next;
        sum += cal_sum(smm->trans_m[indexOf(smm, cur, 2)][indexOf(smm, tmp, 2)], log2);
    }
    free(concat);
    return sum;
}

double simple_genS(smm_t_ptr smm, island_t_ptr island, entry_ptr entry){

    double prob = -1.0;
    switch(smm->order){
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



hmm_t_ptr HiddenMarkovChain(
    int states,
    int e_num,
    char** e_name,
    double** trans_m,
    double** states_m)
{
    return init_hmm(states, e_num, e_name, trans_m, states_m);
}


hmm_t_ptr init_hmm(
    int states,
    int e_num,
    char** e_name,
    double** trans_m,
    double** states_m
    ){

    hmm_t_ptr hmm = (hmm_t_ptr) malloc(sizeof(hmm_t));
    if(trans_m != NULL){
        hmm->trans_m = trans_m;
    }else{
        hmm->trans_m = (double**) malloc(sizeof(double*)*states);
        for(int i=0;i<states;i++)
            hmm->trans_m[i] = (double*) malloc(sizeof(double)*states);
        for(int i=0;i<states;i++){
            for(int j=0;j<states;j++){
                hmm->trans_m[i][j] = 0.0;
            }
        }
    }
    if(states_m != NULL){
        hmm->states_m = states_m;
    }else{
        hmm->states_m = (double**) malloc(sizeof(double*) * states);
        for(int i=0;i<states;i++){
            hmm->states_m[i] = (double*) malloc(sizeof(double)* e_num);
        }
        for(int i=0;i<states;i++){
            for(int j=0;j<e_num;j++){
                hmm->states_m[i][j] = 0.0;
            }
        }
    }
    if(e_name == NULL){
        printf("Not specify the e_name\n");
    }else{
        hmm->e_name = e_name;
    }
    hmm->n_states = states;
    hmm->e_num = e_num;
    return hmm;
}


int indexOf_h(hmm_t_ptr hmm, char* key, int len){

    int value = -1;
    char* tmp = malloc(sizeof(char) * (len+1));
    memset(tmp, '\0', len+1);
    strncpy(tmp, key, len);
    for(int i=0;i<len;i++){
        tmp[i] = tolower(tmp[i]);
    }
    for(int i=0;i<hmm->e_num;i++){
        if(strcmp(hmm->e_name[i], tmp) == 0){
            value = i;
            return value;
        }
    }
    free(tmp);
    return value;
}

double logAdd(double x, double y){
    double temp, diff, z;
    if (x < y)
    {
        temp = x; x = y; y = temp;
    }
    diff = y-x; // notice that diff <= 0
    if (diff < minLogExp)   // if y' is far smaller than x'
        return (x < LSMALL) ? LZERO : x;
    else
    {
        z = exp(diff);
        return x + log2(1.0 + z);
    }
}

double forward_sum(hmm_t_ptr hmm, double* f, int state, char* cur, int *index){
    double sum = 0, f_sum = 0;
    double* f_arr = (double*) malloc(sizeof(double) * hmm->n_states);
    for(int i=0;i<hmm->n_states; i++){
        f_arr[i] = f[i] * hmm->trans_m[i][state];
    }
    for(int i=0;i<hmm->n_states;i++){
        // sum += exp(f_arr[i] - max);
        f_sum += f_arr[i];
    }
    sum = hmm->states_m[state][indexOf_h(hmm, cur, 1)] * f_sum;
    (*index)++;
    free(f_arr);
    return sum;
}

double forward_algorithm(entry_ptr entry, island_t_ptr island, hmm_t_ptr hmm){
    double* f_cur = (double*) malloc(sizeof(double) * hmm->n_states);
    double* f_next = (double*) malloc(sizeof(double) * hmm->n_states);
    double *tmp;
    int* f_iter = (int*) malloc(sizeof(int) * hmm->n_states);
    double s = 0.0f;
    for(int i=0;i<hmm->n_states;i++){
        f_cur[i] = 0;
        f_next[i] = 0;
        f_iter[i] = 0;
    }
    int i;
    double P_x = 0;
    char* row, *cur;
    double p;
    for(row=island->start_addr, i=island->from; i<island->to+1;row=entry->content[++i]){
        for(cur=row; *cur != '\0' && cur != island->close_addr; cur++){
            // printf("========================\n");
            for(int l=0;l<hmm->n_states;l++){
                p = 0;
                for(int k=0;k<hmm->n_states;k++){
                    if(cur == island->start_addr){
                        p = 1*(1.0f/hmm->n_states);
                    }else{
                        p += f_cur[k] * hmm->trans_m[k][l];
                    }
                }
                f_next[l] = p * hmm->states_m[l][indexOf_h(hmm, cur, 1)];
            }
            tmp = f_cur;
            f_cur = f_next;
            f_next = tmp;
        }
    }
    for(int i=0;i<hmm->n_states;i++){
        // * P(x) = sum_{k}f_k(L) * a_{k0}
        P_x += f_cur[i] * 1;
    }
    free(f_cur);
    free(f_iter);
    free(f_next);
    return P_x;
}


void reverse_str(char* dst, const char* src, char* end){
    char *tmp, *cur;
    memset(dst, '\0', strlen(dst));
    cur = dst;
    int len = 0;
    tmp = src;
    while( *tmp != '\0' ){
        if(end != NULL && tmp == end){
            break;
        }
        len++;
        tmp++;
    }
    tmp--;
    int i;
    i = len;
    while(i-- > 0){
        *cur++ = *tmp--;
    }
    *cur = '\0';
}

double backward_algorithm(entry_ptr entry, island_t_ptr island, hmm_t_ptr hmm){
    double* b_cur = (double*) malloc(sizeof(double)* hmm->n_states);
    double* b_next = (double*) malloc(sizeof(double)* hmm->n_states);
    double* tmp;
    for(int i=0;i<hmm->n_states;i++){
        b_cur[i] = 1.0f;
        b_next[i] =  1.0f;
    }
    char *row, *cur;
    char *reverse_buffer = (char*) malloc(sizeof(char)*(MAXLEN+1));
    memset(reverse_buffer, '\0', MAXLEN+1);
    int i;
    for(row=entry->content[island->to], i=island->to; i>island->from-1; row=entry->content[--i]){
        if(i == island->to && island->to == island->from )
            reverse_str(reverse_buffer, island->start_addr, island->close_addr);
        else if(i == island->to && island->to != island->from){
            reverse_str(reverse_buffer, row, island->close_addr);
        }
        else
            reverse_str(reverse_buffer, row, NULL);
        for(cur=reverse_buffer; *(cur+1) != '\0'; cur++){
            for(int k=0;k<hmm->n_states;k++){
                double p = 0;
                for(int l=0;l<hmm->n_states;l++){
                    p += hmm->trans_m[k][l]*hmm->states_m[l][indexOf_h(hmm, cur, 1)]*b_cur[l];
                }
                b_next[k] = p;
            }
            tmp = b_cur;
            b_cur = b_next;
            b_next = tmp;
        }
    }
    double P_x = 0;
    double a_0l = 1.0f / hmm->n_states;
    for(int l=0;l<hmm->n_states;l++){
        // * P(X) = sum_l a_0l * e_l(x_1) * b_l(1)
        P_x += a_0l*hmm->states_m[l][indexOf_h(hmm, island->start_addr, 1)]*b_cur[l];
    }
    free(b_cur);
    free(b_next);
    free(reverse_buffer);
    return P_x;
}









