#ifdef __cplusplus
extern "C" {
#endif
#include "bitarray.h"

#define MAX_SEQ_LEN    16

//#define MAX_TERM  8000
//#define MAX_NTERM 3040
//#define MAX_PROD_SIZE  10000
// for G40 grammar
#define MAX_TERM  4000
#define MAX_NTERM 1248
#define MAX_PROD_SIZE  4200

#define MAX_NTERM_PAIRS MAX_PROD_SIZE
#define MAX_NPROD_SIZE MAX_PROD_SIZE
#define MAX_TPROD_SIZE MAX_PROD_SIZE

// number of gram. productions processed in parallel by the DFE
#define G_PARALLEL 4
// address size of grammar look up tables in bits
#define g_lut 6

typedef bitarray_t bitarray[RESERVE_BITS(MAX_NTERM)+1]; // align size of this BA with the size od struct prod_dfe
typedef bitarray_t bitarray_p[RESERVE_BITS(MAX_NTERM)];

extern int seq[MAX_SEQ_LEN];
extern int seq_len;

struct prod { // holds grammar productions in the first pass
   int lhs, rhs1, rhs2;
   //bitarray lhs_ba; // valid after warshall()
};
extern struct prod tprod[MAX_TPROD_SIZE];
extern struct prod nprod[MAX_NPROD_SIZE];
extern int tprod_size;
extern int nprod_size;

extern int ba_size;

extern bitarray uclosure[MAX_NTERM]; // for transitive closure of unary prods X -> Y and Z -> X

struct prod_cpu { // holds grammar productions in the second pass, for cpu code
   bitarray lhs, rhs1, rhs2;
};

extern struct prod_cpu nprod_cpu[MAX_NTERM_PAIRS]; // nprod_cpu[i] holds describes all productions of the form X -> Y Z
                                                   // for different X, but with the same right-hand side Y Z

struct prod_dfe { // holds grammar productions in fmem of dfe
   bitarray_p lhs;
   uint16_t rhs1;
   uint16_t rhs2;
};

extern int start_neterm;

extern struct prod_dfe nprod_dfe[MAX_NTERM_PAIRS]; // nprod_dfe[i] holds describes all productions of the form X -> Y Z
                                                   // for different X, but with the same right-hand side Y Z

extern int num_nterm_pairs; // number of nonterminal pairs in grammar

int read_grammar(char* filename);

void warshall();

int second_pass();

int read_input_seq(char* filename);

int dfe_mem_write(char* filename);

#ifdef __cplusplus
}
#endif
