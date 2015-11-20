#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Maxfiles.h"
#include "MaxSLiCInterface.h"
#include "helpers.h"

static bitarray *P;

struct timespec start_t1,start_t2,end_t1,end_t2;
double  total_t1,total_t2;

double get_elapsed_time(struct timespec *before, struct timespec *after){
  double deltat_s  = after->tv_sec - before->tv_sec;
  double deltat_ns = after->tv_nsec - before->tv_nsec;
  return deltat_s + deltat_ns*1e-9;
}

void cyk_parse_init(){
    int i, j;
    for ( i=0; i < seq_len; i++ ) { 
        CLEAR(P[0*seq_len+i], MAX_NTERM+32);
    	for ( j=0; j < tprod_size; j++ ) {
    		if ( seq[i] == tprod[j].rhs1 )
                OR(P[0*seq_len+i], uclosure[tprod[j].lhs], MAX_NTERM+32);
        }
    }
}

int cyk_parse_CPU() {
    int i, j, k, p;

    cyk_parse_init();

    for ( i=1; i < seq_len; i++ ) {
    	for ( j=0; j < seq_len-i+1; j ++ ) {
    		for ( k=1; k <= i; k++ ) {
                for (p=0; p < num_nterm_pairs; p++ ) {
                    if (   TEST(P[(k-1)*seq_len+j], nprod_cpu[p].rhs1, MAX_NTERM)
                         && TEST(P[(i-k)*seq_len+j+k], nprod_cpu[p].rhs2, MAX_NTERM) )
                    			OR(P[i*seq_len+j], nprod_cpu[p].lhs, MAX_NTERM);
                }
    		}
        }
    }
    // result check
    return getbit(P[(seq_len-1)*seq_len+0], start_neterm);
}

bitarray *input_seq;

void cyk_parse_init_DFE(){
    // initialization
    int i,j;
    input_seq = (bitarray*) malloc( MAX_SEQ_LEN * sizeof(bitarray));
	//for ( i=0; i < MAX_SEQ_LEN; i++ ) {
	//	for ( j=0; j < ba_size; j++ )
    //		input_seq[ i * ba_size + j ] = tprod_dfe[seq[i]][j];
    //	}
    for ( i=0; i < seq_len; i++ ) {
        CLEAR(input_seq[i], MAX_NTERM+32);
    	for ( j=0; j < tprod_size; j++ ) {
    		if ( seq[i] == tprod[j].rhs1 )
                OR(input_seq[i], uclosure[tprod[j].lhs], MAX_NTERM+32);
        }
    }

}

CykParserV2_actions_t run_s;

int cyk_parse_DFE(max_engine_t *engine) {

    int i;


    uint32_t dfe_result[ba_size];

	cyk_parse_init_DFE();

    int g_max = (num_nterm_pairs + (G_PARALLEL*g_lut) - 1) / (G_PARALLEL*g_lut);

    //  total number of iterations of cyk algorithm:
//  n*1 + (n-1)*2 + ... + (n+1-i) * i +...+ (n+1 - (n-1)) * (n-1)
//  that is
//  sum of (n+1-i)*i = sum of (n+1)*i - i^2 = (n+1) * (sum of i) - sum of i^2 = (n+1)*(n-1)*n/2 - (n-1)*n*(n-1/2)/3
//  where i goes from 1 to n-1
//  multiply this by g_max * delay to get total number of ticks
   	int64_t t1 = (seq_len+1)*seq_len*(seq_len-1)/2;
   	int64_t t2 = (seq_len-1)*seq_len*seq_len-(seq_len-1)*seq_len/2;
    run_s.param_N = MAX_SEQ_LEN + (t1  - t2/3 ) * (g_max + 24);

   	run_s.param_seq_len = seq_len;
   	run_s.outstream_dummy = dfe_result;
    run_s.instream_input_seq = (uint32_t*)input_seq;

	clock_gettime(CLOCK_MONOTONIC, &start_t2);
	CykParserV2_run(engine, &run_s);
	clock_gettime(CLOCK_MONOTONIC, &end_t2);

	for (i=0; i<ba_size;i++) {
		//printf("DFE=%8x CPU=%8x\n", dfe_result[i], ((uint32_t *)P)[((seq_len-1)*seq_len+0)*ba_size+i]);
		if (dfe_result[i] != ((uint32_t *)P)[((seq_len-1)*seq_len+0)*ba_size+i])
			printf("Results differ!\n" );
	}
	free(input_seq);
	return 0;
}

void print_Pcpu(){
	int i,j,k;
	for (i=0; i < seq_len; i++) {
		for (j=0; j < seq_len; j++) {
			printf("%d %d P ", i, j);
			for (k=0; k < ba_size; k++) {
				uint32_t Pcpu = ((uint32_t *)P)[(i*seq_len+j)*ba_size + k];
				printf("%8x ", Pcpu);
			}
			printf("\n");
		}
	}
}


int main(int argc, char * argv[]){

    if (argc < 3) {
       printf("Usage: %s <grammar file> <input seq file> [<dfe mem file>]\n", argv[0]);
       return 1;
    }

    // procitati gramatiku => nprod, tprod
    if ( read_grammar(argv[1]) )
        return 1;

    warshall();

    if ( second_pass() ) // make opt warshall closure, data structures for DFE
        return 1;

    // procitati ulazni niz => seq
    if ( read_input_seq(argv[2]) )
        return 1;

    // write dfe memory contents for a grammar to the third file, if given
    if ( argc == 4 )
    	return dfe_mem_write(argv[3]);

    // cyk algoritam
    P = (bitarray*) calloc( seq_len * seq_len, sizeof(bitarray));

    clock_gettime(CLOCK_MONOTONIC, &start_t1);
    int result = cyk_parse_CPU();
    clock_gettime(CLOCK_MONOTONIC, &end_t1);

    if ( result )
    	printf("The sequence belongs to the language\n");
    else
    	printf("The sequence does not belong to the language\n");



 	printf("Running on DFE.\n");

    max_file_t *maxfile = CykParserV2_init();
    max_engine_t *engine = max_load(maxfile, "*");

	cyk_parse_DFE(engine);

	max_unload(engine);

	printf("Done.\n");

    total_t1 = get_elapsed_time(&start_t1, &end_t1);
    total_t2 = get_elapsed_time(&start_t2, &end_t2);
    printf("Total time CPU: %e, total time DFE: %e\n",total_t1,total_t2);

    //print_Pcpu();

    free(P);
    return 0;
}

