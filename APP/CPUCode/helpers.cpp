#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include "helpers.h"

using namespace std;

static map<string,int> nterm; // maps nonterminals to integer indices
static map<string,int> term; // maps terminals to integer indices

static map<pair<int,int>,int> prod_num; // maps binary productions right hand sizes (nterm pairs) to integers

static int closure_needed;

bitarray uclosure[MAX_NTERM]; // for transitive closure of unary prods X -> Y and Z -> X
// data structures for grammar productions in the first pass, before closure and converting to bitarrays
struct prod tprod[MAX_TPROD_SIZE];
struct prod nprod[MAX_NPROD_SIZE];
int tprod_size;
int nprod_size;
// end first pass data structures

// data structures for a grammar in a final pass, with bitarrays
int ba_size = RESERVE_BITS(MAX_NTERM)+1;

struct prod_cpu nprod_cpu[MAX_NTERM_PAIRS];

//bitarray tprod_dfe[MAX_TERM];
struct prod_dfe nprod_dfe[MAX_NTERM_PAIRS];

int start_neterm;
int num_nterm_pairs;

int seq[MAX_SEQ_LEN];
int seq_len;
// end final pass data structures

int xterm_to_int(map<string,int>& xterm, string nt){
        map<string,int>::iterator it = xterm.find(nt);
        if ( it == xterm.end() ) {
            // nema ga u nterm
            int size = xterm.size() + 1;
            xterm[nt]=size;
            return size;
        }
        else {
            return it->second;
        }

}

int read_grammar(char* filename) {
    ifstream fin(filename);
    string line;
    if ( !fin.is_open() ){
        cout << "Error opening cfg file: " << filename << '\n';
        return 1;
    }

    int err=0;

    tprod_size=0;
    nprod_size=0;

    closure_needed=0;
    for (int i=0; i < MAX_NTERM; i++ )
	CLEAR( uclosure[i], MAX_NTERM);

    while ( getline (fin,line) ) {
        //cout << line << '\n';
        stringstream sin(line);
        string lhs, arrow, rhs1, rhs2;
        sin >> lhs;
        int nlhs = xterm_to_int(nterm,lhs);
        sin >> arrow;
        if (arrow != "->") {
            err=1;
            cout << "Bad format of cfg file.\n";
            break;
        }
        sin >> rhs1;
        if ( rhs1[0]=='\'' ) {
           // terminal production X -> 't'
           int nrhs1 = xterm_to_int(term,rhs1.substr(1, rhs1.size()-2));
           if ( tprod_size >= MAX_TPROD_SIZE ) {
              err=1;
              cout << "Too many terminal productions.\n";
              break;
           }
           tprod[tprod_size].lhs = nlhs;
           tprod[tprod_size++].rhs1 = nrhs1;
	   //cout << nlhs << "->" << nrhs1 << '\n';
        } else {
           sin >> rhs2;
	   if ( rhs2 == "" ) {
               // unary production X -> Y
               closure_needed++;
               int nrhs1 = xterm_to_int(nterm,rhs1);
               putbit(uclosure[nrhs1], nlhs, 1);
               //cout << nlhs << "->" << nrhs1 << '\n';
           } else {
               // binary production X -> Y Z
               int nrhs1 = xterm_to_int(nterm,rhs1);
               int nrhs2 = xterm_to_int(nterm,rhs2);
               if ( nprod_size >= MAX_NPROD_SIZE ) {
                   err=1;
                   cout << "Too many nonterminal productions.\n";
               break;
               }
               nprod[nprod_size].lhs = nlhs;
               nprod[nprod_size].rhs1 = nrhs1;
               nprod[nprod_size++].rhs2 = nrhs2;
               //cout << nlhs << "->" << nrhs1 << ' ' << nrhs2 << '\n';
           }
       }
    } // while
    fin.close();
    if ( nterm.size() > MAX_NTERM-1 ) { // because nterm indices start at 1
        err=1;
        cout << "Too many nonterminals.\n";
    }
    start_neterm = nprod[0].lhs;
    return err;
}

void warshall() {

    int num_nterm = nterm.size();
    for (int i=1; i <= num_nterm; i++ )
        putbit(uclosure[i],i,1);

    if ( !closure_needed )
        return;

     for (int i=1; i <= num_nterm; i++ )
           for (int j=1; j <= num_nterm; j++ )
              for (int k=1; k <= num_nterm; k++ )
                 if ( getbit(uclosure[i], k) && getbit(uclosure[k], j) )
                     putbit(uclosure[i],j,1);
}

int second_pass() {
    int err=0;

    for (int i=0; i < MAX_NTERM_PAIRS; i++ ) {
        CLEAR(nprod_dfe[i].lhs, MAX_NTERM);
        nprod_dfe[i].rhs1=0;
        nprod_dfe[i].rhs2=0;
        CLEAR(nprod_cpu[i].lhs, MAX_NTERM);
        CLEAR(nprod_cpu[i].rhs1, MAX_NTERM);
        CLEAR(nprod_cpu[i].rhs2, MAX_NTERM);
    }

    for (int i=0; i < nprod_size; i++) {
        //putbit(tprod_dfe[tprod[i].rhs1], tprod[i].lhs, 1);
        pair<int,int> npair(nprod[i].rhs1, nprod[i].rhs2);
        map<pair<int,int>,int>::iterator it = prod_num.find(npair);
        int npr;
        if ( it == prod_num.end() ) {
             // nema ga u nterm
             npr = prod_num.size();
             if ( npr >= MAX_NTERM_PAIRS ) {
                err=1;
                printf("Too many different pairs of nonterminals in binary productions.\n");
                return err;
             }
             prod_num[npair]=npr;
             nprod_dfe[npr].rhs1=nprod[i].rhs1;
             nprod_dfe[npr].rhs2=nprod[i].rhs2;
             putbit(nprod_cpu[npr].rhs1, nprod[i].rhs1, 1);
             putbit(nprod_cpu[npr].rhs2, nprod[i].rhs2, 1);
        } else {
             npr = it->second;
        }
        //putbit(nprod_dfe[npr].lhs, nprod[i].lhs, 1);
        OR(nprod_dfe[npr].lhs, uclosure[nprod[i].lhs], 32*(ba_size-1));
        OR(nprod_cpu[npr].lhs, uclosure[nprod[i].lhs], 32*ba_size);
    }
    num_nterm_pairs = prod_num.size();
    cout << "Nonterminals: " << nterm.size() << "\n   Terminals: " << term.size() << '\n';
    cout << "Binary prods: " << nprod_size << "\n Term. prods: " << tprod_size << "\n Unary prods: " << closure_needed << '\n';
    cout << "Unique pairs: " << num_nterm_pairs << '\n';
    return err;
}

int read_input_seq(char* filename) {
    ifstream fin(filename);

    if ( !fin.is_open() ) {
        cout << "Error opening input seq file: " << filename << '\n';
        return 1;
    }

    string word;
    int err=0;
    seq_len = 0;

    while ( fin >> word ){
       if ( seq_len >= MAX_SEQ_LEN ) {
           err=1;
           cout << "Input sequence too long.\n";
           break;
       }
       seq[seq_len++] = xterm_to_int( term, word);
       cout << word << ' ' << seq[seq_len-1] << '\n';
    }

    fin.close();
    return err;
}


int dfe_mem_write(char* filename) {

    int num_tables = (num_nterm_pairs + g_lut - 1) / g_lut;
    while ( num_tables % G_PARALLEL)
    	num_tables++;
    int g_max = num_tables / G_PARALLEL;
    if ( g_max == 1)
    	g_max=2; // because FMEM depth cannot be less than 2
    int g_max_bits = -1;
    while ( 1<<(++g_max_bits) < g_max )
    	;
    const int g_tab_size = (1 << g_lut)*(1<<g_max_bits);
    bitarray* g_tab = new bitarray[g_tab_size];

    ofstream fout(filename);
    if ( num_tables * g_lut > MAX_NTERM_PAIRS) {
    	cout << "Please increase MAX_NTERM_PAIRS to " << (num_tables * g_lut) <<"\n";
    	return 1;
    }
    if ( !fout.is_open() ) {
        cout << "Error opening output dfe mem file: " << filename << '\n';
        return 1;
    }

    // write out rhs1:rhs2 pairs of productions; one line per G_PARALLEL set of look up tables
    int index=0;
    for (int t=0; t < G_PARALLEL; t++) {
        for (int p=0; p < g_lut; p++) {
            for (int i=0; i < g_max; i++) {
            	fout << hex << std::setw(4) << setfill('0') << nprod_dfe[index+(G_PARALLEL * g_lut) * i].rhs2 << hex << std::setw(4) << setfill('0') << nprod_dfe[index+(G_PARALLEL * g_lut) * i].rhs1;
            }
            index++;
            fout << '\n';
        }

    }

    fout << "Proba\n";
    // write out look up tables for productions
    for (int i=0; i < G_PARALLEL; i++) {
    	// prepare table contents in g_tab
    	for (int gg=0; gg < g_max; gg++) {
    		for (int c=0; c< (1 << g_lut); c++) {
    			CLEAR(g_tab[gg * (1<<g_lut) + c],MAX_NTERM);
    			for (int bb=0; bb<g_lut; bb++) {
    				int gc =  i * g_lut + gg * G_PARALLEL * g_lut +bb;
    				if ( c & (1<<bb) )
    					OR(g_tab[gg * (1<<g_lut) + c], nprod_dfe[gc].lhs, MAX_NTERM);
    			}
    		}
    	}
    	// write g_tab contents to the file
        for (int b=0; b<ba_size-1; b++) {
        	for (int c=0; c< g_tab_size; c++) {
        		fout << hex << std::setw(8) << setfill('0') << g_tab[c][b];
        	}
        	fout << "\n";
    	}
    }
    fout.close();
    delete[] g_tab;
	return 0;
}
