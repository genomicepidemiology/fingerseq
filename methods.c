/* Philip T.L.C. Clausen Jan 2020 plan@dtu.dk */

/*
 * Copyright (c) 2017, Philip Clausen, Technical University of Denmark
 * All rights reserved.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *		http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "filebuff.h"
#include "hashfunc.h"
#include "hashmapch.h"
#include "hashmapsc.h"
#include "listsearch.h"
#include "qseqs.h"
#include "seqparse.h"

long unsigned kmer_init(Qseqs *qseq, int kmersize) {
	
	int i;
	long unsigned kmer;
	unsigned char *seq;
	
	i = kmersize < qseq->len ? kmersize : qseq->len;
	
	kmer = 0;
	seq = qseq->seq - 1;
	while(--i) {
		kmer = (kmer << 2) | *++seq;
	}
	
	return kmer;
}

void HashMapSCindex_func(void *src, long unsigned key) {
	HashMapSC_add((HashMapSC *) src, key);
}

void HashMapCHindex_func(void *src, long unsigned key) {
	HashMapCH_add((HashMapCH *) src, key);
}

void linearIndex_func(void *src, long unsigned key) {
	KmerList_push((KmerList *) src, key);
}

void nullIndex_func(void *src, long unsigned key) {
	return;
}

void * indexMethod(void (*func)(void *, long unsigned), FileBuff *inputfile, Qseqs *header, Qseqs *qseq, char *to2bit, int kmersize) {
	
	int i;
	long unsigned kmer, mask, *n;
	unsigned char *seq;
	void *dest;
	
	/* init */
	if(func == &HashMapSCindex_func) {
		dest = HashMapSC_init(1024);
		n = &(((HashMapSC *) dest)->n);
	} else if(func == &HashMapCHindex_func) {
		dest = HashMapCH_init(1024);
		n = &(((HashMapCH *) dest)->n);
	} else {
		dest = KmerList_init(1024);
		n = &(((KmerList *) dest)->n);
	}
	
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * 8 - (kmersize << 1));
	
	/* iterate sequences */
	while(FileBuffgetFsa(inputfile, header, qseq, to2bit)) {
		/* init kmer */
		kmer = kmer_init(qseq, kmersize);
		seq = qseq->seq + (kmersize - 2);
		i = qseq->len - kmersize + 2;
		
		/* add kmers */
		while(--i) {
			kmer = ((kmer << 2) | *++seq) & mask;
			func(dest, kmer);
		}
	}
	
	if(func == &linearIndex_func) {
		/* count and remove duplicates */
		KmerList_countDuplicates((KmerList *) dest);
	} else if(func == &HashMapSCindex_func) {
		/* count collisions */
		mask = HashMapSC_collision((HashMapSC *) dest);
		fprintf(stderr, "# %lu / %lu (%.2f %%) load-factor.\n", *n, ((HashMapSC *) dest)->size + 1, 100.0 * *n / (((HashMapSC *) dest)->size + 1));
		fprintf(stderr, "# %lu / %lu (%.2f %%) collisions.\n", mask, *n, 100.0 * mask / *n);
	} else if(func == &HashMapCHindex_func) {
		/* count collisions */
		mask = HashMapCH_collision((HashMapCH *) dest);
		fprintf(stderr, "# %lu / %lu (%.2f %%) load-factor.\n", *n, ((HashMapCH *) dest)->size + 1, 100.0 * *n / (((HashMapCH *) dest)->size + 1));
		fprintf(stderr, "# %lu / %lu (%.2f %%) collisions.\n", mask, *n, 100.0 * mask / *n);
	}
	
	fprintf(stderr, "# %lu unique kmers saved.\n", *n);
	
	return dest;
}

long unsigned HashMapSCsearch_func(void *src, long unsigned key) {
	return HashMapSC_get((HashMapSC *) src, key);
}

long unsigned HashMapCHsearch_func(void *src, long unsigned key) {
	return HashMapCH_get((HashMapCH *) src, key);
}

long unsigned linearSearch_func(void *src, long unsigned key) {
	return KmerList_linearSearch((KmerList *) src, key);
}

long unsigned binarySearch_func(void *src, long unsigned key) {
	return KmerList_binarySearch((KmerList *) src, key);
}

long unsigned nullSearch_func(void *src, long unsigned key) {
	return 0;
}

void searchMethod(void *src, long unsigned (*func)(void *, long unsigned), FileBuff *inputfile, Qseqs *header, Qseqs *qseq, char *to2bit, int kmersize) {
	
	int i, rc;
	long unsigned kmer, mask, totLookup, match, totmatch, cooccur;
	unsigned char *seq;
	
	/* init */
	totLookup = 0;
	totmatch = 0;
	cooccur = 0;
	mask = 0;
	mask = (~mask) >> (sizeof(long unsigned) * 8 - (kmersize << 1));
	
	/* iterate sequences */
	while(FileBuffgetFsa(inputfile, header, qseq, to2bit)) {
		/* include rc */
		for(rc = 0; rc != 2; ++rc) {
			if(rc) {
				rcQseqs(qseq);
			}
			
			/* init kmer */
			kmer = kmer_init(qseq, kmersize);
			seq = qseq->seq + (kmersize - 2);
			i = qseq->len - kmersize + 2;
			totLookup += qseq->len - kmersize + 1;
			
			/* add kmers */
			while(--i) {
				kmer = ((kmer << 2) | *++seq) & mask;
				if((match = func(src, kmer))) {
					totmatch += match;
					++cooccur;
				}
			}
		}
	}
	
	fprintf(stderr, "# %lu kmer matches towards all templates.\n", totmatch);
	fprintf(stderr, "# %lu kmers co-occured between the query and the templates.\n", cooccur);
	fprintf(stderr, "# Total number of lookups:\t%lu.\n", totLookup);
}

void runMethod(int method, int function, int kmersize, char *queryfilename, char *templatefilename) {
	
	int seqtype;
	long unsigned (*searchFunc)(void *, long unsigned);;
	char to2bit[256];
	void *Method, (*indexFunc)(void *, long unsigned);
	time_t t0, t1;
	FileBuff *inputfile;
	Qseqs *header, *qseq;
	
	/* init */
	memset(to2bit, 8, 256);
	to2bit['\n'] = 16;
	to2bit['A'] = 0;
	to2bit['C'] = 1;
	to2bit['G'] = 2;
	to2bit['T'] = 3;
	to2bit['N'] = 0;
	to2bit['a'] = 0;
	to2bit['c'] = 1;
	to2bit['g'] = 2;
	to2bit['t'] = 3;
	to2bit['n'] = 0;
	to2bit['R'] = 0;
	to2bit['Y'] = 1;
	to2bit['S'] = 2;
	to2bit['W'] = 3;
	to2bit['K'] = 2;
	to2bit['M'] = 0;
	to2bit['B'] = 1;
	to2bit['D'] = 0;
	to2bit['H'] = 3;
	to2bit['V'] = 2;
	to2bit['X'] = 4;
	to2bit['r'] = 0;
	to2bit['y'] = 1;
	to2bit['s'] = 2;
	to2bit['w'] = 3;
	to2bit['k'] = 2;
	to2bit['m'] = 0;
	to2bit['b'] = 1;
	to2bit['d'] = 0;
	to2bit['h'] = 3;
	to2bit['v'] = 2;
	to2bit['x'] = 4;
	to2bit['U'] = 3;
	to2bit['u'] = 3;
	Method = 0;
	if(function == 'k') {
		kmer_hash = &kmer_ashash;
	} else if(function == 's') {
		kmer_hash = &kmer_minimalstandard;
	} else if(function == 'm') {
		kmer_hash = &kmer_murmur;
	} else if(function == 'j') {
		kmer_hash = &kmer_jenkinsoneatatime;
	} else if(function == 'f') {
		kmer_hash = &kmer_FNV1a;
	}
	
	/* allocate */
	header = setQseqs(256);
	qseq = setQseqs(1024);
	inputfile = setFileBuff(1024 * 1024);
	
	/* get template indexes */
	t0 = clock();
	seqtype = openAndDetermine(inputfile, templatefilename);
	if(seqtype & 2) {
		if(method == 'h') {
			indexFunc = &HashMapSCindex_func;
		} else if(method == 'c') {
			indexFunc = &HashMapCHindex_func;
		} else if(method == 'l' || method == 'b') {
			indexFunc = &linearIndex_func;
		} else {
			indexFunc = &nullIndex_func;
		}
		Method = indexMethod(indexFunc, inputfile, header, qseq, to2bit, kmersize);
	} else {
		fprintf(stderr, "File: \"%s\" is fasta format.\n", templatefilename);
		exit(1);
	}
	if(seqtype & 4) {
		gzcloseFileBuff(inputfile);
	} else {
		closeFileBuff(inputfile);
	}
	t1 = clock();
	fprintf(stderr, "#\n# Time to index template with method \"%c\":\t%.2f s.\n", method, difftime(t1, t0) / 1000000);
	
	/* search query */
	t0 = clock();
	seqtype = openAndDetermine(inputfile, queryfilename);
	if(seqtype & 2) {
		if(method == 'h') {
			searchFunc = &HashMapSCsearch_func;
		} else if(method == 'c') {
			searchFunc = &HashMapCHsearch_func;
		} else if(method == 'l') {
			searchFunc = &linearSearch_func;
		} else if(method == 'b') {
			searchFunc = &binarySearch_func;
		} else {
			searchFunc = &nullSearch_func;
		}
		searchMethod(Method, searchFunc, inputfile, header, qseq, to2bit, kmersize);
	} else {
		fprintf(stderr, "File: \"%s\" is fasta format.\n", templatefilename);
		exit(1);
	}
	if(seqtype & 4) {
		gzcloseFileBuff(inputfile);
	} else {
		closeFileBuff(inputfile);
	}
	t1 = clock();
	fprintf(stderr, "#\n# Time to search query against template with method \"%c\":\t%.2f s.\n", method, difftime(t1, t0) / 1000000);
	
	/* clean up */
	destroyQseqs(header);
	destroyQseqs(qseq);
	destroyFileBuff(inputfile);
	if(method == 'h') {
		HashMapSC_destroy((HashMapSC *) Method);
	} else if(method == 'c') {
		HashMapCH_destroy((HashMapCH *) Method);
	} else {
		KmerList_destroy((KmerList *) Method);
	}
}
