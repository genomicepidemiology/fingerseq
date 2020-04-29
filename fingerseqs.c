/* Philip T.L.C. Clausen Mar 2020 plan@dtu.dk */

/*
 * Copyright (c) 2020, Philip Clausen, Technical University of Denmark
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
#define _XOPEN_SOURCE 600
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "filebuff.h"
#include "fingerseqs.h"
#include "pherror.h"
#include "seqparse.h"
#define pairCheck(seqInfo) (seqInfo->pair == 0 && (seqInfo->tech == Illumina || seqInfo->tech == IonTorrent || seqInfo->tech == fastQ))
#define skipLine(src, len)\
	--src;\
	len = 1;\
	while(*++src && *src != '\n') {\
		++len;\
	}\
	if(*src) {\
		++src;\
	}
#define skipNchar(src, N)\
	++N;\
	while(--N && *++src);\

SeqInfo * SeqInfo_init() {
	
	SeqInfo *dest;
	
	if(!(dest = calloc(1, sizeof(SeqInfo)))) {
		ERROR();
	}
	
	return dest;
}

long unsigned minmaxFileBuff(FileBuff *src) {
	
	unsigned seek, len;
	long unsigned max, maxQ;
	unsigned char *buff;
	
	max = 0;
	maxQ = 0;
	buff = src->next;
	while(*buff != 0) {
		seek = 3;
		while(seek && *buff) {
			if(*++buff == '\n') {
				--seek;
			}
		}
		
		len = 0;
		seek = 1;
		while(seek && *buff) {
			if(*++buff == '\n') {
				seek = 0;
			} else {
				++len;
				if(maxQ < *buff) {
					maxQ = *buff;
				}
			}
		}
		
		if(max < len) {
			max = len;
		}
	}
	
	maxQ -= 33;
	
	/* return min in upper 32 bits, and max in the lower 32 bits */
	return (maxQ << 32) | max;
}

int maxFileBuff(FileBuff *src) {
	
	unsigned len, max;
	unsigned char *buff;
	
	max = 0;
	len = 0;
	buff = src->next - 1;
	while(*++buff != 0) {
		if(*buff == '>') {
			if(max < len) {
				max = len;
			}
			len = 0;
		} else if(*buff != '\n') {
			++len;
		}
	}
	if(max < len) {
		max = len;
	}
	
	return max;
}

int matchHead(unsigned char *src1, unsigned char *src2) {
	
	int diff, len;
	
	diff = -1;
	len = 0;
	--src1;
	--src2;
	while(*++src1 && *++src2 && *src1 != '\n' && *src2 != '\n') {
		if(*src1 != *src2) {
			++diff;
		}
		++len;
	}
	if(*src1 && *src2) {
		++len;
	} else {
		diff = 0;
	}
	if(diff <= 4) {
		/* sra reads */
		diff = 0;
	}
	
	return diff ? 0 : len;
}

int isPair(FileBuff *filebuff, FileBuff *filebuff_rc) {
	
	int i, len, fasta;
	unsigned char *buff, *buff_rc;
	
	/* cmp first headers */
	buff = filebuff->next;
	buff_rc = filebuff_rc->next;
	fasta = *buff == '>';
	while(*buff != 0 && *buff_rc != 0) {
		/* cmp headers */
		if((i = matchHead(buff, buff_rc))) {
			buff += i;
			buff_rc += i;
		} else {
			return 0;
		}
		
		/* seek to next header */
		if(fasta) {
			skipLine(buff, len);
			skipLine(buff_rc, len);
		} else {
			skipLine(buff, len);
			skipLine(buff, i);
			skipNchar(buff, len);
			
			skipLine(buff_rc, len);
			skipLine(buff_rc, i);
			skipNchar(buff_rc, len);
		}
	}
	
	return 1;
}

int fingerSeqs(char **filenames, int filenum, int flag) {
	
	char *fastA = "fastA", *fastQ = "fastQ", *Illumina = "Illumina", *IonTorrent = "Ion Torrent", *Nanopore = "Nanopore", *PacBio = "PacBio", *Na = "Na";
	int i, j, seqtype, seqlen;
	long unsigned max, maxQ;
	char to2bit[256];
	FileBuff **seqbuffs, *seqbuff;
	SeqInfo **seqinfos, *seqinfo, *seqinfo_rc;
	
	/*
	for each file
		check phred
		get tech
	
	for each illumina and fasta
		check if paired
	*/
	
	/* rules:
	Phred-64	->	Illumina
	Phred-0		->	PacBio
	*/
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
	
	/* allocate */
	seqinfos = smalloc(filenum * sizeof(SeqInfo*));
	seqbuffs = smalloc(filenum * sizeof(FileBuff*));
	
	/* get phred scale and technology */
	for(i = 0; i < filenum; ++i) {
		seqinfos[i] = (seqinfo = SeqInfo_init());
		seqbuffs[i] = (seqbuff = setFileBuff(1024 * 1024));
		seqtype = openAndDetermine(seqbuff, filenames[i]);
		if(seqtype & 2) {
			seqinfo->phred = -1;
			seqlen = (flag & 1) ? maxFileBuff(seqbuff) : 0;
			if(seqlen <= 0) {
				seqinfo->tech = fastA;
			} else if(seqlen <= 251) {
				seqinfo->tech = Illumina;
			} else if(seqlen <= 500) {
				seqinfo->tech = IonTorrent;
			} else if(seqlen <= 100000) {
				seqinfo->tech = Nanopore;
			} else {
				seqinfo->tech = fastA;
			}
		} else {
			seqinfo->phred = getPhredFileBuff(seqbuff);
			if(seqinfo->phred == 0) {
				seqinfo->tech = PacBio;
			} else if(seqinfo->phred == 64) {
				seqinfo->tech = Illumina;
			} else {
				seqinfo->tech = fastQ;
				maxQ = minmaxFileBuff(seqbuff);
				max = maxQ & 4294967295U;
				maxQ >>= 32;
				
				if(41 < maxQ) {
					if(max <= 500) {
						seqinfo->tech = IonTorrent;
					} else {
						seqinfo->tech = Nanopore;
					}
				} else if(max <= 251) {
					seqinfo->tech = Illumina;
				} else if(max <= 500) {
					seqinfo->tech = IonTorrent;
				} else {
					seqinfo->tech = Nanopore;
				}
			}
		}
		
		/* close file */
		closeFileBuff(seqbuff);
		if(seqtype & 4) {
			seqbuff->strm->avail_out = 0;
		}
	}
	
	/* check pairing */
	for(i = 0; i < filenum; ++i) {
		seqinfo = seqinfos[i];
		seqbuff = seqbuffs[i];
		
		/* check if pairing is possible */
		if(pairCheck(seqinfo)) {
			j = i;
		} else {
			j = filenum;
		}
		while(++j < filenum && !seqinfo->pair) {
			seqinfo_rc = seqinfos[j];
			if(pairCheck(seqinfo_rc)) {
				if(isPair(seqbuff, seqbuffs[j])) {
					seqinfo->tech = Illumina;
					seqinfo->pair = filenames[j];
					seqinfo_rc->tech = Illumina;
					seqinfo_rc->pair = filenames[i];
				}
			}
		}
		if(!seqinfo->pair) {
			seqinfo->pair = Na;
		}
	}
	
	/* print results */
	fprintf(stdout, "#Sample\tTechnology\tPhred-Scale\tPair\n");
	for(i = 0; i < filenum; ++i) {
		seqinfo = seqinfos[i];
		fprintf(stdout, "%s\t%s\t%d\t%s\n", filenames[i], seqinfo->tech, seqinfo->phred, seqinfo->pair);
	}
	
	return 0;
}
