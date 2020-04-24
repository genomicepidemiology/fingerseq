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
#include "filebuff.h"

#ifndef FINGERSEQS
typedef struct seqInfo SeqInfo;
struct seqInfo {
	int phred;
	char *tech;
	char *pair;
};
#define FINGERSEQS 1
#endif

SeqInfo * SeqInfo_init();
long unsigned minmaxFileBuff(FileBuff *src);
int maxFileBuff(FileBuff *src);
int matchHead(unsigned char *src1, unsigned char *src2);
int isPair(FileBuff *filebuff, FileBuff *filebuff_rc);
int fingerSeqs(char **filenames, int filenum);
