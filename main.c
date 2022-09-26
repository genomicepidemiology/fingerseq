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
#include "cmdline.h"
#include "fingerseqs.h"
#include "version.h"
#define missArg(opt) fprintf(stderr, "Missing argument at %s.\n", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

static int helpMessage(FILE *out) {
	
	fprintf(out, "#fingerseq gives standard information from a set of sequence samples. Fasta and fastq format are accepted, gzip compression is allowed.\n");
	fprintf(out, "#   %-24s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'i', "input", "Input file(s)", "stdin");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'f', "flag", "Output flags", "1");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'F', "flag_help", "Help on option \"-f\"", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'b', "buffer_size", "File buffer size", "4194304");
	fprintf(out, "#    -%c, --%-17s\t%-32s\t%s\n", 'v', "version", "Version", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'h', "help", "Shows this helpmessage", "");
	
	return out == stderr;
}

int main(int argc, char *argv[]) {
	
	int args, len, offset, flag, filenum, buffSize;
	char **Arg, *arg, **filenames, opt;
	
	/* set defaults */
	flag = 0;
	filenum = 0;
	filenames = 0;
	buffSize = 4 * 1024 * 1024;
	
	/* parse cmd-line */
	/* parse cmd-line */
	args = argc - 1;
	Arg = argv;
	if(args && **++Arg == '-') {
		len = 1;
		--Arg;
	} else {
		len = 0;
	}
	while(args && len) {
		arg = *++Arg;
		if(*arg++ == '-') {
			if(*arg == '-') {
				/* check if argument is included */
				len = getOptArg(++arg);
				offset = 2 + (arg[len] ? 1 : 0);
				
				/* long option */
				if(*arg == 0) {
					/* terminate cmd-line */
					++Arg;
				} else if(strncmp(arg, "input", len) == 0) {
					filenames = getArgListDie(&Arg, &args, len + offset, "input");
					filenum = getArgListLen(&Arg, &args);
				} else if(strncmp(arg, "flag", len) == 0) {
					flag = getNumArg(&Arg, &args, len + offset, "flag");
				} else if(strncmp(arg, "flag_help", len) == 0) {
					flag = -1;
				} else if(strncmp(arg, "buffer_size", len) == 0) {
					buffSize = getNumArg(&Arg, &args, len + offset, "buffer_size");
				} else if(strncmp(arg, "version", len) == 0) {
					fprintf(stdout, "fingerseq-%s\n", FINGERSEQ_VERSION);
				} else if(strncmp(arg, "help", len) == 0) {
					return helpMessage(stdout);
				} else {
					unknArg(arg - 2);
				}
			} else {
				/* multiple option */
				len = 1;
				opt = *arg;
				while(opt && (opt = *arg++)) {
					++len;
					if(opt == 'i') {
						filenames = getArgListDie(&Arg, &args, len, "i");
						filenum = getArgListLen(&Arg, &args);
						opt = 0;
					} else if(opt == 'f') {
						flag = getNumArg(&Arg, &args, len, "f");
						opt = 0;
					} else if(opt == 'F') {
						flag = -1;
					} else if(opt == 'b') {
						buffSize = getNumArg(&Arg, &args, len, "b");
						opt = 0;
					} else if(opt == 'v') {
						fprintf(stdout, "fingerseq-%s\n", FINGERSEQ_VERSION);
					} else if(opt == 'h') {
						return helpMessage(stdout);
					} else {
						*arg = 0;
						unknArg(arg - 1);
					}
				}
			}
		} else {
			/* terminate cmd-line */
			--arg;
			++args;
			len = 0;
		}
		--args;
	}
	
	/* non-options */
	if(args) {
		filenames = Arg;
		filenum = args;
	}
	
	/* flag help */
	if(flag == -1) {
		fprintf(stdout, "# Format flags output, add them to combine them.\n");
		fprintf(stdout, "#\n");
		fprintf(stdout, "#   1:\tPredict fasta as reads\n");
		fprintf(stdout, "#\n");
		return 0;
	}
	
	/* check input */
	if(filenum == 0) {
		fprintf(stderr, "Missing input.\n");
		return helpMessage(stderr);
	}
	
	/* finger file(s) */
	return fingerSeqs(filenames, filenum, flag, buffSize);
}
