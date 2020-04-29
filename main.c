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
#include "fingerseqs.h"
#include "version.h"
#define missArg(opt) fprintf(stderr, "Missing argument at %s.\n", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);

static int helpMessage(FILE *out) {
	
	fprintf(out, "# fingerseq gives standard information from a set of sequence samples. Fasta and fastq format are accepted, gzip compression is allowed.\n");
	fprintf(out, "# %16s\t%-32s\n", "Options are:", "Desc:");
	fprintf(out, "# %16s\t%-32s\n", "-i", "Input file(s)");
	fprintf(out, "# %16s\t%-32s\n", "-f", "Flag <int>");
	fprintf(out, "# %16s\t%-32s\n", "-fh", "Help on flags");
	fprintf(out, "# %16s\t%-32s\n", "-v", "Version");
	fprintf(out, "# %16s\t%-32s\n", "-h", "Shows this helpmessage");
	
	return out == stderr;
}

int main(int argc, char *argv[]) {
	
	int args, flag, filenum;
	char **arg, **filenames, *errorMsg;
	
	/* set defaults */
	flag = 0;
	filenum = 0;
	filenames = 0;
	
	/* parse cmd-line */
	arg = argv;
	args = 0;
	while(++args < argc) {
		++arg;
		if(*(*arg)++ == '-') {
			if(strcmp(*arg, "i") == 0) {
				/* filenames */
				if(++args < argc) {
					filenames = ++arg;
					filenum = 0;
					while(args < argc && **arg != '-') {
						++filenum;
						++args;
						++arg;
					}
					--args;
					--arg;
				}
			} else if(strcmp(*arg, "f") == 0) {
				if(++args < argc) {
					flag = strtoul(*++arg, &errorMsg, 10);
					if(*errorMsg != 0) {
						invaArg("\"-f\"");
					}
				} else {
					missArg("\"-f\"");
				}
			} else if(strcmp(*arg, "fh") == 0) {
				fprintf(stdout, "# Format flags output, add them to combine them.\n");
				fprintf(stdout, "#\n");
				fprintf(stdout, "#   1:\tPredict fasta as reads\n");
				fprintf(stdout, "#\n");
				return 0;
			} else if(strcmp(*arg, "h") == 0) {
				/* help */
				return helpMessage(stdout);
			} else if(strcmp(*arg, "v") == 0) {
				fprintf(stdout, "fingerseq-%s\n", FINGERSEQ_VERSION);
			} else {
				fprintf(stderr, "Unknown option:\t%s\n", --*arg);
				return helpMessage(stderr);
			}
		} else {
			fprintf(stderr, "Unknown argument:\t%s\n", --*arg);
			return helpMessage(stderr);
		}
	}
	
	/* check input */
	if(filenum == 0) {
		fprintf(stderr, "Missing input.\n");
		return helpMessage(stderr);
	}
	
	/* finger file(s) */
	return fingerSeqs(filenames, filenum, flag);
}
