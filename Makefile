CFLAGS ?= -Wall -O3
CFLAGS += -std=c99
LIBS = filebuff.o fingerseqs.o qseqs.o pherror.o seqparse.o
PROGS = fingerseq

.c .o:
	$(CC) $(CFLAGS) -c -o $@ $<

all: $(PROGS)

fingerseq: main.c libfingerseq.a
	$(CC) $(CFLAGS) -o $@ main.c libfingerseq.a -lz $(LDFLAGS)

libfingerseq.a: $(LIBS)
	$(AR) -csru $@ $(LIBS)

clean:
	$(RM) $(LIBS) $(PROGS) libfingerseq.a


filebuff.o: filebuff.h pherror.h qseqs.h
fingerseqs.o: fingerseqs.h filebuff.h pherror.h seqparse.h
qseqs.o: qseqs.h pherror.h
pherror.o: pherror.h
seqparse.o: seqparse.h filebuff.h qseqs.h
