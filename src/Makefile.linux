BIN    = .
CC     = gcc
CPLP   = -fstrict-aliasing -ffast-math -msse2 -g
#-----------------------------------------------------------------------------
CFLAGS = -O3 -Wall $(CPLP)
#-----------------------------------------------------------------------------
LIBS   = -lm -lpthread
DEPS   = defs.h param.h
PROGS  = $(BIN)/smash-map $(BIN)/smash-reduce $(BIN)/smash-visual
OBJS   = mem.o common.o hash.o rmodel.o msg.o parser.o pos.o paint.o time.o seq.o buffer.o 
all:
	$(MAKE) progs
progs: $(PROGS)
$(BIN)/smash-map: smash-map.c $(DEPS) $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/smash-map smash-map.c $(OBJS) $(LIBS)
$(BIN)/smash-reduce: smash-reduce.c $(DEPS) $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/smash-reduce smash-reduce.c $(OBJS) $(LIBS)
$(BIN)/smash-visual: smash-visual.c $(DEPS) $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/smash-visual smash-visual.c $(OBJS) $(LIBS)
mem.o: mem.c mem.h $(DEPS)
	$(CC) -c $(CFLAGS) mem.c
common.o: common.c common.h $(DEPS)
	$(CC) -c $(CFLAGS) common.c
hash.o: hash.c hash.h $(DEPS)
	$(CC) -c $(CFLAGS) hash.c
rmodel.o: rmodel.c rmodel.h $(DEPS)
	$(CC) -c $(CFLAGS) rmodel.c
msg.o: msg.c msg.h $(DEPS)
	$(CC) -c $(CFLAGS) msg.c
parser.o: parser.c parser.h $(DEPS)
	$(CC) -c $(CFLAGS) parser.c
paint.o: paint.c paint.h $(DEPS)
	$(CC) -c $(CFLAGS) paint.c
pos.o: pos.c pos.h $(DEPS)
	$(CC) -c $(CFLAGS) pos.c
time.o: time.c time.h $(DEPS)
	$(CC) -c $(CFLAGS) time.c
seq.o: seq.c seq.h $(DEPS)
	$(CC) -c $(CFLAGS) seq.c
buffer.o: buffer.c buffer.h $(DEPS)
	$(CC) -c $(CFLAGS) buffer.c
clean:
	/bin/rm -f *.o
cleanall:
	/bin/rm -f *.o $(PROGS)

