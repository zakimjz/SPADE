CC  = g++ -DSGI
#CC = CC -DSGI -no_auto_include
CFLAGS  = -O3
HEADER  = Array.h Itemset.h Lists.h Eqclass.h extl2.h
OBJS	= Itemset.o Array.o Eqclass.o Lists.o extl2.o partition.o
LIBS = -lm -lc
TARGET  = seq calcl2

default:	$(TARGET)

clean:
	rm -rf *~ *.o $(TARGET)

seq: sequence.cc $(OBJS) $(HEADER) 
	$(CC) $(CFLAGS) -o seq sequence.cc $(OBJS) $(LIBS)

Database.o: Database.cc Database.h
	$(CC) $(CFLAGS) -c -o Database.o Database.cc

Lists.o: Lists.cc Lists.h
	$(CC) $(CFLAGS) -c -o Lists.o Lists.cc

Itemset.o: Itemset.cc Itemset.h
	$(CC) $(CFLAGS) -c -o Itemset.o Itemset.cc

Array.o: Array.cc Array.h
	$(CC) $(CFLAGS) -c -o Array.o Array.cc

Eqclass.o: Eqclass.cc Eqclass.h
	$(CC) $(CFLAGS) -c -o Eqclass.o Eqclass.cc

HashTable.o: HashTable.cc HashTable.h
	$(CC) $(CFLAGS) -c -o HashTable.o HashTable.cc

ext.o: ext.cc ext.h
	$(CC) $(CFLAGS) -c -o ext.o ext.cc

extl2.o: extl2.cc extl2.h
	$(CC) $(CFLAGS) -c -o extl2.o extl2.cc

partition.o: partition.cc partition.h
	$(CC) $(CFLAGS) -c -o partition.o partition.cc

calcl2: partition.o calcdb.o calcl2.cc partition.h calcl2.h
	$(CC) $(CFLAGS) -o calcl2 partition.o calcdb.o calcl2.cc $(LIBS)

calcdb.o: calcdb.cc calcdb.h
	$(CC) $(CFLAGS) -c -o calcdb.o calcdb.cc
