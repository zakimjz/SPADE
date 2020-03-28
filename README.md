# SPADE Algorithm

Spade uses the vertical format for mining the set of all frequent sequences from a dataset of many sequences. It also mines the frequent sequences of itemsets.

**Relevant Publications**

*[2001-spade:mlj] Mohammed J. Zaki. SPADE: an efficient algorithm for mining frequent sequences. Machine Learning Journal, 42(1/2):31–60, Jan/Feb 2001. special issue on Unsupervised Learning. doi:10.1023/A:1007652502315.

*[1998-spade] Mohammed J. Zaki. Efficient enumeration of frequent sequences. In 7th ACM International Conference on Information and Knowledge Management. November 1998.

# HOW TO

You also need to download utilities like getconf, exttpose from https://github.com/zakimjz/tposedb

Follow these steps:

1) Generate a data file using the IBM data generator program,
gen. This is the same as the public version from Almaden, but I have
enhanced it to generate a configuration file, and a modified ascii
database file.

For SPADE make sure the DB is in binary format.  The output file
should be named XXX.data (note the 'data' extension). The gen program
also produces a XXX.conf file automatically.

The format of the binary file should be 
<cid> <tid> <numitem> <item list>

run gen seq -help for data generation options.

2) run: exttpose -i XXX -o XXX -l -s LMINSUP

        e.g. exttpose -i XXX -o XXX -x -l -s 0 

note: this produces the files XXX.tpose, and XXX.idx

The XXX.tpose file is the DB in vertical format, and
XXX.idx is an index file specifying where the tid-list for each item
begins.

You can specify a value of LMINSUP to be the same as the one you will use to 
run spade below, in which case you will have to rerun exttpose each time you 
use a new lower MINSUP. Alternatively, you can use a small value for LMINSUP, 
and it will continue to work for all values of MINSUP >= LMINSUP when you
run spade.

The time for inverting is stored in summary.out. The format is:
TPOSE DB_FILENAME X NUMITEMS TOTAL_TIME

(see note one TOTAL_TIME below)

3) Now you can run SPADE as follows:

        spade -i XXX -s MINSUP -e 1 -r

MINSUP is in fractions, i.e., specify 0.5 if you want 50% minsup or
0.01 if you want 1% support.

note that the summary of the run is stored in the summary.out file. The format of this file is as follows:

SPADE DB_FILENAME MINSUP ACTUAL_SUPPORT NUMBER_OF_JOINS JOIN_TIME_FOR_F2
      XXX (F1 F2 F3 .... ) F1_TIME F2_TIME FK_TIME TOTAL_TIME

The values of most interest to you are F1, F2, F3, ... F_k, ... which
give the number of frequent k-sequences and TOTAL_TIME which is the
total Wall-Clock or Elapsed time (and NOT just the CPU time).

Note1: You can use the -o option to print out the actual sequences (to
stdout).  

Note2: -r option does a DFS generation of sequences, which consumes less 
memory, and thus can run even when there are long sequences.
You can try to omit the flag to do BFS generation if you want.

Note3: -e 1 option is a flag indicating spade to compute the support
of 2-sequences from scratch. The number 1 says there is only one DB
partition that will be inverted entirely in main memory. If the
original DB is large then this inversion will obviously take too much
time. So in this case I recommend dividing the DB into chunks of size
roughly 5MB (assuming there is 32MB available to the process). The
exttpose program is equiped to handle this case. If you specify a <-p
NUMPART> flag to exttpose it will divide the DB into NUMPART
chunks. Now you can run spade with -e NUMPART option. You must do this
if the DB is large otherwise the timings for spade will be
skewed. Generally, the more the partitions the better the running time
for spade. For example:        

        exttpose -i XXX -o XXX -l -x -s LMINSUP -p 10
        spade -i XXX -s MINSUP -e 10 -r

Note4: You can use the -m MEMSIZE option to increase the memory
available to the program. MEMSIZE is given in MB. For example if you
have 64 MB available, then use -m 64 (the default is 32MB).

4) In case you have another datafile not generated using the IBM datagen program, but in the same format CID TID #ITEMS ITEM-LIST, then you can use the getconf program to produce a conf file for it and then follow steps 2) and 3).

        Run: getconf -i XXX -o XXX

This assumes that XXX.data is in the required format, and produces XXX.conf

