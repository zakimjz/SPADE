#ifndef __EXT_H_
#define __EXT_H_

#include "partition.h"
#include "Eqclass.h"

#define ITSZ sizeof(int)

class invdb{
public:
   int numcust;
   int **curit;
   int *curcnt;
   int *curitsz;

   invdb(int sz);
   ~invdb();
   void incr(int sz);
   void incr_curit(int midx);
};



extern long AVAILMEM;
extern int DBASE_MAXITEM;
extern float DBASE_AVG_CUST_SZ;
extern float DBASE_AVG_TRANS_SZ;
extern int DBASE_NUM_TRANS;
extern int MINSUPPORT;
extern EqGrNode **eqgraph;
extern int num_partitions;
extern char use_newformat;

extern int make_l1_pass(int bflg);
extern int make_l2_pass();
extern int get_file_l2(char *it2f, char *seqf);

#endif //__EXT_H_
