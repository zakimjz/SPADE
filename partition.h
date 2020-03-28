#ifndef __PARTITION_H_
#define __PARTITION_H_
#include <sys/time.h>

#define ITSZ sizeof(int)
extern struct timeval tp;
#define seconds(tm) gettimeofday(&tp,(struct timezone *)0);\
tm=tp.tv_sec+tp.tv_usec/1000000.0

extern int num_partitions;
 
extern void partition_alloc(char *dataf, char *idxf);
extern void partition_dealloc();
extern void partition_get_blk(int *MAINBUF, int p);
extern int partition_get_blk_sz(int p);
extern int partition_get_max_blksz();
extern int partition_get_idxsup(int it);
extern int partition_get_lidxsup(int idx, int it);
extern int partition_get_idx(int idx, int it);
extern int partition_idxval(int it);
extern int *partition_idx(int idx);
extern void partition_read_item(int *ival, int it);
extern void partition_lclread_item(int *ival, int pnum, int it);
extern void partition_get_minmaxtid(int pnum, int it,
                                       int &minv, int &maxv);
extern void partition_get_minmaxcustid(int *backidx, int numit, int pnum,
                                       int &minv, int &maxv);
#endif// __PARTITION_H_
