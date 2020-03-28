#include <errno.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <malloc.h>
#include <strings.h>
 
#include "extl2.h"

#define seqitcntbufsz 4086

int seqbuf[seqitcntbufsz];
int seqpos=0;

unsigned int **set_sup, **seq_sup;
invdb *invDB;
int EXTBLKSZ;
int *backidx=NULL,  numfreq=0;
extern char print_seq;


invdb::invdb(int sz)
{
   int i;
   numcust = sz;
   curit =(int**) malloc (numcust*sizeof(int *));
   curcnt = (int *) malloc(numcust*ITSZ);
   curitsz = (int *) malloc(numcust*ITSZ);
   curitsz[0] = (int) (DBASE_AVG_CUST_SZ*DBASE_AVG_TRANS_SZ);
   for (i=0; i < numcust; i++){
      curitsz[i] = curitsz[0];
      curit[i] = (int *) malloc (curitsz[i]*ITSZ);
      curcnt[i] = 0;
   }   
}

invdb::~invdb()
{
   int i;
   for (i=0; i < numcust; i++){
      free (curit[i]);
   }
   free(curit);
   free(curcnt);
   free(curitsz);
}

void invdb::incr(int sz)
{
   int oldsz = numcust;
   numcust = sz;
   curit = (int **) realloc(curit, numcust*sizeof(int *));
   curcnt = (int *) realloc(curcnt, numcust*ITSZ);
   curitsz = (int *) realloc(curitsz, numcust*ITSZ);
   if (curit == NULL || curcnt == NULL || curitsz == NULL){
      perror("REALLCO  curit");
      exit(-1);
   }
   
   int i;   
   curitsz[oldsz] = (int) (DBASE_AVG_CUST_SZ*DBASE_AVG_TRANS_SZ);
   for (i=oldsz; i < numcust; i++){
      curitsz[i] = curitsz[oldsz];
      curit[i] = (int *) malloc (curitsz[i]*ITSZ);
      curcnt[i] = 0;      
   }
}

void invdb::incr_curit(int midx)
{
   //cout << "OLD: "  << curitsz[midx] << endl;
   
   curitsz[midx] = (int) (2*curitsz[midx]);
   //cout << "NEW: "  << curitsz[midx] << endl << endl;
   curit[midx] = (int *)realloc(curit[midx], curitsz[midx]*ITSZ);
   if (curit[midx] == NULL){
      perror("REALLCO  curit");
      exit(-1);
   }
}

int cmp2it(const void *a, const void *b)
{
   int **ary = (int **)a;
   int **bry = (int **)b;
   if (ary[0] < bry[0]) return -1;
   else if (ary[0] > bry[0]) return 1;
   else{
      if (ary[1] < bry[1]) return -1;
      else if (ary[1] > bry[1]) return 1;
      else return 0;
   }
}

int make_l1_pass(int bflg)
{
   int i,j;
   int sup,supsz;
   int bsz = 100;

   if (bflg) backidx = (int *) malloc (bsz*ITSZ);

   numfreq = 0;
   int ivalsz=100;
   int *ival = (int *)malloc(ivalsz*ITSZ);
   for (i=0; i < DBASE_MAXITEM; i++){
      supsz = partition_get_idxsup(i);
      if (ivalsz < supsz){
         ivalsz = supsz;
         ival = (int *)realloc (ival, ivalsz*ITSZ);
         if (ival == NULL){
            perror("IVAL NULL");
            exit(-1);
         }
      }
      partition_read_item(ival, i);
      sup = 0;
      if (use_newformat){
         int cid = -1;
         for (j=0; j < supsz; j+= 2){
            if (cid != ival[j]) sup++;
            cid = ival[j];
         }
      }
      else{
         for (j=0; j < supsz;){
            sup++;
            j += ival[j+1]+2;
         }
      }
      if (sup >= MINSUPPORT){
         if (bflg){
            if (numfreq+1 >= bsz){
               bsz = 2*bsz;
               backidx = (int *)realloc (backidx, bsz*ITSZ);
               if (backidx == NULL){
                  perror("BACKIDX NULL");
                  exit(-1);
               }
            }
            backidx[numfreq]  = i;
         } 
         if (print_seq) cout << i << " - " << sup << endl;
         
         numfreq++;
      }
   }
   //cout << "NUMFREQ " << numfreq << endl;

   if (bflg){
      backidx = (int *)realloc(backidx, numfreq*ITSZ);
      if (backidx == NULL){
         perror("BACKIDX NULL");
         exit(-1);
      }
   }
   
   free(ival);
   return numfreq;
}

void add_item_eqgraph(int oldit, char use_seq)
{
   int k;
   //cout << "ADDTOGRAPH " << seqpos << " " << oldit << " " << (int)use_seq << endl;
   if (eqgraph[oldit] == NULL){
      if (use_seq) eqgraph[oldit] = new EqGrNode(0);
      else eqgraph[oldit] = new EqGrNode(seqpos);
   }
   if (use_seq && eqgraph[oldit]->seqelements() == NULL){
      Array *sary = new Array(seqpos);
      eqgraph[oldit]->seqsetelements(sary);
   }
   for (k=0; k < seqpos; k++){
      if (use_seq) eqgraph[oldit]->seqadd_element(seqbuf[k]);
      else eqgraph[oldit]->add_element(seqbuf[k]);
   }   
}

void add_to_eqgraph(int &oldit, char use_seq, int it1, int it2){
   if (oldit == -1) oldit = it1;
   if (oldit != it1){
      add_item_eqgraph(oldit, use_seq);
      oldit = it1;
      seqpos = 0;
   }
   seqbuf[seqpos++] = it2;
}

void process_cust_invert(int curcnt, int *curit){
   int i,j,k,l;
   int nv1, nv2;
   int it1, it2;
   for (i=0; i < curcnt; i=nv1){
      nv1 = i;
      it1 = curit[i];
      while (it1 == curit[nv1] && nv1 < curcnt) nv1+=2;
      for (j=i; j < curcnt; j=nv2){
         nv2 = j;
         it2 = curit[j];
         while (it2 == curit[nv2] && nv2 < curcnt) nv2+=2;
         if (seq_sup[it1] && curit[i+1] < curit[nv2-1]){
            seq_sup[it1][it2]++;
         }
         if (j>i){
            if (seq_sup[it2] && curit[j+1] < curit[nv1-1]){
               seq_sup[it2][it1]++;
            }
            if (set_sup[it1]){
               for (k=i, l=j; k < nv1 && l < nv2;){
                  if (curit[k+1] > curit[l+1]) l+=2;
                  else if (curit[k+1] < curit[l+1]) k+=2;
                  else{
                     set_sup[it1][it2-it1-1]++;
                     break;
                  }
               }
            }
         }
      }
   }
}

void process_invert(int pnum)
{
   int i,k;
   int minv, maxv;
   partition_get_minmaxcustid(backidx, numfreq, pnum, minv, maxv);
   //cout << "MVAL " << minv << " " << maxv << endl;
   if (invDB->numcust < maxv-minv+1)
      invDB->incr(maxv-minv+1);
   
   int supsz;
   int ivalsz=0;
   int *ival = NULL;
   for (i=0; i < numfreq; i++){
      supsz = partition_get_lidxsup(pnum, backidx[i]);
      if (ivalsz < supsz){
         ivalsz = supsz;
         ival = (int *)realloc (ival, ivalsz*ITSZ);
         if (ival == NULL){
            perror("IVAL NULL");
            exit(-1);
         }
      }
      partition_lclread_item(ival, pnum, backidx[i]);
      int cid;
      int midx;
      for (int pos=0; pos < supsz; pos += 2)
      {
         //if (cid != ival[pos]){
         cid = ival[pos];
         midx = cid - minv;
         //if (midx >= maxv-minv+1){
         //   perror("EXCEEDED BOUNDS\n");
         //}
         //cout << "MIDX " << midx << endl;
         if (invDB->curcnt[midx]+2 > invDB->curitsz[midx]){
            invDB->incr_curit(midx);
         }
         invDB->curit[midx][invDB->curcnt[midx]++] = i;
         invDB->curit[midx][invDB->curcnt[midx]++] = ival[pos+1];            
         //}
      }
   }
   for (k=0; k < maxv-minv+1; k++){
      //cout << "MVAL " << mval+k << " " << curcnt[k] << endl << flush;
      if (invDB->curcnt[k] > 0){
         process_cust_invert(invDB->curcnt[k], invDB->curit[k]);
      }
      invDB->curcnt[k] = 0;
   }
}

void get_F2(int &l2cnt)
{
   int j,k;
   int fcnt;

   for (j=0; j < numfreq;j++){
      seqpos = 0;
      if (set_sup[j]){
         for (k=j+1; k < numfreq;k++){
            fcnt = set_sup[j][k-j-1];
            if (fcnt >= MINSUPPORT){
               seqbuf[seqpos++] = backidx[k];
               if (print_seq){
                  cout << backidx[j] << " " << backidx[k] 
                       << " - " << fcnt << endl;
               }
               l2cnt++;
            }
         }
      }
      if (seqpos > 0) add_item_eqgraph(backidx[j], 0);
      seqpos = 0;
      if (seq_sup[j]){
         for (k=0; k < numfreq;k++){
            fcnt = seq_sup[j][k];
            if (fcnt >= MINSUPPORT){
               seqbuf[seqpos++] = backidx[k];
               if (print_seq){
                  cout << backidx[j] << " -> " << backidx[k] 
                       << " - " << fcnt << endl;
               }
               l2cnt++;
            }
         }         
      }
      if (seqpos > 0) add_item_eqgraph(backidx[j], 1);
   }
}


int make_l2_pass()
{
   int i;

   int l2cnt=0;
   int mem_used=0;

   EXTBLKSZ = num_partitions+(DBASE_NUM_TRANS+num_partitions-1)/num_partitions;
   int tsz = (int) (DBASE_AVG_CUST_SZ*DBASE_AVG_TRANS_SZ);
   invDB = new invdb(EXTBLKSZ);
   //mem_used += EXTBLKSZ*2*ITSZ;
   //mem_used += (int) (EXTBLKSZ*tsz*ITSZ);
   //cout << "CURITSZ " << tsz << " " << EXTBLKSZ << " " << mem_used << endl;
   
   set_sup = new unsigned int *[numfreq];        // support for 2-itemsets
   bzero((char *)set_sup, numfreq*sizeof(unsigned int *));   
   mem_used += numfreq*sizeof(unsigned int *);      
   seq_sup = new unsigned int *[numfreq];        // support for 2-itemsets
   bzero((char *)seq_sup, numfreq*sizeof(unsigned int *));   
   mem_used += numfreq*sizeof(unsigned int *);
   
   int low, high;
      
   int itsz = sizeof(unsigned int);
   for (low = 0; low < numfreq; low = high){
      for (high = low; high < numfreq &&
              (mem_used+2*numfreq-high-1)*itsz < AVAILMEM; high++){
         if (numfreq-high-1 > 0){
            set_sup[high] = new unsigned int [numfreq-high-1];
            bzero((char *)set_sup[high], (numfreq-high-1)*itsz);
         }
         seq_sup[high] = new unsigned int [numfreq];
         bzero((char *)seq_sup[high], numfreq*itsz);
         mem_used += (2*numfreq-high-1) * itsz;
         //cout << "MEMUSEDLOOP " << mem_used << " " << endl;
      }
      //cout << "MEMUSEDLOOP " << mem_used << endl;  
      //cout << "LOWHIGH " << high << " " << low << endl;
      for (int p=0; p < num_partitions; p++){
         process_invert(p);
      }
      get_F2(l2cnt);
      // reclaim memory
      for (i = low; i < high; i++)
      {
         if (set_sup[i]) delete [] set_sup[i];
         set_sup[i] = NULL;
         if (seq_sup[i]) delete [] seq_sup[i];
         seq_sup[i] = NULL;
         mem_used -= (2*numfreq-i-1) * itsz;
      }
   }

   delete [] set_sup;
   delete [] seq_sup;
   delete invDB;
   free(backidx);
   return l2cnt;
}


void get_l2file(char *fname, char use_seq, int &l2cnt)
{
   int *cntary;
   int fd = open(fname, O_RDONLY);
   if (fd < 1){
      perror("can't open l2 file");
      exit(errno);
   }   
   int flen = lseek(fd,0,SEEK_END);
   if (flen > 0){
#ifndef DEC
      cntary = (int *) mmap((char *)NULL, flen, PROT_READ,
                             MAP_PRIVATE, fd, 0);
#else
      cntary = (int *) mmap((char *)NULL, flen, PROT_READ,
                             (MAP_FILE|MAP_VARIABLE|MAP_PRIVATE), fd, 0);
#endif
      if (cntary == (int *)-1){
         perror("MMAP ERROR:cntary");  
         exit(errno);
      }
      
      // build eqgraph -- large 2-itemset relations
      int lim = flen/ITSZ;
      int oldit = -1;
      seqpos = 0;
      for (int i=0; i < lim; i += 3){
         if (cntary[i+2] >= MINSUPPORT){
            add_to_eqgraph(oldit, use_seq, cntary[i], cntary[i+1]);
            if (print_seq){
               cout << cntary[i] << ((use_seq)?" -> ":" ") << cntary[i+1]
                    << " - " << cntary[i+2] << endl;
            }
            
            l2cnt++;
         }
      }
      if (seqpos > 0) add_item_eqgraph(oldit, use_seq);
      munmap((caddr_t)cntary, flen);
   }
   close(fd);
}

int get_file_l2(char *it2f, char *seqf)
{
   int l2cnt = 0;

   get_l2file(it2f, 0, l2cnt);
   get_l2file(seqf, 1, l2cnt);

   //cout << "L2 : " << l2cnt << endl;
   return l2cnt;
}
