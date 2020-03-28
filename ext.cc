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

#include "Eqclass.h"
#include "ext.h"
#include "partition.h"



#define TEMPISET "/var/spool/uucppublic/tmpiset"
#define TEMPSEQ "/var/spool/uucppublic/tmpseq"

struct timeval tp;

#define READRATIO 0.75
extern long AVAILMEM;
extern int DBASE_MAXITEM;
extern float DBASE_AVG_CUST_SZ;
extern float DBASE_AVG_TRANS_SZ;
extern int DBASE_NUM_TRANS;
extern int MINSUPPORT;
extern EqGrNode **eqgraph;
Extary **extary;
int **curit;
int *curcnt;
int *curitsz;
int EXTBLKSZ=100;
int *backidx,  numfreq;
int isetfd, seqfd, *eoffsets;
unsigned char *seq2, *itcnt2;
//for writing sequences to external file
#define seqitcntbufsz 1024
int seqbuf[seqitcntbufsz], itcntbuf[seqitcntbufsz];
int seqpos=0, itcntpos=0;

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

void get_ext_cust(int it, int data_fd, int *item_idx)
{
   //if (extary[it]->readall) return;

   if (extary[it]->bufpos &&
       extary[it]->bufpos <= extary[it]->bufsz*READRATIO) return;
   
   int bit = backidx[it];
   int blksz = min (extary[it]->bufsz,
                    extary[it]->endpos-extary[it]->filepos);
   //if (it == 400) cout << "READPOS " << extary[it]->filepos << " " <<
   //                  blksz << " " << extary[it]->endpos << endl;
   lseek(data_fd, (item_idx[bit]+extary[it]->filepos)*sizeof(int), SEEK_SET);
   int off = read(data_fd,(char *)extary[it]->buf,blksz*sizeof(int));
   if (off < blksz*sizeof(int)){
      perror ("ERROR READING");
      exit(errno);
   }

   if (blksz+extary[it]->filepos >= extary[it]->endpos)
      extary[it]->readall = 1;

   extary[it]->bufpos = 0;

   int ntid = extary[it]->buf[extary[it]->bufpos+1];
   if (extary[it]->bufsz < ntid+2){
      cout << "REALLOC " << ntid << endl;
      extary[it]->buf = (int *)realloc(extary[it]->buf, (ntid+2)*sizeof(int));
      if (extary[it]->buf == NULL){
         perror("ERROR: realloc extary in ext.cc");
         exit (-1);
      }
      extary[it]->bufsz = ntid+2;
      lseek(data_fd, (item_idx[bit]+extary[it]->filepos)*sizeof(int),
            SEEK_SET);
      off = read(data_fd,(char *)extary[it]->buf,(ntid+2)*sizeof(int));
   }
}

inline int get_next_cust(int it, int data_fd,int *item_idx)
{
   int ntid = extary[it]->ntid();
   extary[it]->bufpos += (ntid+2);
   extary[it]->filepos += (ntid+2);
   if (!extary[it]->check_empty()){
      ntid = extary[it]->ntid();
      if ((extary[it]->bufpos+ntid+2) > extary[it]->bufsz){
         get_ext_cust(it, data_fd, item_idx);
      }
      return 1;
   }
   return 0;
}

inline void fill_curit(int mval, Extary *ext, int it)
{
   int pos = ext->bufpos;
   int midx;
   while(pos < ext->endpos && ext->buf[pos] < mval+EXTBLKSZ){
      midx = ext->buf[pos] - mval;
      if (curcnt[midx]+1 > curitsz[midx]){
         curitsz[midx] = (int) (1.5*curitsz[midx]);
         curit[midx] = (int *)realloc(curit[midx], curitsz[midx]*sizeof(int));
         if (curit[midx] == NULL){
            perror("REALLCO  curit");
            exit(-1);
         }
      }
      curit[midx][curcnt[midx]++] = it;
      pos += (ext->buf[pos+1]+2);
   }
}

inline void fill_curit_ext(int mval, Extary *ext, int it, int mub)
{
   int pos = ext->bufpos;
   int fpos = ext->filepos;
   int midx;
   //if (it == 400) cout << "FILL " << pos << " " << ext->bufsz << " " <<
   //                  ext->filepos << " " << ext->endpos << " " <<
   //                  mub << endl;
   while(pos < ext->bufsz && fpos < ext->endpos && ext->buf[pos] < mub){
      midx = ext->buf[pos] - mval;
      //if (it == 400) cout << "FILLMIDX " <<  ext->buf[pos] << " " <<
      //                  midx << " " << pos << endl;
      if (curcnt[midx]+1 > curitsz[midx]){
         curitsz[midx] = (int) (1.5*curitsz[midx]);
         curit[midx] = (int *)realloc(curit[midx], curitsz[midx]*sizeof(int));
         if (curit[midx] == NULL){
            perror("REALLCO  curit");
            exit(-1);
         }
      }
      curit[midx][curcnt[midx]++] = it;
      fpos += (ext->buf[pos+1]+2);
      pos += (ext->buf[pos+1]+2);
   }
}

inline void find_matching_cust(int &mval, int mub=0)
{
   int i;

   if (num_partitions){
      for (i=0; i < numfreq; i++){
         if (!extary[i]->isempty()){
            fill_curit(mval, extary[i], i);
         }
      }
   }
   else{
      for (i=0; i < numfreq; i++){
         if (!extary[i]->isempty()){
            fill_curit_ext(mval, extary[i], i, mub);
         }
      }
   }
   mval+=EXTBLKSZ;
}

inline void process_cust(int curcnt, int *curit)
{
   int i,j,k,l;
   int it1, it2;

   for (i=0; i < curcnt; i++){ 
      it1 = curit[i];
      for (j=i; j < curcnt; j++){
         it2 = curit[j];
         if (extary[it1]->fbuf() < extary[it2]->lbuf()){
            if ((++seq2[it1*numfreq+it2]) == 0){
               if (seqpos+2 > seqitcntbufsz){
                  write(seqfd, (char *)seqbuf, seqpos*sizeof(int));
                  seqpos = 0;
               }
               seqbuf[seqpos++] = it1;
               seqbuf[seqpos++] = it2;
               //write(seqfd, (char *)&it1, sizeof(int));
               //write(seqfd, (char *)&it2, sizeof(int));
               //cout << "WROTE " << it1 << " "<< it2 << endl;
            }
         }
         if (j > i){
            if (extary[it2]->fbuf() < extary[it1]->lbuf()){
               if ((++seq2[it2*numfreq+it1]) == 0){
                  if (seqpos+2 > seqitcntbufsz){
                     write(seqfd, (char *)seqbuf, seqpos*sizeof(int));
                     seqpos = 0;
                  }
                  seqbuf[seqpos++] = it2;
                  seqbuf[seqpos++] = it1;
                  //write(seqfd, (char *)&it2, sizeof(int));
                  //write(seqfd, (char *)&it1, sizeof(int));
                  //cout << "WROTE " << it1 << " "<< it2 << endl;
               }
            }
            
            for (k=2, l=2; k < extary[it1]->ntid()+2 && l < extary[it2]->ntid()+2;){
               if (extary[it1]->buff()[k] < extary[it2]->buff()[l]) k++;
               else if (extary[it1]->buff()[k] > extary[it2]->buff()[l]) l++;
               else{
                  if ((++itcnt2[eoffsets[it1]-it1-1+it2]) == 0){                  
                     if (itcntpos+2 > seqitcntbufsz){
                        write(isetfd, (char *)itcntbuf, itcntpos*sizeof(int));
                        itcntpos = 0;
                     }
                     itcntbuf[itcntpos++] = it1;
                     itcntbuf[itcntpos++] = it2;
                     //write(isetfd, (char *)&it1, sizeof(int));
                     //write(isetfd, (char *)&it2, sizeof(int));
                     //cout << "WROTE2 " << it1 << " "<< it2 << endl;
                  }
                  break;
               }
            }
         }
      }
   }
}

void find_equal_cust(int data_fd, int &valid, int *item_idx)
{
   static int mval = 0;
   int i,k;
   int it;
   for (; valid > 1; ){
      //get trans with the same custid
      find_matching_cust(mval);
      
      //advance ext pointer for current items
      //cout << "MVAL " << mval-EXTBLKSZ;
      for (k=0; k < EXTBLKSZ; k++){
         //cout << " " << curcnt[k];
         //only if there are at least two items 
         if (curcnt[k] > 1) process_cust(curcnt[k], curit[k]);
         for (i=0; i < curcnt[k]; i++){
            it = curit[k][i];
            if (num_partitions){
               extary[it]->bufpos += (extary[it]->ntid()+2);
               if (extary[it]->bufpos >= extary[it]->endpos){
                  extary[it]->filepos = -1;
                  valid--;
               }
            }
            else if (!get_next_cust(it, data_fd, item_idx))
               valid--;
         }
         curcnt[k] = 0;
      }
      //cout << " " << valid << endl;
   }
   mval-=EXTBLKSZ;
}

void find_cust(int data_fd, int *item_idx)
{
   int i, k, it;
   int mval = 0;
   int mlb, mub;
   int mmin, mmax;
   int valid;
   while (mval <= DBASE_NUM_TRANS){
      mlb = -1;
      valid = numfreq;
      for (i=0; i < numfreq; i++){
         if (extary[i]->isempty() || extary[i]->readall) valid--;
         if (!extary[i]->isempty()){
            extary[i]->get_min_max(mmin, mmax);
            //cout << " (" << mmin << " " << mmax << ")";
            //cout << i << " MMAX " << mmax << endl;
            if (mmin != -1){
               if (mlb == -1){
                  mlb = mmin;
                  mub = mmax;
               }
               else{
                  if (mmin < mlb) mlb = mmin;
                  if (mmax < mub) mub = mmax;
               }
            }
         }
      }
      if (!valid) mub = DBASE_NUM_TRANS;
      cout << "MAXVAL " << mlb << " " << mub << " " << valid << endl;
      mval = mlb;
      while (mval <= mub){
         find_matching_cust(mval, min(mub+1, mval+EXTBLKSZ));
         
         //advance ext pointer for current items
         //cout << "MVAL " << mval-EXTBLKSZ;
         for (k=0; k < EXTBLKSZ; k++){
            //cout << " " << curcnt[k];
            //only if there are at least two items 
            if (curcnt[k] > 1) process_cust(curcnt[k], curit[k]);
            
            for (i=0; i < curcnt[k]; i++){
               it = curit[k][i];
               //if (mval+k-EXTBLKSZ == 75034)
               //   cout << " (" << it << ")";
               int ntid = extary[it]->ntid();
               extary[it]->bufpos += (ntid+2);
               extary[it]->filepos += (ntid+2);
               if (extary[it]->filepos >= extary[it]->endpos)
                  extary[it]->filepos = -1;
            }
            curcnt[k] = 0;
         }
         //cout << endl;
      }
      for (i=0; i < numfreq; i++){
         if (!extary[i]->isempty() && !extary[i]->readall)
            get_ext_cust(i, data_fd, item_idx);
      }
   }
}

void add_to_eqgraph(int oldit, char use_seq){
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

void write_middle_range(int &begi, int &begj, int endi, int endj,
                        int &l2cnt, int &oldit, char use_seq)
{ 
   if (MINSUPPORT >= 256) return;
   int i, j, ej;
   int lit, idx;
   
   for (i=begi; i <= endi; i++){
      if (i == begi) j = begj;
      else {
         if (use_seq) j = 0;
         else j = i+1;
      }
      if (i == endi) ej = endj;
      else ej = numfreq-1;
      if (!use_seq) idx = eoffsets[i]-i-1;
      for (; j <= ej; j++){
         if (use_seq){
            lit = (int) seq2[i*numfreq+j];
            seq2[i*numfreq+j] = 0;
         }
         else{
            lit = (int) itcnt2[idx+j];
            itcnt2[idx+j] = 0;
         }
         if (lit >= MINSUPPORT){
            if (oldit == -1) oldit = backidx[i];
            if (backidx[i] != oldit){
               add_to_eqgraph(oldit, use_seq);
               oldit = backidx[i];
               seqpos = 0;
            }
            //if (seqpos+3 > seqitcntbufsz){
            //   write(fd, (char *)seqbuf, seqpos*sizeof(int));
            //   seqpos = 0;
            //}
            //seqbuf[seqpos++] = backidx[i];
            seqbuf[seqpos++] = backidx[j];
            //seqbuf[seqpos++] = lit;
            //cout << backidx[i] << " " << backidx[j] << " " << lit << endl;
            l2cnt++;
         }
      }
   }
   if (endj == numfreq-1){
      begi = endi+1;
      if (use_seq) begj=0;
      else begj = begi+1;
   }
   else{
      begi = endi;
      begj = endj+1;
   }
}

void sort_get_l2(int &l2cnt, int fd, char use_seq)
{
   //write 2-itemsets counts to file

   int i, fcnt;
   long sortflen;
   int *sortary;
   int lit, outfd;
   int begi, begj;
   int itbuf[3];
   
   //if ((outfd = open(outfn, (O_WRONLY|O_CREAT|O_TRUNC), 0666)) < 0){
   //   perror("Can't open it2 file");
   //   exit (errno);      
   //}
   sortflen = lseek(fd, 0, SEEK_END);
   if (sortflen < 0){
      perror("SEEK SEQ");
      exit(errno);
   }
   cout << "SORT " << sortflen << endl;
   if (sortflen > 0){
      lseek(fd, 0, SEEK_SET);
#ifdef SGI
      sortary = (int *) mmap((char *)NULL, sortflen,
                             (PROT_READ|PROT_WRITE),
                             MAP_PRIVATE, fd, 0);
#else
      sortary = (int *) mmap((char *)NULL, sortflen,
                             (PROT_READ|PROT_WRITE),
                             (MAP_FILE|MAP_VARIABLE|MAP_PRIVATE), fd, 0);
#endif
      if (sortary == (int *)-1){
         perror("SEQFd MMAP ERROR");
         exit(errno);
      }
      
      qsort(sortary, (sortflen/sizeof(int))/2, 2*sizeof(int), cmp2it);
      
      //cout << "sorted " << (sortflen/sizeof(int))/2 << endl;
      itbuf[0] = itbuf[1] = -1;
      begi = begj = 0;
      seqpos = 0;
      int oldit = -1;
      fcnt = 0;
      for (i=0; i < (sortflen/sizeof(int))/2; i++){
         if (itbuf[0] != sortary[2*i] || itbuf[1] != sortary[2*i+1]){
            if (itbuf[0] != -1){
               write_middle_range(begi, begj, itbuf[0], itbuf[1],
                                  l2cnt, oldit, use_seq);
            }
            if (fcnt >= MINSUPPORT){
               if (oldit == -1) oldit = backidx[itbuf[0]];
               if (backidx[itbuf[0]] != oldit){
                  add_to_eqgraph(oldit, use_seq);
                  oldit = backidx[itbuf[0]];
                  seqpos = 0;
               }
               //if (seqpos+3 > seqitcntbufsz){
               //   write(outfd, (char *)seqbuf, seqpos*sizeof(int));
               //   seqpos = 0;
               //}
               //seqbuf[seqpos++] = backidx[itbuf[0]];
               seqbuf[seqpos++] = backidx[itbuf[1]];
               //seqbuf[seqpos++] = fcnt;
               //cout << itbuf[0] << " " << itbuf[1] << " " << itbuf[2] << endl;
               l2cnt++;
            }
            itbuf[0] = sortary[2*i];
            itbuf[1] = sortary[2*i+1];
            if (use_seq){
               fcnt = (int) seq2[itbuf[0]*numfreq+itbuf[1]];
               seq2[itbuf[0]*numfreq+itbuf[1]] = 0;
            }
            else{
               lit = itbuf[0];
               lit = (eoffsets[lit]-lit-1);
               fcnt = (int) itcnt2[lit+itbuf[1]];
               itcnt2[lit+itbuf[1]] = 0;            
            }
         }
         fcnt+=256;
      }
      write_middle_range(begi, begj, itbuf[0], itbuf[1], l2cnt, oldit, use_seq);
      if (fcnt >= MINSUPPORT){
         if (oldit == -1) oldit = backidx[itbuf[0]];
         if (backidx[itbuf[0]] != oldit){
            add_to_eqgraph(oldit, use_seq);
            oldit = backidx[itbuf[0]];
            seqpos = 0;
         }
         //if (seqpos+3 > seqitcntbufsz){
         //   write(outfd, (char *)seqbuf, seqpos*sizeof(int));
         //   seqpos = 0;
         //}
         //seqbuf[seqpos++] = backidx[itbuf[0]];
         seqbuf[seqpos++] = backidx[itbuf[1]];
         //seqbuf[seqpos++] = fcnt;
         //cout << itbuf[0] << " " << itbuf[1] << " " << itbuf[2] << endl;
         l2cnt++;
      }
      write_middle_range(begi, begj, numfreq-1, numfreq-1, l2cnt, oldit, use_seq);
      if (seqpos > 0){
         add_to_eqgraph(oldit, use_seq);
         seqpos = 0;
      }
      
      //if (seqpos > 0){
      //   write(outfd, (char *)seqbuf, seqpos*sizeof(int));
      //   seqpos = 0;
      //}
      //close(outfd);
      munmap((caddr_t)sortary, sortflen);
   }
   close(fd);
   //cout << "cnt " << l2cnt << endl;
}

int make_l1_pass(int data_fd, int *item_idx, int bflg, char use_newformat)
{
   int i,j;
   int sup,supsz;
   int bsz = 100;
   if (bflg){
      backidx = (int *) malloc (bsz*sizeof(int));
   }
   numfreq = 0;
   int ivalsz=100;
   int *ival = (int *)malloc(ivalsz*sizeof(int));
   for (i=0; i < DBASE_MAXITEM; i++){
      if (num_partitions) supsz = partition_get_idxsup(num_partitions, i);
      else supsz = (item_idx[i+1]-item_idx[i]);

      if (ivalsz < supsz){
         ivalsz = supsz;
         ival = (int *)realloc (ival, ivalsz*sizeof(int));
         if (ival == NULL){
            perror("IVAL NULL");
            exit(-1);
         }
      }
      if (num_partitions){
         partition_read_item(ival, num_partitions, i);
      }
      else{
         lseek(data_fd, item_idx[i]*sizeof(int), SEEK_SET);
         if (read(data_fd, (char *)ival,supsz*sizeof(int)) < 0){
            perror("read item1");
            exit(errno);
         }
      }
      sup = 0;
      if (use_newformat){
         int cid = -1;
         for (j=0; j < supsz; j+= 2){
            if (cid != ival[j]){
               sup++;
               cid = ival[j];
            }
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
               backidx = (int *)realloc (backidx, bsz*sizeof(int));
               if (backidx == NULL){
                  perror("BACKIDX NULL");
                  exit(-1);
               }
            }
            backidx[numfreq] = i;
         }
         numfreq++;
      }
   }
   cout << "NUMFREQ " << numfreq << endl;
   if (bflg){
      backidx = (int *)realloc(backidx, numfreq*sizeof(int));
      if (backidx == NULL){
         perror("BACKIDX NULL");
         exit(-1);
      }
   }

   free(ival);
   return numfreq;
}


void get_join(int *MAINBUF, int pidx, int i1, int i2,
              int ljoin, int ejoin, int mjoin)
{
   int i,j,k,l;
   int nval1, nval2;
   int sup1, sup2;
   int *it1, *it2;
   int itid, jtid;
   
   it1 = &MAINBUF[partition_get_idx(pidx, backidx[i1])];
   it2 = &MAINBUF[partition_get_idx(pidx, backidx[i2])];
   sup1 = partition_get_lidxsup(pidx, backidx[i1]);
   sup2 = partition_get_lidxsup(pidx, backidx[i2]);
   
   for (i=0,j=0; i < sup1 && j < sup2;){
      itid = it1[i];
      jtid = it2[j];
      if (itid > jtid){
         nval2 = it2[j+1];
         j += (nval2+2);
      }
      else if (itid < jtid){
         nval1 = it1[i+1];
         i += (nval1+2);
      }      
      else{
         i++;
         j++;
         nval1 = it1[i++];
         nval2 = it2[j++];
         if (ljoin){
            if (it1[i] < it2[j+nval2-1]){
               if ((++seq2[i1*numfreq+i2]) == 0){
                  if (seqpos+2 > seqitcntbufsz){
                     write(seqfd, (char *)seqbuf, seqpos*sizeof(int));
                     seqpos = 0;
                  }
                  seqbuf[seqpos++] = i1;
                  seqbuf[seqpos++] = i2;
               }
            }
         }
         if (ejoin){
            for (k=i, l=j; k<i+nval1 && l < j+nval2;){
               if (it1[k] < it2[l]) k++;
               else if (it1[k] > it2[l]) l++;
               else{
                  if ((++itcnt2[eoffsets[i1]-i1-1+i2]) == 0){                  
                     if (itcntpos+2 > seqitcntbufsz){
                        write(isetfd, (char *)itcntbuf, itcntpos*sizeof(int));
                        itcntpos = 0;
                     }
                     itcntbuf[itcntpos++] = i1;
                     itcntbuf[itcntpos++] = i2;
                  }
                  break;
               }
            }
         }
         if (mjoin){
            if (it2[j] < it1[i+nval1-1]){
               if ((++seq2[i2*numfreq+i1]) == 0){
                  if (seqpos+2 > seqitcntbufsz){
                     write(seqfd, (char *)seqbuf, seqpos*sizeof(int));
                     seqpos = 0;
                  }
                  seqbuf[seqpos++] = i2;
                  seqbuf[seqpos++] = i1;
               }
            }
         }
         i += nval1;
         j += nval2;
      }
   }
}

void partition_join_L2(int pidx, int *MAINBUF)
{
   int i,j;
   
   for (i=0; i < numfreq; i++){
      cout << "PRE " << i << endl;
      for (j=i; j < numfreq; j++){
         if (j==i) get_join(MAINBUF, pidx, i, j, 1, 0, 0);
         else get_join(MAINBUF, pidx, i, j, 1, 1, 1);
      }
   }
}




//void gen_ext_l2(char *it2f, char *seq2f,
//                int data_fd, int *item_idx)
int gen_ext_l2(int data_fd, int *item_idx)
{
   int i;
   double t1,t2;

   extary = new Extary *[numfreq];
   curit = new int *[EXTBLKSZ];
   curcnt = new int [EXTBLKSZ];
   curitsz = new int [EXTBLKSZ];
   for (i=0; i < EXTBLKSZ; i++){
      curitsz[i] = (int) (DBASE_AVG_CUST_SZ*DBASE_AVG_TRANS_SZ);
      curit[i] = (int *) malloc (curitsz[i]*sizeof(int));
      curcnt[i] = 0;
   }
   
   int extarysz  = AVAILMEM/numfreq;
   extarysz /= sizeof(int);
   cout << "EXTARYSZ " << extarysz << endl;
   if (extarysz < DBASE_AVG_CUST_SZ) extarysz = (int)DBASE_AVG_CUST_SZ;
   for (i=0; i < numfreq; i++){
      if (num_partitions) extary[i] = new Extary();
      else extary[i] = new Extary(extarysz);
   }
   
   if ((seqfd = open(TEMPSEQ, (O_RDWR|O_CREAT|O_TRUNC), 0666)) < 0){
      perror("Can't open out file");
      exit (errno);      
   }
   if ((isetfd = open(TEMPISET, (O_RDWR|O_CREAT|O_TRUNC), 0666)) < 0){
      perror("Can't open out file");
      exit (errno);      
   }

   seq2 = new unsigned char [numfreq*numfreq];
   if (seq2 == NULL){
      perror("SEQ MMAP ERROR");
      exit(errno);
   }
   bzero(seq2, numfreq*numfreq);

   itcnt2  = new unsigned char [(numfreq*(numfreq-1)/2)];
   if (itcnt2 == NULL){
      perror("ITCNT MMAP ERROR");
      exit(errno);
   }
   bzero(itcnt2, ((numfreq*(numfreq-1)/2)));

   eoffsets = new int [numfreq];
   int offt = 0;
   for (i=numfreq-1; i >= 0; i--){
      eoffsets[numfreq-i-1] = offt;
      offt += i;
   }

   seconds(t1);

   int iter = (num_partitions) ? num_partitions:1;
   int *MAINBUF=(int*)malloc(partition_get_max_blksz(num_partitions)*sizeof(int));
   for (int p=0; p < iter; p++){
      int valid = numfreq;
      if (num_partitions){
         partition_get_blk(MAINBUF, p);
         for (i=0; i < numfreq; i++){
            extary[i]->reset();
            extary[i]->endpos = partition_get_lidxsup(p, backidx[i]);
            if (extary[i]->endpos > 0)
               extary[i]->set_buf(&MAINBUF[partition_get_idx(p, backidx[i])]);
            else{
               extary[i]->set_buf(NULL);
               extary[i]->set_filepos(-1);
               valid--;
            }
         }
         find_equal_cust(data_fd, valid, item_idx);
         //cout << "VALID " << valid << endl;
      }
      else{
         for (i=0; i < numfreq; i++){
            //cout << "ITEM " << i << endl << flush;
            extary[i]->endpos = item_idx[backidx[i]+1]-item_idx[backidx[i]];
            get_ext_cust(i, data_fd, item_idx);
         }
         //find_equal_cust(data_fd, valid, item_idx);
         find_cust(data_fd, item_idx);
      }
      //cout << "GOT ALL CUST" << i << endl << flush;
      seconds(t2);
      cout << "PHASE " << p << " " << t2-t1 << endl;
   }
   if (MAINBUF) free(MAINBUF);
   if (seqpos > 0){
      write(seqfd, (char *)seqbuf, seqpos*sizeof(int));
      seqpos = 0;
   }
   if (itcntpos > 0){
      write(isetfd, (char *)itcntbuf, itcntpos*sizeof(int));
      itcntpos = 0;
   }
   //cout << "GOT ALL CUST" << i << endl << flush;

   for (i=0; i < numfreq; i++)
      if (extary[i]){
         if (num_partitions) extary[i]->set_buf(NULL);
         delete extary[i];
      }
   delete [] extary;
   for (i=0; i < EXTBLKSZ; i++) free(curit[i]);
   delete [] curit;
   delete [] curcnt;
   delete [] curitsz;
   
   seconds(t2);
   cout << "FIND_EQUAL " << t2-t1 << endl;
   int l2cnt = 0;
   seconds(t1);
   //iset must happen first
   sort_get_l2(l2cnt, isetfd, 0);
   seconds(t2);
   //cout << "SORT IT2 " << t2-t1 << endl;
   seconds(t1);
   sort_get_l2(l2cnt, seqfd, 1);
   seconds(t2);
   //cout << "SORT SEQ " << t2-t1 << endl;

   //reclaim memory
   free(backidx);
   delete [] eoffsets;

   unlink(TEMPSEQ);
   unlink(TEMPISET);
   return l2cnt;
}


inline void process_mult_cust(int curcnt, int *curit,
                              int **set_sup, int **seq_sup)
{
   int i,j,k,l;
   int it1, it2;

   for (i=0; i < curcnt; i++){ 
      it1 = curit[i];
      for (j=i; j < curcnt; j++){
         it2 = curit[j];
         if (seq_sup[it1] && extary[it1]->fbuf() < extary[it2]->lbuf()){
            seq_sup[it1][it2]++;
         }
         if (j > i){
            if (seq_sup[it2] && extary[it2]->fbuf() < extary[it1]->lbuf()){
               seq_sup[it2][it1]++;              
            }

            if (set_sup[it1]){
               for (k=2, l=2; k < extary[it1]->ntid()+2 &&
                       l < extary[it2]->ntid()+2;){
                  if (extary[it1]->buff()[k] < extary[it2]->buff()[l]) k++;
                  else if (extary[it1]->buff()[k] > extary[it2]->buff()[l]) l++;
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


void find_mult_cust(int data_fd, int &valid, int *item_idx,
                    int **set_sup, int **seq_sup, int &mval)
{
   int i,k;
   int it;
   for (; valid > 1; ){
      //get trans with the same custid
      find_matching_cust(mval);
      
      //advance ext pointer for current items
      //cout << "MVAL " << mval-EXTBLKSZ;
      for (k=0; k < EXTBLKSZ; k++){
         //cout << " " << curcnt[k];
         //only if there are at least two items 
         if (curcnt[k] > 1) process_mult_cust(curcnt[k], curit[k],
                                              set_sup, seq_sup);
         for (i=0; i < curcnt[k]; i++){
            it = curit[k][i];
            if (num_partitions){
               extary[it]->bufpos += (extary[it]->ntid()+2);
               if (extary[it]->bufpos >= extary[it]->endpos){
                  extary[it]->filepos = -1;
                  valid--;
               }
            }
            else if (!get_next_cust(it, data_fd, item_idx))
               valid--;
         }
         curcnt[k] = 0;
      }
      //cout << " " << valid << endl;
   }
   mval-=EXTBLKSZ;
}

void get_mult_L2(int idx, int **set_sup, int **seq_sup, int &l2cnt)
{
   int i;

   int bidx = backidx[idx];
   //2-iset support
   int tcnt = 0;
   for (i=idx+1; i < numfreq; i++){
      //cout << i << " " << idx << " " << i-idx-1 << " " << set_sup[idx][i-idx-1]
      //     << endl << flush;
      if (set_sup[idx][i-idx-1] >= MINSUPPORT)
         tcnt++;
   }
   if (tcnt > 0){
      if (eqgraph[bidx] == NULL) eqgraph[bidx] = new EqGrNode(tcnt);
      for (i=idx+1; i < numfreq; i++){
         if (set_sup[idx][i-idx-1] >= MINSUPPORT)
            eqgraph[bidx]->add_element(backidx[i]);
      }
   }
   l2cnt += tcnt;
   
   //2-seq support
   tcnt = 0;
   for (i=0; i < numfreq; i++){
      if (seq_sup[idx][i] >= MINSUPPORT) tcnt++;
   }
   if (tcnt > 0){
      if (eqgraph[bidx] == NULL) eqgraph[bidx] = new EqGrNode(0);
      if (eqgraph[bidx]->seqelements() == NULL){
         Array *sary = new Array(tcnt);
         eqgraph[bidx]->seqsetelements(sary);
      }
      for (i=0; i < numfreq; i++){
         if (seq_sup[idx][i] >= MINSUPPORT)
            eqgraph[bidx]->seqadd_element(backidx[i]);
      }
   }
   l2cnt += tcnt;
}


int gen_l2_multpass(int data_fd, int *item_idx)
{
   int i;
   double t1,t2;

   //initialize with the size of one partition
   int mem_used = partition_get_max_blksz(num_partitions);
   int *MAINBUF=(int*)malloc(mem_used*sizeof(int));
   extary = new Extary *[numfreq];

   curit = new int *[EXTBLKSZ];
   curcnt = new int [EXTBLKSZ];
   curitsz = new int [EXTBLKSZ];
   cout << "CURITSZ " << DBASE_AVG_CUST_SZ*DBASE_AVG_TRANS_SZ << endl;
   for (i=0; i < EXTBLKSZ; i++){
      curitsz[i] = (int) (DBASE_AVG_CUST_SZ*DBASE_AVG_TRANS_SZ);
      curit[i] = (int *) malloc (curitsz[i]*sizeof(int));
      curcnt[i] = 0;
   }

   int l2cnt=0;
   int extarysz  = AVAILMEM/numfreq;
   extarysz /= sizeof(int);
   //cout << "EXTARYSZ " << extarysz << endl;
   if (extarysz < DBASE_AVG_CUST_SZ) extarysz = (int)DBASE_AVG_CUST_SZ;
   for (i=0; i < numfreq; i++){
      if (num_partitions) extary[i] = new Extary();
      else extary[i] = new Extary(extarysz);
   }
   
   int **set_sup = new int *[numfreq];        // support for 2-itemsets
   bzero((char *)set_sup, numfreq*sizeof(int *));
   int **seq_sup = new int *[numfreq];        // support for 2-sequences
   bzero((char *)seq_sup, numfreq*sizeof(int *));
   int low, high;
   
   seconds(t1);

   int iter = (num_partitions) ? num_partitions:1;
   int mval;
   for (low = 0; low < numfreq; low = high){
      for (high = low; high < numfreq &&
              (mem_used+(2*numfreq-high-1)*sizeof(int)) < AVAILMEM; high++){
         if (numfreq-high-1 > 0){
            set_sup[high] = new int [numfreq-high-1];
            bzero((char *)set_sup[high], (numfreq-high-1)*sizeof(int));
         }
         seq_sup[high] = new int [numfreq];
         bzero((char *)seq_sup[high], numfreq*sizeof(int));
         mem_used += (2*numfreq-high-1) * sizeof(int);
         //cout << "MEMUSEDLOOP " << mem_used << " " <<
         //   (int)seq_sup[high] << " " << (int)set_sup[high] << endl;         
      }
      cout << "MEMUSEDLOOP " << mem_used << endl;  
      cout << "LOWHIGH " << high << " " << low << endl;
      mval = 0;
      for (int p=0; p < iter; p++){
         int valid = numfreq;
         if (num_partitions){
            partition_get_blk(MAINBUF, p);
            for (i=0; i < numfreq; i++){
               extary[i]->reset();
               extary[i]->endpos = partition_get_lidxsup(p, backidx[i]);
               if (extary[i]->endpos > 0)
                  extary[i]->set_buf(&MAINBUF[partition_get_idx(p, backidx[i])]);
               else{
                  extary[i]->set_buf(NULL);
                  extary[i]->set_filepos(-1);
                  valid--;
               }
            }
            cout << "VALID " << valid << endl;
            find_mult_cust(data_fd, valid, item_idx, set_sup, seq_sup, mval);
         }
         
         //cout << "GOT ALL CUST" << i << endl << flush;
         seconds(t2);
         cout << "PHASE " << p << " " << t2-t1 << endl;
      }
      
      // reclaim memory
      for (i = low; i < high; i++)
      {
         get_mult_L2(i, set_sup, seq_sup, l2cnt);
         delete [] set_sup[i];
         set_sup[i] = NULL;
         delete [] seq_sup[i];
         seq_sup[i] = NULL;
         mem_used -= (2*numfreq-i-1) * sizeof(int);
      }      
   }
   
   //reclaim memory
   if (MAINBUF) free(MAINBUF);
   delete [] set_sup;
   delete [] seq_sup;

   for (i=0; i < numfreq; i++)
      if (extary[i]){
         if (num_partitions) extary[i]->set_buf(NULL);
         delete extary[i];
      }
   delete [] extary;
   for (i=0; i < EXTBLKSZ; i++) free(curit[i]);
   delete [] curit;
   delete [] curcnt;
   delete [] curitsz;
   free(backidx);
   
   return l2cnt;
}
