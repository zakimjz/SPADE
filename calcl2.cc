#include <errno.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include "calcl2.h"
#include "partition.h"
#include "calcdb.h"

using namespace std;

#define MBYTE (1024*1024)
#define ITSZ sizeof(int)
#define TEMPISET "/var/spool/uucppublic/tmpiset"
#define TEMPSEQ "/var/spool/uucppublic/tmpseq"

long AVAILMEM = 32*MBYTE;
int DBASE_NUM_TRANS;
int DBASE_MAXITEM;
float DBASE_AVG_TRANS_SZ;
float DBASE_AVG_CUST_SZ; 
int DBASE_TOT_TRANS;

double MINSUP_PER;
int MINSUPPORT;

double L2TIME=0;
char dataf[300];
char idxf[300];
char conf[300];
char input[300];

int **set_sup, **seq_sup;
int **set_cid, **seq_cid;
unsigned char **cset_sup, **cseq_sup;
Extary **extary;
int **curit;
int *curcnt;
int *curitsz;
int EXTBLKSZ=100;
int *backidx,  numfreq;
char new_format=0;
char use_invert=0;
char use_charary = 0;
char use_horz = 0;
char use_apralg = 0;

#define seqitcntbufsz 1024
int seqbuf[seqitcntbufsz], itcntbuf[seqitcntbufsz];
int seqpos=0, itcntpos=0;
int isetfd, seqfd;
ofstream fout;

void parse_args(int argc, char **argv)
{
   extern char * optarg;
   int c;
   
   if (argc < 2)
      cout << "usage: seq -i<infile> -o<outfile> -s<support>\n";
   else{
      while ((c=getopt(argc,argv,"i:m:n:s:fvcha"))!=-1){
         switch(c){
         case 'i': //input file
            sprintf(input, "%s", optarg);
            sprintf(dataf,"%s.tpose", optarg);
            sprintf(idxf,"%s.idx", optarg);
            sprintf(conf,"%s.conf", optarg);
            break;
         case 'm': //amount of mem available
            AVAILMEM = (long) atof(optarg)*MBYTE;
            break;
         case 'n': //#partitions for optimized L2 calculation
            strtok(optarg, " " );
            num_partitions = atoi(optarg);
            if (optarg = strtok(NULL, " "))
               EXTBLKSZ = atoi(optarg);
            break;
         case 's': //min support
            MINSUP_PER = atof(optarg);
            break;
         case 'f':
            new_format = 1;
            break;
         case 'v':
            use_invert = 1;
            new_format = 1;
            //use_charary =1;
            break;
         case 'c':
            use_charary = 1;
            //new_format = 1;
            break;
         case 'h':
            use_horz = 1;
            break;
         case 'a':
            use_apralg = 1;
            use_horz = 1;
            break;
         }
      }
   }
   c= open(conf, O_RDONLY);
   if (c < 0){
      perror("ERROR: invalid conf file\n");
      exit(errno);
   }
   read(c,(char *)&DBASE_NUM_TRANS,ITSZ);
   MINSUPPORT = (int) ceil(MINSUP_PER*DBASE_NUM_TRANS);
   //ensure that support is at least 2
   if (MINSUPPORT < 2) MINSUPPORT = 2;
   cout << "MINSUPPORT " << MINSUPPORT << " " << DBASE_NUM_TRANS << endl;
   read(c,(char *)&DBASE_MAXITEM,ITSZ);
   read(c,(char *)&DBASE_AVG_CUST_SZ,sizeof(float));
   read(c,(char *)&DBASE_AVG_TRANS_SZ,sizeof(float));
   read(c,(char *)&DBASE_TOT_TRANS,ITSZ);
   cout << "CONF " << DBASE_NUM_TRANS << " " << DBASE_MAXITEM << " "
        << DBASE_AVG_TRANS_SZ << endl;
   close(c);
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

int make_l1_pass()
{
   int i,j;
   int sup,supsz;
   int bsz = 100;

   backidx = (int *) malloc (bsz*ITSZ);
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
      if (new_format){
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
         if (numfreq+1 >= bsz){
            bsz = 2*bsz;
            backidx = (int *)realloc (backidx, bsz*ITSZ);
            if (backidx == NULL){
               perror("BACKIDX NULL");
               exit(-1);
            }
         }
         backidx[numfreq] = i;
         numfreq++;
      }
   }
   cout << "NUMFREQ " << numfreq << endl;
   backidx = (int *)realloc(backidx, numfreq*ITSZ);
   if (backidx == NULL){
      perror("BACKIDX NULL");
      exit(-1);
   }
   
   free(ival);
   return numfreq;
}



void add_uniq(int it1, int it2, int *&ary, int &arysz, int &cursz)
{
   int i;
   for (i=0; i < cursz; i+=2){
      if (ary[i] == it1 && ary[i+1] == it2) return;
   }
   if (cursz+2 > arysz){
      arysz = (int)(1.5*arysz);
      ary = (int *) realloc((char *)ary, arysz*sizeof(int));
      cout << "REALLOC " << arysz << endl;
      if (ary == NULL){
         perror("no alloc in add_uniq");
         exit(errno);
      }
   }
   ary[cursz++] = it1;
   ary[cursz++] = it2;
   //cout << "ADDED " << it1 << " " << it2 << " " << cursz <<endl;
}

void process_cust_invert(int curcnt, int *curit)
{
   int i,j;
   static int *tmp;
   static int tmpsz=0;
   static char first=1;
   if (first){
      first = 0;
      tmpsz = seqitcntbufsz;
      tmp = (int *) malloc (tmpsz*sizeof(int));
   }
   int tcnt;
   //cout << "CUST " << curcnt << endl;
   
   //for (i=0; i < curcnt; i+=2)
   //   cout << curit[i] << " " << curit[i+1]  << endl;
   for (i=0; i < curcnt; i+=2){
      //cout << "TRNS " << curit[i] << " " << curit[i+1] << endl;
      tcnt = 0;
      for (j=i; j < curcnt; j+=2){
         if (seq_sup[curit[i]] && curit[i+1] < curit[j+1]){
            add_uniq(curit[i], curit[j], tmp, tmpsz, tcnt);
         }
         if (j>i && curit[i] != curit[j]){
            if (seq_sup[curit[j]] && curit[i+1] > curit[j+1]){
               add_uniq(curit[j], curit[i], tmp, tmpsz, tcnt);
            }
            if (set_sup[curit[i]] && curit[i+1] == curit[j+1]){
               add_uniq(curit[j], -1, tmp, tmpsz, tcnt);
            }
         }
      }
      if (tcnt > 0){
         for (j=0; j < tcnt; j+=2){
            //cout << "TTT "  << tmp[j] << " " << tmp[j+1] << endl;
            if (tmp[j+1] == -1){
               set_sup[curit[i]][tmp[j]-curit[i]-1]++;
            }
            else{
               seq_sup[tmp[j]][tmp[j+1]]++;
            }
         }
      }
      //exit(-1);
   }
}

void process_charcust_invert2(int curcnt, int *curit){
   int i,j,k,l;
   int nv1, nv2;
   int it1, it2;
   for (i=0; i < curcnt; i=nv1){
      nv1 = i;
      it1 = curit[i];
      while (it1 == curit[nv1]) nv1+=2;
      for (j=i; j < curcnt; j=nv2){
         nv2 = j;
         it2 = curit[j];
         while (it2 == curit[nv2]) nv2+=2;
         if (cseq_sup[it1] && curit[i+1] < curit[nv2-1]){
            if ((++cseq_sup[it1][it2]) == 0){
               if (seqpos+2 > seqitcntbufsz){
                  write(seqfd, (char *)seqbuf, seqpos*sizeof(int));
                  seqpos = 0;
               }
               seqbuf[seqpos++] = it1;
               seqbuf[seqpos++] = it2;
            }
         }
         if (j>i){
            if (cseq_sup[it2] && curit[j+1] < curit[nv1-1]){
               if ((++cseq_sup[it2][it1]) == 0){
                  if (seqpos+2 > seqitcntbufsz){
                     write(seqfd, (char *)seqbuf, seqpos*sizeof(int));
                     seqpos = 0;
                  }
                  seqbuf[seqpos++] = it2;
                  seqbuf[seqpos++] = it1;
               }
            }
            if (cset_sup[it1]){
               for (k=i, l=j; k < nv1 && l < nv2;){
                  if (curit[k+1] > curit[l+1]) l+=2;
                  else if (curit[k+1] < curit[l+1]) k+=2;
                  else{
                     if ((++cset_sup[it1][it2-it1-1]) == 0){
                        if (itcntpos+2 > seqitcntbufsz){
                           write(isetfd, (char *)itcntbuf, itcntpos*sizeof(int));
                           itcntpos = 0;
                        }
                        itcntbuf[itcntpos++] = it1;
                        itcntbuf[itcntpos++] = it2;
                     }
                     break;
                  }
               }
            }
         }
      }
   }
}

void process_cust_invert2(int curcnt, int *curit){
   int i,j,k,l;
   int nv1, nv2;
   int it1, it2;
   for (i=0; i < curcnt; i=nv1){
      nv1 = i;
      it1 = curit[i];
      while (it1 == curit[nv1]) nv1+=2;
      for (j=i; j < curcnt; j=nv2){
         nv2 = j;
         it2 = curit[j];
         while (it2 == curit[nv2]) nv2+=2;
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

void process_invert(int pnum, int &mval)
{
   int i,k;

   cout << "MVAL " << mval << endl;
   int supsz;
   int ivalsz=0;
   int *ival = NULL;
   int maxval = 0;
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
      int cid = -1;
      int midx;
      for (int pos=0; pos < supsz; pos += 2)
      {
         //if (cid != ival[pos]){
         cid = ival[pos];
         midx = cid - mval;
         if (maxval < cid ) maxval = cid;
         //cout << "MIDX " << midx << endl;
         if (curcnt[midx]+2 > curitsz[midx]){
            curitsz[midx] = (int) (1.5*curitsz[midx]);
            curit[midx] = (int *)realloc(curit[midx], curitsz[midx]*ITSZ);
            if (curit[midx] == NULL){
               perror("REALLCO  curit");
               exit(-1);
            }
         }
         curit[midx][curcnt[midx]++] = i;
         curit[midx][curcnt[midx]++] = ival[pos+1];            
         //}
      }
   }
   for (k=0; k < EXTBLKSZ; k++){
      //cout << "MVAL " << mval+k << " " << curcnt[k] << endl << flush;
      if (curcnt[k] > 0){
         if (use_charary) process_charcust_invert2(curcnt[k], curit[k]);
         else process_cust_invert2(curcnt[k], curit[k]);
      }
      curcnt[k] = 0;
   }
   mval = maxval+1;
   cout << "MVAL END " << mval << endl;
}

inline void fill_curit(int mval, Extary *ext, int it)
{
   int pos = ext->bufpos;
   int cid = -1;
   int midx;
   if (new_format){
      while (pos < ext->endpos && ext->buf[pos] < mval+EXTBLKSZ){
         if (ext->buf[pos] != cid){
            cid = ext->buf[pos];
            midx = cid - mval;
            if (curcnt[midx]+1 > curitsz[midx]){
               curitsz[midx] = (int) (1.5*curitsz[midx]);
               curit[midx] = (int *)realloc(curit[midx], curitsz[midx]*ITSZ);
               if (curit[midx] == NULL){
                  perror("REALLCO  curit");
                  exit(-1);
               }
            }
            curit[midx][curcnt[midx]++] = it;
         }
         pos += 2;
      }
   }
   else{
      while(pos < ext->endpos && ext->buf[pos] < mval+EXTBLKSZ){
         midx = ext->buf[pos] - mval;
         if (curcnt[midx]+1 > curitsz[midx]){
            curitsz[midx] = (int) (1.5*curitsz[midx]);
            curit[midx] = (int *)realloc(curit[midx], curitsz[midx]*ITSZ);
            if (curit[midx] == NULL){
               perror("REALLCO  curit");
               exit(-1);
            }
         }
         curit[midx][curcnt[midx]++] = it;
         pos += (ext->buf[pos+1]+2);
      }
   }
}

inline void find_matching_cust(int &mval, int mub=0)
{
   int i;

   for (i=0; i < numfreq; i++){
      if (!extary[i]->isempty()){
         fill_curit(mval, extary[i], i);
      }
   }
   mval+=EXTBLKSZ;
}

void process_charmult_cust(int curcnt, int *curit, int mval)
{
   int i,j,k,l;
   int it1, it2;
   int nval1, nval2;
   if (new_format){
      for (i=0; i < curcnt; i++){
         it1 = curit[i];
         nval1 = extary[it1]->bufpos;
         while(mval == extary[it1]->buf[nval1]) nval1 += 2;
         for (j=i; j < curcnt; j++){
            it2 = curit[j];
            nval2 = extary[it2]->bufpos;
            while(mval == extary[it2]->buf[nval2]) nval2 += 2;

            if (cseq_sup[it1] && extary[it1]->ntid() < extary[it2]->buf[nval2-1]){
               if ((++cseq_sup[it1][it2]) == 0){
                  if (seqpos+2 > seqitcntbufsz){
                     write(seqfd, (char *)seqbuf, seqpos*sizeof(int));
                     seqpos = 0;
                  }
                  seqbuf[seqpos++] = it1;
                  seqbuf[seqpos++] = it2;
               }
            }
            if (j > i){
               if (cseq_sup[it2] &&
                   extary[it2]->ntid() < extary[it1]->buf[nval1-1]){
                  if ((++cseq_sup[it2][it1]) == 0){
                     if (seqpos+2 > seqitcntbufsz){
                        write(seqfd, (char *)seqbuf, seqpos*sizeof(int));
                        seqpos = 0;
                     }
                     seqbuf[seqpos++] = it2;
                     seqbuf[seqpos++] = it1;
                  }
               }
               
               if (cset_sup[it1]){
                  for (k=extary[it1]->bufpos, l=extary[it2]->bufpos;
                       k < nval1 && l < nval2;){
                     if (extary[it1]->buf[k+1] < extary[it2]->buf[l+1]) k+=2;
                     else if (extary[it1]->buf[k+1] > extary[it2]->buf[l+1]) l+=2;
                     else{
                        if ((++cset_sup[it1][it2-it1-1]) == 0){
                           if (itcntpos+2 > seqitcntbufsz){
                              write(isetfd, (char *)itcntbuf, itcntpos*sizeof(int));
                              itcntpos = 0;
                           }
                           itcntbuf[itcntpos++] = it1;
                           itcntbuf[itcntpos++] = it2;
                        }
                        break;
                     }
                  }
               }
            }            
         }
      }
   }
   else{
      for (i=0; i < curcnt; i++){
         it1 = curit[i];
         for (j=i; j < curcnt; j++){
            it2 = curit[j];
            if (cseq_sup[it1] && extary[it1]->fbuf() < extary[it2]->lbuf()){
               if ((++cseq_sup[it1][it2]) == 0){
                  if (seqpos+2 > seqitcntbufsz){
                     write(seqfd, (char *)seqbuf, seqpos*sizeof(int));
                     seqpos = 0;
                  }
                  seqbuf[seqpos++] = it1;
                  seqbuf[seqpos++] = it2;
               }
            }
            if (j > i){
               if (cseq_sup[it2] && extary[it2]->fbuf() < extary[it1]->lbuf()){
                  if ((++cseq_sup[it2][it1]) == 0){
                     if (seqpos+2 > seqitcntbufsz){
                        write(seqfd, (char *)seqbuf, seqpos*sizeof(int));
                        seqpos = 0;
                     }
                     seqbuf[seqpos++] = it2;
                     seqbuf[seqpos++] = it1;
                  }
               }
               
               if (cset_sup[it1]){
                  for (k=2, l=2; k < extary[it1]->ntid()+2 &&
                          l < extary[it2]->ntid()+2;){
                     if (extary[it1]->buff()[k] < extary[it2]->buff()[l]) k++;
                     else if (extary[it1]->buff()[k] > extary[it2]->buff()[l]) l++;
                     else{
                        if ((++cset_sup[it1][it2-it1-1]) == 0){
                           if (itcntpos+2 > seqitcntbufsz){
                              write(isetfd, (char *)itcntbuf, itcntpos*sizeof(int));
                              itcntpos = 0;
                           }
                           itcntbuf[itcntpos++] = it1;
                           itcntbuf[itcntpos++] = it2;
                        }
                        break;
                     }
                  }
               }
            }
         }
      }
   }
}

void process_mult_cust(int curcnt, int *curit, int mval)
{
   int i,j,k,l;
   int it1, it2;
   int nval1, nval2;
   if (new_format){
      for (i=0; i < curcnt; i++){
         it1 = curit[i];
         nval1 = extary[it1]->bufpos;
         while(mval == extary[it1]->buf[nval1]) nval1 += 2;
         for (j=i; j < curcnt; j++){
            it2 = curit[j];
            nval2 = extary[it2]->bufpos;
            while(mval == extary[it2]->buf[nval2]) nval2 += 2;

            if (seq_sup[it1] && extary[it1]->ntid() < extary[it2]->buf[nval2-1]){
               seq_sup[it1][it2]++;
            }
            if (j > i){
               if (seq_sup[it2] && extary[it2]->ntid() < extary[it1]->buf[nval1-1]){
                  seq_sup[it2][it1]++;              
               }
               
               if (set_sup[it1]){
                  for (k=extary[it1]->bufpos, l=extary[it2]->bufpos;
                       k < nval1 && l < nval2;){
                     if (extary[it1]->buf[k+1] < extary[it2]->buf[l+1]) k+=2;
                     else if (extary[it1]->buf[k+1] > extary[it2]->buf[l+1]) l+=2;
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
   else{
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
}

void find_mult_cust(int &valid, int &mval)
{
   int i,k;
   int it;
   for (; valid > 1; ){
      //get trans with the same custid
      find_matching_cust(mval);
      
      //advance ext pointer for current items
      //cout << "MVAL " << mval-EXTBLKSZ;
      for (k=0; k < EXTBLKSZ; k++){
         //cout << "MVAL " << mval-EXTBLKSZ+k << " " << curcnt[k] << endl;
         //only if there are at least two items 
         if (curcnt[k] > 1)
            if (use_charary)
               process_charmult_cust(curcnt[k], curit[k], mval-EXTBLKSZ+k);
            else process_mult_cust(curcnt[k], curit[k], mval-EXTBLKSZ+k);
         for (i=0; i < curcnt[k]; i++){
            it = curit[k][i];
            if (new_format){
               while(extary[it]->custid() == k+mval-EXTBLKSZ)
                  extary[it]->bufpos += 2;
            }
            else extary[it]->bufpos += (extary[it]->ntid()+2);
            if (extary[it]->bufpos >= extary[it]->endpos){
               extary[it]->filepos = -1;
               valid--;
            }
         }
         curcnt[k] = 0;
      }
      //cout << " " << valid << endl;
   }
   mval-=EXTBLKSZ;
}

void get_mult_L2(int idx, int &l2cnt)
{
   int i;

   //2-iset support
   int tcnt = 0;
   for (i=idx+1; i < numfreq; i++){
      if (set_sup[idx][i-idx-1] >= MINSUPPORT){
         //cout << idx << " " << i << " " << set_sup[idx][i-idx-1] << endl;
         tcnt++;
      }
   }
   l2cnt += tcnt;
   
   //2-seq support
   tcnt = 0;
   for (i=0; i < numfreq; i++){
      if (seq_sup[idx][i] >= MINSUPPORT){
         //cout << idx << " " << i << " " << seq_sup[idx][i] << endl;
         tcnt++;
      }
   }
   l2cnt += tcnt;
}

void write_middle_range(int &begi, int &begj, int endi, int endj,
                        int &l2cnt, char use_seq)
{ 
   if (MINSUPPORT >= 256) return;
   int i, j, ej;
   int lit;
   
   for (i=begi; i <= endi; i++){
      if (i == begi){
         j = begj;
      }
      else {
         if (use_seq) j = 0;
         else j = i+1;
      }
      if (i == endi) ej = endj;
      else ej = numfreq-1;

      if (!cseq_sup[i] && !cset_sup[i]) continue;
      for (; j <= ej; j++){
         if (use_seq && cseq_sup[i]){
            lit = (int) cseq_sup[i][j];
            cseq_sup[i][j] = 0;
         }
         else if (j>i && cset_sup[i]){
            lit = (int) cset_sup[i][j-i-1];
            cset_sup[i][j-i-1] = 0;
         }
         else lit = 0;
         if (lit >= MINSUPPORT){
            //cout << i << " " << j << " " << lit << endl;
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
   int itbuf[3];
   int begi, begj;
   
   sortflen = lseek(fd, 0, SEEK_CUR);
   if (sortflen < 0){
      perror("SEEK SEQ");
      exit(errno);
   }
   cout << "SORT " << sortflen << endl;
   if (sortflen > 0){
#ifndef DEC
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
      
      itbuf[0] = itbuf[1] = -1;
      begi = begj = 0;
      fcnt = 0;
      for (i=0; i < (sortflen/sizeof(int)); i+=2){
         if (itbuf[0] != -1){
            write_middle_range(begi, begj, itbuf[0], itbuf[1],
                               l2cnt, use_seq);
         }         
         if (itbuf[0] != sortary[i] || itbuf[1] != sortary[i+1]){
            if (fcnt >= MINSUPPORT){
               //cout << itbuf[0] << " " << itbuf[1] << " " << fcnt << endl;
               l2cnt++;
            }
            itbuf[0] = sortary[i];
            itbuf[1] = sortary[i+1];
            if (use_seq){
               fcnt = (int) cseq_sup[itbuf[0]][itbuf[1]];
               cseq_sup[itbuf[0]][itbuf[1]] = 0;
            }
            else{
               fcnt = (int) cset_sup[itbuf[0]][itbuf[1]-itbuf[0]-1];
               cset_sup[itbuf[0]][itbuf[1]-itbuf[0]-1] = 0;
            }
         }
         fcnt += 256;
      }
      write_middle_range(begi, begj, itbuf[0], itbuf[1], l2cnt, use_seq);
      if (fcnt >= MINSUPPORT){
         //cout << itbuf[0] << " " << itbuf[1] << " " << fcnt << endl;
         l2cnt++;
      }
      write_middle_range(begi, begj, numfreq-1, numfreq-1, l2cnt, use_seq);
      munmap((caddr_t)sortary, sortflen);
   }
   else{
      begi = begj = 0;
      write_middle_range(begi, begj, numfreq-1, numfreq-1, l2cnt, use_seq);
   }
}

int gen_l2_multpass()
{
   int i;
   double t1,t2;

   int l2cnt=0;
   int mem_used;

   if (use_invert){
      EXTBLKSZ = num_partitions+(DBASE_NUM_TRANS+num_partitions-1)/num_partitions;
      mem_used = (int) (EXTBLKSZ*DBASE_AVG_CUST_SZ*DBASE_AVG_TRANS_SZ*2*ITSZ);
   }
   else mem_used = partition_get_max_blksz();
   int *MAINBUF=NULL;
   if (!use_invert)
      MAINBUF =(int*)malloc(mem_used*ITSZ);
   
   curit = new int *[EXTBLKSZ];
   curcnt = new int [EXTBLKSZ];
   curitsz = new int [EXTBLKSZ];
   curitsz[0] = (int) (DBASE_AVG_CUST_SZ*DBASE_AVG_TRANS_SZ);
   cout << "CURITSZ " << curitsz[0] << " " << EXTBLKSZ << " " << mem_used << endl;
   for (i=0; i < EXTBLKSZ; i++){
      curitsz[i] = curitsz[0];
      curit[i] = (int *) malloc (curitsz[i]*ITSZ);
      curcnt[i] = 0;
   }
   
   
   if (!use_invert){      
      extary = new Extary *[numfreq];

      int extarysz  = AVAILMEM/numfreq;
      extarysz /= ITSZ;
      //cout << "EXTARYSZ " << extarysz << endl;
      if (extarysz < DBASE_AVG_CUST_SZ) extarysz = (int)DBASE_AVG_CUST_SZ;
      for (i=0; i < numfreq; i++){
         extary[i] = new Extary();
      }
   }

   if (use_charary){
      if ((seqfd = open(TEMPSEQ, (O_RDWR|O_CREAT|O_TRUNC), 0666)) < 0){
         perror("Can't open out file");
         exit (errno);      
      }
      if ((isetfd = open(TEMPISET, (O_RDWR|O_CREAT|O_TRUNC), 0666)) < 0){
         perror("Can't open out file");
         exit (errno);      
      }
   }
   if(use_charary){
      cset_sup = new unsigned char *[numfreq];        // support for 2-itemsets
      bzero((char *)cset_sup, numfreq*sizeof(unsigned char *));
      cseq_sup = new unsigned char *[numfreq];        // support for 2-sequences
      bzero((char *)cseq_sup, numfreq*sizeof(unsigned char *));      
   }
   else{
      set_sup = new int *[numfreq];        // support for 2-itemsets
      bzero((char *)set_sup, numfreq*sizeof(int *));
      seq_sup = new int *[numfreq];        // support for 2-sequences
      bzero((char *)seq_sup, numfreq*sizeof(int *));
   }
   int low, high;
      
   seconds(t1);

   int iter = (num_partitions) ? num_partitions:1;
   int mval;
   int itsz = (use_charary) ? sizeof(unsigned char):ITSZ;
   cout << "ITSZ " << itsz << endl;
   for (low = 0; low < numfreq; low = high){
      if (use_charary){
         lseek(seqfd, 0, SEEK_SET);
         lseek(isetfd, 0, SEEK_SET);
      }
      for (high = low; high < numfreq &&
              (mem_used+(2*numfreq-high-1)*itsz) < AVAILMEM; high++){
         if (numfreq-high-1 > 0){
            if (use_charary){
               cset_sup[high] = new unsigned char [numfreq-high-1];
               bzero((char *)cset_sup[high], (numfreq-high-1)*itsz);
            }
            else{
               set_sup[high] = new int [numfreq-high-1];
               bzero((char *)set_sup[high], (numfreq-high-1)*itsz);
            }
         }
         if (use_charary){
            cseq_sup[high] = new unsigned char [numfreq];
            bzero((char *)cseq_sup[high], numfreq*itsz);
         }
         else{
            seq_sup[high] = new int [numfreq];
            bzero((char *)seq_sup[high], numfreq*itsz);
         }
         mem_used += (2*numfreq-high-1) * itsz;
         //cout << "MEMUSEDLOOP " << mem_used << " " << endl;
         //   (int)seq_sup[high] << " " << (int)set_sup[high] << endl;         
      }
      cout << "MEMUSEDLOOP " << mem_used << endl;  
      cout << "LOWHIGH " << high << " " << low << endl;
      fout << "(" << low << " " << high << ") ";
      mval = 0;
      for (int p=0; p < iter; p++){
         if (use_invert){
            process_invert(p,mval);
         }
         else{
            int valid = numfreq;
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
            find_mult_cust(valid, mval);
         }
         
         //cout << "GOT ALL CUST" << i << endl << flush;
         seconds(t2);
         cout << "PHASE " << p << " " << t2-t1 << endl;
         //cout << "POS " << seqpos << " " << itcntpos << endl;
      }
      
      //cout << "POS " << seqpos << " " << itcntpos << endl;
      if (use_charary){
         if (seqpos > 0){
            write(seqfd, (char *)seqbuf, seqpos*sizeof(int));
            seqpos = 0;
         }
         if (itcntpos > 0){
            write(isetfd, (char *)itcntbuf, itcntpos*sizeof(int));
            itcntpos = 0;
         }
         sort_get_l2(l2cnt, isetfd, 0);
         sort_get_l2(l2cnt, seqfd, 1);
      }

      // reclaim memory
      for (i = low; i < high; i++)
      {
         if (use_charary){
            if (cset_sup[i])delete [] cset_sup[i];
            cset_sup[i] = NULL;
            delete [] cseq_sup[i];
            cseq_sup[i] = NULL;
         }
         else{
            //cout << "DELETE " << i << endl << flush;
            get_mult_L2(i, l2cnt);
            if (set_sup[i]) delete [] set_sup[i];
            set_sup[i] = NULL;
            delete [] seq_sup[i];
            seq_sup[i] = NULL;
         }
         mem_used -= (2*numfreq-i-1) * itsz;
      }
   }

   if (use_charary){
      close(isetfd);
      close(seqfd);
      unlink(TEMPSEQ);
      unlink(TEMPISET);
   }
   
   //reclaim memory
   if (MAINBUF) free(MAINBUF);
   if (use_charary){
      delete [] cset_sup;
      delete [] cseq_sup;
   }
   else{
      delete [] set_sup;
      delete [] seq_sup;
   }
   if (!use_invert){
      for (i=0; i < numfreq; i++)
         if (extary[i]){
            extary[i]->set_buf(NULL);
            delete extary[i];
         }
      delete [] extary;
   }
   for (i=0; i < EXTBLKSZ; i++) free(curit[i]);
   delete [] curit;
   delete [] curcnt;
   delete [] curitsz;
   free(backidx);
   
   return l2cnt;
}

void process_apralg_cust(int fcnt, int *fidx, int custid)
{
   int i, j;
   int ind_i, ind_j;

   if (use_charary){
      for (i = 0; i < fcnt; i++){
         ind_i = fidx[i];
         if (ind_i != -1 && cseq_sup[ind_i]){
            // support for itemsets
            for (j = i+1; j < fcnt && fidx[j] != -1; j++){
               if (set_cid[ind_i]){
                  ind_j = fidx[j];
                  if (set_cid[ind_i][ind_j-ind_i-1] != custid) {
                     if ( (++cset_sup[ind_i][ind_j-ind_i-1]) == 0){
                        if (itcntpos+2 > seqitcntbufsz){
                           write(isetfd, (char *)itcntbuf, itcntpos*sizeof(int));
                           itcntpos = 0;
                        }
                        itcntbuf[itcntpos++] = ind_i;
                        itcntbuf[itcntpos++] = ind_j;
                     }
                     set_cid[ind_i][ind_j-ind_i-1] = custid;
                  }
               }
            }
            
            // support for sequences
            for ( ; j < fcnt; j++){
               ind_j = fidx[j];
               if (ind_j != -1){
                  if (seq_cid[ind_i][ind_j] != custid) {
                     if ( (++cseq_sup[ind_i][ind_j]) == 0){
                        if (seqpos+2 > seqitcntbufsz){
                           write(seqfd, (char *)seqbuf, seqpos*sizeof(int));
                           seqpos = 0;
                        }
                        seqbuf[seqpos++] = ind_i;
                        seqbuf[seqpos++] = ind_j;
                     }
                     seq_cid[ind_i][ind_j] = custid;
                  }
               }
            }
         }
      }
   }
   else{
      for (i = 0; i < fcnt; i++){
         ind_i = fidx[i];
         if (ind_i != -1 && seq_sup[ind_i]){
            // support for itemsets
            for (j = i+1; j < fcnt && fidx[j] != -1; j++){
               if (set_cid[ind_i]){
                  ind_j = fidx[j];
                  if (set_cid[ind_i][ind_j-ind_i-1] != custid) {
                     set_sup[ind_i][ind_j-ind_i-1]++;
                     set_cid[ind_i][ind_j-ind_i-1] = custid;
                  }
               }
            }
            
            // support for sequences
            for ( ; j < fcnt; j++){
               ind_j = fidx[j];
               if (ind_j != -1){
                  if (seq_cid[ind_i][ind_j] != custid) {
                     seq_sup[ind_i][ind_j]++;
                     seq_cid[ind_i][ind_j] = custid;
                  }
               }
            }
         }
      }
   }
}


void process_horz_cust(int fcnt, int *fidx, int *ext[2], char *extc, char **ocust){
   int k, j;
   int ii1, ii2;

   if (use_charary){
      for (k=0; k < fcnt; k++){
         for (j=k; j < fcnt; j++){
            if (cseq_sup[fidx[k]] &&
                ext[0][fidx[k]] < ext[1][fidx[j]]){
               if ( (++cseq_sup[fidx[k]][fidx[j]]) == 0){
                  if (seqpos+2 > seqitcntbufsz){
                  write(seqfd, (char *)seqbuf, seqpos*sizeof(int));
                  seqpos = 0;
                  }
                  seqbuf[seqpos++] = fidx[k];
                  seqbuf[seqpos++] = fidx[j];
               }
            }
            if (j > k){
               if (fidx[k] < fidx[j]){
                  ii1 = fidx[k];
                  ii2 = fidx[j];
               }
               else{
                  ii2 = fidx[k];
                  ii1 = fidx[j];
               }
               if (ocust[ii1] && ocust[ii1][ii2-ii1-1] == 1){
                  if ((++cset_sup[ii1][ii2-ii1-1]) == 0){
                     if (itcntpos+2 > seqitcntbufsz){
                        write(isetfd, (char *)itcntbuf, itcntpos*sizeof(int));
                        itcntpos = 0;
                     }
                     itcntbuf[itcntpos++] = ii1;
                     itcntbuf[itcntpos++] = ii2;
                  }
                  ocust[ii1][ii2-ii1-1] = 0;
               }
               if (cseq_sup[fidx[j]] &&
                   ext[0][fidx[j]] < ext[1][fidx[k]]){
                  if ( (++cseq_sup[fidx[j]][fidx[k]]) == 0){
                     if (seqpos+2 > seqitcntbufsz){
                        write(seqfd, (char *)seqbuf, seqpos*sizeof(int));
                        seqpos = 0;
                     }
                     seqbuf[seqpos++] = fidx[j];
                     seqbuf[seqpos++] = fidx[k];
                  }
               }
            }
         }
         extc[fidx[k]] = 0;
      }
   }
   else{
      for (k=0; k < fcnt; k++){
         for (j=k; j < fcnt; j++){
            if (seq_sup[fidx[k]] && ext[0][fidx[k]] < ext[1][fidx[j]]){
               seq_sup[fidx[k]][fidx[j]]++;
            }
            if (j > k){
               if (fidx[k] < fidx[j]){
                  ii1 = fidx[k];
                  ii2 = fidx[j];
               }
               else{
                  ii2 = fidx[k];
                  ii1 = fidx[j];
               }
               if (ocust[ii1] && ocust[ii1][ii2-ii1-1] == 1){
                  set_sup[ii1][ii2-ii1-1]++;
                  ocust[ii1][ii2-ii1-1] = 0;
               }
               if (seq_sup[fidx[j]] && ext[0][fidx[j]] < ext[1][fidx[k]]){
                  seq_sup[fidx[j]][fidx[k]]++;
               }
            }
         }
         extc[fidx[k]] = 0;
      }
   }
}

void process_horz_data(int &NL1, int &NL2)
{
   double t1, t2;
   int i,j,l;
   int custid, tid, numitem;

   sprintf(dataf, "%s.data", input);
   Dbase_Ctrl_Blk *DCB = new Dbase_Ctrl_Blk(dataf);
   
   
   int *itcnt  = new int [DBASE_MAXITEM];
   int *ocnt = new int [DBASE_MAXITEM];
   bzero((char *)itcnt, ((DBASE_MAXITEM)*ITSZ));
   //bzero((char *)ocnt, ((DBASE_MAXITEM)*ITSZ));
   for (i=0; i < DBASE_MAXITEM; i++) ocnt[i] = -1;
   seconds(t1);
   //count 1 items
   int *buf;
   DCB->get_first_blk();
   DCB->get_next_trans(buf, numitem, tid, custid);
   while (!DCB->eof()){
      for (j=0; j < numitem; j++){
         if (ocnt[buf[j]] != custid) itcnt[buf[j]]++;
         ocnt[buf[j]] = custid;
      }
      DCB->get_next_trans(buf, numitem, tid, custid);
   }
   
   int *freqidx = new int [DBASE_MAXITEM];
   numfreq = 0;
   for (i=0; i < DBASE_MAXITEM; i++){
      if(itcnt[i] >= MINSUPPORT){
         freqidx[i] = numfreq;
         numfreq++;
      }
      else freqidx[i] = -1;
   }
   backidx = new int[numfreq];
   j=0;
   for (i=0; i < DBASE_MAXITEM; i++){
      if (itcnt[i] >= MINSUPPORT)
         backidx[j++] = i;
   }
   seconds(t2);
   cout << "NUMFREQ " << numfreq << " :  " << t2-t1 << endl;
   NL1 = numfreq;

   delete [] ocnt;
   delete [] itcnt;

   int mem_used = 0;

   if (use_charary){
      if ((seqfd = open(TEMPSEQ, (O_RDWR|O_CREAT|O_TRUNC), 0666)) < 0){
         perror("Can't open out file");
         exit (errno);      
      }
      if ((isetfd = open(TEMPISET, (O_RDWR|O_CREAT|O_TRUNC), 0666)) < 0){
         perror("Can't open out file");
         exit (errno);      
      }
   }
   int *ext[2];
   char **ocust, *extc;
   if (!use_apralg){
      ext[0] = new int [numfreq];
      ext[1] = new int [numfreq];
      extc = new char[numfreq];
      bzero(extc, numfreq);
      mem_used += numfreq+2*numfreq*ITSZ;
   
      ocust  = new char *[numfreq];
      bzero((char *)ocust, numfreq*sizeof(char *));
   }
   
   if(use_charary){
      cset_sup = new unsigned char *[numfreq];        // support for 2-itemsets
      bzero((char *)cset_sup, numfreq*sizeof(unsigned char *));
      cseq_sup = new unsigned char *[numfreq];        // support for 2-sequences
      bzero((char *)cseq_sup, numfreq*sizeof(unsigned char *));      
   }
   else{
      set_sup = new int *[numfreq];        // support for 2-itemsets
      bzero((char *)set_sup, numfreq*sizeof(int *));
      seq_sup = new int *[numfreq];        // support for 2-sequences
      bzero((char *)seq_sup, numfreq*sizeof(int *));
   }
   if (use_apralg){
      set_cid = new int *[numfreq];
      bzero((char *)set_cid, numfreq*sizeof(int *));
      seq_cid = new int *[numfreq];
      bzero((char *)seq_cid, numfreq*sizeof(int *));
   }
   
   int low, high;
   int *fidx = (int *) malloc (numfreq*ITSZ);
   int fidxsz = numfreq;
   
   int itsz = (use_charary) ? sizeof(unsigned char):ITSZ;
   seconds(t1);
   int l2cnt = 0;
   for (low = 0; low < numfreq; low = high){
      if (use_charary){
         seqpos = 0;
         itcntpos = 0;
         lseek(seqfd, 0, SEEK_SET);
         lseek(isetfd, 0, SEEK_SET);
      }
      for (high = low; high < numfreq && mem_used < AVAILMEM; high++){
         int newamt = mem_used + (2*numfreq-high-1)*itsz;
         if (use_apralg) newamt += (2*numfreq-high-1)*ITSZ;
         else newamt += numfreq-high-1*sizeof(char);
         if (newamt > AVAILMEM) break;
         
         if (numfreq-high-1 > 0){
            if (!use_apralg){
               ocust[high] = new char [numfreq-high-1];
               bzero((char *)ocust[high], (numfreq-high-1));
            }
            if (use_charary){
               cset_sup[high] = new unsigned char [numfreq-high-1];
               bzero((char *)cset_sup[high], (numfreq-high-1)*itsz);
            }
            else{
               set_sup[high] = new int [numfreq-high-1];
               bzero((char *)set_sup[high], (numfreq-high-1)*itsz);
            }
            if (use_apralg){
               set_cid[high] = new int [numfreq-high-1];
               //bzero((char *)set_cid[high], (numfreq-high-1)*itsz);
               for (i=0; i < numfreq-high-1; i++) set_cid[high][i] = -1;
            }
         }
         if (use_charary){
            cseq_sup[high] = new unsigned char [numfreq];
            bzero((char *)cseq_sup[high], numfreq*itsz);
         }
         else{
            seq_sup[high] = new int [numfreq];
            bzero((char *)seq_sup[high], numfreq*itsz);
         }
         if (use_apralg){
            seq_cid[high] = new int [numfreq];
            //bzero((char *)seq_cid[high], numfreq*itsz);
            for (i=0; i < numfreq; i++) seq_cid[high][i] = -1;
         }
         mem_used = newamt;
      }
      cout << "MEMUSEDLOOP " << mem_used << endl;  
      cout << "LOWHIGH " << high << " " << low << endl;
      fout << "(" << low << " " << high << ") ";
      //cout << "POS " << seqpos << " " << itcntpos << endl;

      int idx;
      int fcnt=0;
      int ocid = -1;
      //count 2-itemsets
      DCB->get_first_blk();
      DCB->get_next_trans(buf, numitem, tid, custid);
      while(!DCB->eof()){
         fcnt=0;
         ocid = custid;
         while(!DCB->eof() && ocid == custid){
            for (j=0; j < numitem; j++){
               idx = freqidx[buf[j]];
               if (idx != -1){
                  if (use_apralg){
                     if (fcnt+1 > fidxsz){
                        fidxsz = (int) (1.5*fidxsz);
                        fidx = (int *)realloc(fidx, fidxsz*ITSZ);
                        if (fidx == NULL){
                           perror("realloc fidx");
                           exit(errno);
                        }
                     }
                     fidx[fcnt++] = idx;
                  }
                  else{
                     if (extc[idx] == 0){
                        extc[idx] = 1;
                        fidx[fcnt] = idx;
                        fcnt++;
                        ext[0][idx] = tid;
                     }
                     ext[1][idx] = tid;
                     
                     if (ocust[idx]){
                        for (l=j+1; l < numitem; l++){
                           if (freqidx[buf[l]] != -1){
                              ocust[idx][freqidx[buf[l]]-idx-1] = 1;
                           }
                        }
                     }
                  }
               }
            }
            if (use_apralg){
               if (fcnt+1 > fidxsz){
                  fidxsz = (int) (1.5*fidxsz);
                  fidx = (int *)realloc(fidx, fidxsz*ITSZ);
                  if (fidx == NULL){
                     perror("realloc fidx");
                     exit(errno);
                  }
               }
               fidx[fcnt++] = -1;
            }
            DCB->get_next_trans(buf, numitem, tid, custid);
         }
         if (use_apralg) process_apralg_cust(fcnt, fidx, ocid);
         else process_horz_cust(fcnt, fidx, ext, extc, ocust);
      }
      seconds(t2);
      cout << "2-IT " << t2-t1 << " " << endl;
      
      //cout << "POS " << seqpos << " " << itcntpos << endl;      
      if (use_charary){
         if (seqpos > 0){
            write(seqfd, (char *)seqbuf, seqpos*sizeof(int));
            seqpos = 0;
         }
         if (itcntpos > 0){
            write(isetfd, (char *)itcntbuf, itcntpos*sizeof(int));
            itcntpos = 0;
         }
         sort_get_l2(l2cnt, isetfd, 0);
         sort_get_l2(l2cnt, seqfd, 1);
      }

      // reclaim memory
      for (i = low; i < high; i++)
      {
         if (!use_apralg){
            if (ocust[i]) delete [] ocust[i];
            ocust[i] = NULL;
            mem_used -= (numfreq-i-1) * sizeof(char);
         }
         if (use_charary){
            //write_charary(i);
            if (cset_sup[i]) delete [] cset_sup[i];
            cset_sup[i] = NULL;
            delete [] cseq_sup[i];
            cseq_sup[i] = NULL;
         }
         else{
            //cout << "DELETE " << i << endl << flush;
            get_mult_L2(i, l2cnt);
            if (set_sup[i]) delete [] set_sup[i];
            set_sup[i] = NULL;
            delete [] seq_sup[i];
            seq_sup[i] = NULL;
         }
         if (use_apralg){
            if (set_cid[i]) delete [] set_cid[i];
            set_cid[i] = NULL;
            delete [] seq_cid[i];
            seq_cid[i] = NULL;
            mem_used -= (2*numfreq-i-1)*ITSZ;
         }
         mem_used -= (2*numfreq-i-1) * itsz;
      }
   }
   NL2 = l2cnt;
   if (use_charary){
      close(isetfd);
      close(seqfd);
      unlink(TEMPSEQ);
      unlink(TEMPISET);
   }
   
   //reclaim memory
   if (use_charary){
      delete [] cset_sup;
      delete [] cseq_sup;
   }
   else{
      delete [] set_sup;
      delete [] seq_sup;
   }
   if (use_apralg){
      delete [] set_cid;
      delete [] seq_cid;
   }
   if (!use_apralg){
      delete [] ext[0];
      delete [] ext[1];
      delete [] extc;
   }
   delete [] backidx;
   delete [] freqidx;
   free(fidx);
   delete DCB;
   seconds(t2);
   L2TIME = t2-t1;
}


int main(int argc, char **argv)
{
   double t1, t2;
   seconds(t1);
   parse_args(argc, argv);

   fout.open("summary.out", ios::app);
   fout << "CALCL2 ";
   if (new_format) fout << "NEWFORMAT ";
   if (use_charary) fout << "CHARARY ";
   if (use_invert) fout << "USEINVERT ";
   if (use_horz) fout << "HORZ ";
   if (use_apralg) fout << "APRALG ";
   fout << dataf << " " << MINSUP_PER << " ";
   int NL1, NL2;
   if (use_horz){
      process_horz_data(NL1, NL2);
   }
   else{
      partition_alloc(dataf, idxf);
      NL1 = make_l1_pass();
      seconds(L2TIME);
      NL2 = gen_l2_multpass();
      seconds(t2);
      L2TIME = t2-L2TIME;
      partition_dealloc();
   }
   seconds(t2);
   cout << "END " << NL1 << " " << NL2 << " " << t2-t1 << " " << L2TIME << endl;

   fout << NL1 << " " << NL2 << " " << t2-t1 << " " << L2TIME << endl;
   fout.close();
}

