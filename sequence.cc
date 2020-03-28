#include <errno.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
//#include <sys/mode.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <sys/types.h>
#include <strings.h>
#include <sys/mman.h>
#include <malloc.h>
#include <strings.h>
#include "Eqclass.h"
#include "Itemset.h"
#include "Lists.h"
#include "extl2.h"
#include "partition.h"

#define min(a, b) ((a) < (b) ? (a) : (b))


#define NONMAXFLG -2
#define MBYTE (1024*1024)

//#define GETBIT(a,b) ((a >> b) & 01)
//#define SETBIT(a,v,b)  ((v != 0) ? (a | (01 << b)): (a & ~(01 << b)))

#define CUSTID(a, b) ((a)[2*(b)])
#define TID(a, b) ((a)[2*(b)+1])
#define LJOIN 0
#define EJOIN 1
#define MJOIN 2

long MEMUSED = 0;
long AVAILMEM = 32*MBYTE;
char fname[300];
char constf[300];               // constraint file
char dataf[300];
char idxf[300];
char conf[300];
char it2f[300];
char seqf[300];

double L2ISECTTIME=0, EXTL1TIME=0, EXTL2TIME=0;

int *NumLargeItemset;
int maxitemsup;
Itemset *item1, *item2; // for use in reading external dbase
Array *interval, *interval2, *interval3;
int maxeqsize = 1;
EqGrNode** eqgraph;
double MAXIMAL_THRESHOLD = 1.0;
int ext_l2_pass = 0;
int use_clique = 0;
char use_constraint =0;
int use_hash = 0;
int num_intersect=0;
int recursive = 0;
FILE *out;
int maxiter = 2;
int min_gap = 0;
int max_gap = (int)pow(2,sizeof(int)*8); //infinity
int window = 0;
char memtrace = 0;
char use_newformat = 1;
int use_ascending = -2;
char use_isetonly = 0;
int L2pruning = 0;
char print_seq = 0;

ofstream mout;
int DBASE_NUM_TRANS;
int DBASE_MAXITEM;
float DBASE_AVG_TRANS_SZ;
float DBASE_AVG_CUST_SZ; 
int DBASE_TOT_TRANS;

double MINSUP_PER;
int MINSUPPORT=-1;


int FreqArraySz = 100;
FreqIt **FreqArray;
int FreqArrayPos = 0;

void process_cluster1(Eqclass *cluster, Lists<Eqclass *> *LargeL, int iter);

void add_freq(Itemset *it, int templ)
{
   FreqIt *freq = new FreqIt(it->itemset()->array(), it->size(), templ);
   if (FreqArrayPos+1 >= FreqArraySz){
      FreqArraySz = (int)(1.5*FreqArraySz);
      FreqArray = (FreqIt **)realloc(FreqArray, FreqArraySz*sizeof(FreqIt*));
      if (FreqArray == NULL){
         perror("no mmeory fro FREqArray ");
         exit(-1);
      }
   }
   FreqArray[FreqArrayPos++] = freq;
}

void print_freqary()
{
   int j=0;
   cout << "FREQARRAY " << FreqArrayPos << ":" << endl;
   for (j=0; j < FreqArrayPos; j++){
      cout << *FreqArray[j];
   }
   cout << "!!!!!!!!!!!!!!!!!!!!" << endl;
}

void parse_args(int argc, char **argv)
{
   extern char * optarg;
   int c;
   
   if (argc < 2)
      cout << "usage: seq -i<infile> -s<support> -e 1\n";
   else{
      while ((c=getopt(argc,argv,"a:bce:fhi:l:m:ors:t:u:v:w:x:"))!=-1){
         switch(c){
         case 'a':
            //if val = -1 then do ascending generation
            //else only generate the eqclass given by the value
            use_ascending = atoi(optarg);
            break;
         case 'b':
            use_isetonly = 1;
            break;
         case 'c': //use constraints file
            use_constraint = 1;
            sprintf(constf,"%s.constraint", fname);            
            break;
         case 'e': //calculate L2 from inverted dbase
            num_partitions = atoi(optarg);
            ext_l2_pass = 1;
            break;
         case 'f':
            use_newformat = 0;
            break;
         case 'h': //use hashing to prune candidates
            use_hash = 1;
            break;
         case 'i': //input file
            sprintf(fname, "%s", optarg);
            sprintf(dataf,"%s.tpose", optarg);
            sprintf(idxf,"%s.idx", optarg);
            sprintf(conf,"%s.conf", optarg);
            sprintf(it2f,"%s.2it", optarg);
            sprintf(seqf,"%s.2seq", optarg);
            break;
         case 'l': //min-gap between items (not implemented)
            min_gap = atoi(optarg);
            break;
         case 'm': //amount of mem available
            AVAILMEM = (long) atof(optarg)*MBYTE;
            break;
         case 'o': //print sequences
            print_seq = 1;
            break;
         case 'r': //use recursive algorithm (doesn't work with pruning)
            recursive = 1;
            break;
         case 's': //min support
            MINSUP_PER = atof(optarg);
            break;
         case 't': //not used
            MAXIMAL_THRESHOLD = atof(optarg);
            break;
         case 'u': //max-gap between items (not implemented)
            max_gap = atoi(optarg);
            break;
         case 'v':
            MINSUPPORT = atoi(optarg);
            break;
         case 'w': //sliding window (not implemented)
            window = atoi(optarg);
            break;
         case 'x':
            memtrace = 1;
            mout.open(optarg, ios::app);
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
   if (MINSUPPORT == -1)
      MINSUPPORT = (int) (MINSUP_PER*DBASE_NUM_TRANS+0.5);
   //ensure that support is at least 2
   if (MINSUPPORT < 2) MINSUPPORT = 2;
   //cout << "MINSUPPORT " << MINSUPPORT << " " << DBASE_NUM_TRANS << endl;
   read(c,(char *)&DBASE_MAXITEM,ITSZ);
   read(c,(char *)&DBASE_AVG_CUST_SZ,sizeof(float));
   read(c,(char *)&DBASE_AVG_TRANS_SZ,sizeof(float));
   read(c,(char *)&DBASE_TOT_TRANS,ITSZ);
   cout << "CONF " << DBASE_NUM_TRANS << " " << DBASE_MAXITEM << " "
        << DBASE_AVG_CUST_SZ << " " << DBASE_AVG_TRANS_SZ << " "
        << DBASE_TOT_TRANS << endl;
   close(c);
}

int choose(int n, int k)
{
   int i;
   int val = 1;

   if (k >= 0 && k <= n){
      for (i=n; i > n-k; i--)
         val *= i;
      for (i=2; i <= k; i++)
         val /= i;
   }
   
   return val;
}

void get_2_intersect(Itemset *ljoin, Itemset *ejoin,
                     int *it1, int *it2, int sup1, int sup2)
{
   int i,j,k,l;
   int nval1, nval2;
   int lflge;
   int olpos, oepos;
   
   num_intersect++;

   int itid, jtid;
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
         //cout << "comm " << i << " " << j << " " << nval1
         //     << " " << nval2 <<endl;
         if (ljoin){
            //lflgl = 0;
            olpos = ljoin->ivalsize();
            if (it1[i] < it2[j+nval2-1]){
               //add tid
               ljoin->add_ival(itid);
               //add dummy value for size
               ljoin->add_ival(0);
               
               for (l=j; l < j+nval2;l++){
                  if (it1[i] < it2[l])
                     ljoin->add_ival(it2[l]);
                  //lflgl=1;
               }
               ljoin->increment_support();
               ljoin->ival()->setitem(olpos+1, ljoin->ivalsize()-olpos-2);
            }
            //else lary->set_size(olpos);
         }
         if (ejoin){
            lflge = 0;
            oepos = ejoin->ivalsize();
            //add tid
            ejoin->add_ival(itid);
            //add dummy value for size
            ejoin->add_ival(0);
            for (k=i, l=j; k<i+nval1 && l < j+nval2;){
               if (it1[k] < it2[l]) k++;
               else if (it1[k] > it2[l]) l++;
               else{
                  ejoin->add_ival(it2[l]);
                  lflge=1;
                  k++;
                  l++;
               }
            }
            if (lflge){
               ejoin->increment_support();
               ejoin->ival()->setitem(oepos+1, ejoin->ivalsize()-oepos-2);
            }
            else ejoin->ival()->set_size(oepos);
         }
         i += nval1;
         j += nval2;
      }
   }
}


void get_2newf_intersect(Itemset *ljoin, Itemset *ejoin,
                         int *it1, int *it2, int sup1, int sup2)
{
   if (!use_newformat){
      get_2_intersect(ljoin, ejoin, it1, it2, sup1, sup2);
      return;
   }
   int i,j,k,l;
   int nval1, nval2;
   int lflge;
   
   num_intersect++;

   int itid, jtid;
   for (i=0,j=0; i < sup1 && j < sup2;){
      itid = it1[i];
      jtid = it2[j];
      if (itid > jtid){
         j += 2;
      }
      else if (itid < jtid){
         i += 2;
      }      
      else{
         nval1 = i;
         nval2 = j;
         while(it1[i] == it1[nval1] && nval1 < sup1) nval1 += 2;
         while(it2[j] == it2[nval2] && nval2 < sup2) nval2 += 2;
         //cout << "comm " << i << " " << j << " " << nval1
         //     << " " << nval2 <<endl;
         if (ljoin){
            if (it1[i+1] < it2[nval2-1]){
               //add tid
               
               for (l=j; l < nval2;l+=2){
                  if (it1[i+1] < it2[l+1]){
                     ljoin->ival()->optadd(itid);
                     ljoin->ival()->optadd(it2[l+1]);
                  }
                  //lflgl=1;
               }
               ljoin->increment_support();
            }
            //else lary->set_size(olpos);
         }
         if (ejoin){
            lflge = 0;
            for (k=i, l=j; k< nval1 && l < nval2;){
               if (it1[k+1] < it2[l+1]) k+=2;
               else if (it1[k+1] > it2[l+1]) l+=2;
               else{
                  ejoin->ival()->optadd(itid);
                  ejoin->ival()->optadd(it2[l+1]);
                  lflge=1;
                  k+=2;
                  l+=2;
               }
            }
            if (lflge){
               ejoin->increment_support();
            }
         }
         i = nval1;
         j = nval2;
      }
   }
}

int interval_comp(int a, int b, int c, int d)
{
   if (a < c) return -1;
   else if (a > c) return 1;
   else{
      if (b < d) return -1;
      else if (b > d) return 1;
      else return 0;
   }
}

void make_itemset(Itemset *it, Array *ary, int cnt)
{
   int i;
   for (i=0; i < ary->size(); i++){
      it->ival()->optadd((*ary)[i]);
   }
   it->set_support(cnt);
}

void get_tmp_intersect(Itemset *&ljoin, Itemset *&ejoin, Itemset *&mjoin,
                       int &lcnt, int &ecnt, int &mcnt,
                       Itemset *it1, Itemset *it2, int iter)
{
   int i,j,k,l;
   int nval1, nval2;
   int lflge;
   Array *lary, *eary, *mary;
   int olpos, oepos, ompos;
   
   num_intersect++;

   lary = interval;
   eary = interval2;
   mary = interval3;
   lary->reset();
   eary->reset();
   mary->reset();

   lcnt = ecnt = mcnt = 0;
   
   int dc1 = it1->support()-MINSUPPORT;
   int dc2 = it2->support()-MINSUPPORT;
   int df1=0;
   int df2=0;
   int itid, jtid;
   for (i=0,j=0; i < it1->ivalsize() && j < it2->ivalsize();){
      if (df1 > dc1 || df2 > dc2) break;
      itid = it1->ival(i);
      jtid = it2->ival(j);
      if (itid > jtid){
         nval2 = it2->ival(j+1);
         j += (nval2+2);
         df2++;
      }
      else if (itid < jtid){
         nval1 = it1->ival(i+1);
         i += (nval1+2);
         df1++;
      }      
      else{
         i++;
         j++;
         nval1 = it1->ival(i++);
         nval2 = it2->ival(j++);
         //cout << "comm " << it1->tid(i) << " " << nval1
         //     << " " << nval2 <<endl;
         if (ljoin){
            //lflgl = 0;
            olpos = lary->size();
            if (it1->ival(i) < it2->ival(j+nval2-1)){
               //add tid
               lary->add(itid);
               //add dummy value for size
               lary->add(0);
               for (l=j; l < j+nval2;l++){
                  if (it1->ival(i) < it2->ival(l)){
                     lary->add(it2->ival(l));
                     //lflgl=1;
                  }
               }
               lcnt++;
               lary->setitem(olpos+1, lary->size()-olpos-2);
            }
            //else lary->set_size(olpos);
         }
         if (ejoin){
            lflge = 0;
            oepos = eary->size();
            //add tid
            eary->add(itid);
            //add dummy value for size
            eary->add(0);
            for (k=i, l=j; k<i+nval1 && l < j+nval2;){
               if (it1->ival(k) < it2->ival(l)) k++;
               else if (it1->ival(k) > it2->ival(l)) l++;
               else{
                  eary->add(it2->ival(l));
                  lflge=1;
                  k++;
                  l++;
               }
            }
            if (lflge){
               ecnt++;
               eary->setitem(oepos+1, eary->size()-oepos-2);
            }
            else eary->set_size(oepos);
         }
         if (mjoin){
            //lflgm = 0;
            ompos = mary->size();
            if (it2->ival(j) < it1->ival(i+nval1-1)){
               //add tid
               mary->add(itid);
               //add dummy value for size
               mary->add(0);            
               for (l=i; l < i+nval1;l++){
                  if (it2->ival(j) < it1->ival(l)){
                     mary->add(it1->ival(l));
                     //lflgm=1;
                  }
               }
               mcnt++;
               mary->setitem(ompos+1, mary->size()-ompos-2);
            }
            //else mary->set_size(ompos);
         }
         i += nval1;
         j += nval2;
      }
   }
   if (ljoin && lcnt >= MINSUPPORT){
      ljoin = new Itemset(iter, lary->size());
      make_itemset(ljoin, lary, lcnt);
   }
   if (ejoin && ecnt >= MINSUPPORT){
      ejoin = new Itemset(iter, eary->size());
      make_itemset(ejoin, eary, ecnt);
   }
   if (mjoin && mcnt >= MINSUPPORT){
      mjoin = new Itemset(iter, mary->size());
      make_itemset(mjoin, mary, mcnt);
   }
}

void get_tmpnewf_intersect(Itemset *&ljoin, Itemset *&ejoin, Itemset *&mjoin,
                       int &lcnt, int &ecnt, int &mcnt,
                       Itemset *it1, Itemset *it2, int iter)
{
   if (!use_newformat){
      get_tmp_intersect(ljoin, ejoin, mjoin, lcnt, ecnt, mcnt, it1, it2, iter);
      return;
   }
   int i,j,k,l;
   int nval1, nval2;
   int lflge;
   Array *lary, *eary, *mary;
   
   num_intersect++;

   lary = interval;
   eary = interval2;
   mary = interval3;
   lary->reset();
   eary->reset();
   mary->reset();

   lcnt = ecnt = mcnt = 0;
   
   int dc1 = it1->support()-MINSUPPORT;
   int dc2 = it2->support()-MINSUPPORT;
   int df1=0;
   int df2=0;
   int itid, jtid;
   for (i=0,j=0; i < it1->ivalsize() && j < it2->ivalsize();){
      if (df1 > dc1 || df2 > dc2) break;
      itid = it1->ival(i);
      jtid = it2->ival(j);
      if (itid > jtid){
         //df must be incremented only once per customer
         while(jtid == it2->ival(j) && j < it2->ivalsize()) j += 2;
         df2++;
      }
      else if (itid < jtid){
         while(itid == it1->ival(i) && i < it1->ivalsize()) i += 2;
         df1++;
      }      
      else{
         nval1 = i;
         nval2 = j;
         while(it1->ival(i) == it1->ival(nval1) && nval1 < it1->ivalsize())
            nval1 += 2;
         while(it2->ival(j) == it2->ival(nval2) && nval2 < it2->ivalsize())
            nval2 += 2;

         if (ljoin){
            if (it1->ival(i+1) < it2->ival(nval2-1)){
               for (l=j; l < nval2;l+=2){
                  if (it1->ival(i+1) < it2->ival(l+1)){
                     lary->optadd(itid);
                     lary->optadd(it2->ival(l+1));
                     //lflgl=1;
                  }
               }
               lcnt++;
            }
         }
         if (ejoin){
            lflge = 0;
            for (k=i, l=j; k< nval1 && l < nval2;){
               if (it1->ival(k+1) < it2->ival(l+1)) k+=2;
               else if (it1->ival(k+1) > it2->ival(l+1)) l+=2;
               else{
                  eary->optadd(itid);
                  eary->optadd(it2->ival(l+1));
                  lflge=1;
                  k+=2;
                  l+=2;
               }
            }
            if (lflge){
               ecnt++;
            }
         }
         if (mjoin){
            if (it2->ival(j+1) < it1->ival(nval1-1)){
               for (l=i; l < nval1;l+=2){
                  if (it2->ival(j+1) < it1->ival(l+1)){
                     mary->optadd(itid);
                     mary->optadd(it1->ival(l+1));
                  }
               }
               mcnt++;
            }
         }
         i = nval1;
         j = nval2;
      }
   }
   if (ljoin && lcnt >= MINSUPPORT){
      ljoin = new Itemset(iter, lary->size());
      make_itemset(ljoin, lary, lcnt);
   }
   if (ejoin && ecnt >= MINSUPPORT){
      ejoin = new Itemset(iter, eary->size());
      make_itemset(ejoin, eary, ecnt);
   }
   if (mjoin && mcnt >= MINSUPPORT){
      mjoin = new Itemset(iter, mary->size());
      make_itemset(mjoin, mary, mcnt);
   }
}


void fill_seq_template(Eqclass *EQ, Eqclass *parent, int LR)
{
   if (LR == 1){
      EQ->set_templ(parent->templ()*2+1);
      EQ->set_templ2(parent->templ()*2);
   }
   else if (LR == 2){
      EQ->set_templ(parent->templ2()*2+1);
      EQ->set_templ2(parent->templ2()*2);
   }
}

int get_valid_el(int it, char *ibvec, char *sbvec)
{
   int i, j;
   int i1, i2;
   int rval = 0;

   // cout << "[" << it << "] = ";
   for (i=0; i < eqgraph[it]->seqnum_elements(); i++){
      sbvec[i] = 0;
      //cout << " " << eqgraph[it]->seqget_element(i);
   }
   //cout << " --- ";
   for (i=0; i < eqgraph[it]->num_elements(); i++){
      ibvec[i] = 0;
      //cout << " " << eqgraph[it]->get_element(i);
   }
   //cout << endl;
   
   for (i=0; i < eqgraph[it]->seqnum_elements(); i++){
      i1 = eqgraph[it]->seqget_element(i);
      for (j=i; j < eqgraph[it]->seqnum_elements(); j++){
         i2 = eqgraph[it]->seqget_element(j);
         if (eqgraph[i1] && eqgraph[i1]->seqfind(i2)){
            sbvec[i] = 1;
            sbvec[j] = 1;
            rval = 1;
         }
         if (j > i){
            if ((eqgraph[i1] && eqgraph[i1]->find(i2)) ||
                (eqgraph[i2] && eqgraph[i2]->seqfind(i1))){
               sbvec[i] = 1;
               sbvec[j] = 1;
               rval = 1;
            }
         }
      }
   }
   

   for (i=0; i < eqgraph[it]->num_elements(); i++){
      i1 = eqgraph[it]->get_element(i);
      if (eqgraph[i1]){
         for (j=i+1; j < eqgraph[it]->num_elements(); j++){
            i2 = eqgraph[it]->get_element(j);
            if (eqgraph[i1]->find(i2)){
               ibvec[i] = 1;
               ibvec[j] = 1;
               rval = 1;
            }
         }
         for (j=0; j < eqgraph[it]->seqnum_elements(); j++){
            i2 = eqgraph[it]->seqget_element(j);
            if (eqgraph[i1]->seqfind(i2)){
               ibvec[i] = 1;
               sbvec[j] = 1;
               rval =1;
            }
         }      
      }
   }
   
   for (i=0; i < eqgraph[it]->seqnum_elements(); i++)
      if (!sbvec[i]) L2pruning++;
   for (i=0; i < eqgraph[it]->num_elements(); i++)
      if (!ibvec[i]) L2pruning++;

   return rval;
}

//construct the next set of eqclasses from external disk
Eqclass* get_ext_eqclass(int it)
{
   double t1, t2;
   seconds(t1);
   //cout << "MEMEXT " << it << " " << MEMUSED << endl;
   int i, k, it2, supsz, supsz2;
   Itemset *ljoin = NULL;
   Itemset *ejoin = NULL;

   char *ibvec, *sbvec;
   ibvec = sbvec = NULL;
   if (eqgraph[it]->num_elements() > 0)
      ibvec = new char[eqgraph[it]->num_elements()];
   if (eqgraph[it]->seqnum_elements() > 0)
      sbvec = new char[eqgraph[it]->seqnum_elements()];

   if (!get_valid_el(it, ibvec, sbvec)) return NULL;
   
   Eqclass *L2 = new Eqclass(1, EQCTYP1);
   if (L2 == NULL)
   {
      perror("memory exceeded : ext_class ");
      exit (errno);
   }
   //init seq pattern templates
   L2->set_templ(1);
   L2->set_templ2(0);
   
   interval->reset();
   interval2->reset();
   
   supsz = partition_get_idxsup(it);
   partition_read_item(interval->array(), it);

   int tmpit;
   for (i=0, k=0; i < eqgraph[it]->num_elements() ||
           k < eqgraph[it]->seqnum_elements();){
      ljoin = NULL;
      ejoin = NULL;
      it2 = DBASE_MAXITEM+1;
      tmpit = DBASE_MAXITEM+1;
      if (i < eqgraph[it]->num_elements() && ibvec[i])
         it2 = eqgraph[it]->get_element(i);
      if (k < eqgraph[it]->seqnum_elements() && sbvec[k])
         tmpit = eqgraph[it]->seqget_element(k);
      if (it2 == tmpit){
         ejoin=(Itemset*)1;
         ljoin=(Itemset*)1;
         k++;
         i++;
         if (it2 == DBASE_MAXITEM+1) continue;
      }
      else if (it2 < tmpit){
         ejoin=(Itemset*)1;
         i++;
      }
      else{
         ljoin=(Itemset*)1;
         k++;
         it2 = tmpit;
      }
      //cout << "JOIN " << it << " " << it2 << " " << ejoin << " " << ljoin << endl << flush;
      supsz2 = partition_get_idxsup(it2);
      
      partition_read_item(interval2->array(), it2);
      
      if (ejoin){
         ejoin = new Itemset(2, min(supsz, supsz2));
         if (ejoin == NULL){
            perror("memory exceeded");
            exit(errno);
         }
      }
      else ejoin = NULL;
      if (ljoin){
         ljoin = new Itemset(2, supsz2);
         if (ljoin == NULL){
            perror("memory exceeded");
            exit(errno);
         }
      }
      else ljoin = NULL;
      //cout << "ljoin " << ljoin << " " << ejoin << " " <<
      //supsz << " " << supsz2 << " " << it << " " << it2 << endl;
      
      get_2newf_intersect(ljoin, ejoin, interval->array(), interval2->array(),
                      supsz, supsz2);
      
      if (ljoin){
         if (ljoin->support() >= MINSUPPORT && !use_isetonly){
            ljoin->reallocival();
            ljoin->add_item(it);
            ljoin->add_item(it2);
            L2->prepend(ljoin);
            //cout << "LARGE ";
            //ljoin->print_seq(L2->templ());
            //NumLargeItemset[1]++;
         }
         else{
            //cout << "DELETED ";
            //ljoin->print_seq(L2->templ());
            delete ljoin;
         }
      }
      if (ejoin){
         if (ejoin->support() >= MINSUPPORT){
            ejoin->reallocival();
            ejoin->add_item(it);
            ejoin->add_item(it2);
            L2->prepend2(ejoin);
            //cout << "LARGE ";
            //ejoin->print_seq(L2->templ2());
            //NumLargeItemset[1]++;
         }
         else{
            //cout << "DELETED ";
            //ejoin->print_seq(L2->templ2());
            delete ejoin;
         }
      }
   }

   //cout << "MEMEXTEND " << it << " " << MEMUSED << endl;
   seconds(t2);
   L2ISECTTIME += t2-t1;
   return L2;
}

void delete_eq_list(Lists<Eqclass *> *eqlist)
{
   ListNodes<Eqclass *> *eqhd = eqlist->head();

   for (;eqhd; eqhd=eqhd->next()){
      delete eqhd->item()->list();
      eqhd->item()->set_list(NULL);
      delete eqhd->item();
   }
   delete eqlist;
}

void fill_join(Itemset *join, Itemset *hdr1, Itemset *hdr2)
{
   int i;

   for (i=0; i < hdr1->size(); i++){
      join->add_item((*hdr1)[i]);
   }
   join->add_item((*hdr2)[hdr2->size()-1]);
}

Itemset *prune_decision(Itemset *it1, Itemset *it2, int ptempl, int jflg, int LR)
{
   int l1 = (*it1)[it1->size()-1];
   int l2 = (*it2)[it2->size()-1];
   if (use_hash && (it1->size() > 2)){
      int i,j,k;
      unsigned int bit, ttpl;
      int nsz;
      FreqIt fit(it1->size(), 0);

      //skip the first two subsets (or omit the last two elements)
      nsz = it1->size()-2;
      
      //cout << "PTEMPL " << ptempl << " " << jflg << endl;
         //cout << *it1;
      //cout << *it2;
      for (i=nsz; i >= 0; i--){
         k=0;
         //form new subset template
         if (i == 0) ttpl = SETBIT(ptempl,0,nsz+1);
         else{
            ttpl = 0;
            for (j=0; j < i-1; j++){
               bit = GETBIT(ptempl,nsz-j+1);
               ttpl = SETBIT(ttpl, bit, nsz-j);
            }
            bit = GETBIT(ptempl, nsz-j+1);
            bit = bit || GETBIT(ptempl, nsz-j);
            ttpl = SETBIT(ttpl, bit, nsz-j);
            j+=2;
            for (; j < nsz+2; j++){
               bit = GETBIT(ptempl,nsz-j+1);
               ttpl = SETBIT(ttpl, bit, nsz-j+1);
            }
         }
         //form new subset by omitting the i-th item
         for (j=0; j < nsz+1; j++){
            if (j != i){
               fit.seq[k++] = (*it1)[j];
            }
         }
         fit.seq[k++] = l1;
         fit.seq[k++] = l2;
         fit.templ = ttpl;

         //if ((*it1)[0] == 101)
         //   cout << "SEARCH " << fit;
         if (fit.seq[0] == (*it1)[0] && !recursive){
            //elements should be in current class
            if (FreqArrayPos > 0){
               if (!EqGrNode::bsearch(0,FreqArrayPos-1,FreqArray,
                                      fit, recursive)){
                  //if ((*it1)[0] == 101) cout << "NOTF" <<endl;
                  //print_freqary();
                  return NULL;
               }
               //if ((*it1)[0] == 101) cout << "FOUND" <<endl;
            }
            else return NULL;
         }
         else if (fit.seq[0] > (*it1)[0]){
            // class must already have been processed, otherwise we can't prune
            if (!eqgraph[fit.seq[0]]->find_freqarray(fit, recursive)){
               //if ((*it1)[0] == 101) cout << "NOTF" <<endl;
               return NULL;
            }
            //if ((*it1)[0] == 101) cout << "FOUND" <<endl;
         }
      }
   }
   else if (it1->size() == 2){
      if (eqgraph[l1]){
         if (jflg == LJOIN || jflg == MJOIN){
            if (!eqgraph[l1]->seqfind(l2))
               return NULL;
         }
         else{
            if (!eqgraph[l1]->find(l2))
               return NULL;
         }
      }
      else return NULL;
      //cout << "FOUND " << endl;
   }
   return (Itemset *)1;
}



void insert_freqarray(Lists<Eqclass *> *LargeL)
{
   //insert frequent itemsets into hash table
   ListNodes<Eqclass *> *chd;
   ListNodes<Itemset *> *hdr1, *hdr2;
   Eqclass *cluster;
   
   chd = LargeL->head();
   for (; chd; chd = chd->next()){
      cluster = chd->item();
      hdr1 = cluster->list()->head();
      for (; hdr1; hdr1=hdr1->next()){
         add_freq(hdr1->item(), cluster->templ());
         //hdr1->item()->print_seq(cluster->templ());
      }
      hdr2 = cluster->list2()->head();
      for (; hdr2; hdr2=hdr2->next()){
         add_freq(hdr2->item(), cluster->templ2());
         //hdr2->item()->print_seq(cluster->templ2());
      }
   }
}

void process_cluster_list1(ListNodes<Itemset *> *hdr1,
                           Lists<Itemset *> *cluster1,
                           Lists<Itemset *> *cluster2, Lists<Eqclass *> *LargeL,
                           int iter, int eqtype, Eqclass *parent)
{
   ListNodes<Itemset *> *hdr2;
   Eqclass *EQ = new Eqclass(iter-1,eqtype);
   if (EQ == NULL){
      perror("memory exceeded");
      exit(errno);
   }
   fill_seq_template(EQ, parent, 2);
   //int first;
   Itemset *ljoin, *ejoin, *mjoin;
   int lsup, esup, msup;
   //cout << "BEG CLUSERT 1 : " << MEMUSED << endl;

   //first = 1;
   hdr2 = cluster2->head();
   for (; hdr2; hdr2=hdr2->next()){
      //ljoin = (Itemset *)1;
      ljoin = prune_decision(hdr1->item(), hdr2->item(), EQ->templ(), LJOIN,2);
      ejoin = NULL;
      mjoin = NULL;
      lsup = esup = msup = 0;
      //cout << "process 1 0 0" << endl;
      if (ljoin || ejoin || mjoin)
         get_tmpnewf_intersect(ljoin, ejoin, mjoin, lsup, esup, msup,
                           hdr1->item(), hdr2->item(), iter);
      if (lsup >= MINSUPPORT){
         NumLargeItemset[iter-1]++;
         fill_join(ljoin, hdr1->item(), hdr2->item());
         //cout << "LARGE ";
         if (print_seq) ljoin->print_seq(EQ->templ());
         EQ->append(ljoin);
      }
   }
   
   hdr2 = cluster1->head();
   for (; hdr2 != hdr1; hdr2=hdr2->next()){
      //ejoin = (Itemset *)1;
      ejoin = prune_decision(hdr1->item(), hdr2->item(), EQ->templ2(), EJOIN,2);
      ljoin = NULL;
      mjoin = NULL;
      lsup = esup = msup = 0;
      //cout << "process 0 1 0" << endl;
      if (ljoin || ejoin || mjoin)
         get_tmpnewf_intersect(ljoin, ejoin, mjoin, lsup, esup, msup,
                           hdr1->item(), hdr2->item(), iter);
      //cout << "AFT JOIN " << MEMUSED << endl;
      if (esup >= MINSUPPORT){
         NumLargeItemset[iter-1]++;
         fill_join(ejoin, hdr1->item(), hdr2->item());
         //cout << "LARGE ";
         if (print_seq) ejoin->print_seq(EQ->templ2());
         EQ->append2(ejoin);
      }
   }
   
   if (EQ){
      if ((EQ->list()->size() > 0) || (EQ->list2()->size() > 0)){
         if (recursive){
            //if (use_hash) insert_freqarray(EQ);
            process_cluster1(EQ, NULL, iter+1);
            delete EQ;
         }
         else LargeL->append(EQ);
      }
      else{
         //   if (use_hash && EQ->list2()->size() == 1)
         //      add_freq(EQ->list2()->head()->item(), EQ->templ2());
         delete EQ;
         EQ = NULL;
      }
   }
   //cout << "END CLUSTER1 : " << MEMUSED << endl;
}

void process_cluster_list2(ListNodes<Itemset *> *hdr1, int i, Eqclass ** EQ,
                           Lists<Itemset *> *cluster, Lists<Eqclass *> *LargeL,
                           int iter, int eqtype, Eqclass *parent)
{
   int j;
   
   ListNodes<Itemset *> *hdr2;
   Itemset *ljoin, *ejoin, *mjoin;
   int lsup, esup, msup;

   //join with sequences
   hdr2 = hdr1;
   for (j=i; hdr2; j++, hdr2=hdr2->next()){
      ljoin = prune_decision(hdr1->item(), hdr2->item(), EQ[i]->templ(), LJOIN,1);
      if (hdr2 == hdr1){
         ejoin = mjoin = NULL;
      }
      else{
         ejoin = prune_decision(hdr2->item(), hdr1->item(), EQ[j]->templ2(), EJOIN,1);
         mjoin = prune_decision(hdr2->item(), hdr1->item(), EQ[j]->templ(), MJOIN,1);
         //ejoin = mjoin = (Itemset *)1;
      }
      //cout << "process 1 1 1" << endl;
      lsup = esup = msup = 0;
      if (ljoin || ejoin || mjoin)
         get_tmpnewf_intersect(ljoin, ejoin, mjoin, lsup, esup, msup,
                           hdr1->item(), hdr2->item(), iter);
      //cout << "SUPPP " << lsup << " " << esup << " " << msup << endl;
      if (lsup >= MINSUPPORT){
         NumLargeItemset[iter-1]++;
         fill_join(ljoin, hdr1->item(), hdr2->item());
         //cout << "LARGE ";
         if (print_seq) ljoin->print_seq(EQ[i]->templ());
         EQ[i]->append(ljoin);
      }
      if (esup >= MINSUPPORT){
         NumLargeItemset[iter-1]++;
         fill_join(ejoin, hdr2->item(), hdr1->item());
         //cout << "LARGE ";
         if (print_seq) ejoin->print_seq(EQ[j]->templ2());
         EQ[j]->append2(ejoin);
      }
      if (msup >= MINSUPPORT){
         NumLargeItemset[iter-1]++;
         fill_join(mjoin, hdr2->item(), hdr1->item());
         //cout << "LARGE ";
         if (print_seq) mjoin->print_seq(EQ[j]->templ());
         EQ[j]->append(mjoin);
      }      
   }
   if ((EQ[i]->list()->size() > 0) || (EQ[i]->list2()->size() > 0)){
      if (recursive){
         //if (use_hash) insert_freqarray(EQ[i]);
         process_cluster1(EQ[i],NULL, iter+1);
         delete EQ[i];
         EQ[i] = NULL;
      }
      else LargeL->append(EQ[i]);
   }
   else{
      //if (use_hash && EQ[i]->list2()->size() == 1)
      //   add_freq(EQ[i]->list2()->head()->item(), EQ[i]->templ2());
      delete EQ[i];
      EQ[i] = NULL;
   }
   
   //cout << "END cluster 2 : " << MEMUSED << endl;
}



void process_cluster1(Eqclass *cluster, Lists<Eqclass *> *LargeL, int iter)
{
   Eqclass **EQ=NULL;
   ListNodes<Itemset *> *hdr1, *hdr2;
   int i;

   if (cluster->list()->head()){
      EQ = new Eqclass *[cluster->list()->size()];
      if (EQ == NULL){
         perror("memory exceeded");
         exit(errno);
      }
      for (i=0; i < cluster->list()->size(); i++){
         EQ[i] = new Eqclass(iter-1,EQCTYP1);
         if (EQ[i] == NULL){
            perror("memory exceeded");
            exit(errno);
         }
         fill_seq_template(EQ[i], cluster, 1);
      }
   }
   
   hdr1 = cluster->list()->head();
   for (i=0; hdr1; hdr1=hdr1->next(), i++){
      //if (use_hash && iter > 3) add_freq(hdr1->item(), cluster->templ());
      process_cluster_list2(hdr1, i, EQ, cluster->list(), LargeL, iter,
                            EQCTYP1, cluster);
   }
   if (EQ) delete [] EQ;
   
   
   hdr2 = cluster->list2()->head();
   for (; hdr2; hdr2=hdr2->next()){
      //if (use_hash && iter > 3) add_freq(hdr2->item(), cluster->templ2());
      process_cluster_list1(hdr2, cluster->list2(), cluster->list(),
                            LargeL, iter, EQCTYP1, cluster);
   }
   
   //if (recursive) delete cluster;
   if (maxiter < iter) maxiter = iter;
   
}


void find_large(Eqclass *cluster, int it)
{
   Lists<Eqclass *> *LargeL, *Candidate;
   ListNodes<Eqclass *> *chd;
   int iter;
   int LargelistSum=0;
   int more;
   
   more = 1;
   Candidate = new Lists<Eqclass *>;
   Candidate->append(cluster);
   //cout << "MEMFIND " << it << " " << MEMUSED << endl;
   for (iter=3; more; iter++){
      LargeL = new Lists<Eqclass *>;
      chd = Candidate->head();
      for (; chd; chd=chd->next()){
         //cout << "EQCLASS ";
         //chd->item()->print_template();
         //chd->item()->print_list(chd->item()->list());
         //cout << "***\n";
         //chd->item()->print_list(chd->item()->list2());
         //cout << "------------------" << endl;
         process_cluster1(chd->item(), LargeL, iter);
         //cout << "BEF MEMFIND " << it << " " << MEMUSED << endl;
         //reclaim memory for this class immediately
         delete chd->item();
         //cout << "AFT MEMFIND " << it << " " << MEMUSED << endl;
         chd->set_item(NULL);
      }
      Candidate->clear();
      delete Candidate;
      //if (maxiter < iter) maxiter = iter;
      
      if (use_hash) insert_freqarray(LargeL);
      chd = LargeL->head();
      LargelistSum = 0;
      for (;chd; chd=chd->next()){
         LargelistSum += chd->item()->list()->size();
         if (chd->item()->list2())
            LargelistSum += chd->item()->list2()->size();
      }
      //print_freqary();
      more = (LargelistSum > 0);
      
      Candidate = LargeL;
      if (memtrace) mout << it << " " << MEMUSED << endl;

      if (!more) {
         LargeL->clear();
         delete LargeL;
      }
      //cout << "AFT DEL " << it << " " << MEMUSED << " " << iter << endl;
   }
   //cout << "MEMLAST " << it << " " << MEMUSED << endl;
}

Eqclass * extract_relevant_items(Eqclass *l2it, Array *cliq)
{
   int i;

   Eqclass *RL2 = new Eqclass(1, EQCTYP1);
   if (RL2 == NULL)
   {
      perror("memory exceeded : ext_class ");
      exit (errno);
   }
   //init seq pattern templates
   RL2->set_templ(1);
   RL2->set_templ2(0);
    
   ListNodes<Itemset *> *ln= l2it->list()->head();
   for (i=0; ln && i < cliq->size()-1; ){
      //cout << "LN " << (*ln->item())[1] << " " << (*cliq)[i+1] << endl;
      if ((*ln->item())[1] == (*cliq)[i+1]){
         RL2->append(ln->item());
         i++;
      }
      ln = ln->next();
   }
   ln= l2it->list2()->head();
   for (i=0; ln && i < cliq->size()-1; ){
      if ((*ln->item())[1] == (*cliq)[i+1]){
         RL2->append2(ln->item());
         i++;
      }
      ln = ln->next();
   }   
   return RL2;
}

void process_class(int it)
{
   //from 2-itemsets from ext disk
   Eqclass *large2it = get_ext_eqclass(it);
   if (large2it == NULL) return;

   if (memtrace) mout << it << " " << MEMUSED << endl;
   Eqclass *l2cliq;
   
   if (use_clique){
      ListNodes<Array *> *clhd = eqgraph[it]->clique()->head();   
      //process each clique
      for (;clhd; clhd=clhd->next()){
         //cout << it << " processing clique " << *clhd->item() << endl;
         //construct large k-itemsets, k > 2
         l2cliq = extract_relevant_items(large2it, clhd->item());
         find_large(l2cliq, it);
      }  
   }
   else{
      if (recursive){
         process_cluster1(large2it, NULL, 3);
         delete large2it;
      }
      else find_large(large2it, it);
   }
   
}

void newSeq()
{
   int i,j;

   if (use_hash)
      FreqArray = (FreqIt **) malloc (FreqArraySz*sizeof(FreqIt*));
   //form large itemsets for each eqclass
   if (use_ascending != -2){
      if (use_ascending == -1){
         for (i=0; i < DBASE_MAXITEM; i++)
            if (eqgraph[i]){
               if (memtrace) mout << i << " " << MEMUSED << endl;
               process_class(i);
               if (memtrace) mout << i << " " << MEMUSED << endl;
            }
      }
      else if (eqgraph[use_ascending])
         process_class(use_ascending);
   }
   else{
      for (i=DBASE_MAXITEM-1; i >= 0; i--){
         if (eqgraph[i]){
            if (memtrace) mout << i << " " << MEMUSED << endl;
            //cout << "PROCESSS ITEM " << i << endl << flush;
            if (use_hash) FreqArrayPos = 0;
            process_class(i);
            if (use_hash){
               if (FreqArrayPos > 0){
                  //cout << "FREQUENT ARRAY3" << endl;
                  FreqIt **fit = new FreqIt *[FreqArrayPos];
                  for (j=0; j < FreqArrayPos; j++){
                     fit[j] = FreqArray[j];
                     //cout << *fit[j];
                  }
                  eqgraph[i]->set_freqarray(fit, FreqArrayPos);
               }
            }
            //cout << " -------- " << endl;
            if (memtrace) mout << i << " " << MEMUSED << endl;
         }
      }
   }
}


void read_files()
{
   int i;
   
   NumLargeItemset = new int [(int)(DBASE_AVG_TRANS_SZ*30)];
   bzero((char *)NumLargeItemset, sizeof(int)*((int)(DBASE_AVG_TRANS_SZ*30)));
   eqgraph = new EqGrNode *[DBASE_MAXITEM];   
   bzero((char *)eqgraph, DBASE_MAXITEM*sizeof(EqGrNode *));   

   double t1,t2;
   if (ext_l2_pass){
      seconds(t1);
      NumLargeItemset[0] = make_l1_pass(1);
      seconds(t2);
      EXTL1TIME = t2-t1;
      NumLargeItemset[1] = make_l2_pass();
      //cout << "L2 " <<  NumLargeItemset[1] <<endl;
      seconds(t1);
      EXTL2TIME = t1-t2;
   }
   else{
      seconds(t1);
      NumLargeItemset[0] = make_l1_pass(0);
      seconds(t2);
      EXTL1TIME = t2-t1;
      NumLargeItemset[1] = get_file_l2(it2f, seqf);
      seconds(t1);
      EXTL2TIME = t1-t2;
   }
   //cout << NumLargeItemset[0] << "LARGE 1 ITEMS\n";
   maxitemsup = 0;
   int sup;
   for (i=0; i < DBASE_MAXITEM; i++) {
      sup = partition_get_idxsup(i);
      if (maxitemsup < sup) maxitemsup  = sup;
   }
   //cout << "MAXITEMSUP " << maxitemsup << endl;
   interval = new Array(maxitemsup);
   interval2 = new Array(maxitemsup);
   interval3 = new Array(maxitemsup);
   //cout << "MAXEQSZIE " << maxeqsize << " " << t2-ts << endl;
}

int main(int argc, char **argv)
{
   int i;

   double ts, te;
   double t1,t2;

   seconds(ts);
   //cout << "BEGIN MEM " << MEMUSED << endl;
   parse_args(argc, argv);

   partition_alloc(dataf, idxf);
   read_files();
   
   //cout << "AFTER READFILE " << MEMUSED << endl;
   seconds(t1);
   newSeq();
   seconds(t2);
   double FKtime = t2-t1;
   //print_freqary();
   //cout << "AFTER SEQ " << MEMUSED << endl;
   seconds(te);
   if ((out = fopen("summary.out", "a+")) == NULL){
      perror("can't open summary file");
      exit(errno);
   }
   fprintf (out, "SPADE ");
   if (use_hash) fprintf (out, "USEHASH ");
   fprintf(out, "%s %f %d %d %f %d (", dataf, MINSUP_PER, MINSUPPORT,
           num_intersect, L2ISECTTIME, L2pruning);
   for (i=0; i < maxiter; i++){
      fprintf(out, "%d ", NumLargeItemset[i]);
      //cout << "ITER " <<  i+1 << " " << NumLargeItemset[i] << endl;
   }
   //cout << "Total elapsed time " << te-ts
   //     << ", NumIntersect " << num_intersect << " L2time "
   //     << L2TIME <<  " " << EXTL2TIME << endl;
   fprintf(out, ") %f %f %f %f\n", EXTL1TIME, EXTL2TIME, FKtime,
           te-ts);
   fclose(out);

   partition_dealloc();

   delete interval;
   delete interval2;
   delete interval3;
   for (i=0; i < DBASE_MAXITEM; i++)
      if (eqgraph[i]) delete eqgraph[i];
   delete [] eqgraph;
   
   if (memtrace){
      mout <<  MEMUSED << endl;
      mout.close();
   }
   //cout << "LAST MEM " << MEMUSED << endl;
   exit(0);
}


