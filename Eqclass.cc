#include <sys/types.h>
#include <unistd.h>
#include <errno.h>

#include "Eqclass.h"


Eqclass::Eqclass(int iset_sz, int eqt){
   Iset_size = iset_sz;
   Eqtype = eqt;
   theList = new Lists<Itemset *>;
   if (theList == NULL){
      perror("memory :: Eqclass");
      exit(errno);
   }
   seqTemplate = seqTemplate2 = 0;
   //seqTemplate = new int[Iset_size];
   theList2 = NULL;
   //seqTemplate2 = NULL;
   if (Eqtype == EQCTYP1){
      theList2 = new Lists<Itemset *>;
      if (theList2 == NULL){
         perror("memory :: Eqclass");
         exit(errno);
      }
      //seqTemplate2 = new int[Iset_size];
   }
   MEMUSED += sizeof(Eqclass);
}

Eqclass::~Eqclass(){
   if (theList){
      theList->clear();
      delete theList;
   }
   theList = NULL;
   //if (seqTemplate) delete [] seqTemplate;
   if (theList2){
      theList2->clear();
      delete theList2;
   }
   theList2 = NULL;
   //if (seqTemplate2) delete [] seqTemplate2;
   MEMUSED -= sizeof(Eqclass);
}
void Eqclass::print_template1()
{
   cout << "TEMPLATE1";
   cout << " " << seqTemplate;
   //for (i=0; i < Iset_size; i++) cout << " " << seqTemplate[i];
   cout << endl;
}
void Eqclass::print_template2()
{
   if (seqTemplate2){
      cout << "TEMPLATE2";
      cout << " " << seqTemplate2;
      //for (i=0; i < Iset_size; i++) cout << " " << seqTemplate2[i];
      cout << endl;
   }
}
void Eqclass::print_template()
{
   cout << "TEMPLATE1";
   cout << " " << seqTemplate;
   //for (i=0; i < Iset_size; i++) cout << " " << seqTemplate[i];
   //cout << endl;
   //if (seqTemplate2){
      cout << " TEMPLATE2";
      cout << " " << seqTemplate2;
      //for (i=0; i < Iset_size; i++) cout << " " << seqTemplate2[i];
      cout << endl;
      //}
}
void Eqclass::print_list(Lists<Itemset *> *ll)
{
   ListNodes<Itemset *> *hd = ll->head();
   for (;hd; hd= hd->next())
      cout << *(hd->item());
   cout << flush;
}

Itemset * Eqclass::uniqsorted(Itemset *it, CMP_FUNC func)
{
   Itemset *rval;
   //ListNodes<Itemset *> * prev = NULL;
   //if (!theList->find_ascend(prev, it, func))
   //   theList->insert(prev, it);
   if (!(rval = theList->find(it, Itemset::Itemcompare))){
      theList->sortedAscend(it, func);
      //  return 0;
      //}
      //else return 1;
   }
   return rval;
}

int Eqclass::subseq(Itemset *it)
{
   ListNodes<Itemset *> *hd = theList->head();
   for (;hd; hd=hd->next()){
      if (it->subsequence(hd->item())){
         return 1;
      }
   }
   return 0;
}


EqGrNode::EqGrNode(int sz)
{
   //totElements = sz;
   //numElements = 0;
   if (sz > 0){
      theElements = new Array(sz);
   }
   else theElements = NULL;
   stheElements = NULL;
   //stheElements = new int[ssz];
   //MEMUSED += sz*sizeof(int);
   freqArray = NULL;
   freqArraySz = 0;
   theCoverL = new Lists<int *>;
   theLargeL = new Lists<Itemset *>;
   theCliqueL = new Lists<Array *>;
   theFlg = 0;
   MEMUSED += sizeof(EqGrNode);
}

EqGrNode::~EqGrNode()
{
   if (theElements) delete theElements;
   if (stheElements) delete stheElements;
   if (freqArray){
      for (int i=0; i < freqArraySz; i++) delete freqArray[i];
      delete [] freqArray;
   }
   delete theCoverL;
   delete theCliqueL;
   delete theLargeL;
   theElements = NULL;
   theCoverL = NULL;
   theFlg = 0;
   MEMUSED -= sizeof(EqGrNode);
}

//assume that elements are sorted in descending order
int EqGrNode::bsearch(int min, int max, FreqIt **freqArray,
                             FreqIt &fit, int recursive)
{
   int mid = (max+min)/2;
   if (max < min) return 0;

   int res = freqArray[mid]->compare(&fit, recursive);
   //if (fit.seq[0] == 101 || fit.seq[0] == 201)
   //   cout << "RES " << res << " " << mid << " *** " << *freqArray[mid];
   if (res == 0) return 1;
   else if (res < 0) return bsearch(min, mid-1, freqArray, fit, recursive);
   else return bsearch(mid+1, max, freqArray, fit, recursive);
}

int EqGrNode::bsearch(int min, int max, int *itary, int it)
{
   int mid = (max+min)/2;
   if (max < min) return 0;

   if (it == itary[mid]) return 1;
   else if (it < itary[mid]) return bsearch(min, mid-1, itary, it);
   else return bsearch(mid+1, max, itary, it);
}


int EqGrNode::find_freqarray(FreqIt &fit, int recursive)
{
   if (freqArraySz > 0)
      return bsearch(0, freqArraySz-1, freqArray, fit, recursive);
   else return 0;
}

//both must be of same length
int FreqIt::compare(Itemset *iset, unsigned int itpl)
{
   for (int i=0; i < seqsz; i++){
      if (seq[i] > (*iset)[i]) return 1;
      else if (seq[i] < (*iset)[i]) return -1;
   }
   if (templ > itpl) return 1;
   else if (templ < itpl) return -1;
   return 0;
}


int FreqIt::compare(FreqIt *fit, int recursive)
{
   int i;
   
   //first compare seqsz, one with larger seqsz is smaller
   if (seqsz > fit->seqsz) return -1;
   else if (seqsz < fit->seqsz) return 1;

   //compare items & template bits
   if (seq[0] < fit->seq[0]) return -1;
   else if (seq[0] > fit->seq[0]) return 1;
   
   int bpos = seqsz-1;
   int b1, b2;
   for (i=1; i < seqsz; i++){
      b1 = GETBIT(templ, bpos-i);
      b2 = GETBIT(fit->templ, bpos-i);
      if (b1 < b2) return -1;
      else if (b1 > b2) return 1;
      
      if (seq[i] < fit->seq[i]) return -1;
      else if (seq[i] > fit->seq[i]) return 1;
   }
   return 0;
}



