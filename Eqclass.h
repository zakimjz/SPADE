#ifndef _EQCLASS_H
#define _EQCLASS_H

#include <iostream>
#include <errno.h>
#include "Lists.h"
#include "Itemset.h"

#define EQCTYP1 1
#define EQCTYP2 2
#define EQCTYP3 3

class EqGrNode;

class Eqclass {
private:
   Lists<Itemset *> *theList;
   int Iset_size;
   unsigned int seqTemplate;
   Lists<Itemset *> *theList2;
   unsigned int seqTemplate2;
   int Eqtype;
 public:
   Eqclass(int iset_sz, int eqt);
   
   ~Eqclass();

   int templ_sz()
   {
      return Iset_size;
   }
   int eqtype()
   {
      return Eqtype;
   }
   Lists<Itemset *> * list()
   {
      return theList;
   }
   Lists<Itemset *> * list2()
   {
      return theList2;
   }
   unsigned int templ()
   {
      return seqTemplate;
   }
   unsigned int templ2()
   {
      return seqTemplate2;
   }
   void set_templ(unsigned int val)
   {
      seqTemplate = val;
   }
   void set_templ2(unsigned int val)
   {
      seqTemplate2 = val;
   }
   void set_list(Lists<Itemset *> * ll)
   {
      theList = ll;
   }
   void append(Itemset *it)
   {
      theList->append(it);
   }
   void append2(Itemset *it)
   {
      theList2->append(it);
   }
   void prepend(Itemset *it)
   {
      theList->prepend(it);
   }
   void prepend2(Itemset *it)
   {
      theList2->prepend(it);
   }
   void print_template1();
   void print_template2();
   void print_template();
   void print_list(Lists<Itemset *> *ll);
   Itemset * uniqsorted(Itemset *it, CMP_FUNC func);
   int subseq(Itemset *it);
};

class FreqIt{
public:
   int *seq;
   int seqsz;
   unsigned int templ;

   FreqIt(int sz, unsigned int tpl)
   {
      templ = tpl;
      seqsz = sz;
      seq = new int[sz];
   }

   FreqIt(int *ary, int sz, unsigned int tpl)
   {
      templ = tpl;
      seqsz = sz;
      seq = new int[sz];
      for (int i=0; i < sz; i++) seq[i] = ary[i];
   }

   ~FreqIt()
   {
      if (seq) delete [] seq;
   }
   
   int compare(Itemset *iset, unsigned int itpl);
   int compare(FreqIt *fit, int recursive);
   friend ostream& operator << (ostream& outputStream, FreqIt& freq)
   {
      outputStream << "FREQ : ";
      for (int i=0; i < freq.seqsz; i++)
         outputStream << " " << freq.seq[i];
      outputStream << " --- " << freq.templ << endl;
      return outputStream;
   }
};

class EqGrNode {
private:
   Array *theElements;
   //int numElements;
   //int totElements;
   Array *stheElements;
   //int snumElements;
   //int stotElements;

   FreqIt **freqArray; //frequent seq from this class
   int freqArraySz;
   
   Lists<int *> *theCoverL;
   Lists<Array *> *theCliqueL;
   Lists<Itemset *> *theLargeL;
   int theFlg; //indicates if class is in memory

public:
   static int bsearch(int min, int max, FreqIt **freqArray,
                      FreqIt &fit, int recursive);
   static int bsearch(int min, int max, int *itary, int it);

   EqGrNode(int sz);
   ~EqGrNode();

   FreqIt **freqarray()
   {
      return freqArray;
   }
   int freqarraysz()
   {
      return freqArraySz;
   }
   void set_freqarray(FreqIt **fit, int sz)
   {
      freqArray = fit;
      freqArraySz = sz;
   }
   int find_freqarray(FreqIt &fit, int recursive);
   
   int getflg()
   {
      return theFlg;
   }
   void setflg(int val)
   {
      theFlg=val;
   }

   Lists<Array *> *clique()
   {
      return theCliqueL;
   }
   
   Lists<Itemset *> *largelist()
   {
      return theLargeL;
   }

   Lists<int *> * cover()
   {
      return theCoverL;
   }
   Array * elements()
   {
      return theElements;
   }
   int num_elements()
   {
      if (theElements) 
         return theElements->size();
      else return 0;
   }
   void add_element(int el)
   {
      //theElements[numElements] = el;
      //numElements++;
      theElements->add(el);
   }
   
   void add_element(int el, int pos)
   {
      theElements->setitem(pos,el);
   }
   int get_element(int pos)
   {
      return (*theElements)[pos];
   }
   //void remove_el(int pos)
   //{
   //   for (int i=pos; i < numElements-1; i++)
   //      theElements[i] = theElements[i+1];
   //   numElements--;
   //}

   void seqsetelements(Array *ary)
   {
      stheElements = ary;
      //stotElements = sz;
      //MEMUSED += sz*sizeof(int);
   }
   
   Array * seqelements()
   {
      return stheElements;
   }
   int seqnum_elements()
   {
      if (stheElements)
         return stheElements->size();
      else return 0;
   }
   void seqadd_element(int el)
   {
      stheElements->add( el);
      //snumElements++;
   }
   
   void seqadd_element(int el, int pos)
   {
      stheElements->setitem(pos,el);
   }
   int seqget_element(int pos)
   {
      return (*stheElements)[pos];
   }

   
   void add_cover(int *eq)
   {
      theCoverL->append(eq);
   }
   int find(int it)
   {
      if (theElements){
         //for (int i=0; i < theElements->size(); i++)
         //   if ((*theElements)[i] == it) return 1;
         return bsearch(0, theElements->size()-1, theElements->array(), it);
      }
      return 0;
   }
   int seqfind(int it)
   {
      if (stheElements){
         //for (int i=0; i < stheElements->size(); i++)
         //   if ((*stheElements)[i] == it) return 1;
         return bsearch(0, stheElements->size()-1, stheElements->array(), it);
      }
      return 0;
   }
   friend ostream& operator << (ostream& outputStream, EqGrNode& EQ);
};

extern void eq_insert(Lists<Eqclass *> &EQC, Itemset *it);
extern void eq_print(Eqclass *LargeEqclass);
#endif

