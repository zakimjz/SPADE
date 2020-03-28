#include <errno.h>
#include "Itemset.h"

Itemset::Itemset(int it_sz, int ival_sz){
   //cout << "ITALLOC " << MEMUSED;   
   theItemset = new Array(it_sz);
   if (theItemset == NULL){
      perror("memory:: Itemset");
      exit(errno);
   }
   theIval = new Array(ival_sz);
   if (theIval == NULL){
      perror("memory:: Ival");
      exit(errno);
   }
   
   //for (int i=0; i < ival_sz; i++)
   //   theIval[i] = NULL;
   theSupport = 0;
   MEMUSED += sizeof(Itemset);
   //cout << " -- " << MEMUSED << endl;   
}

Itemset::~Itemset(){
   //cout << "ITDEL " << MEMUSED;   
   if (theItemset) delete theItemset;
   if (theIval) delete theIval;
   theItemset = NULL;
   theSupport = 0;
   MEMUSED -= sizeof(Itemset);
   //cout << " -- " << MEMUSED << endl;   
}

ostream& operator << (ostream& outputStream, Itemset& itemset){
   outputStream << "ITEM: ";
   outputStream << *itemset.theItemset << " ";
   outputStream << " Support:";
   outputStream << itemset.theSupport;
   outputStream << "\n";
   return outputStream;
}

int Itemset::compare(Itemset& ar2)
{
   int len;
   if (size() <= ar2.size()) len = size();
   else len = ar2.size();
   for(int i=0; i < len; i++){
      if ((*theItemset)[i] > (*ar2.theItemset)[i]) return 1;
      else if ((*theItemset)[i] < (*ar2.theItemset)[i]) return -1;
   }
   if (size() < ar2.size()) return -1;
   else if (size() > ar2.size()) return 1;
   else return 0;
}

//len must be less than length of both Itemsets
int Itemset::compare(Itemset& ar2, int len)
{
   for(int i=0; i < len; i++){
      if ((*theItemset)[i] > (*ar2.theItemset)[i]) return 1;
      else if ((*theItemset)[i] < (*ar2.theItemset)[i]) return -1;
   }
   return 0;
}
int Itemset::compare(Array& ar2, int len)
{
   for(int i=0; i < len; i++){
      if ((*theItemset)[i] > ar2[i]) return 1;
      else if ((*theItemset)[i] < ar2[i]) return -1;
   }
   return 0;
}

int Itemset::compare(Itemset& ar2, int len, unsigned int bvec)
{
   int pos = 0;
   int it;
   for(int i=0; i < len; i++){
      while (!GETBIT(bvec, pos)){
         pos ++;
      }
      it = (*theItemset)[pos++];
      if (it > (*ar2.theItemset)[i]) return 1;
      else if (it < (*ar2.theItemset)[i]) return -1;
   }
   return 0;
}

int Itemset::subsequence(Itemset * ar)
{
   int i,j;
   if (size() > ar->size()) return 0;
   int start = 0;
   for(i=0; i < size(); i++){
      for(j=start; j < ar->size(); j++){
         if ((*theItemset)[i] == (*ar->theItemset)[j]){
            start = j+1;
            break;
         }
      }
      if (j >= ar->size()) return 0;
   }
   return 1;
}

void Itemset::print_seq(int itempl)
{
   int i;
   int sz = size();
   //int templ = itempl;
   //int mask = 1 << (size()-1);
   cout << (*theItemset)[0] << " ";
   for (i=1; i < sz-1; i++){
      //if (templ && mask) cout << "->";
      //templ = templ << 1;
      if (GETBIT(itempl,sz-1-i))
         cout << "-> ";
      cout << (*theItemset)[i] << " ";
   }
   if (GETBIT(itempl,sz-1-i))
      cout << "-> ";
   cout << (*theItemset)[sz-1] << " ";
   cout << "- " << theSupport << endl;      
}



