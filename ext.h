#ifndef __EXT_H_
#define __EXT_H_
#include <sys/time.h>

extern struct timeval tp;
#define seconds(tm) gettimeofday(&tp,(struct timezone *)0);\
tm=tp.tv_sec+tp.tv_usec/1000000.0
#define min(a, b) ((a) < (b) ? (a) : (b))

extern int EXTBLKSZ;
extern int num_partitions;
// extern void gen_ext_l2(char *it2f, char *seq2f,
//                        int data_fd=-1, int *item_idx=NULL);
extern int gen_ext_l2(int data_fd=-1, int *item_idx=NULL);
extern int make_l1_pass(int data_fd, int *item_idx, int bflg, char use_newformat);
extern int gen_l2_multpass(int data_fd=-1, int *item_idx=NULL);

class Extary
{
public:
   int bufsz;
   int bufpos;
   int *buf;
   int filepos;
   int endpos;
   char readall;
   Extary(int sz=0)
   {
      bufsz = sz;
      readall = 0;
      bufpos = 0;
      endpos = 0;
      if (bufsz > 0) buf = (int *) malloc (sizeof(int)*bufsz);
      else buf = NULL;
      filepos = 0;
   }
   ~Extary()
   {
      if (buf) free(buf);
   }
   void set_buf(int *bufadd)
   {
      buf = bufadd;
   }
   void reset()
   {
      bufpos = 0;
      filepos = 0;
      endpos = 0;
   }
   int check_empty()
   {
      if (filepos >= endpos){
         filepos = -1;
         return 1;
      }
      return 0;
   }
   int isempty()
   {
      return (filepos == -1);
   }

   void set_filepos(int fpos)
   {
      filepos = fpos;
   }
   
   int custid()
   {
      return buf[bufpos];
   }

   int get_min_max(int &mmin, int &mmax)
   {
      mmin = -1;
      mmax = -1;
      if (bufpos < bufsz){
         int pos = bufpos;
         int fpos = filepos;
         mmin = buf[pos];
         mmax = buf[pos];
         fpos += (buf[pos+1]+2);
         pos += (buf[pos+1]+2);
         while(pos < bufsz && fpos < endpos){
            //ensure that entire cust fits in buf
            if (pos+buf[pos+1]+2 > bufsz) break;
            
            mmax = buf[pos];
            fpos += (buf[pos+1]+2);
            pos += (buf[pos+1]+2);
         }
         if (fpos >= endpos) return 1;
      }
      return 0;
   }
   //to be used only for partitioned L2
//    int next_custid(int idx=1)
//    {
//       int i, pos;
//       pos = bufpos;
//       for (i=0; i < idx; i++){
//          pos += (pntid(pos)+2);
//          if (pos >= endpos)
//             return -1;
//       }
//       return buf[pos];
//    }
   int next_custid(int pos)
   {
      if (pos >= endpos) return -1;
      return buf[pos];
   }
   int pntid(int pos)
   {
      return buf[pos+1];
   }
   int ntid()
   {
      return buf[bufpos+1];
   }
   int *buff()
   {
      return (&buf[bufpos]);
   }
   int fbuf()
   {
      return buf[bufpos+2];
   }
   int lbuf()
   {
      return buf[bufpos+ntid()+1];
   }
};

#endif //__EXT_H_



