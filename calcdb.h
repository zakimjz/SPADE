#ifndef __DATABASE_H
#define __DATABASE_H

#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>


#define ITSZ sizeof(int)
#define DCBBUFSZ 2048

class Dbase_Ctrl_Blk{
public:
   Dbase_Ctrl_Blk(char *infile, int buf_sz=DCBBUFSZ);
   ~Dbase_Ctrl_Blk();

   void get_next_trans_ext();
   inline void get_first_blk();
   inline void get_next_trans(int *&lbuf, int &numitem, int &tid, int &custid);

   int eof()
   {
      return (readall == 1);
   }
   int fd;     
   int buf_size;
   int * buf;
   int cur_blk_size; 
   int cur_buf_pos;
   int endpos;
   char readall;
};

inline void Dbase_Ctrl_Blk::get_first_blk()
{
   readall=0;
   lseek(fd, 0, SEEK_SET);
   cur_blk_size = (read(fd,(void *)buf, (buf_size*ITSZ)))/ITSZ;
   if (cur_blk_size < 0){
      perror("get_first_blk");
      exit(errno);
   }
   cur_buf_pos = 0;
}

inline  void Dbase_Ctrl_Blk::get_next_trans (int *&lbuf,
                                             int &nitems, int &tid, int &cid)
{
   if (cur_buf_pos+3 >= cur_blk_size ||
       cur_buf_pos+buf[cur_buf_pos+2]+3 > cur_blk_size){
      if (lseek(fd, 0, SEEK_CUR) == endpos) readall = 1;
      if (!readall){
         // Need to get more items from file
         get_next_trans_ext();
      }      
   }
   
   if (!readall){
      cid = buf[cur_buf_pos];
      tid = buf[cur_buf_pos+1];
      nitems = buf[cur_buf_pos+2];
      //if ((cur_buf_pos + nitems + 3) > cur_blk_size)
      // {
      //   if (cur_buf_pos + nitems > cur_blk_size &&
      //       lseek(fd, 0, SEEK_CUR) == endpos) readall = 1;
      //   if (!readall){
      //      // Need to get more items from file
      //      get_next_trans_ext(nitems, tid, cid);
      //   }
      //}
      lbuf = buf + cur_buf_pos + 3;
      cur_buf_pos += nitems + 3;
   }
}
#endif //__DATABASE_H





