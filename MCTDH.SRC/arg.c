

//   BE CAREFUL USING NEXTRAN!!   WITH MPI ALL PROCS MUST BY SYNCED.


#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#ifdef PGFFLAG

    int myiargc_() {         
      extern int    __argc_save;
      return __argc_save - 1;
    } 

    void mygetarg_(int* i, char** buffer) {          
      extern char **__argv_save;

      printf("Go mygetarg! \n");

      strcpy(*buffer, __argv_save[*i+1]);

    //	printf("Arg %i is %50c \n",*i+1,*buffer);
    } 
#endif


unsigned short int bigendian(unsigned short int *littleendian)
{
  int i, ipow,j,ii;

  const int numpow=16;

  i=0;
  j=*littleendian;
  for (ipow=0; ipow<numpow; ipow++) {
    ii = j%2;
    j=(j-ii)/2;
    i=i + ii * pow(2,(ipow+numpow/2)%numpow);
  }

  return i;

}



void writepovray_(int *initwo, double* threedensity, char *filename, long int string_len)
{

  FILE*fileptr;
  int ii, i, isize;
  char mystring[100];
  double max;
  unsigned short int itwo;
  unsigned short int xitwo;
  //  unsigned int iithree[(*initwo)*(*initwo)*(*initwo)];
  unsigned short int iithree[(*initwo)*(*initwo)*(*initwo)];

  unsigned short int imax= 2* 2* 2* 2* 2* 2* 2* 2* 2* 2* 2* 2* 2* 2* 2 *2  - 1 ;


  xitwo=*initwo;
  itwo=bigendian(&xitwo);

  isize=(*initwo)*(*initwo)*(*initwo);

  //   printf("Go %i %i %i \n",isize, itwo, *initwo);


  max=0.0;
  for (i=0; i<isize; i++) {
    if (threedensity[i] > max) max=threedensity[i];
  }


      
  for (i=0; i<isize; i++) {
    if (threedensity[i] > 0.0) {
      iithree[i]=threedensity[i]*imax;
      iithree[i]=bigendian(&(iithree[i]));
    }else{
      iithree[i]=0;
    }
  }


  for (i=0; i<string_len; i++) mystring[i]=filename[i];
  mystring[i]='\0';


  fileptr=fopen(mystring,"wb");


  fwrite(&itwo,sizeof(itwo),1,fileptr);



  fwrite(&itwo,sizeof(itwo),1,fileptr);
  fwrite(&itwo,sizeof(itwo),1,fileptr);
  fwrite(iithree,sizeof(iithree),1,fileptr);
  fclose(fileptr);
  

}

void rand_init_( unsigned int initval )
{

  // TEMP  srand(initval);

  srand(123456);

}


double nextran_()
{

  int iiran;
  
  iiran=rand();

  return 2.0*(double)iiran/(double)RAND_MAX - 1.0;

}


double nextran() {  return nextran_(); }


