/*
 * Nifti and Analyze IO stuff
 *
 * with lots of amendments from Steve Smith, FMRIB
 * and some suggestions from Chris Rorden, AKA MRICRO
 *
 * substantially re-written by Mark Jenkinson, FMRIB
 *
 */


#include "libss.h"
#include "libavw.h"


int fsl_imageexists(const char *filename)
{
  /* tests whether a hdr/img/hdr.gz/img.gz pair exists (in any combination) */
  return FslFileExists(FslMakeBaseName(filename));
}


void avw_read_basic(const char *filename, image_struct *image, int read_data)
{

  FSLIO *fslio;
  short x,y,z,t;
  float xv,yv,zv,tr;
  int size, i;
  void *tmpbuf;
  float slope, intercept;
  int doscaling=0;

  fslio = FslOpen(FslMakeBaseName(filename),"rb");

  image->miscinfo = FslInit();
  FslCloneHeader(image->miscinfo,fslio);

  FslGetDim(fslio,&x,&y,&z,&t);
  image->x=MAX(1,x);
  image->y=MAX(1,y);
  image->z=MAX(1,z);
  image->t=MAX(1,t);

  size=image->x *image->y *image->z *image->t;

  FslGetVoxDim(fslio,&xv,&yv,&zv,&tr);
  image->xv0=xv;
  image->yv0=yv;
  image->zv0=zv;

  image->tr=fabs(tr);

  image->xv=fabs(image->xv0);
  image->yv=fabs(image->yv0);
  image->zv=fabs(image->zv0);

  if (image->xv<0.00001) { image->xv=1; image->xv0=1; }
  if (image->yv<0.00001) { image->yv=1; image->yv0=1; }
  if (image->zv<0.00001) { image->zv=1; image->zv0=1; }

  FslGetIntent(fslio, &(image->intent_code), &(image->intent_p1), &(image->intent_p2), 
	       &(image->intent_p3));


  /* this is the slice ordering field => TODO list */
  /* image->info=(*hdr).dime.funused1; */
  image->info = 0;
  
  FslGetCalMinMax(fslio,&xv,&yv);
  image->min = (FDT) xv;
  image->max = (FDT) yv;

  FslGetAuxFile(fslio,image->lut);
  image->lut[23]=0;
 
  image->bpv=FslGetDataType(fslio,&x)/8;
  image->dt = (int) x;


  /** quit out if we don't need to read the data **/
  if (!read_data) {
    image->i = NULL;
    FslClose(fslio);
    return;
  }

  if ((tmpbuf=malloc(image->bpv * size))==NULL)
    {
      fprintf(stderr,"Aaargh not enough memory!\n");
      exit(1);
    }
  
  if (sizeof(FDT)>image->bpv)
    {
      if ((image->i=(FDT*)malloc(sizeof(FDT) * size))==NULL) 
	/* will need more memory than the file size */
	{
	  fprintf(stderr,"Aaargh not enough memory!\n");
	  exit(1);
	}
    }
  else
    image->i=(FDT*)tmpbuf;

  FslReadVolumes(fslio,(void *) tmpbuf,image->t);
  FslClose(fslio);
   
  doscaling = FslGetIntensityScaling(fslio,&slope,&intercept);

  switch (image->dt) {
  case DT_UNSIGNED_CHAR:
    if (!doscaling) { 
      if ( (strcmp(FDTS,"unsigned char")!=0) && (strcmp(FDTS,"unsignedchar")!=0) )
	{ 
	  for(i=0;i<size;i++) image->i[i] = (FDT) ((unsigned char*)tmpbuf)[i];
	}
    } else { 
      for(i=0;i<size;i++) 
	image->i[i] = (FDT) ((slope * ((unsigned char*)tmpbuf)[i]) + intercept);  
    }
    break;
  case DT_SIGNED_SHORT:
    if (!doscaling) { 
      if ( (strcmp(FDTS,"signed short")!=0) && (strcmp(FDTS,"signedshort")!=0) )
	{ 
	  for(i=0;i<size;i++) image->i[i] = (FDT) ((short*)tmpbuf)[i]; 
	}
    } else { 
      for(i=0;i<size;i++) 
	image->i[i] = (FDT) ((slope * ((short*)tmpbuf)[i]) + intercept);  
    }
    break;
  case DT_SIGNED_INT:
    if (!doscaling) { 
      if ( (strcmp(FDTS,"signed int")!=0) && (strcmp(FDTS,"signedint")!=0) )
	{ 
	  for(i=0;i<size;i++) image->i[i] = (FDT) ((int*)tmpbuf)[i]; 
	}
    } else { 
      for(i=0;i<size;i++) 
	image->i[i] = (FDT) ((slope * ((int*)tmpbuf)[i]) + intercept);  
    }
    break;
  case DT_FLOAT:
    if (!doscaling) { 
      if (strcmp(FDTS,"float")!=0)
	{ 
	  for(i=0;i<size;i++) image->i[i] = (FDT) ((float*)tmpbuf)[i]; 
	}
    } else { 
      for(i=0;i<size;i++) 
	image->i[i] = (FDT) ((slope * ((float*)tmpbuf)[i]) + intercept);  
    }
    break;
  case DT_COMPLEX:
    if (!doscaling) { 
      if (strcmp(FDTS,"double")!=0)
	{ 
	  for(i=0;i<size;i++) image->i[i] = (FDT) ((double*)tmpbuf)[i]; 
	}
    } else { 
      for(i=0;i<size;i++) 
	image->i[i] = (FDT) ((slope * ((double*)tmpbuf)[i]) + intercept);  
    }
    break;
  case DT_DOUBLE:
    if (!doscaling) { 
      if (strcmp(FDTS,"double")!=0)
	{ 
	  for(i=0;i<size;i++) image->i[i] = (FDT) ((double*)tmpbuf)[i]; 
	}
    } else { 
      for(i=0;i<size;i++) 
	image->i[i] = (FDT) ((slope * ((double*)tmpbuf)[i]) + intercept);  
    }
    break;
  default:
    { fprintf(stderr,"Wrong input data type for program\n"); exit(1); }
    break;
  }
  
  /* set correct data type - needed in case type was changed */
  if ( (strcmp(FDTS,"unsigned char")==0) || (strcmp(FDTS,"unsignedchar")==0) )
    { image->dt=DT_UNSIGNED_CHAR; image->dtmin=0;           image->dtmax=255;      }
  if ( (strcmp(FDTS,"signed short")==0) || (strcmp(FDTS,"signedshort")==0) )
    { image->dt=DT_SIGNED_SHORT;  image->dtmin=SHRT_MIN; image->dtmax=SHRT_MAX; }
  if ( (strcmp(FDTS,"signed int")==0) || (strcmp(FDTS,"signedint")==0) )
    {  image->dt=DT_SIGNED_INT;    image->dtmin=INT_MIN;   image->dtmax=INT_MAX;   }
  if (strcmp(FDTS,"float")==0)
    { image->dt=DT_FLOAT;         image->dtmin=-1e10;       image->dtmax=1e10;     }
  if (strcmp(FDTS,"double")==0)
    { image->dt=DT_DOUBLE;    image->dtmin=-1e10;       image->dtmax=1e10;     }
  
  image->bpv=sizeof(FDT);
  
  if (image->tr<1e-5) image->tr=3; /* well it's possible.... */
  image->thresh2=image->min;        /* better than nothing until it gets set properly using find_thresholds */
  image->thresh98=image->max;       /* ditto */
  image->thresh=image->min;         /* can use as a flag for whether find_thresholds has been run */
  image->lthresh=image->dtmin;      /* relevant procs know not to use these unless changed from these defaults */
  image->uthresh=image->dtmax;
  
}


void avw_read(const char *filename, image_struct *image)
{
  avw_read_basic(filename,image,1);
}


void avw_read_hdr(const char *filename, image_struct *image)
{
  avw_read_basic(filename,image,0);
}



int fsloutputtype() 
{
  return FslGetEnvOutputType();
}


void avw_write(const char *filename, image_struct image)
{
  FSLIO *fslio;

  fslio = FslOpen(FslMakeBaseName(filename),"wb");

  FslCloneHeader(fslio,image.miscinfo);
  FslSetIntensityScaling(fslio,1.0,0.0);

  FslSetDim(fslio,image.x,image.y,image.z,image.t);
  FslSetVoxDim(fslio,image.xv0,image.yv0,image.zv0,image.tr);

  FslSetIntent(fslio, image.intent_code, image.intent_p1, image.intent_p2, 
	       image.intent_p3);

  FslSetAuxFile(fslio,image.lut);

  FslSetCalMinMax(fslio,image.min,image.max);
  FslSetDataType(fslio,image.dt);

  FslWriteAllVolumes(fslio,image.i);
  FslClose(fslio);

}


