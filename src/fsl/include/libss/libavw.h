#ifndef _LIBAVW_H_
#define _LIBAVW_H_

#include "fslio/fslio.h"

void avw_read(const char*, image_struct*);
void avw_write(const char*, image_struct);

void avw_read_hdr(const char*, image_struct*);
int fsl_imageexists(const char *filename);
int fsloutputtype();

#endif
