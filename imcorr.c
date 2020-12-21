/****************************************************************************
imcorr - program to do cross-correlation for motion tracking
Put under revision control by Bruce Raup, 2000-01-26

$Id: imcorr.c,v 1.4 2000-03-02 13:05:04-07 braup Exp braup $
*****************************************************************************/

#include <stdio.h>
#include "correlate.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

/*
 *    subimage pulls out a window from one of the main (byte) images and
 *    returns it in a one dimensional float array.  Note that the first
 *    pixel is (0,0).
 */

void subimage(px,ln,winbuf,bnpx,bnls,image,inpx,inls)
int px,ln;
float winbuf[];
int bnpx,bnls;
unsigned char *image;
int inpx,inls;
{
int line,pix;
for(line=0;line<bnls;line++)
	for(pix=0;pix<bnpx;pix++)
		winbuf[line*bnpx+pix]=(float)image[ln*inpx+line*inpx+px+pix];
}


void readimage(infilename,linebuf,imagebuf,np,nl,xoff,yoff,xextent,yextent)
char *infilename;
unsigned char *linebuf,*imagebuf;
int np,nl,xoff,yoff,xextent,yextent;
{
int sfds,i,j,k,l;
unsigned char *inptr2;
fprintf(stderr,"xoff=%d yoff=%d xext=%d yext=%d\n",xoff,yoff,
                xextent,yextent);
fprintf(stderr,"opening %s, size %d by %d\n",infilename,np,nl);
sfds=open(infilename,O_RDONLY);
if(sfds==-1)
        {
        fprintf(stderr,"error opening image: %s\n",infilename);
        exit(-1);
        }

fprintf(stderr,"reading %s\n",infilename);

l=0;

inptr2=imagebuf;

if(yoff>0) lseek(sfds,yoff*np,SEEK_SET);

for(i=0;i<yextent;i++)
{
        if((j=read(sfds,linebuf,(unsigned)np))!=np)
        {
        fprintf(stderr,
		"too few bytes (%d) read from image file %s\n",j,infilename);
        perror("read failed");

        exit(-1);
        }

for(k=0;k<xextent;k++)
       {
       l++;
       inptr2[k]=linebuf[xoff+k];
       }

inptr2+=xextent;
}

fprintf(stderr,"done reading image\n");
fprintf(stderr,"read %d pixels\n",l);
close(sfds);
return;
}


main(argc,argv)
int argc;
char **argv;
{
int np,nl,i,j,imgxoff,imgyoff,imgxextent,imgyextent,numpix,numlines;
unsigned char *simg,*rimg;

char version[] = "1.1";

float ref[16385];			/* Reference subimage (max 128)*/
float srch[65537];			/* Search subimage (max 256)*/

float mincorr;			/* Minimum acceptable correlation strength */
float maxdis;			/* Maximum displacement from nominal allowed */
float disp;			/* Actual displacement from nominal */
float est_err[3];		/* Estimated correlation errors */
float best_fit[2];		/* Best fit offsets from correlation peak */
float npls[2];			/* Size of search subimage */
float nplr[2];			/* Size of reference subimage */
float ioffrq[2];		/* Requested max horiz & vert search offsets */
float nomoff[2];		/* Nominal offsets */
float strength;			/* Strength of correlation */
int fitmeth=1;
int okparam=1;
int gridint;
int lmargin,rmargin,tmargin,bmargin,xoff,yoff;

unsigned char linebuf[10000];

FILE *outfile;


np=atoi(argv[3]);
nl=atoi(argv[4]);

switch(argc)
{
case(11):
	simg=(unsigned char *)malloc(np*nl);
	rimg=(unsigned char *)malloc(np*nl);
	npls[0] = (float)atoi(argv[6]);
	npls[1] = (float)atoi(argv[6]);
	nplr[0] = (float)atoi(argv[7]);
	nplr[1] = (float)atoi(argv[7]);
	mincorr = 2;
	maxdis = -1;
	ioffrq[0] = (float)atoi(argv[6])/2+fabs((float)atoi(argv[9]));
	ioffrq[1] = (float)atoi(argv[6])/2+fabs((float)atoi(argv[10]));
	nomoff[0] = 0.5*(npls[0]-nplr[0])+(float)atoi(argv[9]);
	nomoff[1] = 0.5*(npls[1]-nplr[1])+(float)atoi(argv[10]);
	gridint=atoi(argv[8]);
	xoff=atoi(argv[9]);
	yoff=atoi(argv[10]);
	numpix=np;
	numlines=nl;
	imgxoff=0;
	imgyoff=0;
	break;

case(6):		/* set up defaults */
	simg=(unsigned char *)malloc(np*nl);
	rimg=(unsigned char *)malloc(np*nl);
	npls[0] = 128;
	npls[1] = 128;
	nplr[0] = 32;
	nplr[1] = 32;
	mincorr = 2;
	maxdis = -1;
	ioffrq[0] = 64;
	ioffrq[1] = 64;
	nomoff[0] = 47;
	nomoff[1] = 47;
	gridint=25;
	xoff=0;
	yoff=0;
	numpix=np;
	numlines=nl;
	imgxoff=0;
	imgyoff=0;
	break;

case(15):
	imgxoff=atoi(argv[11]);
	imgyoff=atoi(argv[12]);
	imgxextent=atoi(argv[13]);
	imgyextent=atoi(argv[14]);
	simg=(unsigned char *)malloc(imgxextent*imgyextent);
	rimg=(unsigned char *)malloc(imgxextent*imgyextent);
	npls[0] = (float)atoi(argv[6]);
	npls[1] = (float)atoi(argv[6]);
	nplr[0] = (float)atoi(argv[7]);
	nplr[1] = (float)atoi(argv[7]);
	mincorr = 2;
	maxdis = -1;
	ioffrq[0] = (float)atoi(argv[6])/2+fabs((float)atoi(argv[9]));
	ioffrq[1] = (float)atoi(argv[6])/2+fabs((float)atoi(argv[10]));
	nomoff[0] = 0.5*(npls[0]-nplr[0])+(float)atoi(argv[9]);
	nomoff[1] = 0.5*(npls[1]-nplr[1])+(float)atoi(argv[10]);
	gridint=atoi(argv[8]);
	xoff=atoi(argv[9]);
	yoff=atoi(argv[10]);
	numpix=imgxextent;
	numlines=imgyextent;
	break;

default:
        printf("This is imcorr, release %s\n\n", version);
	printf("Usage:\n imcorr refimg searchimg pixels lines outfile\n");
	printf("or\n imcorr refimg searchimg pixels lines outfile");
	printf("\n search_box_size ref_box_size grid_space x_off y_off\n");
	printf("or\n imcorr refimg searchimg pixels lines outfile");
	printf("\n search_box_size ref_box_size grid_space x_off y_off ");
	printf("\n image_x_off image_y_off image_x_extent image_y_extent\n");
	exit(-1);
	break;
}

outfile=fopen(argv[5],"w");
if(outfile==NULL)
	{
	printf("error: can't open %s\n",argv[5]);
	exit(-1);
	}

#ifdef COMMANDLINE
/* print command line to outfile */
  fprintf( outfile, "# %s", argv[0] );
  for (i = 1; i < argc; ++i)
    fprintf( outfile, " %s", argv[i] );
  fprintf( outfile, "\n" );
#endif

printf("opening %s as ref img and %s as search img, size %d by %d\n",
	argv[1],argv[2],np,nl);
printf("writing to %s\n",argv[5]);

readimage(argv[1],linebuf,rimg,np,nl,imgxoff,imgyoff,numpix,numlines);
readimage(argv[2],linebuf,simg,np,nl,imgxoff,imgyoff,numpix,numlines);



/* Perform grey level (cross) correlation
  --------------------------------------*/
lmargin=(int)(npls[0]/2 + (float)xoff)<nplr[0]/2?
			nplr[0]/2:(int)(npls[0]/2 + (float)xoff);
rmargin=(int)(npls[0]/2 - (float)xoff)<nplr[0]/2?
			nplr[0]/2:(int)(npls[0]/2 - (float)xoff);
tmargin=(int)(npls[1]/2 + (float)yoff)<nplr[1]/2?
			nplr[1]/2:(int)(npls[1]/2 + (float)yoff);
bmargin=(int)(npls[1]/2 - (float)yoff)<nplr[1]/2?
			nplr[1]/2:(int)(npls[1]/2 - (float)yoff);

printf("l %d r %d t %d b %d\n",lmargin,rmargin,tmargin,bmargin);

for(i=lmargin;i<=numpix-rmargin;i+=gridint)
	for(j=tmargin;j<=numlines-bmargin;j+=gridint)
	{

printf("i %d j %d gridint %d\n",i,j,gridint);

	subimage(i-(int)(nplr[0]/2.0),j-(int)(nplr[1]/2.0),ref,
		(int)nplr[0],(int)nplr[1],rimg,numpix,numlines);
	subimage((i-xoff)-(int)(npls[0]/2.0),(j-yoff)-(int)(npls[1]/2.0),srch,
		(int)npls[0],(int)npls[1],simg,numpix,numlines);

	disp=0;
	strength=0;
	best_fit[0]=0;
	best_fit[1]=0;
	est_err[0]=0;
	est_err[1]=0;

	gcorr_(srch, ref, npls, nplr, &mincorr, &fitmeth, &maxdis, ioffrq,
	      nomoff, &okparam, &strength, best_fit, est_err, &disp);

	printf("%4d %4d ",i+imgxoff,j+imgyoff);

	printf("%8.3f %8.3f %4d ",disp,strength,okparam);

	printf("%8.3f %8.3f ",best_fit[0],best_fit[1]);
	printf("%8.3f %8.3f\n",est_err[0],est_err[1]);

	fprintf(outfile,"%4d %4d ",i+imgxoff,j+imgyoff);

	fprintf(outfile,"%8.3f %8.3f %4d ",disp,strength,okparam);

	fprintf(outfile,"%8.3f %8.3f ",best_fit[0],best_fit[1]);
	fprintf(outfile,"%8.3f %8.3f\n",est_err[0],est_err[1]);
	}
}
