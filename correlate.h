#include <stdio.h>

#define STAT_FLAG status
#define CHECK(ARG) status = ARG;

#define GREY_CORR	1
#define EDGE_CORR	2

#define CMLEN 80

#define PARAB	1
#define GAUSS	2
#define RECIP	3
#define INTEGER	4

#define REF_MIN		32
#define REF_MAX	       128
#define REF_BUFSIZE  16384
#define SRCH_MIN	48 
#define SRCH_MAX       256	
#define SRCH_BUFSIZE 65536 

/* Error messages.  Additional error messages may be found in geompak.h
  --------------------------------------------------------------------*/
#define IN_BATCH	201
#define NO_CORR		202
#define CORR_METHOD	203
#define BAD_SIZE	204
#define OUTSIDE_IMAGE	205
#define READ_ERROR	206
#define BAD_COEFS	207
#define BAD_CHIP_DDR	208
#define NO_CHIPS	209
#define BAD_PTID	210
