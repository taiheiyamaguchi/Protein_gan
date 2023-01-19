#include "md.h"

#ifndef _MACROS_H
#define _MAXROS_H

#define   ALLOCATE_INT_ARRAY(x,n) \
	{ if((n) > 0)\
	    x=new int[(n)];\
	  else\
	    x=NULL;\
	}

#define   ALLOCATE_DOUBLE_ARRAY(x,n) \
	{ if((n) > 0)\
	    x=new double[(n)];\
	  else\
	    x=NULL;\
	}

#define   ALLOCATE_INT_MATRIX(x,n,m) \
	{ if((n) > 0 && (m) > 0) {\
	    x=new int*[(n)];\
	    x[0]=new int[(n)*(m)];\
	    for(int i=1;i<(n);i++) {\
	      x[i]=&x[0][(m)*i];\
	    }\
	  } else {\
	    x=NULL;\
	  }\
	}

#define   ALLOCATE_DOUBLE_MATRIX(x,n,m) \
	{ if((n) > 0 && (m) > 0) {\
	    x=new double*[(n)];\
	    x[0]=new double[(n)*(m)];\
	    for(int i=1;i<(n);i++) {\
	      x[i]=&x[0][(m)*i];\
	    }\
	  } else {\
	    x=NULL;\
	  }\
	}

#define   FREE_ARRAY(x) \
	{ if((x) != NULL) delete [] (x); }

#define  FREE_MATRIX(x,n)      \
	{ if(x != NULL) {\
	    delete [] x[0];\
	    delete [] x;\
	  }\
	}

#define   SQR(x)	((x)*(x))

#define   PI		3.1415926535897932384626433832795

#define   BOLTZ		1.98719168260038e-3

#endif
