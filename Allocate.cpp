// Include standard C libraries
#include <cstdlib>
#include <cstdio>

// Include header file
#include "Allocate.h"

void *allocate( unsigned int size )
{
	void *ptr = malloc(size);
	if(ptr == NULL)
	{
		fprintf(stderr,"Memory allocation error, terminating\n");
		exit(EXIT_FAILURE);
	}
	return(ptr);
}
