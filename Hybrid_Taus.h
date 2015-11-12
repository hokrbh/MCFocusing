#ifndef HYBRID_TAUS_H
#define HYBRID_TAUS_H

typedef struct
{
	unsigned int z1, z2, z3, z4;
} TAUS_SEED;

// Function definitions
unsigned int Taus_step( unsigned int &z, int S1, int S2, int S3, unsigned int M );
unsigned int LCG_step( unsigned int &z, unsigned int a, unsigned int c );
double hybrid_Taus( TAUS_SEED &seed );
unsigned int hybrid_Taus_int( TAUS_SEED &seed );

#endif
