#include <float.h>
#include <math.h>

#include "Hybrid_Taus.h"

unsigned int Taus_step( unsigned int &z, int S1, int S2, int S3, unsigned int M )
{
	return z = ( ( ( z & M ) << S3 ) ^ ( ( ( z << S1 ) ^ z ) >> S2 ) );
}

unsigned int LCG_step( unsigned int &z, unsigned int a, unsigned int c )
{
	return z = ( a*z + c );
}

double hybrid_Taus( TAUS_SEED &seed )
{
	return 2.3283064365387e-10*( Taus_step( seed.z1, 13, 19, 12, 4294967294UL )
		^ Taus_step( seed.z2, 2, 25, 4, 4294967288UL )
		^ Taus_step( seed.z3, 3, 11, 17, 4294967280UL )
		^ LCG_step( seed.z4, 1664525, 1013904223UL ) );
}

unsigned int hybrid_Taus_int( TAUS_SEED &seed )
{
	return Taus_step( seed.z1, 13, 19, 12, 4294967294UL )
		^ Taus_step( seed.z2, 2, 25, 4, 4294967288UL )
		^ Taus_step( seed.z3, 3, 11, 17, 4294967280UL )
		^ LCG_step( seed.z4, 1664525, 1013904223UL );
}
