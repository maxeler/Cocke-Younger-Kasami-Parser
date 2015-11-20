#include <inttypes.h> // defines uint32_t

//typedef unsigned int bitarray_t; // if you know that int is 32 bits
typedef uint32_t bitarray_t;

#define RESERVE_BITS(n) (((n)+0x1f)>>5)
#define DW_INDEX(x) ((x)>>5)
#define BIT_INDEX(x) ((x)&0x1f)
#define getbit(array,index) (((array)[DW_INDEX(index)]>>BIT_INDEX(index))&1)
#define putbit(array, index, bit) \
    ((bit)&1 ?  ((array)[DW_INDEX(index)] |= 1<<BIT_INDEX(index)) \
             :  ((array)[DW_INDEX(index)] &= ~(1<<BIT_INDEX(index))) \
             , 0 \
    )

inline int TEST(bitarray_t* arr1, bitarray_t* arr2, int numbits){
    uint32_t result=0;
    int i;
    for ( i=0; i < RESERVE_BITS(numbits); i++ )
	result |= ( arr1[i] & arr2[i] );
    return result;
}

inline void OR(bitarray_t* arr1, bitarray_t* arr2, int numbits){
    int i;
    for ( i=0; i < RESERVE_BITS(numbits); i++ )
	arr1[i] |= arr2[i];
}

inline void COPY(bitarray_t* arr1, bitarray_t* arr2, int numbits){
    int i;
    for ( i=0; i < RESERVE_BITS(numbits); i++ )
	arr1[i] = arr2[i];
}

inline void CLEAR(bitarray_t* arr1, int numbits){
    int i;
    for ( i=0; i < RESERVE_BITS(numbits); i++ )
	arr1[i] = 0;
}

/* Use:
bitarray_t arr[RESERVE_BITS(130)] = {0, 0x12345678,0xabcdef0,0xffff0000,0};
int i = getbit(arr,5);
putbit(arr,6,1);
int x=2;            // the least significant bit is 0
putbit(arr,6,x);    // sets bit 6 to 0 because 2&1 is 0
putbit(arr,6,!!x);  // sets bit 6 to 1 because !!2 is 1
*/
