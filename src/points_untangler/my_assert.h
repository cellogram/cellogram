#ifndef MY_ASSERT_H
#define MY_ASSERT_H

#include<iostream>

#ifdef NDEBUG
#define myAssert(b, msg) {}
#else
#define myAssert(b, msg) {if (!(b)) std::cout<<"ERROR: "<<msg<<"\n"; assert(b);}
#endif

#endif // MY_ASSERT_H
