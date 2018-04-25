#ifndef MY_ASSERT_H
#define MY_ASSERT_H

#include<iostream>

#define myAssert(b, msg) {if (!(b)) std::cout<<"ERROR: "<<msg<<"\n";}
//#define myAssert(b, msg) {}

#endif // MY_ASSERT_H
