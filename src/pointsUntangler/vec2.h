#ifndef VEC2_H
#define VEC2_H

#include<math.h>

namespace cellogram
{
namespace PointsUntangler
{

typedef double scalar;

typedef unsigned int uint;

struct vec2{
    scalar x,y;
    vec2(scalar _x,scalar _y):x(_x),y(_y){}
    vec2():x(0.0f),y(0.0f){}
    vec2(const scalar p[2]):x(p[0]),y(p[1]){}

    vec2 operator +(const vec2& b) const{return vec2(x+b.x,y+b.y);}
    vec2 operator -(const vec2& b) const{return vec2(x-b.x,y-b.y);}
    vec2 operator /(scalar f) const{return vec2(x/f,y/f);}
    vec2 operator *(scalar f) const{return vec2(x*f,y*f);}
    void operator -=(const vec2& b) {x-=b.x;y-=b.y;}
    void operator +=(const vec2& b) {x+=b.x;y+=b.y;}
    void operator /=(scalar f) {x/=f;y/=f;}
    void operator *=(scalar f) {x*=f;y*=f;}
    bool operator ==(const vec2& b) const { return (x==b.x)&&(y==b.y);}
    void normalize(){ scalar n = norm(); if (n!=0) {x/=n; y/=n;} }
    scalar norm() const{ return sqrt(x*x+y*y);}
};

inline static scalar squaredDistance(vec2 a, vec2 b){ return (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y);}
inline static scalar distance(vec2 a, vec2 b){ return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y)); }
inline static scalar cross(vec2 a, vec2 b){ return (a.x)*(b.y) - (a.y)*(b.x); }
inline static scalar dot(vec2 a, vec2 b){ return (a.x)*(b.x) + (a.y)*(b.y); }

inline static vec2 rotateBy60( vec2 a){
    scalar s = 0.86602540378;
    scalar c = 0.5;
    return vec2(a.x*c - a.y*s ,  a.x*s + a.y*c );
}

}
}

#endif // VEC2_H
