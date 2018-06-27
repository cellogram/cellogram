#include<fstream>

#include "my_assert.h"
#include "grid.h"
#include "mesh.h"

namespace cellogram
{
namespace PointsUntangler
{

bool Mesh::importFVFix(const std::string& fnV , const std::string& fnF, const std::string& fnFix){
    std::fstream f;

    std::cout<<"<-- Import Fix from: "<<fnV<<": ";
    f.open(fnV.data(),std::fstream::in);
    assert(f.is_open()&&"cannot open Vert file");
    V.clear();
    while (1){
        scalar x,y,z;
        f>>x>>y>>z;
        assert( (z==0)&&"Z not 0");
        if (f.eof()) break;
        Vert v; v.p = vec2(x,y);
        V.push_back(v);
    }
    f.close();

    f.open(fnF.data(),std::fstream::in);
    assert(f.is_open()&&"cannot open Face file");
    F.clear();
    while (1){
        Face ff;
        ff.dontcare = false;
        f>> ff.vi[0]>>ff.vi[1]>>ff.vi[2];
        if (f.eof()) break;
        F.push_back(ff);
    }
    f.close();

    f.open(fnFix.data(),std::fstream::in);
    assert(f.is_open()&&"cannot open Fix file");
    for (Face& ff:F){
        int i;
        f>>i;
        assert(!f.eof());
        ff.fixed = (i==1);
    }
    f.close();

    buildEdgesFromFaces();
    propagateFixedF2E();
    updateValencies();

    std::cout<<"done ("<<V.size()<<"v/"<< F.size() << "f)\n";
    return true;

}

bool Mesh::importXYZv2(const std::string& filename ){
    std::cout<<"<-- Import XYZ from "<<filename<<": ";
    std::fstream f;
    f.open(filename.data(),std::fstream::in);
    if (!f.is_open()) {
        std::cout<<"Cannot open \""<<filename<<"\"\n";
        return false;
    }
    V.clear();
    while (1){
        scalar x,y;
        f>>x>>y;
        if (f.eof()) break;
        Vert v; v.p = vec2(x,y);
        V.push_back(v);
    }
    f.close();

    std::cout<<"done ("<<V.size()<<"v)\n";
    return true;
}

bool Mesh::importXYZv3(const std::string& filename ){
    std::cout<<"<-- Import XYZ from "<<filename<<": ";
    std::fstream f;
    f.open(filename.data(),std::fstream::in);
    if (!f.is_open()) {
        std::cout<<"Cannot open \""<<filename<<"\"\n";
        return false;
    }
    V.clear();
    while (1){
        scalar x,y;
        scalar dummy;
        f>>x>>y >> dummy;
        if (f.eof()) break;
        Vert v; v.p = vec2(x,y);
        V.push_back(v);
    }
    f.close();

    std::cout<<"done ("<<V.size()<<"v)\n";
    return true;
}
bool Mesh::importXYZ(const std::string& filename ){

    std::cout<<"<-- Import XYZ from "<<filename<<": ";
    std::fstream f;
    f.open(filename.data(),std::fstream::in);
    if (!f.is_open()) {
        std::cout<<"Cannot open \""<<filename<<"\"\n";
        return false;
    }
    int n;
    f>>n;
    V.resize(n);
    //vert.resize(n);

    for (int i=0; i<n; i++) {

        double dummy;
        f>> V[i].p.x >> V[i].p.y >> dummy;
        myAssert(dummy==0,"Error reading???");
    }
    f.close();

    std::cout<<"Done ("<<n<<"v)\n";

    return true;

}

bool Mesh::exportOBJ(const std::string& filename ){
    std::cout<<"--> Exporting: "<<filename<<"\n";

    std::ofstream f;
    f.open(filename.data(),std::fstream::out);
    if (!f.is_open()) return false;

    for (const Vert& v:V) f<<"v "<<v.p.x<<" "<<v.p.y<<" 0\n";
    for (const Face& ff:F) if (!ff.dontcare) f<<"f "<<ff[0]+1<<" "<<ff[1]+1<<" "<<ff[2]+1<<"\n";

    f.close();
    std::cout<<"Done writing OBJ ("<<V.size()<<"v/"<<F.size()<<"f)\n";
    return true;
}


bool Mesh::exportOFF(const std::string& filename ){

    std::cout<<"--> Exporting: "<<filename<<"\n";

    std::ofstream f;
    f.open(filename.data(),std::fstream::out);
    if (!f.is_open()) return false;

    for (const Vert& v:V) f<<"v "<<v.p.x<<" "<<v.p.y<<" 0\n";
    for (const Face& ff:F) if (!ff.dontcare) f<<"f "<<ff[0]+1<<" "<<ff[1]+1<<" "<<ff[2]+1<<"\n";

    f.close();
    std::cout<<"Done writing OBJ ("<<V.size()<<"v/"<<F.size()<<"f)\n";
    return true;
}

static void setColorByValency(int&r, int &g, int &b, int val){
    r = g = b = 255;
    if (val>6) {r/=2; g/=2;}
    if (val>7) {r/=2; g/=2;}
    if (val<6) {b/=2; g/=2;}
    if (val<5) {b/=2; g/=2;}
}

static void setColorByFloodfill(int&r, int &g, int &b, float timeReached, float disputed){
    r = g = b = 255;
    float rr,gg,bb;
    rr = gg = timeReached;
    bb = 1.0;

    gg *= 1.0-disputed;
    bb *= 1.0-disputed;
    r = int(rr*255);
    g = int(gg*255);
    b = int(bb*255);
}

static void setColorByDisputed(int&r, int &g, int &b, float disputed){
    float rr,gg,bb;
    rr = gg = bb = 1.0;
    gg *= 1.0-disputed*2;
    bb *= 1.0-disputed*2;
    r = int(rr*255);
    g = int(gg*255);
    b = int(bb*255);

}

bool Mesh::exportPLY(const std::string& filename, ColMode colMode){
    std::ofstream f;
    f.open((filename+".ply").data(),std::fstream::out);
    std::cout<<"--> Exporting: "<<(filename+".ply")<<": ";
    if (!f.is_open()) return false;

    int nf =0;
    for (const Face& ff:F) if (!ff.dontcare) nf++;
    f
    <<"ply\n"
    <<"format ascii 1.0\n"
    <<"element vertex "<< V.size()<<"\n"
    <<"property float x\n"
    <<"property float y\n"
    <<"property float z\n"
    <<"property uchar red\n"
    <<"property uchar green\n"
    <<"property uchar blue\n"
    <<"property uchar alpha\n"
    <<"element face "<< nf<<"\n"
    <<"property list uchar int vertex_indices\n"
    <<"end_header\n";

    for (const Vert& v:V) {
        int r,g,b;
        r=g=b=255;
        switch (colMode){
        case BY_FLOOD: setColorByFloodfill(r,g,b,v.timeReached,v.disputed); break;
        case BY_VAL: setColorByValency(r,g,b,v.val); break;
        case BY_DISPUTED: setColorByDisputed(r,g,b,v.disputed); break;
        }
        f<<v.p.x<<" "<<v.p.y<<" 0 "<<r<<" "<<g<<" "<<b<<" "<<" 255\n";
    }
    for (const Face& ff:F) if (!ff.dontcare) f<<"3 "<<ff[0]<<" "<<ff[1]<<" "<<ff[2]<<"\n";

    f.close();
    std::cout<<"...done ("<<V.size()<<"/v"<<F.size()<<"/f)\n";
    return true;
}

bool Mesh::exportEdgesPLY(const std::string& filename){
    std::cout<<"--> Exporting (Edges): "<<(filename+".ply")<<": ";
    std::ofstream f;
    f.open(filename.data(),std::fstream::out);

    if (!f.is_open()) return false;

    f
    <<"ply\n"
    <<"format ascii 1.0\n"
    <<"element vertex "<< (E.size()*3)<<"\n"
    <<"property float x\n"
    <<"property float y\n"
    <<"property float z\n"
    <<"property uchar red\n"
    <<"property uchar green\n"
    <<"property uchar blue\n"
    <<"property uchar alpha\n"
    <<"element face "<< E.size() <<"\n"
    <<"property list uchar int vertex_indices\n"
    <<"end_header\n";

    for (const Edge& e:E) {
        int r,g,b;
        if (!e.fixed) { r = b = 0; g = 200;} else r=g=b=255;

        auto pa = V[e.vi[0]].p;
        f<<pa.x<<" "<<pa.y<<" 0 "<<r<<" "<<g<<" "<<b<<" "<<" 255\n";
        f<<pa.x<<" "<<pa.y<<" 0 "<<r<<" "<<g<<" "<<b<<" "<<" 255\n";

        auto pb = V[e.vi[1]].p;
        f<<pb.x<<" "<<pb.y<<" 0 "<<r<<" "<<g<<" "<<b<<" "<<" 255\n";

    }

    for (uint i=0; i<E.size(); i++) f<<"3 "<<(i*3)<<" "<<(i*3+1)<<" "<<(i*3+2)<<"\n";

    f.close();
    std::cout<<"... done\n";
    return true;
}




void setColormap(float d, int &r,int &g, int &b){
    if (d==-1) { r = 0.0; g=200; b=0; return;}
    r = 255;
    d*=.1;
    if (d>1) d = 1;
    if (d<0) d = 0;
    g = (1-d) * 255;
    b = (1-d) * 255;
}

extern Grid g;
bool Mesh::exportEdgesPLY(const std::string& filename, const Grid& grid){
    std::ofstream f;

    std::cout<<"--> Exporting (Edges + Grid): "<<(filename+".ply")<<": ";

    f.open((filename+".ply").data(),std::fstream::out);
    if (!f.is_open()) return false;

    f
    <<"ply\n"
    <<"format ascii 1.0\n"
    <<"element vertex "<< (E.size()*3)<<"\n"
    <<"property float x\n"
    <<"property float y\n"
    <<"property float z\n"
    <<"property uchar red\n"
    <<"property uchar green\n"
    <<"property uchar blue\n"
    <<"property uchar alpha\n"
    <<"element face "<< E.size() <<"\n"
    <<"property list uchar int vertex_indices\n"
    <<"end_header\n";

    for (const Edge& e:E) {
        int r,gg,b;

        setColormap( grid.misalignmentOptimist
                     ( e.vi[0], e.vi[1]) , r,gg,b);
        if (!e.fixed) { r = b = 0; gg = 200;} else r=gg=b=255;

        auto pa = V[e.vi[0]].p;
        f<<pa.x<<" "<<pa.y<<" 0 "<<r<<" "<<gg<<" "<<b<<" "<<" 255\n";
        f<<pa.x<<" "<<pa.y<<" 0 "<<r<<" "<<gg<<" "<<b<<" "<<" 255\n";

        auto pb = V[e.vi[1]].p;
        f<<pb.x<<" "<<pb.y<<" 0 "<<r<<" "<<gg<<" "<<b<<" "<<" 255\n";

    }

    for (uint i=0; i<E.size(); i++) f<<"3 "<<(i*3)<<" "<<(i*3+1)<<" "<<(i*3+2)<<"\n";

    f.close();
    std::cout<<"...done\n";
    return true;
}



}} // namespaces
