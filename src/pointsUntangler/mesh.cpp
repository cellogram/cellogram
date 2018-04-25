#include<fstream>
#include<set>
#include<map>

#include "my_assert.h"
#include "mesh.h"

namespace cellogram
{
namespace PointsUntangler
{

void Mesh::setFaceRegularity(){
    for (Face &f:F) {
        f.regularity =  goodTriangle( f[0],f[1],f[2] )*0.5;
        for (int w=0; w<3; w++)
            f.regularity += V[ f[w] ].distToIrr;
        //std::cout<<"REG = "<<f.regularity<<" = "<<V[ f[0] ].distToIrr<<","<<V[ f[1] ].distToIrr<<","<<V[ f[2] ].distToIrr<<"\n";
    }
}

void Mesh::smoothDisputed(){
    std::vector<scalar> old(V.size());
    std::vector<int> div(V.size());
    for (uint i=0; i<V.size(); i++) { old[i]=V[i].disputed; div[i]=0; V[i].disputed=0;}
    for (const Face& f:F)
    {
        V[f[0]].disputed+=old[f[1]];
        V[f[0]].disputed+=old[f[2]];
        V[f[1]].disputed+=old[f[0]];
        V[f[1]].disputed+=old[f[2]];
        V[f[2]].disputed+=old[f[0]];
        V[f[2]].disputed+=old[f[1]];
        div[f[0]]+=2;
        div[f[1]]+=2;
        div[f[2]]+=2;
    }
    for (uint i=0; i<V.size(); i++) { if (div[i]) V[i].disputed/=div[i];}


}

int Mesh::mostRegularFace()const{
    int res = -1;
    scalar best = -1000000;
    for (int i=0; i<(int)F.size(); i++) {
        //std::cout<<"face "<<i<<": "<<F[i].regularity<<"\n";
        if (F[i].regularity>best) {
            best = F[i].regularity;
            res = i;
        }
    }
    myAssert(res!=-1,"Cannot find face!");
    //std::cout<<"BEST is face "<<res<<" with "<<best<<"("<<goodTriangle(F[res][0],F[res][1],F[res][2])<<"\n";
    return res;
}

void Mesh::setDistanceToIrr(){
    for (Vert &v:V) v.distToIrr = ((v.val!=6) || v.dontcare)? 0 : 1000;


    while (1){
        bool over = true;
        for (const Edge &e:E){
            int minDist = std::min( V[ e[0] ].distToIrr, V[ e[1] ].distToIrr );
            minDist++;
            if (V[ e[0] ].distToIrr > minDist) { V[ e[0] ].distToIrr=minDist;over = false;}
            if (V[ e[1] ].distToIrr > minDist) { V[ e[1] ].distToIrr=minDist;over = false;}
        }
        if (over) break;
    }
}

void Mesh::updateAverageEdge(){
    avgEdge = 0;
    scalar d = 0;
    for (const Edge &e:E) {
        if ( V[e[0]].dontcare ) continue;
        if ( V[e[1]].dontcare ) continue;
        avgEdge += distance(V[e[0]].p , V[e[1]].p);
        d++;
    }
    avgEdge /= d;
}

void Edge::substitute(int fa, int fb){
        if (fi[0]==fa) fi[0]=fb;
        else if (fi[1]==fa) fi[1]=fb;
        else myAssert(false, "SUBSTITUTE NO\n");
}


bool Mesh::importXYZv2(const std::string& filename ){
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

    return true;
}

void Mesh::fromEigen(const Eigen::MatrixXd &mat){
    for(int i = 0; i < mat.rows(); ++i){
        Vert v; v.p = vec2(mat(i,0),mat(i,1));
        V.push_back(v);
    }
}

bool Mesh::importXYZ(const std::string& filename ){

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

    std::cout<<"Done reading "<<n<<" verts\n";

    return true;

}

void Mesh::fuckUp(){
    //Vert v; v.p = vec2(200,600); V.push_back(v);
    V[2420].p = vec2(200,600);

    std::cout<<"Fucked up mesh :D ! ";
}

bool Mesh::exportOBJ(const std::string& filename ){
    std::ofstream f;
    f.open(filename.data(),std::fstream::out);
    if (!f.is_open()) return false;

    for (const Vert& v:V) f<<"v "<<v.p.x<<" "<<v.p.y<<" 0\n";
    for (const Face& ff:F) if (!ff.dontcare) f<<"f "<<ff[0]+1<<" "<<ff[1]+1<<" "<<ff[2]+1<<"\n";

    f.close();
    std::cout<<"Done writing OBJ ("<<V.size()<<" verts, "<<F.size()<<" faces)\n";
    return true;
}


bool Mesh::exportOFF(const std::string& filename ){
    std::ofstream f;
    f.open(filename.data(),std::fstream::out);
    if (!f.is_open()) return false;

    for (const Vert& v:V) f<<"v "<<v.p.x<<" "<<v.p.y<<" 0\n";
    for (const Face& ff:F) if (!ff.dontcare) f<<"f "<<ff[0]+1<<" "<<ff[1]+1<<" "<<ff[2]+1<<"\n";

    f.close();
    std::cout<<"Done writing OBJ ("<<V.size()<<" verts, "<<F.size()<<" faces)\n";
    return true;
}

static void setColorByValency(int&r, int &g, int &b, int val){
    r = g = b = 255;
    if (val>6) {r/=2; g/=2;}
    if (val>7) {r/=2; g/=2;}
    if (val<6) {b/=2; g/=2;}
    if (val<5) {b/=2; g/=2;}
}

static void setColorByFloodfill(int&r, int &g, int &b, float timeReached, int disputed){
    r = g = b = 255;
    if (disputed>0) g=b=std::max(0,255-disputed*75);
    else {
        r = g = timeReached*255;
    }

}


bool Mesh::exportPLY(const std::string& filename, bool colorByFloodfill){
    std::ofstream f;
    f.open(filename.data(),std::fstream::out);
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
        if (colorByFloodfill) setColorByFloodfill(r,g,b,v.timeReached,v.disputed);
        else setColorByValency(r,g,b,v.val);
        f<<v.p.x<<" "<<v.p.y<<" 0 "<<r<<" "<<g<<" "<<b<<" "<<" 255\n";
    }
    for (const Face& ff:F) if (!ff.dontcare) f<<"3 "<<ff[0]<<" "<<ff[1]<<" "<<ff[2]<<"\n";

    f.close();
    std::cout<<"Done writing PLY ("<<V.size()<<" verts, "<<F.size()<<" faces)\n";
    return true;
}


void Mesh::updateValencies(){
    for (Vert& v:V) {v.val=0;}

    for (Edge& e:E) {
        V[e[0]].val++;
        V[e[1]].val++;
    }

    irregulars.clear();
    int nIrregular = 0;
    for (int i=0; i<(int)V.size(); i++) {
        if (V[i].dontcare) continue;
        int irr = V[i].val-6;
        if (irr!=0) irregulars.insert(i);
        nIrregular+=std::abs(irr);
    }
    //std::cout<<nIrregular<<" irregular vertices on "<<V.size()<<"\n";
    nVal = nIrregular;
}

typedef std::pair<int,int> intint;
intint unique(intint p){
    if (p.first>p.second) return intint(p.second,p.first);
    return p;
}

int Face::oppositeVertOfEdge(const Edge &e) const{
    int res = vi[0];
    if (res==e[0] || res==e[1]) res = vi[1];
    if (res==e[0] || res==e[1]) res = vi[2];
    return res;
}

void Mesh::updateIndices(){

    std::map< Edge , int > map;

    E.clear();
    for (int fi=0; fi<(int)F.size(); fi++){
        for (int w=0; w<3; w++){

            Edge e( F[fi][w],F[fi][(w+1)%3] );
            auto f = map.find(e);

            int ei;
            if (f==map.end()) {
                ei = E.size();
                map[e] = ei;
                E.push_back(e);
            }
            else {
                ei = f->second;
                myAssert(ei<(int)E.size(),"Strange");
                myAssert(ei>=0,"Strange");
            }

            F[fi].ei[w] = ei;
            myAssert(ei<(int)E.size(),"Strange");
            myAssert(ei>=0,"Strange");

            if (E[ei].fi[0]==-1) E[ei].fi[0]=fi;
            else if (E[ei].fi[1]==-1) E[ei].fi[1]=fi;
            else myAssert(false,"EDGE "<<ei<<" ("<<E[ei].vi[0]<<","<<E[ei].vi[1]<<") NON MANIFOLD! (face "<<fi<<" w"<<w<<")");

        }
    }

}

void Mesh::dontCareAboutBoundaries(){
    for (Vert& v:V) v.dontcare = false;
    for (const Edge& e:E)
    {
        if (e.fi[1]==-1) V[ e.vi[0] ].dontcare = V[ e.vi[1] ].dontcare = true;
    }
}

void Mesh::propagateDontcareToFaces(){

    for (Face& f:F) {
        f.dontcare = false;
        for (int w=0; w<3; w++)
        if (V[f[w]].dontcare) f.dontcare = true;
    }
}

void Mesh::propagateDontcareToVerts(){

    for (Vert& v:V) v.dontcare = false;
    for (Face& f:F)  {

        for (int w=0; w<3; w++) {
            if (f.dontcare) V[f[w]].dontcare = true;
        }
    }
}

bool Mesh::canFlip(int ei){
    if (ei==14122) return false; // small hack for now
    if (E[ei].fi[1]==-1) return false; // boundary
    if (E[ei].fi[0]==-1) return false; // boundary

    int v0 = E[ei][0];
    int v1 = E[ei][1];

    if (V[v0].val<4) return false;
    if (V[v1].val<4) return false;

    int fa = E[ei].fi[0];
    int fb = E[ei].fi[1];
    int va=-1,vb=-1;
    for (int w=0; w<3; w++) if ((F[fa][w]!=v0) && (F[fa][w]!=v1)) va = F[fa][w];
    for (int w=0; w<3; w++) if ((F[fb][w]!=v0) && (F[fb][w]!=v1)) vb = F[fb][w];
    if (va==vb) return false;

    if (V[v0].dontcare || V[v1].dontcare || V[va].dontcare || V[vb].dontcare) return false;

    return true;
}


// if v1 had to increas valency by delta,
// how much it would be affected by irregular vert v0?
scalar Mesh::deltaEng(int v0, int v1, int delta){
    //int val0 = V[v0].val-6;
    int val1 = V[v1].val-6 - delta;
    //myAssert(val0!=0,"wrong valency");
    float d = distance(V[v0].p,V[v1].p);
    //if (d==0) return 0;
    //return -delta*val0*val1/(10+d*d);
    return -fabs(val1)/(1+d*d);
}

scalar Mesh::goodTriangle(int v0,int v1, int v2) const{
    vec2 e0 = (V[v0].p-V[v1].p);
    vec2 e1 = (V[v1].p-V[v2].p);
    vec2 e2 = (V[v2].p-V[v0].p);
    scalar a0 = e0.norm();
    scalar a1 = e1.norm();
    scalar a2 = e2.norm();

    if (a0<a1) std::swap(a0,a1);
    if (a1<a2) std::swap(a1,a2);
    if (a0<a1) std::swap(a0,a1);
    return a2-a0;
    /*
    scalar avg = (a0+a1+a2)/3.0;
    //return -avg;
    scalar variance = (a0*a0+a1*a1+a2*a2)/3.0-avg*avg;
    return -variance * 10;
    */

    //e0.normalize();
    //e1.normalize();
    //e2.normalize();
    //return 3.0-(dot(e0,e1)+dot(e1,e2)+dot(e2,e0));
}

scalar myPow(scalar b, scalar exp){
    scalar res = 1;
    int i; for (i=0; i<exp; i++) res*=b;
    res*=(exp-i)*b+(1-exp+i);
    return res;
}

scalar Vert::price()const{
    return myPow(0.9,disputed);
    //return (disputed>0)?0.5:1.0;
}

FlipScore Mesh::evaluateFlip(int ei, bool desperate){
    int v0 = E[ei][0];
    int v1 = E[ei][1];
    int fa = E[ei].fi[0];
    int fb = E[ei].fi[1];

    int va=-1,vb=-1;
    for (int w=0; w<3; w++) if ((F[fa][w]!=v0) && (F[fa][w]!=v1)) va = F[fa][w];
    for (int w=0; w<3; w++) if ((F[fb][w]!=v0) && (F[fb][w]!=v1)) vb = F[fb][w];

    /*
    myAssert(v0!=v1,"same vert 01");
    myAssert(v0!=va,"same vert 0a");
    myAssert(v0!=vb,"same vert 0b");
    myAssert(v1!=va,"same vert 1a");
    myAssert(v1!=vb,"same vert 1b");
    myAssert(va!=vb,"same vert ab");
    myAssert(va>=0,"same vert a-1");
    myAssert(vb>=0,"same vert b-1");
    myAssert(v0>=0,"same vert 0-1");
    myAssert(v1>=0,"same vert 1-1");
    */

    scalar valReduction = 0;
    if (!V[v0].dontcare) { if (V[v0].val>6) valReduction+=V[v0].price(); else valReduction-=V[v0].price(); }
    if (!V[v1].dontcare) { if (V[v1].val>6) valReduction+=V[v1].price(); else valReduction-=V[v1].price(); }
    if (!V[va].dontcare) { if (V[va].val<6) valReduction+=V[va].price(); else valReduction-=V[va].price(); }
    if (!V[vb].dontcare) { if (V[vb].val<6) valReduction+=V[vb].price(); else valReduction-=V[vb].price(); }

    /*scalar mult01 = ((V[v0].disputed>0)||(V[v1].disputed>0))? 1.2 : 1;
    scalar multAB = ((V[va].disputed>0)||(V[vb].disputed>0))? 1.2 : 1;*/

    scalar lenReduction=0;
    lenReduction += (distance(V[v1].p,V[v0].p) - distance(V[va].p,V[vb].p))/10.0;

    if (desperate) {
        lenReduction += 10.0*(rand()%1000*0.001);
    }
    myAssert(va!=vb,"Dont!");

    return FlipScore(valReduction,lenReduction);
}


FlipScore lastMove(0,0);
std::vector<int> forbiden;

static bool contains(const std::vector<int>& v, int i){
    for (int j:v) if (j==i) return true; return false;
}

int Mesh::bestFlip(bool force ){
    int winner = -1;

    FlipScore bestScore(0,0);
    if (!force) bestScore = FlipScore(0,0);
    else bestScore = FlipScore(0,-99999);

    //else bestScore.second = -20.0;
    for (int i=0; i<(int)E.size(); i++) {
        if (!canFlip(i)) continue;

        auto score = evaluateFlip(i,force);
        if (force) {
            //if (FlipScore(0,0)<score) continue;
            if (contains(forbiden,i)) continue;
        }
        if (bestScore<score) {
            bestScore = score;
            winner = i;
        }
    }
    forbiden.push_back(winner);
    lastMove = bestScore;
    //std::cout<<"BestScore ("<<bestScore.first<<","<<bestScore.second<<")\n";
    return winner;
}



bool Mesh::sanityCheck(){
    for (int fi=0; fi<(int)F.size(); fi++){
        const Face& f (F[fi]);
        for (int w=0; w<3; w++) {
            int ei = f.ei[w];
            myAssert(ei>=0, "EI>0");
            myAssert(ei<(int)E.size(), "EI<E.size()");
            int v0 = f.vi[w];
            int v1 = f.vi[(w+1)%3];
            myAssert(v0!=v1,"The same F?");
            Edge & e( E[ei] );
            myAssert( ((e.vi[0] == v1) && (e.vi[1] == v0))
                   || ((e.vi[1] == v1) && (e.vi[0] == v0)), "EV");

        }
    }
    for (int ei=0; ei<(int)E.size(); ei++){
        const Edge& e (E[ei]);
        for (int s=0; s<2; s++) {
            int fi = e.fi[s];
            if (s==1) if (fi==-1) continue;
            myAssert(fi>=0, "FI>0");
            myAssert(fi<(int)F.size(), "fi<F.size()");
            int v0 = e.vi[0];
            int v1 = e.vi[1];
            myAssert(v0!=v1,"the same E?");
            Face & f( F[fi] );
            myAssert( ((f.vi[0] == v1) || (f.vi[1] == v1) || (f.vi[2] == v1)), "EF1");
            myAssert( ((f.vi[0] == v0) || (f.vi[1] == v0) || (f.vi[2] == v0)), "EF2");


        }
    }

    return true;

}

void Mesh::updateValence(int vi, int delta){
    V[vi].val+=delta;
    if (V[vi].val == 6) {
        nVal--;
        irregulars.erase(vi);
    } else  {
        nVal++;
        irregulars.insert(vi);
    }

}

void Mesh::applyFlip(int ei){
    //std::cout<<"FLIPPP\n";
    int v0 = E[ei][0];
    int v1 = E[ei][1];
    int fa = E[ei].fi[0];
    int fb = E[ei].fi[1];
    int wa,wb;

    for (int w=0; w<3; w++) if ((F[fa][w]!=v0) && (F[fa][w]!=v1)) wa = w;
    for (int w=0; w<3; w++) if ((F[fb][w]!=v0) && (F[fb][w]!=v1)) wb = w;

    //std::cout<<"ERA: "<<fa<<" ->"<<F[fa].ei[(wa+2)%3]             <<
    //              "  F = "<<E[ F[fa].ei[(wa+2)%3]].fi[0]<<","<<E[ F[fa].ei[(wa+2)%3]].fi[1]<<"\n";

    myAssert(F[fa].ei[(wa+1)%3]==ei,"SENNO!OOOOOOOOOOOOOOOOAAA");
    myAssert(F[fb].ei[(wb+1)%3]==ei,"SENNO!OOOOOOOOOOOOOOOOBBB");
    int va = F[fa][wa];
    int vb = F[fb][wb];

    myAssert(va!=vb,"CAZZO NON MANIFOLD VA E VB\n");

    // fix EV
    E[ei][0] = vb;
    E[ei][1] = va;

    // fix FV
    //myAssert(F[fa].vi[(wa+1)%3]==v1,"CAZZO\n");
    F[fa].vi[(wa+1)%3] = vb;
    F[fb].vi[(wb+1)%3] = va;

    // fix EF
    E[ F[fa].ei[wa] ].substitute(fa,fb);
    E[ F[fb].ei[wb] ].substitute(fb,fa);

    // fix FE
    F[fa].ei[(wa+1)%3] = F[fb].ei[wb];
    F[fb].ei[(wb+1)%3] = F[fa].ei[wa];
    F[fa].ei[wa] = ei;
    F[fb].ei[wb] = ei;

    // fix vert valencies
    updateValence(va,+1);
    updateValence(vb,+1);
    updateValence(v0,-1);
    updateValence(v1,-1);

    //std::cout<<"E': "<<fa<<" ->"<<F[fa].ei[(wa+2)%3]
    //         <<"   F = "<<E[ F[fa].ei[(wa+2)%3]].fi[0]<<","<<E[ F[fa].ei[(wa+2)%3]].fi[1]<<"\n";
    for (int w=0; w<3; w++) myAssert( E[ F[fa].ei[w] ].has(fa), "FAIL CHECK FA "<<w<<" wa="<<wa<<"\n");
    for (int w=0; w<3; w++) myAssert( E[ F[fb].ei[w] ].has(fb), "FAIL CHECK FB "<<w<<" wb="<<wb<<"\n");

    for (int w=0; w<2; w++) myAssert( F[ E[ei].fi[w] ].has(ei), "FAIL CHECK EI\n");

}

void Mesh::removeDontcare(){
    auto oldV = V;
    std::vector<int> oldToNew(V.size(),-1);
    V.clear();
    for(int i=0; i<(int)oldV.size(); i++) {
        const Vert &v(oldV[i]);
        if (!v.dontcare) {
            oldToNew[i] = V.size();
            V.push_back(v);
        }
    }
    auto oldF = F;
    F.clear();
    for(Face& f:oldF) {
        for (int w=0; w<3; w++) f[w] = oldToNew[f[w]];
        if ((f[0]!=-1)&&(f[1]!=-1)&&(f[2]!=-1)) F.push_back(f);
    }
    updateIndices();
    updateValencies();
}

bool Mesh::checkIfBest(){
    if (nVal<bestVal) {
        std::cout<<"Progress: "<<bestVal<<"-->"<<nVal; //<<" (with "<<forcedTurnsStart<<")\n";
        bestVal = nVal;


        // save best conf
        bestConf = F;
        //if (bestVal==0) break;
        return true;
    }
    else {
        std::cout<<"No progress\n";
        return false;
    }
}

void Mesh::regularizeByFlips(int howDeep, int howWide){

    updateIndices();

    dontCareAboutBoundaries();
    for (int i=0; i<2; i++){
        propagateDontcareToFaces();
        propagateDontcareToVerts();
    }

    //updateAverageVecDir();

    updateIndices();
    updateValencies();
    bestVal = nVal;
    //sanityCheck();

    int forcedTurnsStart = 0;
    int forcedTurns = 0;

    int patience = 0;
    for (int i=0; i<10000000; i++) {
        int ei = bestFlip( forcedTurns > 0);

        forcedTurns--;
        if (ei==-1) {
            // no more profitable moves...
            if (checkIfBest()) {
                std::cout<<" (with "<<forcedTurnsStart<<")\n";
                forcedTurnsStart = 0;
                if (nVal==0) break;
            } else {
                // recover last best conf
                //myAssert(forcedTurns<0,"WHY?");
                F = bestConf; updateIndices(); updateValencies();
                forbiden.clear();
            }
            if (patience == 0) {forcedTurnsStart++; patience=howWide;} else patience--;
            //forcedTurnsStart++;

            forcedTurns = forcedTurnsStart;
            if (forcedTurnsStart>howDeep) break;
            continue;
        }

        applyFlip(ei);

        {        // check

            //sanityCheck();
            //updateValencies();
            //updateIndices();
            //int test = updateValencies();
            //myAssert(test==nVal,"Val mismatch "<<test<<"!="<<nVal<<"\n");
            //nVal = test;
        }
    }

    F = bestConf; updateIndices(); updateValencies();
    updateIndices();
    propagateDontcareToFaces();

}
}
}
