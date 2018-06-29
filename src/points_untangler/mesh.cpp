#include "mesh.h"

#include "grid.h"

#include<set>
#include<map>

#include "my_assert.h"


namespace cellogram
{
namespace PointsUntangler
{

/*
void Mesh::setFaceRegularity(){
    for (Face &f:F) {
        f.regularity =  goodTriangle( f[0],f[1],f[2] )*0.5;
        for (int w=0; w<3; w++)
            f.regularity += V[ f[w] ].distToIrr;
    }
}*/

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

int Mesh::bestFace()const{

    int res = -1;
    scalar best = -1000000;
    for (int i=0; i<(int)F.size(); i++) {
        scalar appeal = V[ F[i][0] ].distToIrr + V[ F[i][1] ].distToIrr + V[ F[i][2] ].distToIrr ;

        if (appeal>best) {
            best = appeal;
            res = i;
        }
    }
    myAssert(res!=-1,"Cannot find face!");
    //std::cout<<"BEST is face "<<res<<" with "<<best<<"("<<goodTriangle(F[res][0],F[res][1],F[res][2])<<"\n";
    return res;
}

vec2 Mesh::parellelogramRule(int fi, int ei) const{
    const Face& f( F[fi] );
    int w = f.cornerOfEdge(ei);
    vec2 p0 = V[F[fi].vi[(w+0)%3]].p;
    vec2 p1 = V[F[fi].vi[(w+1)%3]].p;
    vec2 p2 = V[F[fi].vi[(w+2)%3]].p;
    return p1 + p0 - p2;
}

scalar Mesh::parallelogramError( int ei ) const{
    const Edge& e( E[ei] );
    if (e.fi[1]==-1) return 0;
    if (e.fi[0]==-1) return 0;
    const Face& fi= F[e.fi[0]] ;
    const Face& fj= F[e.fi[1]] ;
    int wi = (fi.cornerOfEdge(ei) + 2)%3;
    int wj = (fj.cornerOfEdge(ei) + 2)%3;
    vec2 pi = V[fi.vi[wi]].p;
    vec2 pa = V[e.vi[0]].p;
    vec2 pb = V[e.vi[1]].p;
    vec2 pj = V[fj.vi[wj]].p;

    return squaredDistance(pi,pa+pb-pj);
}

int Mesh::viOnOtherSide(int fa, int ei) const{
    int wa = F[fa].cornerOfEdge(ei);
    int e = F[fa].ei[wa];
    int fb = E[e].fi[0];
    if (fb == fa) fb = E[e].fi[1];
    if (fb==-1) return -1;
    return F[fb].oppositeVertOfEdge( ei );
;

}

void maybeSwap(int &i, int x, int y){
    if (i==x) i=y;
    else
    if (i==y) i=x;
}

void Mesh::swapVert(int vi, int vj){

    //std::swap( V[vi].p , V[vj].p);

    /*
    for (Face& f : F) {
        maybeSwap(f.vi[0],vi,vj);
        maybeSwap(f.vi[1],vi,vj);
        maybeSwap(f.vi[2],vi,vj);
    }
    for (Edge& e : E) {
        maybeSwap(e.vi[0],vi,vj);
        maybeSwap(e.vi[1],vi,vj);
    }
    std::swap( V[vi].val , V[vj].val);*/

}

scalar Mesh::parallelogramError( int vi, int va, int vb, int vj ) const{

    if (vj==-1) return 9e9;
    vec2 pi = V[vi].p;
    vec2 pa = V[va].p;
    vec2 pb = V[vb].p;
    vec2 pj = V[vj].p;

    return squaredDistance(pi,pa+pb-pj);
}


void Mesh::setDistanceToIrr(){
    for (Vert &v:V) v.distToIrr = ((v.val!=6) || v.dontcare)? 0 : 1000;

    // badly shaped faces count as irregular
    for (const Face &f:F) {
        if (distFromEquilateral( f )>0.3)
        V[f[0]].distToIrr = V[f[1]].distToIrr = V[f[2]].distToIrr =0;
    }
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


void Mesh::fromEigen(const Eigen::MatrixXd &mat){
    for(int i = 0; i < mat.rows(); ++i){
        Vert v; v.p = vec2(mat(i,0),mat(i,1));
        V.push_back(v);
    }
}


void Mesh::fuckUp(){
    //Vert v; v.p = vec2(200,600); V.push_back(v);
    V[2420].p = vec2(200,600);

    std::cout<<"Fucked up mesh :D ! ";
}

void Mesh::sanityCheckValencies(){
    //std::cout<<"testingValencies...";
    std::vector<int> val(V.size(),0);
    for (Edge& e:E) {
        val[e[0]]++;
        val[e[1]]++;
    }

    int nerr = 0;
    for (uint i=0; i<V.size(); i++){
        if (val[i] != V[i].val) {
            nerr++;
            std::cout<<"\nVALENCY ERROR: "<<val[i] <<"!="<< V[i].val<<" \n";
        }
    }
    //if (nerr==0) std::cout<<"OK\n";

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

void Mesh::buildEdgesFromFaces(){

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

void Mesh::propagateDontcareV2F(){

    for (Face& f:F) {
        f.dontcare = false;
        for (int w=0; w<3; w++)
        if (V[f[w]].dontcare) f.dontcare = true;
    }
}

void Mesh::propagateFixedF2E(){
    for (Edge& e:E) e.fixed = false;
    for (Face& f:F)  {
        if (f.fixed) for (int w=0; w<3; w++) {
            E[f.ei[w]].fixed = true;
        }
    }

}

void Mesh::propagateDontcareF2V(){

    for (Vert& v:V) v.dontcare = false;
    for (Face& f:F)  {

        for (int w=0; w<3; w++) {
            if (f.dontcare) V[f[w]].dontcare = true;
        }
    }
}

bool Mesh::canFlip(int ei){

    //if (E[ei].fixed) return false;
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

    //if (V[v0].dontcare || V[v1].dontcare || V[va].dontcare || V[vb].dontcare) return false;

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

// 0.0 if equilateral, 1.0 if flat
scalar Mesh::distFromEquilateral(const Face &f) const{

    scalar a0 = (V[f[0]].p-V[f[1]].p).norm();
    scalar a1 = (V[f[1]].p-V[f[2]].p).norm();
    scalar a2 = (V[f[1]].p-V[f[2]].p).norm();

    // sort to a2>a1>a0
    if (a0>a1) std::swap(a0,a1);
    if (a1>a2) std::swap(a1,a2);
    if (a0>a1) std::swap(a0,a1);

    if (a2==0) return 1;
    return (a2-a0)/a2;
}

scalar myPow(scalar b, scalar exp){
    scalar res = 1;
    int i; for (i=0; i<exp; i++) res*=b;
    res*=(exp-i)*b+(1-exp+i);
    return res;
}

scalar Vert::price()const{
    //return myPow(0.9,disputed);
    //return (disputed>0)?0.5:1.0;
    //return 2.0-disputed;
    return 1;
    //return std::max((scalar)0,0.4f-disputed);
}

FlipScore Mesh::evaluateFlip(int ei, Grid &g){
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
    { if (V[v0].val>6) valReduction+=V[v0].price(); else valReduction-=V[v0].price(); }
    { if (V[v1].val>6) valReduction+=V[v1].price(); else valReduction-=V[v1].price(); }
    { if (V[va].val<6) valReduction+=V[va].price(); else valReduction-=V[va].price(); }
    { if (V[vb].val<6) valReduction+=V[vb].price(); else valReduction-=V[vb].price(); }
    /*if (V[v0].val==7) valReduction++; if (V[v0].val==6) valReduction--;
    if (V[v1].val==7) valReduction++; if (V[v1].val==6) valReduction--;
    if (V[va].val==5) valReduction++; if (V[va].val==6) valReduction--;
    if (V[vb].val==5) valReduction++; if (V[vb].val==6) valReduction--;*/

    /*scalar mult01 = ((V[v0].disputed>0)||(V[v1].disputed>0))? 1.2 : 1;
    scalar multAB = ((V[va].disputed>0)||(V[vb].disputed>0))? 1.2 : 1;*/

    scalar misReduction = 0.0;

    if (g.vert.size()) misReduction = g.misalignmentOptimist(v0,v1) - g.misalignmentOptimist(va,vb);

    //lenReduction += (distance(V[v1].p,V[v0].p) - distance(V[va].p,V[vb].p))/10.0;
    //lenReduction += happynessTri( fa ) + happynessTri( fb )

    /*if (desperate) {
        lenReduction += 10.0*(rand()%1000*0.001);
    }*/
    myAssert(va!=vb,"Dont!");

    return FlipScore(valReduction,misReduction);
}


FlipScore lastMove(0,0);
std::vector<int> forbiden;

static bool contains(const std::vector<int>& v, int i){
    for (int j:v) if (j==i) return true; return false;
}

int Mesh::bestFlip(FlipScore &bestScore, Grid &g){
    int winner = -1;

    bestScore = FlipScore(-999,0);

    //else bestScore.second = -20.0;
    for (int i=0; i<(int)E.size(); i++) {
        if (!canFlip(i)) continue;

        auto score = evaluateFlip(i, g);

        if (contains(forbiden,i)) continue;

        if (bestScore<score) {
            bestScore = score;
            winner = i;
        }
    }
    //forbiden.push_back(winner);
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
    int v0 = E[ei][0];
    int v1 = E[ei][1];
    //std::cout<<"FLIPPP "<<v0<<","<<v1<<"\n";
    int fa = E[ei].fi[0];
    int fb = E[ei].fi[1];
    int wa = -1,wb =-1;

    for (int w=0; w<3; w++) if ((F[fa][w]!=v0) && (F[fa][w]!=v1)) wa = w;
    for (int w=0; w<3; w++) if ((F[fb][w]!=v0) && (F[fb][w]!=v1)) wb = w;

    //std::cout<<"ERA: "<<fa<<" ->"<<F[fa].ei[(wa+2)%3]             <<
    //              "  F = "<<E[ F[fa].ei[(wa+2)%3]].fi[0]<<","<<E[ F[fa].ei[(wa+2)%3]].fi[1]<<"\n";

    myAssert(wa!=-1,"NO WA\n");
    myAssert(wb!=-1,"NO WA\n");

    myAssert((F[fa][(wa+1)%3]==v0) || (F[fa][(wa+1)%3]==v1),"NO WA 1\n");
    myAssert((F[fa][(wa+2)%3]==v0) || (F[fa][(wa+2)%3]==v1),"NO WA 2\n");
    myAssert((F[fb][(wb+1)%3]==v0) || (F[fb][(wb+1)%3]==v1),"NO WB 1\n");
    myAssert((F[fb][(wb+2)%3]==v0) || (F[fb][(wb+2)%3]==v1),"NO WB 2\n");

    myAssert(F[fa].ei[(wa+1)%3]==ei,"NO WA3");
    myAssert(F[fb].ei[(wb+1)%3]==ei,"NO WB3");
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

    for (int w=0; w<3; w++) myAssert( E[ F[fa].ei[w] ].has(fa), "FAIL CHECK FA "<<w<<" wa="<<wa<<"\n");
    for (int w=0; w<3; w++) myAssert( E[ F[fb].ei[w] ].has(fb), "FAIL CHECK FB "<<w<<" wb="<<wb<<"\n");

    for (int w=0; w<2; w++) myAssert( F[ E[ei].fi[w] ].has(ei), "FAIL CHECK EI\n");

}

/*
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
    buildEdgesFromFaces();
    updateValencies();
}*/

void Mesh::flipAs(const Grid& g){
    std::cout<<"Flipping to match grid...\n";
    int tot = 0;
    while(1){
        int done = 0;
        for (uint ei=0; ei<E.size(); ei++){
            if (!canFlip(ei)) continue;
            int fa = E[ei].fi[0];
            int fb = E[ei].fi[1];
            int va=-1,vb=-1;
            int v0 = E[ei][0];
            int v1 = E[ei][1];

            for (int w=0; w<3; w++) if ((F[fa][w]!=v0) && (F[fa][w]!=v1)) va = F[fa][w];
            for (int w=0; w<3; w++) if ((F[fb][w]!=v0) && (F[fb][w]!=v1)) vb = F[fb][w];

            //int d01 = g.hopDistanceV(v0,v1);
            //int dab = g.hopDistanceV(va,vb);
            //bool improveSimil = (( d01 > dab ) && ( d01 != -1 ) && ( dab != -1 ) );
            //bool equalSimil = (( d01 == dab ) || ( d01 == -1 ) || ( dab == -1 )  );

            /*int valReduction = 0;
            { if (V[v0].val>6) valReduction++; else valReduction--; }
            { if (V[v1].val>6) valReduction++; else valReduction--; }
            { if (V[va].val<6) valReduction++; else valReduction--; }
            { if (V[vb].val<6) valReduction++; else valReduction--; }
*/
            if ( (!g.areAdjacient(v0,v1)) && (g.areAdjacient(va,vb)) )

            //if  ( dab == 1 && d01!= 1 ) {
            //if  ( dab != -1 && d01!= -1 )
            //if  ( (V[v0].val>=4) || (V[v1].val>=4) )
            /*if ( ( g.posInGrid[va]!=-1 && g.posInGrid[vb]!=-1 &&
                   g.posInGrid[v0]!=-1 && g.posInGrid[v1]!=-1)
                //((V[va].val<4) || (V[vb].val<4))
                 //||
                 && ( //( d01 > dab )
                      //improveSimil
                      g.areAdjacient(va,vb)
                      &&
                      !g.areAdjacient(v0,v1)
                 )
                  )*/
                {
            //if (improveSimil ) { //|| (equalSimil  && (valReduction>0) )) {
                applyFlip(ei);
                done++;
            }
        }
        if (!done) break;
        tot+=done;
    }
    std::cout<<tot<<" flips done!\n";
}

void Mesh::greedyFlips(int howDeep, Grid &g){

    std::cout<<"Greedy mesh regularization by FLIPS...\n";
    //buildEdgesFromFaces();

    //dontCareAboutBoundaries();
    //for (int i=0; i<2; i++){
    //        propagateDontcareV2F();
    //    propagateDontcareF2V();
    //}


    //buildEdgesFromFaces();
    updateValencies();
    //sanityCheck();

    FlipScore totScore;

    bestEver.score=totScore;
    bestEver.F=F;
    bestEver.E=E;

    int patience = 0;
    int ndone=0;
    while (1) {
        FlipScore score;
        int ei = bestFlip( score, g  );

        assert(ei>=0);

        if (score.isPos()) {
            applyFlip(ei);
            totScore += score;
            if (ndone++ > 1000) {
                std::cout<<"ERROR infinte loop in greedy swaps! "<<ei<<" sci=ore move: " << score.valReduction << "," << score.lenReduction << "\n";
                if (ndone>1100) break;
            }
        } else if (howDeep>0) {
            // do a neg move?

            if(bestEver.score < totScore ) {
                bestEver.score=totScore;
                bestEver.F=F;
                bestEver.E=E;
                patience = 0;
                forbiden.clear();
            }

            patience++;

            if (patience>howDeep) {
                // rollback and bail out
                totScore = bestEver.score;
                F = bestEver.F;
                E = bestEver.E;
                updateValencies();
                break;
            }

            applyFlip(ei);
            forbiden.push_back(ei);
            totScore += score;
        }
        //forcedTurns--;
        if (ei==-1) { break; }


        //sanityCheckValencies();
    }

    std::cout<<ndone << " flips done!\n";



}


}} // namespaces
