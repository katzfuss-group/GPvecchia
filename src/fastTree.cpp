#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include <iostream>
#include <set>
#include <queue>
#include <utility>
#include <algorithm>

using namespace arma;
using namespace std;




uvec clusterEqual(mat locs, int K, int dimStart){

  int n = locs.n_rows;
  K = pow(2, ceil(log2(K)));

  map<int, uvec> *regions = new map<int, uvec>();
  map<int, uvec> *newRegions = new map<int, uvec>();

  uvec root = linspace<uvec>(0, n-1, n);
  (*regions)[0] = root;

  for( int power=0; power<log2(K); ++power ){

    //cout << "power: " << power << endl;

    for( auto const& x: *regions ){

      int id = x.first;
      //cout << "id: " << id << endl;
      uvec regInds = x.second;
      mat regLocs = locs.rows(regInds);

      int d = (dimStart + power) % locs.n_cols;
      double cutoff = median(regLocs.col(d));

      uvec r1 = regInds( find(regLocs.col(d)>cutoff) );
      uvec r2 = regInds( find(regLocs.col(d)<cutoff) );

      uvec border = regInds( find(regLocs.col(d)==cutoff) );
      int lengthDiff = r1.size() - r2.size();
      if(lengthDiff>0){
        r2 = join_cols(r2, border.head(lengthDiff));
      } else {
        r1 = join_cols(r1, border.head(abs(lengthDiff)));
      }

      if(border.size() - abs(lengthDiff)>0){
        border = border.tail(border.size() - abs(lengthDiff));
	      int halfLength = border.size()/2;
	      r1 = join_cols(r1, border.head(halfLength));
	      r2 = join_cols(r2, border.tail(border.size() - halfLength));
      }
      (*newRegions)[2*id] = r1;
      (*newRegions)[2*id+1] = r2;

    }

    map<int, uvec> *temp = regions;
    regions = newRegions;
    newRegions = temp;

  }

  uvec clusters = zeros<uvec>(n);

  for(auto const &it: *regions){
    int id = it.first;
    uvec inds = it.second;
    clusters(inds).fill(id);
  }

  return clusters;
}



umat neighbours2NNarray(map<unsigned int,uvec> NNs, int m){

  umat NNarray = zeros<umat>(NNs.size(), m);
  for( auto const &x: NNs ){
    int knotNo = x.first;
    uvec knotNeighb = x.second;
    uvec row = zeros<uvec>(m);
    row.head(knotNeighb.size()) = knotNeighb;
    NNarray.row(knotNo-1) = row.t();
  }

  return NNarray;

}



string parent(string id){

  int last_ = id.find_last_of("_");
  return id.substr(0, last_);

}



umat getNNmatrix(map<string,uvec> knots, int m) {

  map<unsigned int,uvec> neighbours;
  set<string> roots;
  int minLength = -1;
  int effM = 0;

  //find the roots of the tree
  for(auto const &x: knots){
    string key = x.first;
    if(key.length()<minLength){
      minLength = key.length();
      roots.clear();
    }
    if(key.length()<=minLength){
      roots.insert(key);
    }
  }

  for(auto const &id: roots){
    if(knots[id].size()==0) {
      continue;
    }

    vector<unsigned int> condSet;
    for(auto const &knot: knots[id]){
      condSet.insert(condSet.begin(), knot);
      uvec NN = zeros<uvec>(m);
      uvec condS = uvec(condSet);
      NN.head(condSet.size()) = condS;
      neighbours[knot] = condS;
      if(condS.size()>effM){
        effM = condS.size();
      }
    }
  }


  for(auto const &x: knots){
    string nodeId = x.first;
    uvec nodeKnots = knots[nodeId];
    if( roots.count(nodeId)!=0 || nodeKnots.size()==0 ) {
      continue;
    }

    uvec condSet;
    uvec parentKnots = knots[parent(nodeId)];
    unsigned int lastKnot = parentKnots(parentKnots.size()-1);
    condSet = neighbours[lastKnot];

    for( auto const &knot: nodeKnots ){
      uvec knotVec = {knot};
      condSet = join_cols(knotVec, condSet);
      neighbours[knot] = uvec(condSet);
      if(condSet.size()>effM){
        effM = condSet.size();
      }
    }

  }

  umat NNarray = neighbours2NNarray( neighbours, effM );
  return NNarray;

}



int res(string id){
  int m = count(id.begin(), id.end(), '_');
  return m;
}


map<string, uvec> knotTree(mat locs, map<string, uvec> mraParams){

  uvec J = mraParams["J"];
  int M = mraParams["M"](0);
  uvec r = mraParams["r"];
  int N = locs.n_rows;

  map<string, uvec> knots;
  queue<pair<string,uvec>> remaining;

  remaining.push( make_pair( "r", linspace<uvec>(0, N-1, N)) );

  while(!remaining.empty()){

    pair<string, uvec> region = remaining.front();
    string id = region.first;

    //cout << "id: " << id << endl;

    int m = res(id);
    uvec regInds = region.second;
    int n = regInds.size();
    uvec clusters;

    if(m<M){

      int rEff = min(r[m], regInds.size());

      if(rEff>0){
	      knots[id] = regInds.head(rEff)+1; // need to shift ind numbering by one b/c NNarray has to completed with zeros
	      regInds = regInds.tail(regInds.size()-rEff);
      }

      mat regLocs = locs.rows(regInds);

      if(regInds.size()==0) {
	      clusters.reset();
      } else {
	      if( J(m)>regLocs.n_rows ){
	        clusters = linspace<uvec>(0, regLocs.n_rows-1, regLocs.n_rows);
	      } else {
  	      int dimStart = m % 2 + 1;
	        clusters = clusterEqual(regLocs, J(m), dimStart);
	      }
      }

      for( int childNo=0; childNo<J(m); ++childNo ){
	      string childId = id + "_" + to_string(childNo);
	      uvec childInds = regInds( find(clusters==childNo) );
	      remaining.push( make_pair( childId, childInds ) );
      }
    } else {
      knots[id] = regInds+1;
    }
    remaining.pop();
  }
  return knots;
}



// [[Rcpp::export]]
umat generateNNarray(mat locs, uvec J, int M, uvec r, int m){

  map<string, uvec> mraParams;
  mraParams["M"] = M;
  mraParams["J"] = J;
  mraParams["r"] = r;

  map<string, uvec> knots = knotTree(locs, mraParams);
  umat NNarray = getNNmatrix(knots, sum(mraParams["r"]));
  return NNarray;

}


/*

int main(int argc, const char **argv) {

  arma_rng::set_seed(1988);

  int N = 9;
  mat locs = randu(N,2);
  uvec J = {8};
  uvec M = {1};
  uvec r = {1, 1};

  map<string, uvec> mraParams;
  mraParams["M"] = M;
  mraParams["J"] = J;
  mraParams["r"] = r;


  //cout << generateNNarray(locs, mraParams) << endl;

    //for(auto const &x: knots){
    //cout << x.first << " : " << x.second.t() << "\n";
    //}

  return 0;

}*/
