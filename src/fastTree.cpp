#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <set>
#include <queue>
#include <utility>
#include <algorithm>
//#include <experimental/optional>

using namespace Rcpp;
using namespace arma;
using namespace std;




arma::uvec clusterEqual(arma::mat locs, int K, int dimStart){

  int n = locs.n_rows;
  K = pow(2, ceil(log2(K)));

  map<int, arma::uvec> *regions = new map<int, arma::uvec>();
  map<int, arma::uvec> *newRegions = new map<int, arma::uvec>();

  arma::uvec root = linspace<arma::uvec>(0, n-1, n);
  (*regions)[0] = root;

  for( int power=0; power<log2(K); ++power ){

    //cout << "power: " << power << endl;

    for( auto const& x: *regions ){

      int id = x.first;
      //cout << "id: " << id << endl;
      arma::uvec regInds = x.second;
      arma::mat regLocs = locs.rows(regInds);

      int d = (dimStart + power) % locs.n_cols;
      double cutoff = median(regLocs.col(d));

      arma::uvec r1 = regInds( find(regLocs.col(d)>cutoff) );
      arma::uvec r2 = regInds( find(regLocs.col(d)<cutoff) );

      arma::uvec border = regInds( find(regLocs.col(d)==cutoff) );
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

    map<int, arma::uvec> *temp = regions;
    regions = newRegions;
    newRegions = temp;

  }

  arma::uvec clusters = zeros<arma::uvec>(n);

  for(auto const &it: *regions){
    int id = it.first;
    arma::uvec inds = it.second;
    clusters(inds).fill(id);
  }

  return clusters;
}



arma::umat neighbours2NNarray(map<unsigned int,arma::uvec> NNs, int m){

  arma::umat NNarray = zeros<arma::umat>(NNs.size(), m);
  for( auto const &x: NNs ){
    int knotNo = x.first;
    arma::uvec knotNeighb = x.second;
    arma::uvec row = zeros<arma::uvec>(m);
    row.head(knotNeighb.size()) = knotNeighb;
    NNarray.row(knotNo-1) = row.t();
  }

  return NNarray;

}



string parent(string id){

  int last_ = id.find_last_of("_");
  return id.substr(0, last_);

}



arma::umat getNNmatrix(map<string,arma::uvec> knots, int m) {

  map<unsigned int,arma::uvec> neighbours;
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
      arma::uvec NN = zeros<arma::uvec>(m);
      arma::uvec condS = arma::uvec(condSet);
      NN.head(condSet.size()) = condS;
      neighbours[knot] = condS;
      if(condS.size()>effM){
        effM = condS.size();
      }
    }
  }


  for(auto const &x: knots){
    string nodeId = x.first;
    arma::uvec nodeKnots = knots[nodeId];
    if( roots.count(nodeId)!=0 || nodeKnots.size()==0 ) {
      continue;
    }

    arma::uvec condSet;
    arma::uvec parentKnots = knots[parent(nodeId)];
    unsigned int lastKnot = parentKnots(parentKnots.size()-1);
    condSet = neighbours[lastKnot];

    for( auto const &knot: nodeKnots ){
      arma::uvec knotVec = {knot};
      condSet = join_cols(knotVec, condSet);
      neighbours[knot] = arma::uvec(condSet);
      if(condSet.size()>effM){
        effM = condSet.size();
      }
    }

  }

  arma::umat NNarray = neighbours2NNarray( neighbours, effM );
  return NNarray;

}



int res(string id){
  int m = count(id.begin(), id.end(), '_');
  return m;
}


tuple<map<string, arma::uvec>, int, arma::uvec, arma::uvec > knotTree(arma::mat locs, map<string, arma::uvec> mraParams){

  arma::uvec J = mraParams["J"];
  int M = mraParams["M"](0);
  arma::uvec r = mraParams["r"];
  int N = locs.n_rows;

  map<string, arma::uvec> knots;
  queue<pair<string,arma::uvec>> remaining;

  remaining.push( make_pair( "r", linspace<arma::uvec>(0, N-1, N)) );

  int eff_M = 0;
  arma::uvec eff_J = zeros<arma::uvec>(M);
  arma::uvec eff_r = zeros<arma::uvec>(M+1);



  while(!remaining.empty()){

    pair<string, arma::uvec> region = remaining.front();
    string id = region.first;

    int m = res(id);

    eff_M = max(m, eff_M);
    arma::uvec regInds = region.second;
    int n = regInds.size();
    arma::uvec clusters;

    if(m<M){

      int rEff = min(r[m], regInds.size());
      if(eff_r(m)==0){
        eff_r(m) = rEff;
      } else if(eff_r(m)!=rEff) {
        eff_r(m) = 1e8;
      }

      if(rEff>0){
	      knots[id] = regInds.head(rEff)+1; // need to shift ind numbering by one b/c NNarray has to completed with zeros
	      regInds = regInds.tail(regInds.size()-rEff);
      }

      arma::mat regLocs = locs.rows(regInds);

      if(regInds.size()==0) {
	      clusters.reset();
      } else {
	      if( J(m)>regLocs.n_rows ){
	        clusters = linspace<arma::uvec>(0, regLocs.n_rows-1, regLocs.n_rows);
	      } else {
  	      int dimStart = m % 2 + 1;
	        clusters = clusterEqual(regLocs, J(m), dimStart);
	      }
      }

      for( int childNo=0; childNo<J(m); ++childNo ){
	      string childId = id + "_" + to_string(childNo);
	      arma::uvec childInds = regInds( find(clusters==childNo) );
	      remaining.push( make_pair( childId, childInds ) );
      }
    } else {
      knots[id] = regInds+1;
    }
    remaining.pop();
  }

  //tuple<map<string, arma::uvec>, int, arma::uvec, arma::uvec> output = {knots, eff_M, J, eff_r};
  //return output;
  
  return std::make_tuple(knots, eff_M, J, eff_r);
  

}



// [[Rcpp::export]]
List generateNNarray(arma::mat locs, arma::uvec J, int M, arma::uvec r, int m){

  map<string, arma::uvec> mraParams;
  mraParams["M"] = M;
  mraParams["J"] = J;
  mraParams["r"] = r;

  //auto [ knots, Meff, Jeff, reff ] = knotTree(locs, mraParams);

  map<string, uvec> knots;
  int Meff;
  uvec Jeff;
  uvec reff;
  std::tie(knots, Meff, Jeff, reff) = knotTree(locs, mraParams);
  
  arma::umat NNarray = getNNmatrix(knots, sum(mraParams["r"]));

  List output;
  output["NNarray"] = NNarray;

  output["Meff"] = Meff;
  output["Jeff"] = Jeff;
  output["reff"] = reff;

  return output;

}


/*

int main(int argc, const char **argv) {

  arma_rng::set_seed(1988);

  int N = 9;
  mat locs = randu(N,2);
  arma::uvec J = {8};
  arma::uvec M = {1};
  arma::uvec r = {1, 1};

  map<string, arma::uvec> mraParams;
  mraParams["M"] = M;
  mraParams["J"] = J;
  mraParams["r"] = r;


  //cout << generateNNarray(locs, mraParams) << endl;

    //for(auto const &x: knots){
    //cout << x.first << " : " << x.second.t() << "\n";
    //}

  return 0;

}*/
