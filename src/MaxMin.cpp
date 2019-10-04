#include <Rcpp.h>
#include <cmath>
#include <cstring>
#ifdef _WIN32
  #include <malloc.h>
#endif
#include <cstdlib>
using namespace std;
using namespace Rcpp;
 

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// From points.h & points.c

/*A struct containing all the information relevant to the data point*/
typedef struct point {
  int Id;
  int level;
  int d;
  double* x;
  int NChildren;
  int firstChildKey;
  struct point* parent;
  /*Heap Stuff:*/
  double distSQ;
  int hlevel;
  struct point** hnode;
  int hId;
  int leaf;
} point;

/*compute distance between two points*/
//double dist2( const point* x, const point* y ){
//	int i;
//	double result = 0;
//	for (i = 0; i < x->d; ++i) {
//		result += pow(x->x[i] - y->x[i], 2.);
//	}
//	return result;
//}

/*creates a set of n^2 points forming 2-dimensional grid*/
double* create_coords2d(const int n) {
  double h = 1.0 / ((double)(n + 1));
  int k, l;
  int N = pow(n, 2);
  /*allocates a N times 2 array*/
  double* coords = (double*)malloc(sizeof(double) * 2 * N);
  //if (coords == NULL) exit(1);
  /*Fill the allocated array with the right coordinates*/
  for (k = 0; k < n; ++k) {
    for (l = 0; l < n; ++l) {
      coords[2 * (k*n + l) + 0] = h * (double)(l + 1);
      coords[2 * (k*n + l) + 1] = h * (double)(k + 1);
    }
  }
  return coords;
}

/*creates a set of N point uniformly distributed in [0,1]^d
double* create_coordsUnif(const unsigned int d, const unsigned int N) {
  
  double* coords = (double*)malloc(sizeof(double) * d * N);
  if (coords == NULL) exit(1);
  srand(1);
  for (unsigned int k = 0; k < d * N; ++k) {
    coords[k] = ((double)rand()) / ((double)RAND_MAX);
  }
  return coords;
}*/

/*creates point structs from an input of coordinates*/
point* create_Points(double* const x, const int d, const int N) {
  //int k, l;
  int k;
  point* points = (point*)malloc(sizeof(point) * N);
  //if (points == NULL) exit(1);
  for (k = 0; k < N; ++k) {
    points[k].Id = k;
    points[k].d = d;
    points[k].x = &x[d * k];
    points[k].parent = NULL;
    points[k].firstChildKey = 0;
    points[k].NChildren = 0;
    points[k].distSQ = 10000.;
  }
  return points;
}

/*destructor for points*/
void destruct_coords(double* coords) {
  free(coords);
}


/*A max heap to keep track of the next point to be updated*/
typedef struct heap {
  point** elements;
  int N;
} heap;

/*Function comparing squared distances*/
int compareSQ(const void* p1, const void*p2) {
  double diff = ((double)(*((point**)p1))->distSQ - (double)(*((point**)p2))->distSQ);
  if (diff > 0) return -1;
  else return 1;
}

/*Function comparing levels*/
int compareLevel(const void* p1, const void*p2) {
  return ((point*)p1)->level - ((point*)p1)->level;
}

/*Function that initialises a heap*/
void heap_init(heap* h, point* points, const int N) {
  h->elements = (point**)malloc(sizeof(point*) * N);
  h->N = N;
  //if (h->elements == NULL) exit(1);
  int k;
  for (k = 0; k < N; ++k) {
    h->elements[k] = &points[k];
  }
  qsort(h->elements, N, sizeof(point*), compareSQ);
  for (k = 0; k < N; ++k) {
    h->elements[k]->hnode = &(h->elements[k]);
  }
}

/*Function that performs a single step of heap reordering. The output is -1, if no
* swapping has taken place and the new heap id, otherwise.*/
//int heap_sortStep(heap* h, const int k);

/*Function that sorts the heap.
* TODO: Implement using heap structure */
void heap_sort(heap* h) {
  qsort(h->elements, h->N, sizeof(point*), compareSQ);
  int k;
  for (k = 0; k < h->N; ++k) {
    h->elements[k]->hnode = &(h->elements[k]);
    h->elements[k]->hId = k;
  }
}

void heap_destruct(heap* h) {
  free(h->elements);
}

/*A struct that provides storage for the children of the individual nodes*/
typedef struct daycare {
  int size;
  int sizeBuffer;
  point** data;
} daycare;

void daycare_init(daycare* dc, const int N) {
  // int i;
  dc->sizeBuffer = 2 * N;
  dc->size = 0;
  dc->data = (point**)malloc(dc->sizeBuffer * sizeof(point*));
  //if (dc->data == NULL) exit(1);
}

void daycare_add(daycare* dc, point* addendum) {
  if (dc->size == dc->sizeBuffer) {
    dc->sizeBuffer *= 2;
    dc->data = (point**)realloc(dc->data, dc->sizeBuffer * sizeof(point*));
    //if (dc->data == NULL) exit(1);
  }
  ++dc->size;
  dc->data[dc->size - 1] = addendum;
}

void daycare_destruct(daycare* dc) {
  free(dc->data);
}

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// From sortSparse.h & sortSparse.c

typedef struct heapNode {
  /*distance squared to closest point that is already taken out, negative for points that are taken out*/
  double dist;
  struct heapNode** handleHandle;
  struct heapNode* leftChild;
  struct heapNode* rightChild;
  /*might not be needed:*/
  unsigned int Id;
} heapNode;

heapNode* _moveDown(heapNode* const a) {
  /*If the nodes has no children:*/
  if (a->leftChild == NULL) {
    return NULL;
  }
  /*If the node has only one child:*/
  if (a->rightChild == NULL) {
    if (a->dist < a->leftChild->dist) {
      /*swaps a with its left child*/
      const double tempDist = a->dist;
      a->dist = a->leftChild->dist;
      a->leftChild->dist = tempDist;
      *(a->handleHandle) = a->leftChild;
      *(a->leftChild->handleHandle) = a;
      heapNode** const tempHandleHandle = a->handleHandle;
      a->handleHandle = a->leftChild->handleHandle;
      a->leftChild->handleHandle = tempHandleHandle;
      
      const int tempId = a->Id;
      a->Id = a->leftChild->Id;
      a->leftChild->Id = tempId;
      
      return a->leftChild;
    }
    else return NULL;
  }
  /*If the node has two children:*/
  if (a->leftChild->dist > a->rightChild->dist) {
    if (a->dist < a->leftChild->dist) {
      /*swaps a with its left child*/
      const double tempDist = a->dist;
      a->dist = a->leftChild->dist;
      a->leftChild->dist = tempDist;
      *(a->handleHandle) = a->leftChild;
      *(a->leftChild->handleHandle) = a;
      heapNode** const tempHandleHandle = a->handleHandle;
      a->handleHandle = a->leftChild->handleHandle;
      a->leftChild->handleHandle = tempHandleHandle;
      
      const int tempId = a->Id;
      a->Id = a->leftChild->Id;
      a->leftChild->Id = tempId;
      
      return a->leftChild;
    }
    else return NULL;
  }
  else {
    if (a->dist < a->rightChild->dist) {
      /*swaps a with its right child*/
      const double tempDist = a->dist;
      a->dist = a->rightChild->dist;
      a->rightChild->dist = tempDist;
      *(a->handleHandle) = a->rightChild;
      *(a->rightChild->handleHandle) = a;
      heapNode** const tempHandleHandle = a->handleHandle;
      a->handleHandle = a->rightChild->handleHandle;
      a->rightChild->handleHandle = tempHandleHandle;
      
      const int tempId = a->Id;
      //printf( "newId = %d ", a->rightChild->Id);
      a->Id = a->rightChild->Id;
      a->rightChild->Id = tempId;
      
      
      return a->rightChild;
    }
    else return NULL;
  }
}

/*Only works as expected, if the new distance is smaller than the old one*/
void update(heapNode* target, const double newDist) {
  target->dist = newDist;
  while (target != NULL) {
    target = _moveDown(target);
  }
}

/*Initialises array of nodes with proper children and writes it into nodes*/
void heapInit(const unsigned int N, heapNode* const nodes, heapNode** const handles) {
  /*heapNode* ret = (heapNode*) malloc( N * sizeof( heapNode ) );
  if( ret == NULL ) exit(1);*/
  for (unsigned int k = 0; k < N; ++k) {
    if (2 * k + 1 >= N) {
      heapNode newnode = { INFINITY, &handles[k], NULL, NULL };
      memcpy(&(nodes[k]), &newnode, sizeof(heapNode));
    }
    else if (2 * k + 2 >= N) {
      heapNode newnode = { INFINITY, &handles[k], &nodes[2 * k + 1], NULL };
      memcpy(&(nodes[k]), &newnode, sizeof(heapNode));
    }
    else {
      heapNode newnode = { 10000., &handles[k], &nodes[2 * k + 1], &nodes[2 * k + 2] };
      memcpy(&(nodes[k]), &newnode, sizeof(heapNode));
    }
    handles[k] = &nodes[k];
    nodes[k].Id = k;
  }
}

typedef struct ijlookup {
  unsigned int pres_i;
  unsigned int N;
  unsigned int S;
  unsigned int S_Buffer;
  unsigned int* i;
  unsigned int* j;
} ijlookup;

void ijlookup_init(ijlookup* lookup, unsigned int N) {
  lookup->pres_i = 0;
  lookup->N = N;
  lookup->S = 0;
  lookup->S_Buffer = N;
  lookup->i = (unsigned int*)malloc((N + 1) * sizeof(unsigned int));
  //if (lookup->i == NULL) exit(1);
  lookup->j = (unsigned int*)malloc(N * sizeof(unsigned int));
  //if (lookup->j == NULL) exit(1);
  lookup->i[0] = 0;
  lookup->i[1] = 0;
}

void ijlookup_newparent(ijlookup* const lookup) {
  ++lookup->pres_i;
  lookup->i[lookup->pres_i + 1] = lookup->i[lookup->pres_i];
}

void ijlookup_newson(ijlookup* const lookup, const unsigned int Id) {
  ++lookup->S;
  if (lookup->S > lookup->S_Buffer) {
    lookup->S_Buffer *= 2;
    lookup->j = (unsigned int*)realloc(lookup->j, lookup->S_Buffer * sizeof(unsigned int));
    //if (lookup->j == NULL) exit(1);
  }
  lookup->j[lookup->S - 1] = Id;
  ++lookup->i[lookup->pres_i + 1];
}

void ijlookup_destruct(ijlookup* lookup) {
  free(lookup->i);
  free(lookup->j);
}

double dist(const unsigned int i, const unsigned int j, const double* const coords, const unsigned int d) {
  double ret = 0;
  for (int k = 0; k < (int)d; ++k) {
    ret += (coords[d * i + k] - coords[d * j + k]) * (coords[d * i + k] - coords[d * j + k]);
  }
  return sqrt(ret);
}

double dist_2d(const unsigned int i, const unsigned int j, const double* const coords) {
  return sqrt((coords[2 * i] - coords[2 * j]) * (coords[2 * i] - coords[2 * j])
                + (coords[2 * i + 1] - coords[2 * j + 1]) * (coords[2 * i + 1] - coords[2 * j + 1]));
}

double dist_3d(const unsigned int i, const unsigned int j, const double* const coords) {
  return sqrt((coords[3 * i] - coords[3 * j]) * (coords[3 * i] - coords[3 * j])
                + (coords[3 * i + 1] - coords[3 * j + 1]) * (coords[3 * i + 1] - coords[3 * j + 1])
                + (coords[3 * i + 2] - coords[3 * j + 2]) * (coords[3 * i + 2] - coords[3 * j + 2]));
}

double dist2(const unsigned int i, const unsigned int j, const double* const coords, const unsigned int d) {
  double ret = 0;
  for (int k = 0; k < (int)d; ++k) {
    ret += (coords[d * i + k] - coords[d * j + k]) * (coords[d * i + k] - coords[d * j + k]);
  }
  return ret;
}

double dist2_2d(const unsigned int i, const unsigned int j, const double* const coords) {
  return (coords[2 * i] - coords[2 * j]) * (coords[2 * i] - coords[2 * j])
  + (coords[2 * i + 1] - coords[2 * j + 1]) * (coords[2 * i + 1] - coords[2 * j + 1]);
}

double dist2_3d(const unsigned int i, const unsigned int j, const double* const coords) {
  return (coords[3 * i] - coords[3 * j]) * (coords[3 * i] - coords[3 * j])
  + (coords[3 * i + 1] - coords[3 * j + 1]) * (coords[3 * i + 1] - coords[3 * j + 1])
  + (coords[3 * i + 2] - coords[3 * j + 2]) * (coords[3 * i + 2] - coords[3 * j + 2]);
}

//void dist2Blocked(double* results, const unsigned int i, const unsigned int* const j0, const unsigned int n, const double* const coords, const unsigned int d);
//void dist2BlockedVec(double* results, const double* const ivec, const unsigned int* const j0, const unsigned int n, const double* const coords, const unsigned int d);
double dist2fix(const double* x, const unsigned int j, const double* const coords, const unsigned int d);

inline double in_dist2(const unsigned int i, const unsigned int j, const double* const coords, const unsigned int d) {
  double ret = 0;
  for (unsigned int k = 0; k < d; ++k) {
    ret += (coords[d * i + k] - coords[d * j + k]) * (coords[d * i + k] - coords[d * j + k]);
  }
  return ret;
};

void  determineChildren(heapNode* const nodes, heapNode** const handles, ijlookup* const lookup, unsigned int* const parents, const double* const coords, const unsigned int d, const unsigned int N, const unsigned int Id, const unsigned int iter) {
  const double pivotDist = nodes[0].dist;
  /*Need to save the kmin nad kmax beforehand, since otherwise the node might become its own parent, later on */
  const int kmin = lookup->i[parents[Id]];
  const int kmax = lookup->i[parents[Id] + 1];
  ijlookup_newparent(lookup);
  
  
  for (unsigned int k = kmin; (int)k < kmax; ++k) {
    //if( lookup->j[ k ] >= iter ){
    const double tempDist2 = dist2(Id, lookup->j[k], coords, d);
    if (tempDist2 < pivotDist * pivotDist) {
      double jDist = handles[lookup->j[k]]->dist;
      if (tempDist2 < jDist*jDist) {
        update(handles[lookup->j[k]], sqrt(tempDist2));
        jDist = sqrt(tempDist2);
      }
      ijlookup_newson(lookup, lookup->j[k]);
      if (sqrt(tempDist2) + jDist < pivotDist) {
        parents[lookup->j[k]] = iter;
      }
    }
    //}
  }
}

void determineChildren_2d(heapNode* const nodes, heapNode** const handles, ijlookup* const lookup, unsigned int* const parents, const double* const coords, const unsigned int N, const unsigned int Id, const unsigned int iter) {
  const double pivotDist = nodes[0].dist;
  /*Need to save the kmin nad kmax beforehand, since otherwise the node might become its own parent, later on */
  const int kmin = lookup->i[parents[Id]];
  const int kmax = lookup->i[parents[Id] + 1];
  ijlookup_newparent(lookup);
  
  
  for (unsigned int k = kmin; (int)k < kmax; ++k) {
    //if( lookup->j[ k ] >= iter ){
    const double tempDist2 = dist2_2d(Id, lookup->j[k], coords);
    if (tempDist2 < pivotDist * pivotDist) {
      double jDist = handles[lookup->j[k]]->dist;
      if (tempDist2 < jDist*jDist) {
        update(handles[lookup->j[k]], sqrt(tempDist2));
        jDist = sqrt(tempDist2);
      }
      ijlookup_newson(lookup, lookup->j[k]);
      if (sqrt(tempDist2) + jDist < pivotDist) {
        parents[lookup->j[k]] = iter;
      }
    }
    //}
  }
}

void determineChildren_3d(heapNode* const nodes, heapNode** const handles, ijlookup* const lookup, unsigned int* const parents, const double* const coords, const unsigned int N, const unsigned int Id, const unsigned int iter) {
  const double pivotDist = nodes[0].dist;
  /*Need to save the kmin nad kmax beforehand, since otherwise the node might become its own parent, later on */
  const int kmin = lookup->i[parents[Id]];
  const int kmax = lookup->i[parents[Id] + 1];
  ijlookup_newparent(lookup);
  
  
  for (unsigned int k = kmin; (int)k < kmax; ++k) {
    //if( lookup->j[ k ] >= iter ){
    const double tempDist2 = dist2_3d(Id, lookup->j[k], coords);
    if (tempDist2 < pivotDist * pivotDist) {
      double jDist = handles[lookup->j[k]]->dist;
      if (tempDist2 < jDist*jDist) {
        update(handles[lookup->j[k]], sqrt(tempDist2));
        jDist = sqrt(tempDist2);
      }
      ijlookup_newson(lookup, lookup->j[k]);
      if (sqrt(tempDist2) + jDist < pivotDist) {
        parents[lookup->j[k]] = iter;
      }
    }
    //}
  }
}


//void  determineChildrenBlocked(heapNode* const nodes, heapNode** const handles, ijlookup* const lookup, unsigned int* const parents, const double* const coords, const unsigned int d, const unsigned int N, const unsigned int Id, const unsigned int iter);

void create_ordering(unsigned int* P, unsigned int* revP, double* distances, const unsigned int d, const unsigned int N, const double* coords, unsigned int first_node) {
  /*Function to construct the ordering.
  * Inputs:
  *  P:
  *    An N element array containing the hierarchical ordering
  *  revP:
  *    An N element array containing the inverse of the hierarchical ordering
  *  distances:
  *    An N element array containing the distance ( length scale ) of each dof
  *  d:
  *    The number of spatial dimensions
  *  N:
  *    The number of points
  *  coords:
  *    An d*N element array that contains the different points coordinates, with the
  *    coordinates of a given point in contiguous memory
  */
  
  /*Allocate the heap structure:*/
  heapNode* nodes = (heapNode*)malloc(N * sizeof(heapNode));
  //if (nodes == NULL) exit(1);
  heapNode** handles = (heapNode**)malloc(N * sizeof(heapNode*));
  //if (handles == NULL) exit(1);
  /*Initiate the heap structure*/
  heapInit(N, nodes, handles);
  /*initialising lookup*/
  ijlookup lookup;
  ijlookup_init(&lookup, N);
  /*allocate array to store the parents of dof:
  *the i-th entry of parents will contain the number of its parent in the ordering */
  unsigned int* parents = (unsigned int*)malloc(N * sizeof(unsigned int));
  //if (parents == NULL) exit(1);
  
  /* Add the first parent node: */
  /*TODO Make random?*/
  unsigned int rootId = first_node;
  distances[0] = 0.;
  for (unsigned int k = 0; k < N; ++k) {
    ijlookup_newson(&lookup, k);
    if (dist(rootId, k, coords, d) > distances[0]) {
      distances[0] = dist(rootId, k, coords, d);
    }
    update(handles[k], dist(rootId, k, coords, d));
    parents[k] = 0;
  }
  
  for (unsigned int k = 1; k < N; ++k) {
    unsigned int pivotId = nodes[0].handleHandle - handles;
    distances[k] = nodes[0].dist;
    P[k] = pivotId;
    revP[pivotId] = k;
    determineChildren(nodes, handles, &lookup, parents, coords, d, N, pivotId, k);
  }
  
  ijlookup_destruct(&lookup);
  free(parents);
  free(handles);
  free(nodes);
}

void create_ordering_2d(unsigned int* P, unsigned int* revP, double* distances, const unsigned int N, const double* coords, unsigned int first_node) {
  /*Function to construct the ordering.
  * Inputs:
  *  P:
  *    An N element array containing the hierarchical ordering
  *  revP:
  *    An N element array containing the inverse of the hierarchical ordering
  *  distances:
  *    An N element array containing the distance ( length scale ) of each dof
  *  N:
  *    The number of points
  *  coords:
  *    An d*N element array that contains the different points coordinates, with the
  *    coordinates of a given point in contiguous memory
  */
  
  /*Allocate the heap structure:*/
  heapNode* nodes = (heapNode*)malloc(N * sizeof(heapNode));
  //if (nodes == NULL) exit(1);
  heapNode** handles = (heapNode**)malloc(N * sizeof(heapNode*));
  //if (handles == NULL) exit(1);
  /*Initiate the heap structure*/
  heapInit(N, nodes, handles);
  /*initialising lookup*/
  ijlookup lookup;
  ijlookup_init(&lookup, N);
  /*allocate array to store the parents of dof:
  *the i-th entry of parents will contain the number of its parent in the ordering */
  unsigned int* parents = (unsigned int*)malloc(N * sizeof(unsigned int));
  //if (parents == NULL) exit(1);
  
  /* Add the first parent node: */
  /*TODO Make random?*/
  unsigned int rootId = first_node;
  distances[0] = 0.;
  for (unsigned int k = 0; k < N; ++k) {
    ijlookup_newson(&lookup, k);
    if (dist_2d(rootId, k, coords) > distances[0]) {
      distances[0] = dist_2d(rootId, k, coords);
    }
    update(handles[k], dist_2d(rootId, k, coords));
    parents[k] = 0;
  }
  
  for (unsigned int k = 1; k < N; ++k) {
    unsigned int pivotId = nodes[0].handleHandle - handles;
    distances[k] = nodes[0].dist;
    P[k] = pivotId;
    revP[pivotId] = k;
    determineChildren_2d(nodes, handles, &lookup, parents, coords, N, pivotId, k);
  }
  
  ijlookup_destruct(&lookup);
  free(parents);
  free(handles);
  free(nodes);
}


void create_ordering_3d(unsigned int* P, unsigned int* revP, double* distances, const unsigned int N, const double* coords, unsigned int first_node) {
  /*Function to construct the ordering.
  * Inputs:
  *  P:
  *    An N element array containing the hierarchical ordering
  *  revP:
  *    An N element array containing the inverse of the hierarchical ordering
  *  distances:
  *    An N element array containing the distance ( length scale ) of each dof
  *  N:
  *    The number of points
  *  coords:
  *    An d*N element array that contains the different points coordinates, with the
  *    coordinates of a given point in contiguous memory
  */
  
  /*Allocate the heap structure:*/
  heapNode* nodes = (heapNode*)malloc(N * sizeof(heapNode));
  //if (nodes == NULL) exit(1);
  heapNode** handles = (heapNode**)malloc(N * sizeof(heapNode*));
  //if (handles == NULL) exit(1);
  /*Initiate the heap structure*/
  heapInit(N, nodes, handles);
  /*initialising lookup*/
  ijlookup lookup;
  ijlookup_init(&lookup, N);
  /*allocate array to store the parents of dof:
  *the i-th entry of parents will contain the number of its parent in the ordering */
  unsigned int* parents = (unsigned int*)malloc(N * sizeof(unsigned int));
  //if (parents == NULL) exit(1);
  
  /* Add the first parent node: */
  /*TODO Make random?*/
  unsigned int rootId = first_node;
  distances[0] = 0.;
  for (unsigned int k = 0; k < N; ++k) {
    ijlookup_newson(&lookup, k);
    if (dist_3d(rootId, k, coords) > distances[0]) {
      distances[0] = dist_3d(rootId, k, coords);
    }
    update(handles[k], dist_3d(rootId, k, coords));
    parents[k] = 0;
  }
  
  for (unsigned int k = 1; k < N; ++k) {
    unsigned int pivotId = nodes[0].handleHandle - handles;
    distances[k] = nodes[0].dist;
    P[k] = pivotId;
    revP[pivotId] = k;
    determineChildren_3d(nodes, handles, &lookup, parents, coords, N, pivotId, k);
  }
  
  ijlookup_destruct(&lookup);
  free(parents);
  free(handles);
  free(nodes);
}

//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

// New function: MaxMincpp (integrate all c functions)
// Input: a location matrix (nrow = sample size; ncol = dimension)
// Output: the maxmin order. (Indices start from 1.)

//[[Rcpp::export]]
IntegerVector MaxMincpp(NumericMatrix locations)
{
  unsigned int N = locations.nrow();
  int dim = locations.ncol();
  IntegerVector res(N);
  unsigned int* P = (unsigned int*)malloc(N * sizeof(unsigned int));
  if (P == NULL) return res;
  unsigned int* revP = (unsigned int*)malloc(N * sizeof(unsigned int));
  if (revP == NULL) return res;
  double* distances = (double*)malloc(N * sizeof(double));
  if (distances == NULL) return res;
  double *coords = (double*)malloc(dim * N * sizeof(double));
  
  // Find the average point.
  unsigned int first_node;
  double *average_arr, cur_dist2, min_dist2;
  average_arr = new double[dim];
  for(int j = 0; j < dim; j++)
  {
    average_arr[j] = 0.0;
  }
  for (int i = 0; i < (int)N; i++)
  {
    for(int j = 0; j < dim; j++)
    {
      average_arr[j] += (coords[dim * i + j] = locations(i, j));
    }
  }
  for(int j = 0; j < dim; j++)
  {
    average_arr[j] /= N;
  }
  min_dist2 = -1;
  first_node = -1;
  for (int i = 0; i < (int)N; i++)
  {
    cur_dist2 = 0;
    for(int j = 0; j < dim; j++)
    {
      cur_dist2 += (coords[dim * i + j] - average_arr[j]) * (coords[dim * i + j] - average_arr[j]);
    }
    if (min_dist2 < 0 || cur_dist2 < min_dist2)
    {
      min_dist2 = cur_dist2;
      first_node = i;
    }
  }
  
  delete[] average_arr;
  if (dim == 2)
  {
    create_ordering_2d(P, revP, distances, N, coords, first_node);
  }
  else if (dim == 3)
  {
    create_ordering_3d(P, revP, distances, N, coords, first_node);
  }
  else if (dim >= 1)
  {
    create_ordering(P, revP, distances, dim, N, coords, first_node);
  }
  else
  {
    free(P);
    free(revP);
    free(distances);
    destruct_coords(coords);
    return res;
  }
  res[0] = first_node + 1;
  for (int i = 1; i < (int)N; i++)
    res[i] = P[i] + 1;
  free(P);
  free(revP);
  free(distances);
  destruct_coords(coords);
  return res;
}

