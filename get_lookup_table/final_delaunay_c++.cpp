/* Copyright 2017, Gurobi Optimization, Inc. */

/* This example formulates and solves the following simple MIP model:

     maximize    x +   y + 2 z
     subject to  x + 2 y + 3 z <= 4
                 x +   y       >= 1
     x, y, z binary
*/

#include "gurobi_c++.h"
using namespace std;
#include <vector>
#include <limits>
#include <string>
#include <fstream>
#include <sstream>

#include <math.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <list>
#include <set>
#include <queue>
#include <algorithm>
#include <iostream>

#include <time.h>

#include <fstream>
#include <iomanip>
#include <string>

#include <sstream>

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)
#if !defined(MAX)
#define MAX(A, B)   ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B)   ((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654

#define INF 0x3f3f3f3f

//the length of each link in the arm (should be the same as the one used in runtest.m)
#define LINKLENGTH_CELLS 10

double* GetMap(const std::string &filename)
{
  std::ifstream infile;
  infile.open(filename.c_str());
  string STRING;

  vector<vector<double>> Vmap;
  while(!infile.eof())
  {
    getline(infile,STRING);
    std::istringstream iss(STRING);

    vector<double> tmp;
    for(std::string s; iss >> s;)
    {
      tmp.push_back(atof(s.c_str()));
    }

    Vmap.push_back(tmp);
  }

  cout<<"Size of the map would be: "<<Vmap.size()<<"x"<<Vmap[0].size()<<endl;

  return &Vmap[0][0];

}

typedef struct {
  int X1, Y1;
  int X2, Y2;
  int Increment;
  int UsingYIndex;
  int DeltaX, DeltaY;
  int DTerm;
  int IncrE, IncrNE;
  int XIndex, YIndex;
  int Flipped;
} bresenham_param_t;


void ContXY2Cell(double x, double y, short unsigned int* pX, short unsigned int *pY, int x_size, int y_size)
{
    double cellsize = 1.0;
  //take the nearest cell
  *pX = (int)(x/(double)(cellsize));
  if( x < 0) *pX = 0;
  if( *pX >= x_size) *pX = x_size-1;

  *pY = (int)(y/(double)(cellsize));
  if( y < 0) *pY = 0;
  if( *pY >= y_size) *pY = y_size-1;
}


void get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t *params)
{
  params->UsingYIndex = 0;

  if (fabs((double)(p2y-p1y)/(double)(p2x-p1x)) > 1)
    (params->UsingYIndex)++;

  if (params->UsingYIndex)
    {
      params->Y1=p1x;
      params->X1=p1y;
      params->Y2=p2x;
      params->X2=p2y;
    }
  else
    {
      params->X1=p1x;
      params->Y1=p1y;
      params->X2=p2x;
      params->Y2=p2y;
    }

   if ((p2x - p1x) * (p2y - p1y) < 0)
    {
      params->Flipped = 1;
      params->Y1 = -params->Y1;
      params->Y2 = -params->Y2;
    }
  else
    params->Flipped = 0;

  if (params->X2 > params->X1)
    params->Increment = 1;
  else
    params->Increment = -1;

  params->DeltaX=params->X2-params->X1;
  params->DeltaY=params->Y2-params->Y1;

  params->IncrE=2*params->DeltaY*params->Increment;
  params->IncrNE=2*(params->DeltaY-params->DeltaX)*params->Increment;
  params->DTerm=(2*params->DeltaY-params->DeltaX)*params->Increment;

  params->XIndex = params->X1;
  params->YIndex = params->Y1;
}

void get_current_point(bresenham_param_t *params, int *x, int *y)
{
  if (params->UsingYIndex)
    {
      *y = params->XIndex;
      *x = params->YIndex;
      if (params->Flipped)
        *x = -*x;
    }
  else
    {
      *x = params->XIndex;
      *y = params->YIndex;
      if (params->Flipped)
        *y = -*y;
    }
}

int get_next_point(bresenham_param_t *params)
{
  if (params->XIndex == params->X2)
    {
      return 0;
    }
  params->XIndex += params->Increment;
  if (params->DTerm < 0 || (params->Increment < 0 && params->DTerm <= 0))
    params->DTerm += params->IncrE;
  else
    {
      params->DTerm += params->IncrNE;
      params->YIndex += params->Increment;
    }
  return 1;
}


int IsValidLineSegment(double x0, double y0, double x1, double y1, double*  map,
       int x_size,
       int y_size)

{
  bresenham_param_t params;
  int nX, nY; 
    short unsigned int nX0, nY0, nX1, nY1;

    //printf("checking link <%f %f> to <%f %f>\n", x0,y0,x1,y1);
    
  //make sure the line segment is inside the environment
  if(x0 < 0 || x0 >= x_size ||
    x1 < 0 || x1 >= x_size ||
    y0 < 0 || y0 >= y_size ||
    y1 < 0 || y1 >= y_size)
    return 0;

  ContXY2Cell(x0, y0, &nX0, &nY0, x_size, y_size);
  ContXY2Cell(x1, y1, &nX1, &nY1, x_size, y_size);

    //printf("checking link <%d %d> to <%d %d>\n", nX0,nY0,nX1,nY1);

  //iterate through the points on the segment
  get_bresenham_parameters(nX0, nY0, nX1, nY1, &params);
  do {
    get_current_point(&params, &nX, &nY);
    if(map[GETMAPINDEX(nX,nY,x_size,y_size)] == 1)
            return 0;
  } while (get_next_point(&params));

  return 1;
}

int IsValidArmConfiguration(double* angles, int numofDOFs, double*  map,
       int x_size, int y_size)
{
    double x0,y0,x1,y1;
    int i;
    
  //iterate through all the links starting with the base
  x1 = ((double)x_size)/2.0;
  y1 = 0;
  for(i = 0; i < numofDOFs; i++)
  {
    //compute the corresponding line segment
    x0 = x1;
    y0 = y1;
    x1 = x0 + LINKLENGTH_CELLS*cos(2*PI-angles[i]);
    y1 = y0 - LINKLENGTH_CELLS*sin(2*PI-angles[i]);

    //check the validity of the corresponding line segment
    if(!IsValidLineSegment(x0,y0,x1,y1,map,x_size,y_size))
        return 0;
  }   

  return 1; 
}

#define DOF 5

// =============================================KD Tree=====================================================================
// KD tree node
struct kd_node_t{
  int index;
  double x[DOF];
  struct kd_node_t *left, *right;
};

// Initialize KD tree node 
void initialNode(struct kd_node_t *node, double* value){
  for(int i = 0; i < DOF; i++){
    node->x[i] = value[i];
  }
}

// KdTree class
class kdtree{
  public:
  inline double dist(struct kd_node_t *a, struct kd_node_t *b, int dim)
  {
    double t, d = 0;
    while (dim--) {
      t = a->x[dim] - b->x[dim];
      d += t * t;
    }
    return d;
  }

  inline void swap(struct kd_node_t *x, struct kd_node_t *y)
  {
    int index_temp;
    index_temp = x->index;
    x->index = y->index;
    y->index = index_temp;

    double tmp[DOF];
    memcpy(tmp, x->x, sizeof(tmp));
    memcpy(x->x, y->x, sizeof(tmp));
    memcpy(y->x, tmp, sizeof(tmp));
  }

  struct kd_node_t* find_median(struct kd_node_t *start, struct kd_node_t *end, int idx)
  {
    if (end <= start) return NULL;
    if (end == start + 1)
      return start;

    kd_node_t *p, *store, *md = start + (end - start)/2;
    double pivot;
    while (1) {
      pivot = md->x[idx];

      swap(md, end - 1);
      for(store = p = start; p < end; p++){
        if (p->x[idx] < pivot){
          if (p != store)
            swap(p, store);
          store++;
        }
      }
      swap(store, end - 1);

      if (store->x[idx] == md->x[idx])
        return md;

      if (store > md) 
        end = store;
      else
        start = store;
    }
  }

  struct kd_node_t* make_tree(struct kd_node_t *t, int len, int i, int dim)
  {
    struct kd_node_t *n;

    if (!len)
      return 0;
    if ((n = find_median(t, t + len, i))){
      i = (i + 1) % dim;
      n->left = make_tree(t, n - t, i, dim);
      n->right = make_tree(n + 1, t + len - (n + 1), i, dim);
    }
    return n;
  }

  // find nearest kd node 
  void nearest(struct kd_node_t *root, struct kd_node_t *nd, int i, int dim, struct kd_node_t **best, double *best_dist)
  {
    double d, dx, dx2;

    if(!root)
      return;
    d = dist(root, nd, dim);
    dx = root->x[i] - nd->x[i];
    dx2 = dx * dx;

    if(!*best || d < *best_dist){
      *best_dist = d;
      *best = root;
    }

    if (!*best_dist)
      return;

    if (++i >= dim)
      i = 0;

    nearest(dx > 0 ? root->left : root->right, nd, i, dim, best, best_dist);
    if (dx2 >= *best_dist)
      return;
    nearest(dx > 0 ? root->right : root->left, nd, i, dim, best, best_dist);
  }

  // find nodes in raidus of best_dist 
  void near(struct kd_node_t *root, struct kd_node_t *nd, int i, int dim, vector<struct kd_node_t *>& neighbor, double best_dist)
  {
    double d, dx, dx2;

    if(!root)
      return;
    d = sqrt(dist(root, nd, dim));
    dx = root->x[i] - nd->x[i];
    dx2 = fabs(dx);

    if( d < best_dist){
      neighbor.push_back(root);
    }

    if (++i >= dim)
      i = 0;

    near(dx > 0 ? root->left : root->right, nd, i, dim, neighbor, best_dist);
    if (dx2 >= best_dist)
      return;
    near(dx > 0 ? root->right : root->left, nd, i, dim, neighbor, best_dist);
  }

};

//====================================================Probabilistic RoadMap===================================================
// prm node 
struct prmNode{
  int index;
  double x[DOF];
};

// initialize prmNode 
void initialPrmNode(prmNode* node, double *value){
  node->index = 0;
  for(int i = 0; i < DOF; i++)
    node->x[i] = value[i];
}

bool comp(const pair<prmNode*, double> &a, const pair<prmNode*, double> &b)
{
    return a.second < b.second;
}

//======================PRM class==========================
class prm: public kdtree{
public:  
  int V; // number of Vertices 
  vector< list< pair<prmNode*, double> > > graph;  // roadMap
  struct kd_node_t *kdTree;
  struct kd_node_t *root;

public:
  prm(int VIn){
    V = VIn;
    // graph.reserve(V);
    kdTree = (struct kd_node_t*) calloc(V, sizeof(struct kd_node_t));
    root = NULL;
  };

/*  ~prm(){
    if(!graph.empty())
      graph.clear();

    free(kdTree);

    if(root)
      delete root;
  }*/

  // add Vertice to graph 
  void addVertice(prmNode *v){
    list< pair<prmNode*, double> > adj;
    adj.push_back(make_pair(v, 0));
    graph.push_back(adj);
  };

  // add kdNode 
  void addKdNode(struct kd_node_t *node){
    kdTree[node->index] = *node;
  }

  // add Edge to graph 
  void addEdge(prmNode *u, prmNode *v, double weight){
    graph[u->index].push_back(make_pair(v, weight));
    graph[v->index].push_back(make_pair(u, weight));
  };

  // compute dist between two prmNode 
  double dist(prmNode *u, prmNode *v){
    double dist = 0;
    for(int i = 0; i < DOF; i++)
      dist += (u->x[i] - v->x[i])*(u->x[i] - v->x[i]);
    dist = sqrt(dist);
    return dist;
  };

  // check if two prmNode is connected 
  bool checkConnected(int u, int v){
    bool *visited = new bool[V];

    for(int i = 0; i < V; i++)
      visited[i] = false;

    DFSUtil(u, visited);

    if(visited[v]){
      delete [] visited;
      return true;
    }
    else{
      delete [] visited;
      return false;
    }

  };

  // Do Depth First Search in the graph 
  void DFSUtil(int u, bool visited[]){
    if(u > (V - 1) || u < 0)
        cout << " big error ! " << endl; 
    visited[u] = true;
    // printf(" u %d\n",  u);

    list< pair<prmNode*, double> >::iterator iter;
    for(iter = graph[u].begin(); iter != graph[u].end(); iter++){
      if(!visited[(*iter).first->index])
        DFSUtil((*iter).first->index, visited);
    }

  };

  // Construct probabilistic roadmap 
  void process(double radius){

    root = make_tree(kdTree, (V-2), 0, DOF);
    
    for(int number = 0; number < (V-2); number++){

      if(number%1000 == 0)
        cout << "iteration " << number << " times " << endl;

      struct kd_node_t *node = new kd_node_t;
      node->index = number;
      initialNode(node, &(graph[number].front().first->x[0]));      
      
      // find neighbors in the radius
      vector<struct kd_node_t *> neighbor;
      near(root, node, 0, DOF, neighbor, radius);
      delete node;

      // connect to neighbors 
      if(neighbor.size() > 1){
        vector<pair<prmNode*, double> > nearVertice;
        vector<struct kd_node_t *>::iterator iter;
        for(iter = neighbor.begin(); iter != neighbor.end(); iter++){
          if((*iter)->index != number){
            double dist_temp = dist(graph[(*iter)->index].front().first, graph[number].front().first);
            nearVertice.push_back(make_pair(graph[(*iter)->index].front().first, dist_temp));
          }
        }

        std::sort(nearVertice.begin(), nearVertice.end(), comp);
        for(int i = 0; i < nearVertice.size(); i++){
            // bool result = checkConnected(nearVertice[i].first->index, number);
            // if(!result){
            // addEdge(nearVertice[i].first, graph[number].front().first, nearVertice[i].second);
            // }
          bool alreadyIn = false;
          list<pair<prmNode*, double>>::iterator iter;
          for(iter = graph[number].begin(); iter != graph[number].end(); iter++){
            if(nearVertice[i].first->index == (*iter).first->index)
              alreadyIn = true;
          }

          if(!alreadyIn)
              addEdge(nearVertice[i].first, graph[number].front().first, nearVertice[i].second);  
        }  
      }
    }
    
    return;

  };

  double test_delaunay_neighours(double radius)
  {
    root = make_tree(kdTree, V, 0, DOF);
    double best_dist;
    struct kd_node_t *node_near = 0;

    struct kd_node_t *node = new kd_node_t;

    //double center[5] = {PI,PI,PI,PI,PI};

    int rand_index = rand()%V;
    GRBEnv env = GRBEnv();

    //int rand_index = rand() % V;
    //struct kd_node_t *rand_node = kdTree+rand_index;
    node->index = V;
    initialNode(node, &(kdTree[rand_index].x[0]));
    nearest(root, node, 0, DOF, &node_near, &best_dist);

    //cout<<"the nearest node is:"<<node_near->x[0]<<" "<<node_near->x[1]<<" "<<node_near->x[2]<<" "<<node_near->x[3]<<" "<<node_near->x[4]<<endl;

    vector<struct kd_node_t *> neighbors;
    near(root, node_near, 0, DOF, neighbors, radius);

    cout<<"number of neighbors in radius are: "<<neighbors.size()<<endl;

    int num_dn = 0;

    for(auto& neighbor: neighbors)
    {
        if(checkDelaunayNeibour(node_near, neighbor, env))
        {
/*            if(checkDelaunayNeibour(neighbor, node_near, env))
                cout<<"test switch: "<<1<<endl;
            else
                cout<<"test switch: "<<0<<endl;*/
            num_dn++;
        }
    }

    cout<<"number of delaunay neighbors are: "<<num_dn<<endl;

    cout<<"probability for this test is: "<<(double)num_dn/neighbors.size()<<endl;
    return (double)num_dn/neighbors.size();

  };

  double getB(struct kd_node_t *node)
  {
    return node->x[0]*node->x[0]+node->x[1]*node->x[1]+node->x[2]*node->x[2]+node->x[3]*node->x[3]+node->x[4]*node->x[4]; 
  };
  

  bool checkDelaunayNeibour(struct kd_node_t *node1, struct kd_node_t *node2, GRBEnv env)
  {
    int index1 = node1->index;
    int index2 = node2->index;

    struct kd_node_t *node_t = new kd_node_t;

    //cout<<"node1:"<<node1->x[0]<<" "<<node1->x[1]<<" "<<node1->x[2]<<" "<<node1->x[3]<<" "<<node1->x[4]<<endl;
    //cout<<"node2:"<<node2->x[0]<<" "<<node2->x[1]<<" "<<node2->x[2]<<" "<<node2->x[3]<<" "<<node2->x[4]<<endl;
    //cout<<"test node2:"<<(kdTree+index2)->x[0]<<" "<<(kdTree+index2)->x[1]<<" "<<(kdTree+index2)->x[2]<<" "<<(kdTree+index2)->x[3]<<" "<<(kdTree+index2)->x[4]<<endl;

    

    if(index1 == index2) return true;

    try{
        GRBModel model = GRBModel(env);
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);
        
        // Create 
        GRBVar x0 = model.addVar(0.0, 2*PI, 0.0, GRB_CONTINUOUS, "x0");
        GRBVar x1 = model.addVar(0.0, 2*PI, 0.0, GRB_CONTINUOUS, "x1");
        GRBVar x2 = model.addVar(0.0, 2*PI, 0.0, GRB_CONTINUOUS, "x2");
        GRBVar x3 = model.addVar(0.0, 2*PI, 0.0, GRB_CONTINUOUS, "x3");
        GRBVar x4 = model.addVar(0.0, 2*PI, 0.0, GRB_CONTINUOUS, "x4");
        GRBVar x5 = model.addVar(0.0, std::numeric_limits<double>::infinity(), 0.0, GRB_CONTINUOUS, "x5");

        //cout<<"test getB: "<<getB(node1)<<endl;

        model.setObjective(getB(node1) - 2*node1->x[0]*x0 - 2*node1->x[1]*x1 - 2*node1->x[2]*x2 - 2*node1->x[3]*x3 - 2*node1->x[4]*x4+ x5 , GRB_MINIMIZE);
        model.addConstr(getB(node1) + 1 - 2*node1->x[0]*x0 - 2*node1->x[1]*x1 - 2*node1->x[2]*x2 - 2*node1->x[3]*x3 - 2*node1->x[4]*x4+ x5 >=0);


        for(int i = 0; i<V ; i++)
        {
            node_t = kdTree + i;
            if(node_t->index == index1) continue;
            if(node_t->index == index2)
                model.addConstr(getB(node2) - 2*node2->x[0]*x0 - 2*node2->x[1]*x1 - 2*node2->x[2]*x2 - 2*node2->x[3]*x3 - 2*node2->x[4]*x4+ x5 <=0);
            else
            {
                model.addConstr(getB(node_t) - 2*node_t->x[0]*x0 - 2*node_t->x[1]*x1 - 2*node_t->x[2]*x2 - 2*node_t->x[3]*x3 - 2*node_t->x[4]*x4+ x5 >=0);

            }
        }

        model.optimize();
        if(model.get(GRB_DoubleAttr_ObjVal)<0)
            return true;
        else
            return false;

    }catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    return false;
  };

};


//==================================================Test Delaunay Neighours============================================= 

double delaunay_neihour_prob(size_t V, double radius, int N, double *map)
{
  double sum_prob = 0;
  for(int i=0; i<N; i++)
  {
    cout<<"Begin "<<i+1<<"th Test"<<endl;
    int counter = 0;
    prm* myPRM = new prm(V);

    // Constructing PRM graph
    while(counter < V)
    {
        double node_rand[5] = { 2*PI*(rand()/double(RAND_MAX)), \
                                2*PI*(rand()/double(RAND_MAX)), \
                                2*PI*(rand()/double(RAND_MAX)), \
                                2*PI*(rand()/double(RAND_MAX)), \
                                2*PI*(rand()/double(RAND_MAX))};  

        //cout<<node_rand[0]<<" "<<node_rand[1]<<" "<<node_rand[2]<<" "<<node_rand[3]<<" "<<node_rand[4]<<endl;
        if(IsValidArmConfiguration(&node_rand[0], DOF, map, 50, 50))
        {
          struct kd_node_t *kdNode_new = new kd_node_t;
          kdNode_new->index = counter;
          initialNode(kdNode_new, &node_rand[0]);
          myPRM->addKdNode(kdNode_new);

          counter++;
        }
    }
    sum_prob = sum_prob + myPRM->test_delaunay_neighours(radius);

    delete myPRM;

  }

  return sum_prob/N;


}

int main(int argc,char *argv[])
{
  //For saving all the results
    double *map = GetMap("/home/zjcv2012/Planning_Project/prm/map1.txt");

    std::ofstream f;
    string filename = "/home/zjcv2012/Desktop/radius.txt";
    f.open(filename.c_str());
    f<<"Point Size"<<"      "<<"Radius"<<endl;

    srand (time(NULL));
    size_t start_V = 1000;
    size_t end_V = 200000;
    //size_t end_V = 10000;
    size_t incre_V = 1000;

    int N =10;

    double start_radius = 2.5;
    double incre_radius = 0.1;
    double radius;

    double PROB_THRESHOLD = 0.9;

    for(size_t V = start_V; V<=end_V; V+=incre_V)
    {
      radius = start_radius;
      if(radius<PI/10) break; //If radius is too small, then ends the program. 

      cout<<"Start calculating prob for V="<<V<<" and radius="<<radius<<endl;
      while(delaunay_neihour_prob(V,radius,N, map)<PROB_THRESHOLD)
      {
        radius -= incre_radius;
        if(radius<=1.5)
          incre_radius = 0.05;

        cout<<"Start calculating prob for V="<<V<<" and radius="<<radius<<endl;
      }
      start_radius = radius;
      cout<<"\nFor V size: "<<V<<" ,the radius is: "<<radius<<endl<<endl;
      f<<V<<"           "<<radius<<endl;
    
      if(V/incre_V>=10)
        incre_V*=2;
    }
    f.close();

    cout<<"Result is saved!"<<endl;

    //cout<<"probability is: "<<delaunay_neihour_prob(V, radius, N)<<endl;

    return 0; 
}



