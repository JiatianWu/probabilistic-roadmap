/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <vector>
#include <list>
#include <set>
#include <queue>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <chrono>

/* Input Arguments */
#define	MAP_IN      prhs[0]
#define	ARMSTART_IN	prhs[1]
#define	ARMGOAL_IN     prhs[2]
#define	PLANNER_ID_IN     prhs[3]

/* Planner Ids */
#define PRM         0
#define PRMDTC      1

/* Output Arguments */
#define	PLAN_OUT	plhs[0]
#define	PLANLENGTH_OUT	plhs[1]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654

#define INF 0x3f3f3f3f

//the length of each link in the arm (should be the same as the one used in runtest.m)
#define LINKLENGTH_CELLS 10

using namespace std;

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


int IsValidLineSegment(double x0, double y0, double x1, double y1, double*	map,
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

int IsValidArmConfiguration(double* angles, int numofDOFs, double*	map,
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

// Astar Node used to find best path 
class astarNode{
public:
  int index;
  double g;
  double h;
  double f;
  astarNode* parent;

public:
  astarNode(int indexIn, double gIn){
    index = indexIn;
    g = gIn;
    parent = NULL;
  };

  ~astarNode(){
    delete parent;
  };  

  bool operator<(const astarNode &a)const
  {
    return a.g < g;
  };

};

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

    kdTree = (struct kd_node_t*) calloc((V-2), sizeof(struct kd_node_t));
    root = NULL;
  };

  ~prm(){
    if(!graph.empty())
      graph.clear();

    free(kdTree);

    if(root)
      delete root;
  }

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
    printf("break point 0\n");

    for(int i = 0; i < V; i++)
      visited[i] = false;

    int counter = 0;
    DFSUtil(u, visited, counter);

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
  void DFSUtil(int u, bool visited[], int counter){
    printf("before ");
    if(u > (V - 1) || u < 0)
      printf("big error\n");
    printf(" %d\n", visited[u]);
    visited[u] = true;
    counter++;
    printf("counter %d u %d \n", counter, u);

    list< pair<prmNode*, double> >::iterator iter0;
    for(iter0 = graph[u].begin(); iter0 != graph[u].end(); iter0++){
      printf(" %d  %d ", (*iter0).first->index, visited[(*iter0).first->index]);
    }
    printf("\n");

    list< pair<prmNode*, double> >::iterator iter;
    for(iter = graph[u].begin(); iter != graph[u].end(); iter++){
      if(!visited[(*iter).first->index]){
        printf(" dfs %d %d\n", (*iter).first->index, visited[(*iter).first->index]);
        DFSUtil((*iter).first->index, visited, counter);
      }
    }

    return;

  };

  int computeComponent(){

    vector<bool> visited(V, false);
    
    int counter = 0;
    for(int i = 0; i < V; i++){
      if(!visited[i]){
        DFSUtil_backup(i, visited);
        counter++;
      }
    }

    return counter;
  }

    // check if two prmNode is connected 
  bool checkConnected_backup(int u, int v){

    vector<bool> visited(V, false);

    DFSUtil_backup(u, visited);

    if(visited[v]){
      return true;
    }
    else{
      return false;
    }

  };

  // Do Depth First Search in the graph 
  void DFSUtil_backup(int u, vector<bool> &visited){

    visited[u] = true;

    list< pair<prmNode*, double> >::iterator iter;
    for(iter = graph[u].begin(); iter != graph[u].end(); iter++){
      if(!visited[(*iter).first->index]){
        DFSUtil_backup((*iter).first->index, visited);
      }
    }

    return;

  };

  bool checkValid(int u, int v, double ecpilo, double* map, int x_size, int y_size){
    prmNode* start = graph[u].front().first;
    prmNode* goal = graph[v].front().first;
    double dist_temp = dist(start, goal);
    bool result = true;

    double lamda = ecpilo/dist_temp;
    double base = lamda;

    while(result && lamda < 1 ){
      double *node_new = new double[5];
      for(int i = 0; i < DOF; i++){
          node_new[i] = lamda*goal->x[i] + (1-lamda)*start->x[i];
      }
      if(!IsValidArmConfiguration(node_new, DOF, map, x_size, y_size)){
        result = false;
      }

      delete[] node_new;
      lamda = lamda + base;
    }

    return result;
  };

vector<int> generateInteger(int size, int max){
    set<int> s;
    vector<int> result;
    int number;

    while(1){
      int r = rand()%size;
      number = s.count(r);
      if(number == 0){
        s.insert(r);
        result.push_back(r);
      }
      if(result.size() == max){
        return result;
        break;
      }
    }
};

  // Construct probabilistic roadmap 
  void process(double radius, int maxLimit, double ecpilo, double *map, int x_size, int y_size){

    chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
    root = make_tree(kdTree, (V-2), 0, DOF);
    chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    double kdtree_time = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
    cout << " make kdtree time " << kdtree_time << endl;

    bool debugPrint = true;   
    for(int number = 0; number < (V-2); number++){

      if(number%3000 == 0 && debugPrint)
        cout << "processing " << number << " rounds " << endl;

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
        int neighborSize;
        if(neighbor.size() <= maxLimit)
          neighborSize = nearVertice.size();
        else
          neighborSize = maxLimit;

        for(int i = 0; i < neighborSize; i++){
          bool alreadyIn = false;
          list<pair<prmNode*, double>>::iterator iter;
          for(iter = graph[number].begin(); iter != graph[number].end(); iter++){
            if(nearVertice[i].first->index == (*iter).first->index)
              alreadyIn = true;
          }

          if(!alreadyIn){
            bool valid = checkValid(number, nearVertice[i].first->index, ecpilo, map, x_size, y_size);
            if(valid)
              addEdge(nearVertice[i].first, graph[number].front().first, nearVertice[i].second); 
          }
        }  
      }
    }
    
    return;

  };

  // Process Delaunay Triangulation Node Connection Strategy
  void processDTC(double radius, int maxLimit, double ecpilo, double *map, int x_size, int y_size){

    chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
    root = make_tree(kdTree, (V-2), 0, DOF);
    chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    double kdtree_time = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
    cout << " make kdtree time " << kdtree_time << endl;

    bool debugPrint = true;   
    for(int number = 0; number < (V-2); number++){

      if(number%3000 == 0 && debugPrint)
        cout << "processing " << number << " rounds " << endl;

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

        if(nearVertice.size() <= maxLimit){
          for(int i = 0; i < nearVertice.size(); i++){
            bool alreadyIn = false;
            list<pair<prmNode*, double>>::iterator iter;
            for(iter = graph[number].begin(); iter != graph[number].end(); iter++){
              if(nearVertice[i].first->index == (*iter).first->index)
                 alreadyIn = true;
            }

            if(!alreadyIn){
              bool valid = checkValid(number, nearVertice[i].first->index, ecpilo, map, x_size, y_size);
              if(valid)
                addEdge(nearVertice[i].first, graph[number].front().first, nearVertice[i].second); 
            }
          }
        }else{
          vector<int> random = generateInteger(nearVertice.size(), maxLimit);
          for(int i = 0; i < maxLimit; i++){
            int index_v = random[i];

            bool alreadyIn = false;
            list<pair<prmNode*, double>>::iterator iter;
            for(iter = graph[number].begin(); iter != graph[number].end(); iter++){
              if(nearVertice[index_v].first->index == (*iter).first->index)
                 alreadyIn = true;
            }
            
            if(!alreadyIn){
              bool valid = checkValid(number, nearVertice[index_v].first->index, ecpilo, map, x_size, y_size);
              if(valid)
                addEdge(nearVertice[index_v].first, graph[number].front().first, nearVertice[index_v].second); 
            }

          }
        }

      }
    }
    
    return;

  };

  // ADD startNode and goalNode to the graph 
  void addStartGoal(prmNode* node_start, kd_node_t* kdNode_start, prmNode* node_goal, kd_node_t* kdNode_goal, double radius, double ecpilo, double *map, int x_size, int y_size){
    
    assert(node_start->index == graph.size());
    addVertice(node_start);
    assert(node_goal->index == graph.size());
    addVertice(node_goal);

    for(int j = 0; j < 2; j++){

      struct kd_node_t *node = new kd_node_t;
      
      if(j == 0)
        node->index = node_start->index;
      else
        node->index = node_goal->index;

      initialNode(node, &(graph[node->index].front().first->x[0]));      

      vector<struct kd_node_t *> neighbor;
      neighbor.clear();
      near(root, node, 0, DOF, neighbor, radius);

      if(neighbor.size() > 1){
        // cout << " find neighrbor in radius " << endl;
        vector<pair<prmNode*, double> > nearVertice;
        vector<struct kd_node_t *>::iterator iter;
        for(iter = neighbor.begin(); iter != neighbor.end(); iter++){
          if((*iter)->index != node->index){
            double dist_temp = dist(graph[(*iter)->index].front().first, graph[node->index].front().first);
            nearVertice.push_back(make_pair(graph[(*iter)->index].front().first, dist_temp));
          }
        }

        std::sort(nearVertice.begin(), nearVertice.end(), comp);
        for(int i = 0; i < nearVertice.size(); i++){
          bool valid = checkValid(node->index, nearVertice[i].first->index, ecpilo, map, x_size, y_size);
          if(valid)
            addEdge(nearVertice[i].first, graph[node->index].front().first, nearVertice[i].second);
        }    
      }
      else{
        cout << " did not find neighbor in radius " << endl;
        double best_dist;
        struct kd_node_t *node_near = 0;
        nearest(root, node, 0, DOF, &node_near, &best_dist);
        double dist_temp = dist(graph[node_near->index].front().first, graph[node->index].front().first);
        addEdge(graph[node_near->index].front().first, graph[node->index].front().first, dist_temp);
      }

      delete node;

    }

    return;

  };

  astarNode *getLeastNode(list<astarNode *> open)  
  {  
      if(!open.empty())  
      {  
        auto resNode = open.front();  
        for(auto &node:open)  
            if(node->f < resNode->f)  
                resNode = node;  
        return resNode;  
      }  
      return NULL;  
  }; 

  double computeH(astarNode* temp, prmNode* node_goal){
    int index_temp = temp->index;
    double h;
    h = dist(graph[index_temp].front().first, node_goal);
    return h;
  }

  // Return the best path from start to goal if connected, 
  // if start and goal is not connected, show failure and return no path
  void findBestPath(prmNode* node_start, prmNode* node_goal, double*** plan, int* planlength){

    int numComponent = computeComponent();
    cout << " The generated roadmap has " << numComponent << " components."<< endl;

    bool result = checkConnected_backup(node_start->index, node_goal->index);
    // bool result = true;

    if(result){
      double weight = 2;
      double *cost = new double[V];
      bool *visit = new bool[V];   
      for(int i = 0; i < V; i++)
        {
            cost[i] = INF;
            visit[i] = false;
        }
    
      list<astarNode *> open;
      vector<astarNode *> close;
      astarNode *start = new astarNode(node_start->index, 0);
      cost[node_start->index] = 0;
      start->h = computeH(start, node_goal);
      start->f = start->g + weight * start->h;
      open.push_back(start);
        
      int counter = 0;  
      while(!open.empty())
      {
        astarNode *expand = getLeastNode(open);
        open.remove(expand);
        close.push_back(expand);

        if(visit[expand->index])
          continue;
        else 
          visit[expand->index] = true;
 
        if(visit[node_goal->index]){
          break;
        }

        counter++;
        if(counter % 1000 == 0)
           cout << counter << endl;

        list<pair<prmNode*, double>>::iterator iter;
        for(iter = graph[expand->index].begin(), iter++; iter != graph[expand->index].end(); iter++){
          double newG = expand->g + (*iter).second;
          if(!visit[(*iter).first->index] && cost[(*iter).first->index] > newG){  
            astarNode* temp = new astarNode((*iter).first->index, newG);
            temp->parent = expand;         
            cost[temp->index] = temp->g;
            temp->h = computeH(temp, node_goal);
            temp->f = temp->g + weight * temp->h;
            open.push_back(temp);
          }
        }  
      }    

      delete[] visit;
      visit = NULL;
      delete[] cost;
      cost = NULL;


      astarNode* goal = close.back();
      assert(goal->index == node_goal->index);
      list<int> path;
      while(goal != NULL){
        path.push_front(goal->index);
        goal = goal->parent;
      }

      *planlength = path.size();
  
      list<int>::iterator iter = path.begin();
      *plan = (double**) malloc((path.size())*sizeof(double*));
      for(int i = 0; i < *planlength; i++, iter++)
      {
         (*plan)[i] = (double*) malloc(DOF*sizeof(double));
         memcpy((*plan)[i], &(graph[(*iter)].front().first->x[0]), sizeof(graph[(*iter)].front().first->x));
      }

    }
    else{
      cout << " fail " << endl;  
    }
  };

};


//==================================================Probabilistic RoadMap Planner============================================= 
static void plannerPRM(
       double*  map,
       int x_size,
       int y_size,
           double* armstart_anglesV_rad,
           double* armgoal_anglesV_rad,
     int numofDOFs,
     double*** plan,
     int* planlength)
{
  //no plan by default
  *plan = NULL;
  *planlength = 0;


  double radius = 0.8;  // neighborhood radius size
  double V = 30000;  // sampling points number
  double ecpilo = PI/20;
  int maxLimit = 15;
  prm* myPRM = new prm(V);

  int counter = 0;
  while(counter < (V-2)){

    double node_rand[5] = { 2*PI*(rand()/double(RAND_MAX)), \
                            2*PI*(rand()/double(RAND_MAX)), \
                            2*PI*(rand()/double(RAND_MAX)), \
                            2*PI*(rand()/double(RAND_MAX)), \
                            2*PI*(rand()/double(RAND_MAX))};   

    if(IsValidArmConfiguration(&node_rand[0], DOF, map, x_size, y_size)){
        struct prmNode *node_new = new prmNode;
        initialPrmNode(node_new, &node_rand[0]);
        node_new->index = counter;
        myPRM->addVertice(node_new);

        struct kd_node_t *kdNode_new = new kd_node_t;
        kdNode_new->index = node_new->index;
        initialNode(kdNode_new, &node_rand[0]);
        myPRM->addKdNode(kdNode_new);
  
        counter++;
    }
  }

  struct prmNode *node_start = new prmNode; 
  initialPrmNode(node_start, armstart_anglesV_rad);
  node_start->index = myPRM->graph.size();
  struct kd_node_t *kdNode_start = new kd_node_t;
  kdNode_start->index = node_start->index;
  initialNode(kdNode_start, armstart_anglesV_rad);

  struct prmNode *node_goal = new prmNode; 
  initialPrmNode(node_goal, armgoal_anglesV_rad);
  node_goal->index = myPRM->graph.size() + 1;
  struct kd_node_t *kdNode_goal = new kd_node_t;
  kdNode_goal->index = node_goal->index;
  initialNode(kdNode_goal, armgoal_anglesV_rad);

  myPRM->process(radius, maxLimit, ecpilo, map, x_size, y_size);

  int edges = 0;
  vector<int> result;
  for(int i = 0; i < V-2; i++){
    edges += (myPRM->graph[i].size() - 1);
    result.push_back(myPRM->graph[i].size() - 1);
  }

  sort(result.begin(), result.end());
  cout << " The NNC based PRM generated " << edges  << " edges, "<< " the average neighbors per node is " << result[V/2 -1] << endl;

  myPRM->addStartGoal(node_start, kdNode_start, node_goal, kdNode_goal, radius, ecpilo, map, x_size, y_size);
  myPRM->findBestPath(node_start, node_goal, plan, planlength);

  return; 
}

//==================================================Delaunay Triangulation PRM============================================= 
static void plannerPRMDTC(
       double*  map,
       int x_size,
       int y_size,
           double* armstart_anglesV_rad,
           double* armgoal_anglesV_rad,
     int numofDOFs,
     double*** plan,
     int* planlength)
{
  //no plan by default
  *plan = NULL;
  *planlength = 0;

  double radius = 0.8;  // neighborhood radius size
  int V = 30000;  // sampling points number
  double ecpilo = PI/20;
  int maxLimit = 15;
  prm* myPRM = new prm(V);

  int counter = 0;
  while(counter < (V-2)){

    double node_rand[5] = { 2*PI*(rand()/double(RAND_MAX)), \
                            2*PI*(rand()/double(RAND_MAX)), \
                            2*PI*(rand()/double(RAND_MAX)), \
                            2*PI*(rand()/double(RAND_MAX)), \
                            2*PI*(rand()/double(RAND_MAX))};   

    if(IsValidArmConfiguration(&node_rand[0], DOF, map, x_size, y_size)){
        struct prmNode *node_new = new prmNode;
        initialPrmNode(node_new, &node_rand[0]);
        node_new->index = counter;
        myPRM->addVertice(node_new);

        struct kd_node_t *kdNode_new = new kd_node_t;
        kdNode_new->index = node_new->index;
        initialNode(kdNode_new, &node_rand[0]);
        myPRM->addKdNode(kdNode_new);
  
        counter++;
    }
  }

  struct prmNode *node_start = new prmNode; 
  initialPrmNode(node_start, armstart_anglesV_rad);
  node_start->index = myPRM->graph.size();
  struct kd_node_t *kdNode_start = new kd_node_t;
  kdNode_start->index = node_start->index;
  initialNode(kdNode_start, armstart_anglesV_rad);

  struct prmNode *node_goal = new prmNode; 
  initialPrmNode(node_goal, armgoal_anglesV_rad);
  node_goal->index = myPRM->graph.size() + 1;
  struct kd_node_t *kdNode_goal = new kd_node_t;
  kdNode_goal->index = node_goal->index;
  initialNode(kdNode_goal, armgoal_anglesV_rad);

  myPRM->processDTC(radius, maxLimit, ecpilo, map, x_size, y_size);

  int number = 0;
  vector<int> result;
  for(int i = 0; i < V-2; i++){
    number += myPRM->graph[i].size() - 1;
    result.push_back(myPRM->graph[i].size() - 1);
  }
  
  sort(result.begin(), result.end());
  cout << " The DTC based PRM generated " << number  << " edges, "<< " the average neighbors per node is " << result[V/2 -1] << endl;

  myPRM->addStartGoal(node_start, kdNode_start, node_goal, kdNode_goal, radius, ecpilo, map, x_size, y_size);
  myPRM->findBestPath(node_start, node_goal, plan, planlength);

  return; 
};

vector<int> generateInteger_backup(int size, int max){
    set<int> s;
    vector<int> result;
    int number;

    while(1){
      int r = rand()%size;
      number = s.count(r);
      if(number == 0){
        s.insert(r);
        result.push_back(r);
      }
      if(result.size() == max){
        return result;
        break;
      }
    }
};

static void test(
       double*  map,
       int x_size,
       int y_size,
           double* armstart_anglesV_rad,
           double* armgoal_anglesV_rad,
     int numofDOFs,
     double*** plan,
     int* planlength){

    *plan = NULL;
    *planlength = 0;

    int i;
    struct kd_node_t wp[] = {
        {0, {2, 3}}, {1, {5, 4}}, {2, {9, 6}}, {3, {4, 7}}, {4, {8, 1}}, {5, {7, 2}}
    };
    struct kd_node_t testNode = {6, {7, 2}};
    struct kd_node_t *root;
    double best_dist = 1;
    vector<struct kd_node_t *> neighbor;

    kdtree* kdtree;

    root = kdtree->make_tree(wp, sizeof(wp) / sizeof(wp[1]), 0, 2);
 
    kdtree->near(root, &testNode, 0, 2, neighbor, best_dist);

    cout << " neighbor size " << neighbor.size() << endl;

    for(int i = 0; i < neighbor.size(); i++)
      cout << neighbor[i]->index << "  " << neighbor[i]->x[0] << "  " << neighbor[i]->x[1] << endl;

    // srand(time(NULL));
    vector<int> result;
    result = generateInteger_backup(50, 15);
    for(int i = 0; i < result.size(); i++)
      cout << result[i] << endl;

    return;

};

//prhs contains input parameters (3): 
//1st is matrix with all the obstacles
//2nd is a row vector of start angles for the arm 
//3nd is a row vector of goal angles for the arm 
//plhs should contain output parameters (2): 
//1st is a 2D matrix plan when each plan[i][j] is the value of jth angle at the ith step of the plan
//(there are D DoF of the arm (that is, D angles). So, j can take values from 0 to D-1
//2nd is planlength (int)
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[])
     
{ 
    
    /* Check for proper number of arguments */    
    if (nrhs != 4) { 
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidNumInputs",
                "Four input arguments required."); 
    } else if (nlhs != 2) {
	    mexErrMsgIdAndTxt( "MATLAB:planner:maxlhs",
                "One output argument required."); 
    } 
        
    /* get the dimensions of the map and the map matrix itself*/     
    int x_size = (int) mxGetM(MAP_IN);
    int y_size = (int) mxGetN(MAP_IN);
    double* map = mxGetPr(MAP_IN);
    
    /* get the start and goal angles*/     
    int numofDOFs = (int) (MAX(mxGetM(ARMSTART_IN), mxGetN(ARMSTART_IN)));
    if(numofDOFs <= 1){
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidnumofdofs",
                "it should be at least 2");         
    }
    double* armstart_anglesV_rad = mxGetPr(ARMSTART_IN);
    if (numofDOFs != MAX(mxGetM(ARMGOAL_IN), mxGetN(ARMGOAL_IN))){
        	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidnumofdofs",
                "numofDOFs in startangles is different from goalangles");         
    }
    double* armgoal_anglesV_rad = mxGetPr(ARMGOAL_IN);
 
    //get the planner id
    int planner_id = (int)*mxGetPr(PLANNER_ID_IN);
    if(planner_id < 0 || planner_id > 3){
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidplanner_id",
                "planner id should be between 0 and 3 inclusive");         
    }
    
    //call the planner
    double** plan = NULL;
    int planlength = 0;
    
    //you can may be call the corresponding planner function here

    if (planner_id == PRM)
    {
       plannerPRM(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
    }

    if (planner_id == PRMDTC)
    {
       plannerPRMDTC(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
    }

    // test(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);

    printf("planner returned plan of length=%d\n", planlength); 
    
    /* Create return values */
    if(planlength > 0)
    {
        PLAN_OUT = mxCreateNumericMatrix( (mwSize)planlength, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL); 
        double* plan_out = mxGetPr(PLAN_OUT);        
        //copy the values
        int i,j;
        for(i = 0; i < planlength; i++)
        {
            for (j = 0; j < numofDOFs; j++)
            {
                plan_out[j*planlength + i] = plan[i][j];
            }
        }
    }
    else
    {
        PLAN_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL); 
        double* plan_out = mxGetPr(PLAN_OUT);
        //copy the values
        int j;
        for(j = 0; j < numofDOFs; j++)
        {
                plan_out[j] = armstart_anglesV_rad[j];
        }     
    }
    PLANLENGTH_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)1, mxINT8_CLASS, mxREAL); 
    int* planlength_out = (int*)mxGetPr(PLANLENGTH_OUT);
    *planlength_out = planlength;

    return;
}