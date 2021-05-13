/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"
#include <unordered_map>
#include <list>
#include <stdlib.h>
#include <vector>
#include <chrono>

/* Input Arguments */
#define	MAP_IN      prhs[0]
#define	ARMSTART_IN	prhs[1]
#define	ARMGOAL_IN     prhs[2]
#define	PLANNER_ID_IN     prhs[3]

/* Planner Ids */
#define RRT         0
#define RRTCONNECT  1
#define RRTSTAR     2
#define PRM         3

/* Output Arguments */
#define	PLAN_OUT	plhs[0]
#define	PLANLENGTH_OUT	plhs[1]
#define NODES_GENERATED_OUT  plhs[2]
#define PLANTIME_OUT  plhs[3]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)
//#define GETMAPINDEX3D(X, Y, Z, XSIZE, YSIZE, ZSIZE) (Z*(XSIZE*YSIZE) + Y*XSIZE + X)
#define GETMAPINDEX3D(X, Y, Z, XSIZE, YSIZE, ZSIZE) (Z*(XSIZE*YSIZE) + X*YSIZE + Y)

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654

//the length of each link in the arm (should be the same as the one used in runtest.m)
#define LINKLENGTH_CELLS 10

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

int IsValidArmConfiguration(double* angles, int numofDOFs, double* map,
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

// simplified version of the infrastructure above checks if joints are in valid spots
bool IsValidArmConfiguration3D(double* angles, int numofDOFs, double* map, 
							int x_size, int y_size, int z_size)
{
	double l0 = 0.1;
	double l1 = 0.35;
	double l2 = 0.35;
	double l3 = 0.125;

	// this should all be done in a loop
	// sorry it's so gross
	// I was in too deep to change anything, running low on time :(
	double th0 = angles[0];
	double th1 = angles[1];
	double th2 = angles[2];
	double th3 = angles[3];
	double xpos = angles[4];
	double ypos = angles[5];
	double thpos = angles[6];

	// these will correspond to grid coordinates, must multiply by 20
	double m = 20;
	double r0,r1,r2,r3,r4;
	double x0,x1,x2,x3,x4;
	double y0,y1,y2,y3,y4;
	double z0,z1,z2,z3,z4;

	r0 = 0; 
    r1 = 0;
    r2 = (l1*cos(th1));
    r3 = (r2 + l2*cos(th1+th2));
    r4 = (r3 + l3*cos(th1+th2+th3));

    x0 = (0 + xpos);
    x1 = (0 + xpos);
    x2 = ((r2-r1)*cos(th0+thpos) + xpos);
    x3 = ((r3-r1)*cos(th0+thpos) + xpos);
    x4 = ((r4-r1)*cos(th0+thpos) + xpos);

    y0 = (0 + ypos);
    y1 = (0 + ypos);
    y2 = ((r2-r1)*sin(th0+thpos) + ypos);
    y3 = ((r3-r1)*sin(th0+thpos) + ypos);
    y4 = ((r4-r1)*sin(th0+thpos) + ypos);

    z0 = 0;
    z1 = (z0 + l0);
    z2 = (z1 + l1*sin(th1));
    z3 = (z2 + l2*sin(th1 + th2));
    z4 = (z3 + l3*sin(th1 + th2 + th3));

    std::list<int> indices;

    /*mexPrintf("x0,y0,z0: %f,%f,%f\n", m*x0, m*y0, m*z0);
    mexPrintf("x1,y1,z1: %f,%f,%f\n", m*x1, m*y1, m*z1);
	mexPrintf("x2,y2,z2: %f,%f,%f\n", m*x2, m*y2, m*z2);
	mexPrintf("x3,y3,z3: %f,%f,%f\n", m*x3, m*y3, m*z3);
	mexPrintf("x4,y4,z4: %f,%f,%f\n", m*x4, m*y4, m*z4);*/

    int x0i = (int) (m*x0);
    int x1i = (int) (m*x1);
    int x2i = (int) (m*x2);
    int x3i = (int) (m*x3);
    int x4i = (int) (m*x4);

    int y0i = (int) (m*y0);
    int y1i = (int) (m*y1);
    int y2i = (int) (m*y2);
    int y3i = (int) (m*y3);
    int y4i = (int) (m*y4);

    int z0i = (int) (m*z0);
    int z1i = (int) (m*z1);
    int z2i = (int) (m*z2);
    int z3i = (int) (m*z3);
    int z4i = (int) (m*z4);

    // manually check bounds
    if (z0i < 0 || z1i < 0 || z2i < 0 || z3i < 0 || z4i < 0 ||
    	y0i < 0 || y1i < 0 || y2i < 0 || y3i < 0 || y4i < 0 ||
    	x0i < 0 || x1i < 0 || x2i < 0 || x3i < 0 || x4i < 0 ||
    	z0i >= z_size || z1i >= z_size || z2i >= z_size || z3i >= z_size || z4i >= z_size ||
    	y0i >= y_size || y1i >= y_size || y2i >= y_size || z3i >= y_size || z4i >= y_size ||
    	x0i >= x_size || x1i >= x_size || x2i >= x_size || z3i >= x_size || z4i >= x_size ) 
    {
    	return false;
    }

    /*mexPrintf("x0i,y0i,z0i: %d,%d,%d\n", x0i, y0i, z0i);
    mexPrintf("x1i,y1i,z1i: %d,%d,%d\n", x1i, y1i, z1i);
	mexPrintf("x2i,y2i,z2i: %d,%d,%d\n", x2i, y2i, z2i);
	mexPrintf("x3i,y3i,z3i: %d,%d,%d\n", x3i, y3i, z3i);*/
	//mexPrintf("x4i,y4i,z4i: %d,%d,%d\n", x4i, y4i, z4i);

    indices.push_back(GETMAPINDEX3D(x0i, y0i, z0i, x_size, y_size, z_size));
    indices.push_back(GETMAPINDEX3D(x1i, y1i, z1i, x_size, y_size, z_size));
    indices.push_back(GETMAPINDEX3D(x2i, y2i, z2i, x_size, y_size, z_size));
    indices.push_back(GETMAPINDEX3D(x3i, y3i, z3i, x_size, y_size, z_size));
    indices.push_back(GETMAPINDEX3D(x4i, y4i, z4i, x_size, y_size, z_size));

    for (int index : indices) 
    {
    	int val = (int) map[index];
    	if (val == 1) {
    		return false;
    	}
    }
    //mexPrintf("\n");
    return true;
}

// use 2-norm of vector difference for distance
double squaredDistance(double* v1, double* v2, int numofDOFs) {
	double dist = 0;
	for (int i = 0; i < numofDOFs; i++) {
		dist = dist + (v2[i] - v1[i])*(v2[i] - v1[i]);
		//mexPrintf("dist: %f\n", dist);
	}
	dist = sqrt(dist);
	return dist;
}

double* randomConfig(int numofDOFs) {
	double* v = (double*) malloc(numofDOFs*sizeof(double));
	double mag = 2*PI;
	double direction;
	for (int i = 0; i < numofDOFs; i++) {
		int directionInt = rand();
		if (directionInt % 2 == 0) direction = -1.; else direction = 1.; 
		double r = (((double)rand())/((double)RAND_MAX))*mag*direction;
		v[i] = r;
	}
	return v;
}

double* randomConfig3D(int numofDOFs) {
	double* v = (double*) malloc(numofDOFs*sizeof(double));
	double r;
	double mag = 2*PI;
	double magxy = 5;
	double direction;

	for (int i = 0; i < numofDOFs; i++) {
		int directionInt = rand();
		if (directionInt % 2 == 0) direction = -1.; else direction = 1.; 
		if (i == 4 || i == 5) { // hardcode logic to move in x-y plane
			r = (((double)rand())/((double)RAND_MAX))*magxy*direction;
		}
		else { // do the usual arm thing
			r = (((double)rand())/((double)RAND_MAX))*mag*direction;	
		}
		v[i] = r;
	}
	return v;
}


int nearestNode(std::unordered_map <int, double*>* nodes, double* qrand, int tree_id, int numofDOFs) {
	// search current tree for closest node
	int closestNode = -1;
	double minDist = 1000;
	for (int i = 0; i < tree_id; i++) {
		double* q = (*nodes)[i];
		double d = squaredDistance(qrand, q, numofDOFs);
		if (d < minDist) {
			minDist = d;
			closestNode = i;
		}
	}
	return closestNode; 
}

// will return 0 for "trapped", 1 for "advanced", 2 for "reached"
int extend(std::unordered_map <int, double*>* nodes, std::unordered_map <int, int>* tree, 
			double* qrand, int numofDOFs, int* tree_id, double eps, 
			double* map, int x_size, int y_size, int z_size) 
{
	int closestNode = nearestNode(nodes, qrand, *tree_id, numofDOFs);
	double* qclose = (*nodes)[closestNode];

	/*if (IsValidArmConfiguration3D(qrand, numofDOFs, map, x_size, y_size, z_size)) {
		mexPrintf("valid\n");
	} else {
		mexPrintf("invalid\n");
	}*/

	if (squaredDistance(qclose, qrand, numofDOFs) <= eps) {
		if (IsValidArmConfiguration3D(qrand, numofDOFs, map, x_size, y_size, z_size)) {
			double dist = squaredDistance(qclose, qrand, numofDOFs);
			(*nodes)[*tree_id] = qrand; // create node
			(*tree)[*tree_id] = closestNode; // add node to current tree
			return 2;	
		}
	} 
	else {
		// interpolate for each DOF
		int direction;
		double* qnew = (double*) malloc(numofDOFs*sizeof(double));
		for (int i = 0; i < numofDOFs; i++) {
			if (qrand[i] >= qclose[i]) direction = 1; else direction = -1;
			qnew[i] = (qrand[i] - qclose[i])*eps*(1/squaredDistance(qclose, qrand, numofDOFs)) + qclose[i];
		}
		if (IsValidArmConfiguration3D(qnew, numofDOFs, map, x_size, y_size, z_size)) {
			(*nodes)[*tree_id] = qnew;
			(*tree)[*tree_id] = closestNode;
			return 1;
		}
		else {
			// mexPrintf("trapped\n");
			return 0;
		}
	}
}

// connect attempts to extend from begin to target
// will return 0 for "trapped", 2 for "reached"
int connect(std::unordered_map <int, double*>* nodes_target, std::unordered_map <int, double*>* nodes_begin, 
			std::unordered_map <int, int>* tree_target, std::unordered_map <int, int>* tree_begin,
			int* tree_id_target, int* tree_id_begin,
			int* connect_target_node, int* connect_begin_node,
			int numofDOFs, double eps, double* map, int x_size, int y_size, int z_size)
{	
	int stopped = 1;
	int iters = 0;
	while (stopped == 1) {

		double* q = (*nodes_target)[*tree_id_target];
		stopped = extend(nodes_begin, tree_begin, q, numofDOFs, tree_id_begin, eps, map, x_size, y_size, z_size);
		//mexPrintf("stopped: %d\n", stopped);
		//mexPrintf("tree_begin[tree_id_begin]: %d\n", (*tree_begin)[*tree_id_begin]);
		if (stopped == 2) {
			*connect_target_node = *tree_id_target;
			*connect_begin_node = *tree_id_begin;
			return 2;
		}
		if (stopped == 1) (*tree_id_begin)++;
	}
	return 0;
}

double** linearInterp(double* q1, double* q2, int numPoints, int numofDOFs, double* map, int x_size, int y_size, int z_size) {
	double** interp  = (double**) malloc(numPoints*sizeof(double*));
	for (int i = 0; i < numPoints; i++) {
		//mexPrintf("i: %d\n", i);
		double* config = (double*) malloc(numofDOFs*sizeof(double));
		for (int j = 0; j < numofDOFs; j++) {
			//mexPrintf("j: %d\n", j);
			config[j] = (q2[j] - q1[j])*(((double)i)/((double)(numPoints-1.))) + q1[j];
		}
		if (!IsValidArmConfiguration3D(config, numofDOFs, map, x_size, y_size, z_size)) {
			interp = NULL;
			return interp;
		}
		else {
			interp[i] = config;
		}
	}
	return interp;	
}

static void postProcessRRTCONNECT(double* map, int x_size, int y_size, int z_size, double*** plan, int* planlength, int numofDOFs) {
	int beginSearch = 0;
	int endSearch = (*planlength) - 1;
	while (beginSearch != endSearch) {
		double** interp = linearInterp((*plan)[beginSearch], (*plan)[endSearch], endSearch - beginSearch + 1, numofDOFs, map, x_size, y_size, z_size); 
		if (interp != NULL) {
			for (int i = beginSearch; i < endSearch + 1; i++) {
				(*plan)[i] = interp[i - beginSearch];
			}
			beginSearch = endSearch;
			endSearch = *planlength - 1;
		}
		else {
			endSearch--;
		}
	}
}

static void plannerRRTCONNECT(double* map,
			int x_size, 
			int y_size,
			int z_size,
			double* armstart_anglesV_rad,
			double* armgoal_anglesV_rad,
		int numofDOFs,
		double*** plan,
		int* planlength,
		int* nodes_generated,
		double* plantime)
{

	mexPrintf("x size: %d\n", x_size);
	mexPrintf("y size: %d\n", y_size);
	mexPrintf("z size: %d\n", z_size);
	// for (int i = 0; i < x_size*y_size*z_size; i++) {
 //    	double num = map[i];
 //    	mexPrintf("map: %f\n", num);
 //    }

 //    mexPrintf("hi\n");

	//start timer
	auto start = std::chrono::high_resolution_clock::now();

	// check start/goal configs
	bool bad_start = false;
	bool bad_goal = false;
	if (!IsValidArmConfiguration3D(armstart_anglesV_rad, numofDOFs, map, x_size, y_size, z_size)) {
		bad_start = true;
		mexPrintf("Invalid start configuration!\n");
	}
	if (!IsValidArmConfiguration3D(armgoal_anglesV_rad, numofDOFs, map, x_size, y_size, z_size)) {
		bad_goal = true;
		mexPrintf("Invalid goal configuration!\n");
	}
	if (bad_start || bad_goal) {
		mexPrintf("Aborting\n");
		return;
	}

	// assign each node a unique integer id
	std::unordered_map <int, double*> nodes_start; 
	std::unordered_map <int, double*> nodes_goal;
	std::unordered_map <int, int> t_start; // start tree; maps node:parent
	std::unordered_map <int, int> t_goal;  //  goal tree; maps node:parent

	nodes_start[0] = armstart_anglesV_rad;
	nodes_goal[0] = armgoal_anglesV_rad;
	t_start[0] = -1; // "the parent of node 1 is -1"
	t_goal[0] = -1;	// "the parent of node 2 is -2"

	int on_tree = 0; // 0 for t_start, 1 for t_goal;
	int tree_id_start = 1; // for any given time, tree_id should represent the next id to be added to the tree (not in nodes yet --> unsafe to index nodes[tree_id])
	int tree_id_goal = 1;

	int connect_start_node = -1; // t_start node at which trees connect
	int connect_goal_node = -1; // t_goal node at which trees connect

	double eps = 5*PI/180.; // 1 degree
	bool connected = false;

	while (!connected) {

		double* qrand = randomConfig3D(numofDOFs);
		//mexPrintf("qrand: %f,%f,%f,%f,%f,%f,%f\n", qrand[0], qrand[1], qrand[2], qrand[3], qrand[4], qrand[5], qrand[6]);

		// implicitly swap
		if (on_tree == 0) { // t_start
			if (extend(&nodes_start, &t_start, qrand, numofDOFs, &tree_id_start, eps, map, x_size, y_size, z_size) != 0) {
				if (connect(&nodes_start, &nodes_goal, &t_start, &t_goal, &tree_id_start, &tree_id_goal, &connect_start_node, &connect_goal_node, 
							numofDOFs, eps, map, x_size, y_size, z_size) == 2) 
				{
					connected = true;
				}
				tree_id_start++; // connect should increment tree_id_goal
			}
			on_tree = 1;
		}	
		else {
			if (extend(&nodes_goal, &t_goal, qrand, numofDOFs, &tree_id_goal, eps, map, x_size, y_size, z_size) != 0) {
				if (connect(&nodes_goal, &nodes_start, &t_goal, &t_start, &tree_id_goal, &tree_id_start, &connect_goal_node, &connect_start_node, 
							numofDOFs, eps, map, x_size, y_size, z_size) == 2) 
				{
					connected = true;
				}
				tree_id_goal++; // connect should increment tree_id_start
			}
			on_tree = 0;
		}
	}

	// backtrace the path
	int count_start = 1; // this counts the shared mid-path nodes and nodes in the start tree path
	int count_goal = 0; // this counts only nodes in the goal tree path (not the shared node)
	int i = connect_start_node;
	int j = connect_goal_node;
	while (i > 0) {
		i = t_start[i];
		count_start++;
	}
	while (j > 0) {
		j = t_goal[j];
		count_goal++;
	}
	*planlength = count_start + count_goal;
	*nodes_generated = tree_id_goal + tree_id_start;
	*plan = (double**) malloc((*planlength)*sizeof(double*));

	i = connect_start_node;
	j = t_goal[connect_goal_node]; // shared node included in first loop
	int plan_index = count_start - 1;
	while (i != -1) {
		(*plan)[plan_index] = (double*) malloc(numofDOFs*sizeof(double));
		(*plan)[plan_index] = nodes_start[i];
		plan_index--;
		i = t_start[i];
	}
	plan_index = count_start; 
	while  (j != -1) {
		(*plan)[plan_index] = (double*) malloc(numofDOFs*sizeof(double));
		(*plan)[plan_index] = nodes_goal[j];
		plan_index++;
		j = t_goal[j];
	}
	postProcessRRTCONNECT(map, x_size, y_size, z_size, plan, planlength, numofDOFs);

	// end timer
	auto stop = std::chrono::high_resolution_clock::now();
  	double planTime = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
  	*plantime = planTime;
  	mexPrintf("function time: %f\n", planTime);
  	mexPrintf("nodes generated: %d\n", *nodes_generated);
}

static void planner(
		   double*	map,
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
    
    //for now just do straight interpolation between start and goal checking for the validity of samples

    double distance = 0;
    int i,j;
    for (j = 0; j < numofDOFs; j++){
        if(distance < fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]))
            distance = fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]);
    }
    int numofsamples = (int)(distance/(PI/20));
    if(numofsamples < 2){
        printf("the arm is already at the goal\n");
        return;
    }
    *plan = (double**) malloc(numofsamples*sizeof(double*));
    int firstinvalidconf = 1;
    for (i = 0; i < numofsamples; i++){
        (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
        for(j = 0; j < numofDOFs; j++){
            (*plan)[i][j] = armstart_anglesV_rad[j] + ((double)(i)/(numofsamples-1))*(armgoal_anglesV_rad[j] - armstart_anglesV_rad[j]);
        }
        if(!IsValidArmConfiguration((*plan)[i], numofDOFs, map, x_size, y_size) && firstinvalidconf)
        {
            firstinvalidconf = 1;
            printf("ERROR: Invalid arm configuration!!!\n");
        }
    }    
    *planlength = numofsamples;
    
    return;
}

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
    } else if (nlhs != 4) {
	    mexErrMsgIdAndTxt( "MATLAB:planner:maxlhs",
                "One output argument required."); 
    } 
        
    /* get the dimensions of the map and the map matrix itself*/ 
    const mwSize* dims = mxGetDimensions(MAP_IN); // will be of length 3
    int x_size = (int) dims[0];
    int y_size = (int) dims[1];
    int z_size = (int) dims[2];
    double* map = mxGetPr(MAP_IN);

    /*int x_size = (int) mxGetM(MAP_IN);
    int y_size = (int) mxGetN(MAP_IN);
    double* map = mxGetPr(MAP_IN);*/
    
    /* get the start and goal angles*/     
    // int numofDOFs = (int) (MAX(mxGetM(ARMSTART_IN), mxGetN(ARMSTART_IN)));
    int numofDOFs = 7;
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
    int nodes_generated = 0;
    double plantime = 0;
    
    //you can may be call the corresponding planner function here
    if (planner_id == RRTCONNECT)
    {
        plannerRRTCONNECT(map, x_size, y_size, z_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength, &nodes_generated, &plantime);
    }
    else {
    //dummy planner which only computes interpolated path
    planner(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength); 
	}	
    
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
    int* planlength_out = (int*) mxGetPr(PLANLENGTH_OUT);
    *planlength_out = planlength;

    NODES_GENERATED_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)1, mxINT8_CLASS, mxREAL);
    int* nodes_generated_out = (int*) mxGetPr(NODES_GENERATED_OUT);
    *nodes_generated_out = nodes_generated;

    PLANTIME_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)1, mxDOUBLE_CLASS, mxREAL);
    double* plantime_out = (double*) mxGetPr(PLANTIME_OUT);
    *plantime_out = plantime;
    
    return;
    
}





