/**
Written by Albert Gunneström
2019-03-20
*/
#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <omp.h>
#include <limits>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdbool.h>

#include "QuadTreeNode.h"

// #define G 6.67e-11
#define DIMENSION 100
#define TIMESTEP 0.1
#define BOUNDARY 500. // The limits of the barnes hut root square

double G = 6.67e-11;
int num_planets;

using namespace std;

struct point {
  double x; //position
  double y;
  double vx; // velocity
  double vy;
  double fx; // force
  double fy;
  double m; // mass
};

// struct node {
//   int has_particles; // 1 if it contains particle
//   int has_ipoint;
//   point ipoint;
//   int level;
//   double xlow;
//   double xhigh;
//   double ylow;
//   double yhigh;
//
// /* all the children */
//   node *nw;
//   node *ne;
//   node *sw;
//   node *se;
// };



/* benchmark code */
double read_timer() {
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

void calcForces(std::vector<point> &p){ // p contains all points
  double distance, magnitude, direction;
  double xd, yd; // partial directions
  double e = 0.001; // margin to avoid zero division
  for(unsigned i = 0; i < p.size() - 1; i++) {
    for(unsigned j = i + 1; j < p.size(); j++) {
      distance = std::sqrt( std::pow((p[i].x - p[j].x), 2 ) +
  			                    std::pow((p[i].y - p[j].y),2) );
      if(distance < e){ // avoids dividing by zero
        distance = e;
        // printf("too close the distance was %lf\n", distance);
      }

      magnitude = (G * p[i].m * p[j].m) / std::pow(distance, 2);

      xd = p[j].x - p[i].x; // direction with respect to x
      yd = p[j].y - p[i].y; // direction with respect to y
      // printf("did we get here?\n");
      /* update forces */
      p[i].fx = p[i].fx + magnitude * xd / distance; // make not inf
      p[j].fx = p[j].fx - magnitude * xd / distance;
      p[i].fy = p[i].fy + magnitude * yd / distance;
      p[j].fy = p[j].fy - magnitude * yd / distance;
    }
  }
}

void moveBodies(std::vector<point> &p) {
  double dvx, dvy, dpx, dpy; // partial velocities and positions
  for(unsigned i = 0; i < p.size(); i++) {
    dvx = p[i].fx / p[i].m * TIMESTEP;
    dvy = p[i].fy / p[i].m * TIMESTEP;
    dpx = (p[i].vx + dvx/2) * TIMESTEP;
    dpy = (p[i].vy + dvy/2) * TIMESTEP;

    p[i].vx = p[i].vx + dvx; // change velocity
    p[i].vy = p[i].vy + dvy;
    p[i].x = p[i].x + dpx; // change position
    p[i].y = p[i].y + dpy;

    // cout << "node "<< i << " has position [" << p[i].x << "][" << p[i].y << "]" << endl;
    double color = i / (double)num_planets;
    // printf("%lf %lf %lf\n", p[i].x, p[i].y, color);

    p[i].fx = p[i].fy = 0.0; //reset force vector

  }

}

point* createBody(double xp, double yp, double vx, double vy, double fx, double fy,
                  double m, std::vector<point> &bodies) {
  point* newPoint = new point;
  newPoint->x = xp;
  newPoint->y = yp;
  newPoint->vx = vx;
  newPoint->vy = vy;
  newPoint->fx = fx;
  newPoint->fy = fy;
  newPoint->m = m;
  bodies.push_back(*newPoint);
  return newPoint;
}

// used for random numbers
double r(int range) {
  int imax = std::numeric_limits<int>::max();
  double whole = rand() % range;
  // printf("rand = %d\n", rand());
  double fraction = (double) (rand() % 10000) / 10000;
  // printf("fraction %lf\n", fraction);
  double retval = whole + fraction;
  if((int)retval % 2) // allow for negative
    retval = -retval;
  // printf("retval is %lf\n", retval);

  return retval;// random double
}

BoundingBox myboundary = BoundingBoxMake( -500,  -500,  500,  500);
int capacity = 1; // max number of bodies allowed inside a node
QuadTree qtree = QuadTreeMake(myboundary, capacity);
// QuadTreeNode* qtree =

void buildTree(std::vector<point> &p) {
  for(unsigned i = 0; i < p.size(); i++) {
    printf("inserting p with [%lf][%lf]\n", p[i].x, p[i].y);
    QuadTreeNodeDataMake(p[i].x, p[i].y, &p[i]);
  }
}

int main(int argc, char* argv[]){
  num_planets = (argc > 1) ? atoi(argv[1]): 120;
  int num_iterations = (argc > 2) ? atoi(argv[2]): 1000;

  /* create bodies */
  vector<point> bodies;

  G = 1.0; // This greatly increase the gravity



  // args are in form: (xp, yp, vx, vy, fx, fy, m, &bodies)

  srand( 1234567 );//time(NULL) ); // set seed for random

  // create the central and heavier body
  createBody(1, 1, 0, 0, r(10), r(10), abs(r(20)) + 10, bodies); // slightly of middle
  for(int i = 1; i < num_planets; i++) {
    createBody(r(100), r(100), r(4), r(4), r(10), r(10), abs(r(4)) + 1, bodies);
  }



  // QuadTreeNodeDataMake


  double start_time, end_time; /* start and end times */
  start_time = read_timer();

  srand( time(NULL) ); // set seed for random

  /* uses TIMESTEP for making time discrete */
  for(int i = 0; i < num_iterations; i++){
    /* build tree */
    buildTree(bodies);

    /* calculateForces */
    calcForces(bodies);

    /* move bodies */
    moveBodies(bodies);
  }
  end_time = read_timer();

  printf("%g\n", end_time - start_time);

  return 0;
}
