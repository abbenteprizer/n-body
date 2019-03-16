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

struct node {
  int has_particles; // 1 if it contains particle
  int has_innerpoint;
  point innerpoint;
  int level;
  double xlow;
  double xhigh;
  double ylow;
  double yhigh;

/* all the children */
  node *nw;
  node *ne;
  node *sw;
  node *se;
};

node root; // root node

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

void createBody(double xp, double yp, double vx, double vy, double fx, double fy, double m, std::vector<point> &bodies) {
  point* newPoint = new point;
  newPoint->x = xp;
  newPoint->y = yp;
  newPoint->vx = vx;
  newPoint->vy = vy;
  newPoint->fx = fx;
  newPoint->fy = fy;
  newPoint->m = m;
  bodies.push_back(*newPoint);
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
void setBound(node *newNode, double xl, double xh, double yl, double yh){
  newNode->xlow = xl;
  newNode->xhigh = xh;
  newNode->ylow = yl;
  newNode->yhigh = yh;
}
// xlow ylow is middle node division
void insertNode(point p, node *parentNode, int level, double xlow, double ylow) {
  // printf("inserting   point [%lf][%lf]\n", p.x, p.y);
  double lb = (BOUNDARY / (std::pow(2,level))) * 2; // local boundary
  if(parentNode->has_particles > 1) {
    // need to put into child, lets determine which
    if(p.x > xlow + lb) {
      // were in eastern blocks
      if(p.y > ylow + lb) {
        printf("insert into ne\n");
        insertNode(p, parentNode->ne, level + 1, xlow + lb, ylow + lb); // northeast
      } else {
        printf("insert into se\n");
        insertNode(p, parentNode->se, level + 1, xlow + lb, ylow - lb); // southeast
      }
    } else {
      // were in western blocks
      if(p.y > ylow + lb) {
        printf("insert into nw\n");
        insertNode(p, parentNode->nw, level + 1, xlow - lb, ylow + lb); // northwest
      } else {
        printf("insert into sw\n");
        insertNode(p, parentNode->sw, level + 1, xlow - lb, ylow - lb); // southwest
      }
    }

  } else if (parentNode->has_particles == 1){
    printf("splitting node into 4\n");
    parentNode->has_innerpoint = 0;
    parentNode->has_particles = parentNode->has_particles + 1;
    // printf("we have this %d many particle\n", parentNode->has_particles);

    parentNode->nw = (node*) malloc(sizeof(struct node)); // allocate all children
    parentNode->ne = (node*) malloc(sizeof(struct node));
    parentNode->sw = (node*) malloc(sizeof(struct node));
    parentNode->se = (node*) malloc(sizeof(struct node));

    setBound(parentNode->nw, xlow - lb, xlow, ylow, ylow + lb);
    setBound(parentNode->ne, xlow, xlow + lb, ylow, ylow + lb);
    setBound(parentNode->sw, xlow - lb, xlow, ylow - lb, ylow);
    setBound(parentNode->se, xlow, xlow + lb, ylow - lb, ylow);

    parentNode->nw->level = level + 1; // increase level for child
    parentNode->ne->level = level + 1;
    parentNode->sw->level = level + 1;
    parentNode->se->level = level + 1;


    /* insert innerpoint in subnode */
    if(parentNode->innerpoint.x > xlow + lb) {
      // were in eastern blocks
      if(parentNode->innerpoint.y > ylow + lb) {
        printf("insert into ne\n");
        insertNode(parentNode->innerpoint, parentNode->ne, level + 1, xlow + lb, ylow + lb); // northeast
      } else {
        printf("insert into se\n");
        insertNode(parentNode->innerpoint, parentNode->se, level + 1, xlow + lb, ylow - lb); // southeast
      }
    } else {
      // were in western blocks
      if(parentNode->innerpoint.y > ylow + lb) {
        printf("insert into nw\n");
        insertNode(parentNode->innerpoint, parentNode->nw, level + 1, xlow - lb, ylow + lb); // northwest
      } else {
        printf("insert into sw\n");
        insertNode(parentNode->innerpoint, parentNode->sw, level + 1, xlow - lb, ylow - lb); // southwest
      }
    }

    /* insert p in subnode */
    if(p.x > xlow + lb) {
      // were in eastern blocks
      if(p.y > ylow + lb) {
        insertNode(p, parentNode->ne, level + 1, xlow + lb, ylow + lb); // northeast
      } else {
        insertNode(p, parentNode->se, level + 1, xlow + lb, ylow - lb); // southeast
      }
    } else {
      // were in western blocks
      if(p.y > ylow + lb) {
        insertNode(p, parentNode->nw, level + 1, xlow - lb, ylow + lb); // northwest
      } else {
        insertNode(p, parentNode->sw, level + 1, xlow - lb, ylow - lb); // southwest
      }
    }

  } else {
    // printf("parent has %d\n", parentNode->has_particles);
    parentNode->has_particles = parentNode->has_particles + 1;
    parentNode->innerpoint = p;
    parentNode->has_innerpoint = 1;
    printf("added point [%lf,%lf] on level %d %lf %lf\n", p.x, p.y, parentNode->level, parentNode->xlow,parentNode->xhigh);


  }
}

void buildTree(std::vector<point> &p) {
  for(unsigned i = 0; i < p.size(); i++) {
    printf("inserting p with [%lf][%lf]\n", p[i].x, p[i].y);
    insertNode(p[i], &root, root.level, -500, -500);

  }
}


int main(int argc, char* argv[]){
  num_planets = (argc > 1) ? atoi(argv[1]): 120;
  int num_iterations = (argc > 2) ? atoi(argv[2]): 1000;

  /* create bodies */
  vector<point> bodies;

  G = 1.0; // This greatly increase the gravity

  root.level = 1; // remember this
  root.has_particles = 0;
  root.has_innerpoint = 0;
  setBound(&root, -500, 500, -500, 500);

  // args are in form: (xp, yp, vx, vy, fx, fy, m, &bodies)

  srand( 1234567 );//time(NULL) ); // set seed for random

  // create the central and heavier body
  createBody(1, 1, 0, 0, r(10), r(10), abs(r(20)) + 10, bodies); // slightly of middle
  for(int i = 1; i < num_planets; i++) {
    createBody(r(100), r(100), r(4), r(4), r(10), r(10), abs(r(4)) + 1, bodies);
  }

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
