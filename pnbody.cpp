#include <iostream>
#include <cmath>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <limits>
#include <stdio.h>
#include <omp.h>
#include <stdio.h>
#include <time.h>

#include <sys/time.h>


// #define G 6.67e-11
#define DIMENSION 100
#define TIMESTEP 0.1

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
  double distance, magnitude;
  double xd, yd; // partial directions
  double e = 0.001; // margin to avoid zero division
  unsigned j;
  #pragma omp parallel for private(j, distance, magnitude, xd, yd)
  for(unsigned i = 0; i < p.size() - 1; i++) {
    for(j = i + 1; j < p.size(); j++) {
      distance = std::sqrt( std::pow((p[i].x - p[j].x), 2) +
  			                    std::pow((p[i].y - p[j].y), 2) );

      if(distance < e){ // avoids dividing by zero
        // distance = e;
        printf("too close the distance was %lf\n", distance);
      }
      magnitude = (G * p[i].m * p[j].m) / std::pow(distance, 2);

      xd = p[j].x - p[i].x; // direction with respect to x
      yd = p[j].y - p[i].y; // direction with respect to y
      // printf("did we get here?\n");
      /* update forces */
      p[i].fx = p[i].fx + magnitude * xd / distance;
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

    double color = (double) i / (double) num_planets;
    cout << p[i].x << " " << p[i].y << " "  << color  << endl;

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

int main(int argc, char* argv[]){
  int num_iterations = (argc > 1) ? atoi(argv[1]): 100;
  int num_threads = (argc > 2) ? atoi(argv[2]) : 4;
  num_planets = (argc > 3) ? atoi(argv[3]) : 50; // global
  omp_set_num_threads(num_threads);

  G = 1.0; // This greatly increase the gravity

  /* create bodies */
  vector<point> bodies;

  srand( time(NULL) ); // set seed for random

  // create the central and heavier body
  createBody(0, 0, 0, 0, r(10), r(10), abs(r(20)) + 100, bodies);
  for(int i = 1; i < num_planets; i++) {
    // args are in form: (xp, yp, vx, vy, fx, fy, m, &bodies)
    createBody(r(100), r(100), r(4), r(4), r(10), r(10), abs(r(4)) + 1, bodies);
  }

  double start_time, end_time; /* start and end times */
  start_time = omp_get_wtime();

  /* uses TIMESTEP for making time discrete */
  #pragma omp parallel
  {
    #pragma omp single
    {
      for(int i = 0; i < num_iterations; i++){
        /* calculateForces */
        calcForces(bodies);
        /* move bodies */
        moveBodies(bodies);

      }
    }
  }

  end_time = omp_get_wtime();
  printf("%d %g\n", num_threads, end_time - start_time);

  return 0;
}
