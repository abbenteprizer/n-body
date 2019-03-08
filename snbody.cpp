#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <omp.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdbool.h>

// #define G 6.67e-11
#define DIMENSION 100
#define TIMESTEP 0.1

double G = 6.67e-11;

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
  double distance, magnitude, direction;
  double xd, yd; // partial directions
  for(unsigned i = 0; i < p.size() - 1; i++) {
    for(unsigned j = i + 1; j < p.size(); j++) {
      distance = std::sqrt( std::pow(2, (p[i].x - p[j].x)) +
  			                    std::pow(2, (p[i].y - p[j].y)));

      magnitude = (G * p[i].m * p[j].m) / std::pow(2, distance);

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

    // cout << "node "<< i << " has position [" << p[i].x << "][" << p[i].y << "]" << endl;
    cout << p[i].x << " " << p[i].y << " "  << i +0.5 << endl;

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

int main(int argc, char* argv[]){
  int num_iterations = (argc > 1) ? atoi(argv[1]): 100;

  /* create bodies */
  vector<point> bodies;
  // just for fun
  G = 0.1; // This greatly increase the gravity

  // args are in form: (xp, yp, vx, vy, fx, fy, m, &bodies)
  createBody(2, 1.5, 0.3, 0.4, -2,  2, 1, bodies);
  createBody(3, 1, .5, 0.2, -3, 0, 1, bodies);
  // createBody(15, 1, -2, 0.2, 3, 3, 3, bodies);

  double start_time, end_time; /* start and end times */
  start_time = read_timer();

  /* uses TIMESTEP for making time discrete */
  for(int i = 0; i < num_iterations; i++){
    // printf("iteration %d\n", i);
    /* calculateForces */
    // printf("calculating forces\n");
    calcForces(bodies);

    // printf("moving bodies\n");
    /* move bodies */
    moveBodies(bodies);

  }
  end_time = read_timer();

  printf("%g\n", end_time - start_time);

  return 0;
}
