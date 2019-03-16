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

class Quad
{
    // Hold details of the boundary of this node
    Point topLeft;
    Point botRight;

    // Contains details of node
    Node *n;

    // Children of this tree
    Quad *topLeftTree;
    Quad *topRightTree;
    Quad *botLeftTree;
    Quad *botRightTree;

public:
    Quad()
    {
        topLeft = Point(0, 0);
        botRight = Point(0, 0);
        n = NULL;
        topLeftTree  = NULL;
        topRightTree = NULL;
        botLeftTree  = NULL;
        botRightTree = NULL;
    }
    Quad(Point topL, Point botR)
    {
        n = NULL;
        topLeftTree  = NULL;
        topRightTree = NULL;
        botLeftTree  = NULL;
        botRightTree = NULL;
        topLeft = topL;
        botRight = botR;
    }
    void insert(Node*);
    Node* search(Point);
    bool inBoundary(Point);
};

// Insert a node into the quadtree
void Quad::insert(Node *node)
{
    if (node == NULL)
        return;

    // Current quad cannot contain it
    if (!inBoundary(node->pos))
        return;

    // We are at a quad of unit area
    // We cannot subdivide this quad further
    if (abs(topLeft.x - botRight.x) <= 1 &&
        abs(topLeft.y - botRight.y) <= 1)
    {
        if (n == NULL)
            n = node;
        return;
    }

    if ((topLeft.x + botRight.x) / 2 >= node->pos.x)
    {
        // Indicates topLeftTree
        if ((topLeft.y + botRight.y) / 2 >= node->pos.y)
        {
            if (topLeftTree == NULL)
                topLeftTree = new Quad(
                    Point(topLeft.x, topLeft.y),
                    Point((topLeft.x + botRight.x) / 2,
                        (topLeft.y + botRight.y) / 2));
            topLeftTree->insert(node);
        }

        // Indicates botLeftTree
        else
        {
            if (botLeftTree == NULL)
                botLeftTree = new Quad(
                    Point(topLeft.x,
                        (topLeft.y + botRight.y) / 2),
                    Point((topLeft.x + botRight.x) / 2,
                        botRight.y));
            botLeftTree->insert(node);
        }
    }
    else
    {
        // Indicates topRightTree
        if ((topLeft.y + botRight.y) / 2 >= node->pos.y)
        {
            if (topRightTree == NULL)
                topRightTree = new Quad(
                    Point((topLeft.x + botRight.x) / 2,
                        topLeft.y),
                    Point(botRight.x,
                        (topLeft.y + botRight.y) / 2));
            topRightTree->insert(node);
        }

        // Indicates botRightTree
        else
        {
            if (botRightTree == NULL)
                botRightTree = new Quad(
                    Point((topLeft.x + botRight.x) / 2,
                        (topLeft.y + botRight.y) / 2),
                    Point(botRight.x, botRight.y));
            botRightTree->insert(node);
        }
    }
}

// xmid ymid is middle node division
// void insertNode(point p, node *parentNode, int level, double xmid, double ymid) {
//   // printf("inserting   point [%lf][%lf]\n", p.x, p.y);
//   double lb = (BOUNDARY / (level)) / 2; // local boundary
//   if(parentNode->has_particles > 1) {
//     // need to put into child, lets determine which
//     if(p.x > (parentNode->xhigh - parentNode->xlow) / 2) {
//       // were in eastern blocks
//       if(p.y > (parentNode->yhigh - parentNode->ylow) / 2) {
//         insertNode(p, parentNode->ne, level + 1, xmid + lb, ymid + lb); // northeast
//       } else {
//         insertNode(p, parentNode->se, level + 1, xmid + lb, ymid - lb); // southeast
//       }
//     } else {
//       // were in western blocks
//       if(p.y > (parentNode->yhigh - parentNode->ylow) / 2) {
//         insertNode(p, parentNode->nw, level + 1, xmid - lb, ymid + lb); // northwest
//       } else {
//         insertNode(p, parentNode->sw, level + 1, xmid - lb, ymid - lb); // southwest
//       }
//     }
//
//   } else if (parentNode->has_particles == 1){
//     printf("splitting node into 4\n");
//     parentNode->has_innerpoint = 0;
//     parentNode->has_particles = parentNode->has_particles + 1;
//     // printf("we have this %d many particle\n", parentNode->has_particles);
//
//     parentNode->nw = (node*) malloc(sizeof(struct node)); // allocate all children
//     parentNode->ne = (node*) malloc(sizeof(struct node));
//     parentNode->sw = (node*) malloc(sizeof(struct node));
//     parentNode->se = (node*) malloc(sizeof(struct node));
//
//     setBound(parentNode->nw, xmid - lb, xmid, ymid, ymid + lb);
//     setBound(parentNode->ne, xmid, xmid + lb, ymid, ymid + lb);
//     setBound(parentNode->sw, xmid - lb, xmid, ymid - lb, ymid);
//     setBound(parentNode->se, xmid, xmid + lb, ymid - lb, ymid);
//
//     parentNode->nw->level = level + 1; // increase level for child
//     parentNode->ne->level = level + 1;
//     parentNode->sw->level = level + 1;
//     parentNode->se->level = level + 1;
//
//
//     /* insert innerpoint in subnode */
//     if(parentNode->innerpoint.x > (parentNode->xhigh - parentNode->xlow) / 2) {
//       printf("innerpoint.x was %lf, and xmid was %lf\n", parentNode->innerpoint.x, xmid);
//       // were in eastern blocks
//       if(parentNode->innerpoint.y > (parentNode->yhigh - parentNode->ylow) / 2) {
//         insertNode(parentNode->innerpoint, parentNode->ne, level + 1, xmid + lb, ymid + lb); // northeast
//       } else {
//         insertNode(parentNode->innerpoint, parentNode->se, level + 1, xmid + lb, ymid - lb); // southeast
//       }
//     } else {
//       // were in western blocks
//       if(parentNode->innerpoint.y > (parentNode->yhigh - parentNode->ylow) / 2) {
//         insertNode(parentNode->innerpoint, parentNode->nw, level + 1, xmid - lb, ymid + lb); // northwest
//       } else {
//         insertNode(parentNode->innerpoint, parentNode->sw, level + 1, xmid - lb, ymid - lb); // southwest
//       }
//     }
//
//     /* insert p in subnode */
//     if(p.x > (parentNode->xhigh - parentNode->xlow) / 2) {
//       // were in eastern blocks
//       if(p.y > (parentNode->yhigh - parentNode->ylow) / 2) {
//         insertNode(p, parentNode->ne, level + 1, xmid + lb, ymid + lb); // northeast
//       } else {
//         insertNode(p, parentNode->se, level + 1, xmid + lb, ymid - lb); // southeast
//       }
//     } else {
//       // were in western blocks
//       if(p.y > (parentNode->yhigh - parentNode->ylow) / 2) {
//         insertNode(p, parentNode->nw, level + 1, xmid - lb, ymid + lb); // northwest
//       } else {
//         insertNode(p, parentNode->sw, level + 1, xmid - lb, ymid - lb); // southwest
//       }
//     }
//
//   } else {
//     // printf("parent has %d\n", parentNode->has_particles);
//     parentNode->has_particles = parentNode->has_particles + 1;
//     parentNode->innerpoint = p;
//     parentNode->has_innerpoint = 1;
//     printf("added point on level %d, in region [%lf, %lf, %lf][%lf, %lf, %lf]\n ",
//             parentNode->level, xmid - lb, p.x, xmid + lb, ymid - lb, p.y, ymid + lb);
//
//
//   }
// }
//
// void buildTree(std::vector<point> &p) {
//   for(unsigned i = 0; i < p.size(); i++) {
//     printf("inserting p with [%lf][%lf]\n", p[i].x, p[i].y);
//     insertNode(p[i], &root, root.level, 0, 0);
//
//   }
// }


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



  Quad root;
  root.
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
