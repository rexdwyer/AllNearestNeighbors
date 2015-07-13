#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#define REALBIG 9e99  
class EventQueue;
using namespace std;

/****************************************************************
 * Class of points in Euclidean 2-space.
 ****************************************************************/

class Point {
public:
  double x,y;
  // distance squared.
  double distance2(Point * other); 
  // true iff a,b,c are in counterclockwise order.
  // Since everything in Point is public, this declaration is 
  // meaningless, but ccw does naturally fit with the Point class.
  friend bool ccw(Point *a, Point *b, Point *c);
};

/****************************************************************
 * Template for classes of splay trees.  R is the type of the record
 * to be stored in the splay tree.  K is the type of the key on which
 * the tree must be ordered.  R must implement a function that returns
 * a key of type K.
 ****************************************************************/


template <class R, class K>
class SplayTree {
public:
  SplayTree <R,K> *left, *right;
  R * rec;

  SplayTree() { left = right = NULL; rec = NULL; }

  static SplayTree<R,K>* splay(SplayTree<R,K> * t, K splayKey);

  K key() { return rec->key(); }
  void printkey() { rec->printkey(); }
  void showtree(int level);
};

/****************************************************************
 * Class of "sites".  Sites are points plus some labels.
 ****************************************************************/

class Site: public Point {
public:
  int index;
  bool active;  // Is site currently active on sweepline?
  double bestdist;  // How close is the closest neighbor so far?
  Site * neighbor;  // What is the nearest neighbor so far?

  Site(int index) {
    bestdist = REALBIG;
    this->index = index;
    this->neighbor = NULL;
    this->active = false;
  }

  double key() { return(y); } // The key for the splay tree.

  void printkey() { cout << y; } // For tracing.
};


/****************************************************************
 * The instantiated template for the nodes in the splay tree 
 * representing the state of the sweepline.
 ****************************************************************/
typedef SplayTree <Site,double> SweepTree;

/****************************************************************
 * The class representing the sweepline.  Includes a splay tree of 
 * sites, a queue of events to update the splay tree, and the current
 * x coordinate of the sweepline (which moves horizontally.)
 ****************************************************************/
class Sweepline {
public:
  SplayTree<Site,double> * root = NULL;
  EventQueue * eventQueue;
  double sweep_x;

  Sweepline() { sweep_x = -REALBIG;}

  void schedule_deletion_if_needed(SweepTree *p);
  bool now_inactive(Site *p, Site *belowp, Site *abovep);
};

/****************************************************************
 * Events.  There are two subclasses: Insertion and Deletion.
 * Event objects go in the event queue.  They have a "scheduled"
 * x-coordinate and a site to be inserted or removed.
 ****************************************************************/
class Event {
public:
  double event_x; // scheduled "time" of the event. 
  Site *site;     // point to be inserted or deleted 

  // The events are organized in a skew heap implementation of the
  // priority queue abstract data type.  These are the child links.
  Event *left, *right;

  Event() { left = right = NULL; }
  Event(double event_x, Site * site) { 
    left = right = NULL;
    this->event_x = event_x;
    this->site = site;
  }

  // Updates the sweepline.
  virtual void process(Sweepline * sweepline) {};
};

/****************************************************************
 * Deletion events.
 ****************************************************************/
class Deletion: public Event {
public:
  void process(Sweepline * sweepline);
  Deletion(double event_x, Site* site) : Event(event_x,site) {};
};

/****************************************************************
 * Insertion events.
 ****************************************************************/
class Insertion: public Event {
public:
  void process(Sweepline * sweepline);
  Insertion(double event_x, Site* site) : Event(event_x, site) {};
};

/****************************************************************
 * Event queue.  It is a priority queue implemented as a skew heap.
 ****************************************************************/
class EventQueue {
private:
  Event* queue = NULL;
  static Event * merge_queues(Event *q1, Event *q2);
public:
  Sweepline * sweepline;
  bool isEmpty() { return !queue; }
  void enqueue(Event *e); // add a new event to queue.
  Event * dequeue();  // retrieve and remove the next event.
};


/****************************************************************
 * Clear enough, right?
 ****************************************************************/
double Point::distance2(Point * other) {
    double xdiff, ydiff;
    xdiff = x - other->x;
    ydiff = y - other->y;
    return(xdiff*xdiff + ydiff*ydiff);
}

/****************************************************************
 * True iff a,b, and c are in counterclockwise order.
 * This is the sign of a determinant.  The magnitude is proportional
 * to are of triangle abc.
 ****************************************************************/
bool ccw(Point * a, Point *b, Point *c) {
  /* true iff A, B, C form a counterclockwise oriented triangle */
  double cx, cy;
  cx = c->x; cy = c->y;
  if ( ((a->x - cx) * (b->y - cy) - (b->x - cx) * (a->y -cy)) > 0.0 )
	return (true);
  else
	return(false);
};

/****************************************************************
 * Returns the x-coordinate of the center of the circle passing through
 * a, b, and c.  This is used to determine when a site will be 
 * deleted from the sweepline.
 ****************************************************************/
double center_x(Point *a, Point *b, Point *c) {
  /* returns the x-coordinate of the 
     center of the circle passing through A, B & C. */
  double ax, ay, bx, by, cx, cy, numerator, denominator;
  ax = a->x; ay = a->y;
  bx = b->x; by = b->y;
  cx = c->x; cy = c->y;
  numerator = (cx-bx)*(cx-ax) + (cy-by)*(cy-ay);
  denominator = (bx-ax)*(cy-ay) - (by-ay)*(cx-ax);
  return(0.5*(ax + bx - (by-ay)*numerator/denominator));
  /*
	p = V2_sum(V2_times(0.5,V2_sum(a,b)),
	V2_times(V2_dot(V2_sub(c,b), V2_sub(c,a)) /
	(-2 * V2_cprod(V2_sub(b,a), V2_sub(c,a))),
	V2_cross(V2_sub(b,a))));
	*/
};

/****************************************************************
 * Used to display a splay tree for debugging.
 ****************************************************************/
template<class R, class K>
void SplayTree <R,K>::showtree(int level) {
  int i;
  if (left) left->showtree(level+1);
  if (level<10) for (i=0; i<level; i++) cout << "  ";
  else cout << "                    (" << level << ")";
  rec->printkey();
  cout << endl;
  if (right) right->showtree(level+1);
}

/****************************************************************
 * The meat of splay trees: the splaying operation.
 * It rotates the splay key to the root of the tree.
 ****************************************************************/
template<class R, class K>
SplayTree<R,K> * SplayTree <R,K>::splay(SplayTree<R,K> * t, K splayKey) {
  /* Simple top down splay, not requiring i to be in the tree t.  */
  SplayTree<R,K> *n, *l, *r, *q;
  if (t == NULL) return t;
  n = l = r = new SplayTree<R,K>();
  
  for (;;) {
    if (splayKey < t->key()) {
      if (t->left == NULL) break;
      if (splayKey < t->left->key()) {
	q = t->left;                           /* rotate right */
	t->left = q->right;
	q->right = t;
	t = q;
	if (t->left == NULL) break;
      }
      r->left = t;                               /* link right */
      r = t;
      t = t->left;
    } else if (splayKey > t->key()) {
      if (t->right == NULL) break;
      if (splayKey > t->right->key()) {
	q = t->right;                          /* rotate left */
	t->right = q->left;
	q->left = t;
	t = q;
	if (t->right == NULL) break;
      }
      l->right = t;                              /* link left */
      l = t;
      t = t->right;
    } else {
      break;
    }
  }
  l->right = t->left;                                /* assemble */
  r->left = t->right;
  t->left = n->right;
  t->right = n->left;

  return t;
}


/****************************************************************
 * When a site joins or leaves the sweepline, it may be necessary
 * to schedule the deletion of other sites.  This determines the
 * if and when of deletion using the counterclockwise test and 
 * the circle-center calculation.
 ****************************************************************/
void Sweepline::schedule_deletion_if_needed(SweepTree *s) {
  if (!s) return;
  Site * p = s->rec;
  //printf("Entering schedule_deletion_if_needed(%lf,%lf)\n", p->y, sweep_x);

  root = SweepTree::splay(root, p->key());
  root->left = SweepTree::splay(root->left, p->key());
  root->right = SweepTree::splay(root->right, p->key());
  
  if (!(root->left && root->right)) return; 

  double centerx = center_x(root->rec, root->left->rec, root->right->rec) ;
  if (ccw(root->rec, root->left->rec, root->right->rec)
      && (centerx > sweep_x)) {
    eventQueue-> enqueue(new Deletion(centerx, p));
  }
}

/****************************************************************
 * This determines if a site is made immediately inactive by
 * the arrival of a new site on the sweepline.
 ****************************************************************/
bool Sweepline::now_inactive(Site *p, Site *belowp, Site *abovep) {
  /* The points with biggest and smallest y are always active. */
  if (!belowp  || !abovep) return(false);
  if (ccw(p, belowp, abovep)
      && (center_x(p, belowp, abovep) < sweep_x))
    return(true);
  else  return(false);
  
}

/****************************************************************
 * This method implement deletion events. We must splay the site to
 * be deleted to the top of the tree, remove it, and then consider
 * whether its former neighbor need to have there deletion (re-)scheduled.
 ****************************************************************/
void Deletion::process(Sweepline * sweepline) {
  if (!this->site->active) return;  // Previously deleted.
  this->site->active = false;
  SweepTree * swroot = sweepline->root;
  swroot = SweepTree::splay(swroot, this->site->y);
  swroot->left = SweepTree::splay(swroot->left, this->site->y);
  swroot->right = SweepTree::splay(swroot->right, this->site->y);
  if (swroot->right) {
    swroot->right->left = swroot->left;
    swroot = swroot->right;
  }
  else swroot = swroot->left;
  sweepline->root = swroot;
  if (sweepline->root) {
    sweepline->schedule_deletion_if_needed(sweepline->root->right);
    sweepline->schedule_deletion_if_needed(sweepline->root->left);
  }
}

/****************************************************************
 * This method implements insertion events. We must "part" the splay tree
 * around the new site, add the new site to the tree, then repeatedly
 * delete sites above and below the new site if they have become inactive.
 ****************************************************************/
void Insertion::process(Sweepline * sweepline) {
  SweepTree *nu = new SweepTree;
  double bestdist, dist, newy;
  nu->rec  = site;
  newy = site->y;
  site->active = true;
  
  if (site->neighbor != NULL) bestdist = site->distance2(site->neighbor);
  else bestdist = REALBIG;
  
  // Insert new site into tree.
  SweepTree * swroot = sweepline->root;
  swroot = SweepTree::splay(swroot,newy);
  if (swroot == NULL) {
    nu->left = NULL;
    nu->right = NULL;
    swroot = nu;
  } 
  else if (swroot->rec->y == newy) {
    dist = site->distance2(swroot->rec);
    if (dist < bestdist) { 
      bestdist = dist; 
      site->neighbor = swroot->rec;
    }
    nu->left = swroot->left;
    nu->right = swroot->right;
    swroot->rec->active = false;
    swroot = nu;
  }
  else if (swroot->rec->y > site->y) {
    nu->left = swroot->left;
    nu->right = swroot;
    swroot->left = NULL;
    swroot = nu;
  }
  else  /*  (swroot->y < site->y) */ {
    nu->right = swroot->right;
    nu->left = swroot;
    swroot->right = NULL;
    swroot = nu;
  }
  sweepline->root = swroot;

  // Look at sites above the new site to see if they should be deleted.
  if (swroot->right) {
    swroot->right = SweepTree::splay(swroot->right, newy);
    dist = site->distance2(swroot->right->rec);
    if (dist < bestdist) { 
      bestdist = dist; 
      site->neighbor = swroot->right->rec;
    }
    swroot->right->right = SweepTree::splay(swroot->right->right, newy);
    while (swroot->right->right 
	   && sweepline->now_inactive(swroot->right->rec, 
				      swroot->rec, 
				      swroot->right->right->rec)
	   ) {
      swroot->right->rec->active = false;
      swroot->right = swroot->right->right;
      swroot->right->right = SweepTree::splay(swroot->right->right, newy);
      dist = site->distance2(swroot->right->rec);
      if (dist < bestdist) { 
	bestdist = dist; 
	site->neighbor = swroot->right->rec;
      }
    }
  }
  
  // Look at sites below the new site to see if they should be deleted.
  if (swroot->left) {
    swroot->left = SweepTree::splay(swroot->left, newy);
    dist = site->distance2(swroot->left->rec);
    if (dist < bestdist) { 
      bestdist = dist; 
      site->neighbor = swroot->left->rec;
    }
    swroot->left->left = SweepTree::splay(swroot->left->left, newy);
    while (swroot->left->left 
	   && sweepline->now_inactive(swroot->left->rec, 
				      swroot->left->left->rec, 
				      swroot->rec)
	   ) {
      swroot->left->rec->active = false;
      swroot->left = swroot->left->left;
      swroot->left->left = SweepTree::splay(swroot->left->left, newy);
      dist = site->distance2(swroot->left->rec);
      if (dist < bestdist) { 
	bestdist = dist; 
	site->neighbor = swroot->left->rec;
      }
    }
  }
  
  // Possibly schedule deletions of new site and its immediate 
  // neighbors on sweepline.
  sweepline->schedule_deletion_if_needed(swroot);
  sweepline->schedule_deletion_if_needed(swroot->left);
  sweepline->schedule_deletion_if_needed(swroot->right);
  
}

/****************************************************************
 * Adds a new event to the event queue.
 ****************************************************************/
void EventQueue::enqueue(Event * e) {
  e->left = NULL;
  e->right = NULL;
  queue = EventQueue::merge_queues(queue, e);
}

/****************************************************************
 * Removes the next event from the event queue.
 ****************************************************************/
Event * EventQueue::dequeue() {
  Event *p;
  p = queue;
  sweepline->sweep_x = p->event_x;
  queue = EventQueue::merge_queues(queue->left, queue->right);
  p->left = p->right = NULL;
  return p;
}

/****************************************************************
 * Implements the meat of skew heaps.  The basic operation is to
 * merge two trees.  For an insertion, a single-node tree is merged
 * with the existing tree.  For a deletion, the root is removed, and
 * its left and right subtrees are merged.
 * Skew heaps are somewhat like leftist heaps, but don't maintain
 * weights.  Left and right subtrees are swapped every time a node
 * is touched.
 ****************************************************************/
Event * EventQueue::merge_queues(Event *q1, Event *q2) {
  Event *t, result, *q;
  
  result.left = NULL;
  q = &result;
  
  while (q1 && q2) {
    if ((q1->event_x < q2->event_x) ||
	((q1->event_x==q2->event_x) && (q1->site->y < q2->site->y))) {
      q->left = q1;
      q = q1;
      q1 = q->right;
      q->right = q->left;
    }
    else {
      q->left = q2;
      q = q2;
      q2 = q->right;
      q->right = q->left;
    }
  }
  if (q1) q->left = q1; else q->left = q2;
  return(result.left); 
}

/****************************************************************
 * Here is it: the driver.
 ****************************************************************/
main() {
  
  // Read in points from stdin.  They are one per line, represented by
  // two floating point numbers.  Read them into a Site, and collect all
  // sites in a vector.

  vector<Site*> sites;
  std::string line;
  int i = 1;
  while (getline(cin,line)) {
    istringstream iss(line);
    Site * site = new Site(i++);
    iss >> site->x >> site->y;
    sites.push_back(site);
  }
  
  // We have to do two sweeps.  The first finds the nearest neighbor to the
  // left, and the second finds the nearest neighbor to the right.
  // We start by creating a sweepline and an event queue for the first sweep.

  Sweepline * sweepline1 = new Sweepline;
  EventQueue * eventQueue1 = new EventQueue;
  eventQueue1->sweepline = sweepline1;
  sweepline1->eventQueue = eventQueue1;

  // Now create an insertion event for each site and add it to the queue.
  // All insertion events are created at the beginning; deletion events
  // are created as the sweepline progresses.

  for (vector<Site*>::iterator it = sites.begin(); it != sites.end(); ++it) {
    Site *site = *it;
    eventQueue1->enqueue(new Insertion(site->x, site));
  }

  // Now fetch events one at a time and process.

  while (!eventQueue1->isEmpty()) {
    Event * e = eventQueue1->dequeue();
    e->process(eventQueue1->sweepline);
    delete e;
  }
  
  // Now create sweepline and event queue for the second sweep.

  Sweepline * sweepline2 = new Sweepline;
  EventQueue * eventQueue2 = new EventQueue;
  eventQueue2->sweepline = sweepline2;
  sweepline2->eventQueue = eventQueue2;

  // We fool the program into finding neighbors on the right by reflection
  // all the sites over the y-axis.  After reflecting each site, we queue
  // it up for the second sweep.

  for (vector<Site*>::iterator it = sites.begin(); it != sites.end(); ++it) {
    Site *site = *it;
    site->x = -site->x;
    eventQueue2->enqueue(new Insertion(site->x, site));
  }

  // Now work through the event queue for the second sweep.

  while (!eventQueue2->isEmpty()) {
    Event * e = eventQueue2->dequeue();
    e->process(eventQueue2->sweepline);
    delete e;
  }

  // Now un-reflect the sites.  We're not actually referring to the coordinates
  // after this, so this step is superfluous as things stand.

  for (vector<Site*>::iterator it = sites.begin(); it != sites.end(); ++it) {
    Site *site = *it;
    site->x = -site->x;
  }

  // Now march through the vector of sites and print out the index of the
  // nearest neighbor of each site.

  for (vector<Site*>::iterator it = sites.begin(); it != sites.end(); ++it) {
    Site *site = *it;
    if (site->neighbor)
      cout << site->index << " " << site->neighbor->index << endl;
    else  // This would be a bug...
      cout << site->index << " no neighbor" << endl;
  }
}
