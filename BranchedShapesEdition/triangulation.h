/**
 * \file triangulation.h
 * \brief class  to manage Delauny triangulations
 * \author Luis Alvarez \n \n
*/


#ifndef TRIANGULATION_H
#define TRIANGULATION_H
#include <iostream>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <stdio.h>
#include <cstring>
#include <string.h>
#include <time.h>

#define ami_PI 3.14159265358979323846
#define myrand(gamma) ( gamma==1?( (double) rand()/RAND_MAX ):pow((double) rand()/(RAND_MAX+1.),gamma) )




using namespace std;

/**
 * \class  point2d
 * \brief class  to store 2D points and basic methods
 * \author Luis Alvarez
 */
class point2d{
  public :
  float x /** point x coordinate */;
  float y /** point y coordinate */;
  point2d():x( (float) 0), y( (float) 0){};
  ~point2d(){};
  point2d(const float xx , const float yy){x=xx; y=yy;}
  point2d(const float scalar){x=y=scalar;}
  point2d & operator=(const point2d &p){ x=p.x; y=p.y; return *this;}
  point2d & operator=(const float scalar){ x=scalar; y=scalar; return *this;}
  point2d (const point2d &p){x=p.x; y=p.y;}
  point2d operator+(const point2d &p)const {return point2d(x+p.x,y+p.y);}
  point2d operator-(const point2d &p)const {return point2d(x-p.x,y-p.y);}
  point2d operator*(const float &a)const {return point2d(a*x,a*y);}
  point2d operator/(const float &a)const {return point2d(x/a,y/a);}
  float operator*(const point2d &p)const {return ( (float) x*p.x+y*p.y);}
  inline float norm(){float paso=x*x+y*y; return(paso>0.?sqrtf(paso):0.);}
  inline float norm2(){ return(x*x+y*y);}
  void normalize(){float nor=norm(); if(nor>0.){x/=nor; y/=nor;}}
  point2d normal(){float nor=norm(); if(nor>0.) return point2d(-y/nor,x/nor); return point2d(0.,0.); }
  void print(){ std::cout << "point2d : (" << x << "," << y << ")" << std::endl;}

};

float R(point2d &v0,point2d &v1,point2d &v); // return scalar product (v-v0)*(v1-v0)^T

/**
 * \class  affinity
 * \brief class  to manage affine transformations
 * \author Luis Alvarez
 */
class affinity{
  public :
  float A[2][2] /** matriz transform */;
  float cx /** center horizontal coordinate */;
  float cy /** center vertical coordinate */;

  /// BASIC CONSTRUCTORS
  affinity(){A[0][0]=1.; A[0][1]=0.; A[1][0]=0.; A[1][1]=1.; cx=0.; cy=0.;}
  affinity(float B[2][2],float bx,float by){A[0][0]=B[0][0]; A[0][1]=B[0][1]; A[1][0]=B[1][0]; A[1][1]=B[1][1]; cx=bx; cy=by;}
  affinity(float B00,float B01,float B10,float B11,float bx,float by){A[0][0]=B00; A[0][1]=B01; A[1][0]=B10; A[1][1]=B11; cx=bx; cy=by;}
  affinity(float angle,float zoom,float aspec_ratio);
  affinity(float angle1,float angle2,float zoom,float horizontal_scaling);
  /** constructor of A satisfying :  A*p2=p0  && A*(p3-p2) || (p1-p0) */
  affinity(point2d p0,point2d p1,point2d p2,point2d p3);
  /** constructor of A satisfying :  (A*p10=p00+d or A*p20=p00+d)  && A*(p11-p10) || (p01-p00) &&
  the triangle A*(p10,p11,p12) does not intersects the triangle (p00,p01,p02) */
  affinity(point2d p00,point2d p01,point2d p02,point2d p10,point2d p11,point2d p12,float d,int type=0);

  /// BASIC OPERATIONS
  inline point2d operator*(const point2d &p)const {return point2d(A[0][0]*p.x+A[0][1]*p.y+cx,A[1][0]*p.x+A[1][1]*p.y+cy);}
  inline affinity operator+(const point2d p)const {affinity B=*this; B.cx+=p.x; B.cy+=p.y; return B;}
  affinity operator*(const affinity &B)const {
    float D[2][2],dx,dy;
    D[0][0]=A[0][0]*B.A[0][0]+A[0][1]*B.A[1][0];
    D[0][1]=A[0][0]*B.A[0][1]+A[0][1]*B.A[1][1];;
    D[1][0]=A[1][0]*B.A[0][0]+A[1][1]*B.A[1][0];;
    D[1][1]=A[1][0]*B.A[0][1]+A[1][1]*B.A[1][1];;
    dx=cx+A[0][0]*B.cx+A[0][1]*B.cy;
    dy=cy+A[1][0]*B.cx+A[1][1]*B.cy;
    return affinity(D,dx,dy);
  }
};

/**
 * \class  triangulation
 * \brief class  to manage triangulations and basic methods
 * \author Luis Alvarez
 */
class triangle{ public: int v[3]; triangle(){}; triangle(int v0,int v1,int v2){v[0]=v0; v[1]=v1; v[2]=v2;}};
class triangulation{
  public:
  vector<point2d> p; /** array of points representing the vertex */
  vector<int> fr;  /** auxiliary index vector  =0 (the point is inside the triangulation) != the point is in the boundary*/
  vector<triangle> t; /** array of triangles given by 3 indexes to the point2d vector */
  vector<triangle> n; /** triangle neighborhood information. For each triangle segment (t[k],t[(k+1)%3]) represents the index of the neighbor triangle (when negative no neighor triangle exists)*/

  /// BASIC CONSTRUCTORS
  void clear(){p.clear(); fr.clear(); t.clear(); n.clear();}
  triangulation(){p.clear(); fr.clear(); t.clear(); n.clear(); }; /// basic constructor
  triangulation(char *name); /// constructor reading  triangulation from a file
  triangulation(vector<float> &I,int width,int height,float sigma,int radius,float harris_treshold_value,float MinEdgeCos,float MaxEdgeDistance); /// constructor from a silhouette
  triangulation(point2d pmin,point2d pmax); ///constructor of a rectangle using extrema points
  triangulation(point2d pmin,point2d pmax,float min_distance); ///constructor of a rectangle using extrema points and a minimun distance between points
  triangulation(float strip_thick,point2d pmin,point2d pmax); ///constructor of the band of a rectangle with a given strip thick
  triangulation(point2d center,float radius,int Npoints); ///constructor of a circle using Npoints
  triangulation(vector<point2d> &center_line,vector<float> &line_radius ); ///constructor of line structures
  triangulation(vector<point2d> &p); ///constructor from a polygon
  triangulation(point2d center,float y[],float Rx[] ); ///constructor of a decagon symmetric with respect y axis
  triangulation(triangulation &T0, point2d p0,float dt,float a0,float la,float ra,float l0,float ll,float lr,int Niter); /// line structure RANDOM GENERATION
  triangulation(triangulation &T0, vector<point2d> &pV,vector<float> thick_radius); /// constructor of line structures using splines

  /// BASIC METHODS
  int write(const char *name);  /// method to write a triangulation in disk
  void update_n(); /// method to update triangle neighborhood information n
  void erase_point(int i); /// method to erase a point from the triangulation
  int check(); /// method to check if a triangulation is well defined.
  void clean(); /// method to clean up triangulations
  float normalize(); /// method to normalize the point coordinates
  void put_horizontal(); /// TRIANGULATION ROTATION TO OBTAIN THAT THE FIRST AND LAST POINTS ARE HORIZONTAL ALIGNED
  int rebuild(); /// we rebuild triangulation from the collection of points and given triangles
  int rebuild(float max_segment_length); /// we rebuild triangulation from the collection of points and given triangles
  int rebuildOnlyContour(float minArea=1e-3); /// we rebuild triangulation from the collection of points and given triangles using only contour points
  triangulation operator+(const triangulation &Tr); /// addition (joint) of 2 triangulations
  triangulation operator-(triangulation &Tr); /// substraction of 2 triangulations
  triangulation operator*(triangulation &Tr); /// intersection of 2 triangulations
  triangulation operator*(float scale){triangulation T2=(*this); for(int k=0;k<(int) p.size();k++) T2.p[k]=p[k]*scale; return T2;}
  triangulation operator+(point2d d){triangulation T2=(*this); for(int k=0;k<(int) p.size();k++) T2.p[k]=p[k]+d; return T2;}
  triangulation operator-(point2d d){triangulation T2=(*this); for(int k=0;k<(int) p.size();k++) T2.p[k]=p[k]-d; return T2;}
  triangulation symmetry(float a,float b,float c); /// compute symmetry with respect to a line
  triangulation symmetry(point2d c,int Nangles); /// build an Nangles symmetry triangulation from the original one
  triangulation border(float thick); /// build the internal border of a triangulation
  vector<point2d> extrema();/// we compute the extrema (xmin,ymin), (xmax,ymax) of triangulation coordinates
  triangulation deformation(point2d center,float angle0,float par0,float angle1,float par1);
  triangulation deformation2(point2d center,float angle0,float par0,float angle1,float par1);
  triangulation deformation3(point2d center,float angle0,float par0,float angle1,float par1);
  point2d mass_center(){point2d c_=0.; for(int k_=0;k_<(int) p.size();k_++){c_=c_+p[k_];} return(c_/p.size());}
  int inclusion_test(point2d v);  /// test to check if a point is inside a triangulation
  bool inclusion_test(triangulation &T);  /// test to check if a triangulation is inside a triangulation
  bool intersection_test(triangulation &T);  /// test to check if a triangulation intersect another triangulation
  void operator()(affinity A){
    for(int k=0;k<(int) p.size();k++){ p[k]=A*p[k];}
  }
  triangulation operator * (affinity A){ triangulation T2=(*this); for(int k=0;k<(int) T2.p.size();k++){ T2.p[k]=A*T2.p[k];} return T2; }
  vector<triangulation> connected_components(); /// TRIANGULATION SEPARATION IN CONNECTED COMPONENTS
  vector<triangulation> DivisionByLine(float a,float b,float c); /// DIVISION OF A TRIANGULATION IN TWO CUTTING BY A LINE
  double area(); /// triangulation area
  triangulation division(float max_length); /// division of a triangulation such that no triangle edge segment is larger than max_length
  int selectBoundarySegment(int tmin,int tmax,int &tselect); /// select randomly a triangle with a boundary edge segment in a range of triangle index values
  float max_segment_length(int t0); /// return max length of a triangle edge segments.
  void boundary_segments(float min_segment_length,vector<int> &SelectedTriangles,vector<int> &SelectedSegments,vector<int> &SelectedPoints);
  void boundary_segments(float min_segment_length,vector<int> &SelectedTriangles,vector<int> &SelectedSegments,vector<point2d> &SelectedPoints);
  void triangles_containing_a_point(int t_index,int p_index,vector<int> &t_select);
  vector<point2d> point_generation(float max_distance,int type);
  void fix_triangle_orientation();
  int point_index(point2d q); /// return the index of a point. If the point is not in the list it is included.
  triangulation border_stripe_periodicity(); /// fill a border stripe
};

vector<point2d> point_generation_splines(vector<point2d> &p, point2d d0, point2d dN,float segment_length);

int Delauny(triangulation &T,point2d v, int ValorFrontera); // include a point in a triangulation
int Delauny(triangulation &Tr,float a,float b,float c, int ValorFrontera); // include a line in a triangulation
int Delauny(triangulation &Tr,vector<point2d> &p1, int ValorFrontera); // include a collection of segments in a triangulation
void HalfPlaneIntersection(triangulation &Tr,float a,float b,float c); // Intersection with the half plane ax+by+c<0
int Symmetry(triangulation &Tr,float a,float b,float c); // build a new triangulation by symmetry

//affinity estimation from triangle segments
affinity triangles2affinity(triangulation &T1,int tindex1,int index1,triangulation &T2,int tindex2,int index2,float distance,float angle);

//triangulation leaf(float dt,float angle,float dangle,float r,float drl,int Niter,float a0,float a1);
affinity connectT12T0(triangulation &T0,int indexTriangle0,int indexSegment0,triangulation &T1,int indexTriangle1,int indexSegment1,float d,int type=0,point2d p2d=point2d(-1000.,-1000.));

float distance_point2segment(point2d &v0,point2d &v1,point2d &v);

vector<float> ami_harris_value_estimation(vector <float> &I, int width,int height,float sigma);

/// RETURN A RANDOM VALUE IN INTEGER PRECISION BETWEEN 0 AND THE PARAMETER In_
int ami_randI(int In_);


/// CHECKING FOR THE INTERCEPTION OF SEGMENTS v0-v1 and v2-v3
bool check_cross(point2d &v0,point2d &v1,point2d &v2,point2d &v3);

/// CHECKING IF v0 and v1 are in different half-planes of the line v2-v3
bool different_half_planes(point2d &v0,point2d &v1,point2d &v2,point2d &v3);


/// FUNCTION TO COMPUTE THE AVERAGE OF TO ANGLES
double angle_average(double a1,double a2);
/// FUNCTION TO COMPUTE THE DIFFERENCE BETWEEN 2 ANGLES
double angle_dif(double a1,double a2);

/**
 * \class  NodeTree
 * \brief class  to store a tree of nodes to generate and abstract shape
 * \author Luis Alvarez
 */

class NodeTree{ ///AAA
  public :
  vector< point2d > n_; /// node coordinates
  vector< float > r_; /// radius of the circle of each node
  vector< float > s_; /// smoothing factor of the node connection
  vector< vector<int> > i_; /// index of the nodes connected with the current node
  vector< vector<double> > a_; /// angles of the neigbor nodes with respect to the central node
  int shape_type=-1;
  int random_seed=-1;
  float average_node_distance=-1;
  float average_radius=-1;
  float smoothing_factor=-1;
  int branched_finished_in_peak=-1;
  int Xperiods=-1;
  int Yperiods=-1;
  float joint_extrema_probability=-1;
  int MaxNeigbors=-1;
  bool allow_branche_intersection=false;
  float node_radius_reduction_factor=-1;
  int MaxAngles=-1;
  unsigned int MaxNodes=0;
  bool round_termination=true;

  NodeTree(){};
  ~NodeTree(){};

  /// *************************************************************
  ///   CLEAR
  /// *************************************************************
  void clear(){ n_.clear(); r_.clear(); s_.clear(); i_.clear(); a_.clear(); };

  /// *************************************************************
  ///   EXTREMA COMPUTATION
  /// *************************************************************
  void extrema(double &xmin,double &ymin,double &xmax,double &ymax){
    xmin=xmax=n_[0].x;
    ymin=ymax=n_[0].y;
    for(int k=0;k<(int) n_.size();k++){
      if( n_[k].x-r_[k] < xmin) xmin=n_[k].x-r_[k];
      if( n_[k].x+r_[k] > xmax) xmax=n_[k].x+r_[k];
      if( n_[k].y-r_[k] < ymin) ymin=n_[k].y-r_[k];
      if( n_[k].y+r_[k] > ymax) ymax=n_[k].y+r_[k];
    }
  }

  /// *************************************************************
  /// CHECKING IF 2 NODES ARE CONNECTED
  /// *************************************************************
  bool node_connected_checking(int i,int j){
    if(i<0 || i>=(int) i_.size()) return false;
    for(int k=0;k<(int) i_[i].size();k++){
      if(i_[i][k]==j) return true;
    }
    return false;
  }

  /// **************************************************************
  /// DISTANCE FROM A POINT TO TreeNode SEGMENTS
  /// IT RETURNS THE DISTANCE TO THE CLOSEST SEGMENT JOINING 2 CIRCLES
  /// AS PARAMETERS IT ALSO RETURN THE INDEXES OF THE EXTREMA CIRCLES
  /// **************************************************************
  double segment_distance(
  double x,double y,  /// coordinates of point
  int &i,int &j) /// indexes of circles extrema of the segment with lowest distance
  {
    double error=1e50;
    i=j=-1;
    for(int k=0;k<(int) n_.size();k++){
      for(int n=0;n<(int) i_[k].size();n++){
        if(k<i_[k][n]) continue;
        double x1=n_[k].x;
        double y1=n_[k].y;
        double x2=n_[i_[k][n]].x;
        double y2=n_[i_[k][n]].y;
        if( (x2-x1)*(x-x1)+(y2-y1)*(y-y1) < 0  ) continue;
        if( (x1-x2)*(x-x2)+(y1-y2)*(y-y2) < 0  ) continue;
        double norm=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
        double d=fabs( -(y2-y1)*(x-x1)+(x2-x1)*(y-y1) )/norm;
        if(d<error){
          error=d;
          i=k;
          j=i_[k][n];
        }
      }
    }
    return error;
  }


  ///***********************************************************************
  ///                 BASIC PRINT
  ///***********************************************************************
  void print(){
    printf("\nSHAPE PARAMETERS\n");
    for(int k=0;k<(int) n_.size();k++){
      printf("n_[%d]=(%lf,%lf)\n",k,n_[k].x,n_[k].y);
    }
    for(int k=0;k<(int) r_.size();k++){
      printf("r_[%d]=%lf\n",k,r_[k]);
    }
    for(int k=0;k<(int) s_.size();k++){
      printf("s_[%d]=%lf\n",k,s_[k]);
    }
    for(int k=0;k<(int) i_.size();k++){
      printf("i_[%d]= ",k);
      for(int m=0;m<(int) i_[k].size();m++) printf("%d ",i_[k][m]);
      printf("\n");
    }

    for(int k=0;k<(int) a_.size();k++){
      printf("a_[%d]= ",k);
      for(int m=0;m<(int) a_[k].size();m++) printf("%lf ",180.*a_[k][m]/M_PI);
      printf("\n");
    }

    printf("shape_type = %d\n",shape_type);
    printf("random_seed = %d\n",random_seed);
    printf("average_node_distance = %f\n",average_node_distance);
    printf("average_radius = %f\n",average_radius);
    printf("smoothing_factor = %f\n",smoothing_factor);
    printf("branched_finished_in_peak = %d\n",branched_finished_in_peak);
    printf("Xperiods = %d\n",Xperiods);
    printf("Yperiods = %d\n",Yperiods);
    printf("joint_extrema_probability = %f\n",joint_extrema_probability);
    printf("MaxNeigbors = %d\n",MaxNeigbors);
    printf("allow_branche_intersection = %d\n",(int) allow_branche_intersection);
    printf("node_radius_reduction_factor = %f\n",node_radius_reduction_factor);
    printf("MaxAngles = %d\n",MaxAngles);
    printf("MaxNodes = %d\n",(int) MaxNodes);
    printf("round_termination = %d\n\n",(int) round_termination);

  };

  ///***********************************************************************
  ///                 BASIC PRINT
  ///***********************************************************************
  void print_params(){
    printf("\nSHAPE PARAMETERS\n");

    printf("shape_type = %d\n",shape_type);
    printf("random_seed = %d\n",random_seed);
    printf("average_node_distance = %f\n",average_node_distance);
    printf("average_radius = %f\n",average_radius);
    printf("smoothing_factor = %f\n",smoothing_factor);
    printf("branched_finished_in_peak = %d\n",branched_finished_in_peak);
    printf("Xperiods = %d\n",Xperiods);
    printf("Yperiods = %d\n",Yperiods);
    printf("joint_extrema_probability = %f\n",joint_extrema_probability);
    printf("MaxNeigbors = %d\n",MaxNeigbors);
    printf("allow_branche_intersection = %d\n",(int) allow_branche_intersection);
    printf("node_radius_reduction_factor = %f\n",node_radius_reduction_factor);
    printf("MaxAngles = %d\n",MaxAngles);
    printf("MaxNodes = %d\n",(int) MaxNodes);
    printf("round_termination = %d\n\n",(int) round_termination);

  };

  int read(char filename[400]){

    /// OPEN THE FILE
    FILE *f;
    f = fopen (filename, "r");
    if (f == NULL){
      printf("Error opening file %s\n",filename);
      return -1;
    }

    /// READING NodeTree STRUCTURE
    int N;
    double x,y,r,s;
    if (feof(f)) return -1;
    fscanf(f,"%d\n",&N);
    //printf("N=%d\n",N);
    /// UPDATING VECTOR SIZE
    n_.resize(N);
    r_.resize(N);
    s_.resize(N);
    i_.resize(N);

    for(int k=0;k<(int) n_.size();k++){
      if (feof(f)) return -1;
      fscanf(f,"%lf %lf\n",&x,&y);
      n_[k].x=x;
      n_[k].y=y;
      //printf("n_[k].x=%lf n_[k].y=%lf\n",n_[k].x,n_[k].y);
      //system("pause");
    }
    for(int k=0;k<(int) r_.size();k++){
      if (feof(f)) return -1;
      fscanf(f,"%lf\n",&r);
      r_[k]=r;
      //printf("r_[k]=%lf\n",r_[k]);
    }
    for(int k=0;k<(int) s_.size();k++){
      if (feof(f)) return -1;
      fscanf(f,"%lf\n",&s);
      s_[k]=s;
      //printf("s_[k]=%lf\n",s_[k]);
    }
    for(int k=0;k<(int) i_.size();k++){
      if (feof(f)) return -1;
      fscanf(f,"%d\n",&N);
      i_[k].resize(N);
      //printf("i_[k].size()=%d\n",(int) i_[k].size());
      for(int m=0;m<(int) i_[k].size();m++){
        if (feof(f)) return -1;
        fscanf(f,"%d\n",&N);
        i_[k][m]=N;
      }
    }


    if (feof(f)) {return -1;} fscanf(f,"%d\n",&shape_type);
    if (feof(f)) {return -1;} fscanf(f,"%d\n",&random_seed);
    if (feof(f)) {return -1;} fscanf(f,"%f\n",&average_node_distance);
    if (feof(f)) {return -1;} fscanf(f,"%f\n",&average_radius);
    if (feof(f)) {return -1;} fscanf(f,"%f\n",&smoothing_factor);
    if (feof(f)) {return -1;} fscanf(f,"%d\n",&branched_finished_in_peak);
    if (feof(f)) {return -1;} fscanf(f,"%d\n",&Xperiods);
    if (feof(f)) {return -1;} fscanf(f,"%d\n",&Yperiods);
    if (feof(f)) {return -1;} fscanf(f,"%f\n",&joint_extrema_probability);
    if (feof(f)) {return -1;} fscanf(f,"%d\n",&MaxNeigbors);
    int temp; if (feof(f)) {return -1;} fscanf(f,"%d\n",&temp);
    allow_branche_intersection =(bool) temp;
    if (feof(f)) {return -1;} fscanf(f,"%f\n",&node_radius_reduction_factor);
    if (feof(f)) {return -1;} fscanf(f,"%d\n",&MaxAngles);
    if (feof(f)) {return -1;} fscanf(f,"%d\n",&temp);
    MaxNodes=(unsigned int) temp;
    if (feof(f)) {return -1;} fscanf(f,"%d\n",&temp);
    round_termination =(bool) temp;

    fclose(f);
    return 0;
  };

  void write(char filename[400]){

    /// CHECKING VECTOR SIZE
    if(n_.size()!=r_.size()){
      printf("n_.size()!=r_.size()\n");
      return;
    }
    if(n_.size()!=s_.size()){
      printf("n_.size()!=s_.size()\n");
      return;
    }
    if(n_.size()!=i_.size()){
      printf("n_.size()!=i_.size()\n");
      return;
    }

    /// OPEN THE FILE
    FILE *f;
    f = fopen (filename, "w");
    if (f == NULL){
      printf("Error opening file %s\n",filename);
      return;
    }

    /// WRITING NodeTree STRUCTURE
    fprintf(f,"%d\n",(int) n_.size());
    for(int k=0;k<(int) n_.size();k++){
      fprintf(f,"%lf %lf\n",n_[k].x,n_[k].y);
    }
    for(int k=0;k<(int) r_.size();k++){
      fprintf(f,"%lf\n",r_[k]);
    }
    for(int k=0;k<(int) s_.size();k++){
      fprintf(f,"%lf\n",s_[k]);
    }
    for(int k=0;k<(int) i_.size();k++){
      fprintf(f,"%d\n",(int) i_[k].size());
      for(int m=0;m<(int) i_[k].size();m++) fprintf(f,"%d\n",i_[k][m]);
    }
    fprintf(f,"%d\n",shape_type);
    fprintf(f,"%d\n",random_seed);
    fprintf(f,"%lf\n",average_node_distance);
    fprintf(f,"%lf\n",average_radius);
    fprintf(f,"%lf\n",smoothing_factor);
    fprintf(f,"%d\n",branched_finished_in_peak);
    fprintf(f,"%d\n",Xperiods);
    fprintf(f,"%d\n",Yperiods);
    fprintf(f,"%lf\n",joint_extrema_probability);
    fprintf(f,"%d\n",MaxNeigbors);
    fprintf(f,"%d\n",(int) allow_branche_intersection);
    fprintf(f,"%lf\n",node_radius_reduction_factor);
    fprintf(f,"%d\n",MaxAngles);
    fprintf(f,"%d\n",(int) MaxNodes);
    fprintf(f,"%d\n",(int) round_termination);

    fclose(f);
  };

  void erase_point(int j){
    n_.erase(n_.begin()+j,n_.begin()+j+1);
    s_.erase(s_.begin()+j,s_.begin()+j+1);
    r_.erase(r_.begin()+j,r_.begin()+j+1);
    i_.erase(i_.begin()+j,i_.begin()+j+1);

    for(unsigned int m=0;m<i_.size();m++){
      for(unsigned int n=0;n<i_[m].size();n++){
        if(i_[m][n]==j){
          i_[m].erase(i_[m].begin()+n,i_[m].begin()+n+1);
          break;
        }
      }
    }

    for(unsigned int m=0;m<i_.size();m++){
      for(unsigned int n=0;n<i_[m].size();n++){
        if(i_[m][n]>j){ i_[m][n]--; }
      }
    }

  };

  void connect_nodes(int k,int j){
    for(int n=0;n<i_[k].size();n++){
      if(j==i_[k][n]) return;
    }
    i_[k].push_back(j);
    i_[j].push_back(k);
  }

  void disconnect_nodes(int k,int j){
    int pos_kj=-1;
    int pos_jk=-1;
    for(int n=0;n<i_[k].size();n++){
      if(j==i_[k][n]){
        pos_kj=n;
      }
    }
    for(int n=0;n<i_[j].size();n++){
      if(k==i_[j][n]){
        pos_jk=n;
      }
    }
    if(pos_kj==-1 || pos_jk==-1) return;

    i_[k].erase(i_[k].begin()+pos_kj,i_[k].begin()+pos_kj+1);
    i_[j].erase(i_[j].begin()+pos_jk,i_[j].begin()+pos_jk+1);
  }

  void join_nodes(int k,int j){
    vector<int> i2 = i_[j];
    n_.erase(n_.begin()+j,n_.begin()+j+1);
    s_.erase(s_.begin()+j,s_.begin()+j+1);
    r_.erase(r_.begin()+j,r_.begin()+j+1);
    i_.erase(i_.begin()+j,i_.begin()+j+1);

    for(unsigned int n=0;n<i_[k].size();n++){
      if(i_[k][n]==j) i_[k].erase(i_[k].begin()+n,i_[k].begin()+n+1);
       break;
    }

    for(unsigned int n=0;n<i2.size();n++){
      if(i2[n]==k){
        i2.erase(i2.begin()+n,i2.begin()+n+1);
        n--;
        continue;
      }
      if(i2[n]>j) i2[n]--;
    }

    for(unsigned int m=0;m<i_.size();m++){
      for(unsigned int n=0;n<i_[m].size();n++){
        if(i_[m][n]==j){ i_[m][n]=k; }
      }
    }

    for(unsigned int m=0;m<i_.size();m++){
      for(unsigned int n=0;n<i_[m].size();n++){
        if(i_[m][n]>j){ i_[m][n]--; }
      }
    }

    for(unsigned int n=0;n<i2.size();n++){
      unsigned int m;
      for(m=0;m<i_[k].size();m++){
        if(i_[k][m]==i2[n]) break;
      }
      if(m==i_[k].size()) i_[k].push_back(i2[n]);
    }

     for(int k=0;k<(int) i_.size();k++){
       for(int j=0;j<(int) i_[k].size();j++){
         for(int n=j+1;n<(int) i_[k].size();n++){
           if(i_[k][j]==i_[k][n]){
             i_[k].erase(i_[k].begin()+n,i_[k].begin()+n+1);
             n--;
           }
         }
       }
     }


    /// WE CHECK NODES ARE PROPERLY CONNECTED
    for(int k=0;k<(int) i_.size();k++){
      for(int j=0;j<(int) i_[k].size();j++){
        int m=i_[k][j],n;
        for(n=0;n<(int) i_[m].size();n++){
          if(i_[m][n]==k){
            break;
          }
        }
        if(n==(int) i_[m].size()){
           printf("node %d is neighbor of %d but it is not true the other way around\n",k,m);
           //system("pause");
           print();
           return;
           //system("pause");
        }

      }
    }


  };



  ///*********************************************************************************************
  ///   CHECK NodeTree CONFIGURATION. IT RETURNS true IF EVERYTHING IS FINE AND false OTHERWISE
  ///*********************************************************************************************
  bool checking(){
    if(i_.size()==0){
      printf("i_.size()==0\n");
      return false;
    }
    if(n_.size()!=i_.size()){
      printf("n_.size()!=i_.size()\n");
      return false;
    }
    if(r_.size()!=n_.size()){
      printf("r_.size()!=n_.size()\n");
      return false;
    }
    if(s_.size()!=n_.size()){
      printf("s_.size()!=n_.size()\n");
      return false;
    }

    if(n_.size()==1) return true;
    //printf("n_.size()=%d\n",(int) n_.size());

    /// WE CHECK THAT NODES ARE NOT TOO CLOSE. If 2 NODES ARE TOO CLOSE WE REMOVE ONE
    for(unsigned int k=0;k<n_.size();k++){
      for(unsigned int j=k+1;j<n_.size();j++){
        if( (n_[k]-n_[j]).norm2()<1. ){
          //printf("node n_[%d]=(%lf,%lf) is too close to n_[%d]=(%lf,%lf)\n",k,n_[k].x,n_[k].y,j,n_[j].x,n_[j].y);
          join_nodes((int) k, (int) j);
          j--;
        }
      }
    }

    /// WE CHECK NODE INDEX ARE CORRECT
    for(int k=0;k<(int) i_.size();k++){
//      if(i_[k].size()==0){
//        printf("i_[%d].size()==0\n",k);
//        return false;
//      }
      for(int m=0;m<(int) i_[k].size();m++){
        if(i_[k][m]>=(int) n_.size()){
          printf("i_[%d][%d]=%d n_.size()=%d\n",k,m,i_[k][m],(int) n_.size());
          return false;
        }
        if(i_[k][m]==k){
          printf("i_[%d][%d]=%d \n",k,m,i_[k][m]);
          return false;
        }
      }
    }
    /// WE CHECK NODES ARE PROPERLY CONNECTED
    for(int k=0;k<(int) i_.size();k++){
      for(int j=0;j<(int) i_[k].size();j++){
        int m=i_[k][j],n;
        for(n=0;n<(int) i_[m].size();n++){
          if(i_[m][n]==k){
            break;
          }
        }
        if(n==(int) i_[m].size()){
           printf("node %d is neighbor of %d but it is not true the other way around\n",k,m);
           //system("pause");
           print();
           //system("pause");
           return false;
        }

      }
    }

    /// WE CHECK ANGLES IN IT EACH NODE
    if(a_.size()!=n_.size()){
      a_=vector< vector<double> >(i_.size());
      node_angles_computation();
    }
    else{
      for(int k=0;k<(int) a_.size();k++){
        if(i_[k].size()!=a_[k].size()){
          node_angles_computation();
          break;
        }
      }
    }
    return true;
  };

  ///***********************************************************************
  /// COMPUTATION OF THE ANGLES IN EACH NODE CONFIGURATION
  /// NODE INDEX ARE ORDERED ACCORDING TO ANGLE VALUES
  ///***********************************************************************
  void node_angles_computation(){
   a_=vector< vector<double> >(i_.size());
   for(int k=0;k<(int) i_.size();k++){
     a_[k]=vector<double>(i_[k].size());
     for(int m=0;m<(int) i_[k].size();m++){
       point2d p = n_[i_[k][m]]-n_[k];
       a_[k][m] = atan2(p.y,p.x);
     }
     if(i_[k].size()<=1) continue;
     /// ORDER NODE INDEX BY ANGLES
     for(int i=0;i<(int) i_[k].size()-1;i++){
       for(int j=i+1;j<(int) i_[k].size();j++){
         if(a_[k][i]>a_[k][j]){
           double tempf=a_[k][i]; a_[k][i]=a_[k][j]; a_[k][j]=tempf;
           int tempi=i_[k][i]; i_[k][i]=i_[k][j]; i_[k][j]=tempi;
         }
       }
     }
   }
  };

  ///***********************************************************************
  /// COMPUTATION OF THE ANGLES IN ONE NODE CONFIGURATION
  /// NODE INDEX ARE ORDERED ACCORDING TO ANGLE VALUES
  ///***********************************************************************
  void node_angles_computation(int k){
   if((int) a_.size()<k+1) a_.resize(k+1);
   if(a_[k].size()!=i_[k].size()) a_[k]=vector<double>(i_[k].size());
   for(int m=0;m<(int) i_[k].size();m++){
     point2d p = n_[i_[k][m]]-n_[k];
     a_[k][m] = atan2(p.y,p.x);
   }
   /// ORDER NODE INDEX BY ANGLES
   for(int i=0;i<(int) i_[k].size()-1;i++){
     for(int j=i+1;j<(int) i_[k].size();j++){
        if(a_[k][i]>a_[k][j]){
           double tempf=a_[k][i]; a_[k][i]=a_[k][j]; a_[k][j]=tempf;
           int tempi=i_[k][i]; i_[k][i]=i_[k][j]; i_[k][j]=tempi;
        }
     }
   }
  };

  /// UPDATE INDEX AND ANGLE WITH A NEW NODE
  void update_index_angle(int k,int newIndex,double angle){
    if((int) a_.size()<k+1 || (int) i_.size()<k+1) return;
    if(a_[k].size()!=i_[k].size()) return;
    int j=0;
    for(;j<(int) i_[k].size();j++){
      if(a_[k][j]>angle) break;
    }
    i_[k].insert(i_[k].begin()+j,newIndex);
    a_[k].insert(a_[k].begin()+j,angle);
  }

  ///********************************************************
  /// TRIANGULATION GENERATION
  ///********************************************************
  triangulation triangulation_generation(float length_step,bool draw_circle=true){
    if(checking()==false) return triangulation();
    //printf("init triangulation_generation()\n");

    triangulation T;
    for(int k=0;k<(int) n_.size();k++){
      if(i_[k].size()==0){
        T=T+triangulation(n_[k],r_[k],60);
      }
      for(int m=0;m<(int) i_[k].size();m++){
        if(i_[k][m]<k) continue; /// we connect node in increasing order (to avoid duplication)
        /// LOCATION OF THE CORRESPONDING INDEX OF THE NODE
        int ikm=i_[k][m],imk=-1;
        for(int j=0;j<(int) i_[ikm].size();j++){
          if(i_[ikm][j]==k){
            imk=j;
            break;
          }
        }
        if(imk==-1){ /// the nodes are not connected properly
          printf("k=%d,m=%d,ikm=%d,imk=%d\n",k,m,ikm,imk);
          //system("pause");
          return triangulation();
        }

        ///COMPUTATION OF POINTS TO JOIN TO BUILD TRIANGLES
        point2d pk1,pk2,pm1,pm2;
        double ak1,ak2,am1,am2,zmin=0.,zmin0=0.;
        if(a_[k].size()==1){
          ak1 = a_[k][0]+M_PI*0.5;
          pk1 = n_[k]+point2d(cos(ak1),sin(ak1))*r_[k];
          ak2 = a_[k][0]-M_PI*0.5;
          pk2 = n_[k]+point2d(cos(ak2),sin(ak2))*r_[k];

        }
        else{
          int m1=(m+1)%a_[k].size();
          ak1=angle_average(a_[k][m],a_[k][m1]);
          //double z=sqrt(fabs(sin(0.5*angle_dif(a_[k][m],a_[k][m1]))));
          double z=fabs(sin(0.5*angle_dif(a_[k][m],a_[k][m1])));
          if(s_[k]<0.001 && z>zmin0) pk1 = n_[k]+point2d(cos(ak1),sin(ak1))*r_[k]/z;
          //if(true){ pk1 = n_[k]+point2d(cos(ak1),sin(ak1))*r_[k]/z; } //EEE
          else{
            z=sqrt(z);
            pk1 = n_[k]+point2d(cos(ak1),sin(ak1))*r_[k]/(z<zmin?zmin:z);
          }
          m1=m==0?a_[k].size()-1:m-1;
          ak2=angle_average(a_[k][m1],a_[k][m]);
          //z=sqrt(fabs(sin(0.5*angle_dif(a_[k][m1],a_[k][m]))));
          z=fabs(sin(0.5*angle_dif(a_[k][m1],a_[k][m])));
          if(s_[k]<0.001 && z>zmin0)  pk2 = n_[k]+point2d(cos(ak2),sin(ak2))*r_[k]/z;
          //if(true){ pk2 = n_[k]+point2d(cos(ak2),sin(ak2))*r_[k]/z; } //EEE
          else{
            z=sqrt(z);
            pk2 = n_[k]+point2d(cos(ak2),sin(ak2))*r_[k]/(z<zmin?zmin:z);
          }
        }
        if(a_[ikm].size()==1){
          am1 = a_[ikm][0]+M_PI*0.5;
          pm1 = n_[ikm]+point2d(cos(am1),sin(am1))*r_[ikm];
          am2 = a_[ikm][0]-M_PI*0.5;
          pm2 = n_[ikm]+point2d(cos(am2),sin(am2))*r_[ikm];
        }
        else{
          int m1=(imk+1)%a_[ikm].size();
          am1=angle_average(a_[ikm][imk],a_[ikm][m1]);
          //printf("%lf=angle_average(%lf,%lf)\n",angle*180./M_PI,a_[ikm][imk]*180./M_PI,a_[ikm][m1]*180./M_PI);
          //double z=sqrt(fabs(sin(0.5*angle_dif(a_[ikm][imk],a_[ikm][m1]))));
          double z=fabs(sin(0.5*angle_dif(a_[ikm][imk],a_[ikm][m1])));
          if(s_[ikm]<0.001 && z>zmin0) pm1 = n_[ikm]+point2d(cos(am1),sin(am1))*r_[ikm]/z;
          //if(true){ pm1 = n_[ikm]+point2d(cos(am1),sin(am1))*r_[ikm]/z; } //EEE
          else{
            z=sqrt(z);
            pm1 = n_[ikm]+point2d(cos(am1),sin(am1))*r_[ikm]/(z<zmin?zmin:z);
          }
          m1=imk==0?a_[ikm].size()-1:imk-1;
          am2=angle_average(a_[ikm][m1],a_[ikm][imk]);
          //printf("%lf=angle_average(%lf,%lf)\n",angle*180./M_PI,a_[ikm][imk]*180./M_PI,a_[ikm][m1]*180./M_PI);
          //z=sqrt(fabs(sin(0.5*angle_dif(a_[ikm][m1],a_[ikm][imk]))));
          z=fabs(sin(0.5*angle_dif(a_[ikm][m1],a_[ikm][imk])));
          if(s_[ikm]<0.001 && z>zmin0) pm2 = n_[ikm]+point2d(cos(am2),sin(am2))*r_[ikm]/z;
          //if(true){ pm2 = n_[ikm]+point2d(cos(am2),sin(am2))*r_[ikm]/z; } //EEE
          else{
            z=sqrt(z);
            pm2 = n_[ikm]+point2d(cos(am2),sin(am2))*r_[ikm]/(z<zmin?zmin:z);
          }
        }
        if( different_half_planes(pk1,pm1,n_[k],n_[ikm])==true ){
          if( different_half_planes(pk1,pm2,n_[k],n_[ikm])==true ){
            printf("problems in NodeTree::triangulation_generation()\n");
            printf("pk1=(%1.2lf,%1.2lf),pm1=(%1.2lf,%1.2lf),n_[%d]=(%1.2lf,%1.2lf),n_[%d]=(%1.2lf,%1.2lf)\n",
              pk1.x,pk1.y,pm1.x,pm1.y,k,n_[k].x,n_[k].y,ikm,n_[ikm].x,n_[ikm].y);
            //system("pause");
          }
          point2d tempp=pm1; pm1=pm2; pm2=tempp;
          float tempf=am1; am1=am2; am2=tempf;
        }

        /// WE INSERT THE POINTS IN THE TRIANGULATION (IF NECESSARY) AND GET THE POINT INDEX
        int ik1=T.point_index(pk1);
        int ik2=T.point_index(pk2);
        int im1=T.point_index(pm1);
        int im2=T.point_index(pm2);
        int ik=T.point_index(n_[k]);
        int im=T.point_index(n_[ikm]);
        /// WE INSERT BOUNDARY TRIANGLES (IF NECESSARY )
        if(i_[k].size()>2) T.t.push_back( triangle(ik,ik1,ik2) );
        if(i_[ikm].size()>2) T.t.push_back( triangle(im,im1,im2) );
        /// WE COMPUTE THE DERIVATIVE OF THE SPLINES IN THE EXTREMA POINTS OF THE SEGMENTS
        point2d dk1=point2d(-sin(ak1),cos(ak1))*s_[k];
        point2d dm1=point2d(-sin(am1),cos(am1))*s_[ikm];
        point2d dk2=point2d(-sin(ak2),cos(ak2))*s_[k];
        point2d dm2=point2d(-sin(am2),cos(am2))*s_[ikm];
        if( ( dk1.x*(pm1.x-pk1.x)+dk1.y*(pm1.y-pk1.y) ) < 0 ) { dk1.x=-dk1.x; dk1.y=-dk1.y;}
        if( ( dk2.x*(pm2.x-pk2.x)+dk2.y*(pm2.y-pk2.y) ) < 0 ) { dk2.x=-dk2.x; dk2.y=-dk2.y;}
        if( ( dm1.x*(pm1.x-pk1.x)+dm1.y*(pm1.y-pk1.y) ) < 0 ) { dm1.x=-dm1.x; dm1.y=-dm1.y;}
        if( ( dm2.x*(pm2.x-pk2.x)+dm2.y*(pm2.y-pk2.y) ) < 0 ) { dm2.x=-dm2.x; dm2.y=-dm2.y;}

        /// WE COMPUTE THE SPLINES TO CONNECT THE SEGMENTS
        vector<point2d> pV(2);
        pV[0]=pk1; pV[1]=pm1;
        vector<point2d> pV1=point_generation_splines(pV,dk1,dm1,length_step);
        pV[0]=pk2; pV[1]=pm2;
        vector<point2d> pV2=point_generation_splines(pV,dk2,dm2,length_step);

        /// WE INSERT THE POINTS AND TRIANGLES
        int Np = pV1.size()>pV2.size()?pV1.size():pV2.size();
        for(int k=0;k<Np-1;k++){
          int k1,k2,m1,m2;
          if(k<(int) pV1.size()) k1=T.point_index(pV1[k]);
          else k1=T.point_index(pV1[pV1.size()-1]);
          if(k+1<(int) pV1.size()) m1=T.point_index(pV1[k+1]);
          else m1=T.point_index(pV1[pV1.size()-1]);
          if(k<(int) pV2.size()) k2=T.point_index(pV2[k]);
          else k2=T.point_index(pV2[pV2.size()-1]);
          if(k+1<(int) pV2.size()) m2=T.point_index(pV2[k+1]);
          else m2=T.point_index(pV2[pV2.size()-1]);
//          printf("T.p[%d]=(%1.2lf,%1.2lf),T.p[%d]=(%1.2lf,%1.2lf),T.p[%d]=(%1.2lf,%1.2lf),T.p[%d]=(%1.2lf,%1.2lf)\n",
//                 k1,T.p[k1].x,T.p[k1].y,m1,T.p[m1].x,T.p[m1].y,k2,T.p[k2].x,T.p[k2].y,m2,T.p[m2].x,T.p[m2].y);
          //if( T.p[k1].x* T.p[k1].y<1 || T.p[k2].x* T.p[k2].y<1 || T.p[m1].x* T.p[m1].y<1 || T.p[m2].x* T.p[m2].y<1 ) system("pause");
          if(k1!=m1)  T.t.push_back( triangle(k1,m1,k2) );
          if(k2!=m2)  T.t.push_back( triangle(k2,m2,m1) );
        }

        /// WE INSERT A HALF CIRCLE AT THE EXTREMA NODES
        if(draw_circle==true){
          if(a_[k].size()==1){
            double a=a_[k][0]+M_PI*(1.5);
            int k1=T.point_index(n_[k]+point2d(cos(a),sin(a))*r_[k]);
            for(int j=0;j<29;j++){
              double a1=a_[k][0]+M_PI*(0.5+j/30.);
              double a2=a_[k][0]+M_PI*(0.5+(j+1.)/30.);
              int k2=T.point_index(n_[k]+point2d(cos(a1),sin(a1))*r_[k]);
              int k3=T.point_index(n_[k]+point2d(cos(a2),sin(a2))*r_[k]);
              T.t.push_back( triangle(k1,k2,k3) );
            }
          }
          if(a_[ikm].size()==1){
            double a=a_[ikm][0]+M_PI*(1.5);
            int k1=T.point_index(n_[ikm]+point2d(cos(a),sin(a))*r_[ikm]);
            for(int j=0;j<29;j++){
              double a1=a_[ikm][0]+M_PI*(0.5+j/30.);
              double a2=a_[ikm][0]+M_PI*(0.5+(j+1.)/30.);
              int k2=T.point_index(n_[ikm]+point2d(cos(a1),sin(a1))*r_[ikm]);
              int k3=T.point_index(n_[ikm]+point2d(cos(a2),sin(a2))*r_[ikm]);
              T.t.push_back( triangle(k1,k2,k3) );
            }
          }
        }
      }
    }
    T.update_n();

    //printf("T.p.size()=%d\n",(int) T.p.size());

    //printf("end triangulation_generation()\n");


    return T;

  };


  ///********************************************************
  /// nodeTree TRANSFORMATION
  ///********************************************************
  double transformation(double zoom=1.,double zx=0.,double zy=0.,double dx=0.,double dy=0.){
    static time_t timer=0;
    if(zoom!=1 || dx!=0 || dy!=0){
      for(int k=0;k<n_.size();k++){
        n_[k].x+=dx;
        n_[k].y+=dy;
      }

      if(zoom!=1){
        for(int k=0;k<n_.size();k++){
          n_[k].x=zx+zoom*(n_[k].x-zx);
          n_[k].y=zy+zoom*(n_[k].y-zy);
          r_[k]*=zoom;
        }
      }
      time_t timer2;
      time(&timer2);
      double seconds=difftime(timer2,timer);
      timer=timer2;
    }

    return 0;
  }

  ///********************************************************
  /// nodeTree check_equality
  ///********************************************************
  bool check_equality(NodeTree &nT){
    if(n_.size()!=nT.n_.size()) return false;

    for(int k=0;k<(int) n_.size();k++){
      if(n_[k].x!=nT.n_[k].x) return false;
      if(n_[k].y!=nT.n_[k].y) return false;
      if(s_[k]!=nT.s_[k]) return false;
      if(r_[k]!=nT.r_[k]) return false;
      if(i_[k].size()!=nT.i_[k].size()) return false;
      for(int n=0;n<(int) i_[k].size();n++){
        if(i_[k][n]!=nT.i_[k][n]) return(false);
      }
    }

    return(true);
  }

  ///********************************************************
  /// SVG FILE GENERATION
  ///********************************************************
  void svg_generation(char name[200],bool draw_nodes=false){

    node_angles_computation();

    vector < vector<point2d> > eP(n_.size());

    //bool draw_circle = (bool) branched_finished_in_peak;


    FILE *svgFile;
    svgFile = fopen (name , "w");
    if (svgFile == NULL) return;

    //double xmin,ymin,xmax,ymax;
    //extrema(xmin,ymin,xmax,ymax);

    // printf("xmin=%lf,ymin=%lf,xmax=%lf,ymax=%lf\n",xmin,ymin,xmax,ymax);


    //int height=1024*Yperiods;
    //int height=1024;

    //int dis=50;

    //int height = ymax+dis;


    fprintf(svgFile,"<?xml version=\"1.0\" standalone=\"no\"?> \n");
    fprintf(svgFile,"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"  \n");
    fprintf(svgFile,"\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\"> \n");

//    fprintf(svgFile,"\n<!-- SVG BRANCHED SHAPE CREATED WITH THE IPOL DEMO FACILITY --> \n");
//    fprintf(svgFile,"<!-- PARAMETERS USED IN THE DEMO INTERFACE: --> \n");
//    fprintf(svgFile,"   <!-- Branched shape type : %s --> \n",argv[1]);
//    fprintf(svgFile,"   <!-- Initial random seed generator : %s --> \n",argv[2]);
//    fprintf(svgFile,"   <!-- Average node distance : %s --> \n",argv[3]);
//    fprintf(svgFile,"   <!-- Average radius : %s --> \n",argv[4]);
//    fprintf(svgFile,"   <!-- Smoothing factor : %s --> \n",argv[5]);
//    fprintf(svgFile,"   <!-- Branched finished in peak : %s --> \n",argv[6]);
//    fprintf(svgFile,"   <!-- Number of periodic loops in horizontal direction : %s --> \n",argv[7]);
//    fprintf(svgFile,"   <!-- Number of periodic loops in vertical direction : %s --> \n",argv[8]);
//    fprintf(svgFile,"   <!-- Joint extrema probability : %s --> \n",argv[9]);
//    fprintf(svgFile,"   <!-- Allow branche intersection : %s --> \n",argv[10]);
//    fprintf(svgFile,"   <!-- Reduction factor of branche width trough the branches : %s --> \n",argv[11]);
//    fprintf(svgFile,"   <!-- Max number of angles in a node neighborhood  : %s --> \n",argv[12]);
//
//    if(argc==17){
//      fprintf(svgFile,"   <!-- Uploaded silhouette : %s --> \n",argv[11]);
//    }
//    else{
//      fprintf(svgFile,"   <!-- No silhouette uploaded --> \n");
//    }


///    fprintf(svgFile,"<!-- ************************************************* --> \n\n");

    //fprintf(svgFile,"<svg width=\"%dcm\" height=\"%dcm\" viewBox=\"%d %d %d %d\" \n",(int) (xmax-xmin+dis)/85,(int) (ymax-ymin+dis)/85,0,0,(int) xmax+dis,(int) ymax+dis);
    // fprintf(svgFile,"<svg width=\"%dcm\" height=\"%dcm\" viewBox=\"0 0 %d %d\" \n",12*Xperiods,12*Yperiods,1024*Xperiods,1024*Yperiods);
    fprintf(svgFile,"<svg width=\"1024px\" height=\"1024px\" viewBox=\"0 0 1024 1024\" \n");
    fprintf(svgFile,"xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\"> \n");

    //fprintf(svgFile,"<svg width=1024 height=1024 viewBox=\"0 0 1024 1024\"> \n");

//    if(n_.size()==1){
//      fprintf(svgFile,"<circle cx=\"%lf\" cy=\"%lf\" r=\"%lf\" stroke=\"none\" fill=\"red\" /> \n",n_[0].x,n_[0].y,r_[0]);
//      //fprintf(svgFile,"<circle cx=\"%lf\" cy=\"%lf\" r=\"%lf\" stroke=\"red\" stroke-width=\"1\" fill=\"none\" /> \n",n_[k].x,height-n_[k].y,r_[k]);
//      fprintf(svgFile,"</svg>");
//      fclose(svgFile);
//      return;
//    }


    triangulation T;
    for(int k=0;k<(int) n_.size();k++){
      if(i_[k].size()==0){
        fprintf(svgFile,"<circle cx=\"%lf\" cy=\"%lf\" r=\"%lf\" stroke=\"none\" fill=\"red\" /> \n",n_[k].x,n_[k].y,r_[k]);
        continue;
      }
      for(int m=0;m<(int) i_[k].size();m++){
        if(i_[k][m]<k) continue; /// we connect node in increasing order (to avoid duplication)
        #ifdef DEBUG_svg_generation
        printf("LOCATION OF THE CORRESPONDING INDEX OF THE NODE\n");
        #endif
        int ikm=i_[k][m],imk=-1;
        for(int j=0;j<(int) i_[ikm].size();j++){
          if(i_[ikm][j]==k){
            imk=j;
            break;
          }
        }
        if(imk==-1){ /// the nodes are not connected properly
          printf("k=%d,m=%d,ikm=%d,imk=%d\n",k,m,ikm,imk);
          //system("pause");
          return;
        }
        #ifdef DEBUG_svg_generation
        printf("COMPUTATION OF POINTS TO JOIN TO BUILD TRIANGLES\n");
        #endif
        point2d pk1,pk2,pm1,pm2;
        //double ak1,ak2,am1,am2,zmin=0.7,zmin0=0.3;
        double ak1,ak2,am1,am2,zmin=0.,zmin0=0.;
        if(a_[k].size()==1){
          ak1 = a_[k][0]+M_PI*0.5;
          pk1 = n_[k]+point2d(cos(ak1),sin(ak1))*r_[k];
          eP[k].push_back(pk1);
          ak2 = a_[k][0]-M_PI*0.5;
          pk2 = n_[k]+point2d(cos(ak2),sin(ak2))*r_[k];
          eP[k].push_back(pk2);
        }
        else{
          int m1=(m+1)%a_[k].size();
          ak1=angle_average(a_[k][m],a_[k][m1]);
          //double z=sqrt(fabs(sin(0.5*angle_dif(a_[k][m],a_[k][m1]))));
          double z=fabs(sin(0.5*angle_dif(a_[k][m],a_[k][m1])));
          if(s_[k]<0.001 && z>zmin0) pk1 = n_[k]+point2d(cos(ak1),sin(ak1))*r_[k]/z;
          else{
            z=sqrt(z);
            pk1 = n_[k]+point2d(cos(ak1),sin(ak1))*r_[k]/(z<zmin?zmin:z);
          }
          eP[k].push_back(pk1);
          m1=m==0?a_[k].size()-1:m-1;
          ak2=angle_average(a_[k][m1],a_[k][m]);
          //z=sqrt(fabs(sin(0.5*angle_dif(a_[k][m1],a_[k][m]))));
          z=fabs(sin(0.5*angle_dif(a_[k][m1],a_[k][m])));
          if(s_[k]<0.001 && z>zmin0)  pk2 = n_[k]+point2d(cos(ak2),sin(ak2))*r_[k]/z;
          else{
            z=sqrt(z);
            pk2 = n_[k]+point2d(cos(ak2),sin(ak2))*r_[k]/(z<zmin?zmin:z);
          }
          eP[k].push_back(pk2);
        }
        if(a_[ikm].size()==1){
          am1 = a_[ikm][0]+M_PI*0.5;
          pm1 = n_[ikm]+point2d(cos(am1),sin(am1))*r_[ikm];
          eP[ikm].push_back(pm1);
          am2 = a_[ikm][0]-M_PI*0.5;
          pm2 = n_[ikm]+point2d(cos(am2),sin(am2))*r_[ikm];
          eP[ikm].push_back(pm2);
        }
        else{
          int m1=(imk+1)%a_[ikm].size();
          am1=angle_average(a_[ikm][imk],a_[ikm][m1]);
          //printf("%lf=angle_average(%lf,%lf)\n",angle*180./M_PI,a_[ikm][imk]*180./M_PI,a_[ikm][m1]*180./M_PI);
          //double z=sqrt(fabs(sin(0.5*angle_dif(a_[ikm][imk],a_[ikm][m1]))));
          double z=fabs(sin(0.5*angle_dif(a_[ikm][imk],a_[ikm][m1])));
          if(s_[ikm]<0.001 && z>zmin0 ) pm1 = n_[ikm]+point2d(cos(am1),sin(am1))*r_[ikm]/z;
          else{
            z=sqrt(z);
            pm1 = n_[ikm]+point2d(cos(am1),sin(am1))*r_[ikm]/(z<zmin?zmin:z);
          }
          eP[ikm].push_back(pm1);
          m1=imk==0?a_[ikm].size()-1:imk-1;
          am2=angle_average(a_[ikm][m1],a_[ikm][imk]);
          //printf("%lf=angle_average(%lf,%lf)\n",angle*180./M_PI,a_[ikm][imk]*180./M_PI,a_[ikm][m1]*180./M_PI);
          //z=sqrt(fabs(sin(0.5*angle_dif(a_[ikm][m1],a_[ikm][imk]))));
          z=fabs(sin(0.5*angle_dif(a_[ikm][m1],a_[ikm][imk])));
          if(s_[ikm]<0.001 && z>zmin0 ) pm2 = n_[ikm]+point2d(cos(am2),sin(am2))*r_[ikm]/z;
          else{
            z=sqrt(z);
            pm2 = n_[ikm]+point2d(cos(am2),sin(am2))*r_[ikm]/(z<zmin?zmin:z);
          }
          eP[ikm].push_back(pm2);
        }
        if( different_half_planes(pk1,pm1,n_[k],n_[ikm])==true ){
          if( different_half_planes(pk1,pm2,n_[k],n_[ikm])==true ){
            printf("problems in NodeTree::triangulation_generation()\n");
            printf("pk1=(%1.2lf,%1.2lf),pm1=(%1.2lf,%1.2lf),n_[%d]=(%1.2lf,%1.2lf),n_[%d]=(%1.2lf,%1.2lf)\n",
              pk1.x,pk1.y,pm1.x,pm1.y,k,n_[k].x,n_[k].y,ikm,n_[ikm].x,n_[ikm].y);
            //system("pause");
          }
          point2d tempp=pm1; pm1=pm2; pm2=tempp;
          float tempf=am1; am1=am2; am2=tempf;
        }
        #ifdef DEBUG_svg_generation
        printf("WE INSERT THE POINTS IN THE TRIANGULATION (IF NECESSARY) AND GET THE POINT INDEX\n");
        #endif
        int ik1=T.point_index(pk1);
        int ik2=T.point_index(pk2);
        int im1=T.point_index(pm1);
        int im2=T.point_index(pm2);
        int ik=T.point_index(n_[k]);
        int im=T.point_index(n_[ikm]);
        /// WE INSERT BOUNDARY TRIANGLES (IF NECESSARY )
        if(i_[k].size()>2) T.t.push_back( triangle(ik,ik1,ik2) );
        if(i_[ikm].size()>2) T.t.push_back( triangle(im,im1,im2) );
        /// WE COMPUTE THE DERIVATIVE OF THE SPLINES IN THE EXTREMA POINTS OF THE SEGMENTS
        point2d dk1=point2d(-sin(ak1),cos(ak1))*s_[k];
        point2d dm1=point2d(-sin(am1),cos(am1))*s_[ikm];
        point2d dk2=point2d(-sin(ak2),cos(ak2))*s_[k];
        point2d dm2=point2d(-sin(am2),cos(am2))*s_[ikm];
        //if( ( dk1.x*(pm1.x-pk1.x)+dk1.y*(pm1.y-pk1.y) ) < 0 ) { dk1.x=-dk1.x; dk1.y=-dk1.y;}
        //if( ( dk2.x*(pm2.x-pk2.x)+dk2.y*(pm2.y-pk2.y) ) < 0 ) { dk2.x=-dk2.x; dk2.y=-dk2.y;}
        //if( ( dm1.x*(pm1.x-pk1.x)+dm1.y*(pm1.y-pk1.y) ) < 0 ) { dm1.x=-dm1.x; dm1.y=-dm1.y;}
        //if( ( dm2.x*(pm2.x-pk2.x)+dm2.y*(pm2.y-pk2.y) ) < 0 ) { dm2.x=-dm2.x; dm2.y=-dm2.y;}

        /// CAMBIO NUEVO
        if( ( dk1.x*(T.p[im].x-pk1.x)+dk1.y*(T.p[im].y-pk1.y) ) < 0 ) { dk1.x=-dk1.x; dk1.y=-dk1.y;}
        if( ( dk2.x*(T.p[im].x-pk2.x)+dk2.y*(T.p[im].y-pk2.y) ) < 0 ) { dk2.x=-dk2.x; dk2.y=-dk2.y;}
        if( ( dm1.x*(T.p[im].x-pk1.x)+dm1.y*(T.p[im].y-pk1.y) ) < 0 ) { dm1.x=-dm1.x; dm1.y=-dm1.y;}
        if( ( dm2.x*(T.p[im].x-pk2.x)+dm2.y*(T.p[im].y-pk2.y) ) < 0 ) { dm2.x=-dm2.x; dm2.y=-dm2.y;}

        ///
        fprintf(svgFile,"<path fill=\"red\" stroke=\"none\" d=\" \n");
        fprintf(svgFile,"M %lf,%lf\n",pk1.x,pk1.y);
        double norm=(pk1-pm1).norm();
        point2d pk11=pk1+dk1*(norm/3.);
        point2d pm11=pm1-dm1*(norm/3.);
        fprintf(svgFile,"C %lf,%lf %lf,%lf %lf,%lf\n",pk11.x,pk11.y,pm11.x,pm11.y,pm1.x,pm1.y);
        fprintf(svgFile,"L %lf,%lf\n",pm2.x,pm2.y);
        norm=(pk2-pm2).norm();
        point2d pk21=pk2+dk2*(norm/3.);
        point2d pm21=pm2-dm2*(norm/3.);
        fprintf(svgFile,"C %lf,%lf %lf,%lf %lf,%lf\n",pm21.x,pm21.y,pk21.x,pk21.y,pk2.x,pk2.y);
        fprintf(svgFile,"L %lf,%lf\n",pk1.x,pk1.y);

        fprintf(svgFile,"\">\n </path> \n");

        #ifdef DEBUG_svg_generation
        printf("WE INSERT BOUNDARY TRIANGLES (IF NECESSARY )\n");
        #endif
        if(i_[k].size()>2){
          T.t.push_back( triangle(ik,ik1,ik2) );
          fprintf(svgFile,"<path fill=\"red\" stroke=\"none\" d=\" \n");
          fprintf(svgFile,"M %lf,%lf\n",T.p[ik].x,T.p[ik].y);
          fprintf(svgFile,"L %lf,%lf\n",pk1.x,pk1.y);
          fprintf(svgFile,"L %lf,%lf\n",pk2.x,pk2.y);
          fprintf(svgFile,"\">\n </path> \n");
        }
        if(i_[ikm].size()>2){
          T.t.push_back( triangle(im,im1,im2) );
          fprintf(svgFile,"<path fill=\"red\" stroke=\"none\" d=\" \n");
          fprintf(svgFile,"M %lf,%lf\n",T.p[im].x,T.p[im].y);
          fprintf(svgFile,"L %lf,%lf\n",pm1.x,pm1.y);
          fprintf(svgFile,"L %lf,%lf\n",pm2.x,pm2.y);
          fprintf(svgFile,"\">\n </path> \n");
        }




//        /// WE COMPUTE THE SPLINES TO CONNECT THE SEGMENTS
//        vector<point2d> pV(2);
//        pV[0]=pk1; pV[1]=pm1;
//        vector<point2d> pV1=point_generation_splines(pV,dk1,dm1,length_step);
//        pV[0]=pk2; pV[1]=pm2;
//        vector<point2d> pV2=point_generation_splines(pV,dk2,dm2,length_step);
//
//        /// WE INSERT THE POINTS AND TRIANGLES
//        int Np = pV1.size()>pV2.size()?pV1.size():pV2.size();
//        for(int k=0;k<Np-1;k++){
//          int k1,k2,m1,m2;
//          if(k<(int) pV1.size()) k1=T.point_index(pV1[k]);
//          else k1=T.point_index(pV1[pV1.size()-1]);
//          if(k+1<(int) pV1.size()) m1=T.point_index(pV1[k+1]);
//          else m1=T.point_index(pV1[pV1.size()-1]);
//          if(k<(int) pV2.size()) k2=T.point_index(pV2[k]);
//          else k2=T.point_index(pV2[pV2.size()-1]);
//          if(k+1<(int) pV2.size()) m2=T.point_index(pV2[k+1]);
//          else m2=T.point_index(pV2[pV2.size()-1]);
////          printf("T.p[%d]=(%1.2lf,%1.2lf),T.p[%d]=(%1.2lf,%1.2lf),T.p[%d]=(%1.2lf,%1.2lf),T.p[%d]=(%1.2lf,%1.2lf)\n",
////                 k1,T.p[k1].x,T.p[k1].y,m1,T.p[m1].x,T.p[m1].y,k2,T.p[k2].x,T.p[k2].y,m2,T.p[m2].x,T.p[m2].y);
//          //if( T.p[k1].x* T.p[k1].y<1 || T.p[k2].x* T.p[k2].y<1 || T.p[m1].x* T.p[m1].y<1 || T.p[m2].x* T.p[m2].y<1 ) system("pause");
//          if(k1!=m1)  T.t.push_back( triangle(k1,m1,k2) );
//          if(k2!=m2)  T.t.push_back( triangle(k2,m2,m1) );
//        }
//
        #ifdef DEBUG_svg_generation
        printf("WE INSERT A HALF CIRCLE AT THE EXTREMA NODES\n");
        #endif

        if(round_termination==true){
          if(a_[k].size()==1){
            int k1=T.point_index(n_[k]);
            for(int j=0;j<30;j++){
              double a1=a_[k][0]+M_PI*(0.5+j/30.);
              double a2=a_[k][0]+M_PI*(0.5+(j+1.)/30.);
              int k2=T.point_index(n_[k]+point2d(cos(a1),sin(a1))*r_[k]);
              int k3=T.point_index(n_[k]+point2d(cos(a2),sin(a2))*r_[k]);
              T.t.push_back( triangle(k1,k2,k3) );
            }
            double a2=a_[k][0]+M_PI*(0.5);
            double a1=a_[k][0]+M_PI*(1.5);
            point2d p1=n_[k]+point2d(cos(a1),sin(a1))*r_[k];
            point2d p2=n_[k]+point2d(cos(a2),sin(a2))*r_[k];

            fprintf(svgFile,"<path fill=\"red\" stroke=\"none\" d=\" \n");
            fprintf(svgFile,"M %lf,%lf\n",p1.x,p1.y);
            fprintf(svgFile,"A %lf,%lf %d %d,%d %lf,%lf \n",r_[k],r_[k],0,1,0,p2.x,p2.y);
            fprintf(svgFile,"\">\n </path> \n");


          }
          if(a_[ikm].size()==1){
            int k1=T.point_index(n_[ikm]);
            for(int j=0;j<30;j++){
              double a1=a_[ikm][0]+M_PI*(0.5+j/30.);
              double a2=a_[ikm][0]+M_PI*(0.5+(j+1.)/30.);
              int k2=T.point_index(n_[ikm]+point2d(cos(a1),sin(a1))*r_[ikm]);
              int k3=T.point_index(n_[ikm]+point2d(cos(a2),sin(a2))*r_[ikm]);
              T.t.push_back( triangle(k1,k2,k3) );
            }
            double a2=a_[ikm][0]+M_PI*(0.5);
            double a1=a_[ikm][0]+M_PI*(1.5);
            point2d p1=n_[ikm]+point2d(cos(a1),sin(a1))*r_[ikm];
            point2d p2=n_[ikm]+point2d(cos(a2),sin(a2))*r_[ikm];

            fprintf(svgFile,"<path fill=\"red\" stroke=\"none\" d=\" \n");
            fprintf(svgFile,"M %lf,%lf\n",p1.x,p1.y);
            fprintf(svgFile,"A %lf,%lf %d %d,%d %lf,%lf \n",r_[ikm],r_[ikm],0,0,0,p2.x,p2.y);
            fprintf(svgFile,"\">\n </path> \n");

          }
        }
      }
    }



///------------------------------------------------------------------------------------------------------------------------------------------------

  if(draw_nodes==true){
    for(int k=0;k<(int) n_.size();k++){
      if(i_[k].size()==0){
        float r=r_[k]>=10?5:r_[k]/2.;
        fprintf(svgFile,"<circle cx=\"%lf\" cy=\"%lf\" r=\"%lf\" stroke=\"none\" fill=\"black\" /> \n",n_[k].x,n_[k].y,r);
        fprintf(svgFile,"<circle cx=\"%lf\" cy=\"%lf\" r=\"%lf\" stroke=\"black\" stroke-width=\"1\" fill=\"none\" class=\"draggable-circle\"/> \n",n_[k].x,n_[k].y,r_[k]);
        continue;
      }
      for(int m=0;m<(int) i_[k].size();m++){
        if(i_[k][m]<k) continue; /// we connect node in increasing order (to avoid duplication)
        fprintf(svgFile,"<line x1=\"%lf\" y1=\"%lf\" x2=\"%lf\" y2=\"%lf\" style=\"stroke:rgb(0,0,0);stroke-width:1\" />  \n",n_[k].x,n_[k].y,n_[i_[k][m]].x,n_[i_[k][m]].y);
      }
      float r=r_[k]>=10?5:r_[k]/2.;
      fprintf(svgFile,"<circle cx=\"%lf\" cy=\"%lf\" r=\"%lf\" stroke=\"none\" fill=\"black\" /> \n",n_[k].x,n_[k].y,r);
      fprintf(svgFile,"<circle cx=\"%lf\" cy=\"%lf\" r=\"%lf\" stroke=\"black\" stroke-width=\"1\" fill=\"none\" class=\"draggable-circle\"/> \n",n_[k].x,n_[k].y,r_[k]);
    }
    for(int k=0;k<(int) eP.size();k++){
      for(int m=0;m<(int) eP[k].size();m++){
        fprintf(svgFile,"<circle cx=\"%lf\" cy=\"%lf\" r=\"%lf\" stroke=\"none\" fill=\"green\" /> \n",eP[k][m].x,eP[k][m].y,3.);
        fprintf(svgFile,"<line pathLength=\"100\" stroke-dasharray=\"5,10,5,11,5,11,5,10,5,11,5,11,5\" x1=\"%lf\" y1=\"%lf\" x2=\"%lf\" y2=\"%lf\" style=\"stroke:rgb(0,255,0);stroke-width:1\" />  \n",n_[k].x,n_[k].y,eP[k][m].x,eP[k][m].y);
      }
    }
  }
//    T.update_n();
  fprintf(svgFile,"</svg>");

  fclose(svgFile);

  return;

};

  /// *************************************************************
  ///    INSERT A NEW NODE IN THE MIDDLE OF THE SEGMENT
  /// *************************************************************
  void insert_point(point2d q,float radius, float smooth_factor,int i,int j){
    int pos=n_.size();
    if(i<0 || i>=pos || j<0 || j>=pos || i==j) return;

    if(node_connected_checking(i,j)==true){
      disconnect_nodes(i,j);
    }

    n_.push_back(q);
    r_.push_back(radius);
    s_.push_back(smooth_factor);

    i_.push_back(vector<int>(2));
    i_[pos][0]=i;
    i_[pos][1]=j;

    i_[i].push_back(pos);
    i_[j].push_back(pos);


//    for(int k=0;k<(int) i_[i].size();k++){
//      if(i_[i][k]==j){
//        i_[i][k]=pos;
//        break;
//      }
//    }
//    for(int k=0;k<(int) i_[j].size();k++){
//      if(i_[j][k]==i){
//        i_[j][k]=pos;
//        break;
//      }
//    }
    //node_angles_computation();
  };

  /// *************************************************************
  /// INSERT A NEW NODE ASSOCIATED TO A GIVEN ONE
  /// *************************************************************
  void insert_point(point2d q,float radius,float smooth_factor,int i){
    int pos=n_.size();
    if(i<0 || i>=pos) return;

    n_.push_back(q);
    r_.push_back(radius);
    s_.push_back(smooth_factor);

    i_.push_back(vector<int>(1,i));
    i_[i].push_back(pos);

    //node_angles_computation();
  };

  /// *************************************************************
  /// INSERT A NEW ISOLATED NODE
  /// *************************************************************
  void insert_point(point2d q,float radius=-1,float smooth_factor=-1){

    if(radius<0){
      radius=0;
      for(unsigned int k=0;k<r_.size();k++) radius+=r_[k];
      radius/=r_.size();
    }

    if(smooth_factor<0){
      smooth_factor=0;
      for(unsigned int k=0;k<r_.size();k++) smooth_factor+=s_[k];
      smooth_factor/=s_.size();
    }

    n_.push_back(q);
    r_.push_back(radius);
    s_.push_back(smooth_factor);
    i_.push_back(vector<int>());

  };

  float average_length_segment(){
    float d=0;
    int N=0;
    for(unsigned int k=0;k< i_.size();k++){
      for(unsigned int j=0;j< i_[k].size();j++){
        d+=(n_[k]-n_[i_[k][j]]).norm();
        N++;
      }
    }
    return d/N;

  };

  /// ***********************************************************************************
  ///          FUNCTION TO COMPUTE A VECTOR WITH A RANGE OF ANGLES IN [0,2*M_PI]
  /// ***********************************************************************************
  vector<double> angle_range(unsigned int MaxNangles){
    vector<double> a(1,M_PI/2);  //0
//    a.push_back(3*M_PI/4); // 1
//    a.push_back(M_PI/4); // 2
//    a.push_back(5*M_PI/8); // 3
//    a.push_back(3*M_PI/8); // 7
//    a.push_back(7*M_PI/16); // 5
//    a.push_back(9*M_PI/16); // 6
//    a.push_back(7*M_PI/8); // 4
//    a.push_back(M_PI/8); // 8
//    a.push_back(5*M_PI/16); // 9
//    a.push_back(11*M_PI/16); // 10
//    a.push_back(3*M_PI/16); // 11
//    a.push_back(13*M_PI/16); // 12
//    a.push_back(M_PI/16); // 13
//    a.push_back(15*M_PI/16); // 14

    a.push_back(7*M_PI/16); // 5
    a.push_back(9*M_PI/16); // 6
    a.push_back(5*M_PI/8); // 3
    a.push_back(3*M_PI/8); // 7
    a.push_back(5*M_PI/16); // 9
    a.push_back(11*M_PI/16); // 10
    a.push_back(3*M_PI/4); // 1
    a.push_back(M_PI/4); // 2
    a.push_back(3*M_PI/16); // 11
    a.push_back(13*M_PI/16); // 12
    a.push_back(7*M_PI/8); // 4
    a.push_back(M_PI/8); // 8
    a.push_back(M_PI/16); // 13
    a.push_back(15*M_PI/16); // 14

    a.push_back(0.); // 15
    a.push_back(M_PI); // 16
    a.push_back(-M_PI/16); // 17
    a.push_back(-15*M_PI/16); // 18
    a.push_back(-2*M_PI/16); // 19
    a.push_back(-14*M_PI/16); // 20
    a.push_back(-3*M_PI/16); // 21
    a.push_back(-13*M_PI/16); // 22
    a.push_back(-4*M_PI/16); // 23
    a.push_back(-12*M_PI/16); // 24
    a.push_back(-5*M_PI/16); // 25
    a.push_back(-11*M_PI/16); // 26
    a.push_back(-6*M_PI/16); // 27
    a.push_back(-10*M_PI/16); // 28
    a.push_back(-7*M_PI/16); // 29
    a.push_back(-9*M_PI/16); // 30
    a.push_back(-8*M_PI/16); // 31

    if( MaxNangles<a.size() ){
      a.resize(MaxNangles);
    }

//    for(unsigned int k=0;k<a.size();k++){
//        printf("a[%d]=%lf\n",k,a[k]*180./M_PI);
//    }
//    system("pause");



    return a;

//    vector<double> a(1,0.);
//    a.push_back(M_PI);
//    double length=M_PI*0.5;
//    while( length>step){
//      for(double angle=length;angle<(2.*M_PI-1e-10);angle+=2*length){
//        a.push_back(angle);
//      }
//      length/=2;
//    }
////    for(unsigned int k=0;k<a.size();k++){
////        printf("a[%d]=%lf\n",k,a[k]*180./M_PI);
////    }
////    system("pause");
//    return a;
  }

    /// ***********************************************************************************
  ///          FUNCTION TO COMPUTE A VECTOR WITH A RANGE OF ANGLES IN [0,2*M_PI]
  /// ***********************************************************************************
  vector<double> angle_range_rand(unsigned int MaxNangles){
    vector<double> b=angle_range(100);  //0
    if(MaxNangles>=b.size()) return b;

    while(MaxNangles<b.size()){
      int k=rand()%b.size();
      b.erase(b.begin()+k,b.begin()+k+1);
    }

    for(unsigned int k=0;k<b.size();k++){
      unsigned int k2 = rand()%b.size();
      if(k==k2) continue;
      float temp=b[k];
      b[k]=b[k2];
      b[k2]=temp;
    }

    return b;

  }

  /// **********************************************************************
  ///       TEST IF A POINT IS INSIDE THE INITIAL SYMMETRIC REGION
  /// **********************************************************************
  bool InsideSymmetricRegion(float x,float y,float cx,float cy,int SymmetryType){
    //printf("x=%lf y=%lf SymmetryType=%d\n",x,y,SymmetryType);
    if(SymmetryType==0) return true;
    if(SymmetryType==1 || SymmetryType==4){
      if(x<cx) return false;
      //system("pause");
      return true;
    }
    if(SymmetryType==2 || SymmetryType==5){
      if(y<cy || x<cx) return false;
      //system("pause");
      return true;
    }
    if(SymmetryType==3 || SymmetryType==6){
      if(y<cy || x<cx) return false;
      if( (y-cy) > (x-cx) ) return false;
      //system("pause");
      return true;
    }
    return true;
  }

  /// **********************************************************************
  ///  ADD A NodeTree TO THE EXISTING ONE
  /// **********************************************************************
  void add(NodeTree Nt){
    unsigned int size0=n_.size();
    n_.resize(size0+Nt.n_.size());
    r_.resize(size0+Nt.n_.size());
    s_.resize(size0+Nt.n_.size());
    i_.resize(size0+Nt.n_.size());
    for(unsigned int m=0;m<Nt.n_.size();m++){
      n_[size0+m]=Nt.n_[m];
      s_[size0+m]=Nt.s_[m];
      r_[size0+m]=Nt.r_[m];
      for(unsigned int j=0;j<Nt.i_[m].size();j++){
        i_[size0+m].push_back(Nt.i_[m][j]+size0);
      }
    }
  };
  /// **********************************************************************
  ///  TRANSLATION OF A NodeTree
  /// **********************************************************************
  NodeTree translation(point2d p){
    NodeTree Nt=*this;
    for(unsigned int m=0;m<Nt.n_.size();m++){
      Nt.n_[m]=n_[m]+p;
    }
    return Nt;
  };

  /// **********************************************************************
  /// GENERATION OF A SYMMETRY
  /// **********************************************************************
  void SymmetryGeneration(point2d pc, int SymmetryType){
    if(SymmetryType==0) return;
    int Nangles;
    if(SymmetryType==1 || SymmetryType==4) Nangles=2;
    else if(SymmetryType==2 || SymmetryType==5) Nangles=4;
    else if(SymmetryType==3 || SymmetryType==6) Nangles=8;
    else return;

    NodeTree Nt=*this;

    for(int k=1;k<Nangles;k++){
      double alpha=-(2*k)*M_PI/Nangles;
      double a=cos(alpha);
      double b=sin(alpha);
      double c=-pc.x*a-pc.y*b;
      if(SymmetryType==1 || SymmetryType==4){
        for(unsigned int m=0;m<Nt.n_.size();m++){
          Nt.n_[m].x = pc.x-(Nt.n_[m].x-pc.x);
        }
      }
      else{
        for(unsigned int m=0;m<Nt.n_.size();m++){
          float t=2.*(a*Nt.n_[m].x+b*Nt.n_[m].y+c);
          Nt.n_[m] = point2d(Nt.n_[m].x-t*a,Nt.n_[m].y-t*b);
        }
      }
      unsigned int size0=n_.size();
      n_.resize(size0+Nt.n_.size());
      r_.resize(size0+Nt.n_.size());
      s_.resize(size0+Nt.n_.size());
      i_.resize(size0+Nt.n_.size());
      for(unsigned int m=0;m<Nt.n_.size();m++){
        n_[size0+m]=Nt.n_[m];
        s_[size0+m]=Nt.s_[m];
        r_[size0+m]=Nt.r_[m];
        for(unsigned int j=0;j<Nt.i_[m].size();j++){
          i_[size0+m].push_back(Nt.i_[m][j]+size0);
        }
      }

//      vector< point2d > n_; /// node coordinates
//  vector< float > r_; /// radius of the circle of each node
//  vector< float > s_; /// smoothing factor of the node connection
//  vector< vector<int> > i_; /// index of the nodes connected with the current node
//
    }
    //node_angles_computation();

  }

  /// *************************************************************************
  /// PROJECTION OF A POINT TO THE LINES LIMITING THE BASIC SYMMETRY AREA
  /// *************************************************************************
  bool point_line_projection(point2d &p,point2d c,float min_distance,int SymmetryType,float T0_limit=640.){
    if(SymmetryType==0) return true;
    bool projected=false;
    if(SymmetryType==1){
      if((p.x-c.x)<min_distance){ p.x=c.x; projected=true;}
      return projected;
    }

    if(SymmetryType==2 || SymmetryType==5){
      if((p.y-c.y)<min_distance){ p.y=c.y; projected=true;}
      else if( SymmetryType==5 &&  (T0_limit-p.y)<min_distance){p.y=T0_limit; projected=false;}
      if((p.x-c.x)<min_distance){ p.x=c.x; projected=true;}
      else if(SymmetryType==5 && (T0_limit-p.x)<min_distance){p.x=T0_limit; projected=false;}
      return projected;
    }
    if(SymmetryType==3 || SymmetryType==6){
      if(SymmetryType==6){
        if( (p-point2d(T0_limit,T0_limit)).norm2()<(min_distance*min_distance) ){
          p.x=T0_limit;
          p.y=T0_limit;
          return false;
        }
        if( (T0_limit-p.x)<min_distance){p.x=T0_limit; projected=false;}
      }
      point2d q=p-c;
      if(q.norm2()<(min_distance*min_distance)){
        p.x=c.x;
        p.y=c.y;
        return true;
      }
      double norm = 0.7071067812*(q.x-q.y);
      if(norm<min_distance){
        float x=p.x;
        p.x=0.5*(x+c.x+p.y-c.y);
        p.y=0.5*(x-c.x+p.y+c.y);
        if(SymmetryType==6){
          if( (p-point2d(T0_limit,T0_limit)).norm2()<(min_distance*min_distance) ){
            p.x=T0_limit;
            p.y=T0_limit;
            return false;
          }
        }
        return true;
      }
      if(q.y<min_distance){ p.y=c.y; return true;}
    }
    if(SymmetryType==4){
      if((p.x-c.x)<min_distance){ p.x=c.x; projected=true;}
      else if( (T0_limit-p.x)<min_distance){p.x=T0_limit; projected=false;}
      if( (T0_limit-p.y)<min_distance){p.y=T0_limit; projected=false;}
      else if( (p.y-(2.*c.y-T0_limit))<min_distance){p.y=2.*c.y-T0_limit; projected=true;}
    }
    return projected;
  }

  /// ****************************************************
  ///    INSERTION OF NEW NODES
  /// ****************************************************
  bool insert_nodes(
  triangulation &T0,
  float min_segment_length,float max_segment_length,
  float distance_segment_factor,
  float min_radius,float max_radius,
  float distance_node_factor,
  float smooth_factor,
  int type,
  float gamma,
  int Nbranches,
  vector<double> &angles,
  unsigned int MaxNnodes,
  unsigned int MaxNeigbors=100,
  bool allow_branche_intersection=false,
  int SymmetryType=0, /// ==0 no symmetry, ==1->180º, ==2->90º,==3->45º
  float node_radius_reduction_factor=1.) /// reduction of node radius factor
  {
    #ifdef DEBUG
      printf("insert_nodes() init...n_.size()=%d...\n",(int) n_.size());
      int size0=n_.size();
    #endif


    /// WE CHECK IF THE SET OF NODES IS EMPTY
    if(n_.size()==0) return false;

//    /// WE CHECK THE SIZE OF INDEX NEIGHBORHOOD AND ANGLES
//    if(a_.size()!=i_.size()) a_.resize(i_.size());
//    for(int k=0;k<i_.size();k++){
//      if(a_[k].size()!=i_[k].size()){
//        node_angles_computation(k);
//      }
//    }

    //printf("SymmetryType=%d\n",SymmetryType); system("pause");

    distance_node_factor*=max_radius*max_radius;
    distance_segment_factor*=max_radius;

    /// WE DEFINE SOME CONSTRAINT CRITERIUM FOR THE NEW NODES
    bool check_distance_nodes=true; /// the distance between nodes has to be large enough
    bool check_segment_intersection=true; /// the new segment can not intersect any other segment
    bool check_segment_distance=true; /// the distance of the new point to the segments has to be large enough
    if(allow_branche_intersection==true){
      check_segment_intersection=false;
      check_segment_distance=false;
    }

    /// WE DEFINE SOME AUXILIARY VECTOR
    vector<float> Distance2Center(n_.size()); /// squared distance of the node to the first node
    for(unsigned int k=0;k<n_.size();k++){ Distance2Center[k]=(n_[k]-n_[0]).norm2();   }

    /** index of nodes sorted by Distance2Center[] values. In this vector only contains active nodes, that is,
    nodes where we can try to add a new node*/
    vector<int> ActiveNodeIndex;


    /// WE START TO INCLUDE NEW NODES IN THE TREE
    #ifdef DEBUG_insert_nodes
      printf("e ActiveNodeIndex.size()=%d Nbranches=%d MaxNnodes=%d\n",(int) ActiveNodeIndex.size(),Nbranches,MaxNnodes);
    #endif

    unsigned int nsize0=0;
    while( (ActiveNodeIndex.size()>0 || Nbranches >0) && n_.size()< MaxNnodes){

      if(ActiveNodeIndex.size()==0){
        if(nsize0>=n_.size()) break;
        Nbranches--;
        ActiveNodeIndex=vector<int>(n_.size()-nsize0);
        for(unsigned int k=nsize0;k<n_.size();k++){ ActiveNodeIndex[k-nsize0]=k;}
        for(unsigned int k=0;k<ActiveNodeIndex.size();k++){
          for(unsigned int j=k+1;j<ActiveNodeIndex.size();j++){
            if(Distance2Center[ActiveNodeIndex[j]]<Distance2Center[ActiveNodeIndex[k]]){
              unsigned int temp=ActiveNodeIndex[j]; ActiveNodeIndex[j]=ActiveNodeIndex[k]; ActiveNodeIndex[k]=temp;
            }
          }
        }
        nsize0=n_.size();
      }

      /// WE CHOOSE RANDOMLY AN ACTIVE NODE TO ASSOCIATE A NEW NODE
      #ifdef DEBUG_insert_nodes
      printf("f\n");
      #endif
      unsigned int pos0=myrand(1.)*ActiveNodeIndex.size();
      pos0=pos0%ActiveNodeIndex.size();
      unsigned int pos=ActiveNodeIndex[pos0];
      if(i_[pos].size()>=MaxNeigbors){
         ActiveNodeIndex.erase(ActiveNodeIndex.begin()+pos0,ActiveNodeIndex.begin()+pos0+1);
         continue;
      }
      //printf("n_.size()=%d, ActiveNodeIndex.size()=%d, pos=%d\n",(int) n_.size(),(int) ActiveNodeIndex.size(),pos);
      //system("pause");

      /// WE CHOOSE RANDOMLY AN INITIAL ANGLE FOR THE VECTOR BEWTEEN THE ACTIVE NODE AND THE NEW ONE
      #ifdef DEBUG_insert_nodes
      printf("g\n");
      #endif
      double angle0 = -M_PI+ (double) (2.*M_PI)*rand()/RAND_MAX;
      angle0=0;
//      if(type==0 && n_.size()>1){
//        double angle2=atan2( n_[n_.size()-1].y-n_[0].y, n_[n_.size()-1].x-n_[0].x );
//        if(angle2+0.01<M_PI/2) angle0=M_PI/16;
//        else if (angle2-0.01>M_PI/2) angle0=-M_PI/16;
//      }
      /// WE TRY TO INSERT A NEW NODE USING THE ORIENTATION angle0. IF WE CAN'T WE TRY AGAIN USING
      /// NEW ORIENTATIONS IN THE ANGLE RANGE [angle0,angle0+2*M_PI];
      bool new_point_inserted=false; /// flat to control if a new point has been inserted
      //for(double angle=angle0;angle<angle0+2*M_PI;angle+=0.1){
      unsigned int m0=myrand(gamma)*angles.size();
      m0=m0%angles.size();
      //printf("ok1\n");
      #ifdef DEBUG_insert_nodes
      printf("h\n");
      #endif
      for(unsigned int m1,mm=0;mm<angles.size();mm++){
        #ifdef DEBUG_insert_nodes
        //printf("WE COMPUTE THE NEW POINT ORIENTATION WITH RESPECT TO THE ACTIVE NODE\n");
        printf("1\n");
        #endif
        if(mm<=m0) m1=m0-mm;
        else m1=mm;
        //double angle=angle0+angles[(m0+mm)%angles.size()];
        double angle=angle0+angles[m1];
//        if(angle>M_PI) angle-=2.*M_PI;
//        if(angle<-M_PI) angle+=2.*M_PI;
//        /// WE CHECK THE ANGLE DIFFERENCE
//        unsigned int j=0;
//        for(;j<a_[pos].size();j++){
//          //printf("angle=%lf a_[pos][j]=%lf\n",angle,a_[pos][j]);
//          if(fabs(angle_dif(angle,a_[pos][j]))<min_angle_difference || fabs(angle_dif(a_[pos][j],angle))<min_angle_difference ) break;
//        }
//        if(a_[pos].size()>0 && j<a_[pos].size()) continue;



        point2d q=point2d(cos(angle),sin(angle));
        //printf("q=(%lf,%lf)\n",q.x,q.y);

        #ifdef DEBUG_insert_nodes
        //printf("WE CHECK THE DISTANCE OF THE OTHER NEIGBORS OF THE POINT WITH RESPECT TO THE NEW POINT LINE\n");
        printf("2\n");
        #endif
        if(true){
          bool segment_intersection=false;
          for(unsigned int j=0;j< i_[pos].size();j++){
            point2d q1=n_[i_[pos][j]];//-q*radius;
            float angle2=atan2(q1.y-n_[pos].y,q1.x-n_[pos].x);
            if( fabs(angle_dif(angle,angle2))<0.087 || fabs(angle_dif(angle2,angle))<0.087 ){
              //printf("angle=%lf,angle2=%lf\n",angle,angle2);
              segment_intersection=true;
              break;
            }
        }
        if(segment_intersection==true) continue;
        //printf("segment_intersection==false\n");
      }

        #ifdef DEBUG_insert_nodes
        //printf("WE COMPUTE THE NODE RADIUS AND SEGMENT LENGTH\n");
        printf("3\n");
        #endif
        float r0=min_radius+myrand(1.)*(max_radius-min_radius); // (0.95+0.05*myrand(1.));
        if(type==1) r0=r_[pos]*node_radius_reduction_factor;
        float s0=min_segment_length+myrand(1.)*(max_segment_length-min_segment_length);// *(0.5+0.5*myrand(1.));

        #ifdef DEBUG_insert_nodes
        //printf("WE COMPUTE THE NEW NODE\n");
        printf("4\n");
        #endif
        point2d pnew=n_[pos]+q*s0;
        point_line_projection(pnew,point2d(512.,512.),min_segment_length*0.5,SymmetryType);
//printf("m0=%d, angles.size()=%d, pnew=(%lf,%lf), n_[0]=(%lf,%lf)\n",m0,angles.size(),pnew.x,pnew.y,n_[0].x,n_[0].y);
        #ifdef DEBUG_insert_nodes
        //printf("WE CHECK IF THE NEW NODE IS INSIDE T0\n");
        printf("5\n");
        #endif
        if(T0.inclusion_test(pnew)<0 || InsideSymmetricRegion(pnew.x,pnew.y,512.,512.,SymmetryType)==false) continue;

        #ifdef DEBUG_insert_nodes
        //printf("/// WE CHECK IF THE DISTANCE WITH RESPECT TO OTHER NODES IS LARGE ENOUGH\n");
        printf("6\n");
        #endif
        if(check_distance_nodes==true){
          unsigned int k=0;
          for(;k<n_.size();k++){
            if( (n_[k]-pnew).norm2()<distance_node_factor) break; //6*max_radius*max_radius
          }
          if(k<n_.size()) continue;
        }
        #ifdef DEBUG_insert_nodes
        //printf("WE CHECK IF THE NEW SEGMENT INTERSECTS OTHER SEGMENTS\n");
        printf("7\n");
        #endif
        if(check_segment_intersection==true){
          bool segment_intersection=false;
          for(unsigned int k=0;k< i_.size();k++){
            if(k==pos) continue;
            for(unsigned int j=0;j< i_[k].size();j++){
            point2d q1=n_[pos];//-q*radius;
            point2d q2=pnew+q*r0;
            if( check_cross(n_[k],n_[i_[k][j]],q1,q2)==true){
              k=i_.size();
              segment_intersection=true;
              break;
            }
          }
        }
        if(segment_intersection==true) continue;
        //printf("segment_intersection==false\n");
      }
//printf("ok2\n");
      if(check_segment_distance==true){
        bool segment_distance=false;
        for(unsigned int k=0;k< i_.size();k++){
          for(unsigned int j=0;j< i_[k].size();j++){
            if( distance_point2segment(n_[k],n_[i_[k][j]],pnew)<distance_segment_factor){ //2.5*radius
              //printf("pnew=(%1.2lf,%1.2lf), n_[k]=(%1.2lf,%1.2lf), n_[i_[k][j]]=(%1.2lf,%1.2lf),  dist=%1.2lf\n",pnew.x,pnew.y,n_[k].x,n_[k].y,n_[i_[k][j]].x,n_[i_[k][j]].y,distance_point2segment(pnew,n_[k],n_[i_[k][j]]));
              //if(pos==0) system("pause");
              k=i_.size();
              segment_distance=true;
              break;
            }
          }
        }
        if(segment_distance==false){
          for(unsigned int k=0;k<n_.size();k++){
            if(pos==k) continue;
            if( distance_point2segment(n_[pos],pnew,n_[k])<distance_segment_factor){
              segment_distance=true;
              break;
            }
          }
        }
        if(segment_distance==true) continue;
        //printf("distance to segment ok\n");
      }
//printf("ok3\n");
      #ifdef DEBUG_insert_nodes
      //printf("WE INCLUDE THE NEW POINT\n");
      printf("8\n");
      #endif
      new_point_inserted=true;
      //Nattempts=0;
//printf("pnew=(%lf,%lf) d=%lf\n",pnew.x,pnew.y,(n_[pos]-pnew).norm());
//if(pnew.x>1023 || pnew.x<0 || pnew.y>1023 || pnew.y<0 ) system("pause");
      n_.push_back(pnew);
      r_.push_back(r0);
      //update_index_angle(pos,n_.size()-1,angle);
      i_[pos].push_back(n_.size()-1);
      i_.push_back(vector<int>(1,pos));
      float norm2=(pnew-n_[0]).norm2();
      Distance2Center.push_back( norm2 );
      //if(n_.size()%10==1) printf("n_[%d]=(%lf,%lf)\n",(int) n_.size()-1,pnew.x,pnew.y);
//      node_angles_computation(pos);
//      node_angles_computation(n_.size()-1);
      //print(); system("pause");
//      if(type==0){
//        unsigned int k=0;
//        for(;k<ActiveNodeIndex.size();k++){
//          if(Distance2Center[ActiveNodeIndex[k]]>norm2){
//            ActiveNodeIndex.insert(ActiveNodeIndex.begin()+k,n_.size()-1);
//            break;
//          }
//        }
//        if(k==ActiveNodeIndex.size()) ActiveNodeIndex.push_back(n_.size()-1);
//      }

//      for(unsigned int k=0;k<ActiveNodeIndex.size();k++){ printf("%d ",ActiveNodeIndex[k]); } printf("\n");
//      for(unsigned int k=0;k<ActiveNodeIndex.size();k++){ printf("%1.1lf ",Distance2Center[ActiveNodeIndex[k]]); } printf("\n"); system("pause");
      #ifdef DEBUG_insert_nodes
      //printf("NEW POINT INCLUDED\n");
      printf("9\n");
      #endif
      break;
    }
    #ifdef DEBUG_insert_nodes
      printf("0\n");
    #endif
    if(new_point_inserted==false && pos0>=0 && pos0<ActiveNodeIndex.size()){
      //if(type==0) break;
      ActiveNodeIndex.erase(ActiveNodeIndex.begin()+pos0,ActiveNodeIndex.begin()+pos0+1);
    }
    #ifdef DEBUG_insert_nodes
      printf("a\n");
    #endif
  }
  #ifdef DEBUG_insert_nodes
      printf("b\n");
  #endif
  s_ = vector<float>(n_.size(),smooth_factor);
  #ifdef DEBUG_insert_nodes
      printf("c\n");
  #endif

  SymmetryGeneration(point2d(512.,512.),SymmetryType);
  if(checking()==false){
    printf("checking()=%d\n",false);
    n_.clear();
    i_.clear();
    r_.clear();
    a_.clear();
    return false;

    //system("pause");
  }
  #ifdef DEBUG
    printf(" %d nodes included...  end\n",(int) n_.size()-size0);
  #endif

  return true;

  //print();

};


  /// ****************************************************
  ///    INSERTION OF NEW NODES
  /// ****************************************************
  void joint_extrema(float join_probality){

    #ifdef DEBUG
      printf("joint_extrema() init...");
    #endif

    /// WE CHECK IF THE SET OF NODES IS TOO SMALL
    if(n_.size()<3 || join_probality<=0) return;

    /// AUXILARY VECTOR TO MEASURE NODE DISTANCES
    vector<float> D(n_.size());
    int cont=0; /// counter of new segments added

    /// WE GO TROUGH NODES WITH JUST ONE NEIGHBOR
    for(unsigned m=0;m<i_.size();m++){
      if(i_[m].size()!=1 || myrand(1.)>join_probality) continue;
      unsigned int m2=i_[m][0]; /// index of neigbor node

      /// COMPUTATION OF THE VECTOR JOINING THE POINT AND ITS NEIGHBOR
      float co=n_[m2].x-n_[m].x;
      float si=n_[m2].y-n_[m].y;
      float norm=co*co+si*si;
      if(norm<=0) continue;
      norm=sqrtf(norm);
      co/=norm;
      si/=norm;

      /// COMPUTATION OF SQUARED DISTANCE TO THE NODE
      for(unsigned int n=0;n<i_.size();n++){
        if(n==m || n==m2){
            D[n]=1e20;
            continue;
        }
        D[n]=(n_[n]-n_[m]).norm2();
      }

      /// MAX SQUARED DISTANCE TO THE NODE TO LOOK FOR NEW NEIGHBORS
      float minD0=100*norm*norm;
      float minD=0;

      /// LOOKING FOR NEW NEIGHBORS
      while(minD<minD0){
        /// COMPUTATION OF MIN NODE DISTANCE
        minD=1e20;
        unsigned int pos_min=0;
        for(unsigned int n=0;n<i_.size();n++){
          if(D[n]<minD){
            minD=D[n];
            pos_min=n;
          }
        }
        if(minD>minD0) break;
        /// CHECKING THE ORIENTATION OF THE NEW NEIGHBOR
        float co2=n_[pos_min].x-n_[m].x;
        float si2=n_[pos_min].y-n_[m].y;
        float ps=co*co2+si*si2;
        if(ps>0 && (ps*ps)>0.969846310*D[pos_min]){
          D[pos_min]=1e20;
          continue;
        }
        D[pos_min]=1e20;

        /// CHECKING THE INTERSECTION WITH OTHER NEIGHBOR-SEGMENTS IN THE TREE
        bool segment_intersection=false;
        for(unsigned int k=0;k<i_.size();k++){
          if(k==m || k==pos_min) continue;
          for(unsigned int j=0;j< i_[k].size();j++){
            if( i_[k][j]==(int) pos_min || i_[k][j]==(int) m) continue;
            //printf("%d %d %d %d %e\n",(int)m,(int)pos_min,(int)k,i_[k][j],minD);
            if( check_cross(n_[pos_min],n_[m],n_[k],n_[i_[k][j]])==true){
              //system("pause");
              segment_intersection=true;
              k=i_.size();
              break;
            }
          }
        }
        if(segment_intersection==true) continue;
        i_[m].push_back(pos_min);
        i_[pos_min].push_back(m);
        cont++;
        break;

      }
    }
    #ifdef DEBUG_joint_extrema
    printf("number of new segments=%d\n",cont);
    #endif

    node_angles_computation();

    #ifdef DEBUG
      printf("end\n");
    #endif

};

/// ****************************************************
  ///    INSERTION OF NEW NODES
  /// ****************************************************
  void random_generation(
  triangulation &T0,
  float min_segment_length,float max_segment_length,
  float distance_segment_factor,
  float min_radius,float max_radius,
  float distance_node_factor,
  float max_distance_node_factor,
  float min_smooth_factor,float max_smooth_factor,
  float gamma,
  double min_angle_distance,
  unsigned int MaxNnodes,
  unsigned int MaxNeigbors,
  unsigned int MaxAttempts,
  bool allow_branche_intersection=false){

  #ifdef DEBUG
    printf("random_generation() init.....");
  #endif

    /// WE CHECK IF THE SET OF NODES IS EMPTY
    if(T0.p.size()==0) return;

    vector<point2d> ex=T0.extrema();

    clear();

    distance_node_factor*=max_radius*max_radius;
    max_distance_node_factor*=max_radius*max_radius;
    distance_segment_factor*=max_radius;

    /// WE DEFINE SOME CONSTRAINT CRITERIUM FOR THE NEW NODES
    bool check_distance_nodes=true; /// the distance between nodes has to be large enough
    bool check_segment_intersection=true; /// the new segment can not intersect any other segment
    bool check_segment_distance=true; /// the distance of the new point to the segments has to be large enough
    if(allow_branche_intersection==true){
     check_segment_intersection=false;
     check_segment_distance=false;
    }

    unsigned int Attempts=0;
    while(Attempts++<MaxAttempts && n_.size()<MaxNnodes){
      //printf("Attempts=%d\n",Attempts);
      point2d pnew(ex[0].x+myrand(1.)*(ex[1].x-ex[0].x),ex[0].y+myrand(1.)*(ex[1].y-ex[0].y));
      while(T0.inclusion_test(pnew)<0){
        pnew=point2d(ex[0].x+myrand(1.)*(ex[1].x-ex[0].x),ex[0].y+myrand(1.)*(ex[1].y-ex[0].y));
      }

      if(n_.size()==0){
        n_.push_back(pnew);
        s_.push_back(min_smooth_factor+myrand(1.)*(max_smooth_factor-min_smooth_factor));
        r_.push_back(min_radius+myrand(1.)*(max_radius-min_radius));
        i_.resize(1);
        a_.resize(1);
        continue;
      }

      float r0=min_radius+myrand(1.)*(max_radius-min_radius);

      unsigned int pos0=myrand(gamma)*n_.size();
      for(unsigned int m=0;m<n_.size();m++){
        unsigned int pos=(pos0+m)%n_.size();
        if(i_[pos].size()>=MaxNeigbors) continue;
        double angle=atan2(pnew.y-n_[pos].y,pnew.x-n_[pos].x);
        /// WE CHECK ANGLE DISTANCE
        if(a_[pos].size()>0){
          unsigned int k=0;
          for(;k<a_[pos].size();k++){
            if( fabs(a_[pos][k]-angle)<min_angle_distance) break;
            if( 2*M_PI-fabs(a_[pos][k]-angle)<min_angle_distance) break;
            //printf("a_[pos][k]=%lf,angle=%lf,angle_dif(a_[pos][k],angle)=%lf\n",a_[pos][k],angle,angle_dif(a_[pos][k],angle));
          }

          //system("pause");
          if(k<a_[pos].size()) continue;
        }
        point2d q=point2d(cos(angle),sin(angle));

        /// WE CHECK THE DISTANCE WITH RESPECT TO OTHER NODES
        //if( (n_[pos]-pnew).norm2()<distance_node_factor ) continue;
        if(check_distance_nodes==true){
          if( (n_[pos]-pnew).norm2()>max_distance_node_factor) continue; //6*max_radius*max_radius
          unsigned int k=0;
          for(;k<n_.size();k++){
            if( (n_[k]-pnew).norm2()<distance_node_factor) break; //6*max_radius*max_radius

          }
          if(k<n_.size()) continue;
        }

        /// WE CHECK IF THE NEW SEGMENT INTERSECTS OTHER SEGMENTS
        if(check_segment_intersection==true){
          bool segment_intersection=false;
          for(unsigned int k=0;k< i_.size();k++){
            if(k==pos) continue;
            for(unsigned int j=0;j< i_[k].size();j++){
              point2d q1=n_[pos];//-q*radius;
              point2d q2=pnew+q*r0;
              if( check_cross(n_[k],n_[i_[k][j]],q1,q2)==true){
                k=i_.size();
                segment_intersection=true;
                break;
              }
            }
          }
          if(segment_intersection==true) continue;
          //printf("segment_intersection==false\n");
        }

        if(check_segment_distance==true){
          bool segment_distance=false;
          for(unsigned int k=0;k< i_.size();k++){
            for(unsigned int j=0;j< i_[k].size();j++){
              if( distance_point2segment(n_[k],n_[i_[k][j]],pnew)<distance_segment_factor){ //2.5*radius
                //printf("pnew=(%1.2lf,%1.2lf), n_[k]=(%1.2lf,%1.2lf), n_[i_[k][j]]=(%1.2lf,%1.2lf),  dist=%1.2lf\n",pnew.x,pnew.y,n_[k].x,n_[k].y,n_[i_[k][j]].x,n_[i_[k][j]].y,distance_point2segment(pnew,n_[k],n_[i_[k][j]]));
                //if(pos==0) system("pause");
                k=i_.size();
                segment_distance=true;
                break;
              }
            }
          }
          if(segment_distance==false){
            for(unsigned int k=0;k<n_.size();k++){
              if(pos==k) continue;
              if( distance_point2segment(n_[pos],pnew,n_[k])<distance_segment_factor){
                segment_distance=true;
                break;
              }
            }
          }
          if(segment_distance==true) continue;
          //printf("distance to segment ok\n");
        }

        //new_node_inserted=true;
        Attempts=0;
        n_.push_back(pnew);
        r_.push_back(r0);
        s_.push_back(min_smooth_factor+myrand(1.)*(max_smooth_factor-min_smooth_factor));
        i_[pos].push_back(n_.size()-1);
        i_.push_back(vector<int>(1,pos));
        a_[pos].push_back(angle);
        a_.push_back(vector<double>(1,atan2(n_[pos].y-pnew.y,n_[pos].x-pnew.x)));
        //printf("checking()=%d\n",checking());
        break;
      }
      //printf("n_.size()=%d\n",n_.size());
    }
    //printf("n_.size()=%d\n",(int) n_.size());


    //node_angles_computation();

  #ifdef DEBUG
    printf("end\n");
  #endif

};

void normalize(){
  //return;
  if(n_.size()<2) return;
  vector<point2d> eV(2,n_[0]);
  for(unsigned int k=0;k<n_.size();k++){
      //printf("n_[%d]=(%lf,%lf)\n",k,n_[k].x,n_[k].y);
      if( (n_[k].x-2.*r_[k])<eV[0].x ) eV[0].x=n_[k].x-2.*r_[k];
      if( (n_[k].x+2.*r_[k])>eV[1].x ) eV[1].x=n_[k].x+2.*r_[k];
      if( (n_[k].y-2.*r_[k])<eV[0].y ) eV[0].y=n_[k].y-2.*r_[k];
      if( (n_[k].y+2.*r_[k])>eV[1].y ) eV[1].y=n_[k].y+2.*r_[k];
    }
    //printf("\neV=(%lf,%lf),(%lf,%lf)\n",eV[0].x,eV[0].y,eV[1].x,eV[1].y);
    //if(eV[0].x<24 || eV[0].y<24 || eV[1].x >1000 || eV[1].y >1000){
      float xc=(eV[0].x+eV[1].x)/2;
      float yc=(eV[0].y+eV[1].y)/2;
      float scale=(eV[1].x-eV[0].x)/2;
      if(scale<( (eV[1].y-eV[0].y)/2)) scale=(eV[1].y-eV[0].y)/2;
      scale=500./scale;
      //printf("xc=%lf yc=%lf scale=%lf\n",xc,yc,scale);
      for(unsigned int k=0;k<n_.size();k++){
        n_[k].x=512.+(n_[k].x-xc)*scale;
        n_[k].y=512.+(n_[k].y-yc)*scale;
        r_[k]=r_[k]*scale;
     }
  //}

}



};


/// BRANCHED SHAPE GENERATOR
NodeTree BranchedShapeGenerator(
triangulation &Tinit,
int &shape_type,
int &random_seed,
float &average_node_distance,
float &average_radius,
float &smoothing_factor,
int &branched_finished_in_peak,
float &joint_extrema_probability,
int &MaxNeigbors,
bool &allow_branche_intersection,
float &node_radius_reduction_factor,
int &MaxAngles,
unsigned int MaxNodes
);


/// RANDOM BRANCHED SHAPE SELECTION
int get_branched_shape_index();
int get_branched_shape_index2();

/// INITIALIZATION OF NodeTree PARAMETERS
void NodeTree_parameters_initialization(
int argc,char *argv[],int &shape_type, int &random_seed, float &average_node_distance,
float &average_radius, float &smoothing_factor, int &branched_finished_in_peak,
int &Xperiods, int &Yperiods, float &joint_extrema_probability, int &MaxNeigbors,
bool &allow_branche_intersection, float &node_radius_reduction_factor, int &MaxAngles,
unsigned int &MaxNodes,
bool &roud_termination);

/// NodeTree RANDOM GENERATOR
 NodeTree NodeTree_random_generator(int shape_type=-1);

/// GENERATION OF A SIMILAR SHAPE
NodeTree NodeTree_similar_generation(NodeTree &nTo);

#endif

