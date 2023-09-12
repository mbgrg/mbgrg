/**
 * \file triangulation.cpp
 * \brief functions  to manage Delauny triangulations
 * \author Luis Alvarez \n \n
*/


#include <iostream>
#include <vector>
#include <math.h>
#include "triangulation.h"
#include <climits>
#include <time.h>

/// SOME MACROS WE USE
#define ami_PI 3.14159265358979323846
#define CERO 1e-7

#define ami_calloc2d(direccion,tipo,height,width) {int ml,mk; \
          direccion=(tipo **) malloc(sizeof(tipo *)*(height)); \
          direccion[0]=(tipo *)malloc(sizeof(tipo)*(width)*(height)); \
          for(ml=0;ml<(height);ml++) direccion[ml]=&(direccion[0][ml*(width)]); \
          for(ml=0;ml<height;ml++) for(mk=0;mk<width;mk++) direccion[ml][mk]=0; \
}
#define ami_malloc1d(direccion,tipo,size) {direccion=(tipo *) malloc(sizeof(tipo)*(size)); if(direccion==NULL) {printf("fallo malloc\n"); exit(-1);}}
#define ami_free2d(direccion) { free(direccion[0]); free(direccion); }



/******************************************************************
 RETURN DISTANCE FROM A POINT TO A SEGMENT. THE PROJECTION OF THE
 POINT IN THE LINE SEGMENT HAS TO BE INCLUDED IN THE SEGMENT
 ******************************************************************/
float distance_point2segment(point2d &v0,point2d &v1,point2d &v)
{
  float a=v1.x-v0.x;
  float b=v1.y-v0.y;
  float R1=a*(v.x-v0.x)+b*(v.y-v0.y);
  float R2=a*(v.x-v1.x)+b*(v.y-v1.y);
  if(R1*R2>0.) return 1e20;
  float norm=a*a+b*b;
  if(norm==0.){
    printf("problems in distance_point2segment(). norm==0\n");
  }
  return fabs(a*(v.y-v0.y)-b*(v.x-v0.x))/sqrtf(norm);

}


/***************************************************************/
/**         GAUSS METHOD TO SOLVE LINEAR SYSTEMS               */
/***************************************************************/
int ami_gauss(double **A,double *b,int N)
/* DEVUELVE LA SOLUCION EN EL VECTOR b */
     /*double **A,*b; int N; */
{
  double **PASO,max,paso,mul;
  int i,j,i_max,k;
  PASO=(double **)malloc(sizeof(double*)*N);
  for(i=0;i<N;i++)
    PASO[i]=(double *)malloc(sizeof(double)*(N+1));

  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      PASO[i][j]=A[i][j];
    }
    PASO[i][N]=b[i];
  }

  for(i=0;i<N;i++){
    max=fabs(PASO[i][i]);
    i_max=i;
    for(j=i;j<N;j++){
       if(fabs(PASO[j][i])>max){
         i_max=j; max=fabs(PASO[j][i]);
       }
    }
    if(max<10e-30){
      printf("Sistema no tiene Solucion 0\n");
      for(i=0;i<N;i++)
        free(PASO[i]);
      free(PASO);
      return(-1);
    }
    if(i_max>i){
      for(k=0;k<=N;k++){
        paso=PASO[i][k];
        PASO[i][k]=PASO[i_max][k];
        PASO[i_max][k]=paso;
      }
    }
    for(j=i+1;j<N;j++){
      mul=-PASO[j][i]/PASO[i][i];
      for(k=i;k<=N;k++) PASO[j][k]+=mul*PASO[i][k];
    }
  }
  if(fabs(PASO[N-1][N-1])<10e-30){
      printf("Linear system has no solution\n");
      for(i=0;i<N;i++)
       free(PASO[i]);
      free(PASO);
      return(-1);
    }

  for(i=N-1;i>0;i--){
    for(j=i-1;j>=0;j--){
      mul=-PASO[j][i]/PASO[i][i];
      for(k=i;k<=N;k++) PASO[j][k]+=mul*PASO[i][k];
    }
  }
  for(i=0;i<N;i++)
      b[i]=PASO[i][N]/PASO[i][i];

  for(i=0;i<N;i++)
    free(PASO[i]);
  free(PASO);
  return(0);
}

/** FUNCION TO CREATE CUBIC SPLINES */
int ami_crear_splines(vector<float> &x, vector<float> &y, float **s){
    int N=x.size();
	double **A, *b = (double*)malloc(sizeof(double)*N);
	float *h;
	int i;

	/* MEMORY ALLOCATION*/
    ami_calloc2d(A, double, N, N);
	ami_malloc1d(h,float,N);

	/* WE COMPUTE STEPS IN THE x-AXIS */
	for(i=0;i<(N-1);i++){
	  h[i]=x[i+1]-x[i];
	  if(h[i]<=0.){
		  printf("points in the x-axis are not ordered\n");
		  return -1;
	  }
    }

	/* Computation of a */
	for (i=0;i<N;i++) s[i][0] = y[i];

	/* Computation of the linear system */
	A[0][0] = 1.;
	b[0] = 0.;

	for (i=1;i<(N-1);i++)
	{
		A[i][i] = 2.0*(h[i-1]+h[i]);
		A[i][i-1] = h[i-1];
		A[i][(i+1)%N] = h[i];
		b[i] = 3.0*(s[(i+1)%N][0]-s[i][0])/h[i] -
			   3.0*(s[i][0]-s[i-1][0])/h[i-1];
	}

	A[N-1][N-1]=1.;
	b[N-1]=0.;

	/* we solve the linear system */
	if (ami_gauss(A,b,N) < 0) return -1;

	/* computation of b, c y d */
	for (i=0;i<(N-1);i++)
	{
		s[i][2] = b[i];
		s[i][1] = (s[(i+1)%N][0]-s[i][0])/h[i] -
			       h[i]*(2.0*b[i]+b[(i+1)%N])/3.0;
		s[i][3] = (b[(i+1)%N]-b[i])/(3.0*h[i]);
	}

	/* we free memory */
	free(b); free(h);
	ami_free2d(A);

	return 0;
}

/***********************************************************************************************/
/** FUNCTION TO CREATE CUBIS SPLINES ASSIGNING TANGENT DIRECTIONS IN THE FIRST AND LAST POINTS */
/***********************************************************************************************/
int ami_crear_splines(vector<float> &x, vector<float> &y, float **s,
float d0,  /// derivate in x[0]
float dN)  /// derivative in x[x.dim()-1]
{
    int N=x.size();
    if(N!=(int) y.size() || N<2) return -1;

    // BUILDING ANS SOLVING THE LINEAR SYSTEM
	double **A, *b = (double*)malloc(sizeof(double)*N);
	float *h;
	int i;

    ami_calloc2d(A, double, N, N);
	ami_malloc1d(h,float,N);

	for(i=0;i<(N-1);i++){
	  h[i]=x[i+1]-x[i];
	  if(h[i]<=0.){
		  printf("points are not ordered in the x-axis\n");
		  return -1;
	  }
    }

	for (i=0;i<N;i++) s[i][0] = y[i];

	A[0][0] = 2.*h[0]/3.;
	A[0][1] = h[0]/3.;
	b[0] = (y[1]-y[0])/h[0]-d0;

	for (i=1;i<(N-1);i++)
	{
		A[i][i] = 2.0*(h[i-1]+h[i]);
		A[i][i-1] = h[i-1];
		A[i][(i+1)%N] = h[i];
		b[i] = 3.0*(s[(i+1)%N][0]-s[i][0])/h[i] -
			   3.0*(s[i][0]-s[i-1][0])/h[i-1];
	}

	A[N-1][N-1]= 2.*h[N-2]/3.;
	A[N-1][N-2] = h[N-2]/3.;
	b[N-1]=dN-(y[N-1]-y[N-2])/h[N-2];

	/* We solve the linear system*/
	if (ami_gauss(A,b,N) < 0) return -1;

	/* we compute the coefficients b, c y d */
	for (i=0;i<(N-1);i++)
	{
		s[i][2] = b[i];
		s[i][1] = (s[(i+1)%N][0]-s[i][0])/h[i] -
			       h[i]*(2.0*b[i]+b[(i+1)%N])/3.0;
		s[i][3] = (b[(i+1)%N]-b[i])/(3.0*h[i]);
	}

	/* we free the memory */
	free(b); free(h);
	ami_free2d(A);

	return 0;
}

/*************************************************/
/** FUNCTION TO EVALUATE THE SPLINES IN A POINT  */
/*************************************************/
float ami_evaluar_splines(vector<float> &x,float **s,float z){

  int i,j,N=x.size();
  float t;
  /* COMPUTATION OF THE INTERVAL WHERE THE POINT IS */
  j=0;
  i=1;
  while(i<N && x[i]<z){
	i++;
	j++;
  }

  /* EVALUATION OF THE 3 DEGREE POLYNOMIAL */
  t=z-x[j];

  return s[j][0]+t*s[j][1]+t*t*s[j][2]+t*t*t*s[j][3];

}

/******************************************************/
/** FUNCTION TO EVALUATE SPLINE DERIVATIVE IN A POINT */
/******************************************************/
float ami_evaluar_derivada_splines(vector<float> &x,float **s,float z){

  int i,j,N=x.size();
  float t;
  /* COMPUTATION OF THE INTERVAL WHERE THE POINT IS */
  j=0;
  i=1;
  while(i<N && x[i]<z){
	i++;
	j++;
  }

  /* EVALUATION OF THE DERIVATIVE OF THE POLYNOMIAL */
  t=z-x[j];

  return s[j][1]+2*t*s[j][2]+3*t*t*s[j][3];

}

/********************************************************************************/
/**   AFFINITY CONSTRUCTOR USING AN ANGLE, AND  HORIZONTAL AND UNIFORM SCALING  */
/********************************************************************************/
affinity::affinity(float angle,float scale,float horizontal_scaling){
  cx=cy=0.;
  float co=cos(angle);
  float si=sin(angle);
  A[0][0]=horizontal_scaling*scale*co; A[0][1]=scale*si; A[1][0]=-horizontal_scaling*scale*si; A[1][1]=scale*co;
}

/**********************************************************************************/
/**   AFFINITY CONSTRUCTOR USING TWO ANGLES, AND  HORIZONTAL AND UNIFORM SCALING  */
/**********************************************************************************/
affinity::affinity(float angle1,float angle2,float scale,float horizontal_scaling){
  cx=cy=0.;
  float c1=cos(angle1);
  float s1=sin(angle1);
  float c2=cos(angle2);
  float s2=sin(angle2);
  float a=horizontal_scaling;

  A[0][0]=scale*(a*c1*c2-s1*s2);
  A[0][1]=scale*(c1*s2+a*c2*s1);
  A[1][0]=scale*(-c2*s1-a*c1*s2);
  A[1][1]=scale*(c1*c2-a*s1*s2);
}


/**********************************************************************************/
/** constructor of A satisfying :  A*(p2+p3)=(p0+p1)  && A*(p3-p2) || (p1-p0) */
/**********************************************************************************/
affinity::affinity(point2d p0,point2d p1,point2d p2,point2d p3){
    float ux=p1.x-p0.x;
    float uy=p1.y-p0.y;
    float vx=p3.x-p2.x;
    float vy=p3.y-p2.y;;
    float alpha=atan2(-ux*vy+uy*vx,ux*vx+uy*vy);//+_PI;
    float co=cos(alpha);
    float si=sin(alpha);

    cx=p0.x-co*p2.x+si*p2.y;
    cy=p0.y-si*p2.x-co*p2.y;
    A[0][0]=co; A[0][1]=-si;
    A[1][0]=si; A[1][1]=co;

//    printf("p0=(%1.2lf,%1.2lf) p1=(%1.2lf,%1.2lf)\n",p0.x,p0.y,p1.x,p1.y);
//    printf("p2=(%1.2lf,%1.2lf) p3=(%1.2lf,%1.2lf)\n",p2.x,p2.y,p3.x,p3.y);
//    point2d p2b=(*this)*p2;
//    point2d p3b=(*this)*p3;
//    printf("p2b=(%1.2lf,%1.2lf) p3b=(%1.2lf,%1.2lf)\n",p2b.x,p2b.y,p3b.x,p3b.y);
//    system("pause");
}

/**********************************************************************************************/
/** constructor of A satisfying :  (A*p10=p00+d or A*p20=p00+d)  && A*(p11-p10) || (p01-p00) &&
  the triangle A*(p10,p11,p12) does not intersects the triangle (p00,p01,p02) */
/**********************************************************************************************/
affinity::affinity(point2d p00,point2d p01,point2d p02,point2d p10,point2d p11,point2d p12,float d,int type){
    float R0=R(p00,p01,p02);
    affinity A((p00+p01)*0.5,p01,(p10+p11)*0.5,p11);
    point2d Ap12=A*p12;
    //affinity B(p00,p01,p11,p10);
    affinity B;
    if(type==0) B=affinity((p00+p01)*0.5,p01,(p10+p11)*0.5,p10);
    else{
      point2d p=(p01-p00); p.normalize();
      float a=-p.y; float b=p.x; float c=-(a*p00.x+b*p00.y);
      B=affinity(1-2*a*a,-2.*a*b,-2.*a*b,1.-2.*b*b,-2.*a*c,-2.*b*c)*A;
    }
    point2d Bp12=B*p12;
    if( R0*R(p00,p01,Ap12)<R0*R(p00,p01,Bp12)){
       if(d>=0.) (*this)=A;
       else (*this)=B;
    }
    else{
      if(d>=0.) (*this)=B;
      else (*this)=A;
    }
    if(d!=0.){
      float ux=p01.x-p00.x;
      float uy=p01.y-p00.y;
      float norm=sqrt(ux*ux+uy*uy);
      ux/=norm;
      uy/=norm;
      point2d t(-uy,ux);
      point2d p00t=p00+t;
      if( R0*R(p00,p01,p00t)>0.){ cx-=t.x*d; cy-=t.y*d; }
      else {cx+=t.x*d; cy+=t.y*d;}
    }
}


/**********************************************************
ADDITION (JOINT) OF 2 TRIANGULATIONS
**********************************************************/
triangulation triangulation::operator+(const triangulation &Tr)
{
  triangulation T; // new triangulation

  // WE SET T TO THE CURRENT TRIANGULATION
  T.p=p;
  T.fr=fr;
  T.t=t;

  // WE ADD THE SECOND TRIANGULATION
  int Np=p.size();
  int Nt=t.size();
  T.p.insert(T.p.begin()+Np,Tr.p.begin(),Tr.p.end());
  T.fr.insert(T.fr.begin()+Np,Tr.fr.begin(),Tr.fr.end());
  T.t.insert(T.t.begin()+Nt,Tr.t.begin(),Tr.t.end());

  // WE UPDATE THE SECOND TRIANGULATION INDEX
  for(int k=Nt;k<(int) T.t.size();k++){
    for(int l=0;l<3;l++) T.t[k].v[l]+=Np;
  }

  return(T);
}

/**********************************************************
SUBSTRACTION OF 2 TRIANGULATIONS
**********************************************************/
triangulation triangulation::operator-(triangulation &Tr)
{
  triangulation T=(*this); //new triangulation

  // WE COMPUTE TRIANGLE NEIGHBORHOOD INFORMATION
  if(Tr.t.size() !=Tr.n.size()) Tr.update_n();

  vector<point2d> p1; // vector to keep external triangulation edge of the first triangulation
  // WE INCLUDE IN THE FIRST TRIANGULATION THE CONTOUR POINTS OF THE SECOND TRIANGULATION
  for(int k=0;k<(int) Tr.t.size();k++){
    for(int l=0;l<3;l++){
      if(Tr.n[k].v[l]<0){ // contour edge
        Delauny(T,Tr.p[Tr.t[k].v[l]],1);
        Delauny(T,Tr.p[Tr.t[k].v[(l+1)%3]],1);
        p1.push_back(Tr.p[Tr.t[k].v[l]]);
        p1.push_back(Tr.p[Tr.t[k].v[(l+1)%3]]);
      }
    }
  }

  // WE INCLUDE THE COLLECTION OF EXTERNAL EDGES
  Delauny(T,p1,1);

  // WE REMOVE EXTERNAL TRIANGLES IN THE FIRST MODIFIED TRIANGULATION
  for(int k=0;k<(int) T.t.size();k++){
    point2d v=(T.p[T.t[k].v[0]]+T.p[T.t[k].v[1]]+T.p[T.t[k].v[2]])/3.; // triangle center
    for(int m=0;m<(int) Tr.t.size();m++){
      float R1,R2,R3;
      R1=R(Tr.p[Tr.t[m].v[0]],Tr.p[Tr.t[m].v[1]],v);
      R2=R(Tr.p[Tr.t[m].v[1]],Tr.p[Tr.t[m].v[2]],v);
      R3=R(Tr.p[Tr.t[m].v[2]],Tr.p[Tr.t[m].v[0]],v);
      if( (R1>=-1e-10 && R2>=-1e-10 && R3>=-1e-10) || (R1<=1e-10 && R2<=1e-10 && R3<=1e-10) ){ // the triangle center is included in 1 triangle of the initial triangulation
        T.t.erase(T.t.begin()+k);
        k--;
        break;
      }
    }
  }

  return(T);
}

/**********************************************************
INTERSECTION OF 2 TRIANGULATIONS
**********************************************************/
triangulation triangulation::operator*(triangulation &Tr)
{
  triangulation T=(*this); ;

  // WE COMPUTE TRIANGLE NEIGHBORHOOD INFORMATION
  if(Tr.t.size() !=Tr.n.size()) Tr.update_n();

  vector<point2d> p1; // vector to keep external triangulation edge of the first triangulation
  // WE INCLUDE IN THE FIRST TRIANGULATION THE CONTOUR POINTS OF THE SECOND TRIANGULATION
  for(int k=0;k<(int) Tr.t.size();k++){
    for(int l=0;l<3;l++){
      if(Tr.n[k].v[l]<0){
        Delauny(T,Tr.p[Tr.t[k].v[l]],1);
        Delauny(T,Tr.p[Tr.t[k].v[(l+1)%3]],1);
        p1.push_back(Tr.p[Tr.t[k].v[l]]);
        p1.push_back(Tr.p[Tr.t[k].v[(l+1)%3]]);
      }
    }
  }


 // WE INCLUDE IN THE NEW TRIANGULATION THE CONTOUR POINTS OF THE FIRST TRIANGULATION
 if(t.size() !=n.size()) update_n();

  vector<point2d> p2;
  for(int k=0;k<(int) t.size();k++){
    for(int l=0;l<3;l++){
      if(n[k].v[l]<0){
        Delauny(T,p[t[k].v[l]],1);
        Delauny(T,p[t[k].v[(l+1)%3]],1);
        p2.push_back(p[t[k].v[l]]);
        p2.push_back(p[t[k].v[(l+1)%3]]);
      }
    }
  }
  // WE NEED A LOOP TO ENSURE THAT THE EDGE CONTOUR ARE INCLUDED IN THE NEW TRIANGULATION
  int Np=T.p.size()-1;
  while(Np<(int) T.p.size()){
    Np=T.p.size();
    Delauny(T,p2,1);
    Delauny(T,p1,1);
  }


  // WE REMOVE TRIANGLES WHICH ARE NOT IN BOTH TRIANGULATIONS
  for(int k=0;k<(int) T.t.size();k++){ // we check second triangulation
    point2d v=(T.p[T.t[k].v[0]]+T.p[T.t[k].v[1]]+T.p[T.t[k].v[2]])/3.; //triangle center
    int m,l;
    for(m=0;m<(int) Tr.t.size();m++){
      float R1,R2,R3;
      R1=R(Tr.p[Tr.t[m].v[0]],Tr.p[Tr.t[m].v[1]],v);
      R2=R(Tr.p[Tr.t[m].v[1]],Tr.p[Tr.t[m].v[2]],v);
      R3=R(Tr.p[Tr.t[m].v[2]],Tr.p[Tr.t[m].v[0]],v);
      if( (R1>=-1e-10 && R2>=-1e-10 && R3>=-1e-10) || (R1<=1e-10 && R2<=1e-10 && R3<=1e-10) ){ //triangle is included in the second triangulation
        break;
      }
    }
    if(m==(int) Tr.t.size()){ //we remove the triangle if its center is not included in the second triangulation
      T.t.erase(T.t.begin()+k);
      k--;
      continue;
    }

    for(l=0;l<(int) t.size();l++){ // we check first triangulation
      float R1,R2,R3;
      R1=R(p[t[l].v[0]],p[t[l].v[1]],v);
      R2=R(p[t[l].v[1]],p[t[l].v[2]],v);
      R3=R(p[t[l].v[2]],p[t[l].v[0]],v);
      if( (R1>=-1e-10 && R2>=-1e-10 && R3>=-1e-10) || (R1<=1e-10 && R2<=1e-10 && R3<=1e-10) ){
        break;
      }
    }
    if(l==(int) t.size()){
      T.t.erase(T.t.begin()+k);
      k--;
    }
  }

  return(T);
}


/**************************************************
FUNCTION TO READ TRIANGULATION FROM FILE
**************************************************/
triangulation::triangulation(char *name){
    int i,np,nt;
    FILE *f;
    if(f=fopen(name,"r"),!f){
      printf("Problems opening triangulation file %s",name); { int ii_; scanf("%d",&ii_);}//system("pause");
      return ;
    }
    fscanf(f,"%d %d\n",&np,&nt);
    //printf("np=%d,nt=%d\n",np,nt);
    p.resize(np);
    fr.resize(np);
    t.resize(nt);
//    v.resize(nt);
    for(i=0;i<np;i++){
      double x,y;
      int f2;
      //fscanf(f,"%lf %lf %d\n",&(p[i].x),&(p[i].y),&(fr[i]));
      fscanf(f,"%lf %lf %d\n",&x,&y,&f2);
      p[i].x=x; p[i].y=y; fr[i]=f2;
      //printf("p[%d]=(%lf,%lf) fr=%d\n",i,p[i].x,p[i].y,fr[i]);
    }
    for(i=0;i<nt;i++){
      int t0,t1,t2;
      //fscanf(f,"%d %d %d \n",&(t[i].v[0]),&(t[i].v[1]),&(t[i].v[2]));
      fscanf(f,"%d %d %d \n",&t0,&t1,&t2);
      t[i].v[0]=t0;
      t[i].v[1]=t1;
      t[i].v[2]=t2;
      //printf("t[%d]=(%d,%d,%d)\n",i,t[i].v[0],t[i].v[1],t[i].v[2]);
    }
//    for(i=0;i<nt;i++)
//      fscanf(f,"%d %d %d \n",&(v[i].v[0]),&(v[i].v[1]),&(v[i].v[2]));
    fclose(f);
    fix_triangle_orientation();
    update_n();
  }

/******************************************************
FUNCTION TO WRITE A TRIANGULATION IN DISK
******************************************************/
int triangulation::write(const char *name){
  int i;
  FILE *f;
  if(f=fopen(name,"w"),!f) return 1;
  //printf("%d %d\n",p.size(),t.size());
  fprintf(f,"%d %d\n",(int) p.size(),(int) t.size());
  for(i=0;i<(int) p.size();i++){
    //printf("%lf %lf %d\n",p[i].x,p[i].y,fr[i]);
    fprintf(f,"%lf %lf %d\n",p[i].x,p[i].y,fr[i]);
  }
  for(i=0;i<(int) t.size();i++){
    //printf("%d %d %d \n",t[i].v[0],t[i].v[1],t[i].v[2]);
    fprintf(f,"%d %d %d \n",t[i].v[0],t[i].v[1],t[i].v[2]);
  }
  fclose(f);
  return(0);
}

/***********************************************
FUNCTION TO ERASE POINT AND ASSOCIATED TRIANGLES
***********************************************/
void triangulation::erase_point(int i /** index of point to erase */)
{
  if(i>=(int) p.size()) return;
  // we remove the point
  p.erase (p.begin()+i);
  if(i<(int) fr.size()) fr.erase (fr.begin()+i);

  for(int k=0;k<(int) t.size();k++){ // we remove triangles which include the point
    if(t[k].v[0]==i || t[k].v[1]==i || t[k].v[2]==i){ t.erase(t.begin()+k); k--;}
  }

  for(int k=0;k<(int) t.size();k++){ // we update triangle index.
    for(int l=0;l<3;l++) if(t[k].v[l]>i) t[k].v[l]--;
  }

}


/***********************************
 RETURN  TRIANGLE AREA
 **********************************/
double AREA(point2d &v0,point2d &v1,point2d &v2)
{
  double result;
  result= v1.x*v2.y-v2.x*v1.y-v0.x*v2.y+v0.y*v2.x+v0.x*v1.y-v0.y*v1.x;
  result= 0.5*result;
  if(result>0)
    return(result);
  else
    return(-result);
}


/****************************************
 RETURN SCALAR PRODUCT : (v-v0)*(v1-v0)^T
 *****************************************/
float R(point2d &v0,point2d &v1,point2d &v)
{
  float result=(v1.x-v0.x)*(v.y-v0.y)-(v1.y-v0.y)*(v.x-v0.x);
  return(result);
}


/*******************************************************************************
AUXILIARY FUNCTION OF Delauny() FUNCTION TO INCLUDE A POINT IN THE TRIANGULATION
*******************************************************************************/
void CAVIDAD(triangulation &T,int &i0,vector<int> &LisTri,vector<int> &FlatLisTri)
{
  int i,j,k,cont;
  int NLisTri=LisTri.size();
  for(i=0;i<NLisTri;i++){
    cont=0;
    if(FlatLisTri[i]==-1){
      for(j=0;j<3;j++){
        for(k=0;k<3;k++){
          if(T.t[LisTri[i0]].v[j]==T.t[LisTri[i]].v[k]){
             cont++;}
        }
      }
      if(cont==2){
        FlatLisTri[i]=0;
        CAVIDAD(T,i,LisTri,FlatLisTri);
      }
    }
  }
}

/***********************************************************************
FUNCTION TO COMPUTE THE CIRCLE PASSING BY 3 POINTS
***********************************************************************/
int C(point2d v0,point2d v1,point2d v2,float &pxc,float &pyc,float &pR2)
{
  float num1,num2,den;
  den=(v2.x-v0.x)*(v2.y-v1.y)-(v2.y-v0.y)*(v2.x-v1.x);
  if(den==0) return -1;
  den*=2;
  num1=v2.x*v2.x+v2.y*v2.y-v0.x*v0.x-v0.y*v0.y;
  num2=v2.x*v2.x+v2.y*v2.y-v1.x*v1.x-v1.y*v1.y;
  pxc=(num1*(v2.y-v1.y)-num2*(v2.y-v0.y))/den;
  pyc=(num2*(v2.x-v0.x)-num1*(v2.x-v1.x))/den;
  pR2=(v0.x-pxc)*(v0.x-pxc)+(v0.y-pyc)*(v0.y-pyc);
  return 0;

}

/**************************************************
 FUNCTION TO INCLUDE A POINT IN A TRIANGULATION
**************************************************/
int Delauny(triangulation &T,point2d v, int ValorFrontera)
{
  int i,j,k,i0;
  vector<int> LisTri,FlatLisTri;
  float xc,yc,R2;
  vector<triangle> Arista;

  //printf("punto : (%f,%f) \n",v.x,v.y);
  /* COMPROBAMOS SI EL PUNTO ESTA EN LA LISTA DE PUNTOS  */
  for(i=0;i<(int) T.p.size();i++){
    if( ((T.p[i].x-v.x)*(T.p[i].x-v.x)+(T.p[i].y-v.y)*(T.p[i].y-v.y))<1e-6){
	  //printf("punto repetido : (%f,%f) \n",T.p[i].x,T.p[i].y);
	  //printf("THE POINT IS REPEATED IN THE TRIANGULATION\n"); //system("pause");
	  return(-2);
    }
  }

  /* BUSCAMOS LOS TRIANGULOS CUYA CIRCUNFERENCIA CONTIENE AL PUNTO */
  for(i=0;i<(int) T.t.size();i++){
    if( C(T.p[T.t[i].v[0]],T.p[T.t[i].v[1]],T.p[T.t[i].v[2]],xc,yc,R2)==-1){
      printf("Hay un triangulo degenerado en la Triangulacion\n"); //system("pause");
      T.t.erase(T.t.begin()+i);
      //return(-4);
    }
    if( ((v.x-xc)*(v.x-xc)+(v.y-yc)*(v.y-yc))<=R2){
      LisTri.push_back(i);
    }
  }

   /* BUSCAMOS QUE TRIANGULO DE LA LISTA CONTIENE AL PUNTO */
   /* LO ALMACENAMOS EN i0                                 */

  i0=-1;
  //int i1=-1;
  i=0;
  while (i0==-1 && i<(int) LisTri.size()){
   float R1,R2,R3;
   R1=R(T.p[T.t[LisTri[i]].v[0]],T.p[T.t[LisTri[i]].v[1]],v);
   R2=R(T.p[T.t[LisTri[i]].v[1]],T.p[T.t[LisTri[i]].v[2]],v);
   R3=R(T.p[T.t[LisTri[i]].v[2]],T.p[T.t[LisTri[i]].v[0]],v);
   if( (R1>=-CERO && R2>=-CERO && R3>=-CERO) || (R1<=CERO && R2<=CERO && R3<=CERO) ){
     i0=i;
     break;
   }
   i++;
  }

  if(i0==-1) {
    //printf("ERROR IN DELAUNY FUNCTION: THE POINT TO INSERT IS NOT INCLUDED IN THE TRIANGULATION \n");
    //scanf("%d",&i0);
    return(-1);
  }

  /* EL TRIANGULO QUE CONTIENE AL PUNTO ES EL i0 de las lista LisTri */

  /* DETECTAMOS LA CAVIDAD POR VECINDAD */
  /* SI FlatLisTri[i]==0 EL TRIANGULO ESTA EN LA CAVIDAD POR VECINDAD */
  FlatLisTri.resize(LisTri.size());
  for(i=0;i<(int) LisTri.size();i++) FlatLisTri[i]=-1;
  FlatLisTri[i0]=0;
  CAVIDAD(T,i0,LisTri,FlatLisTri);



  /* CREAMOS LA LISTA DE AristaS */
  /* SI Arista[i][2]==0 LA Arista NO ESTA REPETIDA */
  for(i=0;i<(int) LisTri.size();i++){
   if(FlatLisTri[i]==0){
    for(j=0;j<3;j++){
      int i1=T.t[LisTri[i]].v[j];
      int i2=T.t[LisTri[i]].v[(j+1)%3];
      k=0;
      while(k<(int) Arista.size()){
        if((Arista[k].v[0]==i1 && Arista[k].v[1]==i2) ||
	    (Arista[k].v[1]==i1 && Arista[k].v[0]==i2)){
           Arista[k].v[2]=-1;
 	       k=Arista.size()+2;
         }
         k++;
      }
      if(k==(int) Arista.size()){
        Arista.push_back(triangle(i1,i2,0));
      }
    }
   }
   /*else{printf("CAVIDAD POR VECINDAD DISTINTA A CAVIDAD NORMAL\n");}*/
  }

  /* COMPROBAMOS QUE EL AREA QUE QUITAMOS ES IGUAL AL AREA QUE PONEMOS */
 {
    double SumaArea1=0;
    double SumaArea2=0;
    for(i=0;i<(int) LisTri.size();i++){
      if(FlatLisTri[i]==0){
        SumaArea1+=AREA(T.p[T.t[LisTri[i]].v[0]] ,T.p[T.t[LisTri[i]].v[1]] ,T.p[T.t[LisTri[i]].v[2]]);
      }
    }
    for(i=0;i<(int) Arista.size();i++){
      if(Arista[i].v[2]==0){
        SumaArea2+=AREA(v ,T.p[Arista[i].v[0]] ,T.p[Arista[i].v[1]]);
      }
    }
    if(fabs(SumaArea1-SumaArea2)>1e-2*fabs(SumaArea1+1e-10)){
      //printf("Suma Area Triangulos Viejos = %lf\n",SumaArea1);
      //printf("Suma Area Triangulos Nuevos = %lf\n",SumaArea2);
      //system("pause");
    }
 }
  /* ANADIMOS EL PUNTO v A LA LISTA DE NODOS DE LA TRIANGULACION */
  T.p.push_back(v);
  T.fr.push_back(ValorFrontera);

  /* ANADIMOS LOS NUEVOS TRIANGULOS */
  k=0;
  for(i=0;i<(int) Arista.size();i++){
   if(Arista[i].v[2]==0){
    if(fabs(R(T.p[T.p.size()-1],T.p[Arista[i].v[0]],T.p[Arista[i].v[1]]))>CERO ){
      while( k<(int) LisTri.size() && FlatLisTri[k]==-1)  k++;
      if(k<(int) LisTri.size()){
        T.t[LisTri[k]].v[0]=T.p.size()-1;
        T.t[LisTri[k]].v[1]=Arista[i].v[0];
        T.t[LisTri[k]].v[2]=Arista[i].v[1];
        k++;
      }
      else{
        T.t.push_back(triangle(T.p.size()-1,Arista[i].v[0],Arista[i].v[1]));
      }
    }
//    else{
//      printf("Warning: Se intento unir el punto con vertices muy alineados\n");
//      printf("Possible Causa: Introducir un nodo sobre una arista de la frontera del Objeto\n");
//    }
   }
  }
  return 0;

}

/**********************************************************************
COMPUTATION OF THE INTERCEPTION OF A LINE ax+by+c=0 AND A SEGMENT v0-v1
**********************************************************************/
point2d cross(float a,float b,float c,point2d v0,point2d v1){
  float a1=v1.y-v0.y;
  float b1=v0.x-v1.x;
  float c1=-a1*v0.x-b1*v0.y;
  float det=a*b1-b*a1;
  if(fabs(det)<CERO) return(v0);
  float x=(b*c1-c*b1)/det;
  float y=(c*a1-a*c1)/det;
  return point2d(x,y);
}

/**********************************************************************
COMPUTATION OF THE INTERCEPTION OF SEGMENTS v0-v1 and v2-v3
**********************************************************************/
point2d cross(point2d &v0,point2d &v1,point2d &v2,point2d &v3){
  float a00=v1.y-v0.y;
  float a01=v0.x-v1.x;
  float a10=v3.y-v2.y;
  float a11=v2.x-v3.x;
  float b0=a00*v0.x+a01*v0.y;
  float b1=a10*v2.x+a11*v2.y;
  float det=a00*a11-a01*a10;
  if(fabs(det)<CERO) return(v0);
  float x=(b0*a11-b1*a01)/det;
  float y=(a00*b1-a10*b0)/det;
  return point2d(x,y);
}

/**********************************************************************
CHECKING FOR THE INTERCEPTION OF SEGMENTS v0-v1 and v2-v3
**********************************************************************/
bool check_cross(point2d &v0,point2d &v1,point2d &v2,point2d &v3){
  float a00=v1.y-v0.y;
  float a01=v0.x-v1.x;
  if( (a00*(v2.x-v0.x)+a01*(v2.y-v0.y))*(a00*(v3.x-v0.x)+a01*(v3.y-v0.y))>=-1e-20) return false;
  float a10=v3.y-v2.y;
  float a11=v2.x-v3.x;
  if( (a10*(v0.x-v2.x)+a11*(v0.y-v2.y))*(a10*(v1.x-v2.x)+a11*(v1.y-v2.y))>=-1e-20) return false;
  return true;
}

/**********************************************************************
CHECKING IF v0 and v1 are in different half-planes of the line v2-v3
**********************************************************************/
bool different_half_planes(point2d &v0,point2d &v1,point2d &v2,point2d &v3){
  float a10=v3.y-v2.y;
  float a11=v2.x-v3.x;
  if( (a10*(v0.x-v2.x)+a11*(v0.y-v2.y))*(a10*(v1.x-v2.x)+a11*(v1.y-v2.y))>=-1e-20) return false;
  return true;
}

/******************************************************************
INCLUDE A COLLECTION OF SEGMENTS IN A TRIANGULATION
******************************************************************/
int Delauny(triangulation &T,vector<point2d> &p1, int ValorFrontera)
{
  if(p1.size()%2==1) return(-1); // THE NUMBER OF POINTS SHOULD BE EVEN (EACH SEGMENT IS GIVEN BY 2 CONSECUTIVE POINTS)

  triangulation T2=T; // we copy the original triangulation

  // WE INSERT IN THE TRIANGULATION THE INTERSECTION POINTS OF THE TRIANGULATION EDGES AND THE SEGMENTS
  for(int k=0;k<(int) p1.size();k+=2){ // we go through the segments
    for(int m=0;m<(int) T2.t.size();m++){ // we go through the triangles
      for(int l=0;l<3;l++){ // we go through the triangle edges.
        point2d v1=T2.p[T2.t[m].v[l]];
        point2d v2=T2.p[T2.t[m].v[(l+1)%3]];
        float R1=R(p1[k],p1[k+1],v1);
        float R2=R(p1[k],p1[k+1],v2);
        if(R1*R2>-CERO) continue; // the segment and the triangle edge do not intersect.
        float R3=R(v1,v2,p1[k]);
        float R4=R(v1,v2,p1[k+1]);
        if(R3*R4>-CERO) continue; // the segment and the triangle edge do not intersect.
        // WE COMPUTE THE INTERSECTION POINT AND INCLUDE IT IN THE TRIANGULATION
        point2d temp=cross(v1,v2,p1[k],p1[k+1]);
        Delauny(T,temp,ValorFrontera);
      }
    }
  }
  // WE CALL IN A RECURSIVE WAY TO THE FUNCTION TO ENSURE ALL SEGMENTS ARE PROPERLY INCLUDED IN THE TRIANGULATION
  if(T.p.size()>T2.p.size()){
    Delauny(T,p1,ValorFrontera);
  }
  return(0);
}


/**********************************************************************
FUNCTION TO INCLUDE A LINE IN A TRIANGULATION
**********************************************************************/
int Delauny(triangulation &T,float a,float b,float c, int ValorFrontera)
{
  // IF THE POINT ARE TOO CLOSE TO THE LINE WE REPLACE THE POINT BY ITS
  // PROJECTION ON THE LINE
  vector<float> le(T.p.size());
  for(int k=0;k<(int) T.p.size();k++){
    le[k]=a*T.p[k].x+b*T.p[k].y+c;
    if( fabs(le[k])<1e-4 && fabs(le[k])>CERO){
        T.p[k]=point2d(T.p[k].x-le[k]*a,T.p[k].y-le[k]*b);
        le[k]=0.;
    }
  }

  // WE MAKE A COPY OF THE TRIANGULATION
  triangulation T2=T;

  // WE GO THROUGH THE TRIANGLES AND WE INCLUDE NEW POINTS IN THE
  // INTERSECTION OF THE LINE WITH THE TRIANGLES
  for(int k=0;k<(int) T2.t.size();k++){
    point2d v[3]={T2.p[T2.t[k].v[0]],T2.p[T2.t[k].v[1]],T2.p[T2.t[k].v[2]]};
    float s[3]={le[T2.t[k].v[0]],le[T2.t[k].v[1]],le[T2.t[k].v[2]]};
    for(int l=0;l<3;l++){
       if(fabs(s[l])<1e-4 || fabs(s[(l+1)%3])<1e-4) continue; // points are too close to the line
       //printf("l=%d (%lf,%lf) (%lf,%lf) %lf %lf\n",l,v[l].x,v[l].y,v[(l+1)%3].x,v[(l+1)%3].y,s[l],s[(l+1)%3]);
       if(s[l]*s[(l+1)%3]<0){ //line intersects the triangle edge
         point2d temp=cross(a,b,c,v[l],v[(l+1)%3]);
         //if ( (temp-v[l])*(temp-v[(l+1)%3])<0. ){
         //printf("temp.x=%lf temp.y=%lf\n",temp.x,temp.y);
         //Delauny(T,cross(a,b,c,v[l],v[(l+1)%3]),ValorFrontera);
         Delauny(T,temp,ValorFrontera);
      }
    }
  }
  T.clean();
//  for(int k=0;k<T.p.size();k++){
//    printf("p[%d]=(%lf,%lf)\n",k,T.p[k].x,T.p[k].y);
//  }
  // WE CALL IN A RECURSIVE WAY TO THE FUNCTION TO ENSURE THE LINE IS PROPERLY  INCLUDED IN THE TRIANGULATION
  if(T2.p.size()<T.p.size()){
    Delauny(T,a,b,c,ValorFrontera);
  }
  return 0;
}


/********************************************
 Intersection with the half plane ax+by+c<0
 *******************************************/
void HalfPlaneIntersection(triangulation &Tr,float a,float b,float c){
  // WE REMOVE THE POINTS OF 1 SIDE OF THE LINE. LINE ORIENTATION IS GIVEN BY a,b.
  for(int k=0;k<(int) Tr.p.size();k++){
    if( (a*Tr.p[k].x+b*Tr.p[k].y+c)>1e-4){
      //printf("k=%d\n",k); system("pause");
      Tr.erase_point(k);
      k--;
    }
  }
}

/**********************************************************
WE BUILD A NEW TRIANGULATION BY SYMMETRY USING Nangles
LINES DIVIDING  THE CIRCLE OF CENTER C
***********************************************************/
triangulation triangulation::symmetry(point2d c,int Nangles)
{
  if(Nangles<=1) return(*this); // we check the number of angles

  triangulation T=(*this); // we make a copy of the triangulation

  //TO GENERATE A RADIAL CUT WE INCLUDE CUTTING LINES IN THE TRIANGULATION
  float alpha=0.; // we compute the associated angle to Nangles symmetries across the circle
  int Np=T.p.size()-1;
  while((int) T.p.size()>Np){
    Np=T.p.size();
    Delauny(T,cos(alpha),sin(alpha),-c.x*cos(alpha)-c.y*sin(alpha), 1);
    alpha=2.*M_PI/Nangles;
    Delauny(T,cos(-alpha),sin(-alpha),-c.x*cos(-alpha)-c.y*sin(-alpha), 1);
  }

  //WE CUT THE REGION LIMITED BY THE CUTTING LINES
  alpha=0;
  HalfPlaneIntersection(T,-cos(alpha),-sin(alpha),c.x*cos(alpha)+c.y*sin(alpha));
  //alpha=-alpha;
  alpha=-2.*M_PI/Nangles;
  HalfPlaneIntersection(T,cos(alpha),sin(alpha),-c.x*cos(alpha)-c.y*sin(alpha));

  //WE GO THROUGH THE ANGLES COMPUTE THE SYMMETRY AND WE ADD THE TRIANGULATIONS
  triangulation T2=T;
  for(int k=1;k<Nangles;k++){
    alpha=-(2*k)*M_PI/Nangles;
    T2=T2.symmetry(cos(alpha),sin(alpha),-c.x*cos(alpha)-c.y*sin(alpha));
    T=T+T2;
  }
//  T.rebuild();
//  for(int k=0;k<Nangles;k++){
//    alpha=-(2*k)*M_PI/Nangles;
//    Delauny(T,cos(alpha),sin(alpha),-c.x*cos(alpha)-c.y*sin(alpha),1);
//  }

  return(T);
}

/**********************************************************************
WE COMPUTE THE SYMMETRY WITH RESPECT TO A LINE OF A TRIANGULATION
***********************************************************************/
triangulation triangulation::symmetry(float a,float b,float c)
{
  triangulation T=(*this);
  float norm=a*a+b*b;
  for(int k=0;k<(int) p.size();k++){
    float t=2.*(a*p[k].x+b*p[k].y+c);
    t/=norm;
    T.p[k]=point2d(p[k].x-t*a,p[k].y-t*b);
  }
  return(T);
}

/****************************************************************
FUNCTION TO COMPUTE A NEW TRIANGULATION BY SYMMETRY (OLD FUNCTION)
****************************************************************/
int Symmetry(triangulation &Tr,float a,float b,float c)
{
  // WE ADD THE LINE TO THE TRIANGULATION
  Delauny(Tr,a,b,c,1);

  // WE REMOVE THE POINTS OF 1 SIDE OF THE LINE. LINE ORIENTATION IS GIVEN BY a,b.
  for(int k=0;k<(int) Tr.p.size();k++){
    if( (a*Tr.p[k].x+b*Tr.p[k].y+c)>1e-2){
      //printf("k=%d\n",k); system("pause");
      Tr.erase_point(k);
      k--;
    }
  }

  // WE CREATE NEW POINTS AND TRIANGLES BY SYMMETRY WITH RESPECT TO THE LINE
  int Np=Tr.p.size();
  vector<int> index(Tr.p.size()); // index of the new point in the symmetry part
  float norm=a*a+b*b;
  for(int k=0;k<Np;k++){
    float t=2.*(a*Tr.p[k].x+b*Tr.p[k].y+c);
    if(fabs(t)<1e-2){
      index[k]=k;
    }
    else{
      t/=norm;
      Tr.p.push_back(point2d(Tr.p[k].x-t*a,Tr.p[k].y-t*b));
      Tr.fr.push_back(1);
      index[k]=Tr.p.size()-1;
    }
  }
  int Nt=Tr.t.size();
  for(int k=0;k<Nt;k++){
    //printf("%d ",k);
    Tr.t.push_back(triangle(index[Tr.t[k].v[0]],index[Tr.t[k].v[1]],index[Tr.t[k].v[2]]));
  }
  return 0;
}

/************************************************************
FUNCTION TO CHECK IF A TRIANGULATION IS PROPERLY DEFINED
***********************************************************/
int triangulation::check()
{
  // WE CHECK IF ALL VERTEX INDEX ARE IN THE RANGE OF POINT VECTOR RANGE INDEX VALUE
  vector<int> Nv(p.size(),0);
  for(int k=0;k<(int) t.size();k++){
    for(int l=0;l<3;l++){
      if(t[k].v[l]>=(int) p.size()){
        printf("Vertex index : %d  outside range (p.size()=%d)\n",t[k].v[l],(int) p.size());
        { int ii_; scanf("%d",&ii_);} //system("pause");
        return(-1);
      }
      Nv[t[k].v[l]]++;
    }
  }

  // WE CHECK IF ALL POINTS BELONG AT LESS TO A TRIANGLE.
  for(int k=0;k<(int) p.size();k++){
    if(Nv[k]<1){
      printf("point (%lf,%lf) is not included in any vertex\n",p[k].x,p[k].y);
      system("pause");
      return(-2);
    }
  }

  // WE CHECK IF TRIANGLE AREA IS TOO SMALL
  for(int k=0;k<(int) t.size();k++){
    float area=AREA(p[t[k].v[0]] ,p[t[k].v[1]], p[t[k].v[2]]);
    if(area<CERO){
      printf("triangle (%lf,%lf) (%lf,%lf) (%lf,%lf) has a very small area\n",
             p[t[k].v[0]].x,p[t[k].v[0]].y,p[t[k].v[1]].x,p[t[k].v[1]].y,p[t[k].v[2]].x,p[t[k].v[2]].y);
      { int ii_; scanf("%d",&ii_);} //system("pause");
      return(-2);
    }
  }

  // WE CHECK IF 2 POINTS ARE TOO CLOSE
  for(int k=0;k<(int) p.size();k++){
    for(int l=k+1;l<(int) p.size();l++){
      if( (p[k]-p[l]).norm2()<CERO ){
        printf(" p[%d]=(%lf,%lf) and p[%d]=(%lf,%lf) are too close\n",k,p[k].x,p[k].y,l,p[l].x,p[l].y);
        { int ii_; scanf("%d",&ii_);} //system("pause");
        return(-3);
      }
    }
  }

  return(0);
}




/************************************************
TRIANGULATION COORDINATE NORMALIZATION. IT RETURNS
THE SCALE FACTOR APPLIED TO NORMALIZE
************************************************/
float triangulation::normalize(){
  if(p.size()<2) return(-1.);
  float xmin, xmax, ymin, ymax;

  /* we determine the border box */
  xmin = xmax = p[0].x;
  ymin = ymax = p[0].y;
  for(int i=1;i<(int) p.size();i++){
    if(xmin>p[i].x) xmin = p[i].x;
    else if(xmax<p[i].x) xmax = p[i].x;
    if(ymin>p[i].y) ymin = p[i].y;
    else if(ymax<p[i].y) ymax = p[i].y;
  }
  /* we transform the points so that the box becomes a rectangle of minimun size =1 and center 0*/
  float dx = xmax-xmin;
  float dy = ymax-ymin;
  //printf("dx=%lf dy=%lf\n",dx,dy);
  if(dx<1e-20 || dy<1e-20) return(0.);
  float scale=dx>dy?2./dx:2./dy;

  for(int i=0;i<(int) p.size();i++){
    p[i].x = (p[i].x-xmin)*scale-1.;
    p[i].y = (p[i].y-ymin)*scale-1.;
  }
  return(scale);
}

/************************************************
TRIANGULATION ROTATION TO OBTAIN THAT THE FIRST AND
LAST POINTS ARE HORIZONTAL ALIGNED
************************************************/
void triangulation::put_horizontal()
{
  if(p.size()<4) return;
  float y=p[p.size()-1].y-p[0].y+p[p.size()-2].y-p[1].y;
  float x=p[p.size()-1].x-p[0].x+p[p.size()-2].x-p[1].x;

  float a=atan2(y,x);
  float c=cos(a);
  float s=sin(a);
  for(int k=0;k<(int) p.size();k++){
    float t=p[k].x;
    p[k].x=c*p[k].x+s*p[k].y;
    p[k].y=-s*t+c*p[k].y;
  }
}



/*****************************************************************************
FUNCTION TO CLEAN UP TRIANGULATIONS REMOVING A NUMBER OF ERRORS
*****************************************************************************/
void triangulation::clean(){
  // WE LOOK FOR DUPLICATE POINTS
  for(int k=0;k<(int) p.size();k++){
    for(int m=k+1;m<(int) p.size();m++){
      if((p[k]-p[m]).norm2()<1e-6){
        p.erase (p.begin()+m);
        fr.erase (fr.begin()+m);
        for(int n=0;n<(int) t.size();n++){ // we update triangle index.
           for(int l=0;l<3;l++){
             if(t[n].v[l]>m) t[n].v[l]--;
             else if(t[n].v[l]==m) t[n].v[l]=k;
           }
        }
        m--;
      }
    }
  }

  //WE LOOK FOR POINTS NO ASSOCIATED TO TRIANGLES
  vector<int> nt(p.size(),0);
  for(int n=0;n<(int) t.size();n++){ // we update triangle index.
    for(int l=0;l<3;l++){
      nt[t[n].v[l]]++;
    }
  }

  for(int m=0;m<(int) p.size();m++){
    //printf("nt[%d]=%d\n",m,nt[m]);
    if(nt[m]==0){
      //p[m].print(); system("pause");
      p.erase (p.begin()+m);
      fr.erase (fr.begin()+m);
      nt.erase (nt.begin()+m);
      for(int n=0;n<(int) t.size();n++){ // we update triangle index.
        for(int l=0;l<3;l++){
          if(t[n].v[l]>m) t[n].v[l]--;
        }
      }
      m--;
    }
  }
  fix_triangle_orientation();
}

/**********************************************************************************
FUNCTION TO COMPUTE THE EXTREMA POINTS (xmin,ymin), (xmax,ymax)  OF A TRIANGULATION
**********************************************************************************/
// we compute the extrema (xmin,ymin), (xmax,ymax) of triangulation coordinates
vector<point2d> triangulation::extrema(){
  if(p.size()<2) return vector<point2d>();
  vector<point2d> e(2);
  e[0]=e[1]=p[0];
  for(int k=1;k<(int) p.size();k++){
    //printf("p=(%lf,%lf)\n",p[k].x,p[k].y); system("pause");
    if(e[0].x>p[k].x) e[0].x=p[k].x;
    else if(e[1].x<p[k].x) e[1].x=p[k].x;
    if(e[0].y>p[k].y) e[0].y=p[k].y;
    else if(e[1].y<p[k].y) e[1].y=p[k].y;
  }
  return(e);
}


/**********************************************************************************
FUNCTION TO CREATE A RECTANGLE TRIANGULATION FROM 2 EXTREMA POINTS
**********************************************************************************/
triangulation::triangulation(point2d pmin,point2d pmax){
  p.resize(4);
  fr.resize(4);
  t.resize(2);
  p[0]=pmin;
  p[1]=pmax;
  p[2]=point2d(pmax.x,pmin.y);
  p[3]=point2d(pmin.x,pmax.y);
  for(int k=0;k<4;k++) fr[k]=1;
  t[0]=triangle(0,1,2);  t[1]=triangle(1,3,0);

  fix_triangle_orientation();
  update_n();

}

/************************************************************************************
  constructor of the band of a rectangle with a given strip thick
************************************************************************************/
triangulation::triangulation(float strip_thick,point2d pmin,point2d pmax){
  if( (pmax.x-pmin.x)<=2*strip_thick || (pmax.y-pmin.y)<=2*strip_thick ){
    *this=triangulation(pmin,pmax);
    return;
  }
  p.resize(8);
  fr.resize(8);
  t.resize(8);
  int k=0;
  p[k++]=pmin;
  p[k++]=point2d(pmax.x,pmin.y);
  p[k++]=pmax;
  p[k++]=point2d(pmin.x,pmax.y);
  p[k++]=pmin+point2d(strip_thick,strip_thick);
  p[k++]=point2d(pmax.x,pmin.y)+point2d(-strip_thick,strip_thick);
  p[k++]=pmax+point2d(-strip_thick,-strip_thick);
  p[k++]=point2d(pmin.x,pmax.y)+point2d(strip_thick,-strip_thick);

  for(int k=0;k<8;k++) fr[k]=1;
  k=0;
  t[k++]=triangle(0,1,4);  t[k++]=triangle(1,4,5);
  t[k++]=triangle(1,5,2);
  t[k++]=triangle(5,6,2);
  t[k++]=triangle(6,2,3);
  t[k++]=triangle(7,6,3);
  t[k++]=triangle(7,3,0);
  t[k++]=triangle(7,4,0);


  fix_triangle_orientation();
  update_n();
}

/**********************************************************************************
FUNCTION TO CREATE A RECTANGLE TRIANGULATION FROM 2 EXTREMA POINTS and a minimun distance between points
**********************************************************************************/
triangulation::triangulation(point2d pmin,point2d pmax,float min_distance){
  triangulation T(pmin -(pmax-pmin),pmax+pmax-pmin);

  int Npoints=(pmax.x-pmin.x)/min_distance;
  //printf("Npoints=%d\n",Npoints);
  float h=(pmax.x-pmin.x)/Npoints;
  for(int k=0;k<=Npoints;k++){
    Delauny(T,pmin+point2d(h*k,0.),-1);
    Delauny(T,pmin+point2d(h*k,pmax.y-pmin.y),-1);
  }

  Npoints=(pmax.y-pmin.y)/min_distance;
  h=(pmax.y-pmin.y)/Npoints;
  for(int k=1;k<Npoints;k++){
    Delauny(T,pmin+point2d(0.,h*k),-1);
    Delauny(T,pmin+point2d(pmax.x-pmin.x,h*k),-1);
  }

  for(int k=0;k<4;k++) T.erase_point(0);

  T.fix_triangle_orientation();
  T.update_n();

  //system("pause");

  (*this)=T;

}

/******************************************************************************
FUNCTION TO REBUILD A TRIANGULATION
*****************************************************************************/
int triangulation::rebuild(){

  //(*this).clean(); // we first clean the triangulation

  // WE COMPUTE THE EXTREMA OF POINT COORDINATES

  vector<point2d> extrema=(*this).extrema();
  if(extrema.size()<2) return(-1);

  // WE COMPUTE AN INITIAL TRIANGULATION INCLUDING ALL POINTS
  point2d temp=extrema[1]-extrema[0];
  triangulation T(extrema[0]-temp,extrema[1]+temp);

  // WE COMPUTE PROPER INTERIOR POINTS
  vector<int> Nip(p.size(),0);
  for(int k=0;k<(int) t.size();k++){
    for(int l=0;l<3;l++){
      if(n[k].v[l]<0){ // contour edge
        Nip[t[k].v[l]]++;
        Nip[t[k].v[(l+1)%3]]++;
      }
    }
  }

  // WE INCLUDE THE POINTS in the rectangle
  for(int k=0;k<(int) p.size();k++){
    if(Nip[k]==2) Delauny(T,p[k],1);
  }
  //WE REMOVE RECTANGLE CORNERS
  for(int k=0;k<4;k++) T.erase_point(0);

  //WE REMOVE TRIANGLES WHICH ARE NOT IN THE OLD TRIANGULATION
  for(int k=0;k<(int) T.t.size();k++){
    point2d v=(T.p[T.t[k].v[0]]+T.p[T.t[k].v[1]]+T.p[T.t[k].v[2]])/3.;
    if(inclusion_test(v)==-1){
      T.t.erase(T.t.begin()+k);
      k--;
    }
  }
  (*this)=T;
  fix_triangle_orientation();
  for(int k=0;k<(int) t.size();k++){
    for(int m=0;m<3;m++){
      if( t[k].v[m]==1 && t[k].v[(m+1)%3]==0 && ( t[k].v[(m+2)%3]==2 || t[k].v[(m+2)%3]==3) ){
        triangle tT(t[k].v[(m+2)%3],1,0);
        t[k]=t[0];
        t[0]=tT;
      }
    }
  }

  //printf("triangle 0 : %d %d %d\n",t[0].v[0],t[0].v[1],t[0].v[2]);
  update_n();
  return(0);
}

/******************************************************************************
FUNCTION TO REBUILD A TRIANGULATION
*****************************************************************************/
int triangulation::rebuild(float max_segment_length){

  if(max_segment_length<1e-10) return -1;
  //(*this).clean(); // we first clean the triangulation
  rebuild();

  if(n.size()!=t.size()) update_n();

  // WE COMPUTE THE EXTREMA OF POINT COORDINATES
  vector<point2d> extrema=(*this).extrema();
  if(extrema.size()<2) return(-1);

  // WE COMPUTE AN INITIAL RECTANGLE INCLUDING ALL POINTS
  point2d temp=extrema[1]-extrema[0];
  triangulation T(extrema[0]-temp,extrema[1]+temp);

  // WE INCLUDE THE POINTS in the rectangle
  for(int k=0;k<(int) p.size();k++){
    Delauny(T,p[k],1);
  }

  // WE INCLUDE POINTS IN THE BOUNDARY SEGMENTS
  for(int k=0;k<(int) t.size();k++){
    for(int l=0;l<3;l++){
      if(n[k].v[l]==-1){
        float length=(p[t[k].v[l]]-p[t[k].v[(l+1)%3]]).norm();
        if(length>max_segment_length){
          int N=length/max_segment_length+1;
          float d=length/N;
          for(int j=1;j<N;j++){
            point2d q=p[t[k].v[l]]+(p[t[k].v[(l+1)%3]]-p[t[k].v[l]])*(j*d);
            Delauny(T,q,1);
          }
        }
      }
    }

  }

  //WE REMOVE RECTANGLE CORNERS
  for(int k=0;k<4;k++) T.erase_point(0);

  //WE REMOVE TRIANGLES WHICH ARE NOT IN THE OLD TRIANGULATION
  for(int k=0;k<(int) T.t.size();k++){
    point2d v=(T.p[T.t[k].v[0]]+T.p[T.t[k].v[1]]+T.p[T.t[k].v[2]])/3.;
    if(inclusion_test(v)==-1){
      T.t.erase(T.t.begin()+k);
      k--;
    }
  }
  (*this)=T;
  fix_triangle_orientation();
  for(int k=0;k<(int) t.size();k++){
    for(int m=0;m<3;m++){
      if( t[k].v[m]==1 && t[k].v[(m+1)%3]==0 && ( t[k].v[(m+2)%3]==2 || t[k].v[(m+2)%3]==3) ){
        triangle tT(t[k].v[(m+2)%3],1,0);
        t[k]=t[0];
        t[0]=tT;
      }
    }
  }
  update_n();
  return(0);
}

/******************************************************************************
FUNCTION TO REBUILD A TRIANGULATION USING ONLY EXTERNAL CONTOUR TRIANGLE EDGE
NEW CONTOUR TRIANGLE SHOULD HAVE A MINIMUM AREA OF minArea PARAMETER
*****************************************************************************/
int triangulation::rebuildOnlyContour(float minArea){
   minArea*=2;

   (*this).clean(); // we first clean the triangulation

  // WE COMPUTE THE EXTREMA OF POINT COORDINATES
  vector<point2d> extrema=(*this).extrema();
  if(extrema.size()<2) return(-1);

  // WE COMPUTE AN INITIAL TRIANGULATION INCLUDING ALL POINTS
  point2d temp=extrema[1]-extrema[0];
  triangulation T(extrema[0]-temp,extrema[1]+temp);

  // WE STORE THE EXTERNAL TRIANGLE EDGES IN A VECTOR
  vector<int> index;
  if(t.size() !=n.size()) update_n();
  for(int k=0;k<(int) t.size();k++){
    for(int l=0;l<3;l++){
      if(n[k].v[l]<0){
        index.push_back(t[k].v[l]);
        index.push_back(t[k].v[(l+1)%3]);
      }
    }
  }

  // WE JOINT TRIANGLE EDGES WHICH ARE IN THE SAME LINE
  int nA=index.size()+1;
  while(nA>(int) index.size()){
    //printf("index.size()=%d\n",index.size()); system("pause");
    nA=index.size();
    for(int k=0;k<(int) index.size();k+=2){
      for(int m=k+2;m<(int) index.size();m+=2){
        if(m==k) continue;
        // WE COMPARE TRIANGLE EDGE INDEX
        if(index[k]==index[m]){
          float R1=R(p[index[k]],p[index[k+1]],p[index[m+1]]);
          if(fabs(R1)<minArea){ // triangle edge points are aligned
            index[k]=index[m+1];
            index.erase(index.begin()+m);
            index.erase(index.begin()+m);
            m=0;
          }
        }
        else if(index[k]==index[m+1]){
          float R1=R(p[index[k]],p[index[k+1]],p[index[m]]);
          if(fabs(R1)<minArea){ // triangle edge points are aligned
            index[k]=index[m];
            index.erase(index.begin()+m);
            index.erase(index.begin()+m);
            m=0;
          }
        }
        else if(index[k+1]==index[m+1]){
          float R1=R(p[index[k]],p[index[k+1]],p[index[m]]);
          if(fabs(R1)<minArea){ // triangle edge points are aligned
            index[k+1]=index[m];
            index.erase(index.begin()+m);
            index.erase(index.begin()+m);
            m=0;
          }
        }
        else if(index[k+1]==index[m]){
          float R1=R(p[index[k]],p[index[k+1]],p[index[m+1]]);
          if(fabs(R1)<minArea){ // triangle edge points are aligned
            index[k+1]=index[m+1];
            index.erase(index.begin()+m);
            index.erase(index.begin()+m);
            m=0;
          }
        }
      }
    }
  }

  // WE STORE IN A VECTOR THE FINAL TRIANGLE EDGE POINTS WE ARE GOING TO INCLUDE
  // AND WE INCLUDE THE POINTS
  vector<point2d> p1;
  for(int k=0;k<(int) index.size();k+=2){
    p1.push_back(p[index[k]]);
    Delauny(T,p[index[k]],1);
    p1.push_back(p[index[k+1]]);
    Delauny(T,p[index[k+1]],1);
  }

  // WE INCLUDE THE SEGMENTS IN THE RECTANGLE
  Delauny(T,p1,1);

  //WE REMOVE RECTANGLE CORNERS
  for(int k=0;k<4;k++) T.erase_point(0);

  //WE REMOVE TRIANGLES WHICH ARE NOT IN THE OLD TRIANGULATION
  for(int k=0;k<(int) T.t.size();k++){
    point2d v=(T.p[T.t[k].v[0]]+T.p[T.t[k].v[1]]+T.p[T.t[k].v[2]])/3.;
    float R1,R2,R3;
    int m=0;
    for(;m<(int) t.size();m++){
      R1=R(p[t[m].v[0]],p[t[m].v[1]],v);
      R2=R(p[t[m].v[1]],p[t[m].v[2]],v);
      R3=R(p[t[m].v[2]],p[t[m].v[0]],v);
      if( (R1>=-1e-10 && R2>=-1e-10 && R3>=-1e-10) || (R1<=1e-10 && R2<=1e-10 && R3<=1e-10) ){ break; }
    }
    if(m==(int) t.size()){
      T.t.erase(T.t.begin()+k);
      k--;
    }
  }
  (*this)=T;
  fix_triangle_orientation();
  update_n();
  return(0);
}

/****************************************************************************
FUNCTION TO COMPUTE TRIANGLE NEIGBORHOOD INFORMATION OF A TRIANGULATION
IF A TRIANGLE EDGE HAS VALUE -1 MEANS THAT IT IS AN EXTERNAL EDGE OTHERWISE
IT RETURNS THE INDEX OF THE TRIANGLE NEIGBOR.
****************************************************************************/
void triangulation::update_n(){
  n.resize(t.size());
  for(int k=0;k<(int) t.size();k++){ for(int l=0;l<3;l++) n[k].v[l]=-1; }
  for(int k=0;k<(int) t.size();k++){
    for(int l=0;l<3;l++){
      if(n[k].v[l]>-1) continue;
      for(int m=k+1;m<(int) t.size();m++){
        for(int p=0;p<3;p++){
          if( (t[k].v[l]==t[m].v[p] && t[k].v[(l+1)%3]==t[m].v[(p+1)%3]) || (t[k].v[l]==t[m].v[(p+1)%3] && t[k].v[(l+1)%3]==t[m].v[p]) ){ // the 2 triangles share and edge
            n[k].v[l]=m;
            n[m].v[p]=k;
            m=t.size();
            break;
          }
        }
      }
    }
  }
  //UPDATE fr
  if(fr.size()!=p.size()) fr.resize(p.size());
  for(int k=0;k<(int) fr.size();k++) fr[k]=0;
  for(int k=0;k<(int) t.size();k++){
    for(int l=0;l<3;l++){
      if(n[k].v[l]>-1) continue;
      fr[t[k].v[l]]=1;
      fr[t[k].v[(l+1)%3]]=1;
    }
  }

}


//constructor of line structures
triangulation::triangulation(vector<point2d> &center_line,vector<float> &line_radius)
{
  if(center_line.size()!=line_radius.size() || center_line.size()<2){
    (*this)=triangulation();
    return;
  }

  for(int k=0;k<(int) center_line.size();k++){
    point2d pdif;
    if(k==0) pdif=center_line[1]-center_line[0];
    else pdif=center_line[k]-center_line[k-1];
    float norm=pdif.norm();
    if(norm<=0) continue;
    point2d pdifT=point2d(pdif.y,-pdif.x)*line_radius[k]/norm;
    p.push_back(center_line[k]+pdifT);
    p.push_back(center_line[k]-pdifT);
    if(k==0) continue;
    t.push_back(triangle(p.size()-1,p.size()-3,p.size()-2));
    t.push_back(triangle(p.size()-2,p.size()-4,p.size()-3));
  }
  fr=vector<int>(p.size(),1);
  fix_triangle_orientation();
  update_n();
  //clean();

}

/// RANDOM GENERATION OF A LINE TYPE TRIANGULATION BY ITERATION AND USING SPLINES
triangulation::triangulation(
triangulation &T0, /// triangulation where the new triangulation has to be included (if T0 has not triangles, the test is not performed)
point2d p0, /// point to start iterations
float dt, /// distance between consecutive points
float a, /// initial angle orientation of the curve orthogonal vector
float dal,float dar, /// interval of variation of "a" where the variation is randomly selected
float r, /// initial radius of the curve cross-section
float drl,float drr, /// interval of variation of "r" where the variation is randomly selected
int Niter) /// number of iterations to build the curve
{
  vector<point2d> pV;
  //vector<float> aV;
  vector<float> rV;
  vector<float> xV;
  vector<float> yV;

  point2d v_=point2d(sin(a),cos(a));
  point2d p_(0.,0.);
  for(int k=0;k<Niter;k++){
      if(r<0.001) r=0.001;
     //if(r*scale>0.1) r=0.1/scale;
     pV.push_back( p0+p_);
     xV.push_back(pV[pV.size()-1].x);
     yV.push_back(pV[pV.size()-1].y);
//     aV.push_back(a);
     rV.push_back(r);
     //vPr.push_back( p0+p_+point2d(v_.y,-v_.x)*r);
     float w=(float) rand()/RAND_MAX;
     float dr=w*drl+(1-w)*drr;
     r+=dr;
     w=(float) rand()/RAND_MAX;
     float da=w*dal+(1-w)*dar;
     a+=da;
     v_=point2d(sin(a),cos(a));
     p_=p_+v_*dt;
     if(T0.t.size()>0 && T0.inclusion_test(pV[pV.size()-1])==-1 ) break;
  }
//printf("pV.size()=%d\n",pV.size());

  if(pV.size()<2) return;

  // WE CREATE THE SPLINES
  int N=pV.size();
  vector<float> z(N);
  for(int k=0;k<N;k++) z[k]=k;
  float **xS,**yS,**rS; //,**aS;
  ami_calloc2d(xS,float,N,4);
  ami_calloc2d(yS,float,N,4);
//  ami_calloc2d(aS,float,N,4);
  ami_calloc2d(rS,float,N,4);
  ami_crear_splines(z,xV,xS);
  ami_crear_splines(z,yV,yS);
  ami_crear_splines(z,rV,rS);
//  ami_crear_splines(z,aV,aS);


  // WE BUILD THE TRIANGULATION
  for(float q=0.;q<=N-1+1e-6;q+=0.25){
    float xn=ami_evaluar_splines(z,xS,q);
    float yn=ami_evaluar_splines(z,yS,q);
    float rn=ami_evaluar_splines(z,rS,q);
//    float an=ami_evaluar_splines(z,aS,q);
    point2d Grad(ami_evaluar_derivada_splines(z,xS,q),ami_evaluar_derivada_splines(z,yS,q));
    Grad.normalize();
    //if(q==0.){ printf("xn=%lf,yn=%lf,rn=%lf,an=%lf,xn=%lf,yn=%lf,rn=%lf,an=%lf\n",xV[0],yV[0],rV[0],aV[0],xn,yn,rn,an); system("pause"); }
    if(rn<0.01) rn=0.01;
    //point2d disp=point2d(cos(an)*rn,-sin(an)*rn);
    point2d disp=point2d(Grad.y*rn,-Grad.x*rn);
    point2d p1=point2d(xn,yn)+disp;
    point2d p2=point2d(xn,yn)-disp;
//printf("p1=(%lf,%lf) p2=(%lf,%lf)\n",p1.x,p1.y,p2.x,p2.y);
//printf("Inclusion test %d %d \n",T0.inclusion_test(p1),T0.inclusion_test(p2));
    if(T0.t.size()>0 && (T0.inclusion_test(p1)==-1 || T0.inclusion_test(p2)==-1) ) break;
    p.push_back(p1);
    p.push_back(p2);
    if(p.size()<4) continue;
    t.push_back(triangle(p.size()-1,p.size()-3,p.size()-2));
    t.push_back(triangle(p.size()-2,p.size()-4,p.size()-3));
  }
  fr=vector<int>(p.size(),1);
  fix_triangle_orientation();
  update_n();
  clean();

  ami_free2d(xS); ami_free2d(yS); ami_free2d(rS); //ami_free2d(aS);
//  printf("p.size()=%d,t.size()=%d,fr.size()=%d,n.size()=%d\n",p.size(),t.size(),fr.size(),n.size());
}

/// LINE TYPE TRIANGULATION FROM A COLLECTION OF POINTS AND THICK RADIUS USING SPLINES
triangulation::triangulation(
triangulation &T0, /// triangulation where the new triangulation has to be included (if T0 has not triangles, the test is not performed)
vector<point2d> &pV, /// center line points of the triangulation
vector<float> thick_radius) /// cross section radius in each center line point.
{
  int N=pV.size();
  //printf("N=%d p.size()=%d thick_radius.size()=%d\n",N, p.size(),thick_radius.size());
  if( (int) thick_radius.size()!=N || pV.size()<(int) 2 ) return;
  vector<float> xV(N);
  vector<float> yV(N);
  for(int k=0;k<N;k++){
    xV[k]=pV[k].x;
    yV[k]=pV[k].y;
  }

  // WE CREATE THE SPLINES
  vector<float> z(N);
  for(int k=0;k<N;k++) z[k]=k;
  float **xS,**yS,**rS;
  ami_calloc2d(xS,float,N,4);
  ami_calloc2d(yS,float,N,4);
  ami_calloc2d(rS,float,N,4);
  ami_crear_splines(z,xV,xS);
  ami_crear_splines(z,yV,yS);
  ami_crear_splines(z,thick_radius,rS);

  // WE BUILD THE TRIANGULATION
  for(float q=0.;q<=N-1+1e-6;q+=0.25){
    float xn=ami_evaluar_splines(z,xS,q);
    float yn=ami_evaluar_splines(z,yS,q);
    float rn=ami_evaluar_splines(z,rS,q);
    point2d Grad(ami_evaluar_derivada_splines(z,xS,q),ami_evaluar_derivada_splines(z,yS,q));
    Grad.normalize();
    if(rn<0.01) rn=0.01;
    point2d disp=point2d(Grad.y*rn,-Grad.x*rn);
    point2d p1=point2d(xn,yn)+disp;
    point2d p2=point2d(xn,yn)-disp;
    if(T0.t.size()>0 && (T0.inclusion_test(p1)==-1 || T0.inclusion_test(p2)==-1) ) break;
    p.push_back(p1);
    p.push_back(p2);
    if(p.size()<4) continue;
    t.push_back(triangle(p.size()-1,p.size()-3,p.size()-2));
    t.push_back(triangle(p.size()-2,p.size()-4,p.size()-3));
  }
  //printf("N=%d p.size()=%d \n",N, p.size());
  fr=vector<int>(p.size(),1);
  clean();
  fix_triangle_orientation();
  update_n();

  ami_free2d(xS); ami_free2d(yS); ami_free2d(rS); //ami_free2d(aS);
}




///constructor from a polygon
/// points have to be sorted in a clockwise way (first and last points are not equal)
triangulation::triangulation(vector<point2d> &pol)
{
  //printf("pol.size()=%d\n",pol.size());
  if(pol.size()<3) return;

  // WE COMPUTE A RECTANGLE INCLUDING ALL POINTS
  point2d pmin=pol[0];
  point2d pmax=pol[0];
  for(int k=1;k<(int) pol.size();k++){
    if(pol[k].x>pmax.x) pmax.x=pol[k].x;
    else if (pol[k].x<pmin.x) pmin.x=pol[k].x;
    if(pol[k].y>pmax.y) pmax.y=pol[k].y;
    else if (pol[k].y<pmin.y) pmin.y=pol[k].y;
  }
  float dx=(pmax.x-pmin.x)/2+1.;
  float dy=(pmax.y-pmin.y)/2+1.;
  pmin.x-=dx; pmax.x+=dx; pmin.y-=dy; pmax.y+=dy;
  triangulation T,T0(pmin,pmax);

  /// WE INCLUDE ALL POLYGON POINTS IN THE RECTANGLE

  int Nt=0;
  while(Nt++<100){
    T=T0;
    for(unsigned int k=0;k<pol.size();k++){
//printf("pol[%d]=(%lf,%lf)\n",k,pol[k].x,pol[k].y);
      if(Delauny(T,pol[k],1)<0){
        pol.erase(pol.begin()+k);
        k--;
      }
    }
    if( T.p.size()-4!=pol.size() ){
      printf("problems in triangulation::triangulation(vector<point2d> &pol)\n");
      printf("T.t.size()=%d, T.p.size()=%d, pol.size()=%d\n",(int) T.t.size(),(int) T.p.size(),(int) pol.size());
      return;
    }

    /// WE CHECK THAT ALL POLYGON SIDES ARE TRIANGLE EDGE SEGMENTS
    //printf("T.p.size()=%d\n",T.p.size());
    vector<bool> flat(T.p.size(),false);
    for(int k=0;k<(int) T.t.size();k++){
      for(int n=0;n<3;n++){
        int i0=T.t[k].v[n];
        int i1=T.t[k].v[(n+1)%3];
        if(i0<4 || i1<4 ) continue;
        if(i0>i1){i0=i1; i1=T.t[k].v[n];}
        if(i1==i0+1) flat[i0]=true;
        else if(i0==4 && i1==(int) T.p.size()-1) flat[T.p.size()-1]=true;
      }
    }
//for(unsigned int k=4;k<flat.size();k++) {printf("flat[%d]=%d\n",k,(int) flat[k]);}
    unsigned int m=0;
    for(m=4;m<flat.size();m++){
      if(flat[m]==false){
        if( (m-3)<pol.size()) pol.insert(pol.begin()+m-3,(pol[m-4]+pol[m-3])*0.5);
        else pol.push_back((pol[pol.size()-1]+pol[0])*0.5);
        break;
      }
    }
    //printf("m=%d flat.size()=%d\n",m,flat.size()); system("pause");
    if(m==flat.size()) break;
  }


  // WE ERASE THE RECTANGLE CORNERS AND UPDATE TRIANGLE NEIGHBORHOOD INFORMATION
  for(int k=0;k<4;k++) T.erase_point(0);
  T.update_n();
  //(*this)=T; return;

  // WE ERASE TRIANGLES WHICH ARE OUTSIDE OF THE FIGURE ACCORDING TO CLOCKWISE ORIENTATION
  // auxiliary vector to store information about deleted triangles
  vector <int> DEL(T.t.size(),0); // =0(initialization) =1(triangle is keeped) =2(triangle to be erased)

  // WE GO TRHOUGH THE TRIANGLES
  for(int m=0;m<(int) T.t.size();m++){
    // WE SORT THE TRIANGLE INDEXES
    if(T.t[m].v[0]>T.t[m].v[1]){int temp=T.t[m].v[0]; T.t[m].v[0]=T.t[m].v[1]; T.t[m].v[1]=temp; }
    if(T.t[m].v[0]>T.t[m].v[2]){int temp=T.t[m].v[0]; T.t[m].v[0]=T.t[m].v[2]; T.t[m].v[2]=temp; }
    if(T.t[m].v[2]<T.t[m].v[1]){int temp=T.t[m].v[2]; T.t[m].v[2]=T.t[m].v[1]; T.t[m].v[1]=temp; }
    // WE CHECK IF THERE ARE TO CONSECUTIVE POLYGON POINTS IN THE TRIANGLE
    for(int q=0;q<3;q++){
      if(T.t[m].v[(q+1)%3]==T.t[m].v[q]+1 || (T.t[m].v[q]==(int) pol.size()-1 && T.t[m].v[(q+1)%3]==0 ) ){
        DEL[m]=1; // two consecutive polygon points found
        if( R(T.p[T.t[m].v[q]],T.p[T.t[m].v[(q+1)%3]],T.p[T.t[m].v[(q+2)%3]])>0 ){
          DEL[m]=2; // triangle to be erased because is not clockwise oriented
        }
        break;
      }
    }
  }

  // WE MARKED TO BE ERASED TRIANGLES WITH DEL==0 AND WITH A NEIGBOR TRIANGLE TO BE ERASED
  bool del=true;
  while(del==true){ // each time we marked to be erased a new triangle we go again through all the list of triangles to look for more triangles to be erased
    del=false;
    for(int m=0;m<(int) T.t.size();m++){
      if(DEL[m]==0){
        for(int q=0;q<3;q++){
          if(T.n[m].v[q]>=0 && DEL[T.n[m].v[q]]==2){ del=true; DEL[m]=2;}
        }
      }
    }
  }

  // WE ERASE ALL MARKED TRIANGLES AND UPDATE TRIANGLE NEIGHBORHOOD INFORMATION
  for(int m=T.t.size()-1;m>=0;m--){
    if(DEL[m]==2) T.t.erase(T.t.begin()+m);
  }
  T.fix_triangle_orientation();
  T.update_n();

  (*this)=T;
}


float beta(float alpha){
  int t1=alpha/(2.*M_PI);
  float t2=alpha>0.?alpha-t1*2*M_PI:alpha-(t1-1)*2*M_PI;
  t2/=(M_PI/2.);
  if(t2<1.){
    float t=1./5.;
    float gx=t2<t?0.5*t2/t:0.5+0.5*(t2-t)/(1-t);
    float f=4*gx*(1-gx);
    //printf("1: alpha=%lf t2=%lf f=%lf\n",alpha,t2,f);
    return(f);
  }
  else if(t2<2.){
    t2-=1.;
    float t=2./5.;
    float gx=t2<t?0.5*t2/t:0.5+0.5*(t2-t)/(1-t);
    float f=4*gx*(1-gx);
    //printf("2: alpha=%lf t2=%lf f=%lf\n",alpha,t2,f);
    return(f);
  }
  else if(t2<3.){
    t2-=2.;
    float t=3./5.;
    float gx=t2<t?0.5*t2/t:0.5+0.5*(t2-t)/(1-t);
    float f=4*gx*(1-gx);
    //printf("3: alpha=%lf t2=%lf f=%lf\n",alpha,t2,f);
    return(f);
  }
  else{
    t2-=3.;
    float t=4./5.;
    float gx=t2<t?0.5*t2/t:0.5+0.5*(t2-t)/(1-t);
    float f=4*gx*(1-gx);
    //printf("4: alpha=%lf t2=%lf f=%lf\n",alpha,t2,f);
    return(f);
  }

}
/*********************************************************************
DEFORMATION USING OVOID AND TANGENCIAL COMPONENT
******************************************************************/
triangulation triangulation::deformation(point2d center,float angle0,float par0,float angle1,float par1)
{
  //for(float alfa=0;alfa<3.15*2;alfa+=0.05){printf("alfa=%lf beta=%lf\n",180*alfa/3.1416,beta(alfa));}  system("pause");

  triangulation T=(*this);

//for(int k=0;k<=64;k++){ printf("angle : %lf  beta= %lf \n",180.*k*2./64.,beta(k*3.1416*2./64.)); } system("pause");
//par1=0;
  for(int k=0;k<(int) T.p.size();k++) {
    //printf("(%lf,%lf) -> ",T.p[k].x,T.p[k].y);
    point2d temp=T.p[k]-center;
    float rho=temp.norm();
    float alpha=atan2(temp.y,temp.x);
    //float alpha1=alpha+par0*sin(2.*(alpha-angle0));
    float alpha2=par0*beta(alpha-angle0);
    float alpha1=alpha+alpha2;
    float co=cos(alpha1);//*(1-par1*fabs(cos(0.5*(alpha-angle1))*cos(0.5*(alpha-angle1))));
    //float si=sin(alpha1)*(1-par1*fabs(cos(0.5*(alpha-angle1))*cos(0.5*(alpha-angle1))));
    float da=alpha-angle1;
    float temp2=(cos(0.5*da));
    temp2*=1./(1.+0.2*fabs(da));
    float si=sin(alpha1)*(1-par1*temp2);
    //float si=sin(alpha1)*(1-par1*beta(alpha-angle1));
    //float rho2=rho*sqrt(co*co+si*si);
    //printf("(%lf,%lf)  ",T.p[k].x,T.p[k].y);
    float co2=cos(alpha2);
    float si2=sin(alpha2);
    float x=rho*co;
    float y=rho*si;
    T.p[k]=center+point2d(x*co2+y*si2,-x*si2+y*co2);
    //printf("(%lf,%lf) si=%lf \n",T.p[k].x,T.p[k].y,si);
  }
  //system("pause");
  return(T);
}

/*********************************************************************
DEFORMATION USING OVOID AND TANGENCIAL COMPONENT
******************************************************************/
triangulation triangulation::deformation3(point2d center,float angle0,float par0,float angle1,float par1)
{
  //for(float alfa=0;alfa<3.15*2;alfa+=0.05){printf("alfa=%lf beta=%lf\n",180*alfa/3.1416,beta(alfa));}  system("pause");

  triangulation T=(*this);

//for(int k=0;k<=64;k++){ printf("angle : %lf  beta= %lf \n",180.*k*2./64.,beta(k*3.1416*2./64.)); } system("pause");
//par1=0;
  for(int k=0;k<(int) T.p.size();k++) {
    //printf("(%lf,%lf) -> ",T.p[k].x,T.p[k].y);
    point2d temp=T.p[k]-center;
    float alpha=atan2(temp.y,temp.x);
    //float alpha1=alpha+par0*sin(2.*(alpha-angle0));
    float alpha1=angle0+par0*beta(alpha);
    float co=cos(alpha1);
    float si=sin(alpha1);

    //*(1-par1*fabs(cos(0.5*(alpha-angle1))*cos(0.5*(alpha-angle1))));
    //float si=sin(alpha1)*(1-par1*fabs(cos(0.5*(alpha-angle1))*cos(0.5*(alpha-angle1))));
    float da=alpha-angle1;
//    if(k==0) da0=da;
//    if(fabs(da-da0)>fabs(da+2*M_PI-da0)) da=da0-2*M_PI;
//    else if(fabs(da-da0)>fabs(da-2*M_PI-da0)) da=da0+2*M_PI;
//    da0=da;
    float temp2=cos(da);
    temp2*=1./(1.+0.3*fabs(sin(alpha)*(da)));
    float x=co*temp.x+si*temp.y;
    float y=(-si*temp.x+co*temp.y)*(1-par1*temp2);
    //float si=sin(alpha1)*(1-par1*beta(alpha-angle1));
    //float rho2=rho*sqrt(co*co+si*si);
    //printf("(%lf,%lf)  ",T.p[k].x,T.p[k].y);
    T.p[k]=center+point2d(x*co-y*si,x*si+y*co);
    //printf("(%lf,%lf) si=%lf \n",T.p[k].x,T.p[k].y,si);
  }
  //system("pause");
  return(T);
}

triangulation triangulation::deformation2(point2d center,float angle0,float par0,float angle1,float par1)
{
  triangulation T=(*this);

  // WE COMPUTE MAXIMUM DISTANCE OF THE OBJECT POINTS TO THE CENTER
  float rho_max=0;
   for(int k=0;k<(int) T.p.size();k++) {
    float rho=(T.p[k]-center).norm();
    if(rho>rho_max) rho_max=rho;
  }

  for(int k=0;k<(int) T.p.size();k++) {
    //printf("(%lf,%lf) -> ",T.p[k].x,T.p[k].y);
    point2d temp=T.p[k]-center;
    float rho=temp.norm();
    float alpha=atan2(temp.y,temp.x);
    //float alpha1=alpha+par0*sin(2.*(alpha-angle0));
    float alpha2=par0*beta(alpha-angle0);
    float alpha1=alpha+alpha2;
    float co=cos(alpha1);//*(1-par1*fabs(cos(0.5*(alpha-angle1))*cos(0.5*(alpha-angle1))));
    //float si=sin(alpha1)*(1-par1*fabs(cos(0.5*(alpha-angle1))*cos(0.5*(alpha-angle1))));
    float da=alpha-angle1;
    float temp2=(cos(0.5*da));
    temp2*=1./(1.+0.2*fabs(da));
    float si=sin(alpha1)*(1-par1*temp2);
    //float si=sin(alpha1)*(1-par1*beta(alpha-angle1));
    //float rho2=rho*sqrt(co*co+si*si);
    //printf("(%lf,%lf)  ",T.p[k].x,T.p[k].y);
    float co2=cos(alpha2);
    float si2=sin(alpha2);
    //float x=1.5*rho_max*rho*co/(rho_max+0.5*rho);
    //float y=1.5*rho_max*rho*si/(rho_max+0.5*rho);
    float x=1.2*rho_max*rho*co/(rho_max+0.2*rho);
    float y=1.2*rho_max*rho*si/(rho_max+0.2*rho);
    T.p[k]=center+point2d(x*co2+y*si2,-x*si2+y*co2);
    //printf("(%lf,%lf) si=%lf \n",T.p[k].x,T.p[k].y,si);
  }
  //system("pause");
  return(T);
}

/**********************************************************
test to check if a point is inside a triangulation.
it returns the index of the triangle containing de point
or -1 otherwise
**********************************************************/
int triangulation::inclusion_test(point2d v)
{
  for(int i=0;i<(int) t.size();i++){
    float R1=R(p[t[i].v[0]],p[t[i].v[1]],v);
    float R2=R(p[t[i].v[1]],p[t[i].v[2]],v);
    float R3=R(p[t[i].v[2]],p[t[i].v[0]],v);
    if( (R1>=0 && R2>=0 && R3>=0) || (R1<=0 && R2<=0 && R3<=0) ) return(i);
  }
  return(-1);
}

/**********************************************************
test to check if a triangulation is inside a triangulation
**********************************************************/
bool triangulation::inclusion_test(triangulation &T)
{
  // WE CHECK IF ALL THE POINTS OF T ARE INSIDE THE CURRENT TRIANGULATION
  for(int i=0;i<(int) T.p.size();i++){
    if(inclusion_test(T.p[i])<0) return(false);
  }
  //return true; //AAAA
  // WE COMPUTE TRIANGLE NEIGHBORHOOD INFORMATION
  if(t.size() !=n.size()) update_n();

  // WE CHECK IF THE SEGMENTS OF T INTERCEPTS BOUNDARY SEGMENTS OF THE CURRENT TRIANGULATION n.resize(t.size());
  for(int k=0;k<(int) t.size();k++){
    for(int l=0;l<3;l++){
      if(n[k].v[l]>-1) continue;
      for(int m=0;m<(int) T.t.size();m++){
        for(int q=0;q<3;q++){
          if(check_cross( p[t[k].v[l]],p[t[k].v[(l+1)%3]],T.p[T.t[m].v[q]],T.p[T.t[m].v[(q+1)%3]] ) == true ){
//            p[t[k].v[l]].print();
//            p[t[k].v[(l+1)%3]].print();
//            T.p[T.t[m].v[q]].print();
//            T.p[T.t[m].v[(q+1)%3]].print();
//            system("pause");
            return(false);
          }
          //if(check_cross( p[0],p[0],p[0],p[0] ) == true ) return(false);
        }
      }
    }
  }
  return(true);
}

/**********************************************************
test to check if a triangulation intersects  triangulation
**********************************************************/
bool triangulation::intersection_test(triangulation &T)
{
  // WE CHECK IF ANY OF THE POINTS OF T ARE INSIDE THE CURRENT TRIANGULATION
  for(int i=0;i<(int) T.p.size();i++){
    if(inclusion_test(T.p[i])>=0) return(true);
  }

  // WE CHECK IF ANY OF THE POINTS OF THE CURRENT TRIANGULATION ARE INSIDE  T
  for(int i=0;i<(int) p.size();i++){
    if(T.inclusion_test(p[i])>=0) return(true);
  }


  // WE COMPUTE TRIANGLE NEIGHBORHOOD INFORMATION
  if(t.size() !=n.size()) update_n();

  // WE CHECK IF THE SEGMENTS OF T INTERCEPTS BOUNDARY SEGMENTS OF THE CURRENT TRIANGULATION n.resize(t.size());
  for(int k=0;k<(int) t.size();k++){
    for(int l=0;l<3;l++){
      if(n[k].v[l]>-1) continue;
      for(int m=0;m<(int) T.t.size();m++){
        for(int q=0;q<3;q++){
          if(check_cross( p[t[k].v[l]],p[t[k].v[(l+1)%3]],T.p[T.t[m].v[q]],T.p[T.t[m].v[(q+1)%3]] ) == true ) return(true);
          //if(check_cross( p[0],p[0],p[0],p[0] ) == true ) return(false);
        }
      }
    }
  }
  return(false);
}


/** affinity estimation from triangle segments. segments will be oriented at a given distance and orientation
*/
affinity triangles2affinity(
triangulation &T1 /** first triangulation */,
int tindex1 /** triangle index in the first triangulation*/,
int index1 /** segment index in the first triangle*/,
triangulation &T2 /** second triangulation */,
int tindex2 /** triangle index in the second triangulation */,
int index2 /** segment index in the second triangle*/,
float distance /** distance between segment midpoints*/,
float angle /** segment orientation */)
{
  if(tindex1<0 || tindex1>=(int) T1.t.size() || index1<0 || index1>2) return affinity();
  if(tindex2<0 || tindex2>=(int) T2.t.size() || index2<0 || index2>2) return affinity();

  //point coordinates
  point2d p10=T1.p[T1.t[tindex1].v[index1]];
  point2d p11=T1.p[T1.t[tindex1].v[(index1+1)%3]];
  point2d p20=T2.p[T2.t[tindex2].v[index2]];
  point2d p21=T2.p[T2.t[tindex2].v[(index2+1)%3]];

  //segments mid points
  point2d midpoint1=(p10+p11)/2.;
  point2d midpoint2=(p20+p21)/2.;

  //segments normal orientation
  point2d p1n=point2d(p11.y-p10.y,p10.x-p11.x)/(p11-p10).norm();
  point2d p2n=point2d(p21.y-p20.y,p20.x-p21.x)/(p21-p20).norm();
  if(p1n*(midpoint1-T1.p[T1.t[tindex1].v[(index1+2)%3]])<0) p1n=p1n*(-1.);
  if(p2n*(midpoint2-T2.p[T2.t[tindex2].v[(index2+2)%3]])<0) p2n=p2n*(-1.);

//  printf("midpoint1=(%lf,%lf)\n",midpoint1.x,midpoint1.y);
//  printf("p1n=(%lf,%lf)\n",p1n.x,p1n.y);
//  printf("midpoint2=(%lf,%lf)\n",midpoint2.x,midpoint2.y);
//  printf("p2n=(%lf,%lf)\n",p2n.x,p2n.y);

  //we update p1n using parameter angle
  float co=cos(-angle);
  float si=sin(-angle);
  float temp=p1n.x;
  p1n.x=co*p1n.x+si*p1n.y;
  p1n.y=-si*temp+co*p1n.y;

  //angle rotation estimation
  float angle2=acos(-1*(p1n*p2n));
  co=cos(-angle2);
  si=sin(-angle2);
  if( ((co*p2n.x+si*p2n.y)*p1n.x+(-si*p2n.x+co*p2n.y)*p1n.y) > ((co*p2n.x-si*p2n.y)*p1n.x+(si*p2n.x+co*p2n.y)*p1n.y) ) si=-si;

  //end point for the second mid point
  point2d endpoint=midpoint1+p1n*distance;

  //final affinity estimation
  affinity A;
  A.A[0][0]=co; A.A[0][1]=si; A.A[1][0]=-si; A.A[1][1]=co;
  A.cx=endpoint.x-(co*midpoint2.x+si*midpoint2.y);
  A.cy=endpoint.y-(-si*midpoint2.x+co*midpoint2.y);

  // WE CHECK THE RESULT
  point2d c=A*midpoint2;
  //printf("c=(%lf,%lf)\n",c.x,c.y);

  //system("pause");

  return A;

}

void neighborhood_build(triangulation &T,int seed,vector<int> &ConnComp){
  //printf("seed=%d\n",seed);
  //printf("T.t.size()=%d,T.n.size()=%d\n",T.t.size(),T.n.size());
  if(seed<0 || seed>=(int) T.t.size()) return;
  for(int k=0;k<(int) ConnComp.size();k++) if(seed==ConnComp[k]) return;
  //printf("ok1\n");
  ConnComp.push_back(seed);
  //printf("ok2\n");
  for(int k=0;k<3;k++) neighborhood_build(T,T.n[seed].v[k],ConnComp);
  return;
}

vector<triangulation> triangulation::connected_components(){
  vector<triangulation> vT;
  vector<int> flat(t.size(),-1);

  for(int m=0,k=0;k<(int) t.size();k++){
    //printf("m=%d\n",m);
    //for(int l=0;l<t.size();l++) printf("%d ",flat[l]); printf("\n");
    if(flat[k]==-1){
      vector<int> ConnComp;
      neighborhood_build((*this),k,ConnComp);
      for(int j=0;j<(int) ConnComp.size();j++) flat[ConnComp[j]]=m;
      m++;
    }

  }
  vector<int> flatp(p.size(),-1);
  for(int k=0;k<(int) t.size();k++) for(int j=0;j<3;j++) flatp[t[k].v[j]]=flat[k];

  vector<int> new_index(p.size(),-1);
  for(int k=0;k<(int) p.size();k++){
     if(flatp[k]>-1){
       triangulation T2;
       int n=0;
       new_index[k]=0;
       T2.p.push_back(p[k]);
       for(int j=k+1;j<(int) p.size();j++){
         if(flatp[j]==flatp[k]){
           new_index[j]=n++;
           T2.p.push_back(p[j]);
           flatp[j]=-1;
         }
       }
       for(int i=0;i<(int) t.size();i++){
         if(flat[i]==flatp[k]){
           T2.t.push_back(triangle(new_index[t[i].v[0]],new_index[t[i].v[1]],new_index[t[i].v[2]]));
         }
       }
       T2.fr.resize(T2.p.size());
       T2.fix_triangle_orientation();
       T2.update_n();
       vT.push_back(T2);
     }
  }
  return(vT);
}

/// COMPUTE THE INTERSECTION OF A SEGMENT AND A TRIANGULATION BORDER. COMPUTE (a,b), WHERE
/// (px,py)+a*(nx,ny)=(qx,qy)+b*(mx,my). p,n DETERMINES THE SEGMENT AND q,m MOVE ALONG
/// TRIANGULATION BORDER SEGMENTs. b is in [0,1] and a>=0 is as smaller as possible
/// IF THE METHOD FAIL RETURN -1, OTHERWISE IT RETURNS a
float segment_triangulation(point2d p, point2d n, triangulation &T){
  if(T.t.size()!=T.n.size()) T.update_n();
  float a=1e30,a0=1e30,b=-1;
  for(int k=0;k<(int) T.t.size();k++){
    for(int l=0;l<3;l++){
      if(T.n[k].v[l]<0) continue;
      point2d q=T.p[T.t[k].v[l]];
      point2d m=T.p[T.t[k].v[(l+1)%3]]-T.p[T.t[k].v[l]];
      float det=m.y*n.x-m.x*n.y;
      if(fabs(det)<1e-20) continue;
      b=(n.x*p.y-n.y*p.x-n.x*q.y+n.y*q.x)/det;
      if(b<0 || b>1) continue;
      float temp=(m.x*p.y-m.y*p.x-m.x*q.y+m.y*q.x)/det;
      if(temp<=1e-10) continue;
      if(temp<a) a=temp;
    }
  }
  if(a==a0) return -1.;
  return a;
}

triangulation triangulation::border(float thick) // build the external border of a triangulation
{

  if(n.size()!=t.size()) update_n();

  //printf("thick=%lf\n",thick); system("pause");

  triangulation T;

  vector< vector < point2d > > pb(p.size());
  vector< vector < int > > pv(p.size());
  for(int k=0;k<(int) t.size();k++){
    for(int l=0;l<3;l++){
      if(n[k].v[l]<0){ // contour edge
        point2d n1=(p[t[k].v[l]]-p[t[k].v[(l+1)%3]]).normal();
        point2d q=p[t[k].v[l]]+n1*thick;
        if( R(p[t[k].v[l]],p[t[k].v[(l+1)%3]],q)*R(p[t[k].v[l]],p[t[k].v[(l+1)%3]],p[t[k].v[(l+2)%3]]) < 0 ){
        //if(inclusion_test(q)<0){
        //if( (p[t[k].v[l]]+n1-p[t[k].v[(l+2)%3]]).norm2()>(p[t[k].v[l]]-n1-p[t[k].v[(l+2)%3]]).norm2()){
          n1.x=-n1.x; n1.y=-n1.y;
        }
        pb[t[k].v[l]].push_back(p[t[k].v[l]]+n1*fabs(thick));
        pb[t[k].v[l]].push_back(p[t[k].v[(l+1)%3]]+n1*fabs(thick));
        pb[t[k].v[(l+1)%3]].push_back(p[t[k].v[(l+1)%3]]+n1*fabs(thick));
        pb[t[k].v[(l+1)%3]].push_back(p[t[k].v[l]]+n1*fabs(thick));
      }
    }
  }

  vector < point2d > p2(p.size());
  for(int k=0;k<(int) p.size();k++){
    //printf("pb[%d].size()=%d\n",k,pb[k].size()); system("pause");
    if((int) pb[k].size()!=4){
      p2[k]=p[k];
      //printf("thick=%lf\n",thick);
      continue;
    }
    point2d v0=pb[k][1]-pb[k][0];
    point2d v1=pb[k][3]-pb[k][2];
    float ps=v0*v1;
    ps=(ps*ps)/(v0.norm2()*v1.norm2());
    if(ps<0.9) p2[k]=cross(pb[k][0],pb[k][1],pb[k][2],pb[k][3]);
    else p2[k]=pb[k][0];
    int inT=inclusion_test(p2[k]);
    if( (inT<0 && thick>0) || (inT>=0 && thick<0) ){
      float a=segment_triangulation(p[k],p2[k]-p[k],*this);
      //printf("a=%lf",a);
      //system("pause");
      if(a>0 && a<=1){
        inT=inclusion_test(p[k]+(p2[k]-p[k])*a*0.9);
        if( (inT>=0 && thick>0) || (inT<0 && thick<0) ){
          p2[k]=p[k]+(p2[k]-p[k])*a;
          continue;
        }
      }
      inT=inclusion_test(pb[k][0]);
      if( (inT>=0 && thick>0) || (inT<0 && thick<0) ){
        p2[k]=pb[k][0];
        continue;
      }
      inT=inclusion_test(pb[k][2]);
      if( (inT>=0 && thick>0) || (inT<0 && thick<0) ){
        p2[k]=pb[k][2];
        continue;
      }
      for(a=0.99;a>0.02;a=a-0.01){
        inT=inclusion_test(p[k]+(p2[k]-p[k])*a);
        if( (inT>=0 && thick>0) || (inT<0 && thick<0) ){
          p2[k]=p[k]+(p2[k]-p[k])*a;
          a=0.;
          break;
        }
      }
      if(a>0){
        //printf("a=%lf thick=%lf\n",a,thick);
        //system("pause");
        p2[k]=p[k];
      }
    }
  }

//   for(int k=0;k<p.size();k++){
//    if(pb[k].size()==4) continue;
//    int imin;
//    float dmin=1e20;
//    for(int i=0;i<p.size();i++){
//      if(i==k && pb[i].size()!=4) continue;
//      float norm2=(p[k]-p[i]).norm2();
//      if(norm2<dmin){
//         dmin=norm2;
//         imin=i;
//      }
//    }
//    p2[k]=p2[imin];
//   }


  vector< int > index(p.size(),-1);

  for(int k=0;k<(int) t.size();k++){
    for(int l=0;l<3;l++){
      if(n[k].v[l]<0){ // contour edge
        if(index[t[k].v[l]]==-1){
          index[t[k].v[l]]=T.p.size();
          point2d np=p[t[k].v[l]]-p2[t[k].v[l]];
          np.normalize();
          T.p.push_back(p[t[k].v[l]]+np*0.0001);
          T.p.push_back(p2[t[k].v[l]]);
          T.fr.push_back(1);
          T.fr.push_back(1);
        }
        if(index[t[k].v[(l+1)%3]]==-1){
          index[t[k].v[(l+1)%3]]=T.p.size();
          point2d np=p[t[k].v[(l+1)%3]]-p2[t[k].v[(l+1)%3]];
          np.normalize();
          T.p.push_back(p[t[k].v[(l+1)%3]]+np*0.0001);
          T.p.push_back(p2[t[k].v[(l+1)%3]]);
          T.fr.push_back(1);
          T.fr.push_back(1);
        }
        T.t.push_back(triangle(index[t[k].v[l]],index[t[k].v[(l+1)%3]],index[t[k].v[l]]+1));
        T.t.push_back(triangle(index[t[k].v[l]]+1,index[t[k].v[(l+1)%3]],index[t[k].v[(l+1)%3]]+1));
      }
    }
  }
  T.fix_triangle_orientation();
  T.update_n();
  return T;

}


/// DIVISION OF A TRIANGULATION IN TWO CUTTING BY A LINE
vector<triangulation> triangulation::DivisionByLine(float a,float b,float c){
  vector<triangulation> Tresult;
  triangulation Tr1=(*this);
  Delauny(Tr1,a,b,c,0);
//  for(int k=0;k<Tr1.p.size();k++){
//    printf("p[%d]=(%lf,%lf)\n",k,Tr1.p[k].x,Tr1.p[k].y);
//  }
  triangulation Tr2=Tr1;
  for(int k=0;k<(int) Tr1.p.size();k++){
    //printf(" k=%d Tr1.p.size()=%d\n",k,Tr1.p.size());
    if( (a*Tr1.p[k].x+b*Tr1.p[k].y+c)>1e-10 ){
      Tr1.erase_point(k);
      k--;
    }
  }
  Tr1.update_n();
  for(int k=0;k<(int) Tr2.p.size();k++){
    //printf(" k=%d Tr2.p.size()=%d\n",k,Tr2.p.size());
    if( (a*Tr2.p[k].x+b*Tr2.p[k].y+c)<-1e-10 ){
      Tr2.erase_point(k);
      k--;
    }
  }
  Tr2.update_n();
  Tresult.push_back(Tr1);
  Tresult.push_back(Tr2);
  return(Tresult);
}

double triangulation::area(){
  double a=0.;
  for(int i=0;i<(int) t.size();i++)
    a+=AREA(p[t[i].v[0]] ,p[t[i].v[1]] ,p[t[i].v[2]]);

  return a;
}


triangulation triangulation::division(float max_length){
  ///printf("init division()\n");
  // WE COMPUTE THE EXTREMA OF POINT COORDINATES
  vector<point2d> extrema=(*this).extrema();
  if(extrema.size()<2) return triangulation();

  //update_n();

  // WE COMPUTE AN INITIAL TRIANGULATION INCLUDING ALL POINTS
  point2d temp=extrema[1]-extrema[0];
  triangulation T(extrema[0]-temp,extrema[1]+temp);

  ///printf("WE INCLUDE THE POINTS in the rectangle\n");
  for(int k=0;k<(int) p.size();k++){
    Delauny(T,p[k],0);
  }

  update_n();

  ///printf("WE INCLUDE NEW POINTS BY SEGMENT DIVISION IN THE BOUNDARY OF THE ORIGINAL TRIANGULATION\n");
  for(int k=0;k<(int) t.size();k++){
    for(int l=0;l<3;l++){
      if(n[k].v[l]==-1){
        float length=(p[t[k].v[l]]-p[t[k].v[(l+1)%3]]).norm();
        if(length>max_length){
          int N=length/max_length+1;
          for(int n=1;n<N;n++){
            point2d q=p[t[k].v[l]]+(p[t[k].v[(l+1)%3]]-p[t[k].v[l]])*( (float) n/N);
            //printf("p1=(%lf,%lf) q=(%lf,%lf) p2=(%lf,%lf)\n",p[t[k].v[l]].x,p[t[k].v[l]].y,q.x,q.y,p[t[k].v[(l+1)%3]].x,p[t[k].v[(l+1)%3]].y);
            Delauny(T,q,0);
          }
        }
      }
    }
  }

  ///printf("WE INCLUDE NEW POINTS BY SEGMENT DIVISION\n");
  max_length*=max_length;
  float max_length_segment=max_length+1;
  point2d pnew;
  while(max_length_segment>max_length){
    max_length_segment=0;
    //int psize=T.p.size();
    for(int k=0;k<(int) T.t.size();k++){
      for(int l=0;l<3;l++){
        if(T.t[k].v[l]<4 || T.t[k].v[(l+1)%3]<4) continue;
        float length=(T.p[T.t[k].v[l]]-T.p[T.t[k].v[(l+1)%3]]).norm2();
        if(length>max_length_segment){
          max_length_segment=length;
          pnew=(T.p[T.t[k].v[l]]+T.p[T.t[k].v[(l+1)%3]])*0.5;
        }
      }
    }
    //printf("max_length=%lf max_length_segment=%lf T.t.size()=%d \n",max_length,max_length_segment,T.t.size()); //system("pause");
    if(max_length_segment>max_length){
      if(Delauny(T,pnew,0)<0) break;
      //if(Tn.p.size()==psize) break;
    }
  }


  ///printf("WE REMOVE TRIANGLES WHICH ARE NOT IN THE OLD TRIANGULATION\n");
  for(int k=0;k<(int) T.t.size();k++){
    point2d v=(T.p[T.t[k].v[0]]+T.p[T.t[k].v[1]]+T.p[T.t[k].v[2]])/3.;
    if( inclusion_test(v)>=0 ){ continue; }
    T.t.erase(T.t.begin()+k);
    k--;
  }

  ///printf("WE REMOVE RECTANGLE CORNERS\n");
  for(int k=0;k<4;k++) T.erase_point(0);

  //T.update_n();

  ///printf("WE REMOVE POTENTIAL ISOLATED POINTS\n");
  vector<int> flat(T.p.size(),-1);
  for(int k=0;k<(int) T.t.size();k++){
    flat[T.t[k].v[0]]=flat[T.t[k].v[1]]=flat[T.t[k].v[2]]=0;
  }
  for(int k=0;k<(int) T.p.size();k++){
    if(flat[k]==-1){
      T.erase_point(k);
      flat.erase(flat.begin()+k);
      k--;
    }
  }

  ///printf("WE COMPUTE TRIANGLE NEIGHBORHOOD INFORMATION\n");
  T.update_n();

//  printf("p.size()=%d T.p.size()=%d\n",p.size(),T.p.size());
  ///printf("end division()\n");
  return T;
}


triangulation::triangulation(point2d center,float radius,int Npoints) //constructor of a circle using Npoints
{

  vector<point2d> pol;
  float h=2.*3.141617/Npoints;
  for(int k=0;k<Npoints;k++) pol.push_back(center+point2d(cosf(-k*h),sinf(-k*h))*radius);
  *this=triangulation(pol);
}

affinity connectT12T0(triangulation &T0,int indexTriangle0,int indexSegment0,triangulation &T1,int indexTriangle1,int indexSegment1,float d,int type,point2d p2d)
{
  //printf("indexTriangle0=%d indexTriangle1=%d T0.t.size()=%d T1.t.size()=%d\n",indexTriangle0,indexTriangle1,T0.t.size(),T1.t.size());
  if(indexTriangle0>=(int)T0.t.size() || T0.t.size()==0 || T1.t.size()==0) return affinity();
  if(indexTriangle0<0) indexTriangle0=((double) (T0.t.size()-1)*rand()/RAND_MAX);
  if(indexTriangle1<0 || indexTriangle1>=(int) T1.t.size()) indexTriangle1=(int) ((double) (T1.t.size()-1)*rand()/RAND_MAX);
  if(indexSegment0<0 || indexSegment0>2 ) indexSegment0=0;
  if(indexSegment1<0 || indexSegment1>2 ) indexSegment1=0;
//printf("indexTriangle0=%d indexTriangle1=%d T0.t.size()=%d T1.t.size()=%d\n",indexTriangle0,indexTriangle1,T0.t.size(),T1.t.size());
//printf("indexSegment0=%d indexSegment1=%d \n",indexSegment0,indexSegment1);
  if(T0.n.size()!=T0.t.size()) T0.update_n();
  //bool boundary=false;
//  printf("%d %d %d \n", T0.n[indexTriangle0].v[(indexSegment0+0)%3], T0.n[indexTriangle0].v[(indexSegment0+1)%3], T0.n[indexTriangle0].v[(indexSegment0+2)%3]);
  for(int k=0;k<3;k++){
    if (T0.n[indexTriangle0].v[indexSegment0]<0 /*&& T0.n[indexTriangle0].v[(indexSegment0+1)%3]<0 */){
      //boundary=true;
      break;
    }
    indexSegment0=(indexSegment0+1)%3;
  }
  //if(boundary==false) return affinity();
  if(T1.n.size()!=T1.t.size()) T1.update_n();
  //boundary=false;
  for(int k=0;k<3;k++){
    if (T1.n[indexTriangle1].v[indexSegment1]<0 /*&& T1.n[indexTriangle1].v[(indexSegment1+1)%3]<0 */){
      //boundary=true;
      break;
    }
    indexSegment1=(indexSegment1+1)%3;
  }
  //if(boundary==false) return affinity();
  point2d p00=T0.p[T0.t[indexTriangle0].v[indexSegment0]];
  point2d p01=T0.p[T0.t[indexTriangle0].v[(indexSegment0+1)%3]];
  point2d p02=T0.p[T0.t[indexTriangle0].v[(indexSegment0+2)%3]];
  point2d p10=T1.p[T1.t[indexTriangle1].v[indexSegment1]];
  point2d p11=T1.p[T1.t[indexTriangle1].v[(indexSegment1+1)%3]];
  point2d p12=T1.p[T1.t[indexTriangle1].v[(indexSegment1+2)%3]];
  affinity A(p00,p01,p02,p10,p11,p12,d,type);
  if(p2d.x!=-1000. || p2d.y!=-1000.) A=A+(p2d-(p00+p01)*0.5);

  return A;
}

int triangulation::selectBoundarySegment(int tmin,int tmax,int &tselect)
{
  float w=(float) rand()/RAND_MAX;
  tselect=w*tmin+(1.-w)*tmax;
  if(n.size()!=t.size()) update_n();
  bool boundary=false;
  int indexSegment=(100*rand()/RAND_MAX)%3;
  for(int k=0;k<3;k++){
    if (n[tselect].v[indexSegment]<0){
      boundary=true;
      break;
    }
    indexSegment=(indexSegment+1)%3;
  }
  if(boundary==false) return -1;
  return indexSegment;
}

float triangulation::max_segment_length(int t0){
  if(t0<0 || t0>=(int) t.size() ) return 0.;
  float norm2 = (p[t[t0].v[1]]-p[t[t0].v[0]]).norm2();
  float temp = (p[t[t0].v[2]]-p[t[t0].v[1]]).norm2();
  if(temp>norm2) norm2=temp;
  temp = (p[t[t0].v[2]]-p[t[t0].v[0]]).norm2();
  if(temp>norm2) norm2=temp;

  return sqrtf(norm2);
}


void triangulation::triangles_containing_a_point(int t_index,int p_index,vector<int> &t_select){
  for(int k=0;k<(int) t_select.size();k++){
    if(t_select[k]==t_index) return;
  }
  if(t[t_index].v[0]==p_index || t[t_index].v[1]==p_index || t[t_index].v[2]==p_index){
    t_select.push_back(t_index);
    if(n[t_index].v[0]!=-1) triangles_containing_a_point(n[t_index].v[0],p_index,t_select);
    if(n[t_index].v[1]!=-1) triangles_containing_a_point(n[t_index].v[1],p_index,t_select);
    if(n[t_index].v[2]!=-1) triangles_containing_a_point(n[t_index].v[2],p_index,t_select);
  }
  else return;

}

void triangulation::boundary_segments(float min_segment_length,vector<int> &SelectedTriangles,vector<int> &SelectedSegments, vector<int> &SelectedPoints)
{
   if(t.size()==0 || p.size()==0) return;
   if(t.size() !=n.size()) update_n();
   SelectedTriangles.clear();
   SelectedSegments.clear();
   SelectedPoints.clear();
   vector<int> flat(p.size(),-1);

   int ts=-1;
//ts=0; ss=0;
   while(++ts<(int) t.size()){
     int ss;
     for(ss=0;ss<3;ss++){
       if(n[ts].v[ss]==-1) break;
     }
     if(ss==3) continue;
     if(flat[t[ts].v[ss]]==0 /*|| flat[t[ts].v[(ss+1)%3]]==0*/) continue;

   if(ts==(int) t.size()) return;
   bool newt=true;
   float length=0.; //1.+min_segment_length;
   flat[t[ts].v[ss]]=0;
   int ip=t[ts].v[ss];
   SelectedTriangles.push_back(ts);
   SelectedSegments.push_back(ss);
   SelectedPoints.push_back(ip);
   //length=(p[t[ts].v[ss]]-p[t[ts].v[(ss+1)%3]]).norm()*0.5;
   while(newt==true){
     if(length>min_segment_length){
       SelectedTriangles.push_back(ts);
       SelectedSegments.push_back(ss);
       SelectedPoints.push_back(ip);
       length=0.;
     }
     //printf("length=%lf\n",length);
//     bool newt=false;
//     while(newt==false){
        point2d p2=p[ip];
        vector <int> t_select;
        triangles_containing_a_point(ts,ip,t_select);
        //for( int k=0;k<t_select.size();k++) printf("%d ",t_select[k]); printf("\n"); system("pause");
        //if(t_select.size()<1) break;
        newt=false;
        //printf("t_select.size()=%d ts=%d ss=%d ip=%d\n",t_select.size(),ts,ss,ip);
        for(int k=t_select.size()-1;k>-1;k--){
          int it=t_select[k];
          for(int m=0;m<3;m++){
            if(n[it].v[m]==-1){
              if(flat[t[it].v[m]]==-1 && t[it].v[(m+1)%3]==ip){ip=t[it].v[m]; flat[ip]=0;  ts=it; ss=m; newt=true;  break;}
              if(flat[t[it].v[(m+1)%3]]==-1 && t[it].v[m]==ip){ip=t[it].v[(m+1)%3]; flat[ip]=0;  ts=it; ss=m; newt=true;  break;}
            }
          }
          if(newt==true) break;
        }
        if(newt==false) break;
        length+=(p2-p[ip]).norm();
        //printf("length=%lf \n",length);
    // }
   }}
//   printf("(p[ip0]-p[ip]).norm()=%lf\n",(p[ip0]-p[ip]).norm());
//   if((p[ip0]-p[ip]).norm()<(min_segment_length*0.1)){
//      SelectedTriangles.resize(SelectedTriangles.size()-1);
//      SelectedSegments.resize(SelectedSegments.size()-1);
//   }

   //printf("SelectedTriangles.size()=%d\n",SelectedTriangles.size());
}

void triangulation::boundary_segments(float min_segment_length,vector<int> &SelectedTriangles,vector<int> &SelectedSegments, vector<point2d> &SelectedPoints)
{
   //printf("init boundary_segments() .. ");
   if(t.size()==0 || p.size()==0) return;
   if(t.size() !=n.size()) update_n();
   //write("t.txt");
   SelectedTriangles.clear();
   SelectedSegments.clear();
   SelectedPoints.clear();
   vector<int> flat(p.size(),-1);

   vector<triangle> n2=n;
   for(int k=0;k<(int) t.size();k++){
    for(int m=0;m<3;m++){
      if(n2[k].v[m]<0){ //segment in the contour
        float length0=0.5*min_segment_length;
        int k0=k,m0=m;
        while(true){
          //printf("k0=%d m0=%d t[k0].v[m0]=%d \n",k0,m0,t[k0].v[m0]);
          n2[k0].v[m0]=0;
          float norm=(p[t[k0].v[m0]]-p[t[k0].v[(m0+1)%3]]).norm();
          float length=length0+norm;
          //printf("length=%lf\n",length);
          float Np=1.;
          while(length>Np*min_segment_length){
            float d=Np*min_segment_length-length0;
            SelectedPoints.push_back(p[t[k0].v[m0]]+(p[t[k0].v[(m0+1)%3]]-p[t[k0].v[m0]])*(d/norm));
            SelectedTriangles.push_back(k0);
            SelectedSegments.push_back(m0);
            Np+=1.;
            //printf("p[%d]=(%lf,%lf)\n",SelectedPoints.size()-1,SelectedPoints[SelectedPoints.size()-1].x,SelectedPoints[SelectedPoints.size()-1].y);
          }
          length0=length-(Np-1.)*min_segment_length;
          vector <int> t_select;
          triangles_containing_a_point(k0,t[k0].v[(m0+1)%3],t_select);
          int i;
          for(i=0;i<(int) t_select.size();i++){
            int it=t_select[i];
            //printf("t_select[%d]=%d t_select.size()=%d\n",i,it,t_select.size());
            int j;
            for(j=0;j<3;j++){
              //printf("n2[it].v[j]=%d,t[k0].v[(m0+1)%3]=%d,t[it].v[j]=%d\n",n2[it].v[j],t[k0].v[(m0+1)%3],t[it].v[j]);
              if(n2[it].v[j]==-1){
                if( t[k0].v[(m0+1)%3]==t[it].v[j]){
                  k0=it;
                  m0=j;
                  break;
                }
              }
            }
            if(j<3) break;
          }
          if(i==(int) t_select.size()) break;
        }
        //system("pause");
      }
    }
   }
}

triangulation triangulation::border_stripe_periodicity(){
  if((int) p.size()!=8) return triangulation();

  NodeTree nT;
  float dx=p[4].x-p[0].x;
  float dy=p[4].y-p[0].y;

  float zoom=0.02/dx;
  if(zoom>0.5) zoom=0.5;
  int SymmetryType=5+rand()%2;
  float dmin=128.;
  triangulation T0(point2d(512.-dmin,512.-dmin), point2d(512.+dmin,512.+dmin));
  bool labyrinth=myrand(1.)<0.1?true:false;
  float min_segment_length=200.*zoom,max_segment_length=250.*zoom;
  if(labyrinth==true) max_segment_length=min_segment_length;
  while(true){
    point2d p(512.+dmin*myrand(1.),512.+dmin*myrand(1.));
    if(labyrinth==true) p.y=p.x;
    if(nT.InsideSymmetricRegion(p.x,p.y,512.,512.,SymmetryType)==false) continue;
    nT.point_line_projection(p,point2d(512.,512.),min_segment_length*0.5,SymmetryType);
    nT.n_.push_back(p);
    nT.i_.resize(1);
    nT.r_.push_back(50*zoom);
    if(labyrinth==true) nT.s_.push_back(0.);
    else nT.s_.push_back(0.75);
    break;
  }

  float distance_segment_factor=2.5;
  float min_radius=50.*zoom,max_radius=50*zoom;
  float distance_node_factor=5.;
  float smooth_factor=0.75;
  bool allow_branche_intersection=false;
  int type=0;
  float gamma=1.;
  int Nbranches=1000;
  vector<double> angles=nT.angle_range(100);
  unsigned int MaxNodes=10000;
  unsigned int MaxNeigbors=100;

  if(labyrinth==true){
    smooth_factor=0.;
    angles.clear();
    angles.push_back(0.);
    angles.push_back(M_PI/2);
    angles.push_back(M_PI);
    angles.push_back(3*M_PI/2);
    if(myrand(1.)<0.5){
        angles.push_back(M_PI/4);
        angles.push_back(3*M_PI/4);
        angles.push_back(5*M_PI/4);
        angles.push_back(7*M_PI/4);
    }
  }

  nT.insert_nodes(
    T0,
    min_segment_length,max_segment_length,
    distance_segment_factor,
    min_radius,max_radius,
    distance_node_factor,
    smooth_factor,
    type,
    gamma,
    Nbranches,
    angles,
    MaxNodes,
    MaxNeigbors,
    allow_branche_intersection,
    SymmetryType
  );

  float length=256.;
  float Lx=p[1].x-p[0].x;
  float Ly=p[3].y-p[0].y;
  int Nx=round(Lx/dx);
  int Ny=round(Ly/dy);
  if(Nx<2 || Ny<2) return triangulation();

  NodeTree nT0=nT;
//  nT0.print();
//  system("pause");

  for(int m=1;m<Nx;m++){
    nT.add(nT0.translation(point2d((float) m,0.)*length));
    nT.add(nT0.translation(point2d((float) m,(float) (Ny-1))*length));
  }
  nT.add(nT0.translation(point2d((float) 0.,(float) (Ny-1))*length));
  for(int m=1;m<Ny-1;m++){
    nT.add(nT0.translation(point2d(0.,(float) m)*length));
    nT.add(nT0.translation(point2d((float) (Nx-1),(float) m)*length));
  }

  if(nT.checking()==false){
    printf("nT.checking()=%d\n",false);
    return triangulation();
  }

  triangulation T;
  if(labyrinth==false) T=nT.triangulation_generation(10.);
  else T=nT.triangulation_generation(10.,false);

  double sx=dx/length;
  double sy=dy/length;
  for(unsigned int k=0;k<T.p.size();k++){
    T.p[k].x=p[0].x+sx*(T.p[k].x-384.);
    T.p[k].y=p[0].y+sy*(T.p[k].y-384.);
  }
  return T;

}

/// GENERATION OF A COLLECTION OF POINTS INSIDE A TRIANGULATION
vector<point2d> triangulation::point_generation(
float max_distance, ///max distance between points
int type) /*method to generate the points:
            =0 random
            =1 square lattice
            =2 equilateral triangle lattice
            =3 Delauny
            =4 stripe centerline (persian carpets)*/
{
  float max_distance0=max_distance;
  if(max_distance<0.001) type=0;
  else if(max_distance<0.01 && type==3) type=0;
  //printf("max_distance=%lf type=%d\n",max_distance,type);
  if(max_distance<=0. || t.size()==0) return vector<point2d>();
  float Npf=area()/(max_distance*max_distance);
  if(Npf>=INT_MAX) return vector<point2d>();
  int Np=Npf;
  vector<point2d> pV;
  vector<point2d> e=extrema();// we compute the extrema (xmin,ymin), (xmax,ymax) of triangulation coordinates
  if(e.size()!=2) return vector<point2d>();
  int Nx=(e[1].x-e[0].x)/max_distance;
  float dx=(e[1].x-e[0].x)/Nx;
  int Ny=(e[1].y-e[0].y)/max_distance;
  float dy=(e[1].y-e[0].y)/Ny;
  float h=sqrt(3.)/2.;
  max_distance=dx;

//printf("e[0]=(%lf,%lf) e[1]=(%lf,%lf)\n",e[0].x,e[0].y,e[1].x,e[1].y); system("pause") ;
  switch(type%5) {
    case 0 : //random
    {
      vector<float> Area(t.size());
      Area[0]=AREA(p[t[0].v[0]],p[t[0].v[1]],p[t[0].v[2]]);
      for(int k=1;k<(int) t.size();k++) Area[k]=Area[k-1]+AREA(p[t[k].v[0]],p[t[k].v[1]],p[t[k].v[2]]);
      float AreaT=Area[t.size()-1];
      while((int) pV.size()<Np){
        float At=AreaT*rand()/RAND_MAX;
        int it=0;
        while(it<(int) t.size() && Area[it]<At) it++;
        //int it=t.size()*((float) rand()/RAND_MAX);
        //printf("%d ",it);
        //it=it%t.size();
        float w1=1.;
        float w2=1.;
        while((w1+w2)>1.){
          w1=(float) rand()/RAND_MAX;
          w2=(float) rand()/RAND_MAX;
        }
        point2d p2=p[t[it].v[0]]+(p[t[it].v[1]]-p[t[it].v[0]])*w1+(p[t[it].v[2]]-p[t[it].v[0]])*w2;

        pV.push_back(p2);
      }
    }
      break;       // and exits the switch
    case 1 : // square lattice
      for(int x=0;x<Nx;x++){
        for(int y=0;y<Ny;y++){
          point2d p2(e[0].x+x*dx+dx/2,e[0].y+dy/2+y*dy);
          if(inclusion_test(p2)>=0) pV.push_back(p2);
        }
      }
      break;
    case 2 : // equilateral triangle lattice
      {
      int Nx=(e[1].x-e[0].x-max_distance)/(h*max_distance);
      //printf("N=%d\n",N); system("pause");
      if(Nx<=0) break;
      if(Nx%2==1) Nx++;
      float dx=(e[1].x-e[0].x-max_distance)/Nx;
      for(float x=e[0].x+max_distance/2;x<=e[1].x/*-max_distance/2+1e-8*/;x+=2.*dx){
        for(float y=e[0].y+max_distance/2;y<=e[1].y+1e-8/*-max_distance/2*/;y+=dx/h/*max_distance*/){
           point2d p2(x,y);
           if(inclusion_test(p2)>=0) pV.push_back(p2);
        }
      }
      for(float x=e[0].x+max_distance/2+dx;x<=e[1].x/*-max_distance/2+1e-8*/;x+=2.*dx){
        for(float y=e[0].y/*+max_distance/2+max_distance*0.5*/;y<=e[1].y+1e-8/*-max_distance/2*/;y+=dx/h/*max_distance*//*max_distance*/){
           point2d p2(x,y);
           if(inclusion_test(p2)>=0) pV.push_back(p2);
        }
      }
      }
      break;
    case 3 : // Delauny
      {
        triangulation T2=division(max_distance*1.5);
        for(int k=0;k<(int) T2.p.size();k++){
          if(T2.fr[k]==0) pV.push_back(T2.p[k]);
        }
      }
      break;
    case 4 : // stripe centerline (persian carpets)
      {
        if((int) p.size()!=8) return point_generation(max_distance0,3);
        float dx=p[4].x-p[0].x;
        float dy=p[4].y-p[0].y;
        vector<point2d> pT;
        pT.push_back(p[0]+point2d(dx/2.,dy/2.));
        pT.push_back(p[1]+point2d(-dx/2.,dy/2.));
        pT.push_back(p[2]+point2d(-dx/2.,-dy/2.));
        pT.push_back(p[3]+point2d(dx/2.,-dy/2.));
        pT.push_back(p[0]+point2d(dx/2.,dy/2.));
        for(int m=0;m<4;m++){
          point2d p1=pT[m];
          point2d p2=pT[m+1];
          point2d p12=p2-p1;
          float norm=p12.norm();
          p12=p12/norm;
          pV.push_back(p1);
          int Np=norm/max_distance;
          if(Np>1){
            float md2=norm/Np;
            for(int k=1;k<Np;k++){
              pV.push_back(p1+p12*(k*md2));
            }
          }
        }
      }
      break;
  }

  if(pV.size()==0){
    point2d averg(0.,0.);
    for(int k=0;k<(int) p.size();k++) averg=averg+p[k];
    point2d pa=averg/p.size();
    if(inclusion_test(pa)>=0) pV.push_back(pa);
    return pV;
  }

//  /// sorting the point position
//  for(int k=0;k<pV.size();k++){
//    int j=rand()%pV.size();
//    point2d temp=pV[j];
//    pV[j]=pV[k];
//    pV[k]=temp;
//  }

  return pV;
}

void triangulation::fix_triangle_orientation(){
  for(int k=0;k<(int) t.size();k++){
    if(R(p[t[k].v[0]],p[t[k].v[1]],p[t[k].v[2]])<0){
        int temp=t[k].v[0];
        t[k].v[0]=t[k].v[1];
        t[k].v[1]=temp;
    }
  }
}

/// RETURN A RANDOM VALUE IN INTEGER PRECISION BETWEEN 0 AND THE PARAMETER In_
double RAND_MAX2=(double) RAND_MAX*RAND_MAX;
int ami_randI(int In_){
  int I=((double) rand()/RAND_MAX+rand()/RAND_MAX2)*In_;
  return I<In_?I:I-1;
}

/// SPLINE INTERPOLATION
vector<point2d> point_generation_splines(
vector<point2d> &pV, /// control points of the curve
point2d d0, point2d dN, /// tangents to the curve in the first and last points
float segment_length) /// default segment length for discretization of the curve
{
  vector<point2d> pVn;

  /// WE CREATE THE SPLINES
  int N=pV.size();
  vector<float> t(N),z(N);
  t[0]=0;
  for(int k=1;k<N;k++) t[k]=t[k-1]+(pV[k]-pV[k-1]).norm();
  float **xS,**yS;
  ami_calloc2d(xS,float,N,4);
  ami_calloc2d(yS,float,N,4);

  for(int k=0;k<N;k++) z[k]=pV[k].x;
  ami_crear_splines(t,z,xS,d0.x,dN.x);

  for(int k=0;k<N;k++) z[k]=pV[k].y;
  ami_crear_splines(t,z,yS,d0.y,dN.y);

  /// WE BUILD THE CHAINE OF POINTS
  for(float q=0;q<=t[N-1];q+=segment_length){
    float xn=ami_evaluar_splines(t,xS,q);
    float yn=ami_evaluar_splines(t,yS,q);
    //printf("xn=%lf yn=%lf\n",xn,yn);
    pVn.push_back(point2d(xn,yn));
  }
  pVn.push_back(pV[pV.size()-1]);

  /// WE FREE THE MEMORY
  ami_free2d(xS);
  ami_free2d(yS);

  return pVn;
}

/// return the index of a point. If the point is not in the list it is included.
int triangulation::point_index(point2d q){
  for(int k=0;k<(int) p.size();k++){
    if( (p[k]-q).norm2()<1e-10 ) return k;
  }
  p.push_back(q);
  return p.size()-1;
}


///***********************************************************************
///           FUNCTION TO COMPUTE THE AVERAGE OF TO ANGLES
///***********************************************************************
double angle_average(double a1,double a2){
 if(a1<=a2) return (a1+a2)*0.5;
 else return (a1+a2+2*M_PI)*0.5;
};

///***********************************************************************
///        FUNCTION TO COMPUTE THE DIFFERENCE BETWEEN 2 ANGLES
///***********************************************************************
double angle_dif(double a1,double a2){
 double dif;
 if(a1<=a2) dif=a2-a1;
 else dif=a2+2*M_PI-a1;
 //printf("a1=%lf,a2=%lf,dif=%lf\n",a1*180./M_PI,a2*180./M_PI,dif*180./M_PI);
 return dif;
};

///  *********************************************************
///   constructor of a decagon symmetric with respect y axis
///  ***********************************************************
triangulation::triangulation(point2d center,float y[],float Rx[] )
{

  p.resize(11);
  p[0]=center;
  unsigned int pos=1;
  for(unsigned int k=0;k<5;k++){
    p[pos++]=center+point2d(Rx[k],y[k]);
  }
  for(unsigned int k=0;k<5;k++){
    p[pos++]=center+point2d(-Rx[4-k],y[4-k]);
  }
  for(unsigned int k=1;k<10;k++){
    t.push_back(triangle(0,k,k+1));
  }
  t.push_back(triangle(0,1,10));

}


/// AUXILIARY FUNCTION OF BranchedShapeGenerator()
point2d find_point_triangulation_1024x1024(triangulation &Tinit){
  point2d p0;
  while(true){
    p0.x=rand()%1024;
    p0.y=rand()%1024;
    //printf("%d ",Tinit.inclusion_test(p0));
    if( Tinit.inclusion_test(p0)>=0) break;
  }
  return p0;
}


/// **********************************************************************************
/// ******* GENERATOR OF BRANCHED SHAPES **************************************
/// **********************************************************************************
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
)
{
  //printf("init BranchedShapeGenerator()\n");
  //printf("Tinit.p.size()=%d\n",Tinit.p.size());


  srand(random_seed);  rand();

  /// CHECKING PARAMETERS
  if(shape_type<0 || shape_type>10){
    char mes[300];
    sprintf(mes,"The shape type (%d) is out of range [0,10]\n",shape_type+1);
    printf("%s\n",mes);
    //fprintf_demo_failure(mes);
    return NodeTree();
  }

  if(average_node_distance>=0 && average_node_distance<2 ) average_node_distance=2;
  if(average_radius>=0 && average_radius<1 ) average_radius=1;
  if(MaxNeigbors==0 || MaxNeigbors==1) MaxNeigbors=2;
  if(shape_type>2 &&  MaxNeigbors<=2) MaxNeigbors=3;
  if(shape_type!=1 && smoothing_factor<0.001) smoothing_factor=0.001;
  if(shape_type<8 && shape_type>4 && average_node_distance<100) node_radius_reduction_factor=sqrt(node_radius_reduction_factor);
  if(average_node_distance>=0 && average_node_distance<2 ) average_node_distance=2;
  if(average_radius>=0 && average_radius<1 ) average_radius=1;


  triangulation Tinit0=Tinit;
  float average_node_distance0=average_node_distance;
  float average_radius0=average_radius;

  /// NORMALIZATION OF THE NODE DISTANCE
  average_node_distance/=150.;

  /// NODE TREE STRUCTURE
  NodeTree nT;

  /// AUXILIARY VARIABLES
  float node_radius_reduction_factor_sqrt=sqrt(node_radius_reduction_factor);
  srand(random_seed); rand();
  bool round_termination = branched_finished_in_peak==0?true:false;
  float zoom=average_node_distance; //0.1+myrand(1.)*0.9;
  float dmin=64.;
  float R0=512-dmin;
  int type=0;
  int Nbranches=100000;
  float gamma=1.;


  ///-------------------------------------------------------------------------------------------------------
  /// MATISSE 13 GENERATION FROM A CENTRAL POINT
  if(shape_type==0)
  {

    if(Tinit.p.size()==0){
      /// BUILD A DECAGON TO CONSTRAINT THE SHAPE
      float T0y[5]={R0*sinf(2*M_PI/5),R0*sinf(M_PI/5),0,-R0*sinf(M_PI/5),-R0*sinf(2*M_PI/5)};
      float T0R[5]={R0*cosf(2*M_PI/5),R0*cosf(M_PI/5),R0,R0*cosf(M_PI/5),R0*cosf(2*M_PI/5)};
      Tinit=triangulation(point2d(512.,512.),T0y,T0R);
    }

    /// INSERT THE FIRST NODE IN THE TREE
    point2d p0=find_point_triangulation_1024x1024(Tinit);
    //printf("p0=%lf %lf\n",p0.x,p0.y); system("pause");
    nT.n_.push_back(p0);
    nT.i_.resize(1);
    nT.r_.push_back(average_radius*1.5*zoom);
    nT.s_.push_back(0.75);

    /// DEFINE VARIABLES FOR THE TREE GENERATION
    float min_segment_length=200.*zoom,max_segment_length=250.*zoom;
    float distance_segment_factor=2.5;
    float min_radius=average_radius*zoom,max_radius=min_radius*6./5.;
    float distance_node_factor=5.;
    if(allow_branche_intersection==true){
      min_segment_length=175.*zoom;
      max_segment_length=275.*zoom;
      distance_node_factor=4.;
    }

    /// VARIABLES UPDATE ACCORDING TO node_radius_reduction_factor
    float smy=myrand(1.);
    distance_node_factor*=pow(node_radius_reduction_factor,0.5+5.*smy);
    distance_segment_factor*=pow(node_radius_reduction_factor,0.5+1.*smy);

    /// ANGLE COMPUTATION
    vector<double> angles=nT.angle_range(100);
    if(MaxAngles<(int) angles.size()){
      if(rand()%10==0) angles.resize(MaxAngles);
      else  angles=nT.angle_range_rand(MaxAngles);
    }

    /// INCLUSION OF NODES IN THE TREE
    int cont=0;
    while(true){
      cont++;

      /// WE STORE THE CURRENT NUMBER OF NODES
      int Nnodes=nT.n_.size();

      /// UPDATE OF VARIABLES
      min_segment_length*=node_radius_reduction_factor;
      max_segment_length*=node_radius_reduction_factor;
      min_radius*=node_radius_reduction_factor;
      max_radius*=node_radius_reduction_factor;
      if(min_radius<1. || min_segment_length<2.) break;
      distance_segment_factor/=node_radius_reduction_factor_sqrt;
      distance_node_factor/=node_radius_reduction_factor;

      float gamma2=gamma;
      bool allow_branche_intersection2=false;
      float min_segment_length2=min_segment_length;
      float max_segment_length2=max_segment_length;
      float distance_node_factor2=distance_node_factor;
      int MaxNeigbors2=MaxNeigbors;
      if(allow_branche_intersection==true){
        distance_node_factor2*=1.8/pow(node_radius_reduction_factor,1.);
        if(cont%2==0){
          gamma2=1.;
          allow_branche_intersection2=true;
          min_segment_length2=max_segment_length;
          max_segment_length2*=1.2;
          MaxNeigbors2=3;
        }
      }

      /// INSERCIÓN OF NEW NODES
      nT.insert_nodes(
        Tinit,
        min_segment_length2,
        max_segment_length2,
        distance_segment_factor,
        min_radius,max_radius,
        distance_node_factor2,
        smoothing_factor,
        type,
        gamma2,
        1,// Nbranches,
        angles,
        MaxNodes,
        MaxNeigbors2,
        allow_branche_intersection2
      );
      if(Nnodes==(int) nT.n_.size()) break;
    }

#ifdef DEBUG
if(nT.checking()==false){
    printf("checking fails before joint_extrema()\n");
    char filename[400]="nT.txt";
    nT.write(filename);
    system("pause");
}
#endif

    /// WE MANAGE THE JOIN OF NODES EXTREMA OF THE TREE
    nT.joint_extrema(joint_extrema_probability);



#ifdef DEBUG
if(nT.checking()==false){
    printf("checking fails after joint_extrema()\n");
    char filename2[400]="nT2.txt";
    nT.write(filename2);
    system("pause");
}
#endif

    /// WE MANAGE IF THE NODES EXTREMA ARE PICKS
    if(round_termination==false){
      for(unsigned int k=0;k<nT.n_.size();k++){
        if(nT.i_[k].size()==1) {
          nT.r_[k]=1;
          nT.s_[k]=0.1;
          nT.n_[k]=nT.n_[nT.i_[k][0]]+(nT.n_[k]-nT.n_[nT.i_[k][0]])*0.75;
        }
      }
    }

    /// WE NORMALIZE THE NODE COORDINATES TO 1024x1024
    nT.normalize();

  }

  ///-------------------------------------------------------------------------------------------------------
  /// MATISSE 11 LABERINTHIC
  else if(shape_type==1)
  {
    if(Tinit.p.size()==0){
      /// BUILD A DECAGON TO CONSTRAINT THE SHAPE
      float T0y[5]={R0*sinf(2*M_PI/5),R0*sinf(M_PI/5),0,-R0*sinf(M_PI/5),-R0*sinf(2*M_PI/5)};
      float T0R[5]={R0*cosf(2*M_PI/5),R0*cosf(M_PI/5),R0,R0*cosf(M_PI/5),R0*cosf(2*M_PI/5)};
      Tinit=triangulation(point2d(512.,512.),T0y,T0R);
    }

    /// INSERT THE FIRST NODE IN THE TREE
    point2d p0=find_point_triangulation_1024x1024(Tinit);
    nT.n_.push_back(p0);
    nT.i_.resize(1);
    nT.r_.push_back(zoom*average_radius);
    nT.s_.push_back(smoothing_factor);

    /// DEFINE VARIABLES FOR THE TREE GENERATION
    float min_segment_length=200.*zoom,max_segment_length=200.*zoom;
    float distance_segment_factor=2.5;
    float min_radius=zoom*average_radius,max_radius=zoom*average_radius;
    float distance_node_factor=5.;
    if(allow_branche_intersection==true){
      min_segment_length=175.*zoom;
      max_segment_length=275.*zoom;
      distance_node_factor=4.;
    }
    float smy=myrand(1.);
    distance_node_factor*=pow(node_radius_reduction_factor,0.5+2.*smy);
    distance_segment_factor*=pow(node_radius_reduction_factor,0.5+1.*smy);

    /// DEFINE THE INITIAL VECTOR OF DIRECTIONS
    vector<double> angles;
    angles.push_back(0.);
    angles.push_back(M_PI/2);
    angles.push_back(M_PI);
    angles.push_back(3*M_PI/2);
    angles.push_back(M_PI/4);
    angles.push_back(3*M_PI/4);
    angles.push_back(5*M_PI/4);
    angles.push_back(7*M_PI/4);

    /// WE MANAGE 3 COLLECTION OF ANGLES
    int c=rand()%3;
    if(MaxAngles==4 || ( MaxAngles>8 && c==0 ) ){
      angles.resize(4);
    }
    else if(MaxAngles==3 || ( MaxAngles>8 && c==1 ) ){
      angles.resize(3);
      angles[1]=2.*M_PI/3.;
      angles[2]=4.*M_PI/3.;
    }
    else if(MaxAngles<8){
      angles.resize(MaxAngles);
    }

    /// INCLUSION OF NODES IN THE TREE
    int cont=0;
    while(true){
    cont++;
    int Nnodes=nT.n_.size();
    min_segment_length*=node_radius_reduction_factor;
    max_segment_length*=node_radius_reduction_factor;
    min_radius*=node_radius_reduction_factor;
    max_radius*=node_radius_reduction_factor;
    if(min_radius<1. || min_segment_length<2.) break; //min_radius=1;
    distance_segment_factor/=node_radius_reduction_factor_sqrt;
    distance_node_factor/=node_radius_reduction_factor;

    float gamma2=gamma;
    bool allow_branche_intersection2=false;
    float min_segment_length2=min_segment_length;
    float max_segment_length2=max_segment_length;
    float distance_node_factor2=distance_node_factor;
    int MaxNeigbors2=MaxNeigbors;
    if(allow_branche_intersection==true){
      distance_node_factor2*=1.8/pow(node_radius_reduction_factor,1.5);
      if(cont%2==0){
        gamma2=1.;
        allow_branche_intersection2=true;
        min_segment_length2=max_segment_length;
        max_segment_length2*=1.2;
        MaxNeigbors2=3;
      }
    }

    /// INSERCIÓN OF NEW NODES
    nT.insert_nodes(
      Tinit,
      min_segment_length2,
      max_segment_length2,
      distance_segment_factor,
      min_radius,max_radius,
      distance_node_factor2,
      smoothing_factor,
      type,
      gamma2,
      1,// Nbranches,
      angles,
      MaxNodes,
      MaxNeigbors2,
      allow_branche_intersection2
    );
    if(Nnodes==(int) nT.n_.size()) break;
    }

    /// WE MANAGE THE JOIN OF NODES EXTREMA OF THE TREE
    nT.joint_extrema(joint_extrema_probability);

    /// WE MANAGE IF THE NODES EXTREMA ARE PICKS
    if(round_termination==false){
      for(unsigned int k=0;k<nT.n_.size();k++){
        if(nT.i_[k].size()==1) {
          nT.r_[k]=1;
          nT.s_[k]=0.1;
          nT.n_[k]=nT.n_[nT.i_[k][0]]+(nT.n_[k]-nT.n_[nT.i_[k][0]])*0.75;
        }
      }
    }

    /// WE NORMALIZE THE NODE COORDINATES TO 1024x1024
    nT.normalize();
  }

  ///-------------------------------------------------------------------------------------------------------
  /// MATISSE 12 RANDOM DISTRIBUTION
  else if(shape_type==2)
  {
    if(Tinit.p.size()==0){
      /// BUILD A DECAGON TO CONSTRAINT THE SHAPE
      float T0y[5]={R0*sinf(2*M_PI/5),R0*sinf(M_PI/5),0,-R0*sinf(M_PI/5),-R0*sinf(2*M_PI/5)};
      float T0R[5]={R0*cosf(2*M_PI/5),R0*cosf(M_PI/5),R0,R0*cosf(M_PI/5),R0*cosf(2*M_PI/5)};
      Tinit=triangulation(point2d(512.,512.),T0y,T0R);
    }

     /// INSERT THE FIRST NODE IN THE TREE
     point2d p0=find_point_triangulation_1024x1024(Tinit);

     float min_smooth_factor=0.1+myrand(1.)*0.25*smoothing_factor;
     float max_smooth_factor=smoothing_factor;

     nT.n_.push_back(p0);
     nT.i_.resize(1);
     nT.r_.push_back(zoom*average_radius);
     nT.s_.push_back(min_smooth_factor+myrand(1.)*(max_smooth_factor-min_smooth_factor));

     /// VARIABLES UPDATE
     float min_segment_length=200.*zoom,max_segment_length=250.*zoom;
     float distance_segment_factor=2.;
     float min_radius=20.*zoom*average_radius/50.,max_radius=100*zoom*average_radius/50.;
     float distance_node_factor=4.;
     float max_distance_node_factor=32.;
     if(allow_branche_intersection==true){
       min_segment_length=175.*zoom;
       max_segment_length=275.*zoom;
       distance_node_factor=3.;
     }
     double min_angle_distance=M_PI/16.;
     unsigned int MaxAttempts=1000;


     /// INCLUSION OF NODES IN THE TREE
     nT.random_generation(
       Tinit,
       min_segment_length,max_segment_length,
       distance_segment_factor,
       min_radius,max_radius,
       distance_node_factor,
       max_distance_node_factor,
       min_smooth_factor,max_smooth_factor,
       gamma,
       min_angle_distance,
       MaxNodes,
       MaxNeigbors,
       MaxAttempts,
       allow_branche_intersection
     );

    /// WE MANAGE THE JOIN OF NODES EXTREMA OF THE TREE
    nT.joint_extrema(joint_extrema_probability);

    /// WE MANAGE IF THE NODES EXTREMA ARE PICKS
    if(round_termination==false){
      for(unsigned int k=0;k<nT.n_.size();k++){
        if(nT.i_[k].size()==1) {
          nT.r_[k]=1;
          nT.s_[k]=0.1;
          nT.n_[k]=nT.n_[nT.i_[k][0]]+(nT.n_[k]-nT.n_[nT.i_[k][0]])*0.75;
        }
      }
    }

    /// WE NORMALIZE THE NODE COORDINATES TO 1024x1024
    nT.normalize();
  }

  ///-------------------------------------------------------------------------------------------------------
  /// MATISSE 5 SPIDER
  else if(shape_type==3)
  {
    /// BUILD A DECAGON TO CONSTRAINT THE SHAPE
    float T0y[5]={R0*sinf(2*M_PI/5),R0*sinf(M_PI/5),0,-R0*sinf(M_PI/5),-R0*sinf(2*M_PI/5)};
    float T0R[5]={R0*cosf(2*M_PI/5),R0*cosf(M_PI/5),R0,R0*cosf(M_PI/5),R0*cosf(2*M_PI/5)};
    Tinit=triangulation(point2d(512.,512.),T0y,T0R);

    nT.n_.push_back(point2d(512.,512.));
    nT.i_.resize(1);
    nT.r_.push_back(70);
    nT.s_.push_back(smoothing_factor);

    float min_segment_length=200.,max_segment_length=250.;
    float distance_segment_factor=2.;
    float min_radius=30.,max_radius=50.;
    float distance_node_factor=4.;


    Nbranches=1;
    vector<double> angles=nT.angle_range(100);
    if(MaxAngles<(int) angles.size()) angles.resize(MaxAngles);

    if(allow_branche_intersection==true){
      min_segment_length=200.;
      max_segment_length=250.;
      distance_node_factor=4.;
    }

    float smy=myrand(1.);
    distance_node_factor*=pow(node_radius_reduction_factor,0.5+2.*smy);
    distance_segment_factor*=pow(node_radius_reduction_factor,0.5+1.*smy);
    min_radius*=node_radius_reduction_factor;
    max_radius*=node_radius_reduction_factor;

    /// INCLUSION OF NODES IN THE TREE
    nT.insert_nodes(
      Tinit,
      min_segment_length,max_segment_length,
      distance_segment_factor,
      min_radius,max_radius,
      distance_node_factor,
      smoothing_factor,
      type,
      gamma,
      Nbranches,
      angles,
      MaxNodes,
      MaxNeigbors,
      allow_branche_intersection
    );

    min_radius=3.;max_radius=3;
    distance_segment_factor=20.;
    min_segment_length=200.;
    max_segment_length=250.;
    distance_node_factor=1500.;
    Nbranches=1;

    /// INCLUSION OF NODES IN THE TREE
    nT.insert_nodes(
      Tinit,
      min_segment_length,max_segment_length,
      distance_segment_factor,
      min_radius,max_radius,
      distance_node_factor,
      smoothing_factor,
      type,
      gamma,
      Nbranches,
      angles,
      MaxNodes,
      2,
      allow_branche_intersection
    );

    /// WE NORMALIZE THE NODE COORDINATES TO 1024x1024
    nT.normalize();

  }

  ///-------------------------------------------------------------------------------------------------------
  /// MATISSE 6 WATER DROP
  else if(shape_type==4)
  {
    Tinit=triangulation(point2d(0.,0.),point2d(1024.,1024.));
    R0=150+100*myrand(1.);
    nT.clear();

    /// INSERT THE FIRST NODE IN THE TREE
    nT.n_.push_back(point2d(512.,dmin+R0));
    nT.i_.resize(1);
    nT.r_.push_back(R0*average_radius/50);
    nT.s_.push_back(smoothing_factor);
    nT.insert_point( point2d(512.+(-1.+2.*myrand(1.))*30,1024.-dmin) , 1., 0.1, 0);

    /// ANGLE COMPUTATION
    vector<double> angles=nT.angle_range(13);
    int posa=myrand(1.)*angles.size();
    posa=posa%angles.size();
    nT.insert_point( nT.n_[1]-point2d(cos(angles[posa]),sin(angles[posa]))*3 , 5., 1., 0,1);

    /// WE NORMALIZE THE NODE COORDINATES TO 1024x1024
    nT.normalize();
  }

  ///-------------------------------------------------------------------------------------------------------
  /// MATISSE 8 (TREE)
  else if(shape_type==5)
  {
    /// BUILD A DECAGON TO CONSTRAINT THE SHAPE
    float T0y[5]={348.,212.,0.,-238.,-448.};
    float T0R[5]={10.,128.,300.,200.,200};
    if(Tinit.p.size()>2){
      vector<point2d> eV=Tinit.extrema();
      float scalex=(eV[1].x-eV[0].x);
      float scaley=(eV[1].y-eV[0].y);
      float scale=scalex/scaley;
      for(int k=0;k<5;k++) T0R[k]*=scale;
    }


    Tinit=triangulation(point2d(512.,512.),T0y,T0R);

    float min_segment_length=150.*zoom,max_segment_length=200.*zoom;
    float distance_segment_factor=8.;
    float min_radius=average_radius*zoom,max_radius=average_radius*zoom;
    float distance_node_factor=4.;

    gamma=1.5;
    vector<double> angles=nT.angle_range(3);

    /// INSERT THE FIRST NODE IN THE TREE
    nT.clear();
    nT.n_.push_back(point2d(512.,dmin+30.));
    nT.i_.resize(1);
    nT.r_.push_back(average_radius*zoom);
    nT.s_.push_back(smoothing_factor);

    /// INCLUSION OF NODES IN THE MAIN BRANCHE OF THE TREE
    int cont=0;
    while(true){
      cont++;
      int Nnodes=nT.n_.size();
      average_radius*=node_radius_reduction_factor_sqrt;
      min_segment_length*=node_radius_reduction_factor_sqrt;
      max_segment_length*=node_radius_reduction_factor_sqrt;
      min_radius*=node_radius_reduction_factor_sqrt;
      max_radius*=node_radius_reduction_factor_sqrt;
      if(min_radius<1. || min_segment_length<2.) break; //min_radius=1;
      distance_segment_factor/=node_radius_reduction_factor_sqrt;
      distance_node_factor/=node_radius_reduction_factor_sqrt;

      nT.insert_nodes(
        Tinit,
        min_segment_length,max_segment_length,
        distance_segment_factor,
        min_radius*1.3,max_radius*1.3,
        distance_node_factor,
        smoothing_factor,
        type,
        gamma,
        1,// Nbranches,
        angles,
        MaxNodes,
        100,
        false //allow_branche_intersection
      );
      if(Nnodes==(int) nT.n_.size()) break;
    }

    /// INCLUSION OF NODES IN THE LATERAL BRANCHES OF THE TREE
    gamma=3;
    Nbranches=1;
    angles=nT.angle_range(15);
    distance_segment_factor=2.5/pow(node_radius_reduction_factor_sqrt,(float) cont);
    min_segment_length=100.*zoom*pow(node_radius_reduction_factor_sqrt,(float) cont);
    max_segment_length=150.*zoom*pow(node_radius_reduction_factor_sqrt,(float) cont);


    /// INCLUSION OF NODES IN THE LATERAL BRANCHES OF THE TREE
    nT.insert_nodes(
      Tinit,
      min_segment_length,max_segment_length,
      distance_segment_factor/2,
      min_radius,max_radius,
      distance_node_factor/2,
      smoothing_factor,
      type,
      gamma,
      Nbranches,
      angles,
      MaxNodes,
      5,
      false //allow_branche_intersection
   );


    Nbranches=10000;
    angles=nT.angle_range(15);
    if(MaxAngles<(int) angles.size()) angles.resize(MaxAngles);
    distance_segment_factor=2.5/pow(node_radius_reduction_factor_sqrt,(float) cont);
    distance_node_factor=5./pow(node_radius_reduction_factor,(float) cont);

    min_segment_length=100.*zoom;
    max_segment_length=150.*zoom;

//      if(allow_branche_intersection==true){
//        min_segment_length=100.*zoom;
//        max_segment_length=250.*zoom;
//        distance_node_factor=3;
//      }

    float smy=myrand(1.);
    distance_node_factor*=pow(node_radius_reduction_factor,0.5+2.*smy);
    distance_segment_factor*=pow(node_radius_reduction_factor,0.5+1.*smy);

    /// INCLUSION OF NEW NODES IN THE TREE
    cont=0;
    while(true){
      cont++;
      int Nnodes=nT.n_.size();
      min_segment_length*=node_radius_reduction_factor;
      max_segment_length*=node_radius_reduction_factor;
      min_radius*=node_radius_reduction_factor;
      max_radius*=node_radius_reduction_factor;
      if(min_radius<1. || min_segment_length<2.) break; //min_radius=1;
      distance_segment_factor/=node_radius_reduction_factor_sqrt;
      distance_node_factor/=node_radius_reduction_factor;

      float gamma2=gamma;
      bool allow_branche_intersection2=false;
      float min_segment_length2=min_segment_length;
      float max_segment_length2=max_segment_length;
      float distance_node_factor2=distance_node_factor;
      if(cont%2==0 && allow_branche_intersection==true){
        gamma2=1.;
        allow_branche_intersection2=true;
        min_segment_length2=max_segment_length;
        max_segment_length2*=1.2;
        distance_node_factor2*=1.;
      }

      nT.insert_nodes(
        Tinit,
        min_segment_length2,
        max_segment_length2,
        distance_segment_factor/2,
        min_radius,max_radius,
        distance_node_factor2,
        smoothing_factor,
        type,
        gamma2,
        1,// Nbranches,
        angles,
        MaxNodes,
        3,
        allow_branche_intersection2
      );
      if(Nnodes==(int) nT.n_.size()) break;
    }

   /// WE MANAGE THE JOIN OF NODES EXTREMA OF THE TREE
   nT.joint_extrema(joint_extrema_probability);

   /// WE MANAGE IF THE NODES EXTREMA ARE PICKS
   if(round_termination==false){
      for(unsigned int k=0;k<nT.n_.size();k++){
        if(nT.i_[k].size()==1) {
          nT.r_[k]=1;
          nT.s_[k]=0.1;
          nT.n_[k]=nT.n_[nT.i_[k][0]]+(nT.n_[k]-nT.n_[nT.i_[k][0]])*0.75;
        }
      }
    }

    /// WE NORMALIZE THE NODE COORDINATES TO 1024x1024
    nT.normalize();
  }


  ///-------------------------------------------------------------------------------------------------------
  /// MATISSE 10 SEAWOOD 1
  else if(shape_type==6)
  {
    node_radius_reduction_factor=pow(node_radius_reduction_factor,0.75);
    /// BUILD A DECAGON TO CONSTRAINT THE SHAPE
    float T0y[5]={R0*sinf(2*M_PI/5),R0*sinf(M_PI/5),0,-R0*sinf(M_PI/5),-R0*sinf(2*M_PI/5)};
    float T0R[5]={R0*cosf(2*M_PI/5),R0*cosf(M_PI/5),R0,R0*cosf(M_PI/5),R0*cosf(2*M_PI/5)};

   if(Tinit.p.size()>2){
      vector<point2d> eV=Tinit.extrema();
      float scalex=(eV[1].x-eV[0].x);
      float scaley=(eV[1].y-eV[0].y);
      float scale=scalex/scaley;
      for(int k=0;k<5;k++) T0R[k]*=scale;
    }


    Tinit=triangulation(point2d(512.,512.),T0y,T0R);

    nT.clear();
    nT.n_.push_back(point2d(512.,512.-0.9*R0));
    nT.i_.resize(1);
    nT.r_.push_back(zoom*average_radius);
    nT.s_.push_back(smoothing_factor);

    {

      for(unsigned int k=1;(float) k<4./zoom;k++){
        float scale=1.; //=pow(node_radius_reduction_factor,(float) k-1);
        nT.insert_point(nT.n_[nT.n_.size()-1]+point2d(zoom*average_radius*(2.+0.5*myrand(1.)),myrand(1.)*20*zoom*scale),zoom*scale*average_radius,0.75,nT.n_.size()-1);
        if(nT.n_[nT.n_.size()-1].x-1.2*T0R[4]>512.) break;
      }
      nT.insert_point(nT.n_[0]+point2d(-zoom*average_radius*(2.+0.5*myrand(1.)),myrand(1.)*20*zoom),zoom*average_radius,0.75,0);
      for(unsigned int k=2;(float) k<4./zoom;k++){
        float scale=1.; //pow(node_radius_reduction_factor,(float) k-2);
        nT.insert_point(nT.n_[nT.n_.size()-1]+point2d(-zoom*average_radius*(2.+0.5*myrand(1.)),myrand(1.)*20*zoom*scale),zoom*scale*average_radius,0.75,nT.n_.size()-1);
        if(nT.n_[nT.n_.size()-1].x+T0R[4]*1.2<512.) break;
      }
    }

    float min_segment_length=100.*zoom,max_segment_length=150.*zoom;
    float distance_segment_factor=2.;
    //float min_radius=30.*zoom*average_radius/50.,max_radius=40*zoom*average_radius/50.;
    float min_radius=40.*zoom*average_radius/50.,max_radius=zoom*average_radius;
    float distance_node_factor=4.;
    if(allow_branche_intersection==true){
      min_segment_length=100.*zoom;
      max_segment_length=200.*zoom;
      distance_node_factor=4.;
    }

    gamma=3.;
    Nbranches=6/zoom;
    //printf("Nbranches=%d\n",Nbranches);
    vector<double> angles=nT.angle_range(11);
    if(MaxAngles<(int) angles.size()) angles.resize(MaxAngles);

    /// INCLUSION OF NODES IN THE TREE
    int cont=0;
    while(cont++<Nbranches){
      //printf("cont=%d ",cont);
      int Nnodes=nT.n_.size();
      nT.insert_nodes(
        Tinit,
        min_segment_length, //min_segment_length*0.75,
        max_segment_length, //max_segment_length*0.75,
        distance_segment_factor, //cont<(1+Nbranches/2)?0.75*distance_segment_factor:distance_segment_factor,
        min_radius,max_radius,
        distance_node_factor, //cont<(1+Nbranches/2)?0.75*distance_node_factor:distance_node_factor,
        smoothing_factor,
        type,
        gamma,
        1,// Nbranches,
        angles,
        MaxNodes,
        MaxNeigbors,
        allow_branche_intersection
      );
      if(Nnodes==(int) nT.n_.size()) break;
    }


   min_radius=40.*zoom*average_radius/50.;max_radius=average_radius*zoom;
   distance_segment_factor=2./pow(node_radius_reduction_factor,6./zoom);
   distance_node_factor=3./pow(node_radius_reduction_factor,6./zoom);

   float smy=myrand(1.);
   distance_node_factor*=pow(node_radius_reduction_factor,1+5.*smy);
   distance_segment_factor*=pow(node_radius_reduction_factor,0.5+3.*smy);

   Nbranches=100000;
   if(allow_branche_intersection==true){
      min_segment_length=100.*zoom*pow(node_radius_reduction_factor,6./zoom);
      max_segment_length=250.*zoom*pow(node_radius_reduction_factor,6./zoom);
      distance_node_factor=3./pow(node_radius_reduction_factor,6./zoom);
    }
   gamma=4.;
   angles=nT.angle_range(13);
   if(MaxAngles<(int) angles.size()) angles.resize(MaxAngles);

    if(allow_branche_intersection==true){
      distance_node_factor*=2.*pow(node_radius_reduction_factor,4);
    }

    /// INCLUSION OF NODES IN THE TREE
    cont=0;
    while(true){
      cont++;
      int Nnodes=nT.n_.size();


      float gamma2=gamma;
      bool allow_branche_intersection2=false;
      float min_segment_length2=min_segment_length;
      float max_segment_length2=max_segment_length;
      float distance_node_factor2=distance_node_factor;
      if(cont%2==0 && cont>5 && allow_branche_intersection==true){
        gamma2=1.;
        allow_branche_intersection2=true;
        min_segment_length2=max_segment_length;
        max_segment_length2*=1.3;
        distance_node_factor2*=2*pow(node_radius_reduction_factor,4);
      }

      nT.insert_nodes(
        Tinit,
        min_segment_length2,
        max_segment_length2,
        distance_segment_factor,
        min_radius,max_radius,
        distance_node_factor2,
        smoothing_factor,
        type,
        gamma2,
        1,// Nbranches,
        angles,
        MaxNodes,
        3,
        allow_branche_intersection2
      );
      if(Nnodes==(int) nT.n_.size()) break;
      if(cont>5){
        min_segment_length*=node_radius_reduction_factor;
        max_segment_length*=node_radius_reduction_factor;
        min_radius*=node_radius_reduction_factor;
        max_radius*=node_radius_reduction_factor;
        if(min_radius<1. || min_segment_length<2.) break; //min_radius=1;
        distance_segment_factor/=node_radius_reduction_factor_sqrt;
        distance_node_factor/=node_radius_reduction_factor;
      }
    }

   nT.joint_extrema(joint_extrema_probability);

   if(round_termination==false){
      for(unsigned int k=0;k<nT.n_.size();k++){
        if(nT.i_[k].size()==1) {
          nT.r_[k]=1;
          nT.s_[k]=0.1;
          nT.n_[k]=nT.n_[nT.i_[k][0]]+(nT.n_[k]-nT.n_[nT.i_[k][0]])*0.75;
        }
      }
    }

    nT.normalize();

//      char name[400]="Nt.txt";
//      nT.write(name);
  }

    ///-------------------------------------------------------------------------------------------------------
  /// MATISSE 10 BIS SEAWOOD 2
  else if(shape_type==7)
  {
      node_radius_reduction_factor=pow(node_radius_reduction_factor,0.5);

      /// BUILD A DECAGON TO CONSTRAINT THE SHAPE
      float T0y[5]={R0*sinf(2*M_PI/5),R0*sinf(M_PI/5),0,-R0*sinf(M_PI/5),-R0*sinf(2*M_PI/5)};
      float T0R[5]={R0*cosf(2*M_PI/5),R0*cosf(M_PI/5),R0,R0*cosf(M_PI/5),R0*cosf(2*M_PI/5)};

      if(Tinit.p.size()>2){
        vector<point2d> eV=Tinit.extrema();
        float scalex=(eV[1].x-eV[0].x);
        float scaley=(eV[1].y-eV[0].y);
        float scale=scalex/scaley;
        for(int k=0;k<5;k++) T0R[k]*=scale;
      }

      Tinit=triangulation(point2d(512.,512.),T0y,T0R);

      nT.clear();
      nT.n_.push_back(point2d(512.,512.-0.9*R0));
      nT.i_.resize(1);
      nT.r_.push_back(zoom*average_radius);
      nT.s_.push_back(smoothing_factor);

      float min_segment_length=200.*zoom,max_segment_length=200.*zoom;
      float distance_segment_factor=2.;
      float min_radius=40.*zoom*average_radius/50.,max_radius=40*zoom*average_radius/50.;
      float distance_node_factor=4.;

      gamma=2.;
      vector<double> angles=nT.angle_range(7);
      if(MaxAngles<(int) angles.size()) angles.resize(MaxAngles);

      vector<int> posR(1,0);
      vector<int> posL(1,0);
      int direction=1.;
      while(true){
        int pos2=direction==1?posR[posR.size()-1]:posL[posL.size()-1];
        point2d  p;
        if(pos2==0) p=nT.n_[pos2]+point2d(zoom*average_radius*(2.+0.5*myrand(1.))*direction,(0.1+0.6*myrand(1.))*30.*zoom*average_radius/50.);
        else p=nT.n_[pos2]+point2d(zoom*average_radius*(2.+0.5*myrand(1.))*direction,(0.2+1.2*myrand(1.))*30.*zoom*average_radius/50.);
        if(Tinit.inclusion_test(p)<0){
          if(direction==-1) break;
          direction=-1;
          p=nT.n_[0]+point2d(zoom*average_radius*(2.+0.5*myrand(1.))*direction,(0.1+0.6*myrand(1.))*30.*zoom*average_radius/50.);
        }
        if(direction==1){
          nT.insert_point(p,40.*zoom*average_radius/50.,0.75,posR[posR.size()-1]);
          posR.push_back(nT.n_.size()-1);
        }
        else{
          nT.insert_point(p,40.*zoom*average_radius/50.,0.75,posL[posL.size()-1]);
          posL.push_back(nT.n_.size()-1);
        }

      }
      min_segment_length*=node_radius_reduction_factor;
      max_segment_length*=node_radius_reduction_factor;
      min_radius*=node_radius_reduction_factor;
      max_radius*=node_radius_reduction_factor;

      nT.insert_nodes(
          Tinit,
          min_segment_length,max_segment_length,
          distance_segment_factor,
          min_radius,max_radius,
          distance_node_factor,
          smoothing_factor,
          type,
          gamma,
          1,// Nbranches,
          angles,
          MaxNodes,
          3,
          allow_branche_intersection
        );


      float length=0;
      while(true){
        length+=max_segment_length;
        if(length>600.) break;
        int Nnodes=nT.n_.size();
        min_segment_length*=node_radius_reduction_factor;
        max_segment_length*=node_radius_reduction_factor;
        min_radius*=node_radius_reduction_factor;
        max_radius*=node_radius_reduction_factor;
        if(min_radius<1. || min_segment_length<2.) break; //min_radius=1;
        distance_segment_factor/=node_radius_reduction_factor_sqrt;
        distance_node_factor/=node_radius_reduction_factor;

        nT.insert_nodes(
          Tinit,
          min_segment_length,max_segment_length,
          distance_segment_factor,
          min_radius,max_radius,
          distance_node_factor,
          smoothing_factor,
          type,
          gamma,
          1,// Nbranches,
          angles,
          MaxNodes,
          2,
          allow_branche_intersection
        );
        if(Nnodes==(int) nT.n_.size()) break;
      }

      distance_node_factor*=3./2.;
      distance_segment_factor*=4./3.;
      float smy=myrand(1.);
      //smy=0.99;
      distance_node_factor*=pow(node_radius_reduction_factor,1+5.*smy);
      distance_segment_factor*=pow(node_radius_reduction_factor,0.5+3.*smy);

    angles=nT.angle_range(15);
    gamma=4.;
    if(MaxAngles<(int) angles.size()) angles.resize(MaxAngles);


      if(allow_branche_intersection==true){
        distance_node_factor*=2.*pow(node_radius_reduction_factor,4);
      }
      /// INCLUSION OF NODES IN THE TREE
      int cont=0;
      while(true){
        cont++;
        int Nnodes=nT.n_.size();
        min_segment_length*=node_radius_reduction_factor;
        max_segment_length*=node_radius_reduction_factor;
        min_radius*=node_radius_reduction_factor;
        max_radius*=node_radius_reduction_factor;
        if(min_radius<1. || min_segment_length<2.) break; //min_radius=1;
        distance_segment_factor/=node_radius_reduction_factor_sqrt;
        distance_node_factor/=node_radius_reduction_factor;

        float gamma2=gamma;
        bool allow_branche_intersection2=false;
        float min_segment_length2=min_segment_length;
        float max_segment_length2=max_segment_length;
        float distance_node_factor2=distance_node_factor;
        if(cont%2==0 && allow_branche_intersection==true){
          gamma2=1.;
          allow_branche_intersection2=true;
          min_segment_length2=max_segment_length;
          max_segment_length2*=1.3;
          distance_node_factor2*=2*pow(node_radius_reduction_factor,4);
        }

        nT.insert_nodes(
          Tinit,
          min_segment_length2,
          max_segment_length2,
          distance_segment_factor,
          min_radius,max_radius,
          distance_node_factor2,
          smoothing_factor,
          1,//type,
          gamma2,
          1,// Nbranches,
          angles,
          MaxNodes,
          4,
          allow_branche_intersection2
        );
        if(Nnodes==(int) nT.n_.size()) break;
      }



      for(unsigned int k=posL.size()-1;k>1;k--){
        if(nT.i_[posL[k]].size()==1) nT.erase_point(posL[k]);
      }
      for(unsigned int k=posR.size()-1;k>1;k--){
        if(nT.i_[posR[k]].size()==1) nT.erase_point(posR[k]);
      }

     nT.joint_extrema(joint_extrema_probability);

     if(round_termination==false){
        for(unsigned int k=0;k<nT.n_.size();k++){
          if(nT.i_[k].size()==1) {
            nT.r_[k]=1;
            nT.s_[k]=0.1;
            nT.n_[k]=nT.n_[nT.i_[k][0]]+(nT.n_[k]-nT.n_[nT.i_[k][0]])*0.75;
          }
        }
      }
      nT.normalize();
  }

  ///-------------------------------------------------------------------------------------------------------
  /// MATISSE 14 (SYMMETRY)
  else if(shape_type==8)
  {
    int SymmetryType=1+rand()%3;
    /// BUILD A DECAGON TO CONSTRAINT THE SHAPE
    float T0y[5]={R0*sinf(2*M_PI/5),R0*sinf(M_PI/5),0,-R0*sinf(M_PI/5),-R0*sinf(2*M_PI/5)};
    float T0R[5]={R0*cosf(2*M_PI/5),R0*cosf(M_PI/5),R0,R0*cosf(M_PI/5),R0*cosf(2*M_PI/5)};
    Tinit=triangulation(point2d(512.,512.),T0y,T0R);

    float min_segment_length=200.*zoom,max_segment_length=250.*zoom;
    while(true){
      point2d p(512.+(2.*myrand(1.)-1)*R0*0.5,512.+(2.*myrand(1.)-1)*R0*0.5);
      if(nT.InsideSymmetricRegion(p.x,p.y,512.,512.,SymmetryType)==false) continue;
      nT.point_line_projection(p,point2d(512.,512.),min_segment_length*0.5,SymmetryType);
      nT.n_.push_back(p);
      nT.i_.resize(1);
      nT.r_.push_back(50*1.5*zoom*average_radius/50.);
      nT.s_.push_back(smoothing_factor);
      break;
    }


    float distance_segment_factor=2.5;
    float min_radius=average_radius*zoom,max_radius=60*zoom*average_radius/50.;
    float distance_node_factor=5.;

    vector<double> angles=nT.angle_range(100);

    if(min_segment_length<3.5*max_radius){
        double temp=max_radius;
        max_radius=min_segment_length/3.5;
        min_radius*=max_radius/temp;
        nT.r_[0]=max_radius;
    }


    if(allow_branche_intersection==true){
      min_segment_length=200.*zoom;
      max_segment_length=250.*zoom;
      distance_node_factor=4.;
    }

    nT.insert_nodes(
      Tinit,
      min_segment_length,max_segment_length,
      distance_segment_factor,
      min_radius,max_radius,
      distance_node_factor,
      smoothing_factor,
      type,
      gamma,
      Nbranches,
      angles,
      MaxNodes,
      MaxNeigbors,
      allow_branche_intersection,
      SymmetryType
   );

    //nT.print();
    if(myrand(1.)<0.5){
      for(unsigned int k=0;k<nT.n_.size();k++){
        if(nT.i_[k].size()==1) {
          nT.r_[k]=1;
          nT.s_[k]=0.1;
          nT.n_[k]=nT.n_[nT.i_[k][0]]+(nT.n_[k]-nT.n_[nT.i_[k][0]])*0.75;
        }
     }
    }


   if(round_termination==false){
      for(unsigned int k=0;k<nT.n_.size();k++){
        if(nT.i_[k].size()==1) {
          nT.r_[k]=1;
          nT.s_[k]=0.1;
          nT.n_[k]=nT.n_[nT.i_[k][0]]+(nT.n_[k]-nT.n_[nT.i_[k][0]])*0.75;
        }
      }
    }
    nT.normalize();
  }
    ///-------------------------------------------------------------------------------------------------------
  /// MATISSE 15 PERIODICITY
  else if(shape_type==9)
  {
    int SymmetryType=4+rand()%3;
    zoom=average_node_distance*0.8; //0.1+myrand(1.)*0.9;
    dmin=128.;

    Tinit=triangulation(point2d(512.-dmin,512.-dmin), point2d(512.+dmin,512.+dmin));

    float min_segment_length=200.*zoom,max_segment_length=250.*zoom;
    if(smoothing_factor==0) max_segment_length=min_segment_length;
    while(true){
      point2d p(512.+dmin*myrand(1.),512.+dmin*myrand(1.));
      if(smoothing_factor==0.) p.y=p.x;
      if(nT.InsideSymmetricRegion(p.x,p.y,512.,512.,SymmetryType)==false) continue;
      if(nT.point_line_projection(p,point2d(512.,512.),min_segment_length*0.5,SymmetryType)==true) continue;
      nT.n_.push_back(p);
      nT.i_.resize(1);
      nT.r_.push_back(zoom*average_radius);
      nT.s_.push_back(smoothing_factor);
      break;
    }


    float distance_segment_factor=2.5;
    float min_radius=average_radius*zoom,max_radius=average_radius*zoom;
    float distance_node_factor=5.;

    vector<double> angles=nT.angle_range(100);
    if(smoothing_factor==0){
      angles.clear();
      angles.push_back(0.);
      angles.push_back(M_PI/2);
      angles.push_back(M_PI);
      angles.push_back(3*M_PI/2);
      if(myrand(1.)<0.5){
        angles.push_back(M_PI/4);
        angles.push_back(3*M_PI/4);
        angles.push_back(5*M_PI/4);
        angles.push_back(7*M_PI/4);
      }
    }
    if(allow_branche_intersection==true){
      //printf("allow_branche_intersection==true\n");
      min_segment_length=200.*zoom;
      max_segment_length=250.*zoom;
      distance_node_factor=4.;
    }
//printf("min_segment_length=%lf,max_radius=%lf ratio=%lf min_radius=%lf\n",min_segment_length,max_radius,min_segment_length/max_radius,min_radius);
    if(min_segment_length<4*max_radius){
        double temp=max_radius;
        max_radius=min_segment_length/4;
        min_radius*=max_radius/temp;
        nT.r_[0]=max_radius;
    }
//printf("min_segment_length=%lf,max_radius=%lf ratio=%lf min_radius=%lf\n",min_segment_length,max_radius,min_segment_length/max_radius,min_radius);
    nT.insert_nodes(
      Tinit,
      min_segment_length,max_segment_length,
      distance_segment_factor,
      min_radius,max_radius,
      distance_node_factor,
      smoothing_factor,
      type,
      gamma,
      Nbranches,
      angles,
      MaxNodes,
      MaxNeigbors,
      allow_branche_intersection,
      SymmetryType
   );

  for(unsigned int k=0;k<nT.n_.size();k++){
    nT.n_[k].x=512.+(nT.n_[k].x-512.)*4.;
    nT.n_[k].y=512.+(nT.n_[k].y-512.)*4.;
    nT.r_[k]=nT.r_[k]*4.;
  }

   //nT.joint_extrema(joint_extrema_probability);

   if(round_termination==false){
      for(unsigned int k=0;k<nT.n_.size();k++){
        if(nT.i_[k].size()==1) {
          nT.r_[k]=1;
          nT.s_[k]=0.1;
          nT.n_[k]=nT.n_[nT.i_[k][0]]+(nT.n_[k]-nT.n_[nT.i_[k][0]])*0.75;
        }
      }
    }
  }
  ///-------------------------------------------------------------------------------------------------------
  /// SPIRAL
  else if(shape_type==10)
  {

    average_radius*=1.5;

    float step=16.;
    /// BUILD A DECAGON TO CONSTRAINT THE SHAPE
    float T0y[5]={R0*sinf(2*M_PI/5),R0*sinf(M_PI/5),0,-R0*sinf(M_PI/5),-R0*sinf(2*M_PI/5)};
    float T0R[5]={R0*cosf(2*M_PI/5),R0*cosf(M_PI/5),R0,R0*cosf(M_PI/5),R0*cosf(2*M_PI/5)};
    Tinit=triangulation(point2d(512.,512.),T0y,T0R);

    float r0=75.*zoom+(300-75.*zoom)*myrand(2.);  //25.;//2.*200.*zoom;
    while(average_radius*zoom>r0){
      r0+=1.;
    }

    point2d p0(512.+r0,512.);
    nT.n_.push_back(p0);
    nT.i_.resize(1);
    nT.r_.push_back(average_radius*zoom);
    //nT.s_.push_back(smoothing_factor);
    nT.s_.push_back(0.);


    int k;
    float zoom0=zoom;
    node_radius_reduction_factor=sqrt(node_radius_reduction_factor);
    for(k=1;k<=30000;k++){
      r0+=200.*zoom0/step;

      //if(k==1 || k==3) continue;
      //point2d p(512.+(1.+k/16.)*200.*zoom*cos(2.*M_PI*k/16.),512.+(1.+k/16.)*200.*zoom*sin(2.*M_PI*k/16.));
       point2d p(512.+r0*cos(2.*M_PI*k/step),512.+r0*sin(2.*M_PI*k/step));
      if( Tinit.inclusion_test(p)<0) break;
      zoom*=node_radius_reduction_factor;
      float radius=average_radius*zoom;
      if(radius<1.) break;
      nT.insert_point(p,radius,smoothing_factor,nT.n_.size()-1);
    }
    nT.s_[nT.s_.size()-1]=0.01;
    r0-=1.5*200.*zoom0/step;
    point2d p(512.+r0*cos(2.*M_PI*(k-1.5)/step),512.+r0*sin(2.*M_PI*(k-1.5)/step));
    nT.insert_point(p,average_radius*zoom,smoothing_factor,nT.n_.size()-2,nT.n_.size()-1);

   if(round_termination==false){
      for(unsigned int k=0;k<nT.n_.size();k++){
        if(nT.i_[k].size()==1) {
          nT.r_[k]=1;
          nT.s_[k]=0.1;
          nT.n_[k]=nT.n_[nT.i_[k][0]]+(nT.n_[k]-nT.n_[nT.i_[k][0]])*0.75;
        }
      }
    }
    nT.normalize();
  }
  //Tinit=Tinit0; //EEE
  average_node_distance=average_node_distance0;
  average_radius=average_radius0;
  if(average_radius != 50. || MaxAngles!=32 || average_node_distance>91){
   if(nT.n_.size()==1 || ( (shape_type==8 || shape_type==9) && nT.n_.size()<5 ) ){
     average_node_distance = 50.+40*myrand(1.);
     average_radius = 50;
     MaxAngles=32;
     nT=BranchedShapeGenerator(
        Tinit,
        shape_type,
        random_seed,
        average_node_distance,
        average_radius,
        smoothing_factor,
        branched_finished_in_peak,
        joint_extrema_probability,
        MaxNeigbors,
        allow_branche_intersection,
        node_radius_reduction_factor,
        MaxAngles,
        MaxNodes
      );
    }
  }
  //printf("nT.n_.size()=%d\n",(int) nT.n_.size());

  //printf("end BranchedShapeGenerator()\n");
  return nT;
}

/// RANDOM BRANCHED SHAPE SELECTION
int get_branched_shape_index(){
  int Rn[11]={5,5,5,3,1,10,5,5,10,10,3};
  vector<int>  sel;

  for(int k=0;k<11;k++){
    for(int n=0;n<Rn[k];n++){
      sel.push_back(k);
    }
  }

  return( sel[rand()%sel.size()]);

}

int get_branched_shape_index2(){
  int Rn[11]={20,20,20,2,2,20,10,10,20,20,2};
  vector<int>  sel;

  for(int k=0;k<11;k++){
    for(int n=0;n<Rn[k];n++){
      sel.push_back(k);
    }
  }

  return( sel[rand()%sel.size()]);

}

/// INITIALIZATION OF NodeTree PARAMETERS AAA
void NodeTree_parameters_initialization(
int argc,
char *argv[],
int &shape_type, ///argv[1]
int &random_seed, ///argv[2]
float &average_node_distance, ///argv[3]
float &average_radius, ///argv[4]
float &smoothing_factor, ///argv[5]
int &branched_finished_in_peak, ///argv[6]
int &Xperiods, ///argv[7]
int &Yperiods, ///argv[8]
float &joint_extrema_probability, ///argv[9]
int &MaxNeigbors, ///argv[10]
bool &allow_branche_intersection, ///argv[11]
float &node_radius_reduction_factor, ///argv[12]
int &MaxAngles, ///argv[13]
unsigned int &MaxNodes,
bool &round_termination){

  argc=-1;
  if(argc>2 && random_seed<0) random_seed=atoi(argv[2]);
  else{
    srand (time(NULL)); rand();
    random_seed=rand();
  }

  srand (random_seed); rand();



  /// RANDOM SELECTION OF STYLE
  //shape_type=-1;

  if(shape_type<0){
    /// IF THERE EXISTS A SILHOUETTE WE USE STYLE 0,1,2
    if(argc==17){
      shape_type=rand()%3;
    }
    else{
      /// RANDOM SELECTION OF STYLE
      vector<int> sV;
      for(int k=0;k<5;k++) sV.push_back(0);
      for(int k=0;k<5;k++) sV.push_back(1);
      for(int k=0;k<5;k++) sV.push_back(2);
      for(int k=0;k<2;k++) sV.push_back(3);
      for(int k=0;k<1;k++) sV.push_back(4);
      for(int k=0;k<10;k++) sV.push_back(5);
      for(int k=0;k<5;k++) sV.push_back(6);
      for(int k=0;k<5;k++) sV.push_back(7);
      for(int k=0;k<10;k++) sV.push_back(8);
      for(int k=0;k<10;k++) sV.push_back(9);
      for(int k=0;k<2;k++) sV.push_back(10);

      shape_type = sV[rand()%sV.size()];
    }
  }

  average_node_distance=-1;
  average_radius=-1;
  smoothing_factor=-1;
  Xperiods=-1;
  Yperiods=-1;
  joint_extrema_probability=-1;
  MaxNeigbors=-1;
  node_radius_reduction_factor=-1;
  MaxAngles=-1;
  MaxNodes=0;



  if(argc>3) average_node_distance=atof(argv[3]);
  if(average_node_distance<0) average_node_distance=round(40.+myrand(1.5)*160.);

  if(argc>4) average_radius=atof(argv[4]);
  if(average_radius<0) average_radius=round(3.+myrand(1.)*96);

  if(argc>5){
    smoothing_factor=atof(argv[5]);
    if(smoothing_factor<0.) smoothing_factor=round(2.*myrand(1.)*10.)/10.;
    if(atof(argv[5])<0 && shape_type==1 && myrand(1.)<0.5) smoothing_factor=0.;
  }
  else{
    if(shape_type==1 && myrand(1.)<0.5) smoothing_factor=0.;
    else smoothing_factor=round(2.*myrand(1.)*10.)/10.;
  }

  int temp=-1;
  if(argc>6) temp=atoi(argv[6]);
  branched_finished_in_peak = temp<0?(bool) (rand()%2):(bool) temp ;

  if(argc>7) Xperiods=atoi(argv[7]);
  if(Xperiods<=0){
    if(shape_type==9) Xperiods=3; //2+rand()%4;
    else Xperiods=1;
  }

  if(argc>8) Yperiods=atoi(argv[8]);
  if(Yperiods<=0){
    if(shape_type==9) Yperiods=3; //2+rand()%4;
    else Yperiods=1;
  }

  if(argc>9) joint_extrema_probability=atof(argv[9]);
  if(joint_extrema_probability<0.) joint_extrema_probability=round(myrand(1.)*10.)/10.;

  if(argc>10) MaxNeigbors=atoi(argv[10]);
  if(MaxNeigbors<0.) MaxNeigbors=2+rand()%10;

  temp=-1;
  if(argc>11) temp=atoi(argv[11]);
  allow_branche_intersection = temp<0?(bool) (rand()%2):(bool) temp ; //(bool) (atoi(argv[11])%2);
  allow_branche_intersection = false ; //(bool) (atoi(argv[11])%2);

  if(argc>12) node_radius_reduction_factor=atof(argv[12]);
  if(node_radius_reduction_factor<0.) node_radius_reduction_factor=round(80.+20.*myrand(0.5))/100.;

  if(argc>13) MaxAngles=atoi(argv[13]);
  if(MaxAngles<0) MaxAngles=3+rand()%32;

  if(argc>14) MaxNodes=atoi(argv[14]);

  if(MaxNodes==0){
    MaxNodes=(10+rand()%90);
    if(MaxNodes>100) MaxNodes/=(Yperiods*Xperiods);
    if(allow_branche_intersection==true) MaxNodes/=3;
    //if(MaxNodes<50) MaxNodes=50;
  }
  //printf("MaxNodes=%d\n",(int) MaxNodes);

  /************* INIT EXTRA CONDITIONS IMPOSED TO PARAMETERS */
  if( (argc<5 || (argc>4 && atof(argv[4])<0)) && shape_type==6 && average_radius>70) average_radius=70;
  if((argc<5 || (argc>4 && atof(argv[4])<0)) && shape_type==5 && average_radius>60) average_radius=60;
  if((argc<5 || (argc>4 && atof(argv[4])<0)) && shape_type==1 && average_radius>60) average_radius=60;
  if((argc<3 || (argc>2 && atof(argv[2])<0)) && shape_type==10 ) smoothing_factor=1;

  /// PERIODIC PATTERN
  if((argc<4 || (argc>3 && atof(argv[3])<0)) && (Yperiods*Xperiods)>1 /*&& (average_node_distance<40 || average_node_distance>200)*/  ){
      float temp=myrand(1.);
      average_node_distance = round(40*pow((double) Yperiods*Xperiods,0.33)*temp+(1.-temp)*140.);
  }
  if((argc<14 || (argc>13 && atof(argv[13])<0))&& shape_type==9) MaxAngles=32;
  if((argc<11 || (argc>10 && atof(argv[10])<0)) && shape_type==9) MaxNeigbors=32;

  round_termination = rand()%2==0?true:false;

  if(shape_type==5) average_node_distance*=2;

  /************* END EXTRA CONDITIONS IMPOSED TO PARAMETERS */

}

/// NodeTree RANDOM GENERATOR AAA
NodeTree NodeTree_random_generator(int shape_type){
    /// PARAMETERS
    int random_seed;
    float average_node_distance;
    float average_radius;
    float smoothing_factor;
    int branched_finished_in_peak;
    int Xperiods;
    int Yperiods;
    float joint_extrema_probability;
    int MaxNeigbors;
    bool allow_branche_intersection;
    float node_radius_reduction_factor;
    int MaxAngles;
    unsigned int MaxNodes;
    bool round_termination;
    int argc2;
    char *argv2[1];

    NodeTree_parameters_initialization(
        argc2,
        argv2,
        shape_type,
        random_seed,
        average_node_distance,
        average_radius,
        smoothing_factor,
        branched_finished_in_peak,
        Xperiods,
        Yperiods,
        joint_extrema_probability,
        MaxNeigbors,
        allow_branche_intersection,
        node_radius_reduction_factor,
        MaxAngles,
        MaxNodes,
        round_termination
    );

//Xperiods=1;
//Yperiods=1;
//branched_finished_in_peak=false;
//round_termination=false;
//shape_type=5;

    triangulation Tinit;


    NodeTree nT=BranchedShapeGenerator(
      Tinit,
      shape_type,
      random_seed,
      average_node_distance,
      average_radius,
      smoothing_factor,
      branched_finished_in_peak,
      joint_extrema_probability,
      MaxNeigbors,
      allow_branche_intersection,
      node_radius_reduction_factor,
      MaxAngles,
      MaxNodes
    );

    if(Xperiods>1 || Yperiods>1){
      NodeTree nT0=nT;
      for(int i=0;i<Xperiods;i++){
        for(int j=0;j<Yperiods;j++){
          //printf("i=%d,j=%d\n",i,j);
          if(i==0 && j==0) continue;
          if(shape_type==9){
            nT.add( nT0.translation(point2d(1024.*i,1024.*j)) );
            continue;
          }
          //printf("average_radius=%lf, average_node_distance=%lf\n",average_radius,average_node_distance);
          int random_seed2=random_seed+i*Yperiods+j;
          NodeTree nT2=BranchedShapeGenerator(
              Tinit,
              shape_type,
              random_seed2,
              average_node_distance,
              average_radius,
              smoothing_factor,
              branched_finished_in_peak,
              joint_extrema_probability,
              MaxNeigbors,
              allow_branche_intersection,
              node_radius_reduction_factor,
              MaxAngles,
              MaxNodes
          );

          printf("MaxNodes=%d i=%d, j=%d nT2.n_.size()=%d\n",(int) MaxNodes, i,j,(int) nT2.n_.size());
          if(nT2.checking()==false){
              printf("Problems computing the Branched shape. \nMaybe the average node distance is too large for this configuration or \nthe ratio between the average width of branches and the average node distance is too large.\n");
              return NodeTree_random_generator();
          }

          nT.add( nT2.translation(point2d(1024.*i,1024.*j)) );
          if(nT.checking()==false){
              printf("Problems nT.add(). \n");
              return NodeTree_random_generator();
          }
        }
      }
    }

    double xmin,ymin,xmax,ymax,border=100;
    nT.extrema(xmin,ymin,xmax,ymax);
    double scalex = (1024-2*border)/(xmax-xmin);
    double scaley = (1024-2*border)/(ymax-ymin);
    double scale = scalex>scaley?scaley:scalex;

    for(int k=0;k<nT.n_.size();k++){
      nT.n_[k].x = border + (nT.n_[k].x-xmin)*scale;
      nT.n_[k].y = 1024-(border + (nT.n_[k].y-ymin)*scale);
      nT.r_[k]*=scale;
    }

    nT.extrema(xmin,ymin,xmax,ymax);
    double dx=512.-(xmax+xmin)/2.;
    double dy=512.-(ymax+ymin)/2.;

    for(int k=0;k<nT.n_.size();k++){
      nT.n_[k].x += dx;
      nT.n_[k].y += dy;
    }


    nT.shape_type =  shape_type ;
    nT.random_seed =  random_seed;
    nT.average_node_distance =  average_node_distance;
    nT.average_radius =  average_radius;
    nT.smoothing_factor =  smoothing_factor;
    nT.branched_finished_in_peak =  branched_finished_in_peak;
    nT.Xperiods =  Xperiods;
    nT.Yperiods =  Yperiods;
    nT.joint_extrema_probability =  joint_extrema_probability;
    nT.MaxNeigbors =  MaxNeigbors;
    nT.allow_branche_intersection =  allow_branche_intersection;
    nT.node_radius_reduction_factor =  node_radius_reduction_factor;
    nT.MaxAngles =  MaxAngles;
    nT.MaxNodes =  MaxNodes;
    nT.round_termination = round_termination;

   // nT.print_params();

    return nT;


  }



/// GENERATION OF A SIMILAR SHAPE
NodeTree NodeTree_similar_generation(NodeTree &nTo){
    srand (time(NULL)); rand();

    //nTo.print_params();

    /// PARAMETERS
    int shape_type=nTo.shape_type;
    int random_seed=rand();
    random_seed=rand();
    srand(random_seed);
    float average_node_distance = nTo.average_node_distance;
    float average_radius = nTo.average_radius;
    float smoothing_factor = nTo.smoothing_factor;
    int branched_finished_in_peak = nTo.branched_finished_in_peak;
    int Xperiods = nTo.Xperiods;
    int Yperiods = nTo.Yperiods;
    float joint_extrema_probability = nTo.joint_extrema_probability;
    int MaxNeigbors = nTo.MaxNeigbors;
    bool allow_branche_intersection = nTo.allow_branche_intersection;
    float node_radius_reduction_factor = nTo.node_radius_reduction_factor;
    int MaxAngles = nTo.MaxAngles;
    unsigned int MaxNodes = nTo.MaxNodes;
    bool round_termination = nTo.round_termination;

    triangulation Tinit;

    NodeTree nT=BranchedShapeGenerator(
      Tinit,
      shape_type,
      random_seed,
      average_node_distance,
      average_radius,
      smoothing_factor,
      branched_finished_in_peak,
      joint_extrema_probability,
      MaxNeigbors,
      allow_branche_intersection,
      node_radius_reduction_factor,
      MaxAngles,
      MaxNodes
    );

    if(Xperiods>1 || Yperiods>1){
      NodeTree nT0=nT;
      for(int i=0;i<Xperiods;i++){
        for(int j=0;j<Yperiods;j++){
          //printf("i=%d,j=%d\n",i,j);
          if(i==0 && j==0) continue;
          if(shape_type==9){
            nT.add( nT0.translation(point2d(1024.*i,1024.*j)) );
            continue;
          }
          //printf("average_radius=%lf, average_node_distance=%lf\n",average_radius,average_node_distance);
          int random_seed2=random_seed+i*Yperiods+j;
          NodeTree nT2=BranchedShapeGenerator(
              Tinit,
              shape_type,
              random_seed2,
              average_node_distance,
              average_radius,
              smoothing_factor,
              branched_finished_in_peak,
              joint_extrema_probability,
              MaxNeigbors,
              allow_branche_intersection,
              node_radius_reduction_factor,
              MaxAngles,
              MaxNodes
          );

          printf("i=%d, j=%d nT2.n_.size()=%d\n",i,j,(int) nT2.n_.size());
          if(nT2.checking()==false){
              printf("Problems computing the Branched shape. \nMaybe the average node distance is too large for this configuration or \nthe ratio between the average width of branches and the average node distance is too large.\n");
              return NodeTree_random_generator();
          }

          nT.add( nT2.translation(point2d(1024.*i,1024.*j)) );
          if(nT.checking()==false){
              printf("Problems nT.add(). \n");
              return NodeTree_random_generator();
          }
        }
      }
    }

    double xmin,ymin,xmax,ymax,border=100;
    nT.extrema(xmin,ymin,xmax,ymax);
    double scalex = (1024-2*border)/(xmax-xmin);
    double scaley = (1024-2*border)/(ymax-ymin);
    double scale = scalex>scaley?scaley:scalex;

    for(int k=0;k<nT.n_.size();k++){
      nT.n_[k].x = border + (nT.n_[k].x-xmin)*scale;
      nT.n_[k].y = 1024-(border + (nT.n_[k].y-ymin)*scale);
      nT.r_[k]*=scale;
    }

    nT.extrema(xmin,ymin,xmax,ymax);
    double dx=512.-(xmax+xmin)/2.;
    double dy=512.-(ymax+ymin)/2.;

    for(int k=0;k<nT.n_.size();k++){
      nT.n_[k].x += dx;
      nT.n_[k].y += dy;
    }


    nT.shape_type =  shape_type ;
    nT.random_seed =  random_seed;
    nT.average_node_distance =  average_node_distance;
    nT.average_radius =  average_radius;
    nT.smoothing_factor =  smoothing_factor;
    nT.branched_finished_in_peak =  branched_finished_in_peak;
    nT.Xperiods =  Xperiods;
    nT.Yperiods =  Yperiods;
    nT.joint_extrema_probability =  joint_extrema_probability;
    nT.MaxNeigbors =  MaxNeigbors;
    nT.allow_branche_intersection =  allow_branche_intersection;
    nT.node_radius_reduction_factor =  node_radius_reduction_factor;
    nT.MaxAngles =  MaxAngles;
    nT.MaxNodes =  MaxNodes;
    nT.round_termination =  round_termination;

    return nT;


  }
