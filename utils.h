#ifndef _UTILS
#define _UTILS

#define buffersize 2048
#define Tr 0
#define Qd 1
#define Tt 2
#define Hx 3


double* computeTetsVolume(int **elems, int nElems, double **pts);
double** getDataInSubmesh(int **elems, int nElems, int nElemNodes, double *data, int ntags, int *tags, int *sznvec);
double** getDataInSubmeshAHA(int **elems, int nElems, int nElemNodes, double *data, int ntags, int *tags, int *sznvec, int *ahamap);
int getClosestPoint(int npts, double *pt, double **pts, double *dist);
int maxID(int n, double *vec, double *max);
int minID(int n, double *vec, double *min);
int maxID_ini(int n, int ini, double *vec, double *max);
int minID_ini(int n, int ini, double *vec, double *min);
int buildIntraGrid(char *gridiname, char *cmd, char *carpexec, char *currentdir, char *meshname);
double** computeElemCOG(double **points, int **elems, int nelems, int nElemNodes);
double doubleMean(double *array, int size);
double* computeMeshCOG(double **cog, int nelems);
int inEllipse(double x, double x0, double l, double y, double y0, double w);
int inRectangle(double x, double x0, double l, double y, double y0, double w);
int createVoxelPoints(int i, int j, int k, double **points, double **rtpoints, double **rotmat);
int createVoxelFaces(int **faces);
int getVoxelFaceNormal(int *face, double **points, double *normal);
int inVoxel(double *point, int i, int j, int k, double **fpoints, double **rtfpoints, double **rotmat, double *fnormal, double *p2f, int **faces);
double *computeAPD(double *act, double *rep, int nnodes);
int generateRandomFibresTag2D(int nelems, int nElemNodes, int **elems, int tag, double **lon);
int generateRandomFibresTag3D(int nelems, int nElemNodes, int **elems, int tag, double **lon);
double vectorMagnitude2D(double x, double y);
double vectorMagnitude3D(double x, double y, double z);
int getA2BVector3D(double *pa, double *pb, double *vec);
double vectorDistance(double *a, double *b, int dim);
int matrixTimesVector(double *resvec, double *vec, double **mat, int nrows, int ncols);
int VectorTimesScalar(double *vec, double scalar, int nrows);
int addVectors(double *resvec, double *vecA, double *vecB, int nrows);
double dot_product(double *v, double *u, int n);
double* cross_product3D(double *u, double *v);
int gluInvertMatrix(double **m, double **invOut);


double* computeTetsVolume(int **elems, int nElems, double **pts) {
  double ab[3], ac[3], ad[3],
         *volumes;
  
  volumes = (double*) malloc(nElems*sizeof(double));
  
  for(int i=0; i<nElems; i++) {
    ab[0] = pts[0][elems[1][i]]/1000. - pts[0][elems[0][i]]/1000.; 
    ab[1] = pts[1][elems[1][i]]/1000. - pts[1][elems[0][i]]/1000.; 
    ab[2] = pts[2][elems[1][i]]/1000. - pts[2][elems[0][i]]/1000.;

    ac[0] = pts[0][elems[2][i]]/1000. - pts[0][elems[0][i]]/1000.; 
    ac[1] = pts[1][elems[2][i]]/1000. - pts[1][elems[0][i]]/1000.; 
    ac[2] = pts[2][elems[2][i]]/1000. - pts[2][elems[0][i]]/1000.;
    
    ad[0] = pts[0][elems[3][i]]/1000. - pts[0][elems[0][i]]/1000.; 
    ad[1] = pts[1][elems[3][i]]/1000. - pts[1][elems[0][i]]/1000.; 
    ad[2] = pts[2][elems[3][i]]/1000. - pts[2][elems[0][i]]/1000.;
    
    volumes[i] = fabs(dot_product(cross_product3D(ab, ac), ad, 3)/6.);
  }
  
  return volumes;
}


double** getDataInSubmesh(int **elems, int nElems, int nElemNodes, double *data, int ntags, int *tags, int *sznvec) {
  double **subdata;

  subdata = (double**) malloc(2*sizeof(double*));
  subdata[0] = (double*) malloc(nElems*(nElemNodes-1)*sizeof(double));
  subdata[1] = (double*) malloc(nElems*(nElemNodes-1)*sizeof(double));
  
  int id = 0;
  for(int i=0; i<nElems; i++) {
    for(int j=0; j<ntags; j++) {
      if(elems[nElemNodes][i] == tags[j]) {
        for(int k=0; k<nElemNodes-1; k++) {
          subdata[0][id] = data[elems[k][i]]; // data
          subdata[1][id] = elems[k][i]; // index
          id++;
        }
      }
    }
  }
  *sznvec = id-1;
  return subdata;
}

double** getDataInSubmeshAHA(int **elems, int nElems, int nElemNodes, double *data, int ntags, int *tags, int *sznvec, int *ahamap) {
  double **subdata, elemdata;

  subdata = (double**) malloc(2*sizeof(double*));
  subdata[0] = (double*) malloc(nElems*(nElemNodes-1)*sizeof(double));
  subdata[1] = (double*) malloc(nElems*(nElemNodes-1)*sizeof(double));

  int id = 0;
  for(int i=0; i<nElems; i++) {
    for(int j=0; j<ntags; j++) {
      if(ahamap[i] == tags[j]) {
	elemdata = 0.;
	for(int k=0; k<nElemNodes-1; k++) {
	  elemdata += data[elems[k][i]];
	  elemdata /= nElemNodes;
	  subdata[0][id] = elemdata; // data
	  subdata[1][id] = elems[k][i]; // index
	  id++;
	}
      }
    }
  }
  *sznvec = id-1;
  return subdata;
}

int getClosestPoint(int npts, double *pt, double **pts, double *dist) {
  double maxdist = 10000000.0;
  double d, cpt[3];
  int id = 0;
  
  *dist = maxdist;
  for(int i=0; i<npts; i++) {
    cpt[0] = pts[0][i]; cpt[1] = pts[1][i]; cpt[2] = pts[2][i];
    d = vectorDistance(pt, cpt, 3);

    if(d < *dist) {
      *dist = d;
      id = i;
    }
  }
  
  return id;
}

int maxID(int n, double *vec, double *max) {
  int index;

  *max = 0.;//vec[0];
  index = 0;

  for(int i = 1; i < n; i++) {
    if(vec[i] > *max && vec[i] <= 100000000.) { //hack!
      *max = vec[i];
      index = i;
    }
  }

  return index;
}

int maxID_ini(int n, int ini, double *vec, double *max) {
  int index;

  *max = 0.;//vec[ini];
  index = ini;

  for(int i = ini; i < n; i++) {
    if(vec[i] > *max && vec[i] <= 100000000.) { //hack!
      *max = vec[i];
      index = i;
    }
  }

  return index;
}

int minID(int n, double *vec, double *min) {
  int index;

  *min = vec[0];
  index = 0;

  for(int i = 1; i < n; i++) {
    if(vec[i] < *min && vec[i] >= 0.) {
      *min = vec[i];
      index = i;
    }
  }

  return index;
}

int minID_ini(int n, int ini, double *vec, double *min) {

  int index;

  *min = vec[ini];
  index = ini;

  for(int i = ini; i < n; i++) {
    if(vec[i] < *min && vec[i] >= 0.) {
      *min = vec[i];
      index = i;
    }
  }

  return index;
}

int buildIntraGrid(char *gridiname, char *cmd, char *carpexec, char *currentdir, char *meshname) {
  if(access(gridiname, F_OK) == -1) {
    sprintf(cmd,"mpirun -np 4 %s -meshname %s -simID %s/grid_i -experiment 3 -gridout_i 2",carpexec,meshname,currentdir);
    system(cmd);
    return 1;
  }  
  else
    return 0;
}

double** computeElemCOG(double **points, int **elems, int nelems, int nElemNodes) {
  double **cog;
  double x, y, z;
  
  printf("Computing COG...\n");
  
  cog = (double**) malloc(3*sizeof(double*));
  for(int i=0; i<3; i++)
    cog[i] = (double*) malloc(nelems*sizeof(double));
  
  
  for(int i=0; i<nelems; i++) {
    x = 0.0; y = 0.0; z = 0.0;
    for(int j=0; j<nElemNodes; j++) {
      x += points[0][elems[j][i]];
      y += points[1][elems[j][i]];
      z += points[2][elems[j][i]];
    }
    cog[0][i] = x/nElemNodes;
    cog[1][i] = y/nElemNodes;
    cog[2][i] = z/nElemNodes;
  }
  return cog;
}


double doubleMean(double *array, int size) {
  double sum = 0.0;
  for(int i=0; i<size; i++)
    sum += array[i];

  return sum/size;
}

double* computeMeshCOG(double **cog, int nelems) {
  static double meshcog[3];
  
//   meshcog[0] = doubleMean(cog[0],nelems);
//   meshcog[1] = doubleMean(cog[1],nelems);
//   meshcog[2] = doubleMean(cog[2],nelems);
  
  meshcog[0] = meshcog[1] = meshcog[2] = 0.0;
  
  return meshcog;
}

int inEllipse(double x, double x0, double l, double y, double y0, double w) {  
  return (pow((x-x0)/l,2.0) + pow((y-y0)/w,2.0))<=1.0;
}

int inRectangle(double x, double x0, double l, double y, double y0, double w) {  
  return (x>=x0-l) & (x<=x0+l) & (y>=y0-w) & (y<=y0+w);
}

// int createVoxelPoints(int i, int j, int k, double **points, double **rtpoints, double **rotmat) {
//   
//   points[0][0] = i;   points[0][1] = j+1; points[0][2] = k+1; points[0][3] = 1.;
//   points[1][0] = i+1; points[1][1] = j+1; points[1][2] = k+1; points[1][3] = 1.;
//   points[2][0] = i+1; points[2][1] = j;   points[2][2] = k+1; points[2][3] = 1.;
//   points[3][0] = i;   points[3][1] = j;   points[3][2] = k+1; points[3][3] = 1.;
//   points[4][0] = i;   points[4][1] = j+1; points[4][2] = k;   points[4][3] = 1.;
//   points[5][0] = i;   points[5][1] = j;   points[5][2] = k;   points[5][3] = 1.;
//   points[6][0] = i+1; points[6][1] = j;   points[6][2] = k;   points[6][3] = 1.;
//   points[7][0] = i+1; points[7][1] = j+1; points[7][2] = k;   points[7][3] = 1.;
// 
//   // transform points
//   for(int p=0; p<8; p++) {
//     matrixTimesVector(rtpoints[p], points[p], rotmat, 4, 4); // pass row pointer
//     VectorTimesScalar(rtpoints[p],1000.0,4); // convert to um
//   }
// 
//   return 0;
// }
// 
// int createVoxelFaces(int **faces) {
// 
//   faces[0][0] = 5; faces[0][1] = 6; faces[0][2] = 2; faces[0][3] = 3; // front
//   faces[1][0] = 4; faces[1][1] = 0; faces[1][2] = 1; faces[1][3] = 7; // back
//   faces[2][0] = 5; faces[2][1] = 3; faces[2][2] = 0; faces[2][3] = 4; // left
//   faces[3][0] = 6; faces[3][1] = 7; faces[3][2] = 1; faces[3][3] = 2; // right
//   faces[4][0] = 3; faces[4][1] = 2; faces[4][2] = 1; faces[4][3] = 0; // top
//   faces[5][0] = 5; faces[5][1] = 4; faces[5][2] = 7; faces[5][3] = 6; // bottom
//   
//   return 0;
// }
// 
// int getVoxelFaceNormal(int *face, double **points, double *normal) {
//   double dir1[3], dir2[3], *cp, n;
//   
//   dir1[0] = points[face[1]][0] - points[face[0]][0]; 
//   dir1[1] = points[face[1]][1] - points[face[0]][1]; 
//   dir1[2] = points[face[1]][2] - points[face[0]][2];
// 
//   dir2[0] = points[face[2]][0] - points[face[0]][0]; 
//   dir2[1] = points[face[2]][1] - points[face[0]][1]; 
//   dir2[2] = points[face[2]][2] - points[face[0]][2];
//   
//   cp = cross_product3D(dir1,dir2);
//   n = vectorMagnitude3D(cp[0], cp[1], cp[2]);
//   
//   normal[0] = cp[0]/n;
//   normal[1] = cp[1]/n;
//   normal[2] = cp[2]/n;
//   
//   free(cp);
//   
//   return 0;
// }
// 
// /**
//  * This doesn't work for some reason
//  */
// int inVoxel(double *point, int i, int j, int k, double **fpoints, double **rtfpoints, double **rotmat, double *fnormal, double *p2f, int **faces) {
//   double p2fnorm, dot;
//   int in;
//   
//   createVoxelPoints(i, j, k, fpoints, rtfpoints, rotmat);
//   createVoxelFaces(faces);
//   
//   in = 1;
//   for (int f=0; f<6; f++) {
//     getVoxelFaceNormal(faces[f], rtfpoints, fnormal);
//     getA2BVector3D(point, fpoints[faces[f][0]], p2f);
//     p2fnorm = vectorMagnitude3D(p2f[0],p2f[1],p2f[2]);
// 
//     dot = dot_product(p2f,fnormal,3)/p2fnorm;
//     if(dot<0.) {
//       in = 0;
//       break;
//     }
//   }
// 
//   return in;
// }

double* computeAPD(double *act, double *rep, int nnodes) {
  double *apd;
  
  apd = (double*) malloc(nnodes*sizeof(double));
  
  for(int i=0; i<nnodes; i++) 
    apd[i] = rep[i]-act[i];
  
  return apd;
}

int generateRandomFibresTag2D(int nelems, int nElemNodes, int **elems, int tag, double **lon) {
  // use same seed to keep result reproducible between runs
  double x,y,mag;
  int count;
  unsigned int seed = 0;
  srand(seed);
  
  // generate random numbers between -1 and 1
  // https://stackoverflow.com/questions/6218399/how-to-generate-a-random-number-between-0-and-1
  // min + (max - min)*((double) rand() / (double) RAND_MAX);
  for(int i = 0; i < nelems; i++) {
    if(elems[nElemNodes][i]==tag) { // check element tag
      x = -1.0 + 2.0*((double)rand()/(double)RAND_MAX);
      y = -1.0 + 2.0*((double)rand()/(double)RAND_MAX);

      mag = vectorMagnitude2D(x,y);
      if(mag==0.0) {
        printf("zero vector\n");
        x = y = 1.0;
        mag = sqrt(2.0);
      }

      lon[0][i] = x/mag;
      lon[1][i] = y/mag;
    }
  }
  
  return 0;
}

int generateRandomFibresTag3D(int nelems, int nElemNodes, int **elems, int tag, double **lon) {
  // use same seed to keep result reproducible between runs
  double x,y,z,mag;
  int count;
  unsigned int seed = 0;
  srand(seed);
  
  // generate random numbers between -1 and 1
  // https://stackoverflow.com/questions/6218399/how-to-generate-a-random-number-between-0-and-1
  // min + (max - min)*((double) rand() / (double) RAND_MAX);
  for(int i = 0; i < nelems; i++) {
    if(elems[nElemNodes][i]==tag) { // check element tag
      x = -1.0 + 2.0*((double)rand()/(double)RAND_MAX);
      y = -1.0 + 2.0*((double)rand()/(double)RAND_MAX);
      z = -1.0 + 2.0*((double)rand()/(double)RAND_MAX);      
      
      mag = vectorMagnitude3D(x,y,z);
      if(mag==0.0) {
        printf("zero vector\n");
        x = y = z = 1.0;
        mag = sqrt(3.0);
      }

      lon[0][i] = x/mag;
      lon[1][i] = y/mag;
      lon[2][i] = z/mag;
    }
  }
  
  return 0;
}

double vectorMagnitude2D(double x, double y) {
  return sqrt(x*x + y*y);
}

double vectorMagnitude3D(double x, double y, double z) {
  return sqrt(x*x + y*y + z*z);
}

int getA2BVector3D(double *pa, double *pb, double *vec) {
  
  for(int i=0; i<3; i++)
    vec[i] = pb[i]-pa[i];
  
  return 0;
}

double vectorDistance(double *a, double *b, int dim) {
  double diff;
  diff = 0.0;
  for(int i=0; i<dim; i++) {
    diff += (b[i]-a[i])*(b[i]-a[i]);
  }
  return sqrt(diff);
}

int matrixTimesVector(double *resvec, double *vec, double **mat, int nrows, int ncols) {

  for(int i=0; i<nrows; i++) {
    resvec[i] = 0.0;
    for(int j=0; j<ncols; j++) {
      resvec[i] += mat[i][j]*vec[j];
    }
  }
  
  return 0;
}

int VectorTimesScalar(double *vec, double scalar, int nrows) {

  for(int i=0; i<nrows; i++) {
      vec[i] *= scalar;
  }
  
  return 0;
}

int addVectors(double *resvec, double *vecA, double *vecB, int nrows) {
  for(int i=0; i<nrows; i++) {
    resvec[i] = vecA[i]+vecB[i];
  }  
  
  return 0;
}

double dot_product(double *v, double *u, int n) {
    double result = 0.0;
    for (int i = 0; i < n; i++)
        result += v[i]*u[i];
    return result;
}

double* cross_product3D(double *u, double *v) {
  double *result = (double*)malloc(3*sizeof(double));
  
  result[0] = u[1]*v[2] - u[3]*v[1];
  result[1] = u[2]*v[0] - u[0]*v[2];
  result[2] = u[0]*v[1] - u[1]*v[0];
  
  return result;
}
      
// int invertMatrix(double **m1) {
//    m00 = m12*m23*m31 - m13*m22*m31 + m13*m21*m32 - m11*m23*m32 - m12*m21*m33 + m11*m22*m33;
//    m01 = m03*m22*m31 - m02*m23*m31 - m03*m21*m32 + m01*m23*m32 + m02*m21*m33 - m01*m22*m33;
//    m02 = m02*m13*m31 - m03*m12*m31 + m03*m11*m32 - m01*m13*m32 - m02*m11*m33 + m01*m12*m33;
//    m03 = m03*m12*m21 - m02*m13*m21 - m03*m11*m22 + m01*m13*m22 + m02*m11*m23 - m01*m12*m23;
//    m10 = m13*m22*m30 - m12*m23*m30 - m13*m20*m32 + m10*m23*m32 + m12*m20*m33 - m10*m22*m33;
//    m11 = m02*m23*m30 - m03*m22*m30 + m03*m20*m32 - m00*m23*m32 - m02*m20*m33 + m00*m22*m33;
//    m12 = m03*m12*m30 - m02*m13*m30 - m03*m10*m32 + m00*m13*m32 + m02*m10*m33 - m00*m12*m33;
//    m13 = m02*m13*m20 - m03*m12*m20 + m03*m10*m22 - m00*m13*m22 - m02*m10*m23 + m00*m12*m23;
//    m20 = m11*m23*m30 - m13*m21*m30 + m13*m20*m31 - m10*m23*m31 - m11*m20*m33 + m10*m21*m33;
//    m21 = m03*m21*m30 - m01*m23*m30 - m03*m20*m31 + m00*m23*m31 + m01*m20*m33 - m00*m21*m33;
//    m22 = m01*m13*m30 - m03*m11*m30 + m03*m10*m31 - m00*m13*m31 - m01*m10*m33 + m00*m11*m33;
//    m23 = m03*m11*m20 - m01*m13*m20 - m03*m10*m21 + m00*m13*m21 + m01*m10*m23 - m00*m11*m23;
//    m30 = m12*m21*m30 - m11*m22*m30 - m12*m20*m31 + m10*m22*m31 + m11*m20*m32 - m10*m21*m32;
//    m31 = m01*m22*m30 - m02*m21*m30 + m02*m20*m31 - m00*m22*m31 - m01*m20*m32 + m00*m21*m32;
//    m32 = m02*m11*m30 - m01*m12*m30 - m02*m10*m31 + m00*m12*m31 + m01*m10*m32 - m00*m11*m32;
//    m33 = m01*m12*m20 - m02*m11*m20 + m02*m10*m21 - m00*m12*m21 - m01*m10*m22 + m00*m11*m22;
//    scale(1/m1.determinant()); // multiply mat elements by this!!!
// }
// 
// 
// 
// int determinant() {
//    double value;
//    value =
//    m03*m12*m21*m30 - m02*m13*m21*m30 - m03*m11*m22*m30 + m01*m13*m22*m30+
//    m02*m11*m23*m30 - m01*m12*m23*m30 - m03*m12*m20*m31 + m02*m13*m20*m31+
//    m03*m10*m22*m31 - m00*m13*m22*m31 - m02*m10*m23*m31 + m00*m12*m23*m31+
//    m03*m11*m20*m32 - m01*m13*m20*m32 - m03*m10*m21*m32 + m00*m13*m21*m32+
//    m01*m10*m23*m32 - m00*m11*m23*m32 - m02*m11*m20*m33 + m01*m12*m20*m33+
//    m02*m10*m21*m33 - m00*m12*m21*m33 - m01*m10*m22*m33 + m00*m11*m22*m33;
//    return value;
// }
   
#endif
