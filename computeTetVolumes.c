#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include "utils.h"
#include "ioutils.h"

int main(int argc, char *argv[]){

  char *meshname,volsname[buffersize];
       
  double **points, *elemvols;
  
  int **elems, nelems, nelemnodes, etype, nnodes;
  
  FILE *f1;


  // parameters
  meshname  = argv[1];

  // read mesh mesh
  elems = readElemFile(meshname, &nelems, &nelemnodes, &etype);
  points = readPointsFile(meshname, &nnodes);
  
  // compute volume of each tet
  elemvols = computeTetsVolume(elems, nelems, points);
  
  sprintf(volsname,"%s.tvols",meshname);
  f1 = fopen(volsname,"w");
  
  for(int i=0; i<nelems; i++)
    fprintf(f1,"%f\n",elemvols[i]);

  fclose(f1);
}
