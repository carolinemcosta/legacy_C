#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include "utils.h"
#include "ioutils.h"

int main(int argc, char *argv[]){

  char outfile[buffersize], *meshname;
       
  double **points, *elemvols, corevol, bzvol, scarvol, lvvol, bzscarratio, bzlvratio, scarlvratio;
  
  int **elems, nelems, nelemnodes, etype, nnodes;
  
  // parameters
  meshname       = argv[1];
  
  sprintf(outfile,"%s.vols",meshname);
  FILE *f1 = fopen(outfile,"w");

  
  elems = readElemFile(meshname, &nelems, &nelemnodes, &etype);
  points = readPointsFile(meshname, &nnodes);
  
  elemvols = computeTetsVolume(elems, nelems, points);

  corevol = 0.;
  bzvol = 0.;
  lvvol = 0.;
  
  for(int i=0; i<nelems; i++) {
    if(elems[nelemnodes][i] == 3 || elems[nelemnodes][i] == 8)
      corevol += elemvols[i];
    else if(elems[nelemnodes][i] == 4 || elems[nelemnodes][i] == 9)
      bzvol += elemvols[i];
    else
      lvvol += elemvols[i];
  }

  scarvol = corevol + bzvol;
  
  bzscarratio = bzvol/scarvol*100.;
  bzlvratio = bzvol/lvvol*100.;
  scarlvratio = scarvol/lvvol*100.;
   
  fprintf(f1,"%f %f %f\n",bzscarratio,bzlvratio,scarlvratio);
}
