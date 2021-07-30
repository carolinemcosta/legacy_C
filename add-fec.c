#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include "utils.h"
#include "ioutils.h"

int main(int argc, char *argv[]){

  char cmd[buffersize], fecmeshname[buffersize], sfecmeshname[buffersize], *meshname, *laplaceepi, *ofibres;
       
  double **points, **lon, **nslon, *episol, episolmean;
  
  int **elems, nelems, nelemnodes, etype, nnodes, nlon, scartag;
  
  // parameters
  meshname    = argv[1];
  laplaceepi  = argv[2];
  ofibres     = argv[3];
  scartag     = atoi(argv[4]);
  
  // read mesh
  elems  = readElemFile(meshname, &nelems, &nelemnodes, &etype);
  points = readPointsFile(meshname, &nnodes);
  
  // read epi-endo laplace solution
  episol = readActivationFile(laplaceepi, nnodes, 1);

  // read lon with and without scar
  lon = readLonFile(meshname, &nlon, nelems);
  nslon = readLonFile(ofibres, &nlon, nelems);
  
  for(int i=0; i<nelems; i++) {
    // get mean laplace at element
    episolmean = (episol[elems[nelemnodes-1][i]] + episol[elems[nelemnodes-2][i]] + episol[elems[nelemnodes-3][i]] + episol[elems[nelemnodes-4][i]])/4.0;    

    if(episolmean<=0.1) {
      elems[nelemnodes][i] += 200;
    }
    
    // change fibres over scar to original vectors
    if(elems[nelemnodes][i] == scartag+200) {
      lon[0][i] = nslon[0][i];
      lon[1][i] = nslon[1][i];
      lon[2][i] = nslon[2][i];
      lon[3][i] = nslon[3][i];
      lon[4][i] = nslon[4][i];
      lon[5][i] = nslon[5][i];
    }
  }
      
  // no FEC over scar
  // write new element file
  sprintf(fecmeshname,"%s-fec",meshname);
  writeElemFile(fecmeshname, elems, nelems, nelemnodes, etype);
  
  // copy pts and lon files
  sprintf(cmd,"cp %s.pts %s.pts", meshname, fecmeshname);
  system(cmd);
  sprintf(cmd,"cp %s.lon %s.lon", meshname, fecmeshname);
  system(cmd);

  // FEC over scar
  sprintf(sfecmeshname,"%s-fec-scar",meshname);
  // copy pts and elem files
  sprintf(cmd,"cp %s.pts %s.pts", fecmeshname, sfecmeshname);
  system(cmd);
  sprintf(cmd,"cp %s.elem %s.elem", fecmeshname, sfecmeshname);
  system(cmd);  
  // write new lon file
  writeLonFile(sfecmeshname, lon, nlon, nelems);
  
  // free memory
  for(int i=0; i<nelemnodes+1; i++)
    free(elems[i]);
  free(elems);
  
  for(int i=0; i<3; i++)
    free(points[i]);
  free(points);

  for(int i=0; i<3; i++) {
    free(lon[i]);
    free(nslon[i]);
  }
  free(lon);
  free(nslon);

  free(episol);
  
  return 0;
}
