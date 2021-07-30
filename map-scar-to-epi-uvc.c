#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include "utils.h"
#include "ioutils.h"

int printHelp();

int main(int argc, char *argv[]){

  char phivtxfile[buffersize], phictrfile[buffersize], rhovtxfile[buffersize], rhoctrfile[buffersize], phictrfixfile[buffersize], zvtxfile[buffersize], zctrfile[buffersize],
       cmd[buffersize], nmeshname[buffersize], cogfile[buffersize], outfile[buffersize], *meshname, *omeshname, *uvcdir;
       
  double *phictr, *phivtx, *rhoctr, *rhovtx, *zctr, *zvtx,
         maxphiinit, minphiinit, maxrhoinit, minrhoinit, maxzinit, minzinit,
	 maxphi, minphi, maxrho, minrho, maxz, minz, 
	 invrho, **cog, **points,
	 *maxphivec, *minphivec, *maxrhovec, *minrhovec, *maxzvec, *minzvec;
  
  int **elems, **oelems, nelems, nelemnodes, etype, nnodes, label;
  
  if(argc!=4) {
    printHelp();
    exit(1);    
  }
  
  // parameters
  omeshname = argv[1];
  meshname  = argv[2];
  uvcdir    = argv[3];

  // read meshes
  elems  = readElemFile(meshname, &nelems, &nelemnodes, &etype);
  oelems = readElemFile(omeshname, &nelems, &nelemnodes, &etype);
  points = readPointsFile(meshname, &nnodes);

  // ----- Interpolate uvc from nodes to element centers -----
  
  // PHI - circumferential
  sprintf(phivtxfile,"%s/COORDS_PHI_NORM.dat",uvcdir);
  sprintf(phictrfile,"%s/COORDS_PHI_NORM_CTR.dat",uvcdir);
  if(access(phictrfile, F_OK) == -1) {
    sprintf(cmd,"meshtool interpolate node2elem -omsh=%s -idat=%s -odat=%s",meshname,phivtxfile,phictrfile);
    system(cmd);
  }
  
  // read vertex coordinates
  phivtx = readActivationFile(phivtxfile,nelems,1);

  sprintf(phictrfixfile,"%s/COORDS_PHI_NORM_CTR_FIXED.dat",uvcdir);
  if(access(phictrfixfile, F_OK) == -1) { 
    // read interpolated coordinates
    phictr = readActivationFile(phictrfile,nelems,1);
    
    // fix values at junction of PHI coordinate
    for(int i=0; i<nelems; i++) {
      if(phivtx[elems[0][i]] >= 0.99 || phivtx[elems[1][i]] >= 0.99 || phivtx[elems[2][i]] >= 0.99 || phivtx[elems[3][i]] >= 0.99 ||
        phivtx[elems[0][i]] <= 0.001 || phivtx[elems[1][i]] <= 0.001 || phivtx[elems[2][i]] <= 0.001 || phivtx[elems[3][i]] <= 0.001)
        phictr[i] = 1.;    
    }
    // write out file for testing
    writeActivationFile(phictrfixfile, phictr, nelems);
    
    // compute and write out cogs for visualization
    sprintf(cogfile,"%s-cog",meshname);
    if(access(cogfile, F_OK) == -1) {
      cog = computeElemCOG(points, elems, nelems, nelemnodes);
      writePointsFile(cogfile, cog, nelems);
      for(int i=0; i<3; i++)
        free(cog[i]);
      free(cog);  
    }
  }
  else // read fixed interpolated coordinates
    phictr = readActivationFile(phictrfixfile,nelems,1);

  
  // RHO - transmural
  sprintf(rhovtxfile,"%s/COORDS_RHO.dat",uvcdir);
  sprintf(rhoctrfile,"%s/COORDS_RHO_CTR.dat",uvcdir);
  if(access(rhoctrfile, F_OK) == -1) {
    sprintf(cmd,"meshtool interpolate node2elem -omsh=%s -idat=%s -odat=%s",meshname,rhovtxfile,rhoctrfile);
    system(cmd);
  }

  // read vertex coordinates
  rhovtx = readActivationFile(rhovtxfile,nelems,1); 
  // read interpolated coordinates
  rhoctr = readActivationFile(rhoctrfile,nelems,1); 
  
  
  // Z - apico-basal
  sprintf(zvtxfile,"%s/COORDS_Z.dat",uvcdir);
  sprintf(zctrfile,"%s/COORDS_Z_CTR.dat",uvcdir);
  if(access(zctrfile, F_OK) == -1) {
    sprintf(cmd,"meshtool interpolate node2elem -omsh=%s -idat=%s -odat=%s",meshname,zvtxfile,zctrfile);
    system(cmd);
  }

  // read vertex coordinates
  zvtx   = readActivationFile(zvtxfile,nelems,1); 
  // read interpolated coordinates
  zctr   = readActivationFile(zctrfile,nelems,1); 
  
  
  // ---- Map scar/BZ from endo to epi ----
  printf("Finding elements within bounding box with inverted RHO...\n");

  for(int i=0; i<nelems; i++) {
    if(elems[nelemnodes][i]>1) {
      
      // element UVC bounding box - min/max of PHI, RHO, and Z of element nodes
      maxphiinit = -100.; minphiinit = 100.;
      maxrhoinit = -100.; minrhoinit = 100.;
      maxzinit   = -100.; minzinit   = 100.;
      
      maxphi = maxphiinit; minphi = minphiinit;
      maxrho = maxrhoinit; minrho = minrhoinit;
      maxz   = maxzinit;   minz   = minzinit;
      
      // PHI - treat junction
      if(phictr[i] == 1.) {
        maxphi = 1.;
        minphi = 0.99;
      }
      else {
        for(int e=0; e<nelemnodes; e++) {
          if(phivtx[elems[e][i]] < minphi)
            minphi = phivtx[elems[e][i]];
          if(phivtx[elems[e][i]] > maxphi)
            maxphi = phivtx[elems[e][i]];
        }
      }
      
      for(int e=0; e<nelemnodes; e++) {
        // inverted RHO
        invrho = 1-rhovtx[elems[e][i]];
        if(invrho < minrho)
          minrho = invrho;
        if(invrho > maxrho)
          maxrho = invrho;

        // Z
        if(zvtx[elems[e][i]] < minz)
          minz = zvtx[elems[e][i]];
        if(zvtx[elems[e][i]] > maxz)
          maxz = zvtx[elems[e][i]];
      }

      // find elements within bounding box with inverted RHO
      for(int j=0; j<nelems; j++) {
        if(rhoctr[j] >= minrho && rhoctr[j] <= maxrho && phictr[j] >= minphi && phictr[j] <= maxphi && zctr[j] >= minz && zctr[j] <= maxz)
          oelems[nelemnodes][j] = elems[nelemnodes][i];
      }
    }
  }
  
  // ---- Write new files ----
  printf("Writing mesh files...\n");
  
  sprintf(nmeshname,"%s-invscar",meshname);
  writeElemFile(nmeshname, oelems, nelems, nelemnodes, etype);
  // copy points and lon files
  sprintf(cmd,"cp %s.pts %s.pts",meshname,nmeshname);
  system(cmd);
  sprintf(cmd,"cp %s.lon %s.lon",meshname,nmeshname);
  system(cmd);
  
  // ---- Free memory ----
  printf("Memory cleanup...\n");

  for(int i=0; i<nelemnodes+1; i++) {
    free(elems[i]);
    free(oelems[i]);
  }
  free(elems);
  free(oelems);
  
  for(int i=0; i<3; i++)
    free(points[i]);
  free(points);

  free(phivtx);
  free(rhovtx);
  free(zvtx);
  
  free(phictr);
  free(rhoctr);
  free(zctr);
  
 
  return 0;
}

int printHelp() {
  char str[buffersize*5];
  
  strcpy(str, "Usage: map-scar-to-epi [meshname] [uvcdir]\n");
  strcat(str, "Always use the full path!\n");
  strcat(str, "omeshname:         Unlabelled mesh \n");
  strcat(str, "meshname:          Labelled mesh \n");
  strcat(str, "uvcdir:            UVC directory\n");
  
  printf("%s",str);

  return 0;  
}

