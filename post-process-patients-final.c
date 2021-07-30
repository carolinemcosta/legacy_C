#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include "utils.h"
#include "ioutils.h"

int printHelp();
double computeVolumeOfThreshGradients(double *repgrad, double *scaract, double gradthresh, double scarthresh, int nelems, int **elems, double *elemvols);
double computeVolumeOfThreshGradientsNoScar(double *repgrad, double gradthresh, int nelems, int **elems, double *elemvols);
double computeLVvolumes(int nelems, int **elems, double *elemvols, double *bzvol, double *healthyvol);


int main(int argc, char *argv[]){

  char cmd[buffersize], parfile[buffersize], ptsname[buffersize], vname[buffersize],
       actfile[buffersize], repfile[buffersize], apdfile[buffersize], gradrepfile[buffersize], actScarEKfile[buffersize], 
       vtkfile[buffersize], currentdir[buffersize], outfile[buffersize], actEKfile[buffersize], apexvtxname[buffersize], endovtxname[buffersize], epivtxname[buffersize],
       repvtxname[buffersize], apdvtxname[buffersize], repEKsimID[buffersize], apdEKsimID[buffersize], drfilename[buffersize], dafilename[buffersize],
       *carpexec, *imeshname, *simID, *emeshname, *omeshname, *fmeshname, *apexBaseEKsimID, *scarEKsimID, *drdir;
       
  double **points, **ipoints, **fpoints, **cog, *meshcog, *act, *rep, *apd, *gradrep, *scaract, res, *actmax, **bzrep,
         maxrep, minrep, **baserep, transrep, mbaserep, mapexrep,
         maxapd, minapd, **baseapd, transapd, mbaseapd, mapexapd,
	 actmaxrep, actmaxapd, mbaseact, 
	 thickness,
	 totalDR, totalDA, bzDR, bzDA, abDR, abDA, transDR, transDA,
	 *elemvols, *felemvols, scarvol, bzvol, healthyvol;
  
  int **elems, nelems, **felems, fnelems, nelemnodes, etype, nnodes, fnnodes,
      idmaxrep[1], idminrep[1], idmaxapd[1], idminapd[1],
      *ahamap, *apexvtx, **pairs,
      sznvec, scar;
  
  FILE *f1, *f2, *f3;


  if(argc!=7) {
    printHelp();
    exit(1);    
  }
  
  // parameters
  imeshname       = argv[1];
  fmeshname       = argv[2];
  simID           = argv[3];
  drdir           = argv[4];
  scar            = atoi(argv[5]);
  scarEKsimID     = argv[6];  


  if(scar==0) {
    strcpy(imeshname,fmeshname); // use original mesh
  }
  
  sprintf(drfilename,"%s/rep-metrics.dat",drdir);
  f1 = fopen(drfilename,"a");
  
  // files
  sprintf(gradrepfile,"%s/grad-rep.ctr.mag.dat",simID);  
  sprintf(actScarEKfile,"%s/vm_act_seq.dat",scarEKsimID);

  // read original mesh
  felems = readElemFile(fmeshname, &fnelems, &nelemnodes, &etype);
  fpoints = readPointsFile(fmeshname, &fnnodes);
  felemvols = computeTetsVolume(felems, fnelems, fpoints);
      
  if(scar) {
    // read intra mesh
    elems = readElemFile(imeshname, &nelems, &nelemnodes, &etype);
    points = readPointsFile(imeshname, &nnodes);
    
    // repolarisation times and repolarisation gradients
    gradrep = readActivationFile(gradrepfile,nelems,1); // gradient interpolated to element centers

    // scar activation
    scaract = readActivationFile(actScarEKfile,nnodes,1);

    // compute mesh elements volume
    // igrid
    elemvols = computeTetsVolume(elems, nelems, points);

    // compute volume of repolarisation gradients above threshold
    // 5 mm
    sprintf(drfilename,"%s/grad-rep-metrics-3.0-5mm.dat",drdir);
    f2 = fopen(drfilename,"a");

    fprintf(f2,"%f ", computeVolumeOfThreshGradients(gradrep, scaract, 3.0, 6.0, nelems, elems, elemvols));
    fclose(f2);
    
    // 10 mm
    sprintf(drfilename,"%s/grad-rep-metrics-3.0-10mm.dat",drdir);
    f2 = fopen(drfilename,"a");

    fprintf(f2,"%f ", computeVolumeOfThreshGradients(gradrep, scaract, 3.0, 11.0, nelems, elems, elemvols));
    fclose(f2);
    
    // 20 mm
    sprintf(drfilename,"%s/grad-rep-metrics-3.0-20mm.dat",drdir);
    f2 = fopen(drfilename,"a");

    fprintf(f2,"%f ", computeVolumeOfThreshGradients(gradrep, scaract, 3.0, 21.0, nelems, elems, elemvols));
    fclose(f2);    
  
    // all
    sprintf(drfilename,"%s/grad-rep-metrics-3.0-all.dat",drdir);
    f2 = fopen(drfilename,"a");

    fprintf(f2,"%f ", computeVolumeOfThreshGradients(gradrep, scaract, 3.0, 5000.0, nelems, elems, elemvols));
    fclose(f2);
  
    sprintf(drfilename,"%s/lv-volumes.dat",drdir);
    if(~checkOpenFile(drfilename)) {    
      scarvol = 0.;
      bzvol = 0.;
      healthyvol = 0.;
      for(int i=0; i<fnelems; i++) {
        if(felems[nelemnodes][i] == 3)
          scarvol += felemvols[i];
        else if(felems[nelemnodes][i] == 4)
          bzvol += felemvols[i];
        else
          healthyvol += felemvols[i];
      }
      
      f3 = fopen(drfilename,"a");
      fprintf(f3,"%f %f %f\n", scarvol, bzvol, healthyvol);
    }  
    
    // free memory
    for(int i=0; i<nelemnodes+1; i++)
      free(elems[i]);
    free(elems);
    
    for(int i=0; i<3; i++)
      free(points[i]);
    free(points);
    
    free(gradrep);
    free(scaract);
    
    free(elemvols);    
  }
  else { // DR - no scar
    // read grad rep without scar
    gradrep = readActivationFile(gradrepfile,fnelems,1); // gradient interpolated to element centers
    
    // all
    sprintf(drfilename,"%s/grad-rep-metrics-3.0-all.dat",drdir);
    f2 = fopen(drfilename,"a");
    fprintf(f2,"%f ", computeVolumeOfThreshGradientsNoScar(gradrep, 3.0, fnelems, felems, felemvols));
    fclose(f2);
  }    
  
  // free memory
  for(int i=0; i<nelemnodes+1; i++)
    free(felems[i]);
  free(felems);
  
  for(int i=0; i<3; i++)
    free(fpoints[i]);
  free(fpoints);

  free(felemvols);
  
  return 0;
}

int printHelp() {
  char str[buffersize*5];
  
  strcpy(str, "Usage: post-process-patient [imeshname] [fmeshname] [simID] [drdir] [Scar] [scarEKsimID]\n");
  strcat(str, "Always use the full path!\n");
  strcat(str, "imeshname:         Intracellular grid\n");
  strcat(str, "fmeshname:         Labelled smooth mesh\n");
  strcat(str, "simID:             Simulation directory\n");
  strcat(str, "drdir:             Dispersion of repolarisation metrics directory (output) \n");
  strcat(str, "Scar:              Simlation with scar (1) or no scar (2) \n");
  strcat(str, "scarEKsimID:       Scar Eikonal simulation directory \n");
  
  printf("%s",str);

  return 0;  
}

double computeVolumeOfThreshGradients(double *repgrad, double *scaract, double gradthresh, double scarthresh, int nelems, int **elems, double *elemvols) {
  double gradvol;
  
  gradvol  = 0.;
  for(int i=0; i<nelems; i++) {
    if(repgrad[i] >= gradthresh && scaract[elems[0][i]] <= scarthresh && scaract[elems[1][i]] <= scarthresh && scaract[elems[2][i]] <= scarthresh && scaract[elems[3][i]] <= scarthresh)
      gradvol += elemvols[i];
  }
  
  return gradvol;
}

double computeVolumeOfThreshGradientsNoScar(double *repgrad, double gradthresh, int nelems, int **elems, double *elemvols) {
  double gradvol;
  
  gradvol  = 0.;
  for(int i=0; i<nelems; i++) {
    if(repgrad[i] >= gradthresh)
      gradvol += elemvols[i];
  }
  
  return gradvol;
}

double computeLVvolumes(int nelems, int **elems, double *elemvols, double *bzvol, double *healthyvol) {

  *healthyvol = 0.;
  *bzvol = 0.;
  for(int i=0; i<nelems; i++) {
    if(elems[4][i] == 4)
      *bzvol += elemvols[i];
    else
      *healthyvol += elemvols[i];
  }
  
  return *bzvol + *healthyvol;
}
