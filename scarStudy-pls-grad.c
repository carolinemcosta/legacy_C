#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include "utils.h"
#include "ioutils.h"

int labelEllipticalScar(int **elems, int nelems, int nelemnodes, double **lon, double *meshcog, double **cog, double length, double width, double sizebz, int scarshape);
int labelEllipticalScarGrad(int **elems, int nelems, int nelemnodes, double **lon, double *meshcog, double **cog, double length, double width, double sizebz, int scarshape);
int labelRoundRectangularScar(int **elems, int nelems, int nelemnodes, double **lon, double *meshcog, double **cog, double length, double width, double sizebz, int scarshape);
int createBasicParameterFile(char *parfile);
char* createSimulationID(char *meshname, int ionicprop, int conductivity, int stimlocation, int stimdist);
char* createCARPcommandLine(char *carpexec, int nprocs, char *parfile, char *meshname, char *simID, int ionicprop, int conductivity, int stimlocation, int stimdist, double length, int hpc, double meshsize);
int printHelp();
int logSimError(char *currentdir, char *simID, int ionicprop, int conductivity, int stimlocation, int stimdist, int scarshape, int length);
int createTomPBSscript_(int nprocs, char *parfile, char *meshname, char *simID, int ionicprop, int conductivity, int stimlocation, int stimdist, int scarshape, double length, int hpc, double meshsize, int randomfibres);
int createArcherPBSscript_(int nnodes, char *parfile, char *meshname, char *simID, int ionicprop, int conductivity, int stimlocation, int stimdist, int scarshape, double length, int hpc, double meshsize, int randomfibres);

int main(int argc, char *argv[]){

  char cmd[buffersize*10], 
       meshname[buffersize], imeshname[buffersize], imeshptsname[buffersize], omeshname[buffersize], nmeshname[buffersize], ptsname[buffersize], nptsname[buffersize], 
       parfile[buffersize], actfile[buffersize], nactfile[buffersize], repfile[buffersize], nrepfile[buffersize], apdfile[buffersize], vtkfile[buffersize], gradrepname[buffersize],
       simID[buffersize], currentdir[buffersize], outfile[buffersize];
       
  const char *shape;
  
  char *carpexec, *mesherexec, *gradexec, *vtkexec;
       
  double **points, **lon, **cog, *meshcog, *act, *rep, *apd, length, width, bz, meshsize, resolution;
  
  int **elems, nelems, nelemnodes, etype, nnodes, nLon,
      ionicprop, conductivity, stimlocation, stimdist, scarshape,
      built, actOK, apdOK, hpc, nprocs, randomfibres;
  
  struct stat sb;

  if(argc != 15) {
    printHelp();
    exit(1);
  }
  
  // executables
  carpexec     = argv[1];
  mesherexec   = argv[2];
  gradexec     = argv[3];
  
  // simulation parameters
  ionicprop    = atoi(argv[4]); // ionic changes - 0: none, 1: reduced APD, 2: increased APD
  conductivity = atoi(argv[5]); // 0: normal, 1: equally reduced, 2: increased anisotropy, 3: isotropic, 4: slow isotropic
  stimlocation = atoi(argv[6]); // 0: 90 degreee, 1: 45 degrees, 2: 0 degrees
  stimdist     = atoi(argv[7]); // distance to scar center - 0: 15mm, 1: 13mm, 2: 11mm
  scarshape    = atoi(argv[8]); // 0: scar longitudinal to fibers, 1: transvese to fibers
  length       = atof(argv[9])*1000.0; // ellipse major axis (radius) 5,7,9mm (5mm=circle)
  
  // additional parameters
  hpc          = atoi(argv[10]); // run on hpc system - generate PBS script
  randomfibres = atoi(argv[11]);
  
  // mesh setup
  meshsize     = atof(argv[12]);
  resolution   = atof(argv[13]);
  vtkexec      = argv[14];
  
  // set meshname in current directory
  if(hpc)
    currentdir[0] = 0;
  else {
    getcwd(currentdir, sizeof(currentdir));
    strcat(currentdir,"/");
  }
  
  if(length == 0.0) {
    sprintf(omeshname, "%dcm-no-scar",(int)meshsize); // base name
    sprintf(meshname, "%s%s",currentdir,omeshname); // add directory
    strcpy(nmeshname,meshname);
  }
  else {
    shape = scarshape?"trans":"lon";
    if(randomfibres)
      sprintf(omeshname, "%dcm-scar-%s-%dmm-10mm-randomfibres",(int)meshsize,shape,atoi(argv[9])*2);
    else
      sprintf(omeshname, "%dcm-scar-%s-%dmm-10mm",(int)meshsize,shape,atoi(argv[9])*2);
    if(conductivity==8) // bz gradient
      strcat(omeshname, "-BZgrad4"); 
    sprintf(meshname, "%s%s",currentdir,omeshname); 

    
    // labeled mesh
    strcpy(nmeshname,meshname);
    strcat(nmeshname,"-labeled");
    strcat(omeshname,"-labeled");
  }

  // create mesh - if it doesn't exist
  strcpy(ptsname,meshname);
  strcat(ptsname,".pts");
  
  if(!checkOpenFile(ptsname)) {
    // call mesher
    sprintf(cmd,"mesher -mesh %s -size[0] %f -size[1] %f -size[2] 0.0 -resolution[0] %f -resolution[1] %f -resolution[2] %f -tri2D", meshname, meshsize, meshsize, resolution, resolution, resolution);
    printf("%s\n",cmd);
    system(cmd);
  }
    
  // label mesh - if not labeled
  strcpy(nptsname,nmeshname);
  strcat(nptsname,".pts");
  
  if(!checkOpenFile(nptsname) && length > 0.0) {    
    // read mesh
    nnodes = getNnodes(ptsname);
    points = readPointsFile(meshname, &nnodes);
    elems  = readElemFile(meshname, &nelems, &nelemnodes, &etype);
    lon    = readLonFile(meshname, &nLon, nelems);

    // get COGs  
    cog     = computeElemCOG(points, elems, nelems, nelemnodes);
    meshcog = computeMeshCOG(cog, nelems);
    
    // scar dimensions
    width = 5000.0; // core width 10 mm 
    bz = 2000.0; // BZ width 2 mm (total width 14mm)

    // label scar 
//     labelEllipticalScar(elems, nelems, nelemnodes, lon, meshcog, cog, length, width, bz, scarshape);
    labelEllipticalScarGrad(elems, nelems, nelemnodes, lon, meshcog, cog, length, width, bz, scarshape);
//     labelRoundRectangularScar(elems, nelems, nelemnodes, lon, meshcog, cog, length, width, bz, scarshape);

    // generate random fibres at the BZ - alternative isotropy representation
    if(randomfibres)
      generateRandomFibresTag2D(nelems, nelemnodes, elems, 4, lon);

    // write labeled mesh
    writeElemFile(nmeshname, elems, nelems, nelemnodes, etype);
    writeLonFile(nmeshname, lon, nLon, nelems);
    // copy points file    
    sprintf(cmd,"cp %s %s",ptsname,nptsname);
    system(cmd);    
  }
  
  // create basic parameter file
  sprintf(parfile,"%sscarEP-5pls-itensor-dt10-spacedt20-pt.par",currentdir);
  if(access(parfile, F_OK) == -1) {
    createBasicParameterFile(parfile);
  }

  // create simulation directory and check if it already exists
  strcpy(simID, createSimulationID(meshname, ionicprop, conductivity, stimlocation, stimdist));  

  // run simulation locally or create HPC script
  switch(hpc) {
    case 0:
      nprocs = 4;
      
      // original act and rep files
      sprintf(actfile,"%s/vm_activation.dat",simID);
      sprintf(repfile,"%s/vm_repolarisation.dat",simID);
      
      // post-processed act and rep files
      sprintf(nactfile,"%s/vm_activation_nid.dat",simID);
      sprintf(nrepfile,"%s/vm_repolarisation_nid.dat",simID);      
      
      // build labeled model
      sprintf(imeshname,"%s/grid_i/%s_i",currentdir,omeshname);
      printf("%s\n",imeshname);
      sprintf(imeshptsname,"%s.pts",imeshname);      
      buildIntraGrid(imeshptsname, cmd, carpexec, currentdir, nmeshname);
      
      // create command line
      strcpy(cmd,createCARPcommandLine(carpexec, nprocs, parfile, nmeshname, simID, ionicprop, conductivity, stimlocation, stimdist, length, hpc, meshsize));
      printf("%s\n",cmd);
      
      // try to access vm_repolarisation
      if(!checkOpenFile(repfile)) { 
        system(cmd); // call CARP
        // check if simulation was successful
        if(!checkOpenFile(repfile)) { 
          printf("Simulation failed: %s\n", simID);
          logSimError(currentdir, simID, ionicprop, conductivity, stimlocation, stimdist, scarshape, length);
          exit(1);
        }
      }
      else {
        printf("Simulation did not run: %s\n", simID);
      }


      // check if apd and gradients have already been computed
      sprintf(apdfile,"%s/apd.dat",simID);
      
      if (!checkOpenFile(apdfile)) {
        // update number of nodes - now for intra grid 
        nnodes = getNnodes(imeshptsname);
        act = readActivationFilePls(actfile, nnodes, 0, 5);
        rep = readActivationFilePls(repfile, nnodes, 0, 5);

        // compute APD
        apd = computeAPD(act,rep,nnodes);

        // new activation and vm_repolarisation files
        writeActivationFile(nactfile, act, nnodes);
        writeActivationFile(nrepfile, rep, nnodes);
        writeActivationFile(apdfile, apd, nnodes);

        // compute gradients using GlGradient
        sprintf(cmd, "%s -m %s -d %s -c 8 -o %s/grad-apd.dat -s 1000.0", gradexec, imeshname, apdfile, simID);
        system(cmd);
        sprintf(cmd, "%s -m %s -d %s -c 8 -o %s/grad-act.dat -s 1000.0", gradexec, imeshname, nactfile, simID);
        system(cmd);
        sprintf(cmd, "%s -m %s -d %s -c 8 -o %s/grad-rep.dat -s 1000.0", gradexec, imeshname, nrepfile, simID);
        system(cmd);	
      }
      
      // convert to VTK
      sprintf(vtkfile,"%s.vtk",simID);
      if (!checkOpenFile(vtkfile)) {
        sprintf(vtkfile,"%s",simID); // just remove file extension
        sprintf(gradrepname,"%s/grad-rep.dat.vtx.mag.dat",simID);
        sprintf(cmd, "%s -m %s -n %s -n %s -n %s -o %s", vtkexec, imeshname, nactfile, nrepfile, gradrepname, vtkfile);
        system(cmd);     
      }
      
      break;
      
    // tom
    case 1:
      nprocs = 64; //4;
      createTomPBSscript_(nprocs, parfile, nmeshname, simID, ionicprop, conductivity, stimlocation, stimdist, scarshape, length, hpc, meshsize,randomfibres);
      break;
      
    //archer
    case 2:
      nprocs = 24*24; // 24 = # tasks per compute node
      createArcherPBSscript_(nprocs, parfile, nmeshname, simID, ionicprop, conductivity, stimlocation, stimdist, scarshape, length, hpc, meshsize,randomfibres);
      break;
    default:
      printf("unknown platform %d\n", hpc);
      exit(1);      
  }
  
  // free memory
  for(int i=0; i<nelemnodes+1; i++)
    free(elems[i]);
  free(elems);
  
  for(int i=0; i<3; i++)
    free(points[i]);
  free(points);
  
  for(int i=0; i<3; i++)
    free(lon[i]);
  free(lon);
  
  for(int i=0; i<3; i++)
    free(cog[i]);
  free(cog);
  
  free(act);
  free(rep);
  free(apd);
  
  return 0;
}

int labelEllipticalScar(int **elems, int nelems, int nelemnodes, double **lon, double *meshcog, double **cog, double length, double width, double sizebz, int scarshape) {
  int count, inbz, incore;
  double size0, size1;

  if(scarshape) {
    size0 = width;
    size1 = length;
  }
  else { // reverse sizes
    size0 = length;
    size1 = width;
  }
  
  for(int i=0; i<nelems; i++) {
    // check if element is inside core and or bz
    inbz   = inEllipse(cog[i][0], meshcog[0], size0+sizebz, cog[i][1], meshcog[1], size1+sizebz);
    incore = inEllipse(cog[i][0], meshcog[0], size0, cog[i][1], meshcog[0], size1);
    count  = inbz + incore;
    
    if(count==1) // only in BZ
      elems[nelemnodes][i] = 4;
    else if(count==2) { // core
      elems[nelemnodes][i] = 3;
      lon[0][i] = lon[1][i] = lon[2][i] = 0.0; // bath
    }
  }
}

int labelEllipticalScarGrad(int **elems, int nelems, int nelemnodes, double **lon, double *meshcog, double **cog, double length, double width, double sizebz, int scarshape) {
  int count, inbz, intbz, incore;
  double size0, size1, sizesubbz, sizenbz, nsizebz, nsubbz;

  if(scarshape) {
    size0 = width;
    size1 = length;
  }
  else { // reverse sizes
    size0 = length;
    size1 = width;
  }
  
  sizenbz = 1000.0;
  nsubbz = 20.0;
  sizesubbz = sizenbz/nsubbz;
    
  for(int bzcounter=nsubbz; bzcounter>0; bzcounter--) {
    nsizebz = sizebz + sizesubbz*bzcounter; // size of BZ considering gradient
    printf("%f %f %d\n",nsizebz,sizesubbz,bzcounter);

    for(int i=0; i<nelems; i++) {
      // check if element is inside core and or bz
      intbz  = inEllipse(cog[i][0], meshcog[0], size0+nsizebz, cog[i][1], meshcog[1], size1+nsizebz); // transition
      inbz   = inEllipse(cog[i][0], meshcog[0], size0+sizebz,  cog[i][1], meshcog[1], size1+sizebz); // bz
      incore = inEllipse(cog[i][0], meshcog[0], size0,         cog[i][1], meshcog[0], size1); // core
//       count  = inbz + incore;
      
      if(incore) {
        elems[nelemnodes][i] = 3;
        lon[0][i] = lon[1][i] = lon[2][i] = 0.0; // bath
      }
      else if(inbz && ~intbz) {
        elems[nelemnodes][i] = 4;
      }
      else if (intbz){
        elems[nelemnodes][i] = 4+bzcounter;
      }
	
//       if(count==1) { // only in BZ
// 	elems[nelemnodes][i] = 3+bzcounter;
//       }
//       else if(count==2) { // core
// 	elems[nelemnodes][i] = 3;
// 	lon[0][i] = lon[1][i] = lon[2][i] = 0.0; // bath
//       }
    }
  }
}

int labelRoundRectangularScar(int **elems, int nelems, int nelemnodes, double **lon, double *meshcog, double **cog, double length, double width, double sizebz, int scarshape) {
  double corigin1[2], corigin2[2];
  int count, inbz, incore;
  double size0, size1;

  
  for(int i=0; i<nelems; i++) {
    if(scarshape) {
      size0 = width;
      size1 = length;
      corigin1[0] = meshcog[0]; corigin1[1] = meshcog[1]+length; // origin of top circle
      corigin2[0] = meshcog[0]; corigin2[1] = meshcog[1]-length; // origin of bottom circle
      inbz = (inRectangle(cog[i][0], meshcog[0], size0+sizebz, cog[i][1], meshcog[1], size1) || 
              inEllipse(cog[i][0], corigin1[0], width+sizebz, cog[i][1], corigin1[1], width+sizebz) ||
              inEllipse(cog[i][0], corigin2[0], width+sizebz, cog[i][1], corigin2[1], width+sizebz));      
    }
    else { // reverse sizes
      size0 = length;
      size1 = width;
      corigin1[0] = meshcog[0]+length; corigin1[1] = meshcog[1]; // origin of right circle
      corigin2[0] = meshcog[0]-length; corigin2[1] = meshcog[1]; // origin of left circle
      inbz = (inRectangle(cog[i][0], meshcog[0], size0, cog[i][1], meshcog[1], size1+sizebz) || 
              inEllipse(cog[i][0], corigin1[0], width+sizebz, cog[i][1], corigin1[1], width+sizebz) ||
              inEllipse(cog[i][0], corigin2[0], width+sizebz, cog[i][1], corigin2[1], width+sizebz));      
    }

    // check if element is inside core and or bz
    incore = (inRectangle(cog[i][0], meshcog[0], size0, cog[i][1], meshcog[1], size1) || 
              inEllipse(cog[i][0], corigin1[0], width, cog[i][1], corigin1[1], width) ||
              inEllipse(cog[i][0], corigin2[0], width, cog[i][1], corigin2[1], width));
    
    count = inbz + incore;
    
    if(count==1) // only in BZ
      elems[nelemnodes][i] = 4;
    else if(count==2) { // core
      elems[nelemnodes][i] = 3;
      lon[0][i] = lon[1][i] = lon[2][i] = 0.0; // bath
    }
  }
}

int createBasicParameterFile(char *parfile) {
  
  char parstr[buffersize*10];
  
  FILE *f = fopen(parfile,"w");
  if(!f) {
    printf("Could not open file %s\n", parfile);
    exit(1);
  }
  
  // ionicprop region
  strcpy(parstr,"num_imp_regions          = 1\n");
  strcat(parstr,"imp_region[0].num_IDs    = 1\n");
  strcat(parstr,"imp_region[0].ID[0]      = 1\n");
  strcat(parstr,"imp_region[0].im         = TT2\n");
  strcat(parstr,"imp_region[0].im_sv_init = \"TT2_control_100npls_500bcl_INIT.st\"\n\n");
  
  // conductivity region
  strcat(parstr,"num_gregions       = 1\n");
  strcat(parstr,"gregion[0].num_IDs = 1\n");
  strcat(parstr,"gregion[0].ID[0]   = 1\n");
  strcat(parstr,"gregion[0].g_il    = 0.18959\n");
  strcat(parstr,"gregion[0].g_it    = 0.06915\n");
  strcat(parstr,"gregion[0].g_in    = 0.06915\n");
  
  // use intracellular tensor
  strcat(parstr,"bidm_eqv_mono      = 0\n\n");  
  
  // stimulus
  strcat(parstr,"num_stim             = 1\n");
  strcat(parstr,"stimulus[0].stimtype = 0\n");
  strcat(parstr,"stimulus[0].strength = 250.0\n");
  strcat(parstr,"stimulus[0].duration = 2.0\n");
  strcat(parstr,"stimulus[0].start    = 0.0\n");
  strcat(parstr,"stimulus[0].npls     = 5\n\n");
  strcat(parstr,"stimulus[0].bcl      = 500.0\n\n");
  
  // lats
  strcat(parstr,"num_LATs          =  1\n");
  strcat(parstr,"lats[0].measurand =  0\n");
  strcat(parstr,"lats[0].all       =  1\n");
  strcat(parstr,"lats[0].threshold = -20\n");
  strcat(parstr,"lats[0].method    =  1\n\n");
  
  // apd and rep
  strcat(parstr,"compute_APD       =  1\n");
  strcat(parstr,"actthresh         = -20\n");
  strcat(parstr,"recovery_thresh   = -70\n\n");

  // IO
  strcat(parstr,"spacedt   = 20.0\n");
  strcat(parstr,"timedt    = 100.0\n");
  strcat(parstr,"gridout_i = 0\n");
  
  // time step
  strcat(parstr,"dt      = 10\n"); // dt has to be 10us to compute the gradient
  
  // write out parameter file
  fprintf(f,"%s", parstr);
  fclose(f);
}


char* createSimulationID(char *meshname, int ionicprop, int conductivity, int stimlocation, int stimdist) {
  static char simIDstr[buffersize];
  
  strcpy(simIDstr,meshname);

  // ionic properties parameter  
  switch(ionicprop) {
    case 0:
      strcat(simIDstr,"-normalIonProp");
      break;
    case 1:
      strcat(simIDstr,"-reducedAPD");
      break;
    case 2:
      strcat(simIDstr,"-increasedAPD");
      break;
    default:
      printf("unknown ionic properties parameter %d. Option are 0-2\n", ionicprop);
      exit(1);
  }
  
  // conductivity parameter  
  switch(conductivity) {
    case 0:
      strcat(simIDstr,"-normalConductivity");
      break;
    case 1:
      strcat(simIDstr,"-equallyReduced");
      break;
    case 2:
      strcat(simIDstr,"-increasedAnisotropy");
      break;
    case 3:
      strcat(simIDstr,"-isotropicReduced");
      break;      
    case 4:
      strcat(simIDstr,"-isotropic");
      break;
    case 5:
      strcat(simIDstr,"-isotropic-transCond");
      break;      
    case 6:
      strcat(simIDstr,"-isotropic-reducedTransCond");
      break;       
    case 7:
      strcat(simIDstr,"-isotropic-computedCond");
      break;
    case 8:
      strcat(simIDstr,"-increasedAnisotropyGrad4");
      break;      
    default:
      printf("unknown conductivity parameter %d. Option are 0-4\n", conductivity);
      exit(1);
  }
  
  // stimlocation parameter  
  switch(stimlocation) {
    case 0:
      strcat(simIDstr,"-st90");
      break;
    case 1:
      strcat(simIDstr,"-st0");
      break;
    case 2:
      strcat(simIDstr,"-st45");
      break;
    default:
      printf("unknown stimlocation parameter %d. Option are 0-2\n", stimlocation);
      exit(1);
  }

  // stimlocation parameter  
  switch(stimdist) {
    case 0:
      strcat(simIDstr,"-0mm");
      break;
    case 1:
      strcat(simIDstr,"-35mm");
      break;
    case 2:
      strcat(simIDstr,"-40mm");
      break;
    case 3:
      strcat(simIDstr,"-12mm");
      break;     
    default:
      printf("unknown stimlocation parameter %d. Option are 0-2\n", stimdist);
      exit(1);
  }

  return simIDstr; 
}

char* createCARPcommandLine(char *carpexec, int nprocs, char *parfile, char *meshname, char *simID, int ionicprop, int conductivity, int stimlocation, int stimdist, double length, int hpc, double meshsize) {
  
  char static cmd[buffersize*10],mpi[buffersize];
  double x0,y0,z0,dist,midpoint,gincrement;
  
  switch(hpc) {
    case 0:
      sprintf(mpi,"mpirun -np %d",nprocs);
      break;
    case 1: // tom 
      sprintf(mpi,"mpiexec_mpt -np %d",nprocs);
      break;
    case 2: // archer
      sprintf(mpi,"aprun -n %d -N 24",nprocs);
      break;
    default:
      printf("unknown platform %d\n", hpc);
      exit(1);      
  }
  sprintf(cmd,"%s %s +F %s -meshname %s -simID %s ", mpi, carpexec, parfile, meshname, simID);

  // if scar
  if(length>0.0) { 
    
    if(ionicprop>0) {
      // add modified ionic properties  - don't add intermediate region for now
      strcat(cmd,"-num_imp_regions          2 ");
      strcat(cmd,"-imp_region[1].im         TT2 ");      

      strcat(cmd,"-imp_region[1].num_IDs    22 ");
      if(conductivity==8) {
        for(int counter=0; counter<21; counter++) {
          sprintf(cmd,"%s -imp_region[1].ID[%d]      %d ", cmd, counter, counter+4); 
        }
      }
      else {
        strcat(cmd,"-imp_region[1].num_IDs    1 ");
        strcat(cmd,"-imp_region[1].ID[0]      4 ");
      }

      if(ionicprop==1) {
        // increased GKs 200% - decreased APD
        strcat(cmd,"-imp_region[1].im_sv_init \"TT2_gKs200_dt100_100npls_500bcl_INIT.st\" ");
        strcat(cmd,"-imp_region[1].im_param   \"Gks*2.0\" ");    
      }
      else if(ionicprop==2) {
        // reduced GKs 50% - increased APD
        strcat(cmd,"-imp_region[1].im_sv_init \"TT2_gKs50_dt100_100npls_500bcl_INIT.st\" ");
        strcat(cmd,"-imp_region[1].im_param   \"Gks*0.5\" ");    
      }
    }
    else {
      // assign normal properties to imp_region
      if(conductivity==8) {
        strcat(cmd,"-imp_region[0].num_IDs    22 ");
        for(int counter=0; counter<21; counter++) {
          sprintf(cmd,"%s -imp_region[0].ID[%d]      %d ", cmd, counter, counter+4); 
        }
      }
      else {
        strcat(cmd,"-imp_region[0].num_IDs    2 ");
        strcat(cmd,"-imp_region[0].ID[0]      1 ");
        strcat(cmd,"-imp_region[0].ID[1]      4 ");
      }
    }

    // add modified conductivities
    if (conductivity>0) {      
      strcat(cmd,"-num_gregions       2 ");
      strcat(cmd,"-gregion[1].num_IDs 1 ");
      strcat(cmd,"-gregion[1].ID[0]   4 ");      
  
      switch(conductivity) {
        case 1:
          // equally reduced - 3x
          strcat(cmd,"-gregion[1].g_il    0.06319667 ");
          strcat(cmd,"-gregion[1].g_it    0.02305000 ");
          strcat(cmd,"-gregion[1].g_in    0.02305000 ");
          break;	
        case 2:
          // increased anisotropy - 10% Gt
          strcat(cmd,"-gregion[1].g_il    0.18959000 ");
          strcat(cmd,"-gregion[1].g_it    0.00691500 ");
          strcat(cmd,"-gregion[1].g_in    0.00691500 ");
          break;
        case 3:
          // isotropic reduced 3x
          strcat(cmd,"-gregion[1].g_il    0.06319667 ");
          strcat(cmd,"-gregion[1].g_it    0.06319667 ");
          strcat(cmd,"-gregion[1].g_in    0.06319667 ");
          break;
        case 4:
          // isotropic
          strcat(cmd,"-gregion[1].g_il    0.18959000 ");
          strcat(cmd,"-gregion[1].g_it    0.18959000 ");
          strcat(cmd,"-gregion[1].g_in    0.18959000 ");
          break;
        case 5:
          // isotropic using normal transverse conductivity
          strcat(cmd,"-gregion[1].g_il    0.06915 ");
          strcat(cmd,"-gregion[1].g_it    0.06915 ");
          strcat(cmd,"-gregion[1].g_in    0.06915 ");
          break;
        case 6:
          // isotropic using reduced transverse conductivity
          strcat(cmd,"-gregion[1].g_il    0.006915 ");
          strcat(cmd,"-gregion[1].g_it    0.006915 ");
          strcat(cmd,"-gregion[1].g_in    0.006915 ");
          break;
        case 7:
          // isotropic using computed conductivity
          strcat(cmd,"-gregion[1].g_il    0.06897 ");
          strcat(cmd,"-gregion[1].g_it    0.06897 ");
          strcat(cmd,"-gregion[1].g_in    0.06897 ");
          break;
        case 8:        
          // increased anisotropy - 10% Gt
          strcat(cmd,"-gregion[1].g_il    0.18959000 ");
          strcat(cmd,"-gregion[1].g_it    0.00691500 ");
          strcat(cmd,"-gregion[1].g_in    0.00691500 ");

          // transition
          strcat(cmd,"-num_gregions       22 ");
                gincrement = (0.06915 - 0.006915)/21.0;

          for(int counter=2; counter<22; counter++) {
            // increased anisotropy - 100-10% Gt grad
                  sprintf(cmd,"%s -gregion[%d].num_IDs 1 ",       cmd, counter);
                  sprintf(cmd,"%s -gregion[%d].ID[0]   %d ",      cmd, counter, counter+3);      
            sprintf(cmd,"%s -gregion[%d].g_il    0.18959 ", cmd, counter);
            sprintf(cmd,"%s -gregion[%d].g_it    %1.6f ",   cmd, counter, 0.006915+gincrement*(counter-1));
            sprintf(cmd,"%s -gregion[%d].g_in    %1.6f ",   cmd, counter, 0.006915+gincrement*(counter-1));
          }
          break;	  
        default:
          printf("unknown stimlocation parameter %d. Option are 0-4\n", stimlocation);
          exit(1);	
      }    
    }
    else {
      // assign normal conductivities to gregion
      strcat(cmd,"-gregion[0].num_IDs 2 ");
      strcat(cmd,"-gregion[0].ID[0]   1 ");
      strcat(cmd,"-gregion[0].ID[1]   4 ");
    }
  }
  
  midpoint = (meshsize/2.0)*10000.0;
  // add stimulus location
  switch(stimdist) {
    case 0:
      dist = midpoint;// 15000.0;
      break;
    case 1:
      dist = midpoint-35000.0; //2000.0; //13000.0
      break;
    case 2:
      dist = midpoint-40000.0; //4000.0; //11000.0
      break;
    case 3:
      dist = midpoint-12000.0; //4000.0; //11000.0
      break;      
    default:
      printf("unknown stimlocation parameter %d. Option are 0-2\n", stimdist);
      exit(1);
  }  
  switch(stimlocation) {
    case 0:
      x0 = 0.0;
      y0 = -dist;
      break;
    case 1:
      x0 = -dist;
      y0 = 0.0;
      break;
    case 2:
      x0 = -dist;
      y0 = -dist;
      break;
    default:
      printf("unknown stimlocation parameter %d. Option are 0-2\n", stimlocation);
      exit(1);
  }
  
  strcat(cmd,"-num_stim         1 ");
  
  sprintf(cmd,"%s-stimulus[0].x0 %f ",cmd,x0);
  sprintf(cmd,"%s-stimulus[0].y0 %f ",cmd,y0);
  
  strcat(cmd,"-stimulus[0].z0 0.0 ");
  strcat(cmd,"-stimulus[0].xd 250.0 ");
  strcat(cmd,"-stimulus[0].yd 250.0 ");
  strcat(cmd,"-stimulus[0].zd 100.0 ");
  
  strcat(cmd,"-tend 2500.0 ");
/*  
  // add checkpoint
  strcat(cmd,"-chkpt_start 2300.0 ");
  strcat(cmd,"-chkpt_intv  20.0 ");
  strcat(cmd,"-chkpt_stop  2500.0 ");*/

  
  return cmd;
}

int printHelp() {
  char str[buffersize*5];
  
  strcpy(str, "Usage: scarStudy [carpexec] [mesherexec] [gradexec] [ionicprop] [conductivity] [stimlocation] [stimdist] [scarshape] [scarlength]\n");
  strcat(str, "carpexec:     Path to CARP executable\n");
  strcat(str, "mesherexec:   Path to MESHER executable\n");
  strcat(str, "gradexec:     Path to GlGradient executable\n");  
  strcat(str, "ionicprop:    Changes in ionic properties at the BZ - 0: none, 1: reduced APD, 2: increased APD\n");
  strcat(str, "conductivity: Changes in conductivity at the BZ - 0: normal, 1: equally reduced, 2: increased anisotropy, 3: isotropic, 4: isotropic 3x reduced\n");
  strcat(str, "stimlocation: Stimulus location relative to scar major axis - 0: 90 degreee, 1: 0 degrees, 2: 45 degrees\n");
  strcat(str, "stimdist:     Stimulus distance to scar center - 0: 15mm, 1: 13mm, 2: 11mm\n");
  strcat(str, "scarshape:    Scar shape relative to fibre direction - 0: longitudinal, 1: transvese\n");
  strcat(str, "scarlength:   Length of scar major axis (ellipse radius) - 0,5-9mm (0=no scar, 5mm=circle)\n");
  
  printf("%s",str);

  return 0;  
}

int logSimError(char *currentdir, char *simID, int ionicprop, int conductivity, int stimlocation, int stimdist, int scarshape, int length) {
  char outfile[buffersize];
  sprintf(outfile,"%s/errorLog.dat", currentdir);
  
  FILE *f = fopen(outfile,"a");
  if(!f) {
    printf("Could not open file %s\n", outfile);
    exit(1);
  }
  
  fprintf(f,"%s %d %d %d %d %d %d\n", simID, ionicprop, conductivity, stimlocation, stimdist, scarshape, length);  
  fclose(f);
}

int createTomPBSscript_(int nprocs, char *parfile, char *meshname, char *simID, int ionicprop, int conductivity, int stimlocation, int stimdist, int scarshape, double length, int hpc, double meshsize, int randomfibres) {
  char command[buffersize], jobname[buffersize], scriptname[buffersize], carphpc[buffersize], hpcmeshname[buffersize];
  
  sprintf(jobname,"I%dC%dS%dL%d",ionicprop,conductivity,stimlocation,(int)length/1000);
  sprintf(scriptname,"tom-%d-%s-SS%d-SD%d-RF%d.pbs",(int)meshsize,jobname,scarshape,stimdist,randomfibres);
  
  strcpy(carphpc,"/home/cmc16/bin/carp.petsc.pt");
  sprintf(hpcmeshname,"/home/cmc16/meshes/%s",meshname); // add hpc directory

  // create CARP command line  
  strcat(command,createCARPcommandLine(carphpc, nprocs, parfile, hpcmeshname, simID, ionicprop, conductivity, stimlocation, stimdist, length, hpc, meshsize));
  
  // create script
  createTomPBSscript(nprocs, jobname, scriptname, command);
}

int createArcherPBSscript_(int nprocs, char *parfile, char *meshname, char *simID, int ionicprop, int conductivity, int stimlocation, int stimdist, int scarshape, double length, int hpc, double meshsize, int randomfibres) {
  char command[buffersize], jobname[buffersize], scriptname[buffersize], carphpc[buffersize], hpcmeshname[buffersize];
  
  sprintf(jobname,"I%dC%dS%dL%d",ionicprop,conductivity,stimlocation,(int)length/1000);
  sprintf(scriptname,"archer-%d-%s-SS%d-SD%d-RF%d.pbs",(int)meshsize,jobname,scarshape,stimdist,randomfibres);
  
  sprintf(hpcmeshname,"/work/e348/e348/cmcosta/meshes/%s",meshname); // add hpc directory
  strcpy(carphpc,"/work/e348/e348/shared/carp.petsc.pt.archer.r2957");
  
  // create CARP command line
  strcat(command,createCARPcommandLine(carphpc, nprocs, parfile, hpcmeshname, simID, ionicprop, conductivity, stimlocation, stimdist, length, hpc, meshsize));
  
  // create script
  createArcherPBSscript(nprocs, jobname, scriptname, command);
}
