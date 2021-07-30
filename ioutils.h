#ifndef _IOUTILS
#define _IOUTILS

#define buffersize 2048
#define Tr 0
#define Qd 1
#define Tt 2
#define Hx 3

struct NRRDheader_ {
  int dimension;
  int *sizes;
  double *spacedirections;
  double *origin;
  double *resolutions;
} header;

FILE* openFile(char* filename, int read);
int* readMap(char* meshname, int *nelems);
double** readPointsFile(char* meshname, int *nnodes);
double** readVpointsFile(char* meshname, int *npoints);
int** readElemFile(char* meshname, int *nelems, int *nElemNodes, int *etype);
double** readLonFile(char* meshname, int *nLon, int nelems);
int* readVtxFile(char* meshname, int *nnodes);
int writePointsFile(char* meshname, double **points, int nnodes);
int writeElemFile(char* meshname, int **elems, int nelems, int nElemNodes, int etype);
int writeLonFile(char* meshname, double **lon, int nLon, int nelems);
int writeVtxFile(char* meshname, int nnodes, int *nodes);
int writeMap(char* meshname, int nelems, int *map);
double* readActivationFile(char* actfile, int nnodes, int mode);
double* readActivationFilePls(char* actfile, int nnodes, int mode, int npls);
int readWriteActivationFileAllPls(char* actfile, int nnodes, int npls);
int writeActivationFile(char *actfile, double *act, int nnodes);
unsigned short* readNRRD(char *imfile, NRRDheader_ *header);
int skipNRRDcomment(FILE *fid, char *curstr);
int skipNRRDgarbage(FILE *fid, char *curstr);
double** getNRRDtransMat(NRRDheader_ header);
int getNRRDres(NRRDheader_ *header);
int getNnodes(char* mesh_file);
int getNelems(char* mesh_file);
int checkOpenFile(char * filename);
int createArcherPBSscript(int nprocs, char *jobname, char* scriptname, char* command);
int createTomPBSscript(int nprocs, char *jobname, char* scriptname, char* command);

FILE* openFile(char* filename, int read) {
  FILE *f;

  if(read)
    f = fopen(filename,"r");
  else
    f = fopen(filename,"w");
  
  if(!f) {
    printf("Could not open file %s\n", filename);
  }

  return f;
}

int* readMap(char* mapname, int *nelems) {
  int *map;
  
  FILE *f = fopen(mapname,"r");
  if(!f) {
    printf("readMap: Could not open file %s\n", mapname);
    exit(1);
  }
  fscanf(f, "%d\n", nelems);
  
  map = (int*) malloc(*nelems*sizeof(int));
  
  for(int i=0; i<*nelems; i++) {
    fscanf(f, "%d\n", &map[i]);
  }
  fclose(f);
  
  return map;
}

double** readPointsFile(char* meshname, int *nnodes) {
  double **nodes;
  char pts_file[buffersize];
  
  printf("Reading points file...\n");
  
  strcpy(pts_file,meshname);
  strcat(pts_file,".pts");
  
  FILE *f = fopen(pts_file,"r");
  if(!f) {
    printf("readPointsFile: Could not open file %s\n", pts_file);
    exit(1);
  }
  fscanf(f, "%d\n", nnodes);

  nodes = (double**) malloc(3*sizeof(double*));
  for(int i=0;  i<3; i++)
    nodes[i] = (double*) malloc(*nnodes*sizeof(double));

  for(int i=0; i<*nnodes; i++) {
    fscanf(f, "%lf %lf %lf\n", &nodes[0][i], &nodes[1][i], &nodes[2][i]);
  }
  fclose(f);

  return nodes;
}

double** readVpointsFile(char* meshname, int *npoints) {
  double **nodes;
  char pts_file[buffersize];
  
  printf("Reading points file...\n");
  
  strcpy(pts_file,meshname);
  strcat(pts_file,".vpts");
  
  FILE *f = fopen(pts_file,"r");
  if(!f) {
    printf("readPointsFile: Could not open file %s\n", pts_file);
    exit(1);
  }
  fscanf(f, "%d\n", npoints);

  nodes = (double**) malloc(3*sizeof(double*));
  for(int i=0;  i<3; i++)
    nodes[i] = (double*) malloc(*npoints*sizeof(double));

  for(int i=0; i<*npoints; i++) {
    fscanf(f, "%lf %lf %lf\n", &nodes[0][i], &nodes[1][i], &nodes[2][i]);
  }
  fclose(f);

  return nodes;
}

int** readElemFile(char* meshname, int *nelems, int *nElemNodes, int *etype) {
  int **elems;
  char etype_str[3], d[3];
  char elem_file[buffersize];
  
  printf("Reading elements file...\n");
  
  strcpy(elem_file,meshname);
  strcat(elem_file,".elem");
  
  FILE *f = fopen(elem_file,"r");
  if(!f) {
    printf("readElemFile: Could not open file %s\n", elem_file);
    exit(1);
  }
  fscanf(f, "%d\n", nelems);
  fgets(etype_str, 3, f);
  
  if(!strcmp(etype_str,"Tr")) {
    *nElemNodes = 3;
    *etype = Tr;
  }
  else if(!strcmp(etype_str,"Qd")) {    
    *nElemNodes = 4;
    *etype = Qd;
  }
  else if(!strcmp(etype_str,"Tt")) {
    *nElemNodes = 4;
    *etype = Tt;
  }
  else if(!strcmp(etype_str,"Hx")) {
    *nElemNodes = 8;
    *etype = Hx;
  }
  else {
    fprintf(stderr, "Error: Unsupported element type!\n");
    exit(1);
  }
  
  elems = (int**) malloc((*nElemNodes+1)*sizeof(int*));
  for(int i=0; i<*nElemNodes+1; i++)
    elems[i] = (int*) malloc(*nelems*sizeof(int));
  
  for(int i=0; i<*nelems; i++) {
    for(int j=0; j<*nElemNodes+1; j++)
      fscanf(f, "%d", &elems[j][i]);
    fscanf(f, "\n");
    fgets(etype_str, 3, f);
  }
  fclose(f);

  return elems;
}

double** readLonFile(char* meshname, int *nLon, int nelems) {
  double **lon;
  char lon_file[buffersize];
  char c;
  
  printf("Reading fibre file...\n");
  
  strcpy(lon_file,meshname);
  strcat(lon_file,".lon");
  
  FILE *f = fopen(lon_file,"r");
  if(!f) {
    printf("readLonFile: Could not open file %s\n", lon_file);
    exit(1);
  }
  
  // number of vectors
  fscanf(f,"%d\n",nLon);  
  
  if(*nLon < 2) {
    lon = (double**) malloc(3*sizeof(double*));
    for(int i=0; i<3; i++)
      lon[i] = (double*) malloc(nelems*sizeof(double));

    for(int i=0; i<nelems; i++)
      fscanf(f, "%lf %lf %lf\n", &lon[0][i], &lon[1][i], &lon[2][i]);
  }
  else {
    lon = (double**) malloc(6*sizeof(double*));
    for(int i=0; i<6; i++)
      lon[i] = (double*) malloc(nelems*sizeof(double));
    
    for(int i=0; i<nelems; i++)
      fscanf(f, "%lf %lf %lf %lf %lf %lf\n", &lon[0][i], &lon[1][i], &lon[2][i], &lon[3][i], &lon[4][i], &lon[5][i]);
  }
  
  fclose(f);

  return lon;
}

double* readActivationFile(char* actfile, int nnodes, int mode) {
  double tact, *act;
  int ind;
  
  printf("Reading data file...\n");
  
  FILE *f = fopen(actfile,"r");
  if(!f) {
    printf("readActivationFile: Could not open file %s\n", actfile);
    exit(1);
  }

  act = (double*) malloc(nnodes*sizeof(double));
  
  if(mode)
    for(int i=0; i<nnodes; i++)
      fscanf(f, "%lf\n", &act[i]);
  else
    for(int i=0; i<nnodes; i++) {
      fscanf(f, "%d %lf\n", &ind, &tact);
      act[ind] = tact;
    }

  fclose(f);
  return act;
}

double* readActivationFilePls(char* actfile, int nnodes, int mode, int npls) {
  double tact, *act;
  int ind, begin, end;
  
  FILE *f = fopen(actfile,"r");
  if(!f) {
    printf("readActivationFilePls: Could not open file %s\n", actfile);
    exit(1);
  }

  act = (double*) malloc(nnodes*sizeof(double));
  
  begin = nnodes*(npls-1);
  end = nnodes*npls;
  if(mode) {
    for(int i=0; i<begin; i++) // skip previous pulses: http://stackoverflow.com/questions/2799612/how-to-skip-the-first-line-when-fscanning-a-txt-file
      fscanf(f, "%*[^\n]\n");
    for(int i=begin; i<end; i++) // read only last activation
      fscanf(f, "%lf\n", &act[i-begin]);
  }
  else {
    for(int i=0; i<begin; i++) // skip previous pulses
      fscanf(f, "%*[^\n]\n");
    for(int i=begin; i<end; i++) { // read only last activation
      fscanf(f, "%d %lf\n", &ind, &tact);
      act[ind] = tact;
    }
  }

  fclose(f);
  return act;
}

int readWriteActivationFileAllPls(char* actfile, int nnodes, int npls) {
  double tact, **act;
  int ind;
  FILE *f1;
  char nactfile[buffersize];
  
  f1 = fopen(actfile,"r");
  if(!f1) {
    printf("readWriteActivationFileAllPls: Could not open file %s\n", actfile);
    exit(1);
  }

  act = (double**) malloc(npls*sizeof(double));
  for(int i=0; i<npls; i++)
    act[i] = (double*) malloc(nnodes*sizeof(double));
  
  // read all pulses and store act acording to node number
  for(int i=0; i<npls; i++) {
    for(int j=0; j<nnodes; j++) {
      fscanf(f1, "%d %lf\n", &ind, &tact);
      act[i][ind] = tact;
    }
  }

  fclose(f1);
  
  // write one file for each pulse
  for(int i=0; i<npls; i++) {
    sprintf(nactfile,"%s-stim%d.dat",actfile,i);  
    f1 = fopen(nactfile,"w");
    if(!f1) {
      printf("readWriteActivationFileAllPls: Could not open file %s\n", nactfile);
      exit(1);
    }
    for(int j=0; j<nnodes; j++) {
      fprintf(f1,"%f\n",act[i][j]);
    }
    fclose(f1);
  }
  
  // free memory
  for(int i=0; i<npls; i++)
    free(act[i]);
  free(act);  
  
  return 0;
}

int* readVtxFile(char* meshname, int *nnodes) {
  char pts_file[buffersize], type[10];
  int *nodes;
  
  strcpy(pts_file,meshname);
  strcat(pts_file,".vtx");
  
  FILE *f = fopen(pts_file,"r");
  if(!f) {
    printf("readVtxFile: Could not open file %s\n", pts_file);
    exit(1);
  }
  fscanf(f, "%d\n%s\n", nnodes, type);
  
  nodes = (int*) malloc(*nnodes*sizeof(int));

  for(int i=0; i<*nnodes; i++)
      fscanf(f, "%d\n", &nodes[i]);
  
  fclose(f);

  return nodes;
}


int writePointsFile(char* meshname, double **pts, int nnodes) {
  char pts_file[buffersize];
  
  strcpy(pts_file,meshname);
  strcat(pts_file,".pts");
  
  FILE *f = fopen(pts_file,"w");
  if(!f) {
    printf("writePointsFile: Could not open file %s\n", pts_file);
    exit(1);
  }
  fprintf(f, "%d\n", nnodes);

  for(int i=0; i<nnodes; i++)
      fprintf(f, "%lf %lf %lf\n", pts[0][i], pts[1][i], pts[2][i]);
  
  fclose(f);

  return 1;
}

int writeElemFile(char* meshname, int **elems, int nelems, int nElemNodes, int etype) {
  char elem_file[buffersize], etype_str[3];
  
  strcpy(elem_file,meshname);
  strcat(elem_file,".elem");
  
  FILE *f = fopen(elem_file,"w");
  if(!f) {
    printf("writeElemFile: Could not open file %s\n", elem_file);
    exit(1);
  }
  fprintf(f, "%d\n", nelems);
  
  switch(etype) {
    case Tr:
      strcpy(etype_str,"Tr");
      break;
    case Qd:
      strcpy(etype_str,"Qd");
      break;
    case Tt:
      strcpy(etype_str,"Tt");
      break;
    case Hx:
      strcpy(etype_str,"Hx");
      break;
    default: 
      fprintf(stderr, "Error: Unsupported element type!\n");
      exit(1);
  }

  for(int i=0; i<nelems; i++) {
    fprintf(f,"%s", etype_str);
    for(int j=0; j<nElemNodes+1; j++)
      fprintf(f, " %d", elems[j][i]);
    fprintf(f,"\n");
  }
  fclose(f);

  return 1;
}

int writeLonFile(char* meshname, double **lon, int nLon, int nelems) {
  char lon_file[buffersize];
  
  strcpy(lon_file,meshname);
  strcat(lon_file,".lon");
  
  FILE *f = fopen(lon_file,"w");
  if(!f) {
    printf("writeLonFile: Could not open file %s\n", lon_file);
    exit(1);
  }
  fprintf(f, "%d\n", nLon);

  for(int i=0; i<nelems; i++) {
    if(nLon < 2)
      fprintf(f, "%lf %lf %lf\n", lon[0][i], lon[1][i], lon[2][i]);
    else
      fprintf(f, "%lf %lf %lf %lf %lf %lf\n", lon[0][i], lon[1][i], lon[2][i], lon[3][i], lon[4][i], lon[5][i]);
  }
  fclose(f);

  return 1;
}

int writeVtxFile(char* meshname, int nnodes, int *nodes) {
  char pts_file[buffersize];
  
  strcpy(pts_file,meshname);
  strcat(pts_file,".vtx");
  
  FILE *f = fopen(pts_file,"w");
  if(!f) {
    printf("writeVtxFile: Could not open file %s\n", pts_file);
    exit(1);
  }
  fprintf(f, "%d\nextra\n", nnodes);

  for(int i=0; i<nnodes; i++)
      fprintf(f, "%d\n", nodes[i]);
  
  fclose(f);

  return 1;
}

int writeMap(char* mapname, int nelems, int *map) {
  FILE *f = fopen(mapname,"w");
  if(!f) {
    printf("writeMap: Could not open file %s\n", mapname);
    exit(1);
  }
  fprintf(f, "%d\n", nelems);
  
  for(int i=0; i<nelems; i++) {
    fprintf(f, "%d\n", map[i]);
  }
  fclose(f);
  
  return 0;
}

int writeActivationFile(char *actfile, double *act, int nnodes) {
  FILE *f = fopen(actfile,"w");
  if(!f) {
    printf("writeActivationFile: Could not open file %s\n", actfile);
    exit(1);
    return 0;
  }
  
  for(int i=0; i<nnodes; i++)
    fprintf(f,"%lf\n", act[i]);
  
  fclose(f);
  return 1;
}

unsigned short* readNRRD(char *imfile, NRRDheader_ *header) {
  char str[buffersize], *pch;
  int size;
  FILE *f;
    
  // read image file
  f = fopen(imfile,"r");
  if(!f) {
    printf("readNRRD: Could not open file %s\n", imfile);
    exit(1);
    return 0;
  }  
  
  // read dimension, sizes, space directions and origin
  skipNRRDcomment(f,str);
  if(strcmp(str,"NRRD0004")==0) {
    printf("File not in NRRD0004 format\n");
    exit(1);
  }
  
  // type
  skipNRRDcomment(f,str);

  // dimension
  skipNRRDcomment(f,str);
  pch = strtok(str," :"); // discard name of field
  pch = strtok(NULL, " :"); // get dimension
  header->dimension = atoi(pch);
  
  // space
  skipNRRDcomment(f,str);
  
  // sizes
  header->sizes = (int*)malloc(header->dimension*sizeof(int));
  skipNRRDcomment(f,str);
  pch = strtok(str," :"); 
  for(int i=0; i<header->dimension; i++) {
    pch = strtok(NULL, " :"); 
    header->sizes[i] = atoi(pch);
  }

  // space directions
  header->spacedirections = (double*)malloc(header->dimension*header->dimension*sizeof(double));
  skipNRRDcomment(f,str);
  pch = strtok(str," :,()"); 
  pch = strtok(NULL," :,()"); 

  for(int i=0; i<header->dimension*header->dimension; i++) {
    pch = strtok(NULL, " :,()"); 
    header->spacedirections[i] = atof(pch);
  }
  
  // kinds
  skipNRRDcomment(f,str);
  // endian
  skipNRRDcomment(f,str);
  // encoding
  skipNRRDcomment(f,str);
  
  // origin
  header->origin = (double*)malloc(header->dimension*sizeof(double));
  skipNRRDcomment(f,str);
  pch = strtok(str," :,()"); 
  pch = strtok(NULL," :,()");

  for(int i=0; i<header->dimension; i++) {
    pch = strtok(NULL, " :,()"); 
    header->origin[i] = atof(pch);
  }
  
  // skip useless lines
  skipNRRDgarbage(f,str);
  
  // get data size
  size=1;
  for(int i=0; i<header->dimension; i++)
    size *= header->sizes[i];
  
  // needs to be char*!!!
  char *data = (char*)malloc(size*sizeof(unsigned short)); 
  
  // read data
  fread(data, sizeof(unsigned short), size, f);
  
  return (unsigned short*)data;
}

int skipNRRDcomment(FILE *fid, char *curstr) {
  char *str, *pch;

  fgets(curstr, buffersize , fid);
  
  if(curstr != NULL) {
    pch=strchr(curstr,'#');
    while(pch != NULL) {
      fgets(curstr, buffersize , fid);
      pch=strchr(curstr,'#');
    }
  }
}
    
int skipNRRDgarbage(FILE *fid, char *curstr) {
  char *str, *pch;

  fgets(curstr, buffersize , fid);

  while(strcmp(curstr,"\n")!=0) {
    fgets(curstr, buffersize , fid);
  }
}

double** getNRRDtransMat(NRRDheader_ header){
  
  // allocate (dim+1)x(dim+1) matrix
  double **mat = (double**)malloc((header.dimension+1)*sizeof(double*));

  for(int i=0; i<header.dimension+1; i++)
    mat[i] = (double*)malloc((header.dimension+1)*sizeof(double));
  
  for(int i=0; i<header.dimension; i++) {
    for(int j=0; j<header.dimension; j++) {
      mat[j][i] = header.spacedirections[j+i*header.dimension];
    }
  }

  for(int i=0; i<header.dimension; i++)
    mat[header.dimension][i] = 0.0;

  for(int i=0; i<header.dimension; i++)
    mat[i][header.dimension] = header.origin[i];

  mat[header.dimension][header.dimension] = 1.0;
  
//  for(int i=0; i<header.dimension+1; i++)
//    printf("%f %f %f %f\n", mat[i][0], mat[i][1], mat[i][2], mat[i][3]);
  
  return mat;
}
    
int getNRRDres(NRRDheader_ *header) {
  double x,y,z;

  header->resolutions = (double*)malloc(header->dimension*sizeof(double));

  for(int i=0; i<header->dimension; i++) {
    x = header->spacedirections[i*header->dimension];
    y = header->spacedirections[i*header->dimension+1];
    z = header->spacedirections[i*header->dimension+2];

    header->resolutions[i] = vectorMagnitude3D(x,y,z);
  }
  
  return 0;
}  

int getNnodes(char* mesh_file) {
  int nnodes;
  FILE *f = fopen(mesh_file,"r");
  if(!f) {
    printf("getNnodes: Could not open file %s\n", mesh_file);
    exit(1);
  }
  fscanf(f, "%d\n", &nnodes);
  fclose(f);
  
  return nnodes;
}

int getNelems(char* mesh_file) {
  int nelems;
  FILE *f = fopen(mesh_file,"r");
  if(!f) {
    printf("getNnodes: Could not open file %s\n", mesh_file);
    exit(1);
  }
  fscanf(f, "%d\n", &nelems);
  fclose(f);
  
  return nelems;
}

int checkOpenFile(char * filename) {
  int pass = 1;
  FILE *fptr;
  
  if (!(fptr = fopen(filename, "r" ))) {
    printf("checkOpenFile: File %s could not be opened.\n", filename);
    pass = 0;
  }
  else {
    fseek(fptr, 0, SEEK_END);
    unsigned long len = (unsigned long)ftell(fptr);
    if (len == 0) {  //check if the file is empty or not.
      printf("File %s is empty.\n", filename);
      pass = 0;
      fclose(fptr);
    }
  }
  return pass;
}

int createArcherPBSscript(int nprocs, char *jobname, char* scriptname, char* command) {
  char script[buffersize];
  int ncnodes; // # of compute nodes
  
  ncnodes = nprocs/24; //must be multiple of 24
    
  FILE *f = fopen(scriptname,"w");
  if(!f) {
    printf("createArcherPBSscript: Could not open file %s\n", scriptname);
    exit(1);
  }  

  strcat(script,"#!/bin/bash --login\n");
  sprintf(script,"%s#PBS -N %s\n",script,jobname);
  strcat(script,"#PBS -M caroline.mendonca_costa@kcl.ac.uk\n");
  strcat(script,"#PBS -l walltime=24:00:00\n");
  strcat(script,"#PBS -A e348\n");
  sprintf(script,"%s#PBS -l select=%d:bigmem=true\n\n", script, ncnodes);

  strcat(script,"module load cray-petsc\n");
  strcat(script,"module swap cray-mpich\n\n");

  strcat(script,"export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)\n");
  strcat(script,"cd $PBS_O_WORKDIR\n\n");
  
  strcat(script,"export OMP_NUM_THREADS=1\n\n");
  
  // add CARP command line
  strcat(script,command); 
  
  // add PT options
  strcat(script, " -parab_use_pt 1");  
  strcat(script, " -ellip_use_pt 1");
  strcat(script, " -parab_options_file /work/e348/e348/cmcosta/pt_para_amg");
  strcat(script, " -ellip_options_file /work/e348/e348/cmcosta/pt_ell_amg");
 
  fprintf(f,"%s",script);
  fclose(f);  
}

int createTomPBSscript(int nprocs, char *jobname, char* scriptname, char* command) {
  char script[buffersize];
  
  FILE *f = fopen(scriptname,"w");
  if(!f) {
    printf("createTomPBSscript: Could not open file %s\n", scriptname);
    exit(1);
  }
  
  strcat(script,"#!/bin/bash\n");
  sprintf(script,"%s#PBS -N %s\n",script,jobname);
  strcat(script,"#PBS -M caroline.mendonca_costa@kcl.ac.uk\n");
  strcat(script,"#PBS -m bea\n");
  strcat(script,"#PBS -l walltime=150:00:00\n");
  sprintf(script,"%s#PBS -l select=1:ncpus=%d:mpiprocs=%d:mem=128gb\n", script, nprocs, nprocs);
  strcat(script,"#PBS -l place=pack:shared\n\n");

  strcat(script,"#Modules\n");
  strcat(script,"source /etc/profile.d/modules.sh\n");
  strcat(script,"module purge\n");
  strcat(script,"module load mpt/2.15\n");
  strcat(script,"module load intel/composerxe/13.1.1\n");
  strcat(script,"ulimit -s unlimited\n\n");

  strcat(script,"# (Intel compilers)\n");
  strcat(script,"# MPT MPI settings (enabling NUMA mode and logging)\n");
  strcat(script,"export MPI_DSM_DISTRIBUTE=1\n");
  strcat(script,"export MPI_DSM_VERBOSE=1\n\n");

  strcat(script,"# launch MPI job\n");
  strcat(script,"cd $PBS_O_WORKDIR\n\n");
  
  // add CARP command line
  strcat(script,command);

  // add PT options
  strcat(script, " -parab_use_pt 1");  
  strcat(script, " -ellip_use_pt 1");
  strcat(script, " -parab_options_file /home/cmc16/pt_para_amg");
  strcat(script, " -ellip_options_file /home/cmc16/pt_ell_amg");
 
  fprintf(f,"%s",script);
  fclose(f);  
}

#endif
