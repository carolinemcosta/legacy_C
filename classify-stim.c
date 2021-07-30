#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include "utils.h"
#include "ioutils.h"

int main(int argc, char *argv[]){

  char actfilenc[buffersize], actfileprop[buffersize], actfileub[buffersize], nactfileub[buffersize], result[buffersize], ptsname[buffersize], *meshname, *simIDnc, *simIDprop, *simIDub;
       
  double *actnc, *actprop, *actub,
         nactnc, nactprop, nactub,
         minnc, minprop, minub,
	 errnc, errprop, rsmerrnc, rsmerrprop,
	 dicenc, diceprop;
  
  int **elems, nelems, nelemnodes, etype, nnodes, type, nactnodes, intnc, intprop;
  
  FILE *f1, *f2, *f3, *f4;

  // parameters
  meshname  = argv[1];
  simIDnc   = argv[2];
  simIDprop = argv[3];
  simIDub   = argv[4];

  // read number of nodes
  sprintf(ptsname,"%s.pts",meshname);
  nnodes = getNnodes(ptsname);  
  
  // these should be parameters
  // activation file with no-capture
  sprintf(actfilenc,"%s/init_acts_ACTs-thresh.dat",simIDnc);
  // activation file with normal propagation
  sprintf(actfileprop,"%s/init_acts_ACTs-thresh.dat",simIDprop);
  // activation file which we want to classify
  sprintf(actfileub,"%s/init_acts_ACTs-thresh.dat",simIDub);
  
  // open files
  f1 = openFile(actfilenc,1);
  f2 = openFile(actfileprop,1); 
  f3 = openFile(actfileub,1);
  
  // allocate vectors
  actnc   = (double*) malloc(nnodes*sizeof(double));
  actprop = (double*) malloc(nnodes*sizeof(double));
  actub   = (double*) malloc(nnodes*sizeof(double));
  
  intnc = 0; intprop = 0;
  nactnodes = 0;
  // read file and find min act time
  for(int i=0; i<nnodes; i++) {
    fscanf(f1, "%lf\n", &actnc[i]);
    fscanf(f2, "%lf\n", &actprop[i]);
    fscanf(f3, "%lf\n", &actub[i]);
    
    if(actprop[i]>0.) {
      // dice
      if((actnc[i]>0. && actub[i]>0.) || (actnc[i]<0. && actub[i]<0.))
        intnc++;
      nactnodes++;
    }
  }

  fclose(f1); fclose(f2); fclose(f3);

  dicenc = (2.*intnc)/(2.*nactnodes);
  
  // classify type of propagation - 5% error
  if(dicenc > 0.98) // no capture
    type = 0;
  else if(dicenc < 0.40) // normal propagation
    type = 2;
  else // uni-directional block
    type = 1;  
  
  sprintf(result,"%s/type.dat",simIDub);  
  f1 = openFile(result,0);
  fprintf(f1,"%d", type);
  fclose(f1);  
  
 
  return 0;
}
  
  
  
  
