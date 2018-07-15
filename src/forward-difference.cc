#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <mpi.h>

using namespace std;

int main(int argc, char ** argv) {
  int commrank, nproc;
  MPI_Request sendrequest,recvrequest;
  MPI_Status status;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
  double temp;

  int N = 30;
  int nloc = N/nproc;
  double h = 1.0/(N-1);

  vector<vector<double>> u;
  u.resize(nloc+1,vector<double>(nloc+1));

  string name = "data/sol";
  ofstream outfstr[nproc];

  //compute the value
  for(int i=0;i<nproc;i++)
    for(int j=0;j<nloc+1;j++)
      u[j][i] = sin(2*M_PI*(j+(i*nloc))*h);

  //send the ghost values
  if(commrank == 0) {
    temp = u[0][commrank];
    MPI_Isend(&temp, 1, MPI_DOUBLE, commrank+1, 0, MPI_COMM_WORLD,&sendrequest);
    MPI_Wait(&sendrequest,&status);
  }

  else if(commrank == nproc-1) {
    temp = u[nloc][commrank];
    MPI_Irecv(&temp, 1, MPI_DOUBLE, commrank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &recvrequest);
    MPI_Wait(&recvrequest,&status);
    u[nloc][commrank] = temp;
  }

  else {
    temp = u[0][commrank];
    MPI_Isend(&temp, 1, MPI_DOUBLE, commrank+1, 0, MPI_COMM_WORLD,&sendrequest);
    MPI_Wait(&sendrequest,&status);
    MPI_Irecv(&temp, 1, MPI_DOUBLE, commrank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &recvrequest);
    MPI_Wait(&recvrequest,&status);
    u[nloc][commrank] = temp;
  }

  // Write the files
  for(int i=0;i<nproc;i++) {
    outfstr[i].open(name + to_string(i) + ".dat");
    for(int j=0;j<nloc+1;j++)
      outfstr[i] << (j+(i*nloc))*h << " " << u[j][i] << endl;
  }

  MPI_Finalize();
  return 0;
}
