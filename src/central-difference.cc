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
  double sigma = -0.25;

  vector<vector<double>> u_prev;
  u_prev.resize(nloc+2,vector<double>(nloc+2));
  vector<vector<double>> u_next;
  u_next.resize(nloc+2,vector<double>(nloc+2));

  string name = "data/sol";
  ofstream outfstr[nproc][100];

  //compute the value
  for(int i=0;i<nproc;i++)
    for(int j=1;j<nloc+2;j++)
      u_prev[j][i] = sin(2*M_PI*(j-1 +(i*nloc))*h);

  for(int it=0; it < 100;it++)  {
    // Open files for writing
    for(int i=0;i<nproc;i++)
      outfstr[i][it].open(name + to_string(i) + "-" + to_string(it)+ ".dat");
    if(commrank == 0) {
        temp = u_prev[nloc][commrank];
        MPI_Isend(&temp, 1, MPI_DOUBLE, commrank+1, 0, MPI_COMM_WORLD,&sendrequest);
        MPI_Wait(&sendrequest,&status);
        MPI_Irecv(&temp, 1, MPI_DOUBLE, nproc-1 , 0,MPI_COMM_WORLD, &recvrequest);
        MPI_Wait(&recvrequest, &status);
        u_prev[0][commrank] = temp;
        temp = u_prev[1][commrank];
        MPI_Isend(&temp, 1, MPI_DOUBLE, nproc-1, 0, MPI_COMM_WORLD,&sendrequest);
        MPI_Wait(&sendrequest,&status);
        MPI_Irecv(&temp, 1, MPI_DOUBLE, commrank+1 , 0,MPI_COMM_WORLD, &recvrequest);
        MPI_Wait(&recvrequest, &status);
        u_prev[nloc][commrank] = temp;
        for(int i=1; i<nloc+1;i++) {
          u_next[i][commrank] = u_prev[i][commrank] + (sigma/2)*u_prev[i-1][commrank] - (sigma/2)*u_prev[i+1][commrank];
          outfstr[commrank][it] << (i-1 + (commrank*nloc))*h << " " << u_next[i][commrank] << endl;
        }
    }

    else if(commrank == nproc-1) {
        temp = u_prev[nloc-1][commrank];
        MPI_Isend(&temp, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&sendrequest);
        MPI_Wait(&sendrequest,&status);
        MPI_Irecv(&temp, 1, MPI_DOUBLE, commrank-1, 0, MPI_COMM_WORLD, &recvrequest);
        MPI_Wait(&recvrequest,&status);
        u_prev[0][commrank] = temp;
        temp = u_prev[0][commrank];
        MPI_Isend(&temp, 1, MPI_DOUBLE, commrank-1, 0, MPI_COMM_WORLD,&sendrequest);
        MPI_Wait(&sendrequest,&status);
        MPI_Irecv(&temp, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &recvrequest);
        MPI_Wait(&recvrequest,&status);
        u_prev[nloc][commrank] = temp;
        for(int i=1; i<nloc+1;i++) {
          u_next[i][commrank] = u_prev[i][commrank] + (sigma/2)*u_prev[i-1][commrank] - (sigma/2)*u_prev[i+1][commrank];
            outfstr[commrank][it] << (i-1 + (commrank*nloc))*h << " " << u_next[i][commrank] << endl;
        }
    }

    else {
        temp = u_prev[nloc][commrank];
        MPI_Isend(&temp, 1, MPI_DOUBLE, commrank+1, 0, MPI_COMM_WORLD,&sendrequest);
        MPI_Wait(&sendrequest,&status);
        MPI_Irecv(&temp, 1, MPI_DOUBLE, commrank-1, 0, MPI_COMM_WORLD, &recvrequest);
        MPI_Wait(&recvrequest,&status);
        u_prev[0][commrank] = temp;
        temp = u_prev[0][commrank];
        MPI_Isend(&temp, 1, MPI_DOUBLE, commrank-1, 0, MPI_COMM_WORLD,&sendrequest);
        MPI_Wait(&sendrequest,&status);
        MPI_Irecv(&temp, 1, MPI_DOUBLE, commrank+1, 0, MPI_COMM_WORLD, &recvrequest);
        MPI_Wait(&recvrequest,&status);
        u_prev[nloc][commrank] = temp;
        for(int i=1; i<nloc+1;i++) {
            u_next[i][commrank] = u_prev[i][commrank] + (sigma/2)*u_prev[i-1][commrank] - (sigma/2)*u_prev[i+1][commrank];
            outfstr[commrank][it] << (i-1 + (commrank*nloc))*h << " " << u_next[i][commrank] << endl;
        }
    }
    u_prev = u_next;
  }

  MPI_Finalize();
  return 0;
}
