#include <mpi.h>   
#include <stdio.h>
#include <string.h>

int main(int argc, char **argv) {
  
  int rank; 
  int proc;
  int size, i,t;
  char greeting[200];
  MPI_Status status;
  
  MPI_Init(&argc, &argv); /*initialize MPI*/
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); /*RANK OF THIS PROCESSOR*/
  MPI_Comm_size(MPI_COMM_WORLD, &size); /*TOTAL NUMBER OF PROCESSORS*/
  
  
  sprintf(greeting, "Hello Akuma, reporting from processor %d of %d\n", rank, size);
  
  
  if (rank ==0) { /*rank 0 prints greeting*/
    fputs(greeting, stdout);
    for (proc = 1; proc < size; proc++){ /*loop through processes from 0 to size of processes */
      
      MPI_Recv(greeting, sizeof(greeting), MPI_BYTE, proc, 1, MPI_COMM_WORLD, &status); 
      fputs (greeting, stdout); /* rank 0 receives greeting and displays it on console */
      
    }
  }
  else {
    MPI_Send(greeting, strlen(greeting)+1, MPI_BYTE, 0,1,MPI_COMM_WORLD);
  } /* block communication till message is received */
  
  
  
  if (rank == 0) printf("Finish!\n"); /* print this message once all messages have been received */

  MPI_Finalize();  /* Exit */
  
}
