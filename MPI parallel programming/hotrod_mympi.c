#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <gd.h>
#define MAXCOLORS 256
#define OWNER(index) ((nprocs*(index+1)-1)/nx)

/* Global variables */
int nx = NX;              /* number of discrete points including endpoints */
double m = M;             /* initial temperature of rod */
double k = K;             /* D*dt/(dx*dx) */
int nsteps = NSTEPS;      /* number of time steps */
double dx;                /* distance between two grid points: 1/(nx-1) */
double *u;                /* temperature function */
double *u_new;            /* temp. used to update u */
FILE *file;               /* file containing animated GIF */
gdImagePtr im, previm;    /* pointers to GIF images */
int *colors;              /* colors we will use */
int output_par[NX];       // color of parallel 


/* MPI Process in parallel algorithm */

void MPI_Process (int _rank) {
    $scope process_scope = $here;   
    MPI_Comm MPI_COMM_WORLD;

    int nx = -1;              /* number of discrete points including endpoints */
    double k = -1;            /* D*dt/(dx*dx) */
    int nsteps = -1;          /* number of time steps */
    int wstep = -1;           /* write frame every this many time steps */
    double *u;                /* temperature function */
    double *u_new;            /* temp. used to update u */

    int nprocs;    /* number of processes */
    int rank;      /* the rank of this process */
    int left;      /* rank of left neighbor */
    int right;     /* rank of right neighbor on torus */
    int nxl;       /* horizontal extent of one process */
    int first;     /* global index for local index 0 */
    int start;     /* first local index to update */
    int stop;      /* last local index to update */
    double *buf;   /* temp. buffer used on proc 0 only */
    int framecount = 0; //used only on process 0

int firstForProc(int rank) {
        return (rank*nx)/nprocs;	
    }

int countForProc(int rank) {
    int a;
    int b;
	
    a = firstForProc(rank+1);
    b = firstForProc(rank);
    return a-b;
  }

  /* init: initializes global variables. */
  void init() { 
    int i, j;
    int pos = 0; // position in input file

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); //communicator
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  //rank
    MPI_Barrier(MPI_COMM_WORLD); 
      start = MPI_Wtime();
    if (rank == 0) {
      nx = NX;
      k = K;
      nsteps = NSTEPS;
      wstep = WSTEP;
      printf("Diffusion1d (par) with nx=%d, k=%f, nsteps=%d, wstep=%d nprocs=%d\n",
	     nx, k, nsteps, wstep, nprocs);
    }
    
    MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nsteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&wstep, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // assert(nx>=nprocs);

    assert(k>0 && k<.5);
    assert(nx>=2);
    assert(nsteps>=1);
    first = firstForProc(rank);  
    nxl = countForProc(rank);

    if (first == 0 || nxl == 0)
        left = MPI_PROC_NULL;
    else
        left = OWNER(first-1);
    if (first+nxl >= nx || nxl == 0)
         right = MPI_PROC_NULL;
    else
        right = OWNER(first+nxl); 

     u = (double*)$malloc(process_scope, (nxl+2)*sizeof(double));
        assert(u != NULL);
     u_new = (double*)$malloc(process_scope, (nxl+2)*sizeof(double));
        assert(u_new != NULL);
     if (rank == 0) {
        buf = (double*)$malloc(process_scope, (1+nx/nprocs)*sizeof(double));
            for (i=1; i <= nxl; i++)
	            u[i] = u_init[pos++]; // set and increment temp position
            for (i=1; i < nprocs; i++) {
	            int count_i = countForProc(i);

	        for (j=0; j<count_i; j++)
	            buf[j] = u_init[pos++];
	MPI_Send(buf, count_i, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            }
     }
      
    else {
      buf = NULL;
      MPI_Recv(u+1, nxl, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (rank == OWNER(0)) {
      start = 2;
    }
    else {
      start = 1;
    }
    if (rank == OWNER(nx-1)) {
      stop = nxl - 1;
    }
    else {
      stop = nxl;
    }

  /*  if (rank == 0) {
        buf = (double*)malloc((1+nx/nprocs)*sizeof(double));
        assert(buf);
        file = fopen("./parout/out.gif", "wb");
        assert(file);
        colors = (int*)malloc(MAXCOLORS*sizeof(int));
        assert(colors);
    } 
    else {
        buf = NULL;
  } */

  }

 // update color frames 

void write_frame(int time) {
    if (rank != 0) {
      MPI_Send(u+1, nxl, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    } else {
      int source, i, j, count, global_index;
      MPI_Status status;
  
      global_index = 0;
      for (source = 0; source < nprocs; source++) {
	if (source != 0) {
	  MPI_Recv(buf, 1+nx/nprocs, MPI_DOUBLE, source, 0,
		   MPI_COMM_WORLD, &status);
	  MPI_Get_count(&status, MPI_DOUBLE, &count);
	} else {
	  for (i = 1; i <= nxl; i++) buf[i-1] = u[i];
	  count = nxl;
	}
	for (i = 0; i < count; i++) {
	  output_par[global_index] = colorOf(buf[i]);
	  global_index++;
	}
      }
    }
  }

/* exchange_ghost_cells: updates ghost cells using MPI communication */

  void exchange_ghost_cells() { //Prof. Siegel's ghost function
    MPI_Sendrecv(&u[1], 1, MPI_DOUBLE, left, 0,
		 &u[nxl+1], 1, MPI_DOUBLE, right, 0,
		 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&u[nxl], 1, MPI_DOUBLE, right, 0,
		 &u[0], 1, MPI_DOUBLE, left, 0,
		 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

/* update: updates u.  Uses ghost cells.  */
  
  void update() { //Prof Siegel's update function 
    int i;
    
    for (i = start; i <= stop; i++)
      u_new[i] = u[i] + k*(u[i+1] + u[i-1] - 2*u[i]);
    for (i = start; i <= stop; i++)
      u[i] = u_new[i];

      double *tmp = u_new; u_new = u; u = tmp;
  }

  
  void _main() {
    int q;

    MPI_Init(); //initialize all variables
    MPI_COMM_WORLD = MPI_Comm_create(process_scope, MPI_GCOMM_WORLD, _rank);

    init();
    write_frame(0); //begin frame writing 
    for (q = 1; q <= nsteps; q++) {
      exchange_ghost_cells();
      update();
      if (q%wstep==0) write_frame(q);
    }
    MPI_Barrier(MPI_COMM_WORLD); 
      end = MPI_Wtime(); 
    MPI_Finalize();
    MPI_Comm_destroy(MPI_COMM_WORLD);
    free(u);
    free(u_new);
    if (rank == 0)
      free(buf);
  }

  _main();
}

void MPI_par() {
  $proc procs[NPROCS];          // the MPI processes

  MPI_GCOMM_WORLD = CMPI_Gcomm_create($root, NPROCS); //create process instance
  for (int i=0; i<NPROCS; i++) procs[i] = $spawn MPI_Process(i);
  for (int i=0; i<NPROCS; i++) $wait(procs[i]);
  CMPI_Gcomm_destroy(MPI_GCOMM_WORLD);
}  

//run mpi process
void main() {   
  for (int i=0; i<NX; i++) ;
  for (int i=0; i<NSTEPS; i++) ;
  for (int i=0; i<WSTEP; i++) ;
  for (int i=0; i<NPROCS; i++) ;
  MPI_par();
  for (int i=0; i<NX; i++) {
    printf("output[%d] = %d\n", i, output_par[i]);
    
  }
}


