#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <mpi.h>
#include "dbg.h"
#include <unistd.h>

/*

A{sizex X size} * B{size X sizey}

*---->y (size)   *---->y (sizey) 
|                |
|              * |
v                |
x                |
                 v
(sizex)        (size)
*/

void init_matrix(double *matrix, int sizex, int sizey) {
	double **matr = (double**)malloc(sizex * sizeof(double*));
	for (int i = 0; i < sizex; i++) {
		matr[i] = matrix + i * sizey;
	}
	for (int i = 0; i < sizex; i++) {
		for (int j = 0; j < sizey; j++) {
			matr[i][j] = rand() % 10;
		}
	}
	free(matr);
}

void print_m(double *matr, int sizex, int sizey) {
	for (int i = 0; i < sizex; i++) {
		for (int j = 0; j < sizey; j++) {
			printf("%.0f ", matr[i * sizey + j]);
		}
		printf("\n");
	}
}

/*void mul(double *matr1, double *matr2, int sizex, int size, int sizey, double *dst) {
	// fix indexes
	for (int i = 0; i < sizex; i++) {
		for (int j = 0; j < sizey; j++) {
			for (int k = 0; k < size; k++) {
				dst[i][j] = matr1[i][k] * matr2[k][j];
			}
		}
	}
}*/

int size_by_rank(int size, int proc_num, int rank) {
	int part_size = size / proc_num;
	int rest = size % proc_num;
	if (rank < rest) {
		return part_size + 1;
	}
	else {
		return part_size;
	}
}

#define MIN(a, b) ((a) < (b) ? (a) : (b))

int offset_by_rank(int size, int proc_num, int rank) {
	int part_size = size / proc_num;
	int rest = size % proc_num;
	return part_size * rank + MIN(rank, rest);
}

#define X 0
#define Y 1


void create_subgrids(MPI_Comm grid, MPI_Comm *subgridx, MPI_Comm *subgridy) {
	int remain_dims[2];
	remain_dims[0] = true; remain_dims[1] = false;
	MPI_Cart_sub(grid, remain_dims, subgridx);

	remain_dims[0] = false; remain_dims[1] = true;
	MPI_Cart_sub(grid, remain_dims, subgridy);
}

/*int get_xy_leaders(int *subgx_leader, int *subgy_leader, int dimx, int dimy, int coords[2], MPI_Comm subgrid) {
	int *subgx_leaders = (int*)malloc(dimx * sizeof(int));
	int *subgy_leaders = (int*)malloc(dimy * sizeof(int));
	int subg_rank;

	int *tmpx = (int*)malloc(2 * dimx * sizeof(int));
	if (coords[Y] == rooty) { // I'm in leading subgridx
		MPI_Comm_rank(subgridx, &subg_rank);
		int foo[2] = { 0, 0 };
		foo[0] = subg_rank; foo[1] = coords[X];
		
		MPI_Allgather(foo, 2, MPI_INT, tmpx, 2, MPI_INT, subgridx);
		for (int i = 0; i < 2 * dimx; i += 2) {
			subgx_leaders[tmpx[i + 1]] = tmpx[i];
		}
	}
	MPI_Bcast(subgx_leaders, dimx, MPI_INT, 0, MPI_COMM_WORLD);
	if (grid_rank == root) {
		for (int i = 0; i < 2 * dimx; i++) {
			printf("subgridx rank: %d (%d), coord: %d\n", tmp[i], tmp[i + 1]);
			i++;
		}
	}
	int *tmpy = (int*)malloc(2 * dimy * sizeof(int));
	if (coords[X] == rootx) { // I'm in leading subgridy
		MPI_Comm_rank(subgridy, &subg_rank);
		int foo[2] = { 0, 0 };
		foo[0] = subg_rank; foo[1] = coords[Y];

		MPI_Allgather(foo, 2, MPI_INT, tmpy, 2, MPI_INT, subgridy);
		for (int i = 0; i < 2 * dimy; i += 2) {
			subgy_leaders[tmpy[i + 1]]] = tmpy[i];
		}
	}
	MPI_Bcast(subgy_leaders, dimy, MPI_INT, 0, MPI_COMM_WORLD);
}*/

void distribute(double *matrA, double *matrB, int sizex, int size, int sizey, int proc_num, int rank) {
	/*lines_per_proc = (sizex + sizey) / proc_num;
	dimx = sizex / lines_per_proc;
	dimy = sizey / lines_per_proc;*/

	
	int dims[2] = { 0, 0 };
	MPI_Dims_create(proc_num, 2, dims);
	int dimx = dims[0];
	int dimy = dims[1];
	
	int periods[2] = { false, false };
	MPI_Comm grid;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, false, &grid);


	


	int grid_rank;
	MPI_Comm_rank(grid, &grid_rank);

	
	int coords[2] = { 0, 0 };
	MPI_Cart_coords(grid, grid_rank, 2, coords);
	// printf("i'm %d proc (%d in grid) (%d, %d)\n", rank, grid_rank, coords[0], coords[1]);


	// ----------scatter matrix A--------------
	// ----root
	int root = -1;
	int rootx = -1;
	int rooty = -1;
	if (rank == 0) {
		root = grid_rank;
		rootx = coords[X];
		rooty = coords[Y];
	}
	MPI_Bcast(&root, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&rootx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&rooty, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// ----sizes
	// ----offsets
	int *sizesA = (int*)malloc(dimx * sizeof(int));
	int *offsetsA = (int*)malloc(dimx * sizeof(int));

	for (int i = 0; i < dimx; i++) {
		sizesA[i] = size_by_rank(sizex, dimx, i)/* * size*/;
		offsetsA[i] = offset_by_rank(sizex, dimx, i)/* * size*/;
	}

	int *sizesB = (int*)malloc(dimy * sizeof(int));
	int *offsetsB = (int*)malloc(dimy * sizeof(int));

	for (int i = 0; i < dimy; i++) {
		sizesB[i] = size_by_rank(sizey, dimy, i)/* * size*/;
		offsetsB[i] = offset_by_rank(sizey, dimy, i)/* * size*/;
	}

	// ----subgrids
	MPI_Comm subgridx; // same y, different x
	MPI_Comm subgridy; // same x, different y
	create_subgrids(grid, &subgridx, &subgridy);



	
	// ----partA
	double *partA = (double*)malloc(sizesA[coords[X]] * size * sizeof(double));
	double *partB = (double*)malloc(sizesB[coords[Y]] * size * sizeof(double));

	// ----part_size
	int partA_size = sizesA[coords[X]];
	int partB_size = sizesB[coords[Y]];

	
	// ----scattering matrix A----
	if (coords[Y] == rooty) { // if you're the process in root's subgridx
		MPI_Datatype row_type;
		MPI_Type_contiguous(size, MPI_DOUBLE, &row_type);
		MPI_Type_commit(&row_type);

		MPI_Scatterv(matrA, sizesA, offsetsA, row_type, partA, partA_size, row_type, root, subgridx);

		MPI_Type_free(&row_type);
	}
	// ----scattering matrix B----
	if (coords[X] == rootx) {
		MPI_Datatype root_col_type;
		MPI_Type_vector(size, 1, sizey, MPI_DOUBLE, &root_col_type);
		MPI_Type_create_resized(root_col_type, 0, sizeof(double), &root_col_type);
		MPI_Type_commit(&root_col_type);

		MPI_Datatype recver_col_type;
		MPI_Type_vector(size, 1, sizesB[coords[Y]], MPI_DOUBLE, &recver_col_type);
		MPI_Type_create_resized(recver_col_type, 0, sizeof(double), &recver_col_type);
		MPI_Type_commit(&recver_col_type);

		MPI_Scatterv(matrB, sizesB, offsetsB, root_col_type, partB, partB_size, recver_col_type, root, subgridy);

		MPI_Type_free(&root_col_type);
		MPI_Type_free(&recver_col_type);
	}

	
	int *subgy_leaders = (int*)malloc(dimx * sizeof(int));
	int *subgx_leaders = (int*)malloc(dimy * sizeof(int));
	int subg_rank;

	int *tmpy = (int*)malloc(2 * dimx * sizeof(int));
	if (coords[Y] == rooty) { // I'm in leading subgridx
		MPI_Comm_rank(subgridy, &subg_rank);
		int foo[2] = { 0, 0 };
		foo[0] = subg_rank; foo[1] = coords[X];
		
		MPI_Allgather(foo, 2, MPI_INT, tmpy, 2, MPI_INT, subgridx);
		for (int i = 0; i < 2 * dimx; i += 2) {
			subgy_leaders[tmpy[i + 1]] = tmpy[i];
		}
	}
	MPI_Bcast(subgy_leaders, dimx, MPI_INT, 0, MPI_COMM_WORLD);
	/*if (grid_rank == root) {
		for (int i = 0; i < 2 * dimx; i++) {
			printf("subgridx rank: %d (%d), coord: %d\n", tmp[i], tmp[i + 1]);
			i++;
		}
	}*/
	int *tmpx = (int*)malloc(2 * dimy * sizeof(int));
	if (coords[X] == rootx) { // I'm in leading subgridy
		MPI_Comm_rank(subgridx, &subg_rank);
		int foo[2] = { 0, 0 };
		foo[0] = subg_rank; foo[1] = coords[Y];

		MPI_Allgather(foo, 2, MPI_INT, tmpx, 2, MPI_INT, subgridy);
		for (int i = 0; i < 2 * dimy; i += 2) {
			subgx_leaders[tmpx[i + 1]] = tmpx[i];
		}
	}
	MPI_Bcast(subgx_leaders, dimy, MPI_INT, 0, MPI_COMM_WORLD);

	// ranks of subgrids leaders (who placed at the same subgrid with root)
	int subgx_leader = subgx_leaders[coords[Y]];
	int subgy_leader = subgy_leaders[coords[X]];


	

	MPI_Bcast(partA, sizesA[coords[X]] * size, MPI_DOUBLE, subgy_leader, subgridy);
	MPI_Bcast(partB, sizesB[coords[Y]] * size, MPI_DOUBLE, subgx_leader, subgridx);

	/*if (coords[Y] == rooty) {
		printf("I'm in leading subgridx, my rank - %d\n", grid_rank);
	}*/

	if (rank == 0) {
		printf("sizesB: ");
		for (int i = 0; i < dimy; i++) {
			printf("%d ", sizesB[i]);
		}
		printf("\noffsetsB: ");
		for (int i = 0; i < dimy; i++) {
			printf("%d ", offsetsB[i]);
		}
	}
	printf("\n");
	for (int i = 0; i < dimx; i++) {
		sleep(1);
		if (coords[X] == i && coords[Y] == rooty) {
			print_part(partA, size, sizesA[coords[X]]);
			printf("\n");
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	for (int i = 0; i < dimy; i++) {
		sleep(1);
		if (coords[X] == rootx && coords[Y] == i) {
			print_part(partB, sizesB[coords[Y]], size);
			printf("\n");
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	printf("parts\n");
	/*if (coords[X] == 2) {
		for (int i = 0; i < dimy; i++) {
			sleep(1);
			if (i == coords[Y]) {
				// print_part(partB, sizesB[coords[Y]], size);
				print_part(partA, size, sizesA[coords[X]]);
				printf("\n");
			}
			MPI_Barrier(subgridy);
		}
	}*/

	for (int x = 0; x < dimx; x++) {
		for (int y = 0; y < dimy; y++) {
			sleep(1);
			if (x == coords[X] && y == coords[Y]) {
				printf("\nx, y: (%d, %d)\n", x, y);
				print_part(partB, sizesB[coords[Y]], size);
				// print_part(partA, size, sizesA[coords[X]]);
				printf("\n");
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}

	}
	
}

int main(int argc, char *argv[]) {
	int err_code;

	if ((err_code = MPI_Init(&argc, &argv)) != 0) {
		printf("error\n");
		return err_code;
	}

	int sizex = 5;
	int size = 10;
	int sizey = 7;

	int rank = 0, comm_size = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);


	double *matrA = NULL;
	double *matrB = NULL;
	double *partA = NULL;
	double *partB = NULL;
	if (rank == 0) {
		matrA = (double*)malloc(sizex * size * sizeof(double));
		matrB = (double*)malloc(size * sizey * sizeof(double));
		init_matrix(matrA, sizex, size);
		init_matrix(matrB, size, sizey);
		printf("\nA\n");
		print_m(matrA, sizex, size);
		printf("\nB\n");
		print_m(matrB, size, sizey);
		printf("\nlin B\n");
		for (int i = 0; i < size * sizey; i++) {
			printf("%.0f ", matrB[i]);
		}
		printf("\n");
	}

	

	distribute(matrA, matrB, sizex, size, sizey, comm_size, rank);

	/*something_between()??

	combine();???*/

	MPI_Finalize();
	return 0;
}