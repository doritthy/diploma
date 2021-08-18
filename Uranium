#include <stdio.h>
#include <math.h>
#include <mpi.h>
//Globalnie peremennie
#define Nx 80
#define Ny 41
#define Nz 40
//
#define nx 10
#define ny 41
#define nz 40
//
#define nw 3
//
FILE *inf;
//
double H[nx + 2][ny + 2][nz + 2];
double Hn[nx + 2][ny + 2][nz + 2];
double Q[nx + 2][ny + 2][nz + 2];
double Qp[nx + 2][ny + 2][nz + 2];
double Qr_z[nx + 2][ny + 2][nz + 2];
double Qr_o[nx + 2][ny + 2][nz + 2];
double U[nx + 2][ny + 2][nz + 2];
double V[nx + 2][ny + 2][nz + 2];
double W[nx + 2][ny + 2][nz + 2];
double C[nx + 2][ny + 2][nz + 2][3];
double Cn[nx + 2][ny + 2][nz + 2][3];
double D[nx + 2][ny + 2][nz + 2][3];
double buf[nx + 2][ny + 2][nz + 2];
double buf1[nx + 2][ny + 2][nz + 2][3];

/*double Dxx[nx+2][ny+2][nz+2];
double Dyy[nx+2][ny+2][nz+2];
double Dzz[nx+2][ny+2][nz+2];*/

double Lx, Ly, Lz; //razmeri vsei oblasti
double dx, dy, dz, dt, t, tt;
//
int tag_next_x, tag_prev_x;
int tag_next_y, tag_prev_y;
int tag_next_z, tag_prev_z;
//
int rank, size;
int np_x;	// number of processes along x,y,z  Nomer processora v dol' osei x,y,z
int p_x;	// place of current process(cube) along x,y,z Mestopolozhenie dannogo processora vdol' x,y,z
int next_x, prev_x;	// next and previous process by x
					// next and previous process by y
					//
double K, ro, g;
double qs, qs1, H0, Cm0, Cr0;
double h_w; // height of well visota skvazhini
double _D, teta, alpha, v1, v2, beta;
//
double max, max_all, eps, maxV, vel_all, D_all;
double a, ax, ay, az, axc, ayc, azc;
int n = 0;
double start, end;
//
MPI_Status status;
MPI_Request request;
//
// Function prototypes
void pass_D(void);
void pass_H(void);
void pass_C(void);
void Define_Vel(void);
void output_H(void);
int obtain_Hzhup(void);
int obtain_Htak(void);
void output_int(void);
int Define_Q(void);
void Define_D(void);
void obtain_Crmp(void);
void output_C(void);
void pass_Q(void);
void pass_U(void);
//
int main(int argc, char *argv[])
{
	int i, j, k, f;
	// Begin input ********************************************************************
	Lx = 40.;
	Ly = 20.;
	Lz = 20.;
	//
	qs = 100.;
	qs1 = 100.;
	K = 5.;
	ro = 1000.;
	g = 9.8;
	H0 = 10.*100000. / (ro * g);
	Cm0 = 0.15e+01;
	Cr0 = 0.20e+01;
	alpha = 0.1;
	teta = 0.4; //0.7; //0.5; //0.3;
	_D = 0.3; //0.01;
			  //
			  // End input **********************************************************************
	dx = Lx / Nx;
	dy = Ly / (Ny + 1);
	dz = Lz / (Nz + 1);
	dt = dx*dx / 2.;
	tt = 0.54*teta*dx*dx / 1000.;
	np_x = Nx / nx;
	beta = 0.01;
	v1 = 0.3426573; //1.; //98./79.5;
	v2 = 1. + v1; //159.5/79.5;
				  //
	eps = 1.0e-3;
	//
	//
	MPI_Init(&argc, &argv);
	start = MPI_Wtime();
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	//
	p_x = rank;
	//
	next_x = rank + 1;
	prev_x = rank - 1;
	//
	tag_next_x = 0; //p_x * next_x +         (p_y+1)*(p_z+1);
	tag_next_y = 0; //p_y * next_y + (p_x+1)*        (p_z+1);
	tag_prev_x = 0; //p_x * prev_x +         (p_y+1)*(p_z+1);
	tag_prev_y = 0; //p_y * prev_y + (p_x+1)*        (p_z+1);
					//
					//BEGIN DEFINE H=======================================================================================================
	ax = 1. / (dx*dx);
	ay = 1. / (dy*dy);
	az = 1. / (dz*dz);
	axc = 1. / dx;
	ayc = 1. / dy;
	azc = 1. / dz;
	//Initialization-----------------------------------------------------------
	for (i = 0; i <= nx + 1; i++)
		for (j = 0; j <= ny + 1; j++)
			for (k = 0; k <= nz + 1; k++)
			{
				H[i][j][k] = H0;
				Hn[i][j][k] = H0;
			}
	//output_int();
	//----------------------------------------------------------------------------
	Define_Q();
	pass_Q();
	output_int();
	do
	{
		n++;
		obtain_Hzhup();
		MPI_Barrier(MPI_COMM_WORLD);
		pass_H();
		obtain_Htak();
		MPI_Barrier(MPI_COMM_WORLD);
		pass_H();
		for (i = 0; i <= nx + 1; i++)
			for (j = 0; j <= ny + 1; j++)
				for (k = 0; k <= nz + 1; k++)
					H[i][j][k] = Hn[i][j][k];

		MPI_Allreduce(&max, &max_all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		if (rank == 1)
		{
			inf = fopen("inf.txt", "w");
			fprintf(inf, "n = %5d tt = %11.9f dt = %11.9f max = %11.9f\n", n, tt, dt, max_all);
			//printf("%5d %9.6f %9.6f\n",n,max,max_all);
			fclose(inf);
		}
	} while (max_all > eps);
	//END DEFINE H=========================================================================================================
	Define_Vel();
	pass_U();
	Define_D();
	pass_D();
	MPI_Allreduce(&maxV, &vel_all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	tt = 0.5*dx*dx / (3.*alpha*vel_all + teta*_D) / 10.;
	output_H();
	MPI_Barrier(MPI_COMM_WORLD);
	//BEGIN DEFINE C=======================================================================================================
	for (i = 0; i <= nx + 1; i++)
		for (j = 0; j <= ny + 1; j++)
			for (k = 0; k <= nz + 1; k++)
			{
				C[i][j][k][0] = 0.;
				C[i][j][k][1] = Cm0;
				C[i][j][k][2] = 0.;
				Cn[i][j][k][0] = 0.;
				Cn[i][j][k][1] = Cm0;
				Cn[i][j][k][2] = 0.;
			}
	for (n = 0; n <= 100; n++)
	{
		//init=========================

		//initEnd=====================
		obtain_Crmp();
		MPI_Barrier(MPI_COMM_WORLD);
		pass_C();
		//obtain_Ctak();
		//MPI_Barrier(MPI_COMM_WORLD);
		//pass_C();
		MPI_Barrier(MPI_COMM_WORLD);

		for (i = 0; i <= nx + 1; i++)
			for (j = 0; j <= ny + 1; j++)
				for (k = 0; k <= nz + 1; k++)
					C[i][j][k][0] = Cn[i][j][k][0];

		for (i = 0; i <= nx + 1; i++)
			for (j = 0; j <= ny + 1; j++)
				for (k = 0; k <= nz + 1; k++)
					C[i][j][k][1] = Cn[i][j][k][1];

		for (i = 0; i <= nx + 1; i++)
			for (j = 0; j <= ny + 1; j++)
				for (k = 0; k <= nz + 1; k++)
					C[i][j][k][2] = Cn[i][j][k][2];

		if (rank == (int)size / 2)
		{
			inf = fopen("inf_C.txt", "w");
			fprintf(inf, "n = %5d tt = %11.9f C05 = %11.9f\n", n, tt, C[(int)(nx / 2)][(int)(ny / 2)][(int)(nz / 2)][0]);
			//printf("%5d %9.6f %9.6f\n",n,max,max_all);
			fclose(inf);
		}

	}
	output_C();
	//END DEFINE C=========================================================================================================
	end = MPI_Wtime();
	printf("That took %f seconds\n", end - start);
	MPI_Finalize();
	//
	return 1;
}
//FUNCTIONS-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
void pass_C(void)
{
	int i, j, k, f;
	for (i = 0; i <= nx + 1; i++)
		for (k = 1; k <= nz; k++)
		{
			Cn[i][0][k][0] = Cn[i][2][k][0];
			Cn[i][ny + 1][k][0] = Cn[i][ny - 1][k][0];

			Cn[i][0][k][1] = Cn[i][2][k][1];
			Cn[i][ny + 1][k][1] = Cn[i][ny - 1][k][1];

			Cn[i][0][k][2] = Cn[i][2][k][2];
			Cn[i][ny + 1][k][2] = Cn[i][ny - 1][k][2];
		}
	for (i = 0; i <= nx + 1; i++)
		for (j = 0; j <= ny + 1; j++)
		{
			Cn[i][j][0][0] = Cn[i][j][1][0];
			Cn[i][j][nz + 1][0] = Cn[i][j][nz][0];

			Cn[i][j][0][1] = Cn[i][j][1][1];
			Cn[i][j][nz + 1][1] = Cn[i][j][nz][1];

			Cn[i][j][0][2] = Cn[i][j][1][2];
			Cn[i][j][nz + 1][2] = Cn[i][j][nz][2];
		}
	// By x -----------------------------------------------------------------------
	if (p_x != np_x - 1 && p_x != 0)
	{
		if (p_x % 2 == 0)
		{
			MPI_Send(&Cn[nx][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&Cn[nx + 1][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Cn[1][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Cn[0][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
		}
		else
		{
			MPI_Recv(&Cn[0][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Cn[1][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Cn[nx + 1][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Cn[nx][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
		}
	}
	else if (p_x == 0)
	{
		for (j = 0; j <= ny + 1; j++)
			for (k = 1; k <= nz; k++)
			{
				Cn[0][j][k][0] = Cn[2][j][k][0];

				Cn[0][j][k][1] = Cn[2][j][k][1];

				Cn[0][j][k][2] = Cn[2][j][k][2];
			}
		if (p_x != np_x - 1)
		{
			MPI_Send(&Cn[nx][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&Cn[nx + 1][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
		}
		else
			for (j = 0; j <= ny + 1; j++)
				for (k = 1; k <= nz; k++)
				{
					Cn[nx + 1][j][k][0] = Cn[nx - 1][j][k][0];

					Cn[nx + 1][j][k][1] = Cn[nx - 1][j][k][1];

					Cn[nx + 1][j][k][2] = Cn[nx - 1][j][k][2];
				}
	}
	else if (p_x == np_x - 1)
	{
		for (j = 0; j <= ny + 1; j++)
			for (k = 1; k <= nz; k++)
			{
				Cn[nx + 1][j][k][0] = Cn[nx - 1][j][k][0];

				Cn[nx + 1][j][k][1] = Cn[nx - 1][j][k][1];

				Cn[nx + 1][j][k][2] = Cn[nx - 1][j][k][2];
			}
		if (p_x % 2 == 0)
		{
			MPI_Send(&Cn[1][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Cn[0][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
		}
		else
		{
			MPI_Recv(&Cn[0][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Cn[1][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
		}
	}

}
void obtain_Crmp(void)
{
	int i, j, k;
	double Sc = 0.;
	double Dx_p, Dx_m, Dy_p, Dy_m, Dz_p, Dz_m;
	for (i = 1; i <= nx; i++)
		for (j = 1; j <= ny; j++)
			for (k = 1; k <= nz; k++)
			{
				Dx_p = 0.5*(D[i + 1][j][k][0] + D[i][j][k][0]);
				Dx_m = 0.5*(D[i - 1][j][k][0] + D[i][j][k][0]);
				Dy_p = 0.5*(D[i][j + 1][k][1] + D[i][j][k][1]);
				Dy_m = 0.5*(D[i][j - 1][k][1] + D[i][j][k][1]);
				Dz_p = 0.5*(D[i][j][k + 1][2] + D[i][j][k][2]);
				Dz_m = 0.5*(D[i][j][k - 1][2] + D[i][j][k][2]);

				//Reagent=========================
				Cn[i][j][k][0] = C[i][j][k][0] + tt*axc*(axc*Dx_p*(C[i + 1][j][k][0] - C[i][j][k][0]) - axc*Dx_m*(C[i][j][k][0] - C[i - 1][j][k][0])) / teta;
				Cn[i][j][k][0] = Cn[i][j][k][0] + tt*ayc*(ayc*Dy_p*(C[i][j + 1][k][0] - C[i][j][k][0]) - ayc*Dy_m*(C[i][j][k][0] - C[i][j - 1][k][0])) / teta;
				Cn[i][j][k][0] = Cn[i][j][k][0] + tt*azc*(azc*Dz_p*(C[i][j][k + 1][0] - C[i][j][k][0]) - azc*Dz_m*(C[i][j][k][0] - C[i][j][k - 1][0])) / teta;

				//Cn[i][j][k][0] = Cn[i][j][k][0] - tt*axc*0.5*U[i][j][k]*(C[i+1][j][k][0]-C[i-1][j][k][0])/teta;
				//Cn[i][j][k][0] = Cn[i][j][k][0] - tt*ayc*0.5*V[i][j][k]*(C[i][j+1][k][0]-C[i][j-1][k][0])/teta;
				//Cn[i][j][k][0] = Cn[i][j][k][0] - tt*azc*0.5*W[i][j][k]*(C[i][j][k+1][0]-C[i][j][k-1][0])/teta;
				Cn[i][j][k][0] = Cn[i][j][k][0] - tt*(0.5*axc*(U[i][j][k] + fabs(U[i][j][k]))*(C[i][j][k][0] - C[i - 1][j][k][0]) + 0.5*axc*(U[i][j][k] - fabs(U[i][j][k]))*(C[i + 1][j][k][0] - C[i][j][k][0])) / teta;
				Cn[i][j][k][0] = Cn[i][j][k][0] - tt*(0.5*ayc*(V[i][j][k] + fabs(V[i][j][k]))*(C[i][j][k][0] - C[i][j - 1][k][0]) + 0.5*ayc*(V[i][j][k] - fabs(V[i][j][k]))*(C[i][j + 1][k][0] - C[i][j][k][0])) / teta;
				Cn[i][j][k][0] = Cn[i][j][k][0] - tt*(0.5*azc*(W[i][j][k] + fabs(W[i][j][k]))*(C[i][j][k][0] - C[i][j][k - 1][0]) + 0.5*azc*(W[i][j][k] - fabs(W[i][j][k]))*(C[i][j][k + 1][0] - C[i][j][k][0])) / teta;

				Cn[i][j][k][0] = Cn[i][j][k][0] - tt*v1*beta*C[i][j][k][1] * C[i][j][k][0] + tt*Qr_z[i][j][k] * Cr0 / teta - tt*Qr_o[i][j][k] * C[i][j][k][0] / teta;
				//================================
				//Mineral=========================
				Cn[i][j][k][1] = C[i][j][k][1] * (1. - tt*beta*teta*C[i][j][k][0]);
				//================================
				//Poleznoe========================
				Cn[i][j][k][2] = C[i][j][k][2] + tt*axc*(axc*Dx_p*(C[i + 1][j][k][2] - C[i][j][k][2]) - axc*Dx_m*(C[i][j][k][2] - C[i - 1][j][k][2])) / teta;
				Cn[i][j][k][2] = Cn[i][j][k][2] + tt*ayc*(ayc*Dy_p*(C[i][j + 1][k][2] - C[i][j][k][2]) - ayc*Dy_m*(C[i][j][k][2] - C[i][j - 1][k][2])) / teta;
				Cn[i][j][k][2] = Cn[i][j][k][2] + tt*azc*(azc*Dz_p*(C[i][j][k + 1][2] - C[i][j][k][2]) - azc*Dz_m*(C[i][j][k][2] - C[i][j][k - 1][2])) / teta;

				//Cn[i][j][k][0] = Cn[i][j][k][0] - tt*axc*0.5*U[i][j][k]*(C[i+1][j][k][0]-C[i-1][j][k][0])/teta;
				//Cn[i][j][k][0] = Cn[i][j][k][0] - tt*ayc*0.5*V[i][j][k]*(C[i][j+1][k][0]-C[i][j-1][k][0])/teta;
				//Cn[i][j][k][0] = Cn[i][j][k][0] - tt*azc*0.5*W[i][j][k]*(C[i][j][k+1][0]-C[i][j][k-1][0])/teta;
				Cn[i][j][k][2] = Cn[i][j][k][2] - tt*(0.5*axc*(U[i][j][k] + fabs(U[i][j][k]))*(C[i][j][k][2] - C[i - 1][j][k][2]) + 0.5*axc*(U[i][j][k] - fabs(U[i][j][k]))*(C[i + 1][j][k][2] - C[i][j][k][2])) / teta;
				Cn[i][j][k][2] = Cn[i][j][k][2] - tt*(0.5*ayc*(V[i][j][k] + fabs(V[i][j][k]))*(C[i][j][k][2] - C[i][j - 1][k][2]) + 0.5*ayc*(V[i][j][k] - fabs(V[i][j][k]))*(C[i][j + 1][k][2] - C[i][j][k][2])) / teta;
				Cn[i][j][k][2] = Cn[i][j][k][2] - tt*(0.5*azc*(W[i][j][k] + fabs(W[i][j][k]))*(C[i][j][k][2] - C[i][j][k - 1][2]) + 0.5*azc*(W[i][j][k] - fabs(W[i][j][k]))*(C[i][j][k + 1][2] - C[i][j][k][2])) / teta;

				Cn[i][j][k][2] = Cn[i][j][k][2] + tt*v2*beta*C[i][j][k][1] * C[i][j][k][0] - tt*Qp[i][j][k] * C[i][j][k][2] / teta;

			}
}
void pass_D(void)
{
	int i, j, k;
	// By x -----------------------------------------------------------------------
	if (p_x != np_x - 1 && p_x != 0)
	{
		if (p_x % 2 == 0)
		{
			MPI_Send(&D[nx][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&D[nx + 1][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&D[1][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&D[0][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
		}
		else
		{
			MPI_Recv(&D[0][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&D[1][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&D[nx + 1][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&D[nx][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
		}
	}
	else if (p_x == 0)
	{
		if (p_x != np_x - 1)
		{
			MPI_Send(&D[nx][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&D[nx + 1][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
		}

	}
	else if (p_x == np_x - 1)
	{
		if (p_x % 2 == 0)
		{
			MPI_Send(&D[1][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&D[0][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
		}
		else
		{
			MPI_Recv(&D[0][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&D[1][0][0][0], (nz + 2)*(ny + 2) * 3, MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
		}
	}

}
void Define_D(void)
{
	int i, j, k, m;
	double abV, abV1;
	for (i = 0; i <= nx + 1; i++)
		for (j = 0; j <= ny + 1; j++)
			for (k = 0; k <= nz + 1; k++)
			{
				if (U[i][j][k] != 0 && V[i][j][k] != 0 && W[i][j][k] != 0)
				{
					abV = 1. / (sqrt(U[i][j][k] * U[i][j][k] + V[i][j][k] * V[i][j][k] + W[i][j][k] * W[i][j][k]));
					D[i][j][k][0] = alpha*abV*U[i][j][k] * U[i][j][k] + alpha*abV*V[i][j][k] * V[i][j][k] + alpha*abV*W[i][j][k] * W[i][j][k] + teta*_D;
					D[i][j][k][1] = alpha*abV*V[i][j][k] * V[i][j][k] + alpha*abV*U[i][j][k] * U[i][j][k] + alpha*abV*W[i][j][k] * W[i][j][k] + teta*_D;
					D[i][j][k][2] = alpha*abV*W[i][j][k] * W[i][j][k] + alpha*abV*U[i][j][k] * U[i][j][k] + alpha*abV*V[i][j][k] * V[i][j][k] + teta*_D;
				}
				else
				{
					D[i][j][k][0] = teta*_D;
					D[i][j][k][1] = teta*_D;
					D[i][j][k][2] = teta*_D;
				}
			}
	//Correction
}
void Define_Vel(void)
{
	int i, j, k;
	//Initial======================
	for (i = 0; i <= nx + 1; i++)
		for (j = 0; j <= ny + 1; j++)
			for (k = 0; k <= nz + 1; k++)
			{
				U[i][j][k] = 0;
				V[i][j][k] = 0;
				W[i][j][k] = 0;
			}

	//=============================
	for (i = 1; i <= nx; i++)
		for (j = 1; j <= ny; j++)
			for (k = 1; k <= nz; k++)
			{
				U[i][j][k] = -0.5*K*(H[i + 1][j][k] - H[i - 1][j][k]) / dx;
				V[i][j][k] = -0.5*K*(H[i][j + 1][k] - H[i][j - 1][k]) / dy;
				W[i][j][k] = -0.5*K*(H[i][j][k + 1] - H[i][j][k - 1]) / dz;
			}

	//Boundary==============
	//simmetriya
	for (i = 1; i <= nx; i++)
		for (k = 1; k <= nz; k++)
		{
			U[i][0][k] = U[i][2][k];
			U[i][ny + 1][k] = U[i][ny - 1][k];
			V[i][0][k] = -V[i][2][k];
			V[i][ny + 1][k] = V[i][ny - 1][k];
			W[i][0][k] = W[i][2][k];
			W[i][ny + 1][k] = W[i][ny - 1][k];
		}
	//prelipaniya
	for (i = 1; i <= nx; i++)
		for (j = 1; j <= ny; j++)
		{
			U[i][j][0] = 0.;
			U[i][j][nz + 1] = 0.;
			V[i][j][0] = 0.;
			V[i][j][nz + 1] = 0.;
			W[i][j][0] = 0.;
			W[i][j][nz + 1] = 0.;
		}
	//simmerty
	for (j = 1; j <= ny; j++)
		for (k = 1; k <= nz; k++)
		{
			U[0][j][k] = -U[2][j][k];
			U[nx + 1][j][k] = -U[nx - 1][j][k];
			V[0][j][k] = V[2][j][k];
			V[nx + 1][j][k] = V[nx - 1][j][k];
			W[0][j][k] = W[2][j][k];
			W[nx + 1][j][k] = W[nx - 1][j][k];
		}
	//Edges
	for (i = 0; i <= nx + 1; i++)
	{
		U[i][0][0] = 0.;
		V[i][0][0] = 0.;
		W[i][0][0] = 0.;

		U[i][ny + 1][0] = 0.;
		V[i][ny + 1][0] = 0.;
		W[i][ny + 1][0] = 0.;

		U[i][ny + 1][nz + 1] = 0.;
		V[i][ny + 1][nz + 1] = 0.;
		W[i][ny + 1][nz + 1] = 0.;

		U[i][0][nz + 1] = 0.;
		V[i][0][nz + 1] = 0.;
		W[i][0][nz + 1] = 0.;
	}
	for (j = 0; j <= ny + 1; j++)
	{
		U[0][j][0] = 0.;
		V[0][j][0] = 0.;
		W[0][j][0] = 0.;

		U[0][j][nz + 1] = 0.;
		V[0][j][nz + 1] = 0.;
		W[0][j][nz + 1] = 0.;

		U[nx + 1][j][nz + 1] = 0.;
		V[nx + 1][j][nz + 1] = 0.;
		W[nx + 1][j][nz + 1] = 0.;

		U[nx + 1][j][0] = 0.;
		V[nx + 1][j][0] = 0.;
		W[nx + 1][j][0] = 0.;
	}
	for (k = 1; k <= nz; k++)
	{
		U[0][0][k] = -U[2][2][k];
		V[0][0][k] = -V[2][2][k];
		W[0][0][k] = W[2][2][k];

		U[nx + 1][0][k] = -U[nx - 1][2][k];
		V[nx + 1][0][k] = -V[nx - 1][2][k];
		W[nx + 1][0][k] = W[nx - 1][2][k];

		U[nx + 1][ny + 1][k] = -U[nx - 1][ny - 1][k];
		V[nx + 1][ny + 1][k] = -V[nx - 1][ny - 1][k];
		W[nx + 1][ny + 1][k] = W[nx - 1][ny - 1][k];

		U[0][ny + 1][k] = -U[2][ny - 1][k];
		V[0][ny + 1][k] = -V[2][ny - 1][k];
		W[0][ny + 1][k] = W[2][ny - 1][k];
	}
	for (i = 1; i <= nx; i++)
		for (j = 1; j <= ny; j++)
			for (k = 1; k <= nz; k++)
			{
				if (U[i][j][k] > maxV)
					maxV = U[i][j][k];
				if (V[i][j][k] > maxV)
					maxV = V[i][j][k];
				if (W[i][j][k] > maxV)
					maxV = W[i][j][k];
			}
}
void pass_H(void)
{
	int i, j, k;
	for (i = 0; i <= nx + 1; i++)
		for (k = 1; k <= nz; k++)
		{
			Hn[i][0][k] = Hn[i][2][k];
			Hn[i][ny + 1][k] = Hn[i][ny - 1][k];
		}
	for (i = 0; i <= nx + 1; i++)
		for (j = 0; j <= ny + 1; j++)
		{
			Hn[i][j][0] = Hn[i][j][1];
			Hn[i][j][nz + 1] = Hn[i][j][nz];
		}
	// By x -----------------------------------------------------------------------
	if (p_x != np_x - 1 && p_x != 0)
	{
		if (p_x % 2 == 0)
		{
			MPI_Send(&Hn[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&Hn[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Hn[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Hn[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
		}
		else
		{
			MPI_Recv(&Hn[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Hn[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Hn[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Hn[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
		}
	}
	else if (p_x == 0)
	{
		for (j = 0; j <= ny + 1; j++)
			for (k = 1; k <= nz; k++)
			{
				Hn[0][j][k] = Hn[2][j][k];
			}
		if (p_x != np_x - 1)
		{
			MPI_Send(&Hn[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&Hn[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
		}
		else
			for (j = 0; j <= ny + 1; j++)
				for (k = 1; k <= nz; k++)
					Hn[nx + 1][j][k] = Hn[nx - 1][j][k];
	}
	else if (p_x == np_x - 1)
	{
		for (j = 0; j <= ny + 1; j++)
			for (k = 1; k <= nz; k++)
				Hn[nx + 1][j][k] = Hn[nx - 1][j][k];
		if (p_x % 2 == 0)
		{
			MPI_Send(&Hn[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Hn[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
		}
		else
		{
			MPI_Recv(&Hn[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Hn[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
		}
	}

}
int Define_Q(void)
{
	int i, j, k, rank1, rank2, rank3, rank4, rank5;
	double r, r1;
	int w1y, w2y, w3y, w1x, w2x, w3x, w4x, w5x;
	r = 13.;
	r1 = sqrt(3. / 4.*r*r);
	//r = r1*2.;
	//w1y = int((Ly/2.-r1)/dy);
	w2y = int((Ly / 2.) / dy);
	//w3y = int((Ly/2.+r1)/dy);
	//
	//w1x = int((Lx/2.-r   )/dx) - rank*nx;
	//w2x = int((Lx/2.-r/2.)/dx) - rank*nx;
	w3x = int((Lx / 2.) / dx) - rank*nx;
	//w4x = int((Lx/2.+r/2.)/dx) - rank*nx;
	//w5x = int((Lx/2.+r   )/dx) - rank*nx;

	//printf("\n r=%d, %d %d %d %d %d",w1x,w2x,w3x,w4x,w5x);

	for (i = 0; i <= nx + 1; i++)
		for (j = 0; j <= ny + 1; j++)
			for (k = 0; k <= nz + 1; k++)
			{
				Q[i][j][k] = 0.;
				Qp[i][j][k] = 0.;
				Qr_z[i][j][k] = 0.;
				Qr_o[i][j][k] = 0.;
			}

	/*if(0<=w1x && w1x<=nx+1)
	for(k=0;k<=(nz+1)/2;k++)
	{
	Q[w1x][w2y][k] = -qs/6.;
	Qr_z[w1x][w2y][k] = 0.;
	}*/

	/*if(0<=w2x && w2x<=nx+1)
	for(k=0;k<=(nz+1)/2;k++)
	{
	Q[w2x][w1y][k] = -qs/6.;
	Q[w2x][w3y][k] = -qs/6.;
	Qr_z[w2x][w1y][k] = 0.;
	Qr_z[w2x][w3y][k] = 0.;
	}*/
	if (rank == np_x - 1)
		for (k = 3; k <= (nz - 2); k++)
		{
			Q[nx][w2y][k] = qs / 2.;
			Qr_z[nx][w2y][k] = qs / 2.;
		}
	if (rank == 0)
		for (k = 3; k <= (nz - 2); k++)
		{
			Q[1][w2y][k] = qs / 2.;
			Qr_z[1][w2y][k] = qs / 2.;
		}

	if (0 <= w3x && w3x <= nx + 1)
		for (k = 3; k <= (nz - 2); k++)
		{
			Q[w3x][w2y][k] = -qs;
			Qr_o[w3x][w2y][k] = qs;
			Qp[w3x][w2y][k] = qs;
		}
	/*
	if(0<=w4x && w4x<=nx+1)
	for(k=0;k<=(nz+1)/2;k++)
	{
	Q[w4x][w1y][k] = -qs/6.;
	Q[w4x][w3y][k] = -qs/6.;
	}

	if(0<=w5x && w5x<=nx+1)
	for(k=0;k<=(nz+1)/2;k++)
	{
	Q[w5x][w2y][k] = -qs/6.;
	}*/
	return 1;
}
int obtain_Hzhup(void)
{
	int i, j, k;

	for (i = 1; i <= nx; i++)
	{
		for (j = 1; j <= ny; j++)
		{
			for (k = 1; k <= nz; k++)
			{
				if ((i + j + k + n) % 2 == 0)
				{
					Hn[i][j][k] = H[i][j][k]
						+ dt * ax * (H[i + 1][j][k] + H[i - 1][j][k] - 2 * H[i][j][k])
						+ dt * ay * (H[i][j + 1][k] + H[i][j - 1][k] - 2 * H[i][j][k])
						+ dt * az * (H[i][j][k + 1] + H[i][j][k - 1] - 2 * H[i][j][k])
						+ dt * Q[i][j][k] / K;
				}
			}
		}
	}

	return 1;
}
int obtain_Htak(void)
{
	int i, j, k;
	for (i = 1; i <= nx; i++)
	{
		for (j = 1; j <= ny; j++)
		{
			for (k = 1; k <= nz; k++)
			{
				if ((i + j + k + n) % 2 != 0)
				{
					Hn[i][j][k] = (H[i][j][k]
						+ dt * ax * (Hn[i + 1][j][k] + Hn[i - 1][j][k])
						+ dt * ay * (Hn[i][j + 1][k] + Hn[i][j - 1][k])
						+ dt * az * (Hn[i][j][k + 1] + Hn[i][j][k - 1])
						+ dt * Q[i][j][k] / K) / (1. + 2.*dt*(ax + ay + az));
				}
			}
		}
	}
	max = 0.;
	for (i = 1; i <= nx; i++)
		for (j = 1; j <= ny; j++)
			for (k = 1; k <= nz; k++)
			{
				if (fabs(H[i][j][k] - Hn[i][j][k]) > max)
					max = fabs(H[i][j][k] - Hn[i][j][k]);
			}

	// ...
	return 1;
}
void output_H(void)
{
	int i, j, k, idf;
	char a[] = "0123456789";
	char num[] = "out//H_grid_0000.dat";
	//
	//num[15] = a[rank%10];
	//num[14] = a[(rank/10)%10];
	//num[13] = a[(rank/100)%10];
	//num[12] = a[(rank/1000)%10];
	//
	FILE *f;
	double x, y, z; // real domain points
	if (rank != 0)
	{
		MPI_Send(&H[0][0][0], (nx + 2)*(nz + 2)*(ny + 2), MPI_DOUBLE, 0, tag_next_x, MPI_COMM_WORLD);
	}
	else
	{
		f = fopen("Press.tec", "w");
		fprintf(f, "TITLE = \"PRESS:3D\"\n");
		fprintf(f, "VARIABLES = \"X\",\"Y\",\"Z\",\"H\"\n");
		fprintf(f, "ZONE I=%d, J=%d, K=%d, F=POINT\n", Nz, Ny, Nx);
		for (i = 1; i <= nx; i++)
			for (j = 1; j <= ny; j++)
				for (k = 1; k <= nz; k++)
				{
					x = (0.*nx + i)*dx;
					y = j*dy;
					z = k*dz;
					fprintf(f, "%.5f %.5f %.5f %.7f\n", x, y, z, H[i][j][k]);
				}
		for (idf = 1; idf<size; idf++)
		{
			MPI_Recv(&buf[0][0][0], (nx + 2)*(nz + 2)*(ny + 2), MPI_DOUBLE, idf, tag_next_x, MPI_COMM_WORLD, &status);
			for (i = 1; i <= nx; i++)
				for (j = 1; j <= ny; j++)
					for (k = 1; k <= nz; k++)
					{
						x = (idf*nx + i)*dx;
						y = j*dy;
						z = k*dz;
						fprintf(f, "%.5f %.5f %.5f %.7f\n", x, y, z, buf[i][j][k]);
					}
		}
		fclose(f);
	}
}
void output_C(void)
{
	int i, j, k, idf;
	char a[] = "0123456789";
	char num[] = "dat//C_grid_0000.dat";
	//
	//num[15] = a[rank%10];
	//num[14] = a[(rank/10)%10];
	//num[13] = a[(rank/100)%10];
	//num[12] = a[(rank/1000)%10];
	//
	FILE *f;
	double x, y, z; // real domain points
					//CP==============================================================
	if (rank != 0)
	{
		MPI_Send(&C[0][0][0][0], (nx + 2)*(nz + 2)*(ny + 2) * 3, MPI_DOUBLE, 0, tag_next_x, MPI_COMM_WORLD);
	}
	else
	{
		f = fopen("CR.tec", "w");
		fprintf(f, "TITLE = \"CR:3D\"\n");
		fprintf(f, "VARIABLES = \"X\",\"Y\",\"Z\",\"Cr\"\n");
		fprintf(f, "ZONE I=%d, J=%d, K=%d, F=POINT\n", Nz, Ny, Nx);
		for (i = 1; i <= nx; i++)
			for (j = 1; j <= ny; j++)
				for (k = 1; k <= nz; k++)
				{
					x = (0.*nx + i)*dx;
					y = j*dy;
					z = k*dz;
					fprintf(f, "%.5f %.5f %.5f %.7f\n", x, y, z, C[i][j][k][0]);
				}
		for (idf = 1; idf<size; idf++)
		{
			MPI_Recv(&buf1[0][0][0][0], (nx + 2)*(nz + 2)*(ny + 2) * 3, MPI_DOUBLE, idf, tag_next_x, MPI_COMM_WORLD, &status);
			for (i = 1; i <= nx; i++)
				for (j = 1; j <= ny; j++)
					for (k = 1; k <= nz; k++)
					{
						x = (idf*nx + i)*dx;
						y = j*dy;
						z = k*dz;
						fprintf(f, "%.5f %.5f %.5f %.7f\n", x, y, z, buf1[i][j][k][0]);
					}
		}
		fclose(f);
	}
	//CM============================================
	if (rank != 0)
	{
		MPI_Send(&C[0][0][0][0], (nx + 2)*(nz + 2)*(ny + 2) * 3, MPI_DOUBLE, 0, tag_next_x, MPI_COMM_WORLD);
	}
	else
	{
		f = fopen("CM.tec", "w");
		fprintf(f, "TITLE = \"CM:3D\"\n");
		fprintf(f, "VARIABLES = \"X\",\"Y\",\"Z\",\"Cm\"\n");
		fprintf(f, "ZONE I=%d, J=%d, K=%d, F=POINT\n", Nz, Ny, Nx);
		for (i = 1; i <= nx; i++)
			for (j = 1; j <= ny; j++)
				for (k = 1; k <= nz; k++)
				{
					x = (0.*nx + i)*dx;
					y = j*dy;
					z = k*dz;
					fprintf(f, "%.5f %.5f %.5f %.7f\n", x, y, z, C[i][j][k][1]);
				}
		for (idf = 1; idf<size; idf++)
		{
			MPI_Recv(&buf1[0][0][0][0], (nx + 2)*(nz + 2)*(ny + 2) * 3, MPI_DOUBLE, idf, tag_next_x, MPI_COMM_WORLD, &status);
			for (i = 1; i <= nx; i++)
				for (j = 1; j <= ny; j++)
					for (k = 1; k <= nz; k++)
					{
						x = (idf*nx + i)*dx;
						y = j*dy;
						z = k*dz;
						fprintf(f, "%.5f %.5f %.5f %.7f\n", x, y, z, buf1[i][j][k][1]);
					}
		}
		fclose(f);
	}
	//CP===================================================
	if (rank != 0)
	{
		MPI_Send(&C[0][0][0][0], (nx + 2)*(nz + 2)*(ny + 2) * 3, MPI_DOUBLE, 0, tag_next_x, MPI_COMM_WORLD);
	}
	else
	{
		f = fopen("CP.tec", "w");
		fprintf(f, "TITLE = \"CP:3D\"\n");
		fprintf(f, "VARIABLES = \"X\",\"Y\",\"Z\",\"Cp\"\n");
		fprintf(f, "ZONE I=%d, J=%d, K=%d, F=POINT\n", Nz, Ny, Nx);
		for (i = 1; i <= nx; i++)
			for (j = 1; j <= ny; j++)
				for (k = 1; k <= nz; k++)
				{
					x = (0.*nx + i)*dx;
					y = j*dy;
					z = k*dz;
					fprintf(f, "%.5f %.5f %.5f %.7f\n", x, y, z, C[i][j][k][2]);
				}
		for (idf = 1; idf<size; idf++)
		{
			MPI_Recv(&buf1[0][0][0][0], (nx + 2)*(nz + 2)*(ny + 2) * 3, MPI_DOUBLE, idf, tag_next_x, MPI_COMM_WORLD, &status);
			for (i = 1; i <= nx; i++)
				for (j = 1; j <= ny; j++)
					for (k = 1; k <= nz; k++)
					{
						x = (idf*nx + i)*dx;
						y = j*dy;
						z = k*dz;
						fprintf(f, "%.5f %.5f %.5f %.7f\n", x, y, z, buf1[i][j][k][2]);
					}
		}
		fclose(f);
	}
}
void output_int(void)
{
	int i, j, k, idf;
	char a[] = "0123456789";
	char num[] = "int//Q_grid_0000.dat";
	//
	//num[15] = a[rank%10];
	//num[14] = a[(rank/10)%10];
	//num[13] = a[(rank/100)%10];
	//num[12] = a[(rank/1000)%10];
	//
	FILE *f;
	double x, y, z; // real domain points
	if (rank != 0)
	{
		MPI_Send(&Q[0][0][0], (nx + 2)*(nz + 2)*(ny + 2), MPI_DOUBLE, 0, tag_next_x, MPI_COMM_WORLD);
	}
	else
	{
		f = fopen("Wells.dat", "w");
		fprintf(f, "VARIABLES = \"X\",\"Y\",\"Z\",\"Q\"\n");
		fprintf(f, "ZONE I=%d, J=%d, K=%d F=POINT\n", (Nz + 1), (Ny + 1), (Nx + 1));
		for (i = 0; i <= nx; i++)
			for (j = 0; j <= ny; j++)
				for (k = 0; k <= nz; k++)
				{
					x = (0.*nx + i)*dx;
					y = j*dy;
					z = k*dz;
					fprintf(f, "%.5f %.5f %.5f %.7f\n", x, y, z, Q[i][j][k]);
				}
		for (idf = 1; idf<size; idf++)
		{
			MPI_Recv(&buf[0][0][0], (nx + 2)*(nz + 2)*(ny + 2), MPI_DOUBLE, idf, tag_next_x, MPI_COMM_WORLD, &status);
			for (i = 0; i <= nx; i++)
				for (j = 0; j <= ny; j++)
					for (k = 0; k <= nz; k++)
					{
						x = (idf*nx + i)*dx;
						y = j*dy;
						z = k*dz;
						fprintf(f, "%.5f %.5f %.5f %.7f\n", x, y, z, buf[i][j][k]);
					}
		}
		fclose(f);
	}
}
void pass_Q(void)
{

	// By x -----------------------------------------------------------------------
	if (p_x != np_x - 1 && p_x != 0)
	{
		if (p_x % 2 == 0)
		{
			MPI_Send(&Q[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&Q[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Q[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Q[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);

			MPI_Send(&Qp[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&Qp[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Qp[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Qp[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);

			MPI_Send(&Qr_z[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&Qr_z[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Qr_z[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Qr_z[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);

			MPI_Send(&Qr_o[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&Qr_o[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Qr_o[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Qr_o[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
		}
		else
		{
			MPI_Recv(&Q[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Q[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Q[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Q[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);

			MPI_Recv(&Qp[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Qp[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Qp[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Qp[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);

			MPI_Recv(&Qr_z[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Qr_z[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Qr_z[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Qr_z[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);

			MPI_Recv(&Qr_o[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Qr_o[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Qr_o[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Qr_o[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
		}
	}
	else if (p_x == 0)
	{
		if (p_x != np_x - 1)
		{
			MPI_Send(&Q[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&Q[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);

			MPI_Send(&Qp[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&Qp[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);

			MPI_Send(&Qr_z[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&Qr_z[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);

			MPI_Send(&Qr_o[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&Qr_o[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
		}
	}
	else if (p_x == np_x - 1)
	{
		if (p_x % 2 == 0)
		{
			MPI_Send(&Q[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Q[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);

			MPI_Send(&Qp[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Qp[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);

			MPI_Send(&Qr_z[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Qr_z[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);

			MPI_Send(&Qr_o[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&Qr_o[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
		}
		else
		{
			MPI_Recv(&Q[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Q[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);

			MPI_Recv(&Qp[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Qp[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);

			MPI_Recv(&Qr_z[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Qr_z[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);

			MPI_Recv(&Qr_o[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&Qr_o[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
		}
	}

}
void pass_U(void)
{

	// By x -----------------------------------------------------------------------
	if (p_x != np_x - 1 && p_x != 0)
	{
		if (p_x % 2 == 0)
		{
			MPI_Send(&U[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&U[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&U[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&U[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);

			MPI_Send(&V[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&V[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&V[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&V[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);

			MPI_Send(&W[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&W[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&W[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&W[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
		}
		else
		{
			MPI_Recv(&U[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&U[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&U[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&U[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);

			MPI_Recv(&V[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&V[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&V[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&V[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);

			MPI_Recv(&W[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&W[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&W[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
			MPI_Send(&W[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
		}
	}
	else if (p_x == 0)
	{
		if (p_x != np_x - 1)
		{
			MPI_Send(&U[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&U[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);

			MPI_Send(&V[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&V[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);

			MPI_Send(&W[nx][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD);
			MPI_Recv(&W[nx + 1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, next_x, tag_next_x, MPI_COMM_WORLD, &status);
		}
	}
	else if (p_x == np_x - 1)
	{
		if (p_x % 2 == 0)
		{
			MPI_Send(&U[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&U[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);

			MPI_Send(&V[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&V[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);

			MPI_Send(&W[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
			MPI_Recv(&W[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
		}
		else
		{
			MPI_Recv(&U[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&U[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);

			MPI_Recv(&V[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&V[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);

			MPI_Recv(&W[0][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD, &status);
			MPI_Send(&W[1][0][0], (nz + 2)*(ny + 2), MPI_DOUBLE, prev_x, tag_prev_x, MPI_COMM_WORLD);
		}
	}

}
//+-+-+-++--+--++-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-++-+-+-+-+-+-+-+--++-+-=-=-=-=-==-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
