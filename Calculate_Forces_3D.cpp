// Multi-threaded calculation functions
#pragma once

#include <chrono>
#include <random>
#include <cmath>
#include <fstream>
#include <stdexcept>

#include <immintrin.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <vector>

#include "General_Closed_Particle_Simulation.h"

using namespace std;

extern Cell_Structure Cell_Grid;

extern Box_Structure Box;

extern Particle_Structure* Particles;

// simulation variables
extern double Particle_Diameter;

extern double End_Time;
extern double Delta_t;
extern double Epsilon_WCA;

extern int Dimension;
extern int Number_Particles;
extern int Num_Threads;

extern double Region_Size[3];

extern double Sigma_Square;
extern double Sigma_Test_WCA;
extern double Sigma_Test_WCA_Colloid_Particle;
extern double Sigma_Test_LJ;
extern double Epsilon_WCA_x_24;
extern double Delta_t_square;

extern int Neighbor_Count;

extern mutex m;

extern mutex ctrl[4];

extern mutex a;

extern mutex cmp;

extern volatile atomic_bool* Ready;

extern int* Update_Divisions_Start;
extern int* Update_Divisions_End;

extern int* Update_Particles_Start;
extern int* Update_Particles_End;

extern volatile atomic_bool Finished;

extern volatile atomic_int Update_Completed;

extern volatile atomic_int Force_Completed;

extern volatile atomic_int Loop_Completed;

extern int Active_Threads;

extern int Thread_Complete_Test;

extern condition_variable Active;


void Calculate_Forces_Looped_Threaded_3D(int Thread_Num)
{
	a.lock();

	Active_Threads += 1;

	if (Active_Threads == 4)
	{
		Active.notify_all();
	}

	a.unlock();

	//cout << "start " << Thread_Num << endl;

	while (true)
	{
		while (!Ready[Thread_Num]) {}

		Ready[Thread_Num] = false;

		if (Finished) { return; }

		Calculate_Forces_Threaded_3D(Thread_Num);

		// wait until all forces calculated, then update time tick
		Force_Completed += (Thread_Num + 1);

		while (Force_Completed != 10) {}

		for (int i = Update_Particles_Start[Thread_Num]; i < Update_Particles_End[Thread_Num]; i++)
		{
			Particles[i].Position_X += Particles[i].Velocity_X * Delta_t;
			Particles[i].Position_Y += Particles[i].Velocity_Y * Delta_t;
			Particles[i].Position_Z += Particles[i].Velocity_Z * Delta_t;

			Particles[i].Velocity_X += Particles[i].Force_X * Delta_t;
			Particles[i].Velocity_Y += Particles[i].Force_Y * Delta_t;
			Particles[i].Velocity_Z += Particles[i].Force_Z * Delta_t;

			// Debug_Info(0, i);

			Particles[i].Force_X = 0;
			Particles[i].Force_Y = 0;
			Particles[i].Force_Z = 0;
		}

		Loop_Completed += (Thread_Num + 1);
	}
}

void Build_Cell_List_3D()
{
	double X_Divisor, Y_Divisor, Z_Divisor;

	int Local_Cell_X, Local_Cell_Y, Local_Cell_Z;

	int Final_Cell;

	int Div_X = Cell_Grid.Divisions[0];
	int Div_Y = Cell_Grid.Divisions[1];
	int Div_Z = Cell_Grid.Divisions[2];
	int Div_Tot = Cell_Grid.Total_Divisions;

	for (int i = 0; i < Cell_Grid.Total_Divisions; i++)
	{
		Cell_Grid.Particles_Per_Cell[i] = 0;
	}

	int Block_Division = Div_X * Div_Y;

	X_Divisor = Region_Size[0] / Div_X;
	Y_Divisor = Region_Size[1] / Div_Y;
	Z_Divisor = Region_Size[2] / Div_Z;

	for (int i = 0; i < Number_Particles; i++)
	{
		// assign each particle to a cell
		// keeping track of the counts per cell

		Local_Cell_X = (int)floor(Particles[i].Position_X / X_Divisor);
		Local_Cell_Y = (int)floor(Particles[i].Position_Y / Y_Divisor);
		Local_Cell_Z = (int)floor(Particles[i].Position_Z / Z_Divisor);

		Local_Cell_Y *= Div_X;

		Local_Cell_Z *= Block_Division;

		Final_Cell = Local_Cell_Z + Local_Cell_Y + Local_Cell_X;

		Cell_Grid.Particle_Cell_List[Final_Cell][Cell_Grid.Particles_Per_Cell[Final_Cell]] = i;

		Cell_Grid.Particles_Per_Cell[Final_Cell] += 1;
	}
}

void Calculate_Forces_3D()
{
	for (int i = 0; i < Cell_Grid.Total_Divisions; i++)
	{
		if (Cell_Grid.Particles_Per_Cell[i] != 0)
		{
			// note this only runs if more than 1 particle is in the cell
			for (int j = 0; j < Cell_Grid.Particles_Per_Cell[i] - 1; j++)
			{
				// calculate forces on particles in the current cell with each neighbor in this cell
				for (int k = j + 1; k < Cell_Grid.Particles_Per_Cell[i]; k++)
				{
					Calculate_Particle_Forces_3D(Cell_Grid.Particle_Cell_List[i][j], Cell_Grid.Particle_Cell_List[i][k]);
				}
			}

			for (int n = 0; n < Neighbor_Count; n++)
			{
				if (Cell_Grid.Neighbor_List[i][n] != -1)
				{
					if (Cell_Grid.Particles_Per_Cell[Cell_Grid.Neighbor_List[i][n]] != 0)
					{
						for (int j = 0; j < Cell_Grid.Particles_Per_Cell[i]; j++)
						{
							// calculate forces on particles in the current cell with each neighboring cell
							for (int k = 0; k < Cell_Grid.Particles_Per_Cell[Cell_Grid.Neighbor_List[i][n]]; k++)
							{
								Calculate_Particle_Forces_3D(Cell_Grid.Particle_Cell_List[i][j], Cell_Grid.Particle_Cell_List[Cell_Grid.Neighbor_List[i][n]][k]);
							}
						}
					}
				}
			}
		}
	}
}

void Calculate_Forces_Threaded_3D(int Thread_Num)
{
	int start = Cell_Grid.XY_Divisions * Cell_Grid.Z_Div_4 * Thread_Num;

	int end = start + Cell_Grid.XY_Divisions * Cell_Grid.Z_Div_4;

	for (int i = start; i < end; i++)
	{
		if (Cell_Grid.Particles_Per_Cell[i] != 0)
		{
			// note this only runs if more than 1 particle is in the cell
			for (int j = 0; j < Cell_Grid.Particles_Per_Cell[i] - 1; j++)
			{
				// calculate forces on particles in the current cell with each neighbor in this cell
				for (int k = j + 1; k < Cell_Grid.Particles_Per_Cell[i]; k++)
				{
					Calculate_Particle_Forces_3D(Cell_Grid.Particle_Cell_List[i][j], Cell_Grid.Particle_Cell_List[i][k]);
				}
			}

			for (int n = 0; n < Neighbor_Count; n++)
			{
				if (Cell_Grid.Neighbor_List[i][n] != -1)
				{
					if (Cell_Grid.Particles_Per_Cell[Cell_Grid.Neighbor_List[i][n]] != 0)
					{
						for (int j = 0; j < Cell_Grid.Particles_Per_Cell[i]; j++)
						{
							// calculate forces on particles in the current cell with each neighboring cell
							for (int k = 0; k < Cell_Grid.Particles_Per_Cell[Cell_Grid.Neighbor_List[i][n]]; k++)
							{
								Calculate_Particle_Forces_3D(Cell_Grid.Particle_Cell_List[i][j], Cell_Grid.Particle_Cell_List[Cell_Grid.Neighbor_List[i][n]][k]);
							}
						}
					}
				}
			}
		}
	}
}

void Calculate_Particle_Forces_3D(int Index_1, int Index_2)
{
	double Vector_X = Particles[Index_2].Position_X - Particles[Index_1].Position_X;
	double Vector_Y = Particles[Index_2].Position_Y - Particles[Index_1].Position_Y;
	double Vector_Z = Particles[Index_2].Position_Z - Particles[Index_1].Position_Z;

	double r_square = Vector_X * Vector_X + Vector_Y * Vector_Y + Vector_Z * Vector_Z;

	double base_calc;

	double Temp_Force_X;
	double Temp_Force_Y;
	double Temp_Force_Z;

	// WCA potential
	double f_wca = 0;

	if (r_square < Sigma_Test_WCA)
	{
		base_calc = Sigma_Square / r_square;

		base_calc *= base_calc * base_calc;

		f_wca = Epsilon_WCA_x_24 / r_square * (-2.0 * base_calc * base_calc + base_calc);

		Temp_Force_X = Vector_X * f_wca;
		Temp_Force_Y = Vector_Y * f_wca;
		Temp_Force_Z = Vector_Z * f_wca;

		m.lock();
		Particles[Index_1].Force_X += Temp_Force_X;
		Particles[Index_2].Force_X -= Temp_Force_X;

		Particles[Index_1].Force_Y += Temp_Force_Y;
		Particles[Index_2].Force_Y -= Temp_Force_Y;

		Particles[Index_1].Force_Z += Temp_Force_Z;
		Particles[Index_2].Force_Z -= Temp_Force_Z;
		m.unlock();

		//cout << "3D " << Index_1 << " " << Particles[Index_1].Position_X << " " << Particles[Index_1].Position_Y << " " << Particles[Index_1].Position_Z << " " <<
		//	Particles[Index_1].Velocity_X << " " << Particles[Index_1].Velocity_Y << " " << Particles[Index_1].Velocity_Z << " " << Vector_X * f_wca << endl;
	}
}