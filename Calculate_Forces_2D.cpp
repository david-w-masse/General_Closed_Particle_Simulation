// Single threaded calculation functions
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

void Calculate_Forces_Looped_Threaded_2D(int Thread_Num)
{
	a.lock();

	Active_Threads += 1;

	if (Active_Threads == Num_Threads)
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

		Calculate_Forces_Threaded_2D(Thread_Num);

		// wait until all forces calculated, then update time tick
		Force_Completed += (Thread_Num + 1);

		while (Force_Completed != Thread_Complete_Test) {}

		for (int i = Update_Particles_Start[Thread_Num]; i < Update_Particles_End[Thread_Num]; i++)
		{
			Particles[i].Past_Position_X = Particles[i].Position_X;
			Particles[i].Past_Position_Y = Particles[i].Position_Y;

			// Debug_Info(0, i);

			Particles[i].Velocity_X += 0.5 * (Particles[i].Force_X + Particles[i].Past_Force_X) * Delta_t;
			Particles[i].Velocity_Y += 0.5 * (Particles[i].Force_Y + Particles[i].Past_Force_Y) * Delta_t;

			Particles[i].Position_X += Particles[i].Velocity_X * Delta_t + 0.5 * Particles[i].Force_X * Delta_t_square;
			Particles[i].Position_Y += Particles[i].Velocity_Y * Delta_t + 0.5 * Particles[i].Force_Y * Delta_t_square;

			Particles[i].Past_Force_X = Particles[i].Force_X;
			Particles[i].Past_Force_Y = Particles[i].Force_Y;

			Particles[i].Force_X = 0;
			Particles[i].Force_Y = 0;
		}

		Update_Completed += (Thread_Num + 1);
	}
}

void Build_Cell_List_2D()
{
	double Local_Cell_X, Local_Cell_Y;
	double Local_Cell_X_2, Local_Cell_Y_2;
	int Local_Cell_X_3, Local_Cell_Y_3;

	int Final_Cell;

	for (int i = 0; i < Cell_Grid.Total_Divisions; i++)
	{
		Cell_Grid.Particles_Per_Cell[i] = 0;
	}

	for (int i = 0; i < Number_Particles; i++)
	{
		// assign each particle to a cell
		// keeping track of the counts per cell

		Local_Cell_X = Particles[i].Position_X - Box.Center_X;
		Local_Cell_Y = Particles[i].Position_Y - Box.Center_Y;

		Local_Cell_X_2 = (Local_Cell_X * Box.Cos_Angle_Z + Local_Cell_Y * Box.Sin_Angle_Z);
		Local_Cell_Y_2 = (-Local_Cell_X * Box.Sin_Angle_Z + Local_Cell_Y * Box.Cos_Angle_Z);

		Local_Cell_X_2 += Box.Base_Center_X;
		Local_Cell_Y_2 += Box.Base_Center_Y;

		Local_Cell_X_3 = (int)(Local_Cell_X_2 / Cell_Grid.X_Divisor);
		Local_Cell_Y_3 = (int)(Local_Cell_Y_2 / Cell_Grid.Y_Divisor);

		Local_Cell_Y_3 *= Cell_Grid.Divisions[0];

		Final_Cell = Local_Cell_Y_3 + Local_Cell_X_3;

		/*
		if ((Final_Cell >= Cell_Grid.Total_Divisions) || (Final_Cell < 0))
		{
			cout << Local_Cell_X << " " << Local_Cell_Y << " " << Local_Cell_X_2 << " " << Local_Cell_Y_2 << " " << Final_Cell << " " <<
				Box.Cos_Angle_Z << " " << Box.Sin_Angle_Z << " " << Particles[i].Velocity_X << " " << Particles[i].Velocity_Y << endl;

			int b = 0;
			cin >> b;
		}
		*/

		Cell_Grid.Particle_Cell_List[Final_Cell][Cell_Grid.Particles_Per_Cell[Final_Cell]] = i;

		Cell_Grid.Particles_Per_Cell[Final_Cell] += 1;
	}
}

void Calculate_Forces_2D()
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
					Calculate_Particle_Forces_2D(Cell_Grid.Particle_Cell_List[i][j], Cell_Grid.Particle_Cell_List[i][k]);
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
								Calculate_Particle_Forces_2D(Cell_Grid.Particle_Cell_List[i][j], Cell_Grid.Particle_Cell_List[Cell_Grid.Neighbor_List[i][n]][k]);
							}
						}
					}
				}
			}

			if ((Cell_Grid.Segments_Per_Cell[i] != 0) || (Cell_Grid.Points_Per_Cell[i] != 0))
			{
				Calculate_Boundary_Force_2D(i);
			}
		}
	}
}

void Calculate_Forces_Threaded_2D(int Thread_Num)
{
	for (int i = Update_Divisions_Start[Thread_Num]; i < Update_Divisions_End[Thread_Num]; i++)
	{
		if (Cell_Grid.Particles_Per_Cell[i] != 0)
		{
			// note this only runs if more than 1 particle is in the cell
			for (int j = 0; j < Cell_Grid.Particles_Per_Cell[i] - 1; j++)
			{
				// calculate forces on particles in the current cell with each neighbor in this cell
				for (int k = j + 1; k < Cell_Grid.Particles_Per_Cell[i]; k++)
				{
					Calculate_Particle_Forces_2D(Cell_Grid.Particle_Cell_List[i][j], Cell_Grid.Particle_Cell_List[i][k]);
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
								Calculate_Particle_Forces_2D(Cell_Grid.Particle_Cell_List[i][j], Cell_Grid.Particle_Cell_List[Cell_Grid.Neighbor_List[i][n]][k]);
							}
						}
					}
				}
			}

			if ((Cell_Grid.Segments_Per_Cell[i] != 0) || (Cell_Grid.Points_Per_Cell[i] != 0))
			{
				Calculate_Boundary_Force_2D(i);
			}
		}
	}
}

void Calculate_Particle_Forces_2D(int Index_1, int Index_2)
{
	double Vector_X = Particles[Index_2].Position_X - Particles[Index_1].Position_X;
	double Vector_Y = Particles[Index_2].Position_Y - Particles[Index_1].Position_Y;

	double r_square = Vector_X * Vector_X + Vector_Y * Vector_Y;

	double base_calc;

	// WCA potential
	double f_wca = 0;

	if (r_square < Sigma_Test_WCA)
	{
		base_calc = Sigma_Square / r_square;

		base_calc *= base_calc * base_calc;

		f_wca = Epsilon_WCA_x_24 / r_square * (-2.0 * base_calc * base_calc + base_calc);

		// normal

		double Temp_Force_X = Vector_X * f_wca;
		double Temp_Force_Y = Vector_Y * f_wca;

		m.lock();
		Particles[Index_1].Force_X += Temp_Force_X;
		Particles[Index_2].Force_X -= Temp_Force_X;

		Particles[Index_1].Force_Y += Temp_Force_Y;
		Particles[Index_2].Force_Y -= Temp_Force_Y;
		m.unlock();
	}
}

void Calculate_Boundary_Force_2D(int Index_1)
{
	int Index_2;

	double Vector_X;
	double Vector_Y;

	double Check_X;
	double Check_Y;

	double r_square;

	double base_calc;

	// WCA potential
	double f_wca = 0;

	double u_1, u_2, v_1, v_2;
	double x_1, x_2, y_1, y_2;

	double t;

	double x_int, y_int;

	double x_diff, y_diff;

	double Temp_Force_X, Temp_Force_Y;

	double r_x, r_y;

	int seg_index;

	for (int i = 0; i < Cell_Grid.Particles_Per_Cell[Index_1]; i++)
	{
		Index_2 = Cell_Grid.Particle_Cell_List[Index_1][i];

		Vector_X = Particles[Index_2].Position_X - Box.Center_X;
		Vector_Y = Particles[Index_2].Position_Y - Box.Center_Y;

		Check_X = Vector_X * Box.Cos_Angle_Z + Vector_Y * Box.Sin_Angle_Z;
		Check_Y = -Vector_X * Box.Sin_Angle_Z + Vector_Y * Box.Cos_Angle_Z;

		Check_X += Box.Base_Center_X;
		Check_Y += Box.Base_Center_Y;

		for (int j = 0; j < Cell_Grid.Segments_Per_Cell[Index_1]; j++)
		{
			seg_index = Cell_Grid.Segments[Index_1][j];
			
			u_1 = Box.Segments[seg_index].Vertex_0_X;
			u_2 = Box.Segments[seg_index].Vertex_1_X;

			v_1 = Box.Segments[seg_index].Vertex_0_Y;
			v_2 = Box.Segments[seg_index].Vertex_1_Y;

			x_1 = Check_X;
			y_1 = Check_Y;

			x_2 = x_1 + Box.Segments[seg_index].Normal_X;
			y_2 = y_1 + Box.Segments[seg_index].Normal_Y;

			t = -((x_1 - x_2) * (y_1 - v_1) - (y_1 - y_2) * (x_1 - u_1)) / ((x_1 - x_2) * (v_1 - v_2) - (y_1 - y_2) * (u_1 - u_2));

			if ((t >= 0) && (t < 1))
			{
				x_int = u_1 + t * (u_2 - u_1);
				y_int = v_1 + t * (v_2 - v_1);

				x_diff = x_int - x_1;
				y_diff = y_int - y_1;

				r_square = x_diff * x_diff + y_diff * y_diff;

				if (r_square < Sigma_Test_WCA)
				{
					base_calc = Sigma_Square / r_square;

					base_calc *= base_calc * base_calc;

					f_wca = Epsilon_WCA_x_24 / r_square * (-2.0 * base_calc * base_calc + base_calc);

					x_1 = -f_wca * x_diff;
					y_1 = -f_wca * y_diff;

					// calculate torque in rotated coordinates, as its the same
					r_x = x_int - Box.Base_Center_X;
					r_y = y_int - Box.Base_Center_Y;

					// now calculate forces in original coordinates
					Temp_Force_X = x_1 * Box.Cos_Angle_Z - y_1 * Box.Sin_Angle_Z;
					Temp_Force_Y = x_1 * Box.Sin_Angle_Z + y_1 * Box.Cos_Angle_Z;

					m.lock();
					Box.Torque_Z += r_x * y_1 - r_y * x_1;

					Particles[Index_2].Force_X -= Temp_Force_X;
					Particles[Index_2].Force_Y -= Temp_Force_Y;

					Box.Force_X += Temp_Force_X;
					Box.Force_Y += Temp_Force_Y;
					m.unlock();

					//cout << Temp_Force_X << " " << Temp_Force_Y << " " << Particles[Index_2].Position_X << " " << Particles[Index_2].Position_Y << endl;
				}
			}
		}

		for (int j = 0; j < Cell_Grid.Points_Per_Cell[Index_1]; j++)
		{
			seg_index = Cell_Grid.Points_2D[Index_1][j];
			
			u_1 = Box.Points_2D[seg_index].Vertex_X;
			v_1 = Box.Points_2D[seg_index].Vertex_Y;

			x_1 = u_1 - Check_X;
			y_1 = v_1 - Check_Y;

			r_square = x_1 * x_1 + y_1 * y_1;

			if (r_square < Sigma_Test_WCA)
			{
				//cout << "point " << Particles[Index_2].Position_X << " " << Particles[Index_2].Position_Y << endl;

				// cross product with first normal should be positive, second should be negative
				double cross_1 = Box.Points_2D[seg_index].Normal_1_X * y_1 - Box.Points_2D[seg_index].Normal_1_Y * x_1;

				if (cross_1 > 0)
				{
					double cross_2 = Box.Points_2D[seg_index].Normal_0_X * y_1 - Box.Points_2D[seg_index].Normal_0_Y * x_1;

					if (cross_2 <= 0)
					{
						base_calc = Sigma_Square / r_square;

						base_calc *= base_calc * base_calc;

						f_wca = Epsilon_WCA_x_24 / r_square * (-2.0 * base_calc * base_calc + base_calc);

						x_1 = f_wca * x_1;
						y_1 = f_wca * y_1;

						// calculate torque in rotated coordinates, as its the same
						r_x = u_1 - Box.Base_Center_X;
						r_y = v_1 - Box.Base_Center_Y;

						// now calculate forces in original coordinates
						Temp_Force_X = x_1 * Box.Cos_Angle_Z - y_1 * Box.Sin_Angle_Z;
						Temp_Force_Y = x_1 * Box.Sin_Angle_Z + y_1 * Box.Cos_Angle_Z;

						m.lock();
						Box.Torque_Z -= r_x * y_1 - r_y * x_1;

						Particles[Index_2].Force_X += Temp_Force_X;
						Particles[Index_2].Force_Y += Temp_Force_Y;

						Box.Force_X -= Temp_Force_X;
						Box.Force_Y -= Temp_Force_Y;
						m.unlock();
					}
				}
			}
		}
	}
}