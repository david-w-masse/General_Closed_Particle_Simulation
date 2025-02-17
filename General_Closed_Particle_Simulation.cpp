// Particle_Simulation.cpp : Defines the entry point for the application.
//

/*
	NOTES:
	Don't use math.pow for powers of 2 or 3, it is very slow

	multithreading saves ~25% time for 200 particles in 2D with 4 threads
*/
#pragma once

#include <chrono>
#include <random>
#include <cmath>
#include <fstream>
#include <stdexcept>

//#include <emmintrin.h>
#include <immintrin.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <vector>

#include "General_Closed_Particle_Simulation.h"

using namespace std;

Cell_Structure Cell_Grid;

Box_Structure Box;

Particle_Structure* Particles;

// simulation variables
double Particle_Diameter;
double End_Time;
double Delta_t;
double Epsilon_WCA;

int Dimension;
int Number_Particles;
int Num_Threads;

bool Single_Thread;

double Region_Size[3];

double Sigma_Square;
double Sigma_Test_WCA;
double Sigma_Test_LJ;
double Epsilon_WCA_x_24;
double Delta_t_square;
double Half_Delta_t_square;
double Half_Delta_t_square_div_Box_Mass;
double Half_Delta_t;
double Half_Delta_t_div_Box_Mass;
double Half_Delta_t_square_div_Inertia_Z;
double Half_Delta_t_div_Inertia_Z;

double Current_Time;

int Neighbor_Count;

mutex m;

mutex ctrl[4];

mutex a;

mutex cmp;

volatile atomic_bool* Ready;

int* Update_Divisions_Start;
int* Update_Divisions_End;

int* Update_Particles_Start;
int* Update_Particles_End;

volatile atomic_bool Finished;

volatile atomic_int Update_Completed;

volatile atomic_int Force_Completed;

volatile atomic_int Loop_Completed;

int Active_Threads;

int Thread_Complete_Test;

condition_variable Active;

int main(int argc, char* argv[])
{
	int b;
	
	// load simulation parameters from file
	// Also sets up basic simulation data
	bool Data_Load_Success = Load_Simulation_Variables();

	if (!Data_Load_Success)
	{
		cout << "Unable to load simulation data, exiting." << endl;
		return 0;
	}

	Build_Memory();

	// load geometry data from file
	bool Geom_Load_Success = false;
	
	if (Dimension == 2)
	{
		Geom_Load_Success = Load_Geometry_Data_2D();
	}
	else
	{
		Geom_Load_Success = Load_Geometry_Data_3D();
	}

	cout << "good geometry " << Box.Center_X << " " << Box.Center_Y << endl;

	if (!Geom_Load_Success)
	{
		cout << "Unable to load geometry data, exiting." << endl;

		Clear_Memory();
		return 0;
	}

	if (Num_Threads == 1)
	{
		Single_Thread = true;
	}
	else
	{	
		Thread_Complete_Test = 0;
		
		for (int i = 0; i < Num_Threads; i++)
		{
			if (i == 0)
			{
				Update_Particles_Start[i] = 0;
			}
			else
			{
				Update_Particles_Start[i] = Update_Particles_End[i - 1];
			}

			if (i == (Num_Threads - 1))
			{
				Update_Particles_End[i] = Number_Particles;
			}
			else
			{
				Update_Particles_End[i] = Number_Particles / Num_Threads * (i + 1);
			}

			Thread_Complete_Test += (i + 1);
		}

		for (int i = 0; i < Num_Threads; i++)
		{
			if (i == 0)
			{
				Update_Divisions_Start[i] = 0;
			}
			else
			{
				Update_Divisions_Start[i] = Update_Divisions_End[i - 1];
			}

			if (i == (Num_Threads - 1))
			{
				Update_Divisions_End[i] = Cell_Grid.Total_Divisions;
			}
			else
			{
				Update_Divisions_End[i] = Cell_Grid.Divisions[0] * Cell_Grid.Y_Div_T * (i + 1);
			}
		}
	}

	cout << "memory allocated" << endl;

	Assign_Particle_Data();

	Assign_Box_Data();

	Build_Boundary_List();

	cout << "Data Good" << endl;

	Build_Neighbor_List();

	if (Dimension == 2)
	{
		Initialize_Positions_2D(true);
		//Initialize_Positions_Unit_Tests_2D();
	}
	else
	{
		Initialize_Positions_3D(true);
		//Initialize_Positions_Unit_Tests_3D();
	}

	cout << "positions assigned" << endl;

	//cout << Box.Points_2D[0].Use_Point << " " << Box.Points_2D[1].Use_Point << " " << Box.Points_2D[2].Use_Point << " " << Box.Points_2D[3].Use_Point << " " <<
	//	Box.Points_2D[4].Use_Point << " " << Box.Points_2D[5].Use_Point << " " << Box.Points_2D[6].Use_Point << " " << Box.Points_2D[7].Use_Point << endl;
	
	Output_Check_File();
	/*
	Clear_Memory();
	return 0;
	*/

	std::random_device Get_Seed;
	std::mt19937 Random_Normal(Get_Seed());
	std::normal_distribution<double> Random_Force(0, 1);

	std::mt19937 Random_Uniform(Get_Seed());
	std::uniform_real_distribution<double> Random_Reaction(0, 1);

	double chance = 0;

	int steps = 0;

	string iter_val = "";

	if (argc == 2)
	{
		iter_val = argv[1];
	}

	string Position_File_Name = "output/Position" + iter_val + ".txt";
	string Velocity_File_Name = "output/Velocity" + iter_val + ".txt";
	string Box_File_Name = "output/Box" + iter_val + ".txt";

	ofstream Output_File;
	Output_File.open(Position_File_Name, ios::out | ios::trunc);

	ofstream Velocity_File;
	Velocity_File.open(Velocity_File_Name, ios::out | ios::trunc);

	ofstream Box_File;
	Box_File.open(Box_File_Name, ios::out | ios::trunc);

	if (Dimension == 2)
	{
		for (int i = 0; i < Number_Particles; i++)
		{
			Output_File << Particles[i].Position_X << "\t" << Particles[i].Position_Y << "\n";
		}

		for (int i = 0; i < Number_Particles; i++)
		{
			Velocity_File << Particles[i].Velocity_X << "\t" << Particles[i].Velocity_Y << "\n";
		}

		Box_File << Box.Center_X << "\t" << Box.Center_Y << "\t" << Box.Velocity_X << "\t" << Box.Velocity_Y << "\t" 
				 << Box.Angle_Z << "\t" << Box.Angular_Velocity_Z << "\n";
	}
	else
	{
		for (int i = 0; i < Number_Particles; i++)
		{
			Output_File << Particles[i].Position_X << "\t" << Particles[i].Position_Y << "\t" << Particles[i].Position_Z << "\n";
		}

		for (int i = 0; i < Number_Particles; i++)
		{
			Velocity_File << Particles[i].Velocity_X << "\t" << Particles[i].Velocity_Y << "\t" << Particles[i].Velocity_Z << "\n";
		}

		Box_File << Box.Center_X << "\t" << Box.Center_Y << "\t" << Box.Center_Z << "\t" << Box.Velocity_X << "\t" << Box.Velocity_Y << "\t" << Box.Velocity_Z << "\t"
			<< Box.Angle_X << "\t" << Box.Angle_Y << "\t" << Box.Angle_Z << "\t"
			<< Box.Angular_Velocity_X << "\t" << Box.Angular_Velocity_Y << "\t" << Box.Angular_Velocity_Z << "\n";
	}

	double delta_x = 0;
	double delta_y = 0;
	double delta_z = 0;

	cout << "ini complete" << endl;

	Current_Time = 0;

	double check_time = 0;

	for (int i = 0; i < 10000 * 10000; i++)
	{
		check_time += Delta_t;
	}

	auto start = std::chrono::system_clock::now();

	__m256d vec_delta_t = _mm256_set_pd(Delta_t, Delta_t, Delta_t, Delta_t);

	__m256d vec_half_delta_t = _mm256_set_pd(Half_Delta_t, Half_Delta_t, Half_Delta_t, Half_Delta_t);

	__m256d vec_half_delta_t_square = _mm256_set_pd(Half_Delta_t_square, Half_Delta_t_square, Half_Delta_t_square, Half_Delta_t_square);

	__m256d vec_mul_forces = _mm256_set_pd(0, 0, 0, 0);

	__m256d vec_mul_forces_p = _mm256_set_pd(0, 0, 0, 0);

	__m256d vec_mul_vel = _mm256_set_pd(0, 0, 0, 0);

	__m256d vec_mul_step_1 = _mm256_mul_pd(vec_mul_vel, vec_delta_t);

	__m256d vec_mul_step_2 = _mm256_fmadd_pd(vec_mul_forces, vec_half_delta_t_square, vec_mul_step_1);

	double* vec_s_1 = (double*)&vec_mul_step_2;


	if (Single_Thread)
	{
		if (Dimension == 2)
		{
			// velocity verlet integrator
			Build_Cell_List_2D();

			Calculate_Forces_2D();

			while (Current_Time < End_Time)
			{

				for (int i = 0; i < Number_Particles; i++)
				{
					// normal			
					Particles[i].Position_X += Particles[i].Velocity_X * Delta_t + Particles[i].Force_X * Half_Delta_t_square;
					Particles[i].Position_Y += Particles[i].Velocity_Y * Delta_t + Particles[i].Force_Y * Half_Delta_t_square;

					// Debug_Info(0, i);

					Particles[i].Past_Force_X = Particles[i].Force_X;
					Particles[i].Past_Force_Y = Particles[i].Force_Y;

					Particles[i].Force_X = 0;
					Particles[i].Force_Y = 0;
				}
				
				/*
				for (int i = 0; i < Number_Particles; i+=4)
				{
					// AVX
					
					vec_mul_forces = _mm256_set_pd(Particles[i].Force_X, Particles[i + 1].Force_X, Particles[i + 2].Force_X, Particles[i + 3].Force_X);

					vec_mul_vel = _mm256_set_pd(Particles[i].Velocity_X, Particles[i + 1].Velocity_X, Particles[i + 2].Velocity_X, Particles[i + 3].Velocity_X);

					vec_mul_step_1 = _mm256_mul_pd(vec_mul_vel, vec_delta_t);

					vec_mul_step_2 = _mm256_fmadd_pd(vec_mul_forces, vec_half_delta_t_square, vec_mul_step_1);

					vec_s_1 = (double*)&vec_mul_step_2;

					for (int j = 0; j < 4; j++)
					{
						Particles[i + j].Position_X += vec_s_1[3-j];
						Particles[i + j].Past_Force_X = Particles[i + j].Force_X;
						Particles[i + j].Force_X = 0;
					}

					vec_mul_forces = _mm256_set_pd(Particles[i].Force_Y, Particles[i + 1].Force_Y, Particles[i + 2].Force_Y, Particles[i + 3].Force_Y);

					vec_mul_vel = _mm256_set_pd(Particles[i].Velocity_Y, Particles[i + 1].Velocity_Y, Particles[i + 2].Velocity_Y, Particles[i + 3].Velocity_Y);

					vec_mul_step_1 = _mm256_mul_pd(vec_mul_vel, vec_delta_t);

					vec_mul_step_2 = _mm256_fmadd_pd(vec_mul_forces, vec_half_delta_t_square, vec_mul_step_1);

					vec_s_1 = (double*)&vec_mul_step_2;

					for (int j = 0; j < 4; j++)
					{
						Particles[i + j].Position_Y += vec_s_1[3-j];
						Particles[i + j].Past_Force_Y = Particles[i + j].Force_Y;
						Particles[i + j].Force_Y = 0;
					}
				}
				*/
				// update box variables
				Box.Center_X += Box.Velocity_X * Delta_t + Box.Force_X * Half_Delta_t_square_div_Box_Mass;
				Box.Center_Y += Box.Velocity_Y * Delta_t + Box.Force_Y * Half_Delta_t_square_div_Box_Mass;

				Box.Angle_Z += Box.Angular_Velocity_Z * Delta_t + Box.Torque_Z * Half_Delta_t_square_div_Inertia_Z;

				Box.Cos_Angle_Z = cos(Box.Angle_Z);
				Box.Sin_Angle_Z = sin(Box.Angle_Z);

				Box.Past_Force_X = Box.Force_X;
				Box.Past_Force_Y = Box.Force_Y;

				Box.Past_Torque_Z = Box.Torque_Z;

				Box.Force_X = 0;
				Box.Force_Y = 0;

				Box.Torque_Z = 0;

				Build_Cell_List_2D();

				Calculate_Forces_2D();
				
				for (int i = 0; i < Number_Particles; i++)
				{
					Particles[i].Velocity_X += (Particles[i].Force_X + Particles[i].Past_Force_X) * Half_Delta_t;
					Particles[i].Velocity_Y += (Particles[i].Force_Y + Particles[i].Past_Force_Y) * Half_Delta_t;
				}
				
				/*
				for (int i = 0; i < Number_Particles; i += 4)
				{
					// AVX
					vec_mul_forces = _mm256_set_pd(Particles[i].Force_X, Particles[i + 1].Force_X, Particles[i + 2].Force_X, Particles[i + 3].Force_X);

					vec_mul_forces_p = _mm256_set_pd(Particles[i].Past_Force_X, Particles[i + 1].Past_Force_X, Particles[i + 2].Past_Force_X, Particles[i + 3].Past_Force_X);

					vec_mul_vel = _mm256_set_pd(Particles[i].Velocity_X, Particles[i + 1].Velocity_X, Particles[i + 2].Velocity_X, Particles[i + 3].Velocity_X);

					vec_mul_step_1 = _mm256_add_pd(vec_mul_forces, vec_mul_forces_p);

					vec_mul_step_2 = _mm256_mul_pd(vec_mul_step_1, vec_half_delta_t);

					vec_s_1 = (double*)&vec_mul_step_2;

					for (int j = 0; j < 4; j++)
					{
						Particles[i + j].Velocity_X += vec_s_1[3 - j];
					}

					vec_mul_forces = _mm256_set_pd(Particles[i].Force_Y, Particles[i + 1].Force_Y, Particles[i + 2].Force_Y, Particles[i + 3].Force_Y);

					vec_mul_forces_p = _mm256_set_pd(Particles[i].Past_Force_Y, Particles[i + 1].Past_Force_Y, Particles[i + 2].Past_Force_Y, Particles[i + 3].Past_Force_Y);

					vec_mul_vel = _mm256_set_pd(Particles[i].Velocity_Y, Particles[i + 1].Velocity_Y, Particles[i + 2].Velocity_Y, Particles[i + 3].Velocity_Y);

					vec_mul_step_1 = _mm256_add_pd(vec_mul_forces, vec_mul_forces_p);

					vec_mul_step_2 = _mm256_mul_pd(vec_mul_step_1, vec_half_delta_t);

					vec_s_1 = (double*)&vec_mul_step_2;

					for (int j = 0; j < 4; j++)
					{
						Particles[i + j].Velocity_Y += vec_s_1[3 - j];
					}
				}
				*/

				Box.Velocity_X += (Box.Force_X + Box.Past_Force_X) * Half_Delta_t_div_Box_Mass;
				Box.Velocity_Y += (Box.Force_Y + Box.Past_Force_Y) * Half_Delta_t_div_Box_Mass;

				Box.Angular_Velocity_Z += (Box.Torque_Z + Box.Past_Torque_Z) * Half_Delta_t_div_Inertia_Z;

				Current_Time += Delta_t;

				steps += 1;

				if ((steps % 10000) == 0)
				{
					for (int i = 0; i < Number_Particles; i++)
					{
						Output_File << Particles[i].Position_X << "\t" << Particles[i].Position_Y << "\n";
					}

					for (int i = 0; i < Number_Particles; i++)
					{
						Velocity_File << Particles[i].Velocity_X << "\t" << Particles[i].Velocity_Y << "\n";
					}

					Box_File << Box.Center_X << "\t" << Box.Center_Y << "\t" << Box.Velocity_X << "\t" << Box.Velocity_Y << "\t"
						<< Box.Angle_Z << "\t" << Box.Angular_Velocity_Z << "\n";
					/*
					if ((steps % 100000) == 0)
					{
						Debug_Info(1,0);
					}
					*/

					if (Current_Time == check_time)
					{
						Box.Angular_Velocity_Z += 0.002;
						cout << "inpulse" << endl;
					}
				}
			}
		}
		else
		{
			while (Current_Time < End_Time)
			{
				Build_Cell_List_3D();

				Calculate_Forces_3D();

				for (int i = 0; i < Number_Particles; i++)
				{
					// normal
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

				// update box variables
				Box.Center_X += Box.Velocity_X * Delta_t;
				Box.Center_Y += Box.Velocity_Y * Delta_t;
				Box.Center_Z += Box.Velocity_Z * Delta_t;

				Box.Velocity_X += Box.Force_X * Delta_t / Box.Mass;
				Box.Velocity_Y += Box.Force_Y * Delta_t / Box.Mass;
				Box.Velocity_Z += Box.Force_Z * Delta_t / Box.Mass;

				Box.Angle_Z += Box.Angular_Velocity_Z * Delta_t;

				Box.Angular_Velocity_Z += Box.Torque_Z * Delta_t / Box.Inertia_Z;

				Box.Cos_Angle_Z = cos(Box.Angle_Z);
				Box.Sin_Angle_Z = sin(Box.Angle_Z);

				Box.Force_X = 0;
				Box.Force_Y = 0;
				Box.Force_Z = 0;

				Box.Torque_Z = 0;

				Current_Time += Delta_t;


				Current_Time += Delta_t;

				steps += 1;

				if ((steps % 100000) == 0)
				{
					for (int i = 0; i < Number_Particles; i++)
					{
						Output_File << Particles[i].Position_X << "\t" << Particles[i].Position_Y << "\t" << Particles[i].Position_Z << "\n";
					}

					for (int i = 0; i < Number_Particles; i++)
					{
						Velocity_File << Particles[i].Velocity_X << "\t" << Particles[i].Velocity_Y << "\t" << Particles[i].Velocity_Z << "\n";
					}

					Box_File << Box.Center_X << "\t" << Box.Center_Y << "\t" << Box.Center_Z << "\t" << Box.Velocity_X << "\t" << Box.Velocity_Y << "\t" << Box.Velocity_Z << "\t"
						<< Box.Angle_X << "\t" << Box.Angle_Y << "\t" << Box.Angle_Z << "\t"
						<< Box.Angular_Velocity_X << "\t" << Box.Angular_Velocity_Y << "\t" << Box.Angular_Velocity_Z << "\n";
					/*
					if ((steps % 100000) == 0)
					{
						Debug_Info(1,0);
					}
					*/
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < Num_Threads; i++)
		{
			Ready[i] = false;
		}

		Force_Completed = 0;
		Loop_Completed = 0;
		Update_Completed = 0;

		Finished = false;

		void (*thread_func)(int);

		if (Dimension == 2)
		{
			thread_func = Calculate_Forces_Looped_Threaded_2D;
		}
		else
		{
			thread_func = Calculate_Forces_Looped_Threaded_3D;
		}

		vector<std::thread> Threads;

		for (int i = 0; i < Num_Threads; i++)
		{
			Threads.push_back(thread(thread_func, i));
		}

		unique_lock<std::mutex> lk(m);
		Active.wait_for(lk, chrono::milliseconds(10), [] { return Active_Threads == Num_Threads; });
		m.unlock();

		if (Active_Threads != Num_Threads)
		{
			cout << "unable to spawn threads, exitting." << endl;
			Clear_Memory();

			cin >> b;
			return 0;
		}

		cout << Update_Particles_Start[0] << " " << Update_Particles_Start[1] << " " << Update_Particles_Start[2] << " " << Update_Particles_Start[3] << endl;
		cout << Update_Particles_End[0] << " " << Update_Particles_End[1] << " " << Update_Particles_End[2] << " " << Update_Particles_End[3] << endl;

		cout << Cell_Grid.Divisions[0] * Cell_Grid.Y_Div_T * 0 << endl;
		cout << Cell_Grid.Divisions[0] * Cell_Grid.Y_Div_T * 1 << endl;
		cout << Cell_Grid.Divisions[0] * Cell_Grid.Y_Div_T * 2 << endl;
		cout << Cell_Grid.Divisions[0] * Cell_Grid.Y_Div_T * 3 << endl;

		if (Dimension == 2)
		{
			Build_Cell_List_2D();

			Calculate_Forces_2D();

			for (int i = 0; i < Number_Particles; i++)
			{
				Particles[i].Position_X += Particles[i].Velocity_X * Delta_t + Particles[i].Force_X * Half_Delta_t_square;
				Particles[i].Position_Y += Particles[i].Velocity_Y * Delta_t + Particles[i].Force_Y * Half_Delta_t_square;

				// Debug_Info(0, i);

				Particles[i].Past_Force_X = Particles[i].Force_X;
				Particles[i].Past_Force_Y = Particles[i].Force_Y;

				Particles[i].Force_X = 0;
				Particles[i].Force_Y = 0;
			}

			// update box variables
			Box.Center_X += Box.Velocity_X * Delta_t + Box.Force_X * Half_Delta_t_square_div_Box_Mass;
			Box.Center_Y += Box.Velocity_Y * Delta_t + Box.Force_Y * Half_Delta_t_square_div_Box_Mass;

			Box.Angle_Z += Box.Angular_Velocity_Z * Delta_t + Box.Torque_Z * Half_Delta_t_square_div_Inertia_Z;

			Box.Cos_Angle_Z = cos(Box.Angle_Z);
			Box.Sin_Angle_Z = sin(Box.Angle_Z);

			Box.Past_Force_X = Box.Force_X;
			Box.Past_Force_Y = Box.Force_Y;

			Box.Past_Torque_Z = Box.Torque_Z;

			Box.Force_X = 0;
			Box.Force_Y = 0;

			Box.Torque_Z = 0;

			while (Current_Time < End_Time)
			{				
				Build_Cell_List_2D();
				
				for (int i = 0; i < Num_Threads; i++)
				{
					Ready[i] = true;
				}

				while (Update_Completed != Thread_Complete_Test) {}

				Update_Completed = 0;
				Force_Completed = 0;

				//Debug_Info(0, 0);

				// due to update order, write box data before position updates to avoid having to save past data
				Box.Velocity_X += (Box.Force_X + Box.Past_Force_X) * Half_Delta_t_div_Box_Mass;
				Box.Velocity_Y += (Box.Force_Y + Box.Past_Force_Y) * Half_Delta_t_div_Box_Mass;

				Box.Angular_Velocity_Z += (Box.Torque_Z + Box.Past_Torque_Z) * Half_Delta_t_div_Inertia_Z;

				Current_Time += Delta_t;

				steps += 1;

				if ((steps % 10000) == 0)
				{
					for (int i = 0; i < Number_Particles; i++)
					{
						Output_File << Particles[i].Past_Position_X << "\t" << Particles[i].Past_Position_Y << "\n";
					}

					for (int i = 0; i < Number_Particles; i++)
					{
						Velocity_File << Particles[i].Velocity_X << "\t" << Particles[i].Velocity_Y << "\n";
					}

					Box_File << Box.Center_X << "\t" << Box.Center_Y << "\t" << Box.Velocity_X << "\t" << Box.Velocity_Y << "\t" 
							 << Box.Angle_Z << "\t" << Box.Angular_Velocity_Z << "\n";
					/*
					if ((steps % 100000) == 0)
					{
						Debug_Info(1,0);
					}
					*/
				}

				Box.Center_X += Box.Velocity_X * Delta_t + Box.Force_X * Half_Delta_t_square_div_Box_Mass;
				Box.Center_Y += Box.Velocity_Y * Delta_t + Box.Force_Y * Half_Delta_t_square_div_Box_Mass;

				Box.Angle_Z += Box.Angular_Velocity_Z * Delta_t + Box.Torque_Z * Half_Delta_t_square_div_Inertia_Z;

				Box.Cos_Angle_Z = cos(Box.Angle_Z);
				Box.Sin_Angle_Z = sin(Box.Angle_Z);

				Box.Past_Force_X = Box.Force_X;
				Box.Past_Force_Y = Box.Force_Y;

				Box.Past_Torque_Z = Box.Torque_Z;

				Box.Force_X = 0;
				Box.Force_Y = 0;

				Box.Torque_Z = 0;
			}
		}
		else
		{
			while (Current_Time < End_Time)
			{
				Build_Cell_List_3D();

				for (int i = 0; i < Num_Threads; i++)
				{
					Ready[i] = true;
				}

				while (Loop_Completed != 10) {}

				Force_Completed = 0;

				Loop_Completed = 0;

				// Debug_Info(0, i);

				Current_Time += Delta_t;

				steps += 1;

				// update box variables
				Box.Center_X += Box.Velocity_X * Delta_t;
				Box.Center_Y += Box.Velocity_Y * Delta_t;
				Box.Center_Z += Box.Velocity_Z * Delta_t;

				Box.Velocity_X += Box.Force_X * Delta_t / Box.Mass;
				Box.Velocity_Y += Box.Force_Y * Delta_t / Box.Mass;
				Box.Velocity_Z += Box.Force_Z * Delta_t / Box.Mass;

				Box.Angle_Z += Box.Angular_Velocity_Z * Delta_t;

				Box.Angular_Velocity_Z += Box.Torque_Z * Delta_t / Box.Inertia_Z;

				Box.Cos_Angle_Z = cos(Box.Angle_Z);
				Box.Sin_Angle_Z = sin(Box.Angle_Z);

				Box.Force_X = 0;
				Box.Force_Y = 0;
				Box.Force_Z = 0;

				Box.Torque_Z = 0;

				//cout << endl;

				if ((steps % 100000) == 0)
				{
					for (int i = 0; i < Number_Particles; i++)
					{
						Output_File << Particles[i].Position_X << "\t" << Particles[i].Position_Y << "\t" << Particles[i].Position_Z << "\n";
					}

					for (int i = 0; i < Number_Particles; i++)
					{
						Velocity_File << Particles[i].Velocity_X << "\t" << Particles[i].Velocity_Y << "\t" << Particles[i].Velocity_Z << "\n";
					}

					Box_File << Box.Center_X << "\t" << Box.Center_Y << "\t" << Box.Center_Z << "\t" << Box.Velocity_X << "\t" << Box.Velocity_Y << "\t" << Box.Velocity_Z << "\t"
						<< Box.Angle_X << "\t" << Box.Angle_Y << "\t" << Box.Angle_Z << "\t"
						<< Box.Angular_Velocity_X << "\t" << Box.Angular_Velocity_Y << "\t" << Box.Angular_Velocity_Z << "\n";
					/*
					if ((steps % 100000) == 0)
					{
						Debug_Info(1,0);
					}
					*/
				}
			}
		}

		Finished = true;

		for (int i = 0; i < Num_Threads; i++)
		{
			Ready[i] = true;
		}

		for (int i = 0; i < Num_Threads; i++)
		{
			Threads[i].join();
		}
	}

	int a;

	Output_File.close();
	Velocity_File.close();
	Box_File.close();

	auto end = std::chrono::system_clock::now();
	auto elapsed = end - start;
	std::cout << elapsed.count() << '\n';

	Clear_Memory();

	cin >> a;

	return 0;
}

void Debug_Info(int t, int i)
{
	int b;

	double c_m_x, c_m_y;


	if (t == 0)
	{
		if (Particles[i].Position_X >= Region_Size[0])
		{
			cout << "bad X + " << i << " " << Particles[i].Position_X << " " << Particles[i].Position_Y << " " << Particles[i].Position_Z << " " <<
				Particles[i].Velocity_X << " " << Particles[i].Velocity_Y << " " << Particles[i].Velocity_Z << " " << Current_Time << endl;
			cin >> b;
		}

		if (Particles[i].Position_Y >= Region_Size[1])
		{
			cout << "bad Y + " << i << " " << Particles[i].Position_X << " " << Particles[i].Position_Y << " " << Particles[i].Position_Z << " " <<
				Particles[i].Velocity_X << " " << Particles[i].Velocity_Y << " " << Particles[i].Velocity_Z << " " << Current_Time << endl;
			cin >> b;
		}

		if (Particles[i].Position_Z >= Region_Size[2])
		{
			cout << "bad Z + " << i << " " << Particles[i].Position_X << " " << Particles[i].Position_Y << " " << Particles[i].Position_Z << " " <<
				Particles[i].Velocity_X << " " << Particles[i].Velocity_Y << " " << Particles[i].Velocity_Z << " " << Current_Time << endl;
			cin >> b;
		}

		if (Particles[i].Position_X < 0)
		{
			cout << "bad X - " << i << " " << Particles[i].Position_X << " " << Particles[i].Position_Y << " " << Particles[i].Position_Z << " " <<
				Particles[i].Velocity_X << " " << Particles[i].Velocity_Y << " " << Particles[i].Velocity_Z << " " << Current_Time << endl;
			cin >> b;
		}

		if (Particles[i].Position_Y < 0)
		{
			cout << "bad Y - " << i << " " << Particles[i].Position_X << " " << Particles[i].Position_Y << " " << Particles[i].Position_Z << " " <<
				Particles[i].Velocity_X << " " << Particles[i].Velocity_Y << " " << Particles[i].Velocity_Z << " " << Current_Time << endl;
			cin >> b;
		}

		if (Particles[i].Position_Z < 0)
		{
			cout << "bad Z - " << i << " " << Particles[i].Position_X << " " << Particles[i].Position_Y << " " << Particles[i].Position_Z << " " <<
				Particles[i].Velocity_X << " " << Particles[i].Velocity_Y << " " << Particles[i].Velocity_Z << " " << Current_Time << endl;
			cin >> b;
		}

		cout << "loop" << endl;
	}
	else if (t == 1)
	{
		// find center of mass
		c_m_x = 0;
		c_m_y = 0;

		for (int j = 0; j < Number_Particles; j++)
		{
			c_m_x += Particles[j].Position_X;
			c_m_y += Particles[j].Position_Y;
		}

		c_m_x += Box.Center_X * Number_Particles;
		c_m_y += Box.Center_Y * Number_Particles;

		c_m_x = c_m_x / (2 * Number_Particles);
		c_m_y = c_m_y / (2 * Number_Particles);

		cout << c_m_x << " " << c_m_y << endl;
	}
}