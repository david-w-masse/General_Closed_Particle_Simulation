// Initialization code
#pragma once

#include <chrono>
#include <random>
#include <cmath>
#include <fstream>
#include <stdexcept>

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

extern int Neighbor_Count;

void Initialize_Positions_3D(bool zeroize)
{
	std::random_device Get_Seed;
	std::mt19937 Random_Generator(111);
	std::uniform_real_distribution<double> Random_X(2, Region_Size[0] - 2);
	std::uniform_real_distribution<double> Random_Y(2, Region_Size[1] - 2);
	std::uniform_real_distribution<double> Random_Z(2, Region_Size[2] - 2);

	std::uniform_real_distribution<double> Random_VX(-0.1, 0.1);
	std::uniform_real_distribution<double> Random_VY(-0.1, 0.1);
	std::uniform_real_distribution<double> Random_VZ(-0.1, 0.1);

	std::uniform_real_distribution<double> Random_Angle(0, 6.283185307179);
	std::uniform_real_distribution<double> Random_Scale(0.5, Region_Size[0] / 1.6);
	std::uniform_real_distribution<double> Random_Parity(0, 1);
	std::uniform_real_distribution<double> Random_Option(0, 4);

	bool Good_Position;

	double Temp_X, Temp_Y, Temp_Z;
	double Temp_X_2, Temp_Y_2, Temp_Z_2;
	double Temp_X_3, Temp_Y_3, Temp_Z_3;
	double Temp_X_4, Temp_Y_4, Temp_Z_4;

	double Temp_V_X, Temp_V_Y;
	double Temp_V_X_2, Temp_V_Y_2;
	double Temp_V_X_3, Temp_V_Y_3;
	double Temp_V_X_4, Temp_V_Y_4;

	double Temp_Vel;

	double Temp_Angle, Temp_Scale, Temp_Parity;

	double Current_X, Current_Y, Current_Z;

	double Temp_Dist_Min;

	int Tries;

	bool option_2, option_3, option_4;

	bool good_choice;

	double pick_choice;

	int Particles_Remaining = Number_Particles;

	int Particles_Placed = 0;

	double boundary_dist = Particle_Diameter * 1.175;

	if (Particles_Remaining == 1)
	{
		cout << "running without zeroizing" << endl;
		zeroize = false;
	}

	// zeroize means no momentum or angular momentum and mass center is box center
	// uses a combination of asymmetric and symmetric combinations
	if (zeroize)
	{
		for (int i = 0; i < Number_Particles; i++)
		{
			Good_Position = false;

			Tries = 0;

			while (!Good_Position)
			{
				Temp_X = Random_X(Random_Generator);
				Temp_Y = Random_Y(Random_Generator);
				Temp_Z = Random_Y(Random_Generator);

				Good_Position = true;

				for (int j = 0; j < i; j++)
				{
					Temp_Dist_Min = (Particles[i].Size + Particles[j].Size) * 1.175;

					Current_X = Particles[j].Position_X;
					Current_Y = Particles[j].Position_Y;
					Current_Z = Particles[j].Position_Z;

					// check center and all boundaries for every particle
					Good_Position &= Determine_Proper_Seperation_3D(Current_X, Current_Y, Current_Z, Temp_X, Temp_Y, Temp_Z, Temp_Dist_Min);
				}

				if (!Good_Position)
				{
					Tries += 1;

					if ((Tries == 1000000) && !Good_Position)
					{
						throw std::invalid_argument("Error: Cannot place particle.");
					}
				}
			}

			Particles[i].Position_X = Temp_X;
			Particles[i].Position_Y = Temp_Y;
			Particles[i].Position_Z = Temp_Z;

			Particles[i].Velocity_X = Random_VX(Random_Generator);
			Particles[i].Velocity_Y = Random_VY(Random_Generator);
			Particles[i].Velocity_Z = Random_VZ(Random_Generator);
		}
	}
	else
	{
		for (int i = 0; i < Number_Particles; i++)
		{
			Good_Position = false;

			Tries = 0;

			while (!Good_Position)
			{
				Temp_X = Random_X(Random_Generator);
				Temp_Y = Random_Y(Random_Generator);
				Temp_Z = Random_Y(Random_Generator);

				Good_Position = Inside_Boundaries_3D(Temp_X, Temp_Y, Temp_Z, boundary_dist);

				for (int j = 0; j < i; j++)
				{
					Temp_Dist_Min = (Particles[i].Size + Particles[j].Size) * 1.175;

					Current_X = Particles[j].Position_X;
					Current_Y = Particles[j].Position_Y;
					Current_Z = Particles[j].Position_Z;

					Good_Position &= Determine_Proper_Seperation_3D(Current_X, Current_Y, Current_Z, Temp_X, Temp_Y, Temp_Z, Temp_Dist_Min);
				}

				if (!Good_Position)
				{
					Tries += 1;

					if ((Tries == 1000000) && !Good_Position)
					{
						throw std::invalid_argument("Error: Cannot place particle.");
					}
				}
			}

			Particles[i].Position_X = Temp_X;
			Particles[i].Position_Y = Temp_Y;
			Particles[i].Position_Z = Temp_Z;

			Particles[i].Velocity_X = Random_VX(Random_Generator);
			Particles[i].Velocity_Y = Random_VY(Random_Generator);
			Particles[i].Velocity_Z = Random_VZ(Random_Generator);
		}
	}
}

bool Inside_Boundaries_3D(double X_Val, double Y_Val, double Z_Val, double boundary_dist)
{
	return false;
}

void Initialize_Positions_Unit_Tests_3D()
{

	Particles[0].Position_X = 31;
	Particles[0].Position_Y = 35;

	Particles[0].Velocity_X = -0.2;
	Particles[0].Velocity_Y = 0;

	Particles[1].Position_X = 35;
	Particles[1].Position_Y = 36;

	Particles[1].Velocity_X = 0;
	Particles[1].Velocity_Y = 0.05;

	Particles[2].Position_X = 39;
	Particles[2].Position_Y = 34;

	Particles[2].Velocity_X = 0.2;
	Particles[2].Velocity_Y = -0.05;

	/*
	Particles[0].Position_X = 30;
	Particles[0].Position_Y = 40;

	Particles[0].Velocity_X = 0.2;
	Particles[0].Velocity_Y = 0;
	*/
}

bool Determine_Proper_Seperation_3D(double Particle_X, double Particle_Y, double Particle_Z, double Test_X, double Test_Y, double Test_Z, double Min_Distance)
{
	double Diff_X, Diff_Y, Diff_Z, Temp_Dist;

	Diff_X = Particle_X - Test_X;
	Diff_Y = Particle_Y - Test_Y;
	Diff_Z = Particle_Z - Test_Z;

	Temp_Dist = sqrt(pow(Diff_X, 2) + pow(Diff_Y, 2) + pow(Diff_Z, 2));

	if (Temp_Dist < Min_Distance)
	{
		return false;
	}

	return true;
}