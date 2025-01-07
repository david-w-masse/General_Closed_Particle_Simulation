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

void Initialize_Positions_2D(bool zeroize)
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

	double Temp_X, Temp_Y;
	double Temp_X_2, Temp_Y_2;
	double Temp_X_3, Temp_Y_3;
	double Temp_X_4, Temp_Y_4;

	double Temp_V_X, Temp_V_Y;
	double Temp_V_X_2, Temp_V_Y_2;
	double Temp_V_X_3, Temp_V_Y_3;
	double Temp_V_X_4, Temp_V_Y_4;

	double Temp_Vel;

	double Temp_Angle, Temp_Scale, Temp_Parity;

	double Current_X, Current_Y;

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
		while (Particles_Remaining > 0)
		{
			option_2 = true;
			option_3 = true;
			option_4 = true;

			if (Particles_Remaining == 2)
			{
				option_2 = true;
				option_3 = false;
				option_4 = false;
			}

			if (Particles_Remaining == 3)
			{
				option_2 = false;
				option_3 = true;
				option_4 = false;
			}

			if (Particles_Remaining == 4)
			{
				option_2 = true;
				option_3 = false;
				option_4 = true;
			}

			if (Particles_Remaining == 5)
			{
				option_2 = true;
				option_3 = true;
				option_4 = false;
			}
			good_choice = false;

			while (!good_choice)
			{
				pick_choice = Random_Option(Random_Generator);

				// favor 3 particle asymmetric case
				if ((pick_choice < 2) && option_3)
				{
					cout << "placing 3 .. " << endl;

					good_choice = true;

					Good_Position = false;

					Tries = 0;

					while (!Good_Position)
					{
						Temp_Angle = Random_Angle(Random_Generator);
						Temp_Scale = Random_Scale(Random_Generator);
						Temp_Parity = Random_Parity(Random_Generator);

						Temp_Vel = Random_VY(Random_Generator) / 4.0;

						Good_Position = true;

						Temp_X = -sin(Temp_Angle);
						Temp_Y = cos(Temp_Angle);

						Temp_V_X = -sin(Temp_Angle) * Temp_Vel;
						Temp_V_Y = cos(Temp_Angle) * Temp_Vel;

						Temp_X_2 = cos(Temp_Angle) * (-4);
						Temp_Y_2 = sin(Temp_Angle) * (-4);

						Temp_V_X_2 = cos(Temp_Angle) * (-4 * Temp_Vel);
						Temp_V_Y_2 = sin(Temp_Angle) * (-4 * Temp_Vel);

						Temp_X_3 = cos(Temp_Angle) * 4 - sin(Temp_Angle) * (-1);
						Temp_Y_3 = sin(Temp_Angle) * 4 + cos(Temp_Angle) * (-1);

						Temp_V_X_3 = (cos(Temp_Angle) * 4 - sin(Temp_Angle) * (-1)) * Temp_Vel;
						Temp_V_Y_3 = (sin(Temp_Angle) * 4 + cos(Temp_Angle) * (-1)) * Temp_Vel;

						if (Temp_Parity < 0.5)
						{
							Temp_X *= -1;
							Temp_X_2 *= -1;
							Temp_X_3 *= -1;

							Temp_V_X *= -1;
							Temp_V_X_2 *= -1;
							Temp_V_X_3 *= -1;
						}

						Temp_X *= Temp_Scale;
						Temp_X_2 *= Temp_Scale;
						Temp_X_3 *= Temp_Scale;

						Temp_Y *= Temp_Scale;
						Temp_Y_2 *= Temp_Scale;
						Temp_Y_3 *= Temp_Scale;

						Temp_X += Box.Center_X;
						Temp_X_2 += Box.Center_X;
						Temp_X_3 += Box.Center_X;

						Temp_Y += Box.Center_Y;
						Temp_Y_2 += Box.Center_Y;
						Temp_Y_3 += Box.Center_Y;

						//cout << Temp_X << " " << Temp_Y << " " << Temp_X_2 << " " << Temp_Y_2 << " " << Temp_X_3 << " " << Temp_Y_3 << endl;

						Good_Position &= Inside_Boundaries_2D(Temp_X, Temp_Y, boundary_dist);
						Good_Position &= Inside_Boundaries_2D(Temp_X_2, Temp_Y_2, boundary_dist);
						Good_Position &= Inside_Boundaries_2D(Temp_X_3, Temp_Y_3, boundary_dist);

						if (Good_Position)
						{
							for (int j = 0; j < Particles_Placed; j++)
							{
								Temp_Dist_Min = (Particles[Particles_Placed].Size + Particles[j].Size) * 1.175;

								Current_X = Particles[j].Position_X;
								Current_Y = Particles[j].Position_Y;

								Good_Position &= Determine_Proper_Seperation_2D(Current_X, Current_Y, Temp_X, Temp_Y, Temp_Dist_Min);
							}

							for (int j = 0; j < Particles_Placed; j++)
							{
								Temp_Dist_Min = (Particles[Particles_Placed].Size + Particles[j].Size) * 1.175;

								Current_X = Particles[j].Position_X;
								Current_Y = Particles[j].Position_Y;

								Good_Position &= Determine_Proper_Seperation_2D(Current_X, Current_Y, Temp_X_2, Temp_Y_2, Temp_Dist_Min);
							}

							for (int j = 0; j < Particles_Placed; j++)
							{
								Temp_Dist_Min = (Particles[Particles_Placed].Size + Particles[j].Size) * 1.175;

								Current_X = Particles[j].Position_X;
								Current_Y = Particles[j].Position_Y;

								Good_Position &= Determine_Proper_Seperation_2D(Current_X, Current_Y, Temp_X_3, Temp_Y_3, Temp_Dist_Min);
							}
						}

						if (!Good_Position)
						{
							Tries += 1;

							if ((Tries == 1000000) && !Good_Position)
							{
								throw std::invalid_argument("Error: Cannot place particle.");
							}
						}
						else
						{
							cout << "particle_placed " << 3 << endl;
						}
					}

					Particles[Particles_Placed].Position_X = Temp_X;
					Particles[Particles_Placed].Position_Y = Temp_Y;

					Particles[Particles_Placed].Velocity_X = Temp_V_X;
					Particles[Particles_Placed].Velocity_Y = Temp_V_Y;

					Particles_Placed += 1;

					Particles[Particles_Placed].Position_X = Temp_X_2;
					Particles[Particles_Placed].Position_Y = Temp_Y_2;

					Particles[Particles_Placed].Velocity_X = Temp_V_X_2;
					Particles[Particles_Placed].Velocity_Y = Temp_V_Y_2;

					Particles_Placed += 1;

					Particles[Particles_Placed].Position_X = Temp_X_3;
					Particles[Particles_Placed].Position_Y = Temp_Y_3;

					Particles[Particles_Placed].Velocity_X = Temp_V_X_3;
					Particles[Particles_Placed].Velocity_Y = Temp_V_Y_3;

					Particles_Placed += 1;

					Particles_Remaining -= 3;
				}
				// 2 particles
				else if ((pick_choice >= 2) && (pick_choice < 3) && option_2)
				{
					cout << "placing 2 .. " << endl;

					good_choice = true;

					good_choice = true;

					Good_Position = false;

					Tries = 0;

					while (!Good_Position)
					{
						Temp_Angle = Random_Angle(Random_Generator);
						Temp_Scale = 0.5 + Random_Scale(Random_Generator);

						Temp_Vel = Random_VY(Random_Generator);

						Good_Position = true;

						Temp_X = -sin(Temp_Angle);
						Temp_Y = cos(Temp_Angle);

						Temp_V_X = -sin(Temp_Angle) * Temp_Vel;
						Temp_V_Y = cos(Temp_Angle) * Temp_Vel;

						Temp_X_2 = -Temp_X;
						Temp_Y_2 = -Temp_Y;

						Temp_V_X_2 = -Temp_V_X;
						Temp_V_Y_2 = -Temp_V_Y;

						Temp_X *= Temp_Scale;
						Temp_X_2 *= Temp_Scale;

						Temp_Y *= Temp_Scale;
						Temp_Y_2 *= Temp_Scale;

						Temp_X += Box.Center_X;
						Temp_X_2 += Box.Center_X;

						Temp_Y += Box.Center_Y;
						Temp_Y_2 += Box.Center_Y;

						Good_Position &= Inside_Boundaries_2D(Temp_X, Temp_Y, boundary_dist);
						Good_Position &= Inside_Boundaries_2D(Temp_X_2, Temp_Y_2, boundary_dist);

						if (Good_Position)
						{
							for (int j = 0; j < Particles_Placed; j++)
							{
								Temp_Dist_Min = (Particles[Particles_Placed].Size + Particles[j].Size) * 1.175;

								Current_X = Particles[j].Position_X;
								Current_Y = Particles[j].Position_Y;

								Good_Position &= Determine_Proper_Seperation_2D(Current_X, Current_Y, Temp_X, Temp_Y, Temp_Dist_Min);
							}

							for (int j = 0; j < Particles_Placed; j++)
							{
								Temp_Dist_Min = (Particles[Particles_Placed].Size + Particles[j].Size) * 1.175;

								Current_X = Particles[j].Position_X;
								Current_Y = Particles[j].Position_Y;

								Good_Position &= Determine_Proper_Seperation_2D(Current_X, Current_Y, Temp_X_2, Temp_Y_2, Temp_Dist_Min);
							}
						}

						if (!Good_Position)
						{
							Tries += 1;

							if ((Tries == 1000000) && !Good_Position)
							{
								throw std::invalid_argument("Error: Cannot place particle.");
							}
						}
						else
						{
							cout << "particle_placed " << 2 << endl;
						}
					}

					Particles[Particles_Placed].Position_X = Temp_X;
					Particles[Particles_Placed].Position_Y = Temp_Y;

					Particles[Particles_Placed].Velocity_X = Temp_V_X;
					Particles[Particles_Placed].Velocity_Y = Temp_V_Y;

					Particles_Placed += 1;

					Particles[Particles_Placed].Position_X = Temp_X_2;
					Particles[Particles_Placed].Position_Y = Temp_Y_2;

					Particles[Particles_Placed].Velocity_X = Temp_V_X_2;
					Particles[Particles_Placed].Velocity_Y = Temp_V_Y_2;

					Particles_Placed += 1;

					Particles_Remaining -= 2;
				}
				// 4 particles
				else if ((pick_choice >= 3) && (pick_choice < 4) && option_4)
				{
					cout << "placing 4 .. " << endl;

					good_choice = true;

					good_choice = true;

					Good_Position = false;

					Tries = 0;

					while (!Good_Position)
					{
						Temp_Angle = Random_Angle(Random_Generator);
						Temp_Scale = 0.5 + Random_Scale(Random_Generator);

						Temp_Vel = Random_VY(Random_Generator);

						Good_Position = true;

						Temp_X = -sin(Temp_Angle);
						Temp_Y = cos(Temp_Angle);

						Temp_V_X = -sin(Temp_Angle) * Temp_Vel;
						Temp_V_Y = cos(Temp_Angle) * Temp_Vel;

						Temp_X_2 = -Temp_X;
						Temp_Y_2 = -Temp_Y;

						Temp_V_X_2 = -Temp_V_X;
						Temp_V_Y_2 = -Temp_V_Y;

						Temp_X_3 = Temp_X;
						Temp_Y_3 = -Temp_Y;

						Temp_V_X_3 = Temp_V_X;
						Temp_V_Y_3 = -Temp_V_Y;

						Temp_X_4 = -Temp_X;
						Temp_Y_4 = Temp_Y;

						Temp_V_X_4 = -Temp_V_X;
						Temp_V_Y_4 = Temp_V_Y;

						Temp_X *= Temp_Scale;
						Temp_X_2 *= Temp_Scale;
						Temp_X_3 *= Temp_Scale;
						Temp_X_4 *= Temp_Scale;

						Temp_Y *= Temp_Scale;
						Temp_Y_2 *= Temp_Scale;
						Temp_Y_3 *= Temp_Scale;
						Temp_Y_4 *= Temp_Scale;

						Temp_X += Box.Center_X;
						Temp_X_2 += Box.Center_X;
						Temp_X_3 += Box.Center_X;
						Temp_X_4 += Box.Center_X;

						Temp_Y += Box.Center_Y;
						Temp_Y_2 += Box.Center_Y;
						Temp_Y_3 += Box.Center_Y;
						Temp_Y_4 += Box.Center_Y;

						//cout << Temp_X << " " << Temp_Y << " " << Temp_X_2 << " " << Temp_Y_2 << " " << Temp_X_3 << " " << Temp_Y_3 << " " << Temp_X_4 << " " << Temp_Y_4 << endl;

						Good_Position &= Inside_Boundaries_2D(Temp_X, Temp_Y, boundary_dist);
						Good_Position &= Inside_Boundaries_2D(Temp_X_2, Temp_Y_2, boundary_dist);
						Good_Position &= Inside_Boundaries_2D(Temp_X_3, Temp_Y_3, boundary_dist);
						Good_Position &= Inside_Boundaries_2D(Temp_X_4, Temp_Y_4, boundary_dist);

						if (Good_Position)
						{
							for (int j = 0; j < Particles_Placed; j++)
							{
								Temp_Dist_Min = (Particles[Particles_Placed].Size + Particles[j].Size) * 1.175;

								Current_X = Particles[j].Position_X;
								Current_Y = Particles[j].Position_Y;

								Good_Position &= Determine_Proper_Seperation_2D(Current_X, Current_Y, Temp_X, Temp_Y, Temp_Dist_Min);
							}

							for (int j = 0; j < Particles_Placed; j++)
							{
								Temp_Dist_Min = (Particles[Particles_Placed].Size + Particles[j].Size) * 1.175;

								Current_X = Particles[j].Position_X;
								Current_Y = Particles[j].Position_Y;

								Good_Position &= Determine_Proper_Seperation_2D(Current_X, Current_Y, Temp_X_2, Temp_Y_2, Temp_Dist_Min);
							}

							for (int j = 0; j < Particles_Placed; j++)
							{
								Temp_Dist_Min = (Particles[Particles_Placed].Size + Particles[j].Size) * 1.175;

								Current_X = Particles[j].Position_X;
								Current_Y = Particles[j].Position_Y;

								Good_Position &= Determine_Proper_Seperation_2D(Current_X, Current_Y, Temp_X_3, Temp_Y_3, Temp_Dist_Min);
							}

							for (int j = 0; j < Particles_Placed; j++)
							{
								Temp_Dist_Min = (Particles[Particles_Placed].Size + Particles[j].Size) * 1.175;

								Current_X = Particles[j].Position_X;
								Current_Y = Particles[j].Position_Y;

								Good_Position &= Determine_Proper_Seperation_2D(Current_X, Current_Y, Temp_X_4, Temp_Y_4, Temp_Dist_Min);
							}

							// also need to check pairwise distances because small angels might put them very close together
							Temp_Dist_Min = (Particles[Particles_Placed].Size + Particles[Particles_Placed].Size) * 1.175;

							Good_Position &= Determine_Proper_Seperation_2D(Temp_X, Temp_Y, Temp_X_4, Temp_Y_4, Temp_Dist_Min);
							Good_Position &= Determine_Proper_Seperation_2D(Temp_X, Temp_Y, Temp_X_3, Temp_Y_3, Temp_Dist_Min);
						}

						if (!Good_Position)
						{
							Tries += 1;

							if ((Tries == 1000000) && !Good_Position)
							{
								cout << "fail" << endl;
								int b;
								cin >> b;
									
								throw std::invalid_argument("Error: Cannot place particle.");
							}
						}
						else
						{
							cout << "particle_placed " << 4 << endl;
						}
					}

					Particles[Particles_Placed].Position_X = Temp_X;
					Particles[Particles_Placed].Position_Y = Temp_Y;

					Particles[Particles_Placed].Velocity_X = Temp_V_X;
					Particles[Particles_Placed].Velocity_Y = Temp_V_Y;

					Particles_Placed += 1;

					Particles[Particles_Placed].Position_X = Temp_X_2;
					Particles[Particles_Placed].Position_Y = Temp_Y_2;

					Particles[Particles_Placed].Velocity_X = Temp_V_X_2;
					Particles[Particles_Placed].Velocity_Y = Temp_V_Y_2;

					Particles_Placed += 1;

					Particles[Particles_Placed].Position_X = Temp_X_3;
					Particles[Particles_Placed].Position_Y = Temp_Y_3;

					Particles[Particles_Placed].Velocity_X = Temp_V_X_3;
					Particles[Particles_Placed].Velocity_Y = Temp_V_Y_3;

					Particles_Placed += 1;

					Particles[Particles_Placed].Position_X = Temp_X_4;
					Particles[Particles_Placed].Position_Y = Temp_Y_4;

					Particles[Particles_Placed].Velocity_X = Temp_V_X_4;
					Particles[Particles_Placed].Velocity_Y = Temp_V_Y_4;

					Particles_Placed += 1;

					Particles_Remaining -= 4;
				}
			}
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

				Good_Position = Inside_Boundaries_2D(Temp_X, Temp_Y, boundary_dist);

				for (int j = 0; j < i; j++)
				{
					Temp_Dist_Min = (Particles[i].Size + Particles[j].Size) * 1.175;

					Current_X = Particles[j].Position_X;
					Current_Y = Particles[j].Position_Y;

					Good_Position &= Determine_Proper_Seperation_2D(Current_X, Current_Y, Temp_X, Temp_Y, Temp_Dist_Min);
				}

				if (!Good_Position)
				{
					Tries += 1;

					if ((Tries == 1000000) && !Good_Position)
					{
						throw std::invalid_argument("Error: Cannot place particle.");
					}
				}
				else
				{
					//cout << "particle_placed " << i << endl;
				}
			}

			Particles[i].Position_X = Temp_X;
			Particles[i].Position_Y = Temp_Y;

			Particles[i].Velocity_X = Random_VX(Random_Generator);
			Particles[i].Velocity_Y = Random_VY(Random_Generator);
		}
	}
}

bool Inside_Boundaries_2D(double X_Val, double Y_Val, double boundary_dist)
{
	double t;
	
	double x_1, x_2, y_1, y_2, u_1, u_2, v_1, v_2;

	double x_int = 0;
	double y_int = 0;

	double eps = pow(10.0, -15);

	x_1 = 0;
	x_2 = Region_Size[0];
	y_1 = Y_Val;
	y_2 = Y_Val;

	int intersect_counmt = 0;

	// determine how many boundary lines are crossed to get to the proposed point, if it is an odd number, we are inside the region
	for (int i = 0; i < Box.Segment_Count; i++)
	{
		u_1 = Box.Segments[i].Vertex_0_X;
		u_2 = Box.Segments[i].Vertex_1_X;

		v_1 = Box.Segments[i].Vertex_0_Y;
		v_2 = Box.Segments[i].Vertex_1_Y;

		
		if (v_2 == v_1)
		{
			// if horizontal, just check x bounds
			if ((X_Val > u_1) && (X_Val > u_2) && (y_1 == v_1))
			{
				intersect_counmt += 1;
			}
		}
		else if (abs(v_2 - v_1) < eps)
		{
			// if nearly horizontal, return false (outside) if within bounding box of line to avoid numerical errors
			if (v_2 > v_1)
			{
				if ((y_1 >= v_1) && (y_1 < v_2))
				{
					if (u_2 >= u_1)
					{
						if ((X_Val >= u_1) && (X_Val < u_2))
						{
							return false;
						}

						if ((X_Val >= u_1) && (X_Val >= u_2))
						{
							intersect_counmt += 1;
						}
					}
					else
					{
						if ((X_Val < u_2) && (X_Val >= u_1))
						{
							return false;
						}

						if ((X_Val >= u_1) && (X_Val >= u_2))
						{
							intersect_counmt += 1;
						}
					}				
				}
			}
			else
			{
				if ((y_1 >= v_2) && (y_1 < v_1))
				{
					if (u_2 >= u_1)
					{
						if ((X_Val >= u_1) && (X_Val < u_2))
						{
							return false;
						}

						if ((X_Val >= u_1) && (X_Val >= u_2))
						{
							intersect_counmt += 1;
						}
					}
					else
					{
						if ((X_Val < u_2) && (X_Val >= u_1))
						{
							return false;
						}

						if ((X_Val >= u_1) && (X_Val >= u_2))
						{
							intersect_counmt += 1;
						}
					}
				}
			}
		}
		else
		{
			// use general line intersect algorithm
			t = -((x_1 - x_2) * (y_1 - v_1) - (y_1 - y_2) * (x_1 - u_1)) / ((x_1 - x_2) * (v_1 - v_2) - (y_1 - y_2) * (u_1 - u_2));

			if ((t >= 0) && (t < 1))
			{
				x_int = u_1 + t * (u_2 - u_1);

				if (X_Val >= x_int)
				{
					intersect_counmt += 1;
				}
			}
		}
	}

	if ((intersect_counmt % 2) == 0) { return false; }

	// Now we know we are inside the region
	// scan through line segments and points to ensure particles are not interacting with them
	x_1 = X_Val;
	y_1 = Y_Val;

	bool inside_froce_field = false;
	double dist;

	for (int i = 0; i < Box.Segment_Count; i++)
	{
		u_1 = Box.Segments[i].Vertex_0_X;
		u_2 = Box.Segments[i].Vertex_1_X;

		v_1 = Box.Segments[i].Vertex_0_Y;
		v_2 = Box.Segments[i].Vertex_1_Y;

		x_2 = x_1 + Box.Segments[i].Normal_X;
		y_2 = y_1 + Box.Segments[i].Normal_Y;

		t = -((x_1 - x_2) * (y_1 - v_1) - (y_1 - y_2) * (x_1 - u_1)) / ((x_1 - x_2) * (v_1 - v_2) - (y_1 - y_2) * (u_1 - u_2));

		if ((t >= 0) && (t < 1))
		{
			x_int = u_1 + t * (u_2 - u_1);
			y_int = v_1 + t * (v_2 - v_1);

			dist = sqrt(pow(x_int - x_1, 2) + pow(y_int - y_1, 2));

			if (dist < boundary_dist)
			{
				/*
				cout << x_1 << " " << y_1 << " " << x_int << " " << y_int << " " << dist << " " << boundary_dist << endl;

				int b;
				cin >> b;
				*/
				return false;
			}
		}
	}

	// note for now point count is same as segment count
	// we can take a short cut here
	// since we know we are inside the box at this point, and not in range of the segments, we need only check distance to used points
	for (int i = 0; i < Box.Segment_Count; i++)
	{
		if (Box.Points_2D[i].Use_Point)
		{
			u_1 = Box.Points_2D[i].Vertex_X;
			v_1 = Box.Points_2D[i].Vertex_Y;

			dist = sqrt(pow(u_1 - x_1, 2) + pow(v_1 - y_1, 2));

			if (dist < boundary_dist)
			{
				/*
				cout << x_1 << " " << y_1 << " " << x_int << " " << y_int << " " << dist << " " << boundary_dist << endl;

				int b;
				cin >> b;
				*/
				return false;
			}
		}
	}

	return true;
}

void Initialize_Positions_Unit_Tests_2D()
{

	Particles[0].Position_X = Box.Center_X - 4;
	Particles[0].Position_Y = Box.Center_Y;

	Particles[0].Velocity_X = -0.2;
	Particles[0].Velocity_Y = 0;

	Particles[1].Position_X = Box.Center_X;
	Particles[1].Position_Y = Box.Center_Y + 1;

	Particles[1].Velocity_X = 0;
	Particles[1].Velocity_Y = 0.05;

	Particles[2].Position_X = Box.Center_X + 4;
	Particles[2].Position_Y = Box.Center_Y - 1;

	Particles[2].Velocity_X = 0.2;
	Particles[2].Velocity_Y = -0.05;

	/*
	Particles[0].Position_X = 30;
	Particles[0].Position_Y = 40;

	Particles[0].Velocity_X = 0.2;
	Particles[0].Velocity_Y = 0;
	*/
}

bool Determine_Proper_Seperation_2D(double Particle_X, double Particle_Y, double Test_X, double Test_Y, double Min_Distance)
{
	double Diff_X, Diff_Y, Temp_Dist;

	Diff_X = Particle_X - Test_X;
	Diff_Y = Particle_Y - Test_Y;

	Temp_Dist = sqrt(pow((Diff_X), 2) + pow((Diff_Y), 2));

	if (Temp_Dist < Min_Distance)
	{
		return false;
	}

	return true;
}