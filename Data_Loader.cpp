// Load basic simulation data from files
#pragma once

#include <chrono>
#include <random>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <string>

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

extern int* Update_Divisions_Start;
extern int* Update_Divisions_End;

extern int* Update_Particles_Start;
extern int* Update_Particles_End;

extern volatile atomic_bool* Ready;

void Assign_Box_Data()
{
	Box.Velocity_X = 0;
	Box.Velocity_Y = 0;
	Box.Velocity_Z = 0;

	Box.Angular_Velocity_X = 0;
	Box.Angular_Velocity_Y = 0;
	Box.Angular_Velocity_Z = 0;

	Box.Angle_X = 0;
	Box.Angle_Y = 0;
	Box.Angle_Z = 0;

	Box.Torque_X = 0;
	Box.Torque_Y = 0;
	Box.Torque_Z = 0;

	Box.Force_X = 0;
	Box.Force_Y = 0;
	Box.Force_Z = 0;

	Box.Past_Torque_X = 0;
	Box.Past_Torque_Y = 0;
	Box.Past_Torque_Z = 0;

	Box.Past_Force_X = 0;
	Box.Past_Force_Y = 0;
	Box.Past_Force_Z = 0;

	Box.Cos_Angle_X = 1;
	Box.Cos_Angle_Y = 1;
	Box.Cos_Angle_Z = 1;

	Box.Sin_Angle_X = 0;
	Box.Sin_Angle_Y = 0;
	Box.Sin_Angle_Z = 0;
}

void Assign_Particle_Data()
{
	for (int i = 0; i < Number_Particles; i++)
	{
		Particles[i].Position_X = 0;
		Particles[i].Position_Y = 0;
		Particles[i].Position_Z = 0;
		Particles[i].Past_Position_X = 0;
		Particles[i].Past_Position_Y = 0;
		Particles[i].Past_Position_Z = 0;
		Particles[i].Velocity_X = 0;
		Particles[i].Velocity_Y = 0;
		Particles[i].Velocity_Z = 0;
		Particles[i].Force_X = 0;
		Particles[i].Force_Y = 0;
		Particles[i].Force_Z = 0;
		Particles[i].Past_Force_X = 0;
		Particles[i].Past_Force_Y = 0;
		Particles[i].Past_Force_Z = 0;

		Particles[i].Size = Particle_Diameter / 2.0;
	}
}

// Simulation_Data.txt should be in the same directory as executable
bool Load_Simulation_Variables()
{
	ifstream Data_File;
	Data_File.open("Simulation_Data.txt");

	string line;

	string throw_away;
	string string_data;
	string string_data_2;

	int line_count = 0;

	if (Data_File.is_open())
	{
		while (getline(Data_File, line))
		{
			istringstream iss(line);

			if (line.at(0) != '#')
			{
				getline(iss, throw_away, '\t');
				getline(iss, throw_away, '\t');
				getline(iss, string_data, '\n');

				if (line_count == 0) { Dimension = stoi(string_data); }
				if (line_count == 1) { Number_Particles = stoi(string_data); }
				if (line_count == 2) { Particle_Diameter = stod(string_data); }
				if (line_count == 3) { Epsilon_WCA = stod(string_data); }
				if (line_count == 4) { End_Time = stod(string_data); }
				if (line_count == 5) { Delta_t = stod(string_data); }
				if (line_count == 6) { Num_Threads = stoi(string_data); }

				if (line_count == 7) { Region_Size[0] = stod(string_data); }
				if (line_count == 8) { Region_Size[1] = stod(string_data); }
				if (line_count == 9) { Region_Size[2] = stod(string_data); }
				if (line_count == 10) { Cell_Grid.Divisions[0] = stoi(string_data); }
				if (line_count == 11) { Cell_Grid.Divisions[1] = stoi(string_data); }
				if (line_count == 12) { Cell_Grid.Divisions[2] = stoi(string_data); }

				if (line_count == 13) { Box.Mass = stod(string_data); }
				if (line_count == 14) { Box.Movable = (string_data == "true"); }

				line_count++;
			}
		}
	}
	else
	{
		return false;
	}

	Data_File.close();

	if (line_count != 15) { return false; }

	// count geometry file entries
	ifstream Geom_File;

	line_count = 0;

	if (Dimension == 2)
	{
		Geom_File.open("Geometry_Data_2D.txt");
	}
	else
	{
		Geom_File.open("Geometry_Data_3D.txt");
	}

	if (Geom_File.is_open())
	{
		while (getline(Geom_File, line))
		{
			istringstream iss(line);

			if (line.at(0) != '#')
			{
				getline(iss, string_data, '\t');
				getline(iss, string_data_2, '\n');

				line_count++;
			}
		}
	}
	else
	{
		return false;
	}


	Geom_File.close();

	if (line_count == 0)
	{
		return false;
	}
	else
	{
		if (Dimension == 2)
		{
			Box.Segment_Count = line_count;
		}
		else
		{
			Box.Triangle_Count = line_count;
		}
	}

	// set up basic parameters and precalculate reused values based on loaded data
	Delta_t_square = Delta_t * Delta_t;

	Sigma_Square = Particle_Diameter * Particle_Diameter;
	Sigma_Test_WCA = pow(2.0, (2 / 6.0)) * Sigma_Square;
	Sigma_Test_LJ = pow(2.5, 2) * Sigma_Square;

	Epsilon_WCA_x_24 = Epsilon_WCA * 24.0;


	// set up basic cell grid data
	if (Dimension == 2)
	{
		Cell_Grid.Y_Div_T = Cell_Grid.Divisions[1] / Num_Threads;

		Cell_Grid.X_Divisor = Region_Size[0] / Cell_Grid.Divisions[0];
		Cell_Grid.Y_Divisor = Region_Size[1] / Cell_Grid.Divisions[1];

		Cell_Grid.Divisions[2] = 0;
		Cell_Grid.Total_Divisions = Cell_Grid.Divisions[0] * Cell_Grid.Divisions[1];
		Cell_Grid.Total_Boundary_Divisions = (Cell_Grid.Divisions[0] + Cell_Grid.Divisions[1]) * 2;

		Neighbor_Count = 4;
	}
	else
	{
		Cell_Grid.Z_Div_4 = Cell_Grid.Divisions[2] / 4;

		Cell_Grid.X_Divisor = Region_Size[0] / Cell_Grid.Divisions[0];
		Cell_Grid.Y_Divisor = Region_Size[1] / Cell_Grid.Divisions[1];
		Cell_Grid.Z_Divisor = Region_Size[2] / Cell_Grid.Divisions[2];

		Cell_Grid.Total_Divisions = Cell_Grid.Divisions[0] * Cell_Grid.Divisions[1] * Cell_Grid.Divisions[2];

		Cell_Grid.Total_Boundary_Divisions = (Cell_Grid.Divisions[0] * Cell_Grid.Divisions[1]) * 2 +
			(Cell_Grid.Divisions[1] * Cell_Grid.Divisions[2]) * 2 +
			(Cell_Grid.Divisions[0] * Cell_Grid.Divisions[2]) * 2;

		Neighbor_Count = 13;
	}

	Cell_Grid.XY_Divisions = Cell_Grid.Divisions[0] * Cell_Grid.Divisions[1];

	return true;
}

bool Load_Geometry_Data_2D()
{
	ifstream Geom_File;
	Geom_File.open("Geometry_Data_2D.txt");

	string line;

	string string_data;
	string string_data_2;

	int line_count = 0;

	// Fill in starting points of the segments
	if (Geom_File.is_open())
	{
		while (getline(Geom_File, line))
		{
			istringstream iss(line);

			if (line.at(0) != '#')
			{
				getline(iss, string_data, '\t');
				getline(iss, string_data_2, '\n');

				Box.Segments[line_count].Vertex_0_X = stod(string_data);
				Box.Segments[line_count].Vertex_0_Y = stod(string_data_2);

				line_count++;
			}
		}
	}
	else
	{
		return false;
	}

	Geom_File.close();

	// End point of one segment is the starting point of the next
	for (int i = 0; i < Box.Segment_Count-1; i++)
	{
		Box.Segments[i].Vertex_1_X = Box.Segments[i+1].Vertex_0_X;
		Box.Segments[i].Vertex_1_Y = Box.Segments[i+1].Vertex_0_Y;
	}

	Box.Segments[Box.Segment_Count - 1].Vertex_1_X = Box.Segments[0].Vertex_0_X;
	Box.Segments[Box.Segment_Count - 1].Vertex_1_Y = Box.Segments[0].Vertex_0_Y;

	double vec_x = 0;
	double vec_y = 0;
	double norm_len = 0;

	// generate normals assuming travelling counterclockwise around perimeter
	for (int i = 0; i < Box.Segment_Count; i++)
	{
		vec_x = Box.Segments[i].Vertex_1_X - Box.Segments[i].Vertex_0_X;
		vec_y = Box.Segments[i].Vertex_1_Y - Box.Segments[i].Vertex_0_Y;

		norm_len = pow(vec_x * vec_x + vec_y * vec_y, 0.5);

		Box.Segments[i].Normal_X = vec_y / norm_len;
		Box.Segments[i].Normal_Y = -vec_x / norm_len;

		Box.Points_2D[i].Use_Point = false;
	}

	// add in points necessary for angles larger than 180
	// whenever cross product of two successive normals is negative
	for (int i = 0; i < Box.Segment_Count - 1; i++)
	{
		if ((Box.Segments[i].Normal_X * Box.Segments[i + 1].Normal_Y - Box.Segments[i + 1].Normal_X * Box.Segments[i].Normal_Y) < 0)
		{
			Box.Points_2D[i + 1].Vertex_X = Box.Segments[i].Vertex_1_X;
			Box.Points_2D[i + 1].Vertex_Y = Box.Segments[i].Vertex_1_Y;

			Box.Points_2D[i + 1].Normal_0_X = Box.Segments[i].Normal_X;
			Box.Points_2D[i + 1].Normal_0_Y = Box.Segments[i].Normal_Y;

			Box.Points_2D[i + 1].Normal_1_X = Box.Segments[i + 1].Normal_X;
			Box.Points_2D[i + 1].Normal_1_Y = Box.Segments[i + 1].Normal_Y;

			Box.Points_2D[i + 1].Use_Point = true;
		}
	}

	if ((Box.Segments[Box.Segment_Count - 1].Normal_X * Box.Segments[0].Normal_Y - Box.Segments[0].Normal_X * Box.Segments[Box.Segment_Count - 1].Normal_Y) < 0)
	{
		Box.Points_2D[0].Vertex_X = Box.Segments[Box.Segment_Count - 1].Vertex_1_X;
		Box.Points_2D[0].Vertex_Y = Box.Segments[Box.Segment_Count - 1].Vertex_1_Y;

		Box.Points_2D[0].Normal_0_X = Box.Segments[Box.Segment_Count - 1].Normal_X;
		Box.Points_2D[0].Normal_0_Y = Box.Segments[Box.Segment_Count - 1].Normal_Y;

		Box.Points_2D[0].Normal_1_X = Box.Segments[0].Normal_X;
		Box.Points_2D[0].Normal_1_Y = Box.Segments[0].Normal_Y;

		Box.Points_2D[0].Use_Point = true;
	}

	// now calculate the center of mass, assuming uniform density in segments
	double m = 0;
	double cm_x = 0;
	double cm_y = 0;

	double loc_cm_x = 0;
	double loc_cm_y = 0;

	double len = 0;

	for (int i = 0; i < Box.Segment_Count; i++)
	{
		loc_cm_x = (Box.Segments[i].Vertex_0_X + Box.Segments[i].Vertex_1_X) / 2.0;
		loc_cm_y = (Box.Segments[i].Vertex_0_Y + Box.Segments[i].Vertex_1_Y) / 2.0;

		vec_x = Box.Segments[i].Vertex_1_X - Box.Segments[i].Vertex_0_X;
		vec_y = Box.Segments[i].Vertex_1_Y - Box.Segments[i].Vertex_0_Y;

		len = pow(vec_x * vec_x + vec_y * vec_y, 0.5);

		m += len;

		cm_x += loc_cm_x * len;
		cm_y += loc_cm_y * len;
	}

	Box.Center_X = cm_x / m;
	Box.Center_Y = cm_y / m;
	Box.Center_Z = 0;

	Box.Base_Center_X = Box.Center_X;
	Box.Base_Center_Y = Box.Center_Y;
	Box.Base_Center_Z = Box.Center_Z;

	// Now calculate inertia
	// Thin rods always have inertia equal to d^2/12 * m + m * C^2
	// where d is length and C is distance to the rod center
	double c = 0;

	Box.Inertia_X = 0;
	Box.Inertia_Y = 0;
	Box.Inertia_Z = 0;

	for (int i = 0; i < Box.Segment_Count; i++)
	{
		loc_cm_x = (Box.Segments[i].Vertex_0_X + Box.Segments[i].Vertex_1_X) / 2.0;
		loc_cm_y = (Box.Segments[i].Vertex_0_Y + Box.Segments[i].Vertex_1_Y) / 2.0;

		vec_x = loc_cm_x - Box.Center_X;
		vec_y = loc_cm_y - Box.Center_Y;

		// keep as c^2
		c = vec_x * vec_x + vec_y * vec_y;

		vec_x = Box.Segments[i].Vertex_1_X - Box.Segments[i].Vertex_0_X;
		vec_y = Box.Segments[i].Vertex_1_Y - Box.Segments[i].Vertex_0_Y;

		len = pow(vec_x * vec_x + vec_y * vec_y, 0.5);

		Box.Inertia_Z += Box.Mass * (len / m) * (c + len * len / 12);
	}

	cout << Box.Inertia_Z << endl;

	return true;
}

bool Load_Geometry_Data_3D()
{
	return false;
}

void Build_Boundary_List()
{
	// first scan through every cell and determine how many points and segemnt are in each.
	// build the cell grid memory accordingly
	// then scan through points and lines and assign segments accordingly
	if (Dimension == 2)
	{
		bool** temp_bounds;
		temp_bounds = new bool* [Cell_Grid.Total_Divisions];

		for (int i = 0; i < Cell_Grid.Total_Divisions; i++)
		{
			temp_bounds[i] = new bool[Box.Segment_Count];
		}

		for (int i = 0; i < Cell_Grid.Total_Divisions; i++)
		{
			for (int j = 0; j < Box.Segment_Count; j++)
			{
				temp_bounds[i][j] = false;
			}
		}

		double start_x, end_x;
		double start_y, end_y;

		double current_x, current_y;

		int Div_X = Cell_Grid.Divisions[0];
		int Div_Y = Cell_Grid.Divisions[1];

		double X_Divisor = Region_Size[0] / Div_X;
		double Y_Divisor = Region_Size[1] / Div_Y;

		double divs = 0;
		int int_divs = 0;
		int current_cell;

		int index_i, index_j;

		int count = 0;

		if (Div_X >= Div_Y)
		{
			divs = Div_X * 4.0;
			int_divs = Div_X * 4;
		}
		else
		{
			divs = Div_Y * 4.0;
			int_divs = Div_Y * 4;
		}

		for (int i = 0; i < Box.Segment_Count; i++)
		{
			start_x = Box.Segments[i].Vertex_0_X;
			end_x = Box.Segments[i].Vertex_1_X;

			start_y = Box.Segments[i].Vertex_0_Y;
			end_y = Box.Segments[i].Vertex_1_Y;

			for (int j = 0; j < int_divs; j++)
			{
				current_x = start_x + ((double)j / (divs - 1)) * (end_x - start_x);
				current_y = start_y + ((double)j / (divs - 1)) * (end_y - start_y);

				index_i = (int)floor(current_x / X_Divisor);
				index_j = (int)floor(current_y / Y_Divisor);

				// put this boundary in the grid
				if (index_i == Div_X) { index_i -= 1; }
				if (index_j == Div_Y) { index_j -= 1; }

				current_cell = index_j * Div_X + index_i;

				temp_bounds[current_cell][i] = true;

				if ((index_i + 1) < Div_X)
				{
					current_cell = index_j * Div_X + index_i + 1;

					temp_bounds[current_cell][i] = true;

					if ((index_j + 1) < Div_Y)
					{
						current_cell = (index_j + 1) * Div_X + index_i + 1;

						temp_bounds[current_cell][i] = true;
					}

					if ((index_j - 1) >= 0)
					{
						current_cell = (index_j - 1) * Div_X + index_i + 1;

						temp_bounds[current_cell][i] = true;
					}
				}

				if ((index_i - 1) >= 0)
				{
					current_cell = index_j * Div_X + index_i - 1;

					temp_bounds[current_cell][i] = true;

					if ((index_j + 1) < Div_Y)
					{
						current_cell = (index_j + 1) * Div_X + index_i - 1;

						temp_bounds[current_cell][i] = true;
					}

					if ((index_j - 1) >= 0)
					{
						current_cell = (index_j - 1) * Div_X + index_i - 1;

						temp_bounds[current_cell][i] = true;
					}
				}

				if ((index_i + 2) < Div_X)
				{
					current_cell = index_j * Div_X + index_i + 2;

					temp_bounds[current_cell][i] = true;
				}

				if ((index_i - 2) >= 0)
				{
					current_cell = index_j * Div_X + index_i - 2;

					temp_bounds[current_cell][i] = true;
				}

				if ((index_j + 1) < Div_Y)
				{
					current_cell = (index_j + 1) * Div_X + index_i;

					temp_bounds[current_cell][i] = true;
				}

				if ((index_j - 1) >= 0)
				{
					current_cell = (index_j - 1) * Div_X + index_i;

					temp_bounds[current_cell][i] = true;
				}

				if ((index_j + 2) < Div_Y)
				{
					current_cell = (index_j + 2) * Div_X + index_i;

					temp_bounds[current_cell][i] = true;
				}

				if ((index_j - 2) >= 0)
				{
					current_cell = (index_j - 2) * Div_X + index_i;

					temp_bounds[current_cell][i] = true;
				}
			}
		}

		// now we have all the data for the cell grid, build the memory
		for (int i = 0; i < Cell_Grid.Total_Divisions; i++)
		{
			for (int j = 0; j < Box.Segment_Count; j++)
			{
				if (temp_bounds[i][j])
				{
					count += 1;
				}
			}

			Cell_Grid.Segments[i] = new int[count];

			Cell_Grid.Segments_Per_Cell[i] = count;

			count = 0;

			// copy segment data to grid
			for (int j = 0; j < Box.Segment_Count; j++)
			{
				if (temp_bounds[i][j])
				{
					Cell_Grid.Segments[i][count] = j;

					count += 1;
				}
			}

			count = 0;
		}

		// repeat the same process for points (Note segment count is the same as point count)
		for (int i = 0; i < Cell_Grid.Total_Divisions; i++)
		{
			for (int j = 0; j < Box.Segment_Count; j++)
			{
				temp_bounds[i][j] = false;
			}
		}

		for (int i = 0; i < Box.Segment_Count; i++)
		{
			start_x = Box.Points_2D[i].Vertex_X;
			start_y = Box.Points_2D[i].Vertex_Y;

			index_i = (int)floor(start_x / X_Divisor);
			index_j = (int)floor(start_y / Y_Divisor);

			// put this boundary in the grid
			if (index_i == Div_X) { index_i -= 1; }
			if (index_j == Div_Y) { index_j -= 1; }

			current_cell = index_j * Div_X + index_i;

			temp_bounds[current_cell][i] = Box.Points_2D[i].Use_Point;

			if ((index_i + 1) < Div_X)
			{
				current_cell = index_j * Div_X + index_i + 1;

				temp_bounds[current_cell][i] = Box.Points_2D[i].Use_Point;

				if ((index_j + 1) < Div_Y)
				{
					current_cell = (index_j + 1) * Div_X + index_i + 1;

					temp_bounds[current_cell][i] = Box.Points_2D[i].Use_Point;
				}

				if ((index_j - 1) >= 0)
				{
					current_cell = (index_j - 1) * Div_X + index_i + 1;

					temp_bounds[current_cell][i] = Box.Points_2D[i].Use_Point;
				}
			}

			if ((index_i - 1) >= 0)
			{
				current_cell = index_j * Div_X + index_i - 1;

				temp_bounds[current_cell][i] = Box.Points_2D[i].Use_Point;

				if ((index_j + 1) < Div_Y)
				{
					current_cell = (index_j + 1) * Div_X + index_i - 1;

					temp_bounds[current_cell][i] = Box.Points_2D[i].Use_Point;
				}

				if ((index_j - 1) >= 0)
				{
					current_cell = (index_j - 1) * Div_X + index_i - 1;

					temp_bounds[current_cell][i] = Box.Points_2D[i].Use_Point;
				}
			}

			if ((index_j + 1) < Div_Y)
			{
				current_cell = (index_j + 1) * Div_X + index_i;

				temp_bounds[current_cell][i] = Box.Points_2D[i].Use_Point;
			}

			if ((index_j - 1) >= 0)
			{
				current_cell = (index_j - 1) * Div_X + index_i;

				temp_bounds[current_cell][i] = Box.Points_2D[i].Use_Point;
			}
		}

		// now we have all the data for the cell grid, build the memory
		for (int i = 0; i < Cell_Grid.Total_Divisions; i++)
		{
			for (int j = 0; j < Box.Segment_Count; j++)
			{
				if (temp_bounds[i][j])
				{
					count += 1;
				}
			}

			Cell_Grid.Points_2D[i] = new int[count];

			Cell_Grid.Points_Per_Cell[i] = count;

			count = 0;

			// copy point data to grid
			for (int j = 0; j < Box.Segment_Count; j++)
			{
				if (temp_bounds[i][j])
				{
					Cell_Grid.Points_2D[i][count] = j;

					count += 1;
				}
			}

			count = 0;
		}

		Cell_Grid.Memory_Built = true;

		for (int i = 0; i < Cell_Grid.Total_Divisions; i++)
		{
			delete[] temp_bounds[i];
		}

		delete[] temp_bounds;
	}
	else
	{

	}
}

void Build_Neighbor_List()
{
	int Base_Block, Neighbor_Index;

	int Div_X = Cell_Grid.Divisions[0];
	int Div_Y = Cell_Grid.Divisions[1];
	int Div_Z = Cell_Grid.Divisions[2];
	int Div_Tot = Cell_Grid.Total_Divisions;

	int Block_Size = Div_X * Div_Y;

	if (Dimension == 2)
	{
		for (int i = 0; i < Div_X; i++)
		{
			for (int j = 0; j < Div_Y; j++)
			{
				Neighbor_Index = j * Div_X + i;

				// Upper Left Neighbor
				if ((i == 0) || (j == (Div_Y - 1)))
				{
					Cell_Grid.Neighbor_List[Neighbor_Index][3] = -1;
				}
				else
				{
					Cell_Grid.Neighbor_List[Neighbor_Index][3] = (j + 1) * Div_X + (i - 1);
				}

				// Upper Neighbor
				if (j == (Div_Y - 1))
				{
					Cell_Grid.Neighbor_List[Neighbor_Index][2] = -1;
				}
				else
				{
					Cell_Grid.Neighbor_List[Neighbor_Index][2] = (j + 1) * Div_X + i;
				}

				// Upper Right Neighbor
				if ((i == (Div_X - 1)) || (j == (Div_Y - 1)))
				{
					Cell_Grid.Neighbor_List[Neighbor_Index][1] = -1;
				}
				else
				{
					Cell_Grid.Neighbor_List[Neighbor_Index][1] = (j + 1) * Div_X + (i + 1);
				}

				// Right Neighbor
				if (i == (Div_X - 1))
				{
					Cell_Grid.Neighbor_List[Neighbor_Index][0] = -1;
				}
				else
				{
					Cell_Grid.Neighbor_List[Neighbor_Index][0] = j * Div_X + (i + 1);
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < Div_X; i++)
		{
			for (int j = 0; j < Div_Y; j++)
			{
				for (int k = 0; k < Div_Z; k++)
				{
					Base_Block = k * Block_Size;

					Neighbor_Index = Base_Block + j * Div_X + i;

					// Right Neighbor
					if (i == (Div_X - 1))
					{
						Cell_Grid.Neighbor_List[Neighbor_Index][0] = -1;
					}
					else
					{
						Cell_Grid.Neighbor_List[Neighbor_Index][0] = Base_Block + j * Div_X + (i + 1);
					}

					// Upper Right Neighbor
					if ((i == (Div_X - 1)) || (j == (Div_Y - 1)))
					{
						Cell_Grid.Neighbor_List[Neighbor_Index][1] = -1;
					}
					else
					{
						Cell_Grid.Neighbor_List[Neighbor_Index][1] = Base_Block + (j + 1) * Div_X + (i + 1);
					}

					// Upper Neighbor
					if (j == (Div_Y - 1))
					{
						Cell_Grid.Neighbor_List[Neighbor_Index][2] = -1;
					}
					else
					{
						Cell_Grid.Neighbor_List[Neighbor_Index][2] = Base_Block + (j + 1) * Div_X + i;
					}

					// Upper Left Neighbor
					if ((i == 0) || (j == (Div_Y - 1)))
					{
						Cell_Grid.Neighbor_List[Neighbor_Index][3] = -1;
					}
					else
					{
						Cell_Grid.Neighbor_List[Neighbor_Index][3] = Base_Block + (j + 1) * Div_X + (i - 1);
					}

					if (k != (Div_Z - 1))
					{
						// Forward Right Neighbor
						if (i == (Div_X - 1))
						{
							Cell_Grid.Neighbor_List[Neighbor_Index][4] = -1;
						}
						else
						{
							Cell_Grid.Neighbor_List[Neighbor_Index][4] = Base_Block + Block_Size + j * Div_X + (i + 1);
						}

						// Forward Upper Right Neighbor
						if ((i == (Div_X - 1)) || (j == (Div_Y - 1)))
						{
							Cell_Grid.Neighbor_List[Neighbor_Index][5] = -1;
						}
						else
						{
							Cell_Grid.Neighbor_List[Neighbor_Index][5] = Base_Block + Block_Size + (j + 1) * Div_X + (i + 1);
						}

						// Forward Upper Neighbor
						if (j == (Div_Y - 1))
						{
							Cell_Grid.Neighbor_List[Neighbor_Index][6] = -1;
						}
						else
						{
							Cell_Grid.Neighbor_List[Neighbor_Index][6] = Base_Block + Block_Size + (j + 1) * Div_X + i;
						}

						// Forward Upper Left Neighbor
						if ((i == 0) || (j == (Div_Y - 1)))
						{
							Cell_Grid.Neighbor_List[Neighbor_Index][7] = -1;
						}
						else
						{
							Cell_Grid.Neighbor_List[Neighbor_Index][7] = Base_Block + Block_Size + (j + 1) * Div_X + (i - 1);
						}

						// Forward Left Neighbor
						if (i == 0)
						{
							Cell_Grid.Neighbor_List[Neighbor_Index][8] = -1;
						}
						else
						{
							Cell_Grid.Neighbor_List[Neighbor_Index][8] = Base_Block + Block_Size + j * Div_X + (i - 1);
						}

						if (j == 0)
						{
							Cell_Grid.Neighbor_List[Neighbor_Index][9] = -1;
							Cell_Grid.Neighbor_List[Neighbor_Index][10] = -1;
							Cell_Grid.Neighbor_List[Neighbor_Index][11] = -1;
						}
						else
						{
							// Forward Lower Left Neighbor
							if (i == 0)
							{
								Cell_Grid.Neighbor_List[Neighbor_Index][9] = -1;
							}
							else
							{
								Cell_Grid.Neighbor_List[Neighbor_Index][9] = Base_Block + Block_Size + (j - 1) * Div_X + (i - 1);
							}

							// Forward Lower Neighbor
							Cell_Grid.Neighbor_List[Neighbor_Index][10] = Base_Block + Block_Size + (j - 1) * Div_X + i;

							// Forward Lower Right Neighbor
							if (i == (Div_X - 1))
							{
								Cell_Grid.Neighbor_List[Neighbor_Index][4] = -1;
							}
							else
							{
								Cell_Grid.Neighbor_List[Neighbor_Index][11] = Base_Block + Block_Size + (j - 1) * Div_X + i + 1;
							}
						}

						// Forward Center Neighbor
						Cell_Grid.Neighbor_List[Neighbor_Index][12] = Base_Block + Block_Size + j * Div_X + i;
					}
					else
					{
						Cell_Grid.Neighbor_List[Neighbor_Index][4] = -1;
						Cell_Grid.Neighbor_List[Neighbor_Index][5] = -1;
						Cell_Grid.Neighbor_List[Neighbor_Index][6] = -1;
						Cell_Grid.Neighbor_List[Neighbor_Index][7] = -1;
						Cell_Grid.Neighbor_List[Neighbor_Index][8] = -1;
						Cell_Grid.Neighbor_List[Neighbor_Index][9] = -1;
						Cell_Grid.Neighbor_List[Neighbor_Index][10] = -1;
						Cell_Grid.Neighbor_List[Neighbor_Index][11] = -1;
						Cell_Grid.Neighbor_List[Neighbor_Index][12] = -1;
					}
				}
			}
		}
	}
}

void Build_Memory()
{
	Particles = new Particle_Structure[Number_Particles];

	double Cell_Volume;
	double Particle_Volume;

	// calculate how many particles can fit in the smallest grid
	// depending on dimension, then only allocate that many spots (with some to spare)
	// assume particles are smallest size
	if (Dimension == 2)
	{
		Cell_Volume = (Region_Size[0] * Region_Size[1]) / (Cell_Grid.Divisions[0] * Cell_Grid.Divisions[1]);
		Particle_Volume = Particle_Diameter * Particle_Diameter * 0.75;

		Cell_Grid.Max_Per_Cell = (int)floor(Cell_Volume / Particle_Volume * 2);
	}
	else
	{
		Cell_Volume = (Region_Size[0] * Region_Size[1] * Region_Size[2]) / (Cell_Grid.Divisions[0] * Cell_Grid.Divisions[1] * Cell_Grid.Divisions[2]);
		Particle_Volume = Particle_Diameter * Particle_Diameter * Particle_Diameter * 0.5;

		Cell_Grid.Max_Per_Cell = (int)floor(Cell_Volume / Particle_Volume * 2);
	}

	// Don't allocate more than the max number of particles
	if (Cell_Grid.Max_Per_Cell > Number_Particles)
	{
		Cell_Grid.Max_Per_Cell = Number_Particles;
	}

	Cell_Grid.Particle_Cell_List = new int* [Cell_Grid.Total_Divisions];

	Cell_Grid.Neighbor_List = new int* [Cell_Grid.Total_Divisions];

	Cell_Grid.Particles_Per_Cell = new int[Cell_Grid.Total_Divisions];

	for (int j = 0; j < Cell_Grid.Total_Divisions; j++)
	{
		Cell_Grid.Particle_Cell_List[j] = new int[Cell_Grid.Max_Per_Cell];

		Cell_Grid.Neighbor_List[j] = new int[Neighbor_Count];

		Cell_Grid.Particles_Per_Cell[j] = 0;
	}

	Update_Divisions_Start = new int[Num_Threads];
	Update_Divisions_End = new int[Num_Threads];

	Update_Particles_Start = new int[Num_Threads];
	Update_Particles_End = new int[Num_Threads];

	Ready = new volatile atomic_bool[Num_Threads];

	if (Dimension == 2)
	{
		Box.Segments = new Segment_Structure[Box.Segment_Count];
		Box.Points_2D = new Boundary_Point_Structure_2D[Box.Segment_Count];

		Cell_Grid.Segments_Per_Cell = new int[Cell_Grid.Total_Divisions];
		Cell_Grid.Points_Per_Cell = new int[Cell_Grid.Total_Divisions];

		Cell_Grid.Segments = new int * [Cell_Grid.Total_Divisions];
		Cell_Grid.Points_2D = new int * [Cell_Grid.Total_Divisions];
	}
	else
	{
		Box.Triangles = new Triangle_Structure[Box.Triangle_Count];
	}
}

void Clear_Memory()
{
	delete[] Particles;
	delete[] Cell_Grid.Particles_Per_Cell;

	for (int j = 0; j < Cell_Grid.Total_Divisions; j++)
	{
		delete[] Cell_Grid.Particle_Cell_List[j];
		delete[] Cell_Grid.Neighbor_List[j];
	}

	delete[] Cell_Grid.Particle_Cell_List;
	delete[] Cell_Grid.Neighbor_List;

	delete[] Update_Divisions_Start;
	delete[] Update_Divisions_End;

	delete[] Update_Particles_Start;
	delete[] Update_Particles_End;

	delete[] Ready;

	if (Dimension == 2)
	{
		delete[] Box.Segments;
		delete[] Box.Points_2D;

		delete[] Cell_Grid.Segments_Per_Cell;
		delete[] Cell_Grid.Points_Per_Cell;

		if (Cell_Grid.Memory_Built)
		{
			for (int j = 0; j < Cell_Grid.Total_Divisions; j++)
			{
				delete[] Cell_Grid.Segments[j];
				delete[] Cell_Grid.Points_2D[j];
			}
		}

		delete[] Cell_Grid.Segments;
		delete[] Cell_Grid.Points_2D;
	}
	else
	{
		delete[] Box.Triangles;
	}
}

void Output_Check_File()
{
	ofstream Check_File;
	Check_File.open("Check.txt", ios::out | ios::trunc);

	Check_File << Cell_Grid.Total_Divisions << "\t" << Cell_Grid.Divisions[0] << "\t" << Cell_Grid.Divisions[1] << "\n";

	Check_File << Box.Segment_Count << "\n";

	Check_File << "\n";

	for (int i = 0; i < Box.Segment_Count; i++)
	{
		Check_File << Box.Segments[i].Vertex_0_X << "\t" << Box.Segments[i].Vertex_0_Y << "\t" << Box.Segments[i].Vertex_1_X << "\t" << Box.Segments[i].Vertex_1_Y << "\n";
	}
	Check_File << "\n";

	for (int i = 0; i < Cell_Grid.Total_Divisions; i++)
	{
		Check_File << Cell_Grid.Segments_Per_Cell[i] << "\t";
	}
	Check_File << "\n";

	for (int i = 0; i < Cell_Grid.Total_Divisions; i++)
	{
		if (Cell_Grid.Segments_Per_Cell[i] > 0)
		{
			for (int j = 0; j < Cell_Grid.Segments_Per_Cell[i]; j++)
			{
				Check_File << Cell_Grid.Segments[i][j] << "\t";
			}
			Check_File << "\n";
		}
	}

	for (int i = 0; i < Cell_Grid.Total_Divisions; i++)
	{
		Check_File << Cell_Grid.Points_Per_Cell[i] << "\t";
	}

	Check_File << "\n";

	for (int i = 0; i < Cell_Grid.Total_Divisions; i++)
	{
		if (Cell_Grid.Points_Per_Cell[i] > 0)
		{
			for (int j = 0; j < Cell_Grid.Points_Per_Cell[i]; j++)
			{
				Check_File << Cell_Grid.Points_2D[i][j] << "\t";
			}
			Check_File << "\n";
		}
	}

	Check_File.close();

	ofstream Particle_File;
	Particle_File.open("Check_2.txt", ios::out | ios::trunc);

	for (int i = 0; i < Number_Particles; i++)
	{
		Particle_File << Particles[i].Position_X << "\t" << Particles[i].Position_Y << "\t" << Particles[i].Position_Z << "\n";
	}

	Particle_File.close();
}