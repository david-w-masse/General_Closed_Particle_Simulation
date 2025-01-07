// Particle_Simulation.h : Include file for standard system include files,
// or project specific include files.

#pragma once

#include <iostream>

struct Cell_Structure
{
	int Total_Divisions;
	int Total_Boundary_Divisions;
	int XY_Divisions;

	int Z_Div_4;
	int Y_Div_T;

	int Max_Per_Cell;

	int Divisions[3];

	double X_Divisor;
	double Y_Divisor;
	double Z_Divisor;

	int** Particle_Cell_List;

	int** Neighbor_List;

	int* Particles_Per_Cell;

	int* Segments_Per_Cell;

	int* Points_Per_Cell;

	int** Segments;

	int** Points_2D;

	bool Memory_Built = false;
};

struct Segment_Structure
{
	double Vertex_0_X;
	double Vertex_0_Y;

	double Vertex_1_X;
	double Vertex_1_Y;

	double Normal_X;
	double Normal_Y;
};

struct Particle_Structure
{
	double Position_X;
	double Position_Y;
	double Position_Z;

	double Past_Position_X;
	double Past_Position_Y;
	double Past_Position_Z;

	double Velocity_X;
	double Velocity_Y;
	double Velocity_Z;

	double Force_X;
	double Force_Y;
	double Force_Z;

	double Past_Force_X;
	double Past_Force_Y;
	double Past_Force_Z;

	double Size;
};

struct Boundary_Point_Structure_2D
{
	double Vertex_X;
	double Vertex_Y;

	double Normal_0_X;
	double Normal_0_Y;

	double Normal_1_X;
	double Normal_1_Y;

	bool Use_Point;
};

struct Triangle_Structure
{
	double Vertex_0_X;
	double Vertex_0_Y;
	double Vertex_0_Z;

	double Vertex_1_X;
	double Vertex_1_Y;
	double Vertex_1_Z;

	double Vertex_2_X;
	double Vertex_2_Y;
	double Vertex_2_Z;

	double Normal_0_X;
	double Normal_0_Y;
	double Normal_0_Z;

	double Normal_1_X;
	double Normal_1_Y;
	double Normal_1_Z;

};

struct Box_Structure
{
	double Center_X;
	double Center_Y;
	double Center_Z;

	double Base_Center_X;
	double Base_Center_Y;
	double Base_Center_Z;

	double Velocity_X;
	double Velocity_Y;
	double Velocity_Z;

	double Angular_Velocity_X;
	double Angular_Velocity_Y;
	double Angular_Velocity_Z;

	double Angle_X;
	double Angle_Y;
	double Angle_Z;

	double Cos_Angle_X;
	double Cos_Angle_Y;
	double Cos_Angle_Z;

	double Sin_Angle_X;
	double Sin_Angle_Y;
	double Sin_Angle_Z;

	double Inertia_X;
	double Inertia_Y;
	double Inertia_Z;

	double Torque_X;
	double Torque_Y;
	double Torque_Z;

	double Force_X;
	double Force_Y;
	double Force_Z;

	double Past_Torque_X;
	double Past_Torque_Y;
	double Past_Torque_Z;

	double Past_Force_X;
	double Past_Force_Y;
	double Past_Force_Z;

	double Mass;

	bool Movable;

	int Segment_Count;
	int Triangle_Count;

	Segment_Structure* Segments;
	Triangle_Structure* Triangles;

	Boundary_Point_Structure_2D* Points_2D;
};

bool Load_Simulation_Variables();
bool Load_Geometry_Data_2D();
bool Load_Geometry_Data_3D();
void Output_Check_File();

void Build_Memory();
void Clear_Memory();

void Assign_Particle_Data();
void Assign_Box_Data();

void Build_Neighbor_List();
void Build_Boundary_List();

void Initialize_Positions_2D(bool);
void Initialize_Positions_3D(bool);

void Initialize_Positions_Unit_Tests_2D();
void Initialize_Positions_Unit_Tests_3D();

bool Inside_Boundaries_2D(double, double, double);
bool Inside_Boundaries_3D(double, double, double, double);

void Build_Cell_List_2D();
void Build_Cell_List_3D();

void Calculate_Forces_2D();
void Calculate_Forces_3D();

void Calculate_Forces_Looped_Threaded_2D(int);
void Calculate_Forces_Looped_Threaded_3D(int);

void Calculate_Forces_Threaded_2D(int);
void Calculate_Forces_Threaded_3D(int);

void Calculate_Particle_Forces_2D(int, int);
void Calculate_Particle_Forces_3D(int, int);

void Calculate_Boundary_Force_2D(int);

bool Determine_Proper_Seperation_2D(double, double, double, double, double);
bool Determine_Proper_Seperation_3D(double, double, double, double, double, double, double);

void Debug_Info(int, int);
