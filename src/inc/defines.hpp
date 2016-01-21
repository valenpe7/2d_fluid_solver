#ifndef DEFINES_H
#define DEFINES_H

// definitions for x and y directions
const int X = 0;
const int Y = 1;

// definitions of "compass rose" directions
const int D_EAST		= 0;	// direction +X
const int D_WEST		= 1;	// direction -X
const int D_NORTH		= 2;	// direction +Y
const int D_SOUTH		= 3;	// direction -Y

// definitions of boundary conditions
const int BC_FREE_SLIP	= 1;
const int BC_NO_SLIP	= 2;
const int BC_OUTFLOW	= 3;
const int BC_INFLOW		= 4;   

// definitions of flag field macros for obstacles and boundary cell definitions
#define C_B       0x0000      	// obstacle cell
#define B_N       0x0001
#define B_S       0x0002
#define B_W       0x0004
#define B_E       0x0008
#define B_NW      0x0005
#define B_SW      0x0006
#define B_NE      0x0009
#define B_SE      0x000A
#define C_F       0x0010     	// fluid cell

// Macros for Poisson equation, denoting whether there is an obstacle cell adjacent to some direction
#define eps_E  !(dat->flag[i+1][j] < C_F)
#define eps_W  !(dat->flag[i-1][j] < C_F)
#define eps_N  !(dat->flag[i][j+1] < C_F)
#define eps_S  !(dat->flag[i][j-1] < C_F)

// definitions for bounding box
const int BB_MIN_X		= 0;
const int BB_MIN_Y		= 1;
const int BB_MAX_X		= 2;
const int BB_MAX_Y		= 3;

#endif //DEFINES_H
