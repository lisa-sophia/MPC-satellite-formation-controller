/* This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright (C) 2014-2018 by Tobias Englert, Knut Graichen, Felix Mesmer,
 * Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * Developed at the Institute of Measurement, Control, and Microtechnology,
 * Ulm University. All rights reserved.
 *
 * GRAMPC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * GRAMPC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GRAMPC. If not, see <http://www.gnu.org/licenses/>
 *
 *
 *
 *
 * This probfct file provides an interface to GRAMPC. The most general
 * formulation of the optimal control problem (OCP) that can be solved
 * by GRAMPC has the following structure
 *                                           _T
 *                                          /
 *      min    J(u,p,T;x0) = V(T,x(T),p) + / l(t,x(t),u(t),p) dt
 *   u(.),p,T                            _/
 *                                      0
 *              .
 *      s.t.   Mx(t) = f(t0+t,x(t),u(t),p), x(0) = x0
 *             g(t,x(t),u(t),p)  = 0,   gT(T,x(T),p)  = 0
 *             h(t,x(t),u(t),p) <= 0,   hT(T,x(t),p) <= 0
 *             u_min <= u(t) <= u_max
 *             p_min <=  p   <= p_max
 *             T_min <=  T   <= T_max
 *
 *
 */


#include "probfct.h"

 /* square macro */
#define POW2(a) ((a)*(a))

// W is the mean motion w = sqrt(mu/r³) = 0.001106815901, 
// where mu = G*M and r = earth's radius + orbit height = r + a
#define W (typeRNum)(sqrt((6.6743e-11 * 5.9723e24)/(pow(6378000+500000,3))))

// define setpoint here (reference array for the states) 
typeRNum xRef[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

// matrix P from the Lyapunov function V(x) = x^T*P*x
typeRNum P[6][6] = {
    {1.07786996565576,        0.0,       0.464598059590278, 262.442240474333,         0.0,       387.485312310805},
    {       0.0,       0.795061209829976,        0.0,               0.0,       292.134658374445,        0.0      },
    {0.464598059590278,        0.0,       1.56939820749626, -40.0352497803358,        0.0,       689.473441363119},
    {262.442240474333,        0.0,       -40.0352497803358, 116295.882947240,        0.0,       37586.4139960403},
    {       0.0,       292.134658374445,        0.0,               0.0,       296169.690670132,        0.0      },
    {387.485312310805,        0.0,       689.473441363119,  37586.4139960403,        0.0,       397874.483646452}
};

// matrix (A-B*K) with feedback control matrix K to estimate Lyapunov state x_dot = (A-BK)x
typeRNum AminusBK[3][6] = {
    {-1.55466767119236e-07,           0.0,         3.83184526594936e-06,	-0.00180498819166418,           0.0,           0.00302973414983929},
    {           0.0,        -1.71153947560417e-06,           0.0,                    0.0,         -0.000986377632746424,           0.0        },
    {-1.18507447567354e-06,           0.0,         -5.02689584348447e-07,	-0.00139758065297048,           0.0,           -0.00201910606594911}
};    

/** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng), 
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
{
	*Nx = 12;   // dimension of states is 12: x[1..6] = x_mpc and x[7..12] = x_lyap
	*Nu = 3;    // dimension of controls is 3: u = [u_x, u_y, u_z]
    *Np = 3;    // dimension of params is 3 the store u values of the previous time step for u_dot constraints
	*Nh = 19;   /* dimension of inequality constraints is 19: 
                 * min and max for each x_mpc = x_1 .. x_6 = 2 * 6 = 12 
                 * min and max for each   u   = u_x .. u_z = 2 * 3 = 6  
                 * Lyapunov constraint 0 >= V(x_mpc) - V(x_lyap) = 1 
                 * 12 + 6 + 1 = 19 constraints total */
	*Ng = 0;
	*NgT = 0;
	*NhT = 0;
}


/** System function f(t,x,u,p,userparam) 
    ------------------------------------ **/
void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{   
    // MPC states, calculated by the dynamics model x_dot = A*x - B*u
    out[0] = x[3];
    out[1] = x[4];
    out[2] = x[5];
    out[3] = (typeRNum)(2)*W*x[5] + u[0];
    out[4] = -POW2(W)*x[1] + u[1];
    out[5] = (typeRNum)(3)*POW2(W)*x[2] - (typeRNum)(2)*W*x[3] + u[2];

    // Lyapunov states, calculated by the dynamics model x_dot = (A-BK)x
    out[6] = x[9]-xRef[3];
    out[7] = x[10]-xRef[4];
    out[8] = x[11]-xRef[5];
    out[9]  = AminusBK[0][0]*(x[6]-xRef[0]) + AminusBK[0][2]*(x[8]-xRef[2]) + AminusBK[0][3]*(x[9]-xRef[3])  + AminusBK[0][5]*(x[11]-xRef[5]);
    out[10] = AminusBK[1][1]*(x[7]-xRef[1]) + AminusBK[1][4]*(x[10]-xRef[4]);
    out[11] = AminusBK[2][0]*(x[6]-xRef[0]) + AminusBK[2][2]*(x[8]-xRef[2]) + AminusBK[2][3]*(x[9]-xRef[3])  + AminusBK[2][5]*(x[11]-xRef[5]);
}

/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{   
    // df/dx_mpc
    out[0] = 0;
    out[1] = -POW2(W)*vec[4];
    out[2] = (typeRNum)(3)*POW2(W)*vec[5];
    out[3] = vec[0] - (typeRNum)(2)*W*vec[5];
    out[4] = vec[1];
    out[5] = vec[2] + (typeRNum)(2)*W*vec[3];
    
   // df/dx_lyap
    out[6]  = AminusBK[0][0]*vec[9] + AminusBK[2][0]*vec[11];
    out[7]  = AminusBK[1][1]*vec[10];
    out[8]  = AminusBK[0][2]*vec[9] + AminusBK[2][2]*vec[11];
    out[9]  = vec[6] + AminusBK[0][3]*vec[9] + AminusBK[2][3]*vec[11];
    out[10] = vec[7] + AminusBK[1][4]*vec[10];
    out[11] = vec[8] + AminusBK[0][5]*vec[9] + AminusBK[2][5]*vec[11];
}

/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	out[0] = vec[0];
    out[1] = vec[1];
    out[2] = vec[2];
}
/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
    out[0] = 0;
    out[1] = 0;
    out[2] = 0;
}


/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam)  = 1/2 (x^T * Q * x + u^T * R * u)
    -------------------------------------------------- **/
void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
    // weighting values of the weight matrices Q and R are stored in the "param" and can be changed in the "initData.m" file
	ctypeRNum* param = (ctypeRNum*)userparam;

	out[0] = (param[0] * POW2(x[0] - xdes[0])
            + param[1] * POW2(x[1] - xdes[1])
            + param[2] * POW2(x[2] - xdes[2])
            + param[3] * POW2(x[3] - xdes[3])
            + param[4] * POW2(x[4] - xdes[4])
            + param[5] * POW2(x[5] - xdes[5])
            + param[6] * POW2(u[0] - udes[0])
            + param[7] * POW2(u[1] - udes[1])
            + param[8] * POW2(u[2] - udes[2])) / 2;
}
/** Gradient dl/dx **/
void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
    // values of the weight matrix Q are stored in the "param[0]..[5]"
	ctypeRNum* param = (ctypeRNum*)userparam;

    // dl/dx_mpc
	out[0] = param[0] * (x[0] - xdes[0]);
	out[1] = param[1] * (x[1] - xdes[1]);
    out[2] = param[2] * (x[2] - xdes[2]);
    out[3] = param[3] * (x[3] - xdes[3]);
    out[4] = param[4] * (x[4] - xdes[4]);
    out[5] = param[5] * (x[5] - xdes[5]);

    // dl/dx_lyap
    out[6] = 0;
    out[7] = 0;
    out[8] = 0;
    out[9] = 0;
    out[10] = 0;
    out[11] = 0;
}
/** Gradient dl/du **/
void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
    // values of the weight matrix R are stored in the "param[6]..[8]"
	ctypeRNum* param = (ctypeRNum*)userparam;

	out[0] = param[6] * (u[0] - udes[0]);
    out[1] = param[7] * (u[1] - udes[1]);
    out[2] = param[8] * (u[2] - udes[2]);
}
/** Gradient dl/dp **/
void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
    out[0] = 0;
    out[1] = 0;
    out[2] = 0;
}


/** Terminal cost V(T,x(T),p,xdes,userparam) 
    ---------------------------------------- **/
void Vfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
}
/** Gradient dV/dx **/
void dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
}
/** Gradient dV/dp **/
void dVdp(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
}
/** Gradient dV/dT **/
void dVdT(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
}


/** Equality constraints g(t,x(t),u(t),p,userparam) = 0 
    --------------------------------------------------- **/
void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian dg/dx multiplied by vector vec, i.e. (dg/dx)^T*vec or vec^T*(dg/dx) **/
void dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dg/du multiplied by vector vec, i.e. (dg/du)^T*vec or vec^T*(dg/du) **/
void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dg/dp multiplied by vector vec, i.e. (dg/dp)^T*vec or vec^T*(dg/dp) **/
void dgdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Inequality constraints 0 >= h(t,x(t),u(t),p,userparam)
    ------------------------------------------------------ **/
void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{   
    /***** max and min values on the states and controls *******/
    
    // "params" stores the max and min values for the states and can be changed in "initData.m"
	ctypeRNum* param = (ctypeRNum*)userparam;

    // limitations on x
	out[0] =  param[9] - x[0];
	out[1] = -param[10] + x[0];
    
	out[2] =  param[11] - x[1];
	out[3] = -param[12] + x[1];
    
    out[4] =  param[13] - x[2];
	out[5] = -param[14] + x[2];
    
    // limitations on x_dot
    out[6] =  param[15] - x[3];
	out[7] = -param[16] + x[3];
    
    out[8] =  param[17] - x[4];
	out[9] = -param[18] + x[4];
    
    out[10] =  param[19] - x[5];
	out[11] = -param[20] + x[5];   
    
    
    //limitations on u_dot: ( |du/dt| = |u(t) - u(t-1)| <= 3.07692307692307e-06 )
    
    // in "startMPC.m" it is set vec.p(:,i) = vec.u(:,i-1); --> therefore p(t) = u(t-1)
    out[12] = -3.07692307692307e-06 - (u[0] - p[0]);
    out[13] = -3.07692307692307e-06 + (u[0] - p[0]);
    
    out[14] = -3.07692307692307e-06 - (u[1] - p[1]);
    out[15] = -3.07692307692307e-06 + (u[1] - p[1]);
    
    out[16] = -3.07692307692307e-06 - (u[2] - p[2]);
    out[17] = -3.07692307692307e-06 + (u[2] - p[2]);
    
    
    
    /***** Lyapunov constraints *******/
    
    // 0 >= V(x_mpc) - V(x_lyap)  ;   V(x) = x^T * P * x
    
    // temporary variables to store (interim) results
    typeRNum xMPC[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    typeRNum xLyap[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    typeRNum Vmpc = 0.0;
    typeRNum Vlyap = 0.0;
    typeRNum res1[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    typeRNum res2[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   
    // first initalize states
    for(int i=0; i<6; i++) {
        xMPC[i] = x[i];;
        xLyap[i] = x[i+6];
    }
     
    // multiply x^T * P
    for(int j=0; j<6; j++) {
        for(int k=0; k<6; k++) {
            res1[j] += xMPC[k]*P[k][j];
            res2[j] += xLyap[k]*P[k][j];
        }
    }
    
    // multiply (x^T * P) * x and store result
    for(int l=0; l<6; l++) {
        Vmpc += res1[l]*xMPC[l];
        Vlyap += res2[l]*xLyap[l];  
    }
 
    // Constraint V_lyap >= V_mpc   <-->   0 >= V_mpc - V_lyap
    out[18] = Vmpc - Vlyap;
}

/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{      
    // dh/dx_mpc
    out[0] = -vec[0]  + vec[1]  + 2*P[0][0]*vec[18]*x[0]  + 2*P[1][0]*vec[18]*x[1]  + 2*P[2][0]*vec[18]*x[2]  + 2*P[3][0]*vec[18]*x[3] 	+ 2*P[4][0]*vec[18]*x[4]  + 2*P[5][0]*vec[18]*x[5];
	out[1] = -vec[2]  + vec[3]  + 2*P[1][0]*vec[18]*x[0]  + 2*P[1][1]*vec[18]*x[1]  + 2*P[1][2]*vec[18]*x[2]  + 2*P[1][3]*vec[18]*x[3] 	+ 2*P[1][4]*vec[18]*x[4]  + 2*P[1][5]*vec[18]*x[5];
    out[2] = -vec[4]  + vec[5]  + 2*P[2][0]*vec[18]*x[0]  + 2*P[2][1]*vec[18]*x[1]  + 2*P[2][2]*vec[18]*x[2]  + 2*P[2][3]*vec[18]*x[3] 	+ 2*P[2][4]*vec[18]*x[4]  + 2*P[2][5]*vec[18]*x[5];
    out[3] = -vec[6]  + vec[7]  + 2*P[3][0]*vec[18]*x[0]  + 2*P[3][1]*vec[18]*x[1]  + 2*P[3][2]*vec[18]*x[2]  + 2*P[3][3]*vec[18]*x[3] 	+ 2*P[3][4]*vec[18]*x[4]  + 2*P[3][5]*vec[18]*x[5];
    out[4] = -vec[8]  + vec[9]  + 2*P[4][0]*vec[18]*x[0]  + 2*P[4][1]*vec[18]*x[1]  + 2*P[4][2]*vec[18]*x[2]  + 2*P[4][3]*vec[18]*x[3] 	+ 2*P[4][4]*vec[18]*x[4]  + 2*P[4][5]*vec[18]*x[5];
    out[5] = -vec[10] + vec[11] + 2*P[5][0]*vec[18]*x[0]  + 2*P[5][1]*vec[18]*x[1]  + 2*P[5][2]*vec[18]*x[2]  + 2*P[5][3]*vec[18]*x[3] 	+ 2*P[5][4]*vec[18]*x[4]  + 2*P[5][5]*vec[18]*x[5];
    
    // dh/dx_lyap
    out[6] =  -(2*P[0][0]*vec[18]*x[6]  + 2*P[1][0]*vec[18]*x[7] 	+ 2*P[2][0]*vec[18]*x[8]  + 2*P[3][0]*vec[18]*x[9]  + 2*P[4][0]*vec[18]*x[10]  + 2*P[5][0]*vec[18]*x[11]);
    out[7] =  -(2*P[1][0]*vec[18]*x[6]  + 2*P[1][1]*vec[18]*x[7] 	+ 2*P[1][2]*vec[18]*x[8]  + 2*P[1][3]*vec[18]*x[9]  + 2*P[1][4]*vec[18]*x[10]  + 2*P[1][5]*vec[18]*x[11]);
    out[8] =  -(2*P[2][0]*vec[18]*x[6]  + 2*P[2][1]*vec[18]*x[7] 	+ 2*P[2][2]*vec[18]*x[8]  + 2*P[2][3]*vec[18]*x[9]  + 2*P[2][4]*vec[18]*x[10]  + 2*P[2][5]*vec[18]*x[11]); 
    out[9] =  -(2*P[3][0]*vec[18]*x[6]  + 2*P[3][1]*vec[18]*x[7] 	+ 2*P[3][2]*vec[18]*x[8]  + 2*P[3][3]*vec[18]*x[9]  + 2*P[3][4]*vec[18]*x[10]  + 2*P[3][5]*vec[18]*x[11]);
    out[10] = -(2*P[4][0]*vec[18]*x[6]  + 2*P[4][1]*vec[18]*x[7] 	+ 2*P[4][2]*vec[18]*x[8]  + 2*P[4][3]*vec[18]*x[9]  + 2*P[4][4]*vec[18]*x[10]  + 2*P[4][5]*vec[18]*x[11]);
    out[11] = -(2*P[5][0]*vec[18]*x[6]  + 2*P[5][1]*vec[18]*x[7] 	+ 2*P[5][2]*vec[18]*x[8]  + 2*P[5][3]*vec[18]*x[9]  + 2*P[5][4]*vec[18]*x[10]  + 2*P[5][5]*vec[18]*x[11]);
    
}

/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
	out[0] = -vec[12] + vec[13];
    out[1] = -vec[14] + vec[15];
    out[2] = -vec[16] + vec[17];
    
}
/** Jacobian dh/dp multiplied by vector vec, i.e. (dh/dp)^T*vec or vec^T*(dg/dp) **/
void dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
    out[0] = vec[12] - vec[13];
    out[1] = vec[14] - vec[15];
    out[2] = vec[16] - vec[17];
}


/** Terminal equality constraints gT(T,x(T),p,userparam) = 0 
    -------------------------------------------------------- **/
void gTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian dgT/dx multiplied by vector vec, i.e. (dgT/dx)^T*vec or vec^T*(dgT/dx) **/
void dgTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dgT/dp multiplied by vector vec, i.e. (dgT/dp)^T*vec or vec^T*(dgT/dp) **/
void dgTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dgT/dT multiplied by vector vec, i.e. (dgT/dT)^T*vec or vec^T*(dgT/dT) **/
void dgTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Terminal inequality constraints hT(T,x(T),p,userparam) <= 0 
    ----------------------------------------------------------- **/
void hTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
void dhTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dhT/dp multiplied by vector vec, i.e. (dhT/dp)^T*vec or vec^T*(dhT/dp) **/
void dhTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
void dhTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Additional functions required for semi-implicit systems 
    M*dx/dt(t) = f(t0+t,x(t),u(t),p) using the solver RODAS 
    ------------------------------------------------------- **/
/** Jacobian df/dx in vector form (column-wise) **/
void dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian df/dx in vector form (column-wise) **/
void dfdxtrans(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian df/dt **/
void dfdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian d(dH/dx)/dt  **/
void dHdxdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *vec, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Mass matrix in vector form (column-wise, either banded or full matrix) **/
void Mfct(typeRNum *out, typeUSERPARAM *userparam)
{
}
/** Transposed mass matrix in vector form (column-wise, either banded or full matrix) **/
void Mtrans(typeRNum *out, typeUSERPARAM *userparam)
{
}
