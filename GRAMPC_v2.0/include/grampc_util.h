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
 */


#ifndef GRAMPC_UTIL_H
#define GRAMPC_UTIL_H

#include "grampc_init.h"
#include "grampc_mess.h"
#include "grampc_run.h"
#include "grampc_setopt.h"
#include "probfct.h"

/* Scaling */
#define PTR_SCALE_X   grampc->rws->rwsScale
#define PTR_SCALE_ADJ grampc->rws->rwsScale + grampc->param->Nx
#define PTR_SCALE_U   grampc->rws->rwsScale + 2 * grampc->param->Nx
#define PTR_SCALE_P   grampc->rws->rwsScale + 2 * (grampc->param->Nx + grampc->param->Nu)

#define ASSIGN_X(x_,x,grampc)       unscale_states(PTR_SCALE_X, x, grampc); x_ = PTR_SCALE_X;
#define ASSIGN_ADJ(adj_,adj,grampc) unscale_adjoints(PTR_SCALE_ADJ, adj, grampc); adj_ = PTR_SCALE_ADJ;
#define ASSIGN_U(u_,u,grampc)       unscale_controls(PTR_SCALE_U, u, grampc); u_ = PTR_SCALE_U;
#define ASSIGN_P(p_,p,grampc)       unscale_parameters(PTR_SCALE_P, p, grampc); p_ = PTR_SCALE_P;

void createNumMatrix(typeRNum **cs, const size_t size);
void createIntMatrix(typeInt **cs, const size_t size);
void resizeNumMatrix(typeRNum **cs, const size_t size);
void resizeIntMatrix(typeInt **cs, const size_t size);
void resize_rwsGeneral(const typeGRAMPC *grampc);
void resize_rwsRodas(const typeGRAMPC *grampc);
void setLWorkRodas(const typeGRAMPC *grampc);

typeInt CastDvec2Intvec(typeInt** Intvec, const double* Numvec, const size_t size);
typeInt CastDvec2Numvec(typeRNum** Realvec, const double* Numvec, const size_t size);

void discretize_time(typeRNum *tvec, typeRNum T, const typeGRAMPC *grampc);

void check_ControlLimits(const typeGRAMPC* grampc);

char* IntegratorInt2Str(ctypeInt INT_Integrator);
char* LineSearchTypeInt2Str(ctypeInt INT_LineSearchType);

void lsearch_fit2(typeRNum *kfit, typeRNum *Jfit, ctypeRNum *k, ctypeRNum *J);
void interplin(typeRNum *varint, ctypeRNum *tvec, ctypeRNum *varvec, ctypeRNum tint,
	ctypeInt Nvar, ctypeInt Nvec, ctypeInt searchdir);

void unscale_states(typeRNum *out, ctypeRNum *x, const typeGRAMPC *grampc);
void unscale_adjoints(typeRNum *out, ctypeRNum *adj, const typeGRAMPC *grampc);
void unscale_controls(typeRNum *out, ctypeRNum *u, const typeGRAMPC *grampc);
void unscale_parameters(typeRNum *out, ctypeRNum *p, const typeGRAMPC *grampc);
void unscale_time(typeRNum *out, ctypeRNum T, const typeGRAMPC *grampc);

void scale_states(typeRNum *out, ctypeRNum *x, const typeGRAMPC *grampc);
void scale_adjoints(typeRNum *out, ctypeRNum *adj, const typeGRAMPC *grampc);
void scale_controls(typeRNum *out, ctypeRNum *u, const typeGRAMPC *grampc);
void scale_parameters(typeRNum *out, ctypeRNum *p, const typeGRAMPC *grampc);
void scale_time(typeRNum *out, typeRNum T, const typeGRAMPC *grampc);

void scale_constraints(typeRNum *c, ctypeRNum *cScale, ctypeInt Ncon);
void scale_cost(typeRNum *J, ctypeRNum JScale, ctypeInt Ncost);

/* C = a */
void MatSetScalar(typeRNum *C, ctypeRNum a, ctypeInt n1, ctypeInt n2);
/* C = A */
void MatCopy(typeRNum *C, ctypeRNum *A, ctypeInt n1, ctypeInt n2);
/* C = A + B */
void MatAdd(typeRNum *C, ctypeRNum *A, ctypeRNum *B, ctypeInt n1, ctypeInt n2);
/* C = A * B */
void MatMult(typeRNum *C, ctypeRNum *A, ctypeRNum *B, ctypeInt n1, ctypeInt n2, ctypeInt n3);
/* norm = ||A||^2 */
void MatNorm(typeRNum *norm, ctypeRNum *A, ctypeInt n1, ctypeInt n2);
/* norm = ||A - B||^2 */
void MatDiffNorm(typeRNum *norm, ctypeRNum *A, ctypeRNum *B, ctypeInt n1, ctypeInt n2);

/* estimates PenaltyMin on basis of the first MPC iteration */
typeInt grampc_estim_penmin(typeGRAMPC *grampc, ctypeInt rungrampc);

#endif /* GRAMPC_UTIL_H */
