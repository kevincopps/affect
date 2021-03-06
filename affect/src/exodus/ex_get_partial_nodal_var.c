/*
 * Copyright (c) 1998 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
/*****************************************************************************
*
* exgnnv - ex_get_partial_nodal_var
*
* environment - UNIX
*
* entry conditions -
*   input parameters:
*	int	exoid			exodus file id
*	int	time_step		whole time step number
*	int	nodeal_var_index	index of desired nodal variable
*       int     start_node		starting location for read
*	int	num_nodes		number of nodal points
*
* exit conditions -
*	float*	var_vals		array of nodal variable values
*
* revision history -
*
*  $Id: ne_gnnv.c,v 1.16 2008/01/25 15:47:35 gdsjaar Exp $
*
*****************************************************************************/

#include "netcdf.h"       // for nc_inq_varid, NC_NOERR, etc
#include <exodusII.h>     // for exerrval, ex_err, etc
#include <exodusII_int.h> // for EX_WARN, ex_comp_ws, etc
#include <stddef.h>       // for size_t
#include <stdio.h>
#include <sys/types.h> // for int64_t

/*
 * reads the values of a single nodal variable for a single time step from
 * the database; assume the first time step and nodal variable index is 1
 */

int ex_get_partial_nodal_var(int exoid, int time_step, int nodal_var_index, int64_t start_node,
                             int64_t num_nodes, void *var_vals)
{
  int    varid;
  int    status;
  size_t start[3], count[3];
  char   errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  /* Verify that time_step is within bounds */
  {
    int num_time_steps = ex_inquire_int(exoid, EX_INQ_TIME);
    if (time_step <= 0 || time_step > num_time_steps) {
      snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: time_step is out-of-range. Value = %d, valid "
                                       "range is 1 to %d in file id %d",
               time_step, num_time_steps, exoid);
      ex_err("ex_get_partial_nodal_var", errmsg, EX_BADPARAM);
      return (EX_FATAL);
    }
  }

  if (ex_large_model(exoid) == 0) {
    /* read values of the nodal variable */
    if ((status = nc_inq_varid(exoid, VAR_NOD_VAR, &varid)) != NC_NOERR) {
      exerrval = status;
      snprintf(errmsg, MAX_ERR_LENGTH, "Warning: could not find nodal variables in file id %d",
               exoid);
      ex_err("ex_get_partial_nodal_var", errmsg, exerrval);
      return (EX_WARN);
    }

    start[0] = --time_step;
    start[1] = --nodal_var_index;
    start[2] = --start_node;

    count[0] = 1;
    count[1] = 1;
    count[2] = num_nodes;
  }
  else {
    /* read values of the nodal variable  -- stored as separate variables... */
    /* Get the varid.... */
    if ((status = nc_inq_varid(exoid, VAR_NOD_VAR_NEW(nodal_var_index), &varid)) != NC_NOERR) {
      exerrval = status;
      snprintf(errmsg, MAX_ERR_LENGTH, "Warning: could not find nodal variable %d in file id %d",
               nodal_var_index, exoid);
      ex_err("ex_get_partial_nodal_var", errmsg, exerrval);
      return (EX_WARN);
    }

    start[0] = --time_step;
    start[1] = --start_node;

    count[0] = 1;
    count[1] = num_nodes;
  }

  if (ex_comp_ws(exoid) == 4) {
    status = nc_get_vara_float(exoid, varid, start, count, var_vals);
  }
  else {
    status = nc_get_vara_double(exoid, varid, start, count, var_vals);
  }

  if (status != NC_NOERR) {
    exerrval = status;
    snprintf(errmsg, MAX_ERR_LENGTH, "ERROR: failed to get nodal variables in file id %d", exoid);
    ex_err("ex_get_partial_nodal_var", errmsg, exerrval);
    return (EX_FATAL);
  }
  return (EX_NOERR);
}
