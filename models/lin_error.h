/*
 *  lin_error.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef LIN_ERROR_H
#define LIN_ERROR_H

// Includes from models:
#include "error_transformer_node.h"
#include "error_transformer_node_impl.h"

namespace nest
{

/** @BeginDocumentation
@ingroup Neurons
@ingroup rate

Name: lin_error - Linear rate model

Description:

lin_error is an implementation of a linear rate model with
input function \f$ input(h) = g * h \f$.

The model supports connections to other rate models with either zero or
non-zero delay, and uses the secondary_event concept introduced with
the gap-junction framework.

Parameters:

The following parameters can be set in the status dictionary.
\verbatim embed:rst
===============  ======= ==================================================
 rate            real    Rate (unitless)
 tau             ms      Time constant of rate dynamics
 lambda          real    Passive decay rate
 mu              real    Mean input
 sigma           real    Noise parameter
 g               real    Gain parameter
 rectify_output  boolean Switch to restrict rate to values >= 0
===============  ======= ==================================================
\endverbatim

References:

\verbatim embed:rst
.. [1] Hahne J, Dahmen D, Schuecker J, Frommer A, Bolten M, Helias M, Diesmann
       M (2017). Integration of continuous-time dynamics in a spiking neural
       network simulator. Frontiers in Neuroinformatics, 11:34.
       DOI: https://doi.org/10.3389/fninf.2017.00034
.. [2] Hahne J, Helias M, Kunkel S, Igarashi J, Bolten M, Frommer A, Diesmann M
       (2015). A unified framework for spiking and gap-junction interactions
       in distributed neuronal network simulations.
       Frontiers Neuroinformatics, 9:22.
       DOI: https://doi.org/10.3389/fninf.2015.00022
\endverbatim

Sends: InstantaneousRateConnectionEvent, DelayedRateConnectionEvent

Receives: InstantaneousRateConnectionEvent, DelayedRateConnectionEvent,
DataLoggingRequest

Author: David Dahmen, Jan Hahne, Jannis Schuecker

SeeAlso: rate_connection_instantaneous, rate_connection_delayed
*/

class nonlinearities_lin_error
{
private:
  /** gain factor of gain function */
  double g_;
  /** linear factor in multiplicative excitatory coupling*/

public:
  /** sets default parameters */
  nonlinearities_lin_error()
    : g_( 1.0 )
  {
  }

  void get( DictionaryDatum& ) const; //!< Store current values in dictionary
  void set( const DictionaryDatum& ); //!< Set values from dicitonary

  double input( double h );               // non-linearity on input
};

inline double
nonlinearities_lin_error::input( double h )
{
  return g_ * h;
}

typedef error_transformer_node< nest::nonlinearities_lin_error >
  rate_transformer_error;

template <>
void RecordablesMap< rate_transformer_error >::create();


} // namespace nest


#endif /* #ifndef LIN_ERROR_H */
