/*
 *  histentry.h
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

/**
 * \file histentry.h
 * Part of definition of Archiving_Node which is capable of
 * recording and managing a spike history.
 * \author Moritz Helias, Abigail Morrison
 * \note moved to separate file to avoid circular inclusion in node.h
 * \date april 2006
 */

#ifndef HISTENTRY_H
#define HISTENTRY_H

// Includes from nestkernel:
#include "nest_types.h"

namespace nest
{

// entry in the spiking history
class histentry
{
public:
  histentry( double t, double Kminus, double Kminus_triplet, size_t access_counter );

  double t_;              //!< point in time when spike occurred (in ms)
  double Kminus_;         //!< value of Kminus at that time
  double Kminus_triplet_; //!< value of triplet STDP Kminus at that time
  size_t access_counter_; //!< access counter to enable removal of the entry, once all neurons read it
};

// entry in the history of plasticity rules which consider additional factors
class histentry_extended
{
public:
  histentry_extended( double t, double dw, size_t access_counter );

  double t_;              //!< point in time when spike occurred (in ms)
  double dw_;             //!< value dependend on the additional factor
  size_t access_counter_; //!< access counter to enable removal of the entry, once all neurons read it
};

// entry in the history of plasticity rules which consider additional factors
class histentry_eextended
{
public:
  histentry_eextended( double t, double I1_L, double I1_s, double I2_L, double I2_s, size_t access_counter );

  double t_;              //!< point in time when spike occurred (in ms)
  double I1_L_;              //!< point in time when spike occurred (in ms)
  double I1_s_;              //!< point in time when spike occurred (in ms)
  double I2_L_;             //!< value dependend on the additional factor
  double I2_s_;             //!< value dependend on the additional factor
  size_t access_counter_; //!< access counter to enable removal of the entry, once all neurons read it
};
}

#endif
