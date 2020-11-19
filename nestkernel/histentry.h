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
 * Part of definition of ArchivingNode which is capable of
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

  double t_; //!< point in time when spike occurred (in ms)
  double dw_;
  //! how often this entry was accessed (to enable removal, once read by all
  //! neurons which need it)
  size_t access_counter_;

  friend bool operator<( const histentry_extended he, double  t );
};

// entry for eprop with two entries: one for the learning signal and one for the membrane potential
class histentry_eprop
{
public:
  histentry_eprop( double t, double V_m, double readout_signal, double normalization, double target_signal, size_t access_counter );

  double t_; //!< point in time when spike occurred (in ms)
  double V_m_;
  double readout_signal_;
  double normalization_;
  double target_signal_;
  //! how often this entry was accessed (to enable removal, once read by all
  //! neurons which need it)
  size_t access_counter_;

  friend bool operator<( const histentry_eprop he, double  t );
};

inline bool operator<( const histentry_extended he, double t )
{
    return ( he.t_ ) < t;
}

inline bool operator<( const histentry_eprop he, double t )
{
    return ( he.t_ ) < t;
}
}

#endif
