/*
 *  eprop_archiving_node.h
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

#ifndef EPROP_ARCHIVING_NODE_H
#define EPROP_ARCHIVING_NODE_H

// C++ includes:
#include <deque>

// Includes from nestkernel:
#include "histentry.h"
#include "nest_time.h"
#include "nest_types.h"
#include "archiving_node.h"
#include "synaptic_element.h"

// Includes from sli:
#include "dictdatum.h"

namespace nest
{

/**
 * \class Eprop_Archiving_Node
 * a archiving node which additionally archives parameters
 * needed for the Clopath plasticity rule
 */
class Eprop_Archiving_Node : public Archiving_Node
{

public:
  /**
   * \fn Eprop_Archiving_Node()
   * Constructor.
   */
  Eprop_Archiving_Node();

  /**
   * \fn Eprop_Archiving_Node()
   * Copy Constructor.
   */
  Eprop_Archiving_Node( const Eprop_Archiving_Node& );

  /**
   * \fn void get_LTP_history(long t1, long t2,
   * std::deque<Archiver::histentry>::iterator* start,
   * std::deque<Archiver::histentry>::iterator* finish)
   * Sets pointer start (finish) to the first (last) entry in LTP_history
   * whose time argument is between t1 and t2
   */
  void get_eprop_history( double t1,
    double t2,
    std::deque< histentry_cl >::iterator* start,
    std::deque< histentry_cl >::iterator* finish );

  /**
   * \fn double get_theta_plus()
   * Returns threshold theta_plus_
   */
  double get_theta_plus() const;

protected:
  /**
   * \fn void write_eprop_history( Time const& t_sp,
   * double u, double u_bar_plus, double u_bar_minus, double u_bar_bar )
   * Writes and reads the delayed_u_bar_[plus/minus] buffers and
   * calls write_LTD_history and write_LTP_history if
   * the corresponding Heaviside functions yield 1.
   */
  void write_eprop_history( Time const& t_sp,
    double learning_signal );

  void get_status( DictionaryDatum& d ) const;
  void set_status( const DictionaryDatum& d );

private:
  std::deque< histentry_cl > eprop_history_;

  double theta_plus_;
};

inline double
Eprop_Archiving_Node::get_theta_plus() const
{
  return theta_plus_;
}

} // of namespace
#endif
