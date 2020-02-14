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
    std::deque< histentry_eprop >::iterator* start,
    std::deque< histentry_eprop >::iterator* finish,
    bool decrease_access_counter );

  void get_spike_history( double t1,
    double t2,
    std::deque< double >::iterator* start,
    std::deque< double >::iterator* finish);

  void tidy_eprop_history( double t1 );
  void tidy_spike_history( double t1 );

  double get_eprop_history_len() const;
  double get_spike_history_len() const;
  double get_ls_per_syn_len() const;

  //TODO: make history private again!
  std::deque< histentry_eprop > eprop_history_;
  std::deque< double > spike_history_;

protected:
  void write_readout_history( Time const& t_sp,
    double learning_signal );

  /**
   * \fn void write_eprop_history( Time const& t_sp,
   * double u, double u_bar_plus, double u_bar_minus, double u_bar_bar )
   * Writes and reads the delayed_u_bar_[plus/minus] buffers and
   * calls write_LTD_history and write_LTP_history if
   * the corresponding Heaviside functions yield 1.
   */
  void write_eprop_history( Time const& t_sp,
    double V_m,
    double V_th );

  void write_spike_history( Time const& t_sp );

  void add_learning_to_hist( DelayedRateConnectionEvent& e );

  double pseudo_deriv( double V_m, double V_th ) const;
  // TODO: V_th as variable of archiving node? Or archiving node as template class?

  void get_status( DictionaryDatum& d ) const;
  void set_status( const DictionaryDatum& d );

  void init_eprop_buffers();

  double get_update_interval();
  //TODO: propagate information from readout neuron

private:

  double dampening_factor_; // called gamma in paper
  double update_interval_;
  std::vector< histentry_extended > last_spike_per_synapse_;

};

inline double
Eprop_Archiving_Node::get_eprop_history_len() const
{
  return eprop_history_.size();
}

inline double
Eprop_Archiving_Node::get_ls_per_syn_len() const
{
  return last_spike_per_synapse_.size();
}

} // of namespace
#endif