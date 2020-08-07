/*
 *  eprop_archiving_node.cpp
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

#include "eprop_archiving_node.h"

// Includes from nestkernel:
#include "kernel_manager.h"

// Includes from sli:
#include "dictutils.h"

namespace nest
{

// member functions for Eprop_Archiving_Node

nest::Eprop_Archiving_Node::Eprop_Archiving_Node()
  : Archiving_Node()
  , dampening_factor_( 0.3 )
  , update_interval_( 100.0 )
{
}

nest::Eprop_Archiving_Node::Eprop_Archiving_Node(
  const Eprop_Archiving_Node& n )
  : Archiving_Node( n )
  , dampening_factor_( n.dampening_factor_ )
  , update_interval_( n.update_interval_ )
{
}

void
nest::Eprop_Archiving_Node::get_status( DictionaryDatum& d ) const
{
  Archiving_Node::get_status( d );

  def< double >( d, names::dampening_factor, dampening_factor_ );
  def< double >( d, names::update_interval, update_interval_ );
}

void
nest::Eprop_Archiving_Node::set_status( const DictionaryDatum& d )
{
  Archiving_Node::set_status( d );

  // We need to preserve values in case invalid values are set
  double new_dampening_factor = dampening_factor_;
  double new_update_interval = update_interval_;
  updateValue< double >( d, names::dampening_factor, new_dampening_factor );
  updateValue< double >( d, names::update_interval, new_update_interval );

  dampening_factor_ = new_dampening_factor;
  update_interval_ = new_update_interval;
}

void
nest::Eprop_Archiving_Node::init_eprop_buffers( double delay )
{
  // register first etry for every synapse. If it is already in the list increase access counter.
  std::vector< histentry_extended >::iterator it_reg = std::lower_bound(
      last_spike_per_synapse_.begin(),
      last_spike_per_synapse_.end(),
      delay - kernel().connection_manager.get_stdp_eps() );
  if ( it_reg == last_spike_per_synapse_.end() ||
      fabs( delay - it_reg->t_ ) > kernel().connection_manager.get_stdp_eps() )
  {
    last_spike_per_synapse_.insert( it_reg, histentry_extended( delay, 0.0, 1 ) );
  }
  else
  {
    it_reg->access_counter_++;
  }
}

double
nest::Eprop_Archiving_Node::get_update_interval()
{
  return update_interval_;
}

int
nest::Eprop_Archiving_Node::get_update_interval_steps()
{
  return Time( Time::ms( update_interval_ ) ).get_steps();
}

double
nest:: Eprop_Archiving_Node::get_spike_history_len() const
{
  return spike_history_.size();
}

void
nest::Eprop_Archiving_Node::print_t_ls_per_syn()
{
  std::cout << "t_ls per syn:" << std::endl;
  for ( std::vector< histentry_extended >::iterator it = last_spike_per_synapse_.begin();
      it != last_spike_per_synapse_.end(); it++)
  {
    std::cout << it->t_ << "  " << it->access_counter_ << ",  ";
  }
  std::cout << std::endl;
}

void
nest::Eprop_Archiving_Node::print_eprop_history()
{
  std::cout << "eprop hist t, h, readout, target, normalization:" << std::endl;
  std::deque< histentry_eprop >::iterator runner = eprop_history_.begin();
  while( runner != eprop_history_.end() )
  {
    std::cout << runner->t_ << " " << runner->V_m_ << " " << runner->readout_signal_
      << " " << runner->target_signal_ << " " << runner->normalization_ << "|  ";
    runner++;
  }
  std::cout << std::endl;
}

void
nest::Eprop_Archiving_Node::find_eprop_hist_entries( double t1,
  double t2,
  std::deque< histentry_eprop >::iterator* start,
  std::deque< histentry_eprop >::iterator* finish )
{
  // set pointer to entries of eprop history hist that correspond to the times t1 and t2.
  *finish = eprop_history_.end();
  if ( eprop_history_.empty() )
  {
    *start = *finish;
    return;
  }
  else
  {
    // compute the position of the pointers that point to the successor of the hist entries with
    // times t1 and t2. This is straight forward because there are no time steps missing in the
    // eprop history. We just have to take care that *start points at least to hist.begin() and
    // *finish at most to hist.end().
    // DEBUG: set pointers to one step earlier (removed + 1)
    double t_first = eprop_history_.begin()->t_;
    int pos_t1 = std::max( 0,
        ( (int) std::round( ( t1 - t_first ) / Time::get_resolution().get_ms() ) ) + 0*1 );
    int pos_t2 = std::min( (int)( eprop_history_.size() ),
        ( (int) std::round( ( t2 - t_first ) / Time::get_resolution().get_ms() ) ) + 0*1 );

    std::deque< histentry_eprop >::iterator it_first = eprop_history_.begin();
    *start = it_first + std::max( 0, pos_t1);
    *finish = it_first + std::max( 0, pos_t2);
  }
}

void
nest::Eprop_Archiving_Node::register_update( double t_lastupdate,
   double t_update )
{
  // register spike time if it is not in the list, otherwise increase access counter.
  std::vector< histentry_extended >::iterator it_reg = std::lower_bound(
      last_spike_per_synapse_.begin(),
      last_spike_per_synapse_.end(),
      t_update - kernel().connection_manager.get_stdp_eps() );
  if ( it_reg == last_spike_per_synapse_.end() ||
      fabs( t_update - it_reg->t_ ) > kernel().connection_manager.get_stdp_eps() )
  {
    last_spike_per_synapse_.insert( it_reg, histentry_extended( t_update, 0.0, 1 ) );
  }
  else
  {
    it_reg->access_counter_++;
  }
  // search for old entry and decrease access counter and delete entry if the access counter
  // equals zero
  it_reg = std::lower_bound(
      last_spike_per_synapse_.begin(),
      last_spike_per_synapse_.end(),
      t_lastupdate - kernel().connection_manager.get_stdp_eps() );
  if ( it_reg == last_spike_per_synapse_.end() ||
      fabs( t_lastupdate - it_reg->t_ ) > kernel().connection_manager.get_stdp_eps() )
  {
    std::cout << "found nothing, searched for:" << t_lastupdate << std::endl;
  }
  else
  {
    it_reg->access_counter_--;
    // delete old entry
    if ( it_reg->access_counter_ == 0 )
    {
      it_reg = last_spike_per_synapse_.erase( it_reg );
    }
  }
}

void
nest::Eprop_Archiving_Node::get_eprop_history( double t1,
  double t2,
  double t3,
  std::deque< histentry_eprop >::iterator* start,
  std::deque< histentry_eprop >::iterator* finish )
{
  register_update( t1, t3 );
  nest::Eprop_Archiving_Node::find_eprop_hist_entries( t1, t2, start, finish );
}

void
nest::Eprop_Archiving_Node::get_spike_history( double t1,
  double t2,
  std::deque< double >::iterator* start,
  std::deque< double >::iterator* finish)
{
  // set pointer to entries of LTP hist that correspond to the times t1 and t2.
  *finish = spike_history_.end();
  if ( spike_history_.empty() )
  {
    *start = *finish;
    return;
  }
  else
  {
    std::deque< double >::iterator runner1 = std::lower_bound(
        spike_history_.begin(),
        spike_history_.end(),
        t1 + kernel().connection_manager.get_stdp_eps() );
    *start = runner1;

    std::deque< double >::iterator runner2 = std::lower_bound(
        runner1,
        spike_history_.end(),
        t2 + kernel().connection_manager.get_stdp_eps() );
    *finish = runner2;
  } //else
}

void
nest::Eprop_Archiving_Node::tidy_eprop_history( double t1 )
{
  double smallest_time_to_keep = ( last_spike_per_synapse_.begin() )->t_;
  if ( !eprop_history_.empty() )
  {
    // erase history for times smaller than the smallest last spike time.
    // search for coresponding hist entry
    std::deque< histentry_eprop >::iterator start;
    std::deque< histentry_eprop >::iterator finish;
    nest::Eprop_Archiving_Node::find_eprop_hist_entries(
       0.0, smallest_time_to_keep, &start, &finish );
    // erase entries that are no longer used
    eprop_history_.erase( eprop_history_.begin(), finish );
  }
  // tidy spike history
  for ( std::deque< double >::iterator runner = spike_history_.begin();
      runner != spike_history_.end(); runner++ )
  {
    if ( *runner + 1.0e-6 < smallest_time_to_keep )
    {
      spike_history_.erase( runner );
    }
    else
    {
      break;
    }
  }
}

void
nest::Eprop_Archiving_Node::write_eprop_history( Time const& t_sp,
  double diff_V_m_V_th,
  double V_th )
{
  const double t_ms = t_sp.get_ms();

  if ( n_incoming_ )
  {
    // create new entry in history
    // DEBUG: additional factor 1 / V_th
    double h = pseudo_deriv( diff_V_m_V_th, V_th ) / V_th;
    eprop_history_.push_back( histentry_eprop( t_ms, h, 0.0, 0.0, 0.0, 0 ) );
  }
}


void
nest::Eprop_Archiving_Node::write_spike_history( Time const& t_sp )
{
  const double t_ms = t_sp.get_ms();
  spike_history_.push_back( t_ms );
}

void
nest::Eprop_Archiving_Node::add_learning_to_hist( LearningSignalConnectionEvent& e )
{
  const double weight = e.get_weight();
  const long delay = e.get_delay_steps();
  const Time stamp = e.get_stamp();

  // TODO: Do we need to sutract the resolution? Examine delays in the network.
  double t_ms = stamp.get_ms() - Time::get_resolution().get_ms();

  std::deque< histentry_eprop >::iterator start;
  std::deque< histentry_eprop >::iterator finish;

  // Get part of history to which the learning signal is added
  // This increases the access counter which is undone below
  nest::Eprop_Archiving_Node::find_eprop_hist_entries(
     t_ms, t_ms + Time::delay_steps_to_ms(delay), &start, &finish );
  std::vector< unsigned int >::iterator it = e.begin();
  /*
  double ls = e.get_coeffvalue( it );
  if ( ls != 0.0 && t_ms < 3510.0 )
  {
    std::cout << "add learning_signal at t = " << t_ms << std::endl;
  }
  */
  while ( start != finish && it != e.end() )
  {
    // Add learning signal and reduce access counter
    double regression = e.get_coeffvalue(it);
    double readout_signal = e.get_coeffvalue(it);
    double target_signal = e.get_coeffvalue(it);
    double recall = e.get_coeffvalue(it);
    if (recall == 1.)
    {
      double old_norm = start->normalization_;
      if (regression == 1.)
      {
          start->readout_signal_ += weight * readout_signal;
          start->target_signal_ += weight * target_signal;
          start->normalization_ = 1.;
      }
      else if (regression == 0.)
      {
          start->readout_signal_ += weight * std::exp(readout_signal);
          start->target_signal_ += weight * target_signal;
          start->normalization_ += std::exp(readout_signal);
      }
      //std::cout << "added normalization: " << start->t_ << ", " << stamp.get_ms() << ", " << old_norm << " -> " <<
      //  start->normalization_ << std::endl;
    }
    else
    {
       start->readout_signal_ += 0.;
       start->target_signal_ += 0.;
       start->normalization_ += 1.;
    }
    start++;
  }
}

double
nest::Eprop_Archiving_Node::pseudo_deriv( double diff_V_m_V_th, double V_th_const ) const
{
  // DEBUG: v_scaled = (Vm - adaptive_thr) / V_th,
  // where adaptive_thr is the spiking threshold including the adaptive part and
  // V_th is the constant part of the threshold. In the normal LIF neuron
  // adaptive_thr = V_th
  double norm_diff_threshold = 1.0 - std::fabs( ( diff_V_m_V_th ) / V_th_const );
  return dampening_factor_ * ( ( norm_diff_threshold > 0.0 ) ? norm_diff_threshold : 0.0 );
}

} // of namespace nest
