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
nest::Eprop_Archiving_Node::init_eprop_buffers()
{
  last_spike_per_synapse_.push_back( histentry_extended( -1000.0, 0.0, n_incoming_ ) );
}

double
nest::Eprop_Archiving_Node::get_update_interval()
{
  return update_interval_;
}


double
nest:: Eprop_Archiving_Node::get_spike_history_len() const
{
  return spike_history_.size();
}


void
nest::Eprop_Archiving_Node::get_eprop_history( double t1,
  double t2,
  std::deque< histentry_eprop >::iterator* start,
  std::deque< histentry_eprop >::iterator* finish,
  bool decrease_access_counter = true )
{
  // t1 = t_last_spike_ equals -1000.0 - dendritic delay at the beginning of the simulation. To find
  // the correct entry in last_spikes_per_synapse we use the the max function.
  t1 = std::max( -1000.0, t1 );
  if ( decrease_access_counter )
  {
    // register spike time if it is not in the list, otherwise increase access counter.
    std::vector< histentry_extended >::iterator it_reg = std::lower_bound(
        last_spike_per_synapse_.begin(),
        last_spike_per_synapse_.end(),
        t2 - kernel().connection_manager.get_stdp_eps() );
    if ( it_reg == last_spike_per_synapse_.end() ||
        fabs( t2 - it_reg->t_ ) > kernel().connection_manager.get_stdp_eps() )
    {
      last_spike_per_synapse_.insert( it_reg, histentry_extended( t2, 0.0, 1 ) );
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
        t1 - kernel().connection_manager.get_stdp_eps() );
    if ( it_reg == last_spike_per_synapse_.end() ||
        fabs( t1 - it_reg->t_ ) > kernel().connection_manager.get_stdp_eps() )
    {
      std::cout << "found nothing, searched for:" << t1 << std::endl;
      /*
      std::cout << "hist:" << std::endl;
      for ( std::vector< histentry_extended >::iterator it = last_spike_per_synapse_.begin();
          it != last_spike_per_synapse_.end(); it++)
      {
        std::cout << it->t_ << "  " << it->access_counter_ << ",  ";
      }
      std::cout << std::endl;
      */
    }
    else
    {
      it_reg->access_counter_--;
    }
    // delete old entry
    if ( decrease_access_counter && it_reg->access_counter_ == 0 )
    {
      it_reg = last_spike_per_synapse_.erase( it_reg );
    }
  }

  // set pointer to entries of LTP hist that correspond to the times t1 and t2.
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
    double t_first = eprop_history_.begin()->t_;
    int pos_t1 = std::max( 0,
        ( (int) std::round( ( t1 - t_first ) / Time::get_resolution().get_ms() ) ) + 1 );
    int pos_t2 = std::min( (int)( eprop_history_.size() ),
        ( (int) std::round( ( t2 - t_first ) / Time::get_resolution().get_ms() ) ) + 1 );

    std::deque< histentry_eprop >::iterator it_first = eprop_history_.begin();
    *start = it_first + std::max( 0, pos_t1);
    *finish = it_first + std::max( 0, pos_t2);
  }
}

void
nest::Eprop_Archiving_Node::get_spike_history( double t1,
  double t2,
  std::deque< double >::iterator* start,
  std::deque< double >::iterator* finish)
{
  t1 = std::max( -1000.0, t1 );

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
  if ( !eprop_history_.empty() )
  {
    // erase history for times smaller than the smallest last spike time.
    // search for coresponding hist entry
    t1 = std::max( -1000.0, t1 );

    std::deque< histentry_eprop >::iterator start;
    std::deque< histentry_eprop >::iterator finish;
    
    // Get part of history to which the learning signal is added
    // This increases the access counter which is undone below
    nest::Eprop_Archiving_Node::get_eprop_history(
       1000.0, ( last_spike_per_synapse_.begin() )->t_, &start, &finish, false );

    /*
    std::deque< histentry_eprop >::iterator it_del_upper = std::lower_bound(
        eprop_history_.begin(),
        eprop_history_.end(),
        ( last_spike_per_synapse_.begin() )->t_ + kernel().connection_manager.get_stdp_eps() );

    if ( finish != it_del_upper )
    {
      std::cout << "got ";
      if ( finish != eprop_history_.end() )
      {
        std::cout << "here";
        std::cout << finish->t_ << " and ";
      }
      else
      {
        std::cout << "there";
        std::cout << "end() and ";
      }
      std::cout << it_del_upper->t_ << " at: " << t1 << std::endl;
      throw std::exception();
    }
    */
    // erase entries that are no longer used
    eprop_history_.erase( eprop_history_.begin(), finish );
  }
}


void
nest::Eprop_Archiving_Node::tidy_spike_history( double t1_spk )
{
  if ( !spike_history_.empty() )
  {
    t1_spk = std::max( -1000.0, t1_spk );
    std::deque< double >::iterator start_spk;
    std::deque< double >::iterator finish_spk;
    nest::Eprop_Archiving_Node::get_spike_history(
       1000.0, ( last_spike_per_synapse_.begin() )->t_, &start_spk, &finish_spk );
    spike_history_.erase( spike_history_.begin(), finish_spk );
  }
}

void
nest::Eprop_Archiving_Node::write_readout_history( Time const& t_sp,
  double learning_signal )
{
  const double t_ms = t_sp.get_ms();
  //std::cout << "learning_signal: " << learning_signal << std::endl;

  if ( n_incoming_ )
  {
    /* old code
    // prune all entries from history which are no longer needed
    // except the penultimate one. we might still need it.
    while ( eprop_history_.size() > 1 )
    {
      if ( eprop_history_.front().access_counter_ >= n_incoming_ )
      {
        eprop_history_.pop_front();
      }
      else
      {
        break;
      }
    }
    */
    // create new entry in history
    eprop_history_.push_back( histentry_eprop( t_ms, 0.0, learning_signal, 0 ) );
  }
}

void
nest::Eprop_Archiving_Node::write_eprop_history( Time const& t_sp,
  double V_m,
  double V_th )
{
  const double t_ms = t_sp.get_ms();

  if ( n_incoming_ )
  {
    // create new entry in history
    double h = pseudo_deriv( V_m, V_th );
    eprop_history_.push_back( histentry_eprop( t_ms, h, 0.0, 0 ) );
  }
  /*
  std::cout << "learning hisotry: " << std::endl;
  for ( std::deque< histentry_eprop >::iterator runner = eprop_history_.begin();
      runner != eprop_history_.end(); runner++ )
  {
    std::cout << runner->V_m_ << " ";
  }
  std::cout << std::endl;
  */
}


void
nest::Eprop_Archiving_Node::write_spike_history( Time const& t_sp )
{
  const double t_ms = t_sp.get_ms();
  spike_history_.push_back( t_ms );
}

void
nest::Eprop_Archiving_Node::add_learning_to_hist( DelayedRateConnectionEvent& e )
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
  nest::Eprop_Archiving_Node::get_eprop_history(
     t_ms, t_ms + Time::delay_steps_to_ms(delay), &start, &finish, false );

  std::vector< unsigned int >::iterator it = e.begin();

  //std::cout << "t_ms = " << t_ms << std::endl;
  //std::cout << "learning to hist: " << std::endl;

  // The call to get_coeffvalue( it ) in this loop also advances the iterator it
  while ( start != finish && it != e.end() )
  {
    // Add learning signal and reduce access counter
    start->learning_signal_ += weight * e.get_coeffvalue( it );
    ( start->access_counter_ )--;
    //std::cout << start->t_ << ", " << start->V_m_ << ", " << start->learning_signal_ << "; ";
    start++;
  }
  //std::cout << std::endl;
}

double
nest::Eprop_Archiving_Node::pseudo_deriv( double V_m, double V_th ) const
{
  double norm_diff_threshold = 1.0 - std::fabs( ( V_m - V_th ) / V_th );
  return dampening_factor_ * ( ( norm_diff_threshold > 0.0 ) ? norm_diff_threshold : 0.0 );
}

} // of namespace nest
