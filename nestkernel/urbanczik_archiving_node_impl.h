/*
 *  urbanczik_archiving_node_impl.h
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

#include "urbanczik_archiving_node.h"

// Includes from nestkernel:
#include "kernel_manager.h"

// Includes from sli:
#include "dictutils.h"

namespace nest
{

// member functions for Urbanczik_Archiving_Node
template < class urbanczik_parameters >
nest::Urbanczik_Archiving_Node< urbanczik_parameters >::Urbanczik_Archiving_Node()
  : Archiving_Node()
{
}

template < class urbanczik_parameters >
nest::Urbanczik_Archiving_Node< urbanczik_parameters >::Urbanczik_Archiving_Node( const Urbanczik_Archiving_Node& n )
  : Archiving_Node( n )
{
}

template < class urbanczik_parameters >
void
nest::Urbanczik_Archiving_Node< urbanczik_parameters >::get_status( DictionaryDatum& d ) const
{
  Archiving_Node::get_status( d );
}

template < class urbanczik_parameters >
void
nest::Urbanczik_Archiving_Node< urbanczik_parameters >::set_status( const DictionaryDatum& d )
{
  Archiving_Node::set_status( d );
}

template < class urbanczik_parameters >
void
nest::Urbanczik_Archiving_Node< urbanczik_parameters >::init_urbanczik_buffers( size_t comp )
{
  hist_block_steps_ = std::max( (int) ( 2 * kernel().connection_manager.get_max_delay() + 1 ), 100 );
  hist_block_ms_ = hist_block_steps_ * Time::get_resolution().get_ms();
  block_margin_ = ( kernel().connection_manager.get_max_delay() + 1 ) *
    Time::get_resolution().get_ms();
  // the first entry in the histentry contains the time at which the block ends
  // this last time is not included in the block
  last_spike_per_synapse_[ comp - 1 ].push_back(
      histentry_extended( -1000.0 - block_margin_ + hist_block_ms_, 0.0, n_incoming_ ) );
}

template < class urbanczik_parameters >
void
nest::Urbanczik_Archiving_Node< urbanczik_parameters >::get_urbanczik_history( double t1,
  double t2,
  std::deque< histentry_extended >::iterator* start,
  std::deque< histentry_extended >::iterator* finish,
  int comp )
{
  // register spike time if it is not in the list, otherwise increase access counter.
  if ( t2 + kernel().connection_manager.get_stdp_eps() >
      ( last_spike_per_synapse_[ comp - 1 ].rbegin() )->t_ )
  {
    // make sure that the new block does not overlap with the last one
    double block_end_time = std::max(
        ( last_spike_per_synapse_[ comp -1 ].rbegin() )->t_ + hist_block_ms_,
        t2 - block_margin_ + hist_block_ms_ );
    last_spike_per_synapse_[ comp - 1 ].push_back( histentry_extended( block_end_time, 0.0, 1 ) );
  }
  else
  {
    // since the block size is larger than min delay, the spike has to be registered either in the
    // last or the second last block
    if ( t2 + kernel().connection_manager.get_stdp_eps() >
        ( ( last_spike_per_synapse_[ comp - 1 ].rbegin() )++ )->t_ )
    {
      // register spike in last block
      ( last_spike_per_synapse_[ comp - 1 ].rbegin() )->access_counter_++;
    }
    else
    {
      // register spike in second last block
      ( ( last_spike_per_synapse_[ comp - 1 ].rbegin() )++ )->access_counter_++;
    }
  }
  // search for old entry and decrease access counter and delete entry if the access counter
  // equals zero
  std::vector< histentry_extended >::iterator it_reg = std::lower_bound(
      last_spike_per_synapse_[ comp - 1 ].begin(),
      last_spike_per_synapse_[ comp - 1 ].end(),
      t1 + kernel().connection_manager.get_stdp_eps() );
  if ( it_reg == last_spike_per_synapse_[ comp - 1 ].end() )
  {
    std::cout << "found nothing, searched for:" << t1 << std::endl;
    throw std::exception();
  }
  else
  {
    it_reg->access_counter_--;
  }
  // delete old entry
  if ( it_reg->access_counter_ == 0 )
  {
    it_reg = last_spike_per_synapse_[ comp - 1 ].erase( it_reg );
  }

  // set pointer to entries of LTP hist that correspond to the times t1 and t2.
  *finish = urbanczik_history_[ comp - 1 ].end();
  if ( urbanczik_history_[ comp - 1 ].empty() )
  {
    *start = *finish;
    return;
  }
  else
  {
    // compute the position of the pointers that point to the successor of the hist entries with
    // times t1 and t2. This is straight forward because there are no time steps missing in the
    // urbanczik history. We just have to take care that *start points at least to hist.begin() and
    // *finish at most to hist.end().
    double t_first = urbanczik_history_[ comp - 1 ].begin()->t_;
    int pos_t1 = std::max( 0,
        ( (int) std::round( ( t1 - t_first ) / Time::get_resolution().get_ms() ) ) + 1 );
    int pos_t2 = std::min( (int)( urbanczik_history_[ comp - 1 ].size() ),
        ( (int) std::round( ( t2 - t_first ) / Time::get_resolution().get_ms() ) ) + 1 );

    std::deque< histentry_extended >::iterator it_first = urbanczik_history_[ comp - 1 ].begin();
    *start = it_first + std::max( 0, pos_t1 );
    *finish = it_first + std::max( 0, pos_t2 );
  }
}

template < class urbanczik_parameters >
void
nest::Urbanczik_Archiving_Node< urbanczik_parameters >::tidy_urbanczik_history( double t1, int comp )
{
  if ( last_spike_per_synapse_[ comp - 1 ].empty() )
  {
    throw std::exception();
  }
  if ( !urbanczik_history_[ comp - 1 ].empty() )
  {
    // erase history for times smaller than the smallest last spike time.
    // search for coresponding hist entry

    /* old code, we can compute the corresponding entry directly
    t1 = std::max( -1000.0, t1 );
    std::deque< histentry_extended >::iterator it_del_upper = std::lower_bound(
        urbanczik_history_[ comp - 1 ].begin(),
        urbanczik_history_[ comp - 1 ].end(),
        ( last_spike_per_synapse_[ comp - 1 ].begin() )->t_ + kernel().connection_manager.get_stdp_eps() );
    */

    std::deque< histentry_extended >::iterator finish = urbanczik_history_[ comp - 1 ].end();
    double t_first = urbanczik_history_[ comp - 1 ].begin()->t_;
    double t_oldest_spike = ( last_spike_per_synapse_[ comp - 1 ].begin() )->t_;
    int pos_t2 = std::min( (int)( urbanczik_history_[ comp - 1 ].size() ),
        ( (int) std::round( ( t_oldest_spike - t_first ) / Time::get_resolution().get_ms() ) ) + 1 );

    finish = urbanczik_history_[ comp - 1 ].begin() + std::max( 0, pos_t2 );
    // erase entries that are no longer used
    urbanczik_history_[ comp - 1 ].erase( urbanczik_history_[ comp - 1 ].begin(), finish );
  }
}

template < class urbanczik_parameters >
void
nest::Urbanczik_Archiving_Node< urbanczik_parameters >::write_urbanczik_history( Time const& t_sp,
  double V_W,
  int n_spikes,
  int comp )
{
  const double t_ms = t_sp.get_ms();

  const double g_D = urbanczik_params->g_conn[ urbanczik_parameters::SOMA ];
  const double g_L = urbanczik_params->g_L[ urbanczik_parameters::SOMA ];
  const double E_L = urbanczik_params->E_L[ urbanczik_parameters::SOMA ];
  const double V_W_star = ( ( E_L * g_L + V_W * g_D ) / ( g_D + g_L ) );

  if ( n_incoming_ )
  {
    double dPI = ( n_spikes - urbanczik_params->phi( V_W_star ) * Time::get_resolution().get_ms() )
      * urbanczik_params->h( V_W_star );
    urbanczik_history_[ comp - 1 ].push_back( histentry_extended( t_ms, dPI, 0 ) );
  }
}

} // of namespace nest
