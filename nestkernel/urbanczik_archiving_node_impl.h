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
nest::Urbanczik_Archiving_Node< urbanczik_parameters >::get_urbanczik_history( double t1,
  double t2,
  std::deque< histentry_extended >::iterator* start,
  std::deque< histentry_extended >::iterator* finish,
  int comp )
{
  *finish = urbanczik_history_[ comp - 1 ].end();
  if ( urbanczik_history_[ comp - 1 ].empty() )
  {
    *start = *finish;
    return;
  }
  else
  {
    std::deque< histentry_extended >::iterator runner = urbanczik_history_[ comp - 1 ].begin();
    // To have a well defined discretization of the integral, we make sure
    // that we exclude the entry at t1 but include the one at t2 by subtracting
    // a small number so that runner->t_ is never equal to t1 or t2.
    while ( ( runner != urbanczik_history_[ comp - 1 ].end() ) && ( runner->t_ - 1.0e-6 < t1 ) )
    {
      ++runner;
    }
    *start = runner;
    while ( ( runner != urbanczik_history_[ comp - 1 ].end() ) && ( runner->t_ - 1.0e-6 < t2 ) )
    {
      ( runner->access_counter_ )++;
      ++runner;
    }
    *finish = runner;
  }
}

template < class urbanczik_parameters >
void
nest::Urbanczik_Archiving_Node< urbanczik_parameters >::get_urbanczik_value( double t_lastspike,
    double& I1_L, double& I1_s, double& I2_L, double& I2_s, int comp )
{
  if ( urbanczik_history_compr_[ comp - 1 ].empty() )
  {
    I1_L = 0.0;
    I1_s = 0.0;
    I2_L = 0.0;
    I2_s = 0.0;
  }
  else
  {
    /*
    std::cout << "compressed history: " << std::endl;
    for ( std::deque< histentry_eextended >::iterator runner = urbanczik_history_compr_[ comp - 1 ].begin();
        runner != urbanczik_history_compr_[ comp - 1 ].end(); runner++)
    {
      std::cout << runner->t_ << "  " << runner->I1_L_ << "  " << runner->I1_s_ << "; ";
    }
    std::cout << std::endl;
    */

    std::deque< histentry_eextended >::iterator runner = urbanczik_history_compr_[ comp - 1 ].begin();
    while ( std::abs( t_lastspike - runner->t_ ) > 1.0e-6
        && runner != urbanczik_history_compr_[ comp - 1 ].end())
    {
      runner++;
    }
    if ( runner == urbanczik_history_compr_[ comp - 1 ].end() )
    {
      I1_L = 0.0;
      I1_s = 0.0;
      I2_L = 0.0;
      I2_s = 0.0;
    }
    else
    {
      runner->access_counter_ -= 1;
      I1_L = runner->I1_L_;
      I1_s = runner->I1_s_;
      I2_L = runner->I2_L_;
      I2_s = runner->I2_s_;
    }
  }
}

template < class urbanczik_parameters >
void
nest::Urbanczik_Archiving_Node< urbanczik_parameters >::compress_urbanczik_history(
    double t_compr_end, double tau_Delta, int comp )
{
  /* For this procedure to work we have to assume that: 1) incoming spikes are processed in the
   * order of their time stamps and 2) that each presynaptic neuron sends at most one spike event per
   * delta_t (resolution of the simulation). */
  // t_compr_end = t_spike - dendritic_delay
  if ( n_incoming_ )
  {
    double t_last_update = -1000.0;

    // prune all entries from history which are no longer needed
    // except the penultimate one. we might still need it.
    if ( !urbanczik_history_compr_[ comp - 1 ].empty() )
    {
      // if hist is not empty, the time of the last entry is the time of the last update
      t_last_update = urbanczik_history_compr_[ comp - 1 ].rbegin()->t_;
      std::deque< histentry_eextended >::iterator runner = urbanczik_history_compr_[ comp - 1 ].begin();
      while ( runner != urbanczik_history_compr_[ comp - 1 ].end() )
      {
        if ( runner->access_counter_ == 0 )
        {
          runner = urbanczik_history_compr_[ comp - 1 ].erase( runner );
        }
        else
        {
          runner++;
        }
      }
    }

    // if this is not the first spike in this time step, it must not process the history
    if ( std::abs( t_last_update - t_compr_end ) < 1.0e-6 )
    {
      urbanczik_history_compr_[ comp - 1 ].rbegin()->access_counter_ += 1;
      //std::cout << "I'm second" << std::endl;
    }
    else
    {
      // first update history
      double const tau_L = get_tau_L( comp );
      double const tau_s = get_tau_syn_ex( comp );
      double PI_integral_L = 0.0;
      double PI_integral_s = 0.0;
      double dPI_exp_integral_L = 0.0;
      double dPI_exp_integral_s = 0.0;
      while ( ( !urbanczik_history_[ comp - 1].empty() )
          && ( urbanczik_history_[ comp - 1].begin()->t_ - t_compr_end < 1.0e-6 ) )
      {
        std::deque< histentry_extended >::iterator he = urbanczik_history_[ comp - 1].begin();

        double const t_up = he->t_;     // from t_lastspike to t_spike
        double const minus_delta_t_up = t_last_update - t_up; // from 0 to -delta t
        double const minus_t_down = t_up - t_compr_end;          // from -t_spike to 0
        double const PI_L = exp( minus_delta_t_up / tau_L ) * he->dw_;
        double const PI_s = exp( minus_delta_t_up / tau_s ) * he->dw_;
        PI_integral_L += PI_L;
        PI_integral_s += PI_s;
        dPI_exp_integral_L += exp( minus_t_down / tau_Delta ) * PI_L;
        dPI_exp_integral_s += exp( minus_t_down / tau_Delta ) * PI_s;
        urbanczik_history_[ comp - 1].pop_front();
      }
      /*
      while ( ( runner != urbanczik_history_compr_[ comp - 1 ].end() ) && ( runner->t_ < t_compr_end - 5.0*tau_x
            ) )
      {
        runner++;
        std::cout << runner->t_ << "  " << t_compr_end << "  " << tau_x << "; ";
      }
      */
      std::deque< histentry_eextended >::iterator runner = urbanczik_history_compr_[ comp - 1 ].begin();
      while ( runner != urbanczik_history_compr_[ comp - 1 ].end() )
      {
        double const decay_L = std::exp( ( runner->t_ - t_last_update ) / tau_L );
        double const decay_s = std::exp( ( runner->t_ - t_last_update ) / tau_s );
        runner->I1_L_ += decay_L * PI_integral_L;
        runner->I1_s_ += decay_s * PI_integral_s;
        runner->I2_L_ = runner->I2_L_ * exp( ( t_last_update - t_compr_end ) / tau_Delta ) +
          decay_L * dPI_exp_integral_L;
        runner->I2_s_ = runner->I2_s_ * exp( ( t_last_update - t_compr_end ) / tau_Delta ) +
          decay_s * dPI_exp_integral_s;
        runner++;
      }
      // secondly, create new entry for current spike
      urbanczik_history_compr_[ comp - 1 ].push_back( histentry_eextended( t_compr_end, 0.0, 0.0, 0.0, 0.0, 1 ) );
    }
  }
}

template < class urbanczik_parameters >
void
nest::Urbanczik_Archiving_Node< urbanczik_parameters >::write_urbanczik_history( Time const& t_sp,
  double V_W,
  double U,
  int comp )
{
  const double t_ms = t_sp.get_ms();

  const double g_D = urbanczik_params->g_conn[ urbanczik_parameters::SOMA ];
  const double g_L = urbanczik_params->g_L[ urbanczik_parameters::SOMA ];
  const double E_L = urbanczik_params->E_L[ urbanczik_parameters::SOMA ];
  const double V_W_star = ( ( E_L * g_L + V_W * g_D ) / ( g_D + g_L ) );

  if ( n_incoming_ )
  {
    // prune all entries from history which are no longer needed
    // except the penultimate one. we might still need it.
    while ( urbanczik_history_[ comp - 1 ].size() > 1 )
    {
      if ( urbanczik_history_[ comp - 1 ].front().access_counter_ >= n_incoming_ )
      {
        urbanczik_history_[ comp - 1 ].pop_front();
      }
      else
      {
        break;
      }
    }

    double dPI = ( ( urbanczik_params->phi( U ) - urbanczik_params->phi( V_W_star ) )
        * Time::get_resolution().get_ms() ) * urbanczik_params->h( V_W_star );
    urbanczik_history_[ comp - 1 ].push_back( histentry_extended( t_ms, dPI, 0 ) );
  }
}

} // of namespace nest
