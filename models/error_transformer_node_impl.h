/*
 *  error_transformer_node_impl.h
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

#ifndef ERROR_TRANSFORMER_NODE_IMPL_H
#define ERROR_TRANSFORMER_NODE_IMPL_H

#include "error_transformer_node.h"

// C++ includes:
#include <cmath> // in case we need isnan() // fabs
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"
#include "gauss_rate.h"

namespace nest
{

/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

template < class TNonlinearities >
RecordablesMap< error_transformer_node< TNonlinearities > >
  error_transformer_node< TNonlinearities >::recordablesMap_;


/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

template < class TNonlinearities >
nest::error_transformer_node< TNonlinearities >::Parameters_::Parameters_()
  : linear_summation_( true )
  , tau_m_( 10.0 )                                  // ms
  , c_m_( 250.0 )                                   // pF
  , E_L_( -70.0 )                                   // mV
  , I_e_( 0.0 )                                     // pA
  , V_min_( -std::numeric_limits< double >::max() ) // relative E_L_-55.0-E_L_

{
}

template < class TNonlinearities >
nest::error_transformer_node< TNonlinearities >::State_::State_()
  : rate_( 0.0 )
  , y0_( 0.0 )
  , y3_( 0.0 )
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

template < class TNonlinearities >
void
nest::error_transformer_node< TNonlinearities >::Parameters_::get(
  DictionaryDatum& d ) const
{
  def< bool >( d, names::linear_summation, linear_summation_ );
  def< double >( d, names::E_L, E_L_ ); // Resting potential
  def< double >( d, names::I_e, I_e_ );
  def< double >( d, names::V_min, V_min_ + E_L_ );
  def< double >( d, names::C_m, c_m_ );
  def< double >( d, names::tau_m, tau_m_ );
}

template < class TNonlinearities >
double
nest::error_transformer_node< TNonlinearities >::Parameters_::set(
  const DictionaryDatum& d )
{
  updateValue< bool >( d, names::linear_summation, linear_summation_ );

  const double ELold = E_L_;
  updateValue< double >( d, names::E_L, E_L_ );
  const double delta_EL = E_L_ - ELold;
  if ( updateValue< double >( d, names::V_min, V_min_ ) )
  {
    V_min_ -= E_L_;
  }

  updateValue< double >( d, names::I_e, I_e_ );
  updateValue< double >( d, names::C_m, c_m_ );
  updateValue< double >( d, names::tau_m, tau_m_ );

  if ( c_m_ <= 0 )
  {
    throw BadProperty( "Capacitance must be >0." );
  }

  if ( tau_m_ <= 0 )
  {
    throw BadProperty( "Membrane time constant must be > 0." );
  }
  return delta_EL;
}

template < class TNonlinearities >
void
nest::error_transformer_node< TNonlinearities >::State_::get(
  DictionaryDatum& d, const Parameters_& p ) const
{
  def< double >( d, names::rate, rate_ ); // Rate
  def< double >( d, names::V_m, y3_ + p.E_L_ ); // Membrane potential
}

template < class TNonlinearities >
void
nest::error_transformer_node< TNonlinearities >::State_::set(
  const DictionaryDatum& d, const Parameters_& p, double delta_EL)
{
  updateValue< double >( d, names::rate, rate_ ); // Rate

  if ( updateValue< double >( d, names::V_m, y3_ ) )
  {
    y3_ -= p.E_L_;
  }
  else
  {
    y3_ -= delta_EL;
  }
}

template < class TNonlinearities >
nest::error_transformer_node< TNonlinearities >::Buffers_::Buffers_(
  error_transformer_node< TNonlinearities >& n )
  : logger_( n )
{
}

template < class TNonlinearities >
nest::error_transformer_node< TNonlinearities >::Buffers_::Buffers_(
  const Buffers_&,
  error_transformer_node< TNonlinearities >& n )
  : logger_( n )
{
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

template < class TNonlinearities >
nest::error_transformer_node< TNonlinearities >::error_transformer_node()
  : Eprop_Archiving_Node()
  , S_()
  , B_( *this )
{
  recordablesMap_.create();
  Node::set_node_uses_wfr( kernel().simulation_manager.use_wfr() );
}

template < class TNonlinearities >
nest::error_transformer_node< TNonlinearities >::error_transformer_node(
  const error_transformer_node& n )
  : Eprop_Archiving_Node( n )
  , nonlinearities_( n.nonlinearities_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
  Node::set_node_uses_wfr( kernel().simulation_manager.use_wfr() );
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

template < class TNonlinearities >
void
nest::error_transformer_node< TNonlinearities >::init_state_( const Node& proto )
{
  const error_transformer_node& pr = downcast< error_transformer_node >( proto );
  S_ = pr.S_;
}

template < class TNonlinearities >
void
nest::error_transformer_node< TNonlinearities >::init_buffers_()
{
  B_.delayed_rates_.clear(); // includes resize

  // resize buffers
  const size_t buffer_size = kernel().connection_manager.get_min_delay();
  B_.instant_rates_.resize( buffer_size, 0.0 );
  B_.last_y_values.resize( buffer_size, 0.0 );

  B_.logger_.reset(); // includes resize
  B_.spikes_.clear();   // includes resize
  B_.currents_.clear(); // includes resize
  Archiving_Node::clear_history();
}

template < class TNonlinearities >
void
nest::error_transformer_node< TNonlinearities >::calibrate()
{
  B_.logger_
    .init(); // ensures initialization in case mm connected after Simulate

  const double h = Time::get_resolution().get_ms();
  V_.P33_ = std::exp( -h / P_.tau_m_ );
  V_.P30_ = 1 / P_.c_m_ * ( 1 - V_.P33_ ) * P_.tau_m_;
}

/* ----------------------------------------------------------------
 * Update and event handling functions
 */

template < class TNonlinearities >
bool
nest::error_transformer_node< TNonlinearities >::update_( Time const& origin,
  const long from,
  const long to,
  const bool called_from_wfr_update )
{
  assert(
    to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  const size_t buffer_size = kernel().connection_manager.get_min_delay();
  const double wfr_tol = kernel().simulation_manager.get_wfr_tol();
  bool wfr_tol_exceeded = false;

  // allocate memory to store rates to be sent by rate events
  std::vector< double > new_rates( buffer_size, 0.0 );
  std::vector< double > new_learning_signals( buffer_size, 0.0 );

  for ( long lag = from; lag < to; ++lag )
  {

      S_.y3_ = V_.P30_ * ( S_.y0_ + P_.I_e_ ) + V_.P33_ * S_.y3_
        + B_.spikes_.get_value( lag );

      S_.y3_ = ( S_.y3_ < P_.V_min_ ? P_.V_min_ : S_.y3_ );

    // set new input current
    S_.y0_ = B_.currents_.get_value( lag );

    // store rate
    new_rates[ lag ] = S_.rate_;
    double new_learning_signal = std::abs(S_.y3_ - S_.rate_);
    new_learning_signals [ lag ] = new_learning_signal;
    // reinitialize output rate
    S_.rate_ = 0.0;

    double delayed_rates = 0;
    if ( called_from_wfr_update )
    {
      // use get_value_wfr_update to keep values in buffer
      delayed_rates = B_.delayed_rates_.get_value_wfr_update( lag );
    }
    else
    {
      // use get_value to clear values in buffer after reading
      delayed_rates = B_.delayed_rates_.get_value( lag );
    }

    if ( P_.linear_summation_ )
    {
      S_.rate_ +=
        nonlinearities_.input( delayed_rates + B_.instant_rates_[ lag ] );

    }
    else
    {
      S_.rate_ += delayed_rates + B_.instant_rates_[ lag ];
    }

    if ( called_from_wfr_update )
    {
      // check if deviation from last iteration exceeds wfr_tol
      wfr_tol_exceeded = wfr_tol_exceeded
        or fabs( S_.rate_ - B_.last_y_values[ lag ] ) > wfr_tol;
      // update last_y_values for next wfr iteration
      B_.last_y_values[ lag ] = S_.rate_;
    }
    else
    {
      // rate logging
      B_.logger_.record_data( origin.get_steps() + lag );
    }
    write_readout_history( Time::step( origin.get_steps() + lag + 1), new_learning_signal);

  }

  if ( not called_from_wfr_update )
  {
    // Send delay-rate-neuron-event. This only happens in the final iteration
    // to avoid accumulation in the buffers of the receiving neurons.
    DelayedRateConnectionEvent drve;
    // drve.set_coeffarray( new_rates )

    drve.set_coeffarray( new_learning_signals );
    kernel().event_delivery_manager.send_secondary( *this, drve );

    // clear last_y_values
    std::vector< double >( buffer_size, 0.0 ).swap( B_.last_y_values );

    // modifiy new_rates for rate-neuron-event as proxy for next min_delay
    for ( long temp = from; temp < to; ++temp )
    {
      new_rates[ temp ] = S_.rate_;
      new_learning_signals[ temp ] = std::abs(S_.y3_ - S_.rate_);
    }
  }

  // Send rate-neuron-event
  InstantaneousRateConnectionEvent rve;
  rve.set_coeffarray( new_rates );
  kernel().event_delivery_manager.send_secondary( *this, rve );

  // Reset variables
  std::vector< double >( buffer_size, 0.0 ).swap( B_.instant_rates_ );

  return wfr_tol_exceeded;
}

template < class TNonlinearities >
bool
nest::error_transformer_node< TNonlinearities >::is_eprop_readout()
    {
        return true;
    }


template < class TNonlinearities >
void
nest::error_transformer_node< TNonlinearities >::handle(
  InstantaneousRateConnectionEvent& e )
{
  const double weight = e.get_weight();

  size_t i = 0;
  std::vector< unsigned int >::iterator it = e.begin();
  // The call to get_coeffvalue( it ) in this loop also advances the iterator it
  while ( it != e.end() )
  {
    if ( P_.linear_summation_ )
    {
      B_.instant_rates_[ i ] += weight * e.get_coeffvalue( it );
    }
    else
    {
      B_.instant_rates_[ i ] +=
        weight * nonlinearities_.input( e.get_coeffvalue( it ) );
    }
    ++i;
  }
}

template < class TNonlinearities >
void
nest::error_transformer_node< TNonlinearities >::handle(
  DelayedRateConnectionEvent& e )
{
  const double weight = e.get_weight();
  const long delay = e.get_delay_steps();

  size_t i = 0;
  std::vector< unsigned int >::iterator it = e.begin();
  // The call to get_coeffvalue( it ) in this loop also advances the iterator it
  while ( it != e.end() )
  {
    if ( P_.linear_summation_ )
    {
      B_.delayed_rates_.add_value( delay + i, weight * e.get_coeffvalue( it ) );
    }
    else
    {
      B_.delayed_rates_.add_value(
        delay + i, weight * nonlinearities_.input( e.get_coeffvalue( it ) ) );
    }
    ++i;
  }
}


template < class TNonlinearities >
void
nest::error_transformer_node< TNonlinearities >::handle( SpikeEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  // EX: We must compute the arrival time of the incoming spike
  //     explicity, since it depends on delay and offset within
  //     the update cycle.  The way it is done here works, but
  //     is clumsy and should be improved.
  B_.spikes_.add_value(
    e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
    e.get_weight() * e.get_multiplicity() );
}


template < class TNonlinearities >
void
nest::error_transformer_node< TNonlinearities >::handle( CurrentEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  const double c = e.get_current();
  const double w = e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(
    e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
    w * c );
}


template < class TNonlinearities >
void
nest::error_transformer_node< TNonlinearities >::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

} // namespace

#endif /* #ifndef ERROR_TRANSFORMER_NODE_IMPL_H */
