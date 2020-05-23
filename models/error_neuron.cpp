/*
 *  error_neuron.cpp
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

#include "error_neuron.h"

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

RecordablesMap< error_neuron >error_neuron::recordablesMap_;

template <>
void
RecordablesMap< error_neuron >::create()
{
  insert_( names::target_rate, &error_neuron::get_target_rate_ );
  insert_( names::learning_signal, &error_neuron::get_learning_signal_ );
  insert_( names::V_m, &error_neuron::get_V_m_ );
  insert_( names::len_eprop_hist, &error_neuron::get_eprop_history_len );
  insert_( names::len_ls_per_syn, &error_neuron::get_ls_per_syn_len );
}
/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

nest::error_neuron::Parameters_::Parameters_()
  : tau_m_( 10.0 )                                  // ms
  , c_m_( 250.0 )                                   // pF
  , E_L_( -70.0 )                                   // mV
  , I_e_( 0.0 )                                     // pA
  , V_min_( -std::numeric_limits< double >::max() ) // relative E_L_-55.0-E_L_
  , t_start_ls_( 0.0 )                               // ms
{
}

nest::error_neuron::State_::State_()
  : target_rate_( 0.0 )
  , learning_signal_( 0.0 )
  , y0_( 0.0 )
  , y3_( 0.0 )
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
nest::error_neuron::Parameters_::get(
  DictionaryDatum& d ) const
{
  def< double >( d, names::E_L, E_L_ ); // Resting potential
  def< double >( d, names::I_e, I_e_ );
  def< double >( d, names::V_min, V_min_ + E_L_ );
  def< double >( d, names::C_m, c_m_ );
  def< double >( d, names::tau_m, tau_m_ );
  def< double >( d, names::start, t_start_ls_ );
}

double
nest::error_neuron::Parameters_::set(
  const DictionaryDatum& d )
{
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
  updateValue< double >( d, names::start, t_start_ls_ );

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

void
nest::error_neuron::State_::get(
  DictionaryDatum& d, const Parameters_& p ) const
{
  def< double >( d, names::target_rate, target_rate_); // target_rate
  def< double >( d, names::learning_signal, learning_signal_ );
  def< double >( d, names::V_m, y3_ + p.E_L_ ); // Membrane potential
}

void
nest::error_neuron::State_::set(
  const DictionaryDatum& d, const Parameters_& p, double delta_EL)
{
  updateValue< double >( d, names::target_rate, target_rate_ ); // target_rate
  updateValue< double >( d, names::learning_signal, learning_signal_ );

  if ( updateValue< double >( d, names::V_m, y3_ ) )
  {
    y3_ -= p.E_L_;
  }
  else
  {
    y3_ -= delta_EL;
  }
}


nest::error_neuron::Buffers_::Buffers_(
  error_neuron& n )
  : logger_( n )
{
}

nest::error_neuron::Buffers_::Buffers_(
  const Buffers_&,
  error_neuron& n )
  : logger_( n )
{
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::error_neuron::error_neuron()
  : Eprop_Archiving_Node()
  , P_()
  , S_()
  , B_( *this )
{
  recordablesMap_.create();
}

nest::error_neuron::error_neuron(
  const error_neuron& n )
  : Eprop_Archiving_Node( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
nest::error_neuron::init_state_( const Node& proto )
{
  const error_neuron& pr = downcast< error_neuron >( proto );
  S_ = pr.S_;
}

void
nest::error_neuron::init_buffers_()
{
  B_.delayed_rates_.clear(); // includes resize

  B_.logger_.reset(); // includes resize
  B_.spikes_.clear();   // includes resize
  B_.currents_.clear(); // includes resize
  Archiving_Node::clear_history();
}

void
nest::error_neuron::calibrate()
{
  B_.logger_
    .init(); // ensures initialization in case mm connected after Simulate

  const double h = Time::get_resolution().get_ms();
  V_.P33_ = std::exp( -h / P_.tau_m_ );
  V_.P30_ = 1 / P_.c_m_ * ( 1 - V_.P33_ ) * P_.tau_m_;
  V_.step_start_ls_ = Time( Time::ms( P_.t_start_ls_ ) ).get_steps();
  V_.prnt = true;
  //std::cout << "start = " << V_.step_start_ls_ << std::endl;
}

/* ----------------------------------------------------------------
 * Update and event handling functions
 */

void
nest::error_neuron::update_( Time const& origin,
  const long from,
  const long to)
{
  assert(
    to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  const double h = Time::get_resolution().get_ms();
  const size_t buffer_size = kernel().connection_manager.get_min_delay();

  // allocate memory to store rates to be sent by rate events
  std::vector< double > readout_and_target_signals( 3*buffer_size, 0.0 );

  for ( long lag = from; lag < to; ++lag )
  {
    // DEBUG: added reset after each T to be compatible with tf code
    int t_mod_T = ( origin.get_steps() + lag - 2 ) % get_update_interval_steps();
    if ( t_mod_T == 0 )
    {
      S_.y3_ = 0.0;
      B_.spikes_.clear();   // includes resize
      V_.prnt = true;
    }
    // DEBUG: introduced factor ( 1 - exp( -dt / tau_m ) ) for campatibility wit tf code
      S_.y3_ = V_.P30_ * ( S_.y0_ + P_.I_e_ ) + V_.P33_ * S_.y3_ + ( 1 - V_.P33_ ) * B_.spikes_.get_value( lag );
      S_.y3_ = ( S_.y3_ < P_.V_min_ ? P_.V_min_ : S_.y3_ );

      // DEBUG: changed sign (see tf code) (maybe this is not true any more)
      // S_.learning_signal_ = ( S_.target_rate_ - (S_.y3_ + P_.E_L_) );
      // 0: regression 1: classification
      readout_and_target_signals [ 3*lag ] = 0.;  // TODO replace by settable parameter
      if ( t_mod_T > V_.step_start_ls_ )
      {
        // TODO: replace -1 by ls
        //std::cout << "error neuron send: " << S_.learning_signal_ << std::endl;
        // new_learning_signals [ lag ] = S_.learning_signal_;
        readout_and_target_signals [ 3*lag + 1 ] =  S_.y3_ + P_.E_L_;
        readout_and_target_signals [ 3*lag + 2 ] = S_.target_rate_;
        // std::cout << "rs error " << readout_and_target_signals [ lag + 1 ] << std::endl;
        // std::cout << "ts error " << readout_and_target_signals [ lag + 2 ] << std::endl;
        // std::cout << "...." << std::endl;
      }
      else
      {
        // new_learning_signals[ lag ] = 0.0;
        readout_and_target_signals [ 3*lag + 1 ] = 0.;
        readout_and_target_signals [ 3*lag + 2 ] = 0.;
      }

      S_.y0_ = B_.currents_.get_value( lag ); // set new input current
      S_.target_rate_ =  1. * B_.delayed_rates_.get_value( lag );

      B_.logger_.record_data( origin.get_steps() + lag );
      write_readout_history( Time::step( origin.get_steps() + lag + 1), S_.y3_ + P_.E_L_,
          S_.target_rate_);
  }

  // Send delay-rate-neuron-event. This only happens in the final iteration
  // to avoid accumulation in the buffers of the receiving neurons.
  LearningSignalConnectionEvent drve;
  drve.set_coeffarray( readout_and_target_signals );
  kernel().event_delivery_manager.send_secondary( *this, drve );

  return;
}

bool
nest::error_neuron::is_eprop_readout()
    {
        return true;
    }

void
nest::error_neuron::handle(
  DelayedRateConnectionEvent& e )
{
  const double weight = e.get_weight();
  const long delay = e.get_delay_steps();

  size_t i = 0;
  std::vector< unsigned int >::iterator it = e.begin();
  // The call to get_coeffvalue( it ) in this loop also advances the iterator it
  while ( it != e.end() )
  {
    B_.delayed_rates_.add_value(delay + i, weight * 1. * e.get_coeffvalue( it ) ) ;
    ++i;
  }
}

void
nest::error_neuron::handle( SpikeEvent& e )
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

void
nest::error_neuron::handle( CurrentEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  const double c = e.get_current();
  const double w = e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(
    e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
    w * c );
}

void
nest::error_neuron::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

} // namespace
