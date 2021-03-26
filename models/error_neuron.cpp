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
  // DEBUG II: the commented line below was used in the pattern generation task
  //insert_( names::learning_signal, &error_neuron::get_learning_signal_ );
  insert_( names::learning_signal, &error_neuron::get_last_ls_ );
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
  , regression_( true )
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
  def< bool >( d, names::regression, regression_ );
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
  updateValue< bool >( d, names::regression, regression_ );

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
  : EpropArchivingNode()
  , P_()
  , S_()
  , B_( *this )
{
  recordablesMap_.create();
}

nest::error_neuron::error_neuron(
  const error_neuron& n )
  : EpropArchivingNode( n )
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
  ArchivingNode::clear_history();
}

void
nest::error_neuron::calibrate()
{
  B_.logger_
    .init(); // ensures initialization in case mm connected after Simulate

  const double h = Time::get_resolution().get_ms();
  V_.P33_ = std::exp( -h / P_.tau_m_ );
  V_.P30_ = 1 / P_.c_m_ * ( 1 - V_.P33_ ) * P_.tau_m_;
  V_.step_start_ls_ = Time( Time::ms( std::max( P_.t_start_ls_, 0.0 ) + h ) ).get_steps();
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

  const size_t buffer_size = kernel().connection_manager.get_min_delay();

  // allocate memory to store learning signal to be sent by learning signal events
  size_t n_entries = 2;
  std::vector< double > readout_and_target_signals( n_entries*buffer_size, 0.0 );

  // allocate memory to store readout signal to be sent by rate events
  std::vector< double > readout_signal_buffer( buffer_size, 0.0 );

  for ( long lag = from; lag < to; ++lag )
  {
    // DEBUG: added reset after each T to be compatible with tf code
    int t_mod_T = ( origin.get_steps() + lag - 2 ) % get_update_interval_steps();
    if ( t_mod_T == 0 )
    {
      S_.y3_ = 0.0;
      B_.spikes_.clear();   // includes resize
    }
    // DEBUG: introduced factor ( 1 - exp( -dt / tau_m ) ) for campatibility wit tf code
    S_.y3_ = V_.P30_ * ( S_.y0_ + P_.I_e_ ) + V_.P33_ * S_.y3_ + ( 1 - V_.P33_ ) * B_.spikes_.get_value( lag );
    S_.y3_ = ( S_.y3_ < P_.V_min_ ? P_.V_min_ : S_.y3_ );

    // compute the readout signal
    double readout_signal = S_.y3_ + P_.E_L_;
    // write exp( readout signal ) into the buffer which is used to send it to the other error neurons
    // in case of a regression task we don't need this and therefore set it to zero
    readout_signal_buffer[ lag ] = P_.regression_ ? 0.0 : std::exp( readout_signal );

    // DEBUG: changed sign (see tf code)
    readout_and_target_signals[ n_entries*lag ] = origin.get_steps() + lag + 1;

    // compute normalized learning signal from values stored in state_buffer_
    // which now contains the correct normalization because in the meantime
    // the other readout neurons have sent their membrane potential
    // the entries of the state_buffer_ are
    // 0: readout_signal, 1: target_signal, 2: normalization
    double normalized_learning_signal;
    if ( t_mod_T >= V_.step_start_ls_ )
    {
      // if recall is active, compute normalized learning signal
      if ( P_.regression_ )
      {
        // if this is a regression task, use the bare membrane potential
        normalized_learning_signal = V_.state_buffer_ [ 0 ] /
        V_.state_buffer_ [ 2 ] - V_.state_buffer_ [ 1 ];
      }
      else
      {
        // if this is a classification task, use exp( membrane potential )
        normalized_learning_signal = std::exp( V_.state_buffer_ [ 0 ] ) /
        V_.state_buffer_ [ 2 ] - V_.state_buffer_ [ 1 ];
      }
    }
    else
    {
      // if recall is inactive, set normalized learning signal to zero
      normalized_learning_signal = 0.0;
    }

    // fill the state buffer with new values
    if ( t_mod_T >= V_.step_start_ls_ - 1 )
    {
      // if the recall is active, fill state_buffer_ with the current state
      V_.state_buffer_ [ 0 ] = readout_signal;
      V_.state_buffer_ [ 1 ] = S_.target_rate_;
      V_.state_buffer_ [ 2 ] = 1.0;
    }
    else
    {
      // if the recall is inactive, fill state_buffer_ with zeros
      V_.state_buffer_ [ 0 ] = 0.0;
      V_.state_buffer_ [ 1 ] = 0.0;
      V_.state_buffer_ [ 2 ] = 1.0;
    }

    // write normalized learning signal into history. Use the previous time step:
    // origin.get_steps() + lag (without + 1) because of the buffering in
    // readout_signal_buffer
    // previously, the function write_readout_history was used for this purpose
    Time const& t_norm_ls = Time::step( origin.get_steps() + lag );
    const double t_norm_ls_ms = t_norm_ls.get_ms();
    eprop_history_.push_back( histentry_eprop( t_norm_ls_ms,
          0.0, normalized_learning_signal, 0 ) );

    // store the normalized learning signal in the buffer which is send to
    // the recurrent neurons via the learning signal connection
    readout_and_target_signals[ n_entries*lag + 1 ] = normalized_learning_signal;
    S_.y0_ = B_.currents_.get_value( lag ); // set new input current
    S_.target_rate_ =  B_.delayed_rates_.get_value( lag );

    B_.logger_.record_data( origin.get_steps() + lag );
  }

  // time as it is in the last iteration of the for loop modulo update interval
  int t_mod_T_final = ( origin.get_steps() + to - 3 ) % get_update_interval_steps();
  // send learning signal and readout signal only if recall is active
  if ( t_mod_T_final >= V_.step_start_ls_ )
  {
    // send learning signal
    // TODO: it would be much more efficient to send this in larger batches
    LearningSignalConnectionEvent drve;
    drve.set_coeffarray( readout_and_target_signals );
    kernel().event_delivery_manager.send_secondary( *this, drve );
  }
  // time one time step larger than t_mod_T_final because the readout has to
  // be sent one time step in advance so that the normalization can be computed
  // and the learning signal is ready as soon as the recall starts.
  if ( !P_.regression_ )
  {
    int t_mod_T_final_p1 = ( origin.get_steps() + to - 2 ) % get_update_interval_steps();
    if ( t_mod_T_final_p1 >= V_.step_start_ls_ )
    {
      // send readout signal only if this is a classification task
      // rate connection to connect to other readout neurons
      DelayedRateConnectionEvent readout_event;
      readout_event.set_coeffarray( readout_signal_buffer );
      kernel().event_delivery_manager.send_secondary( *this, readout_event );
    }
  }

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
  assert( 0 <= e.get_rport() && e.get_rport() < 3 );
  const double weight = e.get_weight();
  const long delay = e.get_delay_steps();
  size_t i = 0;
  std::vector< unsigned int >::iterator it = e.begin();
  if ( e.get_rport() == READOUT_SIG - MIN_RATE_RECEPTOR )
  {
    // handle port for readout signal
    // The call to get_coeffvalue( it ) in this loop also advances the iterator it
    while ( it != e.end() )
    {
      double readout_signal = e.get_coeffvalue( it );
      V_.state_buffer_ [ 2 ] += readout_signal;
      ++i;
    }
  }
  else  if ( e.get_rport() == TARGET_SIG - MIN_RATE_RECEPTOR )
  {
    // handle port for target signal
    while ( it != e.end() )
    {
      // The call to get_coeffvalue( it ) in this loop also advances the iterator it
      B_.delayed_rates_.add_value(delay + i, weight * 1. * e.get_coeffvalue( it ) ) ;
      ++i;
    }
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
