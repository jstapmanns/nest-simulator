/*
 *  clopath_stdp_connection.h
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

#ifndef CLOPATH_STDP_CONNECTION_H
#define CLOPATH_STDP_CONNECTION_H

/* BeginDocumentation
  Name: clopath_stdp_connection - TODO add documentation

  Description:


  Parameters:
   tau_plus   double - Time constant of STDP window, potentiation in ms
                       (tau_minus defined in post-synaptic neuron)
   lambda     double - Step size
   alpha      double - Asymmetry parameter (scales depressing increments as
                       alpha*lambda)
   mu_plus    double - Weight dependence exponent, potentiation
   mu_minus   double - Weight dependence exponent, depression
   Wmax       double - Maximum allowed weight

  Transmits: SpikeEvent

  References:
   [1] Clopath et al. (2010) Connectivity reflects coding:
       a model of voltage-based STDP with homeostasis.
       Nature Neuroscience 13:3,344--352

  Author: Jonas Stapmanns, David Dahmen, Jan Hahne

  SeeAlso: synapsedict, tsodyks_synapse, stdp_synapse
*/

// C++ includes:
#include <cmath>

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connection.h"
#include "connector_model.h"
#include "event.h"
#include "ring_buffer.h"

// Includes from sli:
#include "dictdatum.h"
#include "dictutils.h"


namespace nest
{

// connections are templates of target identifier type (used for pointer /
// target index addressing) derived from generic connection template
template < typename targetidentifierT >
class Clopath_STDPConnection : public Connection< targetidentifierT >
{

public:
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  Clopath_STDPConnection();


  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  Clopath_STDPConnection( const Clopath_STDPConnection& );

  // Explicitly declare all methods inherited from the dependent base
  // ConnectionBase. This avoids explicit name prefixes in all places these
  // functions are used. Since ConnectionBase depends on the template parameter,
  // they are not automatically found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_delay;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set properties of this connection from the values given in dictionary.
   */
  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param cp common properties of all synapses (empty).
   */
  void send( Event& e, thread t, const CommonSynapseProperties& cp );


  class ConnTestDummyNode : public ConnTestDummyNodeBase
  {
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using ConnTestDummyNodeBase::handles_test_event;
    port
    handles_test_event( SpikeEvent&, rport )
    {
      return invalid_port_;
    }
  };

  void
  check_connection( Node& s,
    Node& t,
    rport receptor_type,
    const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;

    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    t.register_stdp_connection( t_lastspike_ - get_delay() );
  }

  void
  set_weight( double w )
  {
    weight_ = w;
  }

private:
  double
  facilitate_( double w, double dw, double x_bar )
  {
    w += dw * x_bar;
    return w < Wmax_ ? w : Wmax_;
  }

  double
  depress_( double w, double dw )
  {
    w -= dw;
    return w > Wmin_ ? w : Wmin_;
  }

  // data members of each connection
  double weight_;
  double x_bar_;
  double tau_x_; // TODO: save tau_x in synapse?
  double Wmax_;
  double Wmin_;
  std::vector< double > delayed_u_bar_plus_;
  std::vector< double > delayed_u_bar_minus_;
  size_t read_idx_;
  size_t delay_length_;

  double t_lastspike_;
};


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param t The thread on which this connection is stored.
 * \param cp Common properties object, containing the stdp parameters.
 */
template < typename targetidentifierT >
inline void
Clopath_STDPConnection< targetidentifierT >::send( Event& e,
  thread t,
  const CommonSynapseProperties& )
{
  // const double old_w = weight_;
  double t_spike = e.get_stamp().get_ms();

  // use accessor functions (inherited from Connection< >) to obtain delay and
  // target
  Node* target = get_target( t );
  double dendritic_delay = get_delay();

  // get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque< histentry_cl >::iterator start;
  std::deque< histentry_cl >::iterator finish;

  // For a new synapse, t_lastspike_ contains the point in time of the last
  // spike. So we initially read the
  // history(t_last_spike - dendritic_delay, ..., T_spike-dendritic_delay]
  // which increases the access counter for these entries.
  // At registration, all entries' access counters of
  // history[0, ..., t_last_spike - dendritic_delay] have been
  // incremented by Archiving_Node::register_stdp_connection(). See bug #218 for
  // details.
  target->get_LTP_history( t_lastspike_ - dendritic_delay,
    t_spike - dendritic_delay,
    &start,
    &finish );
  // facilitation due to post-synaptic spikes since last pre-synaptic spike
  double minus_dt;
  while ( start != finish )
  {
    minus_dt = t_lastspike_ - ( start->t_ + dendritic_delay );
    if ( minus_dt == 0 )
    {
      continue;
    }
    weight_ =
      facilitate_( weight_, start->dw_, x_bar_ * exp( minus_dt / tau_x_ ) );
    ++start;
  }

  // depression due to new pre-synaptic spike
  weight_ =
    depress_( weight_, target->get_LTD_value( t_spike - dendritic_delay ) );

  e.set_receiver( *target );
  e.set_weight( weight_ );
  // use accessor functions (inherited from Connection< >) to obtain delay in
  // steps and rport
  e.set_delay_steps( get_delay_steps() );
  e.set_rport( get_rport() );
  e();

  x_bar_ =
    x_bar_ * std::exp( ( t_lastspike_ - t_spike ) / tau_x_ ) + 1.0 / tau_x_;

  t_lastspike_ = t_spike;
}


template < typename targetidentifierT >
Clopath_STDPConnection< targetidentifierT >::Clopath_STDPConnection()
  : ConnectionBase()
  , weight_( 1.0 )
  , x_bar_( 0.0 )
  , tau_x_( 15.0 )
  , Wmin_( 0.0 )
  , Wmax_( 100.0 )
  , read_idx_( 0 )
  , delay_length_( 100 )
  , t_lastspike_( 0.0 )
{
  delayed_u_bar_plus_.resize( delay_length_ );
  delayed_u_bar_minus_.resize( delay_length_ );
}

template < typename targetidentifierT >
Clopath_STDPConnection< targetidentifierT >::Clopath_STDPConnection(
  const Clopath_STDPConnection< targetidentifierT >& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , x_bar_( rhs.x_bar_ )
  , tau_x_( rhs.tau_x_ )
  , Wmin_( rhs.Wmin_ )
  , Wmax_( rhs.Wmax_ )
  , t_lastspike_( rhs.t_lastspike_ )
{
}

template < typename targetidentifierT >
void
Clopath_STDPConnection< targetidentifierT >::get_status(
  DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< double >( d, names::x_bar, x_bar_ );
  def< double >( d, names::tau_x, tau_x_ );
  def< double >( d, names::Wmin, Wmin_ );
  def< double >( d, names::Wmax, Wmax_ );
  def< long >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
Clopath_STDPConnection< targetidentifierT >::set_status(
  const DictionaryDatum& d,
  ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
  updateValue< double >( d, names::x_bar, x_bar_ );
  updateValue< double >( d, names::tau_x, tau_x_ );
  updateValue< double >( d, names::Wmin, Wmin_ );
  updateValue< double >( d, names::Wmax, Wmax_ );

  // check if weight_ and Wmax_ has the same sign
  if ( not( ( ( weight_ >= 0 ) - ( weight_ < 0 ) )
         == ( ( Wmax_ >= 0 ) - ( Wmax_ < 0 ) ) ) )
  {
    throw BadProperty( "Weight and Wmax must have same sign." );
  }

  // check if weight_ and Wmin_ has the same sign
  if ( not( ( ( weight_ >= 0 ) - ( weight_ < 0 ) )
         == ( ( Wmin_ >= 0 ) - ( Wmin_ < 0 ) ) ) )
  {
    throw BadProperty( "Weight and Wmin must have same sign." );
  }
}

} // of namespace nest

#endif // of #ifndef CLOPATH_STDP_CONNECTION_H