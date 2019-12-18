/*
 *  urbanczik_connection.h
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

#ifndef URBANCZIK_CONNECTION_BC_H
#define URBANCZIK_CONNECTION_BC_H

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

/** @BeginDocumentation
@ingroup Synapses
@ingroup stdp

Name: urbanczik_synapse - Synapse type for a plastic synapse after Urbanczik and Senn.

Description:

urbanczik_synapse is a connector to create Urbanczik synapses as defined in [1] that can connect suitable
multicompartment models. In contrast to usual STDP, the change of the synaptic weight does not only depend on the pre-
and postsynaptic spike timing but also on the postsynaptic dendritic potential.

Urbanczik synapses require archiving of continuous quantities. Therefore they can only be connected to neuron models
that are capable of doing this archiving. So far, the only compatible model is pp_cond_exp_mc_urbanczik.

Parameters:

\verbatim embed:rst
=======   ======  ==========================================================
eta       real    Learning rate
tau_Delta real    Time constant of low pass filtering of the weight change
Wmax      real    Maximum allowed weight
Wmin      real    Minimum allowed weight
=======   ======  ==========================================================
\endverbatim

All other parameters are stored in in the neuron models that are compatible with the Urbanczik synapse.

Remarks:

So far the implementation of the urbanczik_synapse only supports two-compartment neurons.

Transmits: SpikeEvent

References:

\verbatim embed:rst
.. [1] R. Urbanczik, W. Senn (2014). Learning by the Dendritic Prediction of Somatic Spiking. Neuron, 81, 521 - 528.
\endverbatim
Authors: Jonas Stapmanns, David Dahmen, Jan Hahne

SeeAlso: stdp_synapse, clopath_synapse, pp_cond_exp_mc_urbanczik
*/
// connections are templates of target identifier type (used for pointer /
// target index addressing) derived from generic connection template
template < typename targetidentifierT >
class UrbanczikConnectionBC : public Connection< targetidentifierT >
{

public:
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  UrbanczikConnectionBC();


  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  UrbanczikConnectionBC( const UrbanczikConnectionBC& );

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
  check_connection( Node& s, Node& t, rport receptor_type, const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;

    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    t.register_stdp_connection( t_lastspike_ - get_delay(), get_delay() );
  }

  void
  set_weight( double w )
  {
    weight_ = w;
  }

private:
  // data members of each connection
  double weight_;
  double init_weight_;
  double tau_Delta_;
  double eta_;
  double Wmin_;
  double Wmax_;
  double PI_integral_L_;
  double PI_integral_s_;
  double PI_exp_integral_L_;
  double PI_exp_integral_s_;
  double tau_L_trace_;
  double tau_s_trace_;

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
UrbanczikConnectionBC< targetidentifierT >::send( Event& e, thread t, const CommonSynapseProperties& )
{
  double t_spike = e.get_stamp().get_ms();
  // use accessor functions (inherited from Connection< >) to obtain delay and target
  Node* target = get_target( t );
  double dendritic_delay = get_delay();

  // get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque< histentry_extended >::iterator start;
  std::deque< histentry_extended >::iterator finish;

  // for now we only support two-compartment neurons
  // in this case the dendritic compartment has index 1
  const int comp = 1;

  // compress history
  target->compress_urbanczik_history( t_spike - dendritic_delay, tau_Delta_, comp );
  // get weight change
  double const g_L = target->get_g_L( comp );
  double const tau_L = target->get_tau_L( comp );
  double const C_m = target->get_C_m( comp );
  double const tau_s = weight_ > 0.0 ? target->get_tau_syn_ex( comp ) : target->get_tau_syn_in( comp );
  double I1_L;
  double I1_s;
  double I2_L;
  double I2_s;
  target->get_urbanczik_value( t_lastspike_ - dendritic_delay, I1_L, I1_s, I2_L, I2_s, comp );
  PI_integral_L_ += tau_L_trace_ * I1_L;
  PI_integral_s_ += tau_s_trace_ * I1_s;
  PI_exp_integral_L_ = ( exp( ( t_lastspike_ - t_spike ) / tau_Delta_ ) * PI_exp_integral_L_
      + tau_L_trace_ * I2_L );
  PI_exp_integral_s_ = ( exp( ( t_lastspike_ - t_spike ) / tau_Delta_ ) * PI_exp_integral_s_
      + tau_s_trace_ * I2_s );
  //std::cout << "summed integrals: " << PI_integral_L_ - PI_integral_s_ - PI_exp_integral_L_ + PI_exp_integral_s_;
  //std::cout << "  PI: " << PI_integral_L_ - PI_integral_s_ << "  PI_exp: " << PI_exp_integral_L_ - PI_exp_integral_s_ << std::endl;
  weight_ = ( init_weight_ + ( PI_integral_L_ - PI_integral_s_ - PI_exp_integral_L_ + PI_exp_integral_s_ )
      * 15.0 * C_m * tau_s * eta_ / ( g_L * ( tau_L - tau_s ) ) );

  if ( weight_ > Wmax_ )
  {
    weight_ = Wmax_;
  }
  else if ( weight_ < Wmin_ )
  {
    weight_ = Wmin_;
  }

  e.set_receiver( *target );
  e.set_weight( weight_ );
  // use accessor functions (inherited from Connection< >) to obtain delay in steps and rport
  e.set_delay_steps( get_delay_steps() );
  e.set_rport( get_rport() );
  e();

  // compute the trace of the presynaptic spike train
  tau_L_trace_ = tau_L_trace_ * std::exp( ( t_lastspike_ - t_spike ) / tau_L ) + 1.0;
  tau_s_trace_ = tau_s_trace_ * std::exp( ( t_lastspike_ - t_spike ) / tau_s ) + 1.0;

  t_lastspike_ = t_spike;
}


template < typename targetidentifierT >
UrbanczikConnectionBC< targetidentifierT >::UrbanczikConnectionBC()
  : ConnectionBase()
  , weight_( 1.0 )
  , init_weight_( 1.0 )
  , tau_Delta_( 100.0 )
  , eta_( 0.07 )
  , Wmin_( 0.0 )
  , Wmax_( 100.0 )
  , PI_integral_L_( 0.0 )
  , PI_integral_s_( 0.0 )
  , PI_exp_integral_L_( 0.0 )
  , PI_exp_integral_s_( 0.0 )
  , tau_L_trace_( 0.0 )
  , tau_s_trace_( 0.0 )
  , t_lastspike_( -1.0 )
{
}

template < typename targetidentifierT >
UrbanczikConnectionBC< targetidentifierT >::UrbanczikConnectionBC( const UrbanczikConnectionBC< targetidentifierT >& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , init_weight_( rhs.init_weight_ )
  , tau_Delta_( rhs.tau_Delta_ )
  , eta_( rhs.eta_ )
  , Wmin_( rhs.Wmin_ )
  , Wmax_( rhs.Wmax_ )
  , PI_integral_L_( rhs.PI_integral_L_ )
  , PI_integral_s_( rhs.PI_integral_s_ )
  , PI_exp_integral_L_( rhs.PI_exp_integral_L_ )
  , PI_exp_integral_s_( rhs.PI_exp_integral_s_ )
  , tau_L_trace_( rhs.tau_L_trace_ )
  , tau_s_trace_( rhs.tau_s_trace_ )
  , t_lastspike_( rhs.t_lastspike_ )
{
}

template < typename targetidentifierT >
void
UrbanczikConnectionBC< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< double >( d, names::tau_Delta, tau_Delta_ );
  def< double >( d, names::eta, eta_ );
  def< double >( d, names::Wmin, Wmin_ );
  def< double >( d, names::Wmax, Wmax_ );
  def< long >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
UrbanczikConnectionBC< targetidentifierT >::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
  updateValue< double >( d, names::tau_Delta, tau_Delta_ );
  updateValue< double >( d, names::eta, eta_ );
  updateValue< double >( d, names::Wmin, Wmin_ );
  updateValue< double >( d, names::Wmax, Wmax_ );

  init_weight_ = weight_;
  // check if weight_ and Wmin_ has the same sign
  if ( not( ( ( weight_ >= 0 ) - ( weight_ < 0 ) ) == ( ( Wmin_ >= 0 ) - ( Wmin_ < 0 ) ) ) )
  {
    throw BadProperty( "Weight and Wmin must have same sign." );
  }

  // check if weight_ and Wmax_ has the same sign
  if ( not( ( ( weight_ >= 0 ) - ( weight_ < 0 ) ) == ( ( Wmax_ > 0 ) - ( Wmax_ <= 0 ) ) ) )
  {
    throw BadProperty( "Weight and Wmax must have same sign." );
  }
}

} // of namespace nest

#endif // of #ifndef URBANCZIK_CONNECTION_BC_H
