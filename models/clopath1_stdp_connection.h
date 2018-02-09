/*
 *  stdp_connection.h
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

#ifndef CLOPATH1_STDP_CONNECTION_H
#define CLOPATH1_STDP_CONNECTION_H

/* BeginDocumentation
  Name: stdp_synapse - Synapse type for spike-timing dependent
   plasticity.

  Description:
   stdp_synapse is a connector to create synapses with spike time
   dependent plasticity (as defined in [1]). Here the weight dependence
   exponent can be set separately for potentiation and depression.

  Examples:
   multiplicative STDP [2]  mu_plus = mu_minus = 1.0
   additive STDP       [3]  mu_plus = mu_minus = 0.0
   Guetig STDP         [1]  mu_plus = mu_minus = [0.0,1.0]
   van Rossum STDP     [4]  mu_plus = 0.0 mu_minus = 1.0

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
   [1] Guetig et al. (2003) Learning Input Correlations through Nonlinear
       Temporally Asymmetric Hebbian Plasticity. Journal of Neuroscience

   [2] Rubin, J., Lee, D. and Sompolinsky, H. (2001). Equilibrium
       properties of temporally asymmetric Hebbian plasticity, PRL
       86,364-367

   [3] Song, S., Miller, K. D. and Abbott, L. F. (2000). Competitive
       Hebbian learning through spike-timing-dependent synaptic
       plasticity,Nature Neuroscience 3:9,919--926

   [4] van Rossum, M. C. W., Bi, G-Q and Turrigiano, G. G. (2000).
       Stable Hebbian learning from spike timing-dependent
       plasticity, Journal of Neuroscience, 20:23,8812--8821

  FirstVersion: March 2006
  Author: Moritz Helias, Abigail Morrison
  Adapted by: Philipp Weidel
  SeeAlso: synapsedict, tsodyks_synapse, static_synapse
*/

// C++ includes:
#include <cmath>

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connection.h"
#include "connector_model.h"
#include "event.h"

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
   * \param t_lastspike Point in time of last spike sent.
   * \param cp common properties of all synapses (empty).
   */
  void send( Event& e,
    thread t,
    double t_lastspike,
    const CommonSynapseProperties& cp );


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
    double t_lastspike,
    const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;

    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    t.register_stdp_connection( t_lastspike - get_delay() );
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
	//std::cout << "facilitate!" << std::endl;
    w += dw * x_bar;
    return w < Wmax_ ? w : Wmax_;
  }

  double
  depress_( double w, double dw )
  {
    w -= dw;
    return w > 0.0 ? w : 0.0;
  }

  // data members of each connection
  double weight_;
  double x_bar_;
  double tau_x_;  // TO DO: save tau_x in synapse?
  double Wmax_;
};


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param t The thread on which this connection is stored.
 * \param t_lastspike Time point of last spike emitted
 * \param cp Common properties object, containing the stdp parameters.
 */
template < typename targetidentifierT >
inline void
Clopath_STDPConnection< targetidentifierT >::send( Event& e,
  thread t,
  double t_lastspike,
  const CommonSynapseProperties& )
{
  const double old_w = weight_;
  // synapse STDP depressing/facilitation dynamics
  //   if(t_lastspike >0) {std::cout << "last spike " << t_lastspike <<
  //   std::endl ;}
  double t_spike = e.get_stamp().get_ms();
  // t_lastspike_ = 0 initially

  // use accessor functions (inherited from Connection< >) to obtain delay and
  // target
  Node* target = get_target( t );
  double dendritic_delay = get_delay();

  // get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque< histentry_cl >::iterator start;
  std::deque< histentry_cl >::iterator finish;

  // For a new synapse, t_lastspike contains the point in time of the last
  // spike. So we initially read the
  // history(t_last_spike - dendritic_delay, ..., T_spike-dendritic_delay]
  // which increases the access counter for these entries.
  // At registration, all entries' access counters of
  // history[0, ..., t_last_spike - dendritic_delay] have been
  // incremented by Archiving_Node::register_stdp_connection(). See bug #218 for
  // details.
  //std::cout << "t_lastspike = " << t_lastspike << "  t_spike = " << t_spike << std::endl;
  target->get_LTP_history(
    t_lastspike - dendritic_delay, t_spike - dendritic_delay, &start, &finish );
  // facilitation due to post-synaptic spikes since last pre-synaptic spike
  double minus_dt;
  while ( start != finish )
  {
    minus_dt = t_lastspike - ( start->t_ + dendritic_delay );
    if ( minus_dt == 0 )
    {
      continue;
    }
    weight_ = facilitate_( weight_, start->dw_,  
       x_bar_ * exp( minus_dt / tau_x_ ) );
    ++start;
  }

  const double fa_weight = weight_;
  // depression due to new pre-synaptic spike
  weight_ =
    depress_( weight_, target->get_LTD_value( t_spike - dendritic_delay) );

  e.set_receiver( *target );
  std::cout << "facilitation:  " << fa_weight - old_w << "   depression:  " << weight_ - fa_weight 
    << "   delta:  " << weight_ - old_w << "   synapse weight:  " << weight_ << std::endl;
  e.set_weight( weight_ );
  // use accessor functions (inherited from Connection< >) to obtain delay in
  // steps and rport
  e.set_delay( get_delay_steps() );
  e.set_rport( get_rport() );
  e();

  x_bar_ = x_bar_ * std::exp( ( t_lastspike - t_spike ) / tau_x_ ) + 1.0;
}


template < typename targetidentifierT >
Clopath_STDPConnection< targetidentifierT >::Clopath_STDPConnection()
  : ConnectionBase()
  , weight_( 1.0 )
  , x_bar_( 0.0 )
  , tau_x_( 15.0 )
  , Wmax_( 100.0 )
{
}

template < typename targetidentifierT >
Clopath_STDPConnection< targetidentifierT >::Clopath_STDPConnection(
  const Clopath_STDPConnection< targetidentifierT >& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , x_bar_( rhs.x_bar_ )
  , tau_x_( rhs.tau_x_ )
  , Wmax_( rhs.Wmax_ )
{
}

template < typename targetidentifierT >
void
Clopath_STDPConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< double >( d, names::x_bar, x_bar_ );
  def< double >( d, names::tau_x, tau_x_ );
  def< double >( d, names::Wmax, Wmax_ );
  def< long >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
Clopath_STDPConnection< targetidentifierT >::set_status( const DictionaryDatum& d,
  ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
  updateValue< double >( d, names::x_bar, x_bar_ );
  updateValue< double >( d, names::tau_x, tau_x_ );
  updateValue< double >( d, names::Wmax, Wmax_ );

  // check if weight_ and Wmax_ has the same sign
  if ( not( ( ( weight_ >= 0 ) - ( weight_ < 0 ) )
         == ( ( Wmax_ >= 0 ) - ( Wmax_ < 0 ) ) ) )
  {
    throw BadProperty( "Weight and Wmax must have same sign." );
  }
}

} // of namespace nest

#endif // of #ifndef CLOPATH1_STDP_CONNECTION_H
