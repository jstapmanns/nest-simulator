/*
 *  clopath_connection_bc.h
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

#ifndef CLOPATH_CONNECTION_BC_H
#define CLOPATH_CONNECTION_BC_H

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
@ingroup clopath_s

Name: clopath_synapse - Synapse type for voltage-based STDP after Clopath.

Description:

clopath_synapse is a connector to create Clopath synapses as defined
in [1]. In contrast to usual STDP, the change of the synaptic weight does
not only depend on the pre- and postsynaptic spike timing but also on the
postsynaptic membrane potential.

Clopath synapses require archiving of continuous quantities. Therefore Clopath
synapses can only be connected to neuron models that are capable of doing this
archiving. So far, compatible models are aeif_psc_delta_clopath and
hh_psc_alpha_clopath.

Parameters:

\verbatim embed:rst
=======  ======  ==========================================================
tau_x    ms      Time constant of the trace of the presynaptic spike train
Wmax     real    Maximum allowed weight
Wmin     real    Minimum allowed weight
=======  ======  ==========================================================
\endverbatim

Other parameters like the amplitudes for long-term potentiation (LTP) and
depression (LTD) are stored in in the neuron models that are compatible with the
Clopath synapse.

Transmits: SpikeEvent

References:

\verbatim embed:rst
.. [1] Clopath et al. (2010). Connectivity reflects coding:
       a model of voltage-based STDP with homeostasis.
       Nature Neuroscience 13:3, 344--352. DOI: https://doi.org/10.1038/nn.2479
.. [2] Clopath and Gerstner (2010). Voltage and spike timing interact
       in STDP – a unified model. Frontiers in Synaptic Neuroscience 2:25.
       DOI: https://doi.org/10.3389/fnsyn.2010.00025
.. [3] Voltage-based STDP synapse (Clopath et al. 2010) on ModelDB
       https://senselab.med.yale.edu/ModelDB/showmodel.cshtml?model=144566
\endverbatim
Authors: Jonas Stapmanns, David Dahmen, Jan Hahne

SeeAlso: stdp_synapse, aeif_psc_delta_clopath, hh_psc_alpha_clopath
*/
// connections are templates of target identifier type (used for pointer /
// target index addressing) derived from generic connection template
template < typename targetidentifierT >
class ClopathConnectionBC : public Connection< targetidentifierT >
{

public:
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  ClopathConnectionBC();


  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  ClopathConnectionBC( const ClopathConnectionBC& );

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
  double
  depress_( double w, double dw )
  {
    w -= eta_ * dw;
    return w > Wmin_ ? w : Wmin_;
  }

  double
  facilitate_( double w, double dw, double x_bar )
  {
    w += eta_ * dw * x_bar;
    return w < Wmax_ ? w : Wmax_;
  }

  // data members of each connection
  double weight_;
  double x_bar_;
  double tau_x_;
  double Wmin_;
  double Wmax_;
  double eta_;

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
ClopathConnectionBC< targetidentifierT >::send( Event& e, thread t, const CommonSynapseProperties& )
{
  double t_spike = e.get_stamp().get_ms();
  // use accessor functions (inherited from Connection< >) to obtain delay and
  // target
  Node* target = get_target( t );
  double dendritic_delay = get_delay();

  target->compress_LTP_history( tau_x_, t_spike - dendritic_delay );
  // facilitation due to post-synaptic activity since last pre-synaptic spike
  //double old_w = weight_;
  weight_ = facilitate_( weight_, target->get_LTP_value( t_lastspike_ - dendritic_delay ), x_bar_ );
  //std::cout << "facilitation: " << weight_ - old_w << "  x_bar = " << x_bar_ << std::endl;

  // depression due to new pre-synaptic spike
  weight_ = depress_( weight_, target->get_LTD_value( t_spike - dendritic_delay ) );

  e.set_receiver( *target );
  e.set_weight( weight_ );
  // use accessor functions (inherited from Connection< >) to obtain delay in
  // steps and rport
  e.set_delay_steps( get_delay_steps() );
  e.set_rport( get_rport() );
  e();

  // compute the trace of the presynaptic spike train
  x_bar_ = x_bar_ * std::exp( ( t_lastspike_ - t_spike ) / tau_x_ ) + 1.0 / tau_x_;

  t_lastspike_ = t_spike;
}


template < typename targetidentifierT >
ClopathConnectionBC< targetidentifierT >::ClopathConnectionBC()
  : ConnectionBase()
  , weight_( 1.0 )
  , x_bar_( 0.0 )
  , tau_x_( 15.0 )
  , Wmin_( 0.0 )
  , Wmax_( 100.0 )
  , eta_( 1.0 )
  , t_lastspike_( 0.0 )
{
}

template < typename targetidentifierT >
ClopathConnectionBC< targetidentifierT >::ClopathConnectionBC( const ClopathConnectionBC< targetidentifierT >& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , x_bar_( rhs.x_bar_ )
  , tau_x_( rhs.tau_x_ )
  , Wmin_( rhs.Wmin_ )
  , Wmax_( rhs.Wmax_ )
  , eta_( rhs.eta_ )
  , t_lastspike_( rhs.t_lastspike_ )
{
}

template < typename targetidentifierT >
void
ClopathConnectionBC< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< double >( d, names::x_bar, x_bar_ );
  def< double >( d, names::tau_x, tau_x_ );
  def< double >( d, names::Wmin, Wmin_ );
  def< double >( d, names::Wmax, Wmax_ );
  def< double >( d, names::eta, eta_ );
  def< long >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
ClopathConnectionBC< targetidentifierT >::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
  updateValue< double >( d, names::x_bar, x_bar_ );
  updateValue< double >( d, names::tau_x, tau_x_ );
  updateValue< double >( d, names::Wmin, Wmin_ );
  updateValue< double >( d, names::Wmax, Wmax_ );
  updateValue< double >( d, names::eta, eta_ );

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

#endif // of #ifndef CLOPATH_CONNECTION_BC_H
