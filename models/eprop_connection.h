/*
 *  eprop_connection.h
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

#ifndef EPROP_CONNECTION_H
#define EPROP_CONNECTION_H

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
Name: eprop_connection - Synapse type for ... TODO

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

tau_x    double - Time constant of the trace of the presynaptic spike train.
Wmax     double - Maximum allowed weight.
Wmin     double - Minimum allowed weight.

Other parameters like the amplitudes for depression and facilitation are
stored in in the neuron models that are compatible with the Clopath synapse.

Transmits: SpikeEvent

References:

[1] Clopath et al. (2010) Connectivity reflects coding:
    a model of voltage-based STDP with homeostasis.
    Nature Neuroscience 13:3, 344--352
[2] Clopath and Gerstner (2010) Voltage and spike timing interact
    in STDP – a unified model. Front. Synaptic Neurosci. 2:25
    doi: 10.3389/fnsyn.2010.00025
[3] Voltage-based STDP synapse (Clopath et al. 2010) on ModelDB
    https://senselab.med.yale.edu/ModelDB/showmodel.cshtml?model=144566

Authors: Jonas Stapmanns, David Dahmen, Jan Hahne

SeeAlso: stdp_synapse, aeif_psc_delta_clopath, hh_psc_alpha_clopath
*/
// connections are templates of target identifier type (used for pointer /
// target index addressing) derived from generic connection template
template < typename targetidentifierT >
class EpropConnection : public Connection< targetidentifierT >
{

public:
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  EpropConnection();


  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  EpropConnection( const EpropConnection& );

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

    if ( Eprop_Archiving_Node* t_eprop = dynamic_cast< Eprop_Archiving_Node* >( &t ) )
    {
      //std::cout << "register eprop node" << std::endl;
      if (t_eprop->is_eprop_readout() )  // if target is a readout neuron
      {
        t_eprop->init_eprop_buffers( 3.0 * get_delay() );
      }
      else
      {
        t_eprop->init_eprop_buffers( 2.0 * get_delay() );
      }
    }
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
  // TODO: tau_alpha_ can be determined from neuron params
  double tau_alpha_; // time constant corresponding to leak term of rec neurons
  double tau_kappa_; // time constant corresponding to leak term of output neurons
  double learning_rate_;
  double update_interval_;
  double Wmin_;
  double Wmax_;

  double t_lastspike_;
  double t_lastupdate_;
  double t_nextupdate_;
  double last_e_trace_;
  double t_prime_int_trace_;
  double keep_traces_;
  double rate_reg_;
  double target_firing_rate_;

  std::vector< double > pre_syn_spike_times_;
};


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param t The thread on which this connection is stored.
 * \param cp Common properties object, containing the stdp parameters.
 */
template < typename targetidentifierT >
inline void
EpropConnection< targetidentifierT >::send( Event& e,
  thread t,
  const CommonSynapseProperties& )
{
  double t_spike = e.get_stamp().get_ms();
  // use accessor functions (inherited from Connection< >) to obtain delay and
  // target
  Node* target = get_target( t );
  double dendritic_delay = get_delay();

  // store times of incoming spikes to enable computation of eligibility trace
  pre_syn_spike_times_.push_back( t_spike );

  // do update only 
  if ( t_spike > t_nextupdate_ )
  {
    if ( keep_traces_ < 1.0 )
    {
      last_e_trace_ = 0.0;
      t_prime_int_trace_ = 0.0;
    }
    // get spike history in relevant range (t1, t2] from post-synaptic neuron
    std::deque< histentry_eprop >::iterator start;
    std::deque< histentry_eprop >::iterator finish;

    std::deque< double >::iterator start_spk;
    std::deque< double >::iterator finish_spk;
    // For a new synapse, t_lastspike_ contains the point in time of the last
    // spike. So we initially read the
    // history(t_last_spike - dendritic_delay, ..., T_spike-dendritic_delay]
    // which increases the access counter for these entries.
    // At registration, all entries' access counters of
    // history[0, ..., t_last_spike - dendritic_delay] have been
    // incremented by Archiving_Node::register_stdp_connection(). See bug #218 for
    // details.
    //DEBUG: added 2*delay to be in sync with TF code
    double t_update_ = ( floor( t_spike / update_interval_ ) ) * update_interval_ + 2.0 *
      dendritic_delay;
    //std::cout << "t_update = " << t_update_ << " update interval = " << update_interval_ << std::endl;
    // double t1 = std::max( 0.0, ( floor( t_lastspike_ / update_interval_ ) ) * update_interval_ );
    // double t2 = t1 + update_interval_;

    //std::cout << "in synapse at time: " << t_spike << ", t_lu: " << t_lastupdate_
     //<< ", t_u: " << t_update_ << std::endl;
    if ( start != finish )
    {
      ++start;
    }
    // target->get_eprop_history( t1 - dendritic_delay,
    //     t2 - dendritic_delay,
    //     &start,
    //     &finish );

    target->get_spike_history( t_lastupdate_ - dendritic_delay,
        t_update_ - dendritic_delay,
        &start_spk,
        &finish_spk );
    int nspikes = std::distance(start_spk, finish_spk);
    // std::cout << "nspikes " << nspikes << std::endl;
    double const dt = Time::get_resolution().get_ms();
    double kappa = std::exp( -dt / tau_kappa_ );
    std::vector< double >::iterator t_pre_spike = pre_syn_spike_times_.begin();
    double dw = 0.0;
    if (target->is_eprop_readout() )  // if target is a readout neuron
    {
      target->get_eprop_history( t_lastupdate_ + dendritic_delay,
          t_update_ + dendritic_delay,
          &start,
          &finish );

      /*
      if ( Eprop_Archiving_Node* t_eprop = dynamic_cast< Eprop_Archiving_Node* >( target ) )
      {
        std::cout << "update readout" << std::endl;
        t_eprop->print_eprop_history();
      }
      */
      //std::cout << "I'm a readout neuron" << std::endl;
      //std::cout << "lu + delay = " << t_lastupdate_ + dendritic_delay << "  t_update + delay = " <<
        //t_update_ + dendritic_delay << std::endl;
      //std::cout << "trace: t, z_hat, ls:" << std::endl;
      while ( start != finish )
      {
        last_e_trace_ *= kappa;
        if ( std::fabs( *t_pre_spike - start->t_ + dendritic_delay ) < 1.0e-6 )
        {
          // DEBUG: inserted factor ( 1 - dacay )
          last_e_trace_ += ( 1.0 - kappa );
          t_pre_spike++;
        }
        //std::cout << start->t_ << " " << last_e_trace_ << " " << start->learning_signal_ << "; ";
        dw += start->learning_signal_ * last_e_trace_;
        start++;
      }
     dw *= learning_rate_ * dt;
     //std::cout << std::endl << "dw out = " << dw << std::endl;
    }
    else  // if target is a neuron of the recurrent network
    {
      target->get_eprop_history( t_lastupdate_ - 0.0 * dendritic_delay,
          t_update_ - 0.0 * dendritic_delay,
          &start,
          &finish );

      std::vector< double > elegibility_trace;
      double alpha = target->get_leak_propagator();
      // compute the sum of the elegibility trace because it is used for the firing rate
      // regularization
      double sum_eleg_tr = 0.0;
      if ( target->is_eprop_adaptive() )
      {
        // if the target is of type aif_psc_delta_eprop (adaptive threshold)
        double beta = target->get_beta();
        double rho = target->get_adapt_propagator();
        double epsilon = 0.0;
        for ( std::deque< histentry_eprop >::iterator runner = start; runner != finish; runner++ )
        {
          double pseudo_deriv = runner->V_m_;
          // Eq.(22)
          last_e_trace_ *= alpha;
          // Eq.(27)
          epsilon = pseudo_deriv * last_e_trace_ + ( rho - beta * pseudo_deriv ) * epsilon;
          if ( std::fabs( *t_pre_spike - runner->t_ ) < 1.0e-6 )
          {
            // DEBUG: inserted factor ( 1 - dacay )
            last_e_trace_ += ( 1.0 - alpha );
            t_pre_spike++;
          }
          double eleg_tr = pseudo_deriv * ( last_e_trace_  - beta * epsilon );
          sum_eleg_tr += eleg_tr;
          // Eq.(28)
          elegibility_trace.push_back( eleg_tr );
        }
      }
      else
      {
        // if the target is of type iaf_psc_delta_eprop
        /*
        std::cout << "pre syn spike times:" << std::endl;
        for ( std::vector< double >::iterator sp_it = pre_syn_spike_times_.begin(); sp_it !=
            pre_syn_spike_times_.end(); sp_it++ )
        {
          std::cout << *sp_it << ", ";
        }
        std::cout << std::endl;
        */
        /*
        if ( Eprop_Archiving_Node* t_eprop = dynamic_cast< Eprop_Archiving_Node* >( target ) )
        {
          std::cout << "update recurrent" << std::endl;
          t_eprop->print_eprop_history();
        }
        */
        //std::cout << "t, ls, h, z_hat: " << std::endl;
        for ( std::deque< histentry_eprop >::iterator runner = start; runner != finish; runner++ )
        {
          // Eq.(22)
          last_e_trace_ *= alpha;
          //DEBUG: added dendritic delay for sync with TF code
          if ( std::fabs( *t_pre_spike - runner->t_ + dendritic_delay ) < 1.0e-6 )
          {
            // DEBUG: inserted factor ( 1 - dacay )
            last_e_trace_ += ( 1.0 - alpha );
            t_pre_spike++;
          }
          //std::cout << runner->t_ << " " << runner->learning_signal_ << " " << runner->V_m_ << " " << last_e_trace_ << ";  ";
          double eleg_tr = runner->V_m_ * last_e_trace_;
          sum_eleg_tr += eleg_tr;
          // Eq.(23)
          elegibility_trace.push_back( eleg_tr );
        }
        //std::cout << std::endl;
      }

      /*
      std::cout << "elegibility trace (" << elegibility_trace.size() << "): " << std::endl;
      for ( std::vector< double >::iterator it = elegibility_trace.begin(); it !=
      elegibility_trace.end(); it++)
      {
      std::cout << *it << ", ";
      }
      std::cout << std::endl;
      */
      /*
      std::cout  << "learning_signal: " << std::endl;;
      for ( std::deque< histentry_eprop >::iterator runner = start; runner != finish; runner++ )
      {
        std::cout << runner->learning_signal_ << " ";
      }
      std::cout << std::endl;
      */

      int t_prime = 0;
      //double sum_t_prime = 0.0;
      double sum_t_prime_new = 0.0;
      // std::cout << "dw(t): " << std::endl;
      //std::cout << "learning signal:" << std::endl;
      while ( start != finish )
      {
        // DEBUG: inserted factor ( 1 - decay )
        sum_t_prime_new = kappa * sum_t_prime_new + ( 1.0 - kappa ) * elegibility_trace[ t_prime ];
        dw += ( sum_t_prime_new * dt + std::pow( kappa, t_prime ) * t_prime_int_trace_ ) * start->learning_signal_;
        //std::cout << start->learning_signal_ << ",  ";
        // std::cout << dw*dt + weight_ << " ";
        //std::cout << start->V_m_ << ", ";
        t_prime++;
        start++;
      }
      //std::cout << std::endl;
      // firing rate regularization
      // compute average firing rate since last update. factor 1000 to convert into Hz
      double av_firing_rate = 1000.0 * nspikes / (t_update_ - t_lastupdate_);
      // Eq.(56)
      dw += rate_reg_ * ( target_firing_rate_ - av_firing_rate ) * sum_eleg_tr;
      /*
      std::cout << "target f_rate = " << target_firing_rate_ << "  actual f_rate = " <<
        av_firing_rate << "  dw_f_rate = " << rate_reg_ * ( target_firing_rate_ - av_firing_rate ) *
        sum_eleg_tr << std::endl;
        */
      dw *= dt*learning_rate_;
      t_prime_int_trace_ += sum_t_prime_new * dt;
      //std::cout << "dw rec/in: " << dw << std::endl;
    }

    weight_ += dw;
    //std::cout << "dw: " << dw << std::endl;
    /*
    if ( ! target->is_eprop_readout() )  // if target is a readout neuron
    {
      std::cout << "weight: " << weight_ << std::endl;
    }
    */

    // TODO: keep this in final implementation
    /*
    if ( weight_ > Wmax_ )
    {
      weight_ = Wmax_;
    }
    else if ( weight_ < Wmin_ )
    {
      weight_ = Wmin_;
    }
    */
    t_lastupdate_ = t_update_;
    t_nextupdate_ += ( floor( ( t_spike - t_nextupdate_ ) / update_interval_ ) + 1 ) *
      update_interval_;
    // clear history of presynaptic spike because we don't need them any more
    pre_syn_spike_times_.clear();
    pre_syn_spike_times_.push_back( t_spike );
    target->tidy_eprop_history( t_lastupdate_ - dendritic_delay );
    target->tidy_spike_history( t_lastupdate_ - dendritic_delay );
  }

  e.set_receiver( *target );
  e.set_weight( weight_ );
  // use accessor functions (inherited from Connection< >) to obtain delay in
  // steps and rport
  e.set_delay_steps( get_delay_steps() );
  e.set_rport( get_rport() );
  e();

  t_lastspike_ = t_spike;
}


template < typename targetidentifierT >
EpropConnection< targetidentifierT >::EpropConnection()
  : ConnectionBase()
  , weight_( 1.0 )
  , tau_alpha_( 10.0 )
  , tau_kappa_( 10.0 )
  , learning_rate_( 0.0001 )
  , update_interval_( 100.0 )
  , Wmin_( 0.0 )
  , Wmax_( 100.0 )
  , t_lastspike_( 0.0 )
  , t_lastupdate_( 0.0 )
  , t_nextupdate_( 100.0 )
  , last_e_trace_( 0.0 )
  , t_prime_int_trace_( 0.0 )
  , keep_traces_( true )
  , rate_reg_( 0 )
  , target_firing_rate_( 0.01 )
{
}

template < typename targetidentifierT >
EpropConnection< targetidentifierT >::EpropConnection(
  const EpropConnection< targetidentifierT >& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , tau_alpha_( rhs.tau_alpha_ )
  , tau_kappa_( rhs.tau_kappa_ )
  , learning_rate_( rhs.learning_rate_ )
  , update_interval_( rhs.update_interval_ )
  , Wmin_( rhs.Wmin_ )
  , Wmax_( rhs.Wmax_ )
  , t_lastspike_( rhs.t_lastspike_ )
  , t_lastupdate_( rhs.t_lastupdate_ )
  , t_nextupdate_( rhs.t_nextupdate_ )
  , last_e_trace_( rhs.last_e_trace_ )
  , t_prime_int_trace_( rhs.t_prime_int_trace_ )
  , keep_traces_( rhs.keep_traces_ )
  , rate_reg_( rhs.rate_reg_ )
  , target_firing_rate_( rhs.target_firing_rate_ )
{
}

template < typename targetidentifierT >
void
EpropConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< double >( d, names::tau_alpha, tau_alpha_ );
  def< double >( d, names::tau_kappa, tau_kappa_ );
  def< double >( d, names::learning_rate, learning_rate_ );
  def< double >( d, names::update_interval, update_interval_ );
  def< double >( d, names::Wmin, Wmin_ );
  def< double >( d, names::Wmax, Wmax_ );
  def< double >( d, names::keep_traces, keep_traces_ );
  def< double >( d, names::rate_reg, rate_reg_ );
  def< double >( d, names::target_firing_rate, target_firing_rate_ );
  def< long >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
EpropConnection< targetidentifierT >::set_status( const DictionaryDatum& d,
  ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
  updateValue< double >( d, names::tau_alpha, tau_alpha_ );
  updateValue< double >( d, names::tau_kappa, tau_kappa_ );
  updateValue< double >( d, names::learning_rate, learning_rate_ );
  updateValue< double >( d, names::update_interval, update_interval_ );
  updateValue< double >( d, names::Wmin, Wmin_ );
  updateValue< double >( d, names::Wmax, Wmax_ );
  updateValue< double >( d, names::keep_traces, keep_traces_ );
  updateValue< double >( d, names::rate_reg, rate_reg_ );
  updateValue< double >( d, names::target_firing_rate, target_firing_rate_ );

  // TODO: t_nextupdate and t_lastupdate should be initialized even if set_status is not called
  // DEBUG: added + delay to correct for the delay of the learning signal
  t_nextupdate_ = update_interval_ + 2.0 * get_delay(); // TODO: is this waht we want?
  //DEBUG: shifted initial value of t_lastupdate to be in sync with TF code
  t_lastupdate_ = 2.0 * get_delay();

  // check if weight_ and Wmin_ has the same sign
  if ( not( ( ( weight_ >= 0 ) - ( weight_ < 0 ) )
         == ( ( Wmin_ >= 0 ) - ( Wmin_ < 0 ) ) ) )
  {
    // throw BadProperty( "Weight and Wmin must have same sign." );
  }

  // check if weight_ and Wmax_ has the same sign
  if ( not( ( ( weight_ >= 0 ) - ( weight_ < 0 ) )
         == ( ( Wmax_ > 0 ) - ( Wmax_ <= 0 ) ) ) )
  {
    // throw BadProperty( "Weight and Wmax must have same sign." );
  }

  if ( tau_alpha_ <= 0.0 )
  {
    throw BadProperty( "The synaptic time constant tau must be greater than zero." );
  }

  if ( update_interval_ <= 0.0 )
  {
    throw BadProperty( "The synaptic update interval must be greater than zero." );
  }
}

} // of namespace nest

#endif // of #ifndef URBANCZIK_CONNECTION_H
