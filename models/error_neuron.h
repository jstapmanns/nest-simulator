/*
 *  error_neuron.h
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

#ifndef ERROR_NEURON_H
#define ERROR_NEURON_H

// Generated includes:
#include "config.h"

// C++ includes:
#include <string>

// Includes from nestkernel:
#include "eprop_archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "node.h"
#include "normal_randomdev.h"
#include "poisson_randomdev.h"
#include "ring_buffer.h"
// #include "recordables_map.h"
#include "universal_data_logger.h"

namespace nest
{

/** @BeginDocumentation
@ingroup Neurons
@ingroup rate

Name: error_neuron - Rate neuron that sums up incoming rates
                and applies a nonlinearity specified via the template.

Description:

The rate transformer node simply applies the nonlinearity specified in the
input-function of the template class to all incoming inputs.
An important application is to provide the possibility to
apply different nonlinearities to different incoming connections of the
same rate neuron by connecting the sending rate neurons to the
rate transformer node and connecting the rate transformer node to the
receiving rate neuron instead of using a direct connection.

Remarks:

- Weights on connections from and to the error_neuron_
  are handled as usual.
- Delays are honored on incoming and outgoing connections.

Receives: DelayedRateConnectionEvent

Sends: DelayedRateConnectionEvent

Parameters:

Only the parameters from the class Nonlinearities can be set in the
status dictionary.

Author: Mario Senden, Jan Hahne, Jannis Schuecker

FirstVersion: November 2017
*/
class error_neuron : public Eprop_Archiving_Node
{

public:
  typedef Node base;

  error_neuron();
  error_neuron( const error_neuron& );

  /**
   * Import sets of overloaded virtual functions.
   * We need to explicitly include sets of overloaded
   * virtual functions into the current scope.
   * According to the SUN C++ FAQ, this is the correct
   * way of doing things, although all other compilers
   * happily live without.
   */

  using Node::handle;
  using Node::sends_secondary_event;

  using Node::handles_test_event;

  void handle( DelayedRateConnectionEvent& );
  void handle( SpikeEvent& );
  void handle( CurrentEvent& );
  void handle( LearningSignalConnectionEvent& );
  void handle( DataLoggingRequest& );

  port handles_test_event( DelayedRateConnectionEvent&, rport );
  port handles_test_event( SpikeEvent&, rport );
  port handles_test_event( CurrentEvent&, rport );
  port handles_test_event( LearningSignalConnectionEvent&, rport );
  port handles_test_event( DataLoggingRequest&, rport );

  void
  sends_secondary_event( LearningSignalConnectionEvent& )
  {
  }


  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );
  bool is_eprop_readout();

private:
  void init_state_( const Node& proto );
  void init_buffers_();
  void calibrate();

  void update_( Time const&, const long, const long );

  void update( Time const&, const long, const long );
  double phi( double );

  // The next two classes need to be friends to access the State_ class/member
  friend class RecordablesMap< error_neuron >;
  friend class UniversalDataLogger< error_neuron >;

  // ----------------------------------------------------------------

  /**
   * Independent parameters of the model.
   */
  struct Parameters_
  {
    /** Target of non-linearity.
        True (default): Gain function applied to linearly summed input.
        False: Gain function applied to each input before summation.
    **/
    double phi_max_;    //!< Parameter of the rate function
    double rate_slope_; //!< Parameter of the rate function
    double beta_;       //!< Parameter of the rate function
    double theta_;      //!< Parameter of the rate function
    double tau_m_;  // Membrane time constant in ms.
    double c_m_;  // Membrane capacitance in pF.
    double E_L_;  // Resting potential in mV.
    double I_e_;  // External DC current.
    double V_min_;  // Lower bound relative to resting potential.
    double t_start_ls_; // time after which a learning signal is sent to the recurrent neurons
    bool regression_; // regression if true else classification
    Parameters_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary

    double set( const DictionaryDatum& );
  };

  // ----------------------------------------------------------------

  /**
   * State variables of the model.
   */
  struct State_
  {
    double target_rate_; //!< Rate
    double learning_signal_; //!< Rate
    double y0_;
    double y3_; // Membrane potential relative to resting potential.
    State_(); //!< Default initialization

    void get( DictionaryDatum&, const Parameters_&) const;

    /** Set values from dictionary.
     * @param dictionary to take data from
     * @param current parameters
     * @param Change in reversal potential E_L specified by this dict
     */
    void set( const DictionaryDatum&, const Parameters_&, double );
  };

  // ----------------------------------------------------------------

  /**
   * Buffers of the model.
   */
  struct Buffers_
  {
    Buffers_( error_neuron& );
    Buffers_( const Buffers_&, error_neuron& );

    // buffer for rate vector received by DelayRateConnection
    RingBuffer delayed_rates_;

    /** buffers and summs up incoming spikes/currents */
    RingBuffer spikes_;
    RingBuffer currents_;

    //! Logger for all analog data
    UniversalDataLogger< error_neuron > logger_;
  };

  // ----------------------------------------------------------------

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {

    double P30_;
    double P33_;
    int step_start_ls_; // step after which a learning signal is sent to the recurrent neurons
    int T_steps_; // length of one period T in steps
    bool prnt;
    librandom::RngPtr rng_;                   //!< random number generator of my own thread
    librandom::PoissonRandomDev poisson_dev_; //!< random deviate generator

   };

  void write_readout_history( Time const& t_sp,
    double readout_signal, double target_signal, double norm );

  void add_learning_to_hist( LearningSignalConnectionEvent& e );

  // DEBUG II:use this function to read learning signal in case of evidence accumulation task
  double
  get_last_ls_() const
  {
    if ( eprop_history_.size() > 2 )
    {
      return ( ( eprop_history_.rbegin() ) + 2 )->target_signal_
        - ( ( eprop_history_.rbegin() ) + 2 )->readout_signal_
        / ( ( eprop_history_.rbegin() ) + 2 )->normalization_ ;
    }
    return 0.0;
  }

  // DEBUG II: use this function to read learning signal in case of pattern generation task
  double
  get_learning_signal_() const
  {
    return S_.target_rate_ - (S_.y3_ + P_.E_L_);
  }

  //! Read out the rate
  double
  get_target_rate_() const
  {
    return S_.target_rate_;
  }

  //! Read out the real membrane potential
  double
  get_V_m_() const
  {
    return S_.y3_ + P_.E_L_;
  }

  // ----------------------------------------------------------------

  Parameters_ P_;
  State_ S_;
  Variables_ V_;
  Buffers_ B_;

  //! Mapping of recordables names to access functions
  static RecordablesMap< error_neuron > recordablesMap_;
};

inline void
error_neuron::update( Time const& origin,
  const long from,
  const long to )
{
  update_( origin, from, to);
}

inline double
error_neuron::phi( double u )
{
  return P_.phi_max_ / ( 1.0 + P_.rate_slope_ * exp( P_.beta_ * ( P_.theta_ - u ) ) );
}

inline port
error_neuron::handles_test_event(
  DelayedRateConnectionEvent&,
  rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}


inline port
error_neuron::handles_test_event(
        SpikeEvent&, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline port
error_neuron::handles_test_event(
        CurrentEvent&, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline port
error_neuron::handles_test_event( LearningSignalConnectionEvent&,
  rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline port
error_neuron::handles_test_event(
  DataLoggingRequest& dlr,
  rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}

inline void
error_neuron::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d, P_ );
  Eprop_Archiving_Node::get_status( d );
  ( *d )[ names::recordables ] = recordablesMap_.get_list();

}

inline void
error_neuron::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_; // temporary copy in case of errors
  const double delta_EL = ptmp.set( d ); // throws if BadProperty
  State_ stmp = S_;      // temporary copy in case of errors
  stmp.set( d, ptmp, delta_EL );         // throws if BadProperty

  // We now know that (stmp) is consistent. We do not
  // write it back to (S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  Eprop_Archiving_Node::set_status( d );

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
  S_ = stmp;

}

} // namespace

#endif /* #ifndef ERROR_NEURON_H */
