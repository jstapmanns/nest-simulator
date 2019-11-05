/*
 *  error_transformer_node.h
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

#ifndef ERROR_TRANSFORMER_NODE_H
#define ERROR_TRANSFORMER_NODE_H

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
#include "recordables_map.h"
#include "universal_data_logger.h"

namespace nest
{

/** @BeginDocumentation
@ingroup Neurons
@ingroup rate

Name: error_transformer_node - Rate neuron that sums up incoming rates
                and applies a nonlinearity specified via the template.

Description:

The rate transformer node simply applies the nonlinearity specified in the
input-function of the template class to all incoming inputs. The boolean
parameter linear_summation determines whether the input function is applied to
the summed up incoming connections (True, default value) or to each input
individually (False).
An important application is to provide the possibility to
apply different nonlinearities to different incoming connections of the
same rate neuron by connecting the sending rate neurons to the
rate transformer node and connecting the rate transformer node to the
receiving rate neuron instead of using a direct connection.

Remarks:

- Weights on connections from and to the error_transformer_node
  are handled as usual.
- Delays are honored on incoming and outgoing connections.

Receives: DelayedRateConnectionEvent

Sends: DelayedRateConnectionEvent

Parameters:

Only the parameter
- linear_summation
and the parameters from the class Nonlinearities can be set in the
status dictionary.

Author: Mario Senden, Jan Hahne, Jannis Schuecker

FirstVersion: November 2017
*/
template < class TNonlinearities >
class error_transformer_node : public Eprop_Archiving_Node
{

public:
  typedef Node base;

  error_transformer_node();
  error_transformer_node( const error_transformer_node& );

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
  void handle( DataLoggingRequest& );

  port handles_test_event( DelayedRateConnectionEvent&, rport );
  port handles_test_event( SpikeEvent&, rport );
  port handles_test_event( CurrentEvent&, rport );
  port handles_test_event( DataLoggingRequest&, rport );

  void
  sends_secondary_event( DelayedRateConnectionEvent& )
  {
  }


  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );
  bool is_eprop_readout();

private:
  void init_state_( const Node& proto );
  void init_buffers_();
  void calibrate();

  TNonlinearities nonlinearities_;

  void update_( Time const&, const long, const long );

  void update( Time const&, const long, const long );

  // The next two classes need to be friends to access the State_ class/member
  friend class RecordablesMap< error_transformer_node< TNonlinearities > >;
  friend class UniversalDataLogger< error_transformer_node< TNonlinearities > >;

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
    bool linear_summation_;

    double tau_m_;  // Membrane time constant in ms.
    double c_m_;  // Membrane capacitance in pF.
    double E_L_;  // Resting potential in mV.
    double I_e_;  // External DC current.
    double V_min_;  // Lower bound relative to resting potential.

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
    double rate_; //!< Rate
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
    Buffers_( error_transformer_node& );
    Buffers_( const Buffers_&, error_transformer_node& );

    // buffer for rate vector received by DelayRateConnection
    RingBuffer delayed_rates_;

    /** buffers and summs up incoming spikes/currents */
    RingBuffer spikes_;
    RingBuffer currents_;

    //! Logger for all analog data
    UniversalDataLogger< error_transformer_node > logger_;
  };

  // ----------------------------------------------------------------

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {

    double P30_;
    double P33_;

   };

  //! Read out the rate
  double
  get_rate_() const
  {
    return S_.rate_;
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
  static RecordablesMap< error_transformer_node< TNonlinearities > >
    recordablesMap_;
};

template < class TNonlinearities >
inline void
error_transformer_node< TNonlinearities >::update( Time const& origin,
  const long from,
  const long to )
{
  update_( origin, from, to);
}

template < class TNonlinearities >
inline port
error_transformer_node< TNonlinearities >::handles_test_event(
  DelayedRateConnectionEvent&,
  rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}


template < class TNonlinearities >
inline port
error_transformer_node< TNonlinearities >::handles_test_event(
        SpikeEvent&, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

template < class TNonlinearities >
inline port
error_transformer_node< TNonlinearities >::handles_test_event(
        CurrentEvent&, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}


template < class TNonlinearities >
inline port
error_transformer_node< TNonlinearities >::handles_test_event(
  DataLoggingRequest& dlr,
  rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}

template < class TNonlinearities >
inline void
error_transformer_node< TNonlinearities >::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d, P_ );
  Archiving_Node::get_status( d );
  ( *d )[ names::recordables ] = recordablesMap_.get_list();

  nonlinearities_.get( d );
}

template < class TNonlinearities >
inline void
error_transformer_node< TNonlinearities >::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_; // temporary copy in case of errors
  const double delta_EL = ptmp.set( d ); // throws if BadProperty
  State_ stmp = S_;      // temporary copy in case of errors
  stmp.set( d, ptmp, delta_EL );         // throws if BadProperty

  // We now know that (stmp) is consistent. We do not
  // write it back to (S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  Archiving_Node::set_status( d );

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
  S_ = stmp;

  nonlinearities_.set( d );
}

} // namespace

#endif /* #ifndef ERROR_TRANSFORMER_NODE_H */
