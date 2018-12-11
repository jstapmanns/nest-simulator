# -*- coding: utf-8 -*-
#
# test_clopath_stdp_synapse.py
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.

"""
Test functionality of the Clopath stdp synapse
"""

import unittest
import nest
import numpy as np

HAVE_GSL = nest.sli_func("statusdict/have_gsl ::")


@nest.check_stack
@unittest.skipIf(not HAVE_GSL, 'GSL is not available')
class ClopathSynapseTestCase(unittest.TestCase):
    """Test Clopath stdp synapse"""

    def test_ConnectNeuronsWithClopathSynapse(self):
        """Ensures that the Clopath stdp synapse can only be used with
        supported neuron models."""

        supported_models = [
            'aeif_cbvg_2010',
            'hh_psc_alpha_clopath',
        ]

        for nm in supported_models:
            nest.ResetKernel()

            n = nest.Create(nm,2)

            nest.Connect(n, n, {"rule": "all_to_all"},
                    {"model": "clopath_stdp_synapse"})

        not_supported_models = [n for n in nest.Models(mtype='nodes') 
                               if n not in supported_models]

        for nm in not_supported_models:
            nest.ResetKernel()

            n = nest.Create(nm,2)

            # try to connect with clopath_rule
            with self.assertRaises(nest.NESTError):
                nest.Connect(n, n, {"rule": "all_to_all"},
                         {"model": "clopath_stdp_synapse"})



    def test_SynapseFunctionWithAeifModel(self):
        """Ensure that spikes are properly processed by the
        aeif_cbvg_2010 model and that spikes send through the
        clopath_stdp_synapse are received by the target neuron"""

        nest.ResetKernel()

        # Create neurons and devices
        nrns = nest.Create('aeif_cbvg_2010', 2, {'V_m':-70.6})
        prrt_nrn = nest.Create('parrot_neuron', 1)

        spike_times = [10.0]
        sg = nest.Create('spike_generator', 1, {'spike_times':spike_times})

        mm = nest.Create('multimeter', params={'record_from' : ['V_m'], 'interval':1.0})

        # Connect spike generator to parrot neurons
        nest.Connect(sg, prrt_nrn)

        # static connection
        conn_dict = {'rule':'all_to_all'}
        static_syn_dict = {'model':'static_synapse', 'weight':2.0, 'delay': 1.0}
        nest.Connect(prrt_nrn, nrns[0:1], conn_dict, static_syn_dict)

        # Clopath stdp connection
        cl_stdp_syn_dict = {'model':'clopath_stdp_synapse', 'weight':2.0, 'delay': 1.0}
        nest.Connect(prrt_nrn, nrns[1:2], conn_dict, cl_stdp_syn_dict)

        # Connect multimeter
        nest.Connect(mm, nrns)

        # Simulation using nest
        nest.Simulate(20.)
        
        # Read out
        data = nest.GetStatus(mm)
        senders = data[0]['events']['senders']
        voltages = data[0]['events']['V_m']
        
        vm1 = voltages[np.where(senders==1)]
        vm2 = voltages[np.where(senders==1)]
        
        self.assertTrue(np.allclose(vm1,vm2,rtol=1e-5))
        self.assertTrue(np.isclose(vm2[11]-vm2[10],2,rtol=1e-5))


def suite():
    suite = unittest.makeSuite(ClopathSynapseTestCase, 'test')
    return suite


def run():
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())


if __name__ == "__main__":
    run()
