import os
from unittest import TestCase
from floes.yank_binding import job as floe


class YankBindingTestCase(TestCase):
    """ Invocation: python -m pytest floes/floe_tests/ -s """

    def test_simple_run(self):
        receptor_path = os.path.join(os.path.dirname(__file__), '..', '..', 'examples', 'data', 'T4-protein.pdb')
        ligand_path = os.path.join(os.path.dirname(__file__), '..', '..', 'examples', 'data', 'p-xylene.mol2')
        # run_args should be a list like sys.argv
        run_args = [
            '--receptor', receptor_path,
            '--molecules', ligand_path,
            '--success', 'success.sdf',
            '--failure', 'failure.sdf',
            '--simulation_time', '0.0001',
            '--nsteps_per_iteration', '5',
            '--minimize', '0',
        ]
        self.assertFalse(floe.run(args=run_args))
