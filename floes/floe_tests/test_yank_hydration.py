import os
from unittest import TestCase
from floes.yank_hydration import job as floe


class YankHydrationTestCase(TestCase):
    """ Invocation: python -m pytest floes/floe_tests/ -s """

    def test_simple_run(self):
        molecules_path = os.path.join(os.path.dirname(__file__), '..', '..', 'examples', 'data', 'freesolv_mini.oeb.gz')
        # run_args should be a list like sys.argv
        run_args = [
            '--molecules', molecules_path,
            '--success', 'success.sdf',
            '--failure', 'failure.sdf',
            '--simulation_time', '0.001',
            '--nsteps_per_iteration', '50',
        ]
        self.assertFalse(floe.run(args=run_args))
