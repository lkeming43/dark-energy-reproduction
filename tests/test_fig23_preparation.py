"""Independent checks of conversion, source covariance conventions and row order."""
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from src.data_loading import load_scenario_config, load_supernova, prepare_combined_dataset
from src.distance_conversion import convert_supernova_to_dm_over_rd, propagate_covariance

ROOT = Path(__file__).resolve().parents[1]


class PreparationTests(unittest.TestCase):
    def test_distance_and_jacobian(self):
        # mu=35 -> luminosity distance 100 Mpc; at z=1 and rd=50 -> y=1.
        sn = {'z': np.array([1.,1.]), 'observable':np.array([35.,35.]),
              'observable_kind':'distance_modulus'}
        cal = {'rd_mpc':50.}
        np.testing.assert_allclose(convert_supernova_to_dm_over_rd(sn,cal),1.)
        cov = np.array([[.04,.01],[.01,.09]])
        eps = 1e-5
        plus = dict(sn,observable=sn['observable']+eps)
        minus = dict(sn,observable=sn['observable']-eps)
        derivative = (convert_supernova_to_dm_over_rd(plus,cal)-convert_supernova_to_dm_over_rd(minus,cal))/(2*eps)
        np.testing.assert_allclose(propagate_covariance(sn,cov,cal),
                                   cov*derivative[:,None]*derivative[None,:],rtol=1e-8)

    def test_release_covariance_and_selection(self):
        for name,expected in [('fig2_pantheon',1590),('fig3_desy5',1754)]:
            config = load_scenario_config(ROOT/'configs'/f'{name}.yaml')
            cfg = config['supernova']
            frame = pd.read_csv(ROOT/cfg['data_file'],sep=cfg['separator'])
            original_z = frame.zHD.to_numpy().copy()
            sn = load_supernova(config)
            self.assertEqual(len(sn['z']),expected)
            np.testing.assert_array_equal(sn['z'],original_z[sn['row_index']])
            raw = np.loadtxt(ROOT/cfg['covariance_file'],skiprows=1).reshape(len(frame),len(frame))
            expected_diag = np.diag(raw).copy()
            if name == 'fig3_desy5':
                expected_diag += frame.MUERR_FINAL.to_numpy()**2
            np.testing.assert_allclose(np.diag(sn['covariance']),expected_diag[sn['row_index']])
            np.testing.assert_allclose(sn['covariance'],((raw+raw.T)/2)[np.ix_(sn['row_index'],sn['row_index'])]
                                       + (np.diag(frame.MUERR_FINAL.to_numpy()[sn['row_index']]**2)
                                          if name == 'fig3_desy5' else 0))

    def test_combined_order_and_raw_data_integrity(self):
        for name,total in [('fig2_pantheon',1596),('fig3_desy5',1760)]:
            config = load_scenario_config(ROOT/'configs'/f'{name}.yaml')
            sn = load_supernova(config)
            before = sn['z'].copy()
            data = prepare_combined_dataset(config)
            self.assertEqual(len(data['z']),total)
            np.testing.assert_array_equal(data['z'][6:],before)
            self.assertTrue(np.all(data['covariance'][:6,6:]==0))
            np.testing.assert_array_equal(data['sn_observable'],sn['observable'])


if __name__ == '__main__':
    unittest.main()
