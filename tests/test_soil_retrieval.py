import unittest
import numpy as np
from navasar import retrieve_smi, event_response, closure_phase_from_ifgrams
from navasar.dielectric_inversion import dielectric_covariance, invert_closure_moisture, MODEL_NOTES


class SoilRetrievalTests(unittest.TestCase):
    def test_public_event_api(self):
        t=np.arange(60.)*12
        closure,_=event_response(t,[8.],[40.],[120.])
        result=retrieve_smi(t,closure,median_window=1)
        self.assertEqual(result['event_count'],1)
        self.assertLess(result['closure_rmse_deg'],1e-5)
        self.assertTrue(np.isnan(retrieve_smi(t,t*0)['event_smi']).all())

    def test_network_api(self):
        phase=np.linspace(0,1,10)
        cov=np.exp(1j*(phase[None,:]-phase[:,None]))
        np.testing.assert_allclose(closure_phase_from_ifgrams(cov)[0],0,atol=1e-10)

    def test_opt_in_and_frequency(self):
        kw=dict(model='hallikainen',frequency_ghz=5.,sand=.51,clay=.13)
        with self.assertRaises(ValueError): dielectric_covariance([.1,.2,.15],**kw)
        kw['frequency_ghz']=5.405
        self.assertTrue(np.isfinite(dielectric_covariance([.1,.2,.15],allow_experimental=True,**kw)).all())
        kw['frequency_ghz']=19.
        with self.assertRaises(ValueError): dielectric_covariance([.1,.2,.15],allow_experimental=True,**kw)

    def test_all_local_models_covariance_and_inverse(self):
        moisture=np.array([.10,.14,.22,.18,.15,.12,.11,.10])
        for model in MODEL_NOTES:
            with self.subTest(model=model):
                kw=dict(model=model,frequency_ghz=5.,sand=.51,clay=.13,allow_experimental=True)
                cov=dielectric_covariance(moisture,**kw)
                np.testing.assert_allclose(cov,cov.conj().T,atol=1e-12)
                np.testing.assert_allclose(np.diag(cov),1,atol=1e-12)
                self.assertGreater(np.linalg.eigvalsh(cov).min(),-1e-10)
                result=invert_closure_moisture(cov,np.arange(8.)*12,reference_moisture=.1,
                                              starts=12,max_nfev=300,**kw)
                self.assertEqual(result['moisture'][0],.1)
                self.assertTrue(np.all(np.isfinite(result['moisture'])))
                self.assertTrue(np.all((result['moisture']>=.005)&(result['moisture']<=.45)))
                self.assertLess(result['closure_rmse_deg'],1.)
                self.assertFalse(result['absolute_moisture_validated'])
                self.assertIn('multistart_max_spread',result)


if __name__=='__main__': unittest.main()
