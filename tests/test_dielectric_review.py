import unittest
import numpy as np
from navasar.dielectric import hallikainen_dielectric
from navasar.dobson_reviewed import dobson_ulaby

class ReviewedDielectricTests(unittest.TestCase):
    def test_hand_calculated_table_rows(self):
        # Independently evaluated Table E.1 at S=51%, C=13%, mv=0.2.
        self.assertAlmostEqual(hallikainen_dielectric(.2,4.,.51,.13).real,10.86344)
        self.assertAlmostEqual(hallikainen_dielectric(.2,4.,.51,.13).imag,1.50904)
        self.assertAlmostEqual(hallikainen_dielectric(.2,6.,.51,.13).real,10.24868)
        self.assertAlmostEqual(hallikainen_dielectric(.2,6.,.51,.13).imag,1.94732)

    def test_adjacent_interpolation_and_texture_maps(self):
        lo=hallikainen_dielectric(.2,4.);hi=hallikainen_dielectric(.2,6.)
        np.testing.assert_allclose(hallikainen_dielectric(.2,5.405),lo+(hi-lo)*1.405/2)
        s=np.array([[.51,.4],[.3,.2]])
        v=hallikainen_dielectric(.2,6.,s,.13)
        for ij in np.ndindex(s.shape):
            self.assertAlmostEqual(v[ij],hallikainen_dielectric(.2,6.,s[ij],.13))
        with self.assertRaises(ValueError): hallikainen_dielectric(.2,19.)
        with self.assertRaises(ValueError): hallikainen_dielectric(.2,6.,51,13)

    def test_dobson_texture_and_dry_limit(self):
        a=dobson_ulaby(.2,5.405,.51,.13)
        b=dobson_ulaby(.2,5.405,.3,.3)
        self.assertGreater(abs(a-b),.1)
        self.assertEqual(dobson_ulaby(0.,5.405).imag,0)
        self.assertAlmostEqual(dobson_ulaby(0.,5.405).real,(1+.66*1.3)**(1/.65))
        with self.assertRaises(ValueError): dobson_ulaby(.6,5.405)

if __name__=='__main__': unittest.main()
