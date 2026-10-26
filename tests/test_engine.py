"""Scientific invariants and FDR-family checks for the consolidated implementation."""
from pathlib import Path
from unittest.mock import patch
import subprocess,sys,unittest
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'common/02_oinformation'))
import compute_hoi_oinfo as engine
class EngineChecks(unittest.TestCase):
    def test_xor_and_redundancy_sign_and_local_identity(self):
        pairs=np.tile(np.array([[0,0],[0,1],[1,0],[1,1]]),(32,1))
        xor=np.column_stack([pairs,pairs[:,0]^pairs[:,1]])
        red=np.tile(np.array([[0,0,0],[1,1,1]]),(64,1))
        for matrix,expected in [(xor,-1.),(red,1.)]:
            for n_states in [2,3,5]:
                result=engine.information_from_codes(matrix,n_states)
                local=engine.local_oinfo_from_codes(matrix,n_states)
                self.assertAlmostEqual(result['hoi'],expected,places=8)
                self.assertAlmostEqual(local.mean(),result['hoi'],places=8)
    def test_bh_known_family(self):
        np.testing.assert_allclose(engine.benjamini_hochberg([.01,.04,.03,.9]),[.04,.16/3,.16/3,.9])
        with self.assertRaises(ValueError):engine.benjamini_hochberg([.2,np.nan])
    def test_fdr_before_sign_selection_and_family_boundaries(self):
        frame=pd.DataFrame(dict(sample=['a']*4,block=['x','x','y','y'],order=[3]*4,
            p_down=[.01,.9,.03,.8],observed_hoi=[-1,1,-1,1]))
        result=engine.correct_triplet_fdr(frame,'sample-block')
        np.testing.assert_allclose(result.q_down_BH,[.02,.9,.06,.8])
        self.assertEqual(result.significant_synergy.tolist(),[True,False,False,False])
        self.assertEqual(result.fdr_family_size.tolist(),[2]*4)
        np.testing.assert_allclose(engine.correct_triplet_fdr(frame,'sample').q_down_BH,[.04,.9,.06,.9])
    def test_seven_signal_sparse_state_local_global_identity(self):
        x=np.random.default_rng(42).integers(0,5,size=(5000,7))
        global_value=engine.information_from_codes(x,5)["hoi"]
        self.assertAlmostEqual(engine.local_oinfo_from_codes(x,5).mean(),global_value,places=12)

    def test_invalid_discrete_states(self):
        with self.assertRaises(ValueError):engine.local_oinfo_from_codes([[0,1,2.1]],5)
        with self.assertRaises(ValueError):engine.information_from_codes([[0,1,5]],5)
    def test_new_run_defaults(self):
        with patch.object(sys, "argv", ["compute_hoi_oinfo.py", "--input-file", "unused", "--output-dir", "unused"]):
            args=engine.parse_args()
        self.assertEqual(args.n_permutations,10000)
        self.assertEqual(args.permutation_level,"all")
        self.assertEqual((args.permutation_min_order,args.permutation_max_order),(3,3))

    def test_cli_rejects_no_permutation_mode(self):
        p=subprocess.run([sys.executable,str(ROOT/'common/02_oinformation/compute_hoi_oinfo.py'),
            '--input-file','unused','--output-dir','unused','--n-permutations','0'],capture_output=True,text=True)
        self.assertNotEqual(p.returncode,0)
        self.assertIn('must be positive',p.stderr)
if __name__=='__main__':unittest.main()
