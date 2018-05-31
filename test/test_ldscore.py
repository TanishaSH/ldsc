import ldscore.ldscore as ld
import unittest
import bitarray as ba
import numpy as np
import os
import nose
import ldscore.parse as ps


##########################################################################
#                                  MISC FUNCTIONS                                        #
##########################################################################


def test_getBlockLefts():
    coords = (np.arange(1, 6), np.arange(1, 6), (1, 4, 6, 7, 7, 8))
    max_dist = (5, 0, 2)
    correct = (np.zeros(5), np.arange(0, 5), (0, 1, 1, 2, 2, 2))
    assert np.all([np.all(ld.getBlockLefts(coor, max_d) == ca) for coor, max_d, ca in zip(coords, max_dist, correct)])


def test_block_left_to_right():
    block_left = ((0, 0, 0, 0, 0),(0, 1, 2, 3, 4, 5), (0, 0, 2, 2))
    correct_answer = ((5, 5, 5, 5, 5), (1, 2, 3, 4, 5, 6),(2, 2, 4, 4) )
    assert np.all(([np.all(ld.block_left_to_right(bl) == ca) for bl, ca in zip(block_left, correct_answer)]))


##########################################################################
#                                    BED PARSER                          #
##########################################################################


PROJECT_PATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
PLINK_TEST_FILES_DIR = os.path.join(PROJECT_PATH, 'test/plink_test/')

class test_bed(unittest.TestCase):

    def setUp(self):
        self.M = 8
        self.N = 5
        self.bim = ps.PlinkBIMFile(os.path.join(PLINK_TEST_FILES_DIR, 'plink.bim'))

    def test_bed(self):
        bed = ld.PlinkBEDFile(os.path.join(PLINK_TEST_FILES_DIR, 'plink.bed'), self.N, self.bim)
        # remove three monomorphic SNPs
        print(bed.geno)
        print(bed.m)
        assert bed.m == 4
        # no individuals removed
        print(bed.n)
        assert self.N == bed.n
        # 5 indivs * 4 polymorphic SNPs
        print(len(bed.geno))
        assert len(bed.geno) == 64
        print(bed.freq)
        correct = np.array(
            [0.59999999999999998, 0.59999999999999998, 0.625, 0.625])
        assert np.all(bed.freq == correct)

    def test_filter_snps(self):
        keep_snps = [1, 4]
        bed = ld.PlinkBEDFile(os.path.join(PLINK_TEST_FILES_DIR,'plink.bed'), self.N, self.bim,
                              keep_snps=keep_snps)
        assert bed.m == 1
        assert bed.n == 5
    # pad bits are initialized with random memory --> can't test them
        assert bed.geno[0:10] == ba.bitarray('0001011111')

    def test_filter_indivs(self):
        keep_indivs = [0, 1]
        bed = ld.PlinkBEDFile(os.path.join(PLINK_TEST_FILES_DIR, 'plink.bed'), self.N, self.bim,
                              keep_indivs=keep_indivs)
        assert bed.m == 2
        assert bed.n == 2
        # pad bits are initialized with random memory --> can't test them
        assert bed.geno[0:4] == ba.bitarray('0001')
        assert bed.geno[8:12] == ba.bitarray('0001')

    def test_filter_indivs_and_snps(self):
        keep_indivs = [0, 1]
        keep_snps = [1, 5]
        bed = ld.PlinkBEDFile(os.path.join(PLINK_TEST_FILES_DIR, 'plink.bed'), self.N, self.bim,
                              keep_snps=keep_snps, keep_indivs=keep_indivs)
        self.assertEqual(bed.m, 1)
        self.assertEqual(bed.n, 2)
        print(bed.geno)
        self.assertTrue(bed.geno[0:4] == ba.bitarray('0001'))

    @nose.tools.raises(ValueError)
    def test_bad_filename(self):
        bed = ld.PlinkBEDFile(os.path.join(PLINK_TEST_FILES_DIR, 'plink.bim'), 9, self.bim)

    @nose.tools.raises(ValueError)
    def test_nextSNPs_errors1(self):
        bed = ld.PlinkBEDFile(os.path.join(PLINK_TEST_FILES_DIR,'plink.bed'), self.N, self.bim)
        bed.nextSNPs(0)

    @nose.tools.raises(ValueError)
    def test_nextSNPs_errors2(self):
        bed = ld.PlinkBEDFile(os.path.join(PLINK_TEST_FILES_DIR, 'plink.bed'), self.N, self.bim)
        bed.nextSNPs(5)

    # @param.expand([([1]), ([2]), ([3])])
    def test_nextSNPs(self):
        b_num = [3, 2, 1]
        for b in b_num:
            bed = ld.PlinkBEDFile(os.path.join(PLINK_TEST_FILES_DIR, 'plink.bed'), self.N, self.bim)
            x = bed.nextSNPs(b)
            print("x", x)
            print(np.std(x, axis=0))
            self.assertEqual(x.shape, (5, b))
            self.assertTrue(np.all(np.abs(np.mean(x, axis=0)) < 0.01))
            self.assertTrue(np.all(np.abs(np.std(x, axis=0) - 1) < 0.01))

    def test_nextSNPs_maf_ref(self):
        b = 4
        bed = ld.PlinkBEDFile(os.path.join(PLINK_TEST_FILES_DIR, 'plink.bed'), self.N, self.bim)
        x = bed.nextSNPs(b)
        bed._currentSNP -= b
        y = bed.nextSNPs(b, minorRef=True)
        assert np.all(x == -y)
