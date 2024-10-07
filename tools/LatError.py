import array
import ctypes
import math
import pprint
import sys

import ROOT

import PlotUtils

sys.path.insert(0, '/exp/minerva/app/users/jyhan/cmtuser/anal_flux/reading_file_inc_spline')
from .KinematicsCalculator import KinematicsCalculator


class LatShiftError(PlotUtils.MnvRecoShifter):

    def __init__(self,  shifts, n_universes=10):        
        super(LatShiftError, self).__init__()
        #self.recoil_shift = 0.0        
        self._rndm_seeds = dict( [(s, 4826354113533 + i*3) for i, s in enumerate(shifts)] )
        self.shifts = shifts
        print("what is shifts number ??? ", self.shifts, self._rndm_seeds)
        self._shift_sigmas = {}
        for shift, rndm_seed in self._rndm_seeds.items():
            #print "shift and rndm_seed = ", shift, rndm_seed
            shifts = getattr(ROOT, "std::vector<double>")()
            #print "shifts = ", shifts, n_universes
            self.GetRandomShiftVector(shifts, rndm_seed, n_universes)
            #print "found shift sigma = ", shifts
            self._shift_sigmas[shift] = shifts            
            #print "shift test = ", shift, self._shift_sigmas[shift]

        self.kin_cal = KinematicsCalculator(correct_beam_angle=True, correct_MC_energy_scale=False, calc_true=True)         

    def recoil_sys(self, intree, name, universe):
        
        recoilError = 0.0
        rename = name + '_Response'

        if rename=='Proton_Response':
            recoilError = intree.part_response_recoil_proton_id_err
        elif rename=='Low_Neutron_Response':
            recoilError = intree.part_response_recoil_low_neutron_id_err
        elif rename=='Mid_Neutron_Response':
            recoilError = intree.part_response_recoil_high_neutron_id_err
        elif rename=='Pion_Response':
            recoilError = intree.part_response_recoil_meson_id_err
        elif rename=='Muon_Response':
            recoilError = intree.part_response_recoil_muon_id_err
        elif rename=='EM_Response':
            recoilError = intree.part_response_recoil_em_id_err
        elif rename=='Other_Response':
            recoilError = intree.part_response_recoil_other_id_err
            
        if recoilError < 0.0:
            recoilError = 0.0

        #recoil_shift = recoilError*self._shift_sigmas[0]
        #recoil_test = self.recoil_shift
        #PlotUtils.MnvRecoShifter.GetShift_ccqe_recoilE( recoil_test, recoilError, universe+1, 1.0)
        #print "This is inside of function => ", PlotUtils.MnvRecoShifter.m_rShiftHadron[universe]
        #print recoil_test, recoilError
        recoil_shift = []
        for i in range(universe):
            recoil_shift.append(recoilError*self._shift_sigmas[name][i])

        return recoil_shift


    def recal_kin_recoilsys(self, Enu, Ee, Thetae, shift):
        
        recal_Enu = Enu + shift/1000.0
        recal_Erecoil = Enu - Ee + shift/1000.0

        recal_q2 = self.kin_cal.Q2_QE(Ee, Thetae, recal_Enu)
        recal_q0 = recal_Erecoil
        recal_q3 = math.sqrt(recal_q2 + math.pow(recal_q0,2))

        recal_kin = {'Enu':recal_Enu, 'q0':recal_q0, 'q3':recal_q3}

        return recal_kin
