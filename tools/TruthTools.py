import ROOT
import math
from operator import itemgetter
import random

def GetIndexlistFromPDG(event,PDG,include_anti_particle):
    cond = (lambda x,y : abs(x)==abs(y)) if include_anti_particle else (lambda x,y :x==y)
    return [i for i,x in enumerate(event.mc_FSPartPDG) if cond(x,PDG)]

def MostEnergeticParticle(event,PDG,anti_particle=True):
    particle_index = GetIndexlistFromPDG(event,PDG,anti_particle)
    if len(particle_index) == 0 :
        return None
    elif len(particle_index) > 1:
        m = max(event.mc_FSPartE[i] for i in particle_index) 
        return event.mc_FSPartE.index(m)
    else : 
        return particle_index[0]

def LeadingParticleEnergy(event,PDG,anti_particle=True):
    particle_index = GetIndexlistFromPDG(event,PDG,anti_particle)
    if len(particle_index) == 0 :
        return None
    else:
        return max(event.mc_FSPartE[i] for i in particle_index)

def CalcApothem(x,y):
    x=abs(x)
    y=abs(y)
    if ( x == 0 or y/x > 1/math.sqrt(3)):
        return (y+x/math.sqrt(3))/2*math.sqrt(3)
    else:
        return x

def EnergyofParticle(event,PDG,anti_particle=True):
    if anti_particle:
        particle_index = [i for i,x in enumerate(event.mc_FSPartPDG) if abs(x) == abs(PDG)]
    else:
        particle_index = [i for i,x in enumerate(event.mc_FSPartPDG) if x == PDG]
    if len(particle_index) == 0 :
        return None
    elif len(particle_index) > 1:
        en = max(event.mc_FSPartE[i] for i in particle_index)
        m = event.mc_FSPartE.index(en) 
        #print event.mc_FSPartPDG[m], event.mc_FSPartE[m]/1e3
        if event.mc_FSPartPDG[m] == 111:
            return event.mc_FSPartE[m]/1e3
        elif event.mc_FSPartPDG[m] == 211:
            return event.mc_FSPartE[m]/1e3 - 0.939
    else :
        return event.mc_FSPartE[particle_index[0]]

def FindThetaFromIndex(event,index):
    if index is None:
        return None
    p = ROOT.Math.RotationX(math.radians(-3.3))*ROOT.Math.PxPyPzEVector(event.mc_FSPartPx[index],event.mc_FSPartPy[index],event.mc_FSPartPz[index],event.mc_FSPartE[index])
    return math.degrees(p.Theta())

def TransverseMomentum(event, PDG,anti_particle=True): 
    if anti_particle:
        particle_index = [i for i,x in enumerate(event.mc_FSPartPDG) if abs(x) == abs(PDG)]
    else:
        particle_index = [i for i,x in enumerate(event.mc_FSPartPDG) if x == PDG]
   
    if len(particle_index) == 0:
        return None
    elif len(particle_index) == 1: 
        index = particle_index[0]
        p = ROOT.Math.RotationX(math.radians(-3.3))*ROOT.Math.PxPyPzEVector(event.mc_FSPartPx[index],event.mc_FSPartPy[index],event.mc_FSPartPz[index],event.mc_FSPartE[index])
        Pt = math.sqrt((p.Px())**2 + (p.Py())**2)/1e3 #assuming beam is in z for beam coordinates 
        return Pt

def TotalKEProtons(event,PDG,anti_particle=True):
    kine = 0
    if anti_particle:
        particle_index = [i for i,x in enumerate(event.mc_FSPartPDG) if abs(x) == abs(PDG)]
    else:
        particle_index = [i for i,x in enumerate(event.mc_FSPartPDG) if x == PDG] 
    if len(particle_index) == 0:
        return None
    elif len(particle_index) > 1:
        for i in particle_index:
            kine = kine + (event.mc_FSPartE[i]/1e3)# - 0.938)
        return kine
    else:
        return (event.mc_FSPartE[0]/1e3)# - 0.938)

def TotalKENeutrons(event,PDG,anti_particle=True):
    kine = 0
    if anti_particle:
        particle_index = [i for i,x in enumerate(event.mc_FSPartPDG) if abs(x) == abs(PDG)]
    else:
        particle_index = [i for i,x in enumerate(event.mc_FSPartPDG) if x == PDG]
    if len(particle_index) == 0:
        return None
    elif len(particle_index) > 1:
        for i in particle_index:
            kine = kine + (event.mc_FSPartE[i]/1e3 - 0.939)
        return kine
    else:
        return (event.mc_FSPartE[0]/1e3 - 0.939)


def FindSubrun(event):
    for i in range(len(event.mc_FSPartPDG)):
        print(event.mc_FSPartPDG[i], event.mc_FSPartE[i], event.prong_dEdXMeanFrontTracker[0])
    

def VertexDifferenceX(event):
    truth_vertexX = event.mc_vtx[0]
    reco_vertexX = event.NuECCQE_vtx[0]
    diffX = abs(reco_vertexX) - abs(truth_vertexX) 
    return diffX

def VertexDifferenceY(event):
    truth_vertexY = event.mc_vtx[1]
    reco_vertexY = event.NuECCQE_vtx[1]
    diffY = abs(reco_vertexY) - abs(truth_vertexY)
    return diffY

def VertexDifferenceZ(event):
    truth_vertexZ = event.mc_vtx[2]
    reco_vertexZ = event.NuECCQE_vtx[2]
    diffZ = abs(reco_vertexZ) - abs(truth_vertexZ)  
    return diffZ
 
def ParticlePDG(event, index): #for high energy pi0s  
    PDG = []
    if index != None:
        if event.mc_FSPartE[index]/1e3 > 4.0: 
            for i in range(len(event.mc_FSPartPDG)):
                PDG.append(event.mc_FSPartPDG[i])   
            return event.mc_incomingE/1e3 

def HighEnergyPi0(event, index, metric):
    if index != None:
        if event.mc_FSPartE[index]/1e3 > 4.0: #if greater than 4GeV
            if metric == 3: #if DIS
                return event.mc_w/1e3 #return W
               
def FractionalPi0E(event, index):
    total = 0
    pi0totE = 0
    if index != None:
        if event.mc_FSPartE[index]/1e3 > 4.0:
            for i in range(len(event.mc_FSPartPDG)): 
                if event.mc_FSPartPDG[i] == 111:
                    pi0totE = pi0totE + event.mc_FSPartE[i]
                else:
                    total = total + event.mc_FSPartE[i]
            fraction = pi0totE/1e3 / (total/1e3 +  pi0totE/1e3) 
            return fraction

def LeadingPi0Fraction(event, index):
    total = 0
    leadingpi0 = 0
    if index != None and event.mc_intType == 3 and event.mc_FSPartE[index]/1e3>4.0: #if DIS event
        leadingpi0 = event.mc_FSPartE[index]
        for i in range(len(event.mc_FSPartPDG)): 
            if i == index: 
                continue
            else:
                total = total + event.mc_FSPartE[i]
        fraction = leadingpi0/1e3 / (total/1e3 + leadingpi0/1e3)
        return fraction

def clusterdEdX(event):
    counter = 0
    total = 0
    fraction = 0
    energy = 0
    totalE = 0
    energypdg = 0 
    for i in range(len(event.eProngClusters_cal_energy)):
        pdg = event.GetVecElem("prong_ClusterCalPDG",0,i)  
        if abs(pdg) == 11: 
             counter += 1.0
             energypdg += event.GetVecElem("prong_ClusterCalEnergy",0,i) 
        elif pdg > 1:
             total += 1.0
             energy += event.GetVecElem("prong_ClusterCalEnergy",0,i)
        elif abs(pdg) == 11:
             total += 1.0        
    if (energypdg + energy) > 0 and event.prong_dEdXMeanFrontTracker[0] > 2.4: 
        fractionalE = energypdg / (energypdg + energy) 
        return fractionalE

def clusterdEdXneutron(event):
    counter = 0
    total = 0
    fraction = 0
    energy = 0
    totalE = 0
    energypdg = 0 
    for i in range(len(event.eProngClusters_cal_energy)):
        pdg = event.GetVecElem("prong_ClusterCalPDG",0,i)  
        if abs(pdg) == 211:
             counter += 1.0
             energypdg += event.GetVecElem("prong_ClusterCalEnergy",0,i)   
        elif pdg > 1:
             total += 1.0 
             energy += event.GetVecElem("prong_ClusterCalEnergy",0,i)
        elif abs(pdg) == 11:
             total += 1.0
             energy += event.GetVecElem("prong_ClusterCalEnergy",0,i)
    if (energypdg + energy) > 0 and event.prong_dEdXMeanFrontTracker[0] > 2.4:
        fraction = counter / (total+counter)
        fractionalE = energypdg / (energypdg + energy) 
        return fraction

def dE(event):
    b = []
    sumde = 0
    counter = 0
    totcounter = 0
    fraction = 0
    sumEde = 0
    energyfraction = 0

    firstbin_Cal = event.GetVecElem("prong_ClusterCal_bin",0,0)

    for i in range(0,25):
       if event.prong_dEdXMeanFrontPositionTracker[0] == 25*i and firstbin_Cal != -1:
           b = [i,i+1,i+2,i+3]
       elif event.prong_dEdXMeanFrontPositionTracker[0] == 25*i and firstbin_Cal == -1:
           b = [i-1,i,i+1,i+2]

    for j in range(0,35):
        ClusterCal_bin = event.GetVecElem("prong_ClusterCal_bin",0,j)
        ClusterCal_de = event.GetVecElem("prong_ClusterCal_de",0,j)
        ClusterPDG = event.GetVecElem("prong_ClusterCalPDG",0,j) 

        bin_number = event.GetVecElem("prong_binned_energy_bin_indices",0,j)
        binned_dE = event.GetVecElem("prong_binned_energy_bin_contents",0,j)
    
        if ClusterCal_bin in b: 
            sumde = sumde + ClusterCal_de  
            if ClusterPDG != -999:
                totcounter = totcounter + 1.0
            if abs(ClusterPDG) == 211:
                counter = counter + 1.0  
                sumEde = sumEde + ClusterCal_de
            #print sumEde
            #print ClusterPDG, ClusterCal_de, sumde, event.prong_dEdXMeanFrontTracker[0]*100,  'pdg de sumde'
    if totcounter != 0:
        fraction = counter / totcounter #returns fraction of pdg seen for the de dx
        energyfraction = sumEde / sumde 
    else:
        pass
        #print "return zero here"
   
    return energyfraction
    #print sumde/10, event.prong_dEdXMeanFrontTracker[0]

def vtxdiff(event):
    reco_vtx = event.vtx
    truth_vtx = event.mc_vtx
    
    diff_x = reco_vtx[0] - truth_vtx[0]
    diff_y = reco_vtx[1] - truth_vtx[1]
    diff_z = reco_vtx[2] - truth_vtx[2] 
    return diff_z

def Leakage(event):
    OutsideConeE = event.kin_cal.reco_OutsideCone_E
    trueE = event.kin_cal.true_E_e 
    recoE = event.kin_cal.reco_E_e
    trueTh = event.kin_cal.true_theta_e
 
    fraction = OutsideConeE / trueE 
    difference = OutsideConeE - trueE*0.008
    #fraction = OutsideConeE/recoE
    angle = event.kin_cal.reco_theta_e 
    recoq0 = event.kin_cal.reco_q0
    true_pt = TransverseMomentum(event, 11)
    
    if 3.5 < true_pt <= 5.0:
   
        return difference

def Zcoordinate(event):
    reco_vtx = event.vtx
    z = reco_vtx[2]/10 #cm   
    return z

def VertexR(event):
    reco_vtx = event.vtx
    x = reco_vtx[0]/10
    y = reco_vtx[1]/10
    R2 = x*x + y*y
    return R2/1e3

def Eth2(event):
    Ee = event.kin_cal.reco_E_e
    theta = math.radians(event.kin_cal.reco_theta_e)
    Eth2 = Ee * (theta * theta) 
    return Eth2

def Eth(event):
    Ee = event.kin_cal.reco_E_e
    theta = math.radians(event.kin_cal.reco_theta_e)
    Eth = Ee * (theta)   
    return Eth

def tEth(event):
    Ee = event.kin_cal.true_E_e
    theta = event.kin_cal.true_theta_e_rad
    tEth = Ee * theta
    return tEth

def binEe(event):
    Ee = event.kin_cal.reco_E_e 
    if len(event.prong_MedianPlaneShowerWidth) != 0:
        width = event.prong_MedianPlaneShowerWidth[0] 
        if 4.0 <= Ee < 10.0:   
            return width 
 
def TransverseShower(event):
    Ee = event.kin_cal.reco_E_e
    asymm = {"Numerator": {}, "Denominator": {}}
    for view in ("X","U","V"):
        for frac_part in asymm:
            br_name = "prong_TransverseShowerAsymmetry%s%s" % (frac_part, view)
            if hasattr(event, br_name):
                asymm[frac_part][view] = getattr(event, br_name)[0] 
    if len(asymm["Numerator"]) == len(asymm["Denominator"]) and len(asymm["Numerator"]) > 0:
        # this is error-weighted mean for X, (U+V).  see DocDB 11319, p. 9
        asymm_xhat = 4 * asymm["Numerator"]["X"] + asymm["Numerator"]["U"] + asymm["Numerator"]["V"]
        asymm_xhat /= 5 * sum(asymm["Denominator"].values())          
   
        asymm_yhat = asymm["Numerator"]["U"] - asymm["Numerator"]["V"]
        asymm_yhat /= math.sqrt(3) * (asymm["Denominator"]["U"] + asymm["Denominator"]["V"])  # uhat - vhat = sqrt(3) * yhat
        
        A =  math.sqrt(asymm_xhat**2 + asymm_yhat**2) 
        return A

def ExcessSort(event):
    recoq0 = event.kin_cal.reco_q0 
    recoq3 = event.kin_cal.reco_q3 

    if 0.6 < recoq3 < 1.0 and recoq0 < 0.25 and event.prong_dEdXMeanFrontTracker[0] > 2.4:
        print(event.mc_run,  event.mc_subrun, event.mc_nthEvtInFile+1, event.mc_intType, event.mc_incomingE, event.mc_incoming,event.prong_dEdXMeanFrontTracker[0]) 
        #return event.kin_cal.reco_E_e
    return event.mc_incoming

def Trans(event):
    cluster = []
    if len(event.prong_axis_vector) != 0:
            theta = event.GetVecElem("prong_axis_vector",0,0)
            phi = event.GetVecElem("prong_axis_vector",0,1)
            
            sx = math.sin(theta)*math.cos(phi)
            sy = math.sin(theta)*math.sin(phi)
            sz = math.cos(theta)
            v2 = ROOT.TVector3(sx, sy, sz)
            v2.RotateX(math.radians(3.3))

            recoX = event.NuECCQE_vtx[0]/1e3
            recoY = event.NuECCQE_vtx[1]/1e3
            recoZ = event.NuECCQE_vtx[2]/1e3

            if recoX == 0 or recoY == 0 or recoZ ==0:
                recoX = event.mc_vtx[0]
                recoY = event.mc_vtx[1]
                recoZ = event.mc_vtx[2]
            vVertex = ROOT.TVector3(recoX, recoY, recoZ)
            for i in range(0,event.GetVecDouble("ExtraEnergyClusters_X").size()):
                ExtraX = event.ExtraEnergyClusters_X[i]/1e3
                ExtraY = event.ExtraEnergyClusters_Y[i]/1e3
                ExtraZ = event.ExtraEnergyClusters_Z[i]/1e3

                vExtra = ROOT.TVector3(ExtraX, ExtraY, ExtraZ)
                v1 = vExtra - vVertex
                #v1 = ROOT.TVector3(abs(ExtraX-recoX),abs(ExtraY-recoY),abs(recoZ-ExtraZ)) 
                longitudinal = v1.Dot(v2)
                transverse = math.sqrt(v1.Dot(v1) - longitudinal*longitudinal)

                energy = event.ExtraEnergyClusters_energy[i]
                totenergy = sum(event.ExtraEnergyClusters_energy)
                transWeighted = (transverse*energy)/totenergy
                if ExtraX == 0 or ExtraY == 0 or ExtraZ == 0:
                    cluster.append(-1000)
                else:
                    cluster.append(transverse)  
            return cluster

def Long(event):

    cluster = []
    if len(event.prong_axis_vector) != 0:
            theta = event.GetVecElem("prong_axis_vector",0,0)
            phi = event.GetVecElem("prong_axis_vector",0,1)

            sx = math.sin(theta)*math.cos(phi)
            sy = math.sin(theta)*math.sin(phi)
            sz = math.cos(theta)
            v2 = ROOT.TVector3(sx, sy, sz)
            v2.RotateX(math.radians(3.3))

            recoX = event.NuECCQE_vtx[0]/1e3
            recoY = event.NuECCQE_vtx[1]/1e3
            recoZ = event.NuECCQE_vtx[2]/1e3

            if recoX == 0 or recoY == 0 or recoZ ==0:
                recoX = event.mc_vtx[0]
                recoY = event.mc_vtx[1]
                recoZ = event.mc_vtx[2]
            vVertex = ROOT.TVector3(recoX, recoY, recoZ)

            for i in range(0,event.GetVecDouble("ExtraEnergyClusters_X").size()):
                ExtraX = event.ExtraEnergyClusters_X[i]/1e3
                ExtraY = event.ExtraEnergyClusters_Y[i]/1e3
                ExtraZ = event.ExtraEnergyClusters_Z[i]/1e3

                vExtra = ROOT.TVector3(ExtraX, ExtraY, ExtraZ)
                v1 = vExtra - vVertex
                #v1 = ROOT.TVector3(abs(ExtraX-recoX),abs(ExtraY-recoY),abs(recoZ-ExtraZ))
                longitudinal = v1.Dot(v2)
                transverse = math.sqrt(v1.Dot(v1) - longitudinal*longitudinal)

                energy = event.ExtraEnergyClusters_energy[i]
                totenergy = sum(event.ExtraEnergyClusters_energy)
                longWeighted = (longitudinal*energy)/totenergy
                if ExtraX == 0 or ExtraY == 0 or ExtraZ == 0:
                    cluster.append(-1000)
                else:
                    cluster.append(longitudinal)
            #print cluster, 'long', len(cluster)
            return cluster
 

def GetWeight(event):
    cluster = []
    if len(event.prong_axis_vector) != 0:
        for i in range(0,event.GetVecDouble("ExtraEnergyClusters_X").size()):
            ExtraX = event.ExtraEnergyClusters_X[i]/1e3
            ExtraY = event.ExtraEnergyClusters_Y[i]/1e3
            ExtraZ = event.ExtraEnergyClusters_Z[i]/1e3
            energy = event.ExtraEnergyClusters_energy[i]
            totenergy = sum(event.ExtraEnergyClusters_energy)
            
            #if ExtraX == 0 or ExtraY == 0 or ExtraZ == 0:
            #    cluster.append(-1000)
            #else:
            cluster.append(energy)
        #print cluster , 'energy', len(cluster)
        return cluster

def InlineVertexStudyE(event): 
    
    if event.prong_dEdXMeanFrontTracker[0] > 2.4 and event.Psi < 0.1 and abs(event.mc_incoming) == 12:
        diffZ = VertexDifferenceZ(event)
        #print event.UpstreamInlineEnergy, diffZ, event.mc_run,  event.mc_subrun, event.mc_nthEvtInFile+1, event.mc_intType, event.mc_incomingE 
        return event.UpstreamInlineEnergy

def InlineVertexStudyVtx(event): 
    if event.prong_dEdXMeanFrontTracker[0] > 2.4 and event.Psi < 0.1 and event.UpstreamInlineEnergy > 10 and abs(event.mc_incoming) == 12:
        diffZ = VertexDifferenceZ(event) 
        return diffZ

def test(event):
    IsCC = event.mc_current == 1
    IsNC = event.mc_current == 2
    IsNotNue = abs(event.mc_incoming) != 12
    Is2p2h = event.mc_intType == 8
    IsElastic = event.mc_intType == 7
    IsCoherent = event.mc_intType == 4
    IsPC = event.mc_processType ==5
    
    IsNuE = abs(event.mc_incoming) == 12
    IsAntiNu = event.mc_incoming < 0
   
    particle = []
    for i in range(0, len(event.mc_FSPartPDG)):
        particle.append(event.mc_FSPartPDG[i])  
  
    IsPi0InFinalState = 111 in particle
    IsHeavyBaryon = (3112 or 3122 or 3212 or 3222 or 4112 or 4122 or 4212 or 4222) in particle
    IsMeson = (211 or -211 or 321 or -321 or 323 or -323 or 111 or 130 or 310 or 311) in particle
    if len(event.mc_FSPartPDG) > 0:
        IsDeexcitationPhoton = event.mc_FSPartPDG[0] == 22 and event.mc_FSPartE[0] < 10
    IsPhoton = 22 in particle and event.mc_FSPartPDG[0] != 22   

    CCQElike = IsCC and IsNuE and not IsPi0InFinalState and not IsMeson and not IsHeavyBaryon and not IsPhoton
    notCCQElike = IsNuE and (IsPi0InFinalState or IsMeson or IsPhoton or IsHeavyBaryon)
    NuEElastic = IsElastic 
    CCPi0 = IsCC and IsPi0InFinalState and IsNotNue 
    NCCohPi0 = IsCoherent and IsNC and IsPi0InFinalState
    
    if CCQElike or notCCQElike or NuEElastic or CCPi0:
        return 1.0
    elif NCCohPi0:   
        return 1.0
    else:
        return 1.0

def Rebin(event): 
    if event.UpstreamInlineEnergy < 0.4:
        return 0.4
    else: 
        return event.UpstreamInlineEnergy

def Print(event):
    if event.prong_dEdXMeanFrontTracker[0] > 2.4 and event.Psi < 0.1 and event.UpstreamInlineEnergy < 10:
        print(event.UpstreamInlineEnergy, 'inline', event.kin_cal.reco_E_e) 
        Long(event)
        VertexDifferenceX(event)
        VertexDifferenceY(event)
        VertexDifferenceZ(event)
        print(event.mc_run,  event.mc_subrun, event.mc_nthEvtInFile+1, event.mc_intType, event.mc_incomingE, event.mc_incoming,event.prong_dEdXMeanFrontTracker[0])
        print(" ") 
    return None

def EvailStudy(event):
    if event.prong_dEdXMeanFrontTracker[0] < 2.4 and abs(event.mc_incoming) == 12: 
        true = event.kin_cal.true_visE
        reco = event.kin_cal.reco_visE
        #if true > 1.0:
            #print true, reco, event.mc_run,  event.mc_subrun, event.mc_nthEvtInFile+1, event.mc_incomingE
    return reco

def OverlappingEnergy(event):
    if event.kin_cal.true_visE > 0.8:
        coneE = event.kin_cal.reco_InsideCone_E
        for i in range(0,len(event.mc_FSPartPDG)):
            pdg = event.mc_FSPartPDG[i]
            energy = event.mc_FSPartE[i]/1e3
            print(pdg, energy, coneE, coneE-energy)
        print(" ")



