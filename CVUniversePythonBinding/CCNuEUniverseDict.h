#ifndef CVUNIVERSEPYTHONBINDING_H
#define CVUNIVERSEPYTHONBINDING_H

#include <vector>
#include "PythonMinervaUniverse.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "PlotUtils/MichelSystematics.h"
#include "PlotUtils/AngleSystematics.h"
#include "PlotUtils/ResponseSystematics.h"
#include "PlotUtils/MuonSystematics.h"
#include "PlotUtils/Cut.h"
#include "PlotUtils/CCInclusiveCuts.h"
#include "PlotUtils/Model.h"
#include "PlotUtils/GENIEReweighter.h"
#include "PlotUtils/FSIReweighter.h"
#include "PlotUtils/FluxAndCVReweighter.h"
#include "PlotUtils/LowRecoil2p2hReweighter.h"
#include "PlotUtils/LowQ2PiReweighter.h"
#include "PlotUtils/MKReweighter.h"
#include "PlotUtils/MINOSEfficiencyReweighter.h"
#include "PlotUtils/RPAReweighter.h"
#include "PlotUtils/SuSAFromValencia2p2hReweighter.h"
#include "PlotUtils/TargetMassSystematics.h"
#include "PlotUtils/TargetUtils.h"

template class PlotUtils::Model<PythonMinervaUniverse>;
template class PlotUtils::Reweighter<PythonMinervaUniverse>;
template class PlotUtils::GENIEReweighter<PythonMinervaUniverse>;
template class PlotUtils::FSIReweighter<PythonMinervaUniverse>;
template class PlotUtils::FluxAndCVReweighter<PythonMinervaUniverse>;
template class PlotUtils::LowRecoil2p2hReweighter<PythonMinervaUniverse>;
template class PlotUtils::LowQ2PiReweighter<PythonMinervaUniverse>;
template class PlotUtils::MKReweighter<PythonMinervaUniverse>;
template class PlotUtils::MINOSEfficiencyReweighter<PythonMinervaUniverse>;
template class PlotUtils::RPAReweighter<PythonMinervaUniverse>;
template class PlotUtils::SuSAFromValencia2p2hReweighter<PythonMinervaUniverse>;
template class PlotUtils::BeamAngleXUniverse<PythonMinervaUniverse>;
template class PlotUtils::BeamAngleYUniverse<PythonMinervaUniverse>;
template class PlotUtils::FluxUniverse<PythonMinervaUniverse>;
template class PlotUtils::GenieUniverse<PythonMinervaUniverse>;
template class PlotUtils::GenieRvx1piUniverse<PythonMinervaUniverse>;
template class PlotUtils::GenieFaCCQEUniverse<PythonMinervaUniverse>;
template class PlotUtils::ResponseUniverse<PythonMinervaUniverse>;
template class PlotUtils::MichelEfficiencyUniverse<PythonMinervaUniverse>;
template class PlotUtils::LowQ2PionUniverse<PythonMinervaUniverse>;
template class PlotUtils::RPAUniverse<PythonMinervaUniverse>;
template class PlotUtils::Universe2p2h<PythonMinervaUniverse>;
template class PlotUtils::MuonUniverseMinerva<PythonMinervaUniverse>;
template class PlotUtils::MuonUniverseMinos<PythonMinervaUniverse>;
template class PlotUtils::TargetMassScintillatorUniverse<PythonMinervaUniverse>;
template class PlotUtils::TargetMassWaterUniverse<PythonMinervaUniverse>;
template class PlotUtils::TargetMassCarbonUniverse<PythonMinervaUniverse>;
template class PlotUtils::TargetMassLeadUniverse<PythonMinervaUniverse>;
template class PlotUtils::TargetMassIronUniverse<PythonMinervaUniverse>;
template class PlotUtils::Cut<PythonMinervaUniverse>;
template class reco::HasMINOSMatch<PythonMinervaUniverse>;
template class reco::ZRange<PythonMinervaUniverse>;
template class reco::Apothem<PythonMinervaUniverse>;

#endif /* CVUNIVERSEPYTHONBINDING_H */


