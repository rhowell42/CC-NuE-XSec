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
template class PlotUtils::Model<PythonMinervaUniverse>;
template class PlotUtils::Reweighter<PythonMinervaUniverse>;
template class PlotUtils::GENIEReweighter<PythonMinervaUniverse>;
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
template class PlotUtils::Cut<PythonMinervaUniverse>;
template class reco::HasMINOSMatch<PythonMinervaUniverse>;
template class reco::ZRange<PythonMinervaUniverse>;
template class reco::Apothem<PythonMinervaUniverse>;
#endif /* CVUNIVERSEPYTHONBINDING_H */
