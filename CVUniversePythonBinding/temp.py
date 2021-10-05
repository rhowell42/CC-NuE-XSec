l = ["LowQ2Pi","MK","MINOSEfficiency","MnvHadron","RPA","SuSAFromValencia2p2h",]

for i in l:
    print ('#include "PlotUtils/{}Reweighter.h"'.format(i))

for i in l:
    print('template class PlotUtils::{}Reweighter<PythonMinervaUniverse>;'.format(i))

for i in l:
    print(' <class name="PlotUtils::{}Reweighter<PythonMinervaUniverse,PlotUtils::detail::empty>" />'.format(i))
