import ROOT,PlotUtils
import os
import sys

bins = [
    (9,2),
    (11,3),
    (12,4),
    (13,5)
]

def check(filepath):
    tfile = ROOT.TFile.Open(filepath)
    h = tfile.Get("Eavail_q3_true_signal")
    for b in bins:
        try:
            if h.GetBinContent(b[0],b[1])>0:
                return True
        except AttributeError as e:
            print(e)
            return True
    return False

if __name__=="__main__":
    for i in os.listdir(sys.argv[1]):
        if check("{}/{}".format(sys.argv[1],i)):
            print(i)
