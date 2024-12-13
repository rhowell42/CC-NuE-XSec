import os
import numpy as np

def MergeAsimovs():
    while True:
        print("Type the path of the directory that contains the Asimov delta chi2 files")
        path = input()
        if len(path) == 0:
            print("Type a valid path")
            continue

        if os.path.isdir(path):
            results = []
            for f in os.listdir(path):
                f_path = path+'/'+f
                if os.path.isdir(f_path):
                    for f_ in os.listdir(f_path):
                        file = f_path+'/'+f_
                        res = np.loadtxt(file,delimiter=',')
                        if np.isnan(res).any() or np.isinf(res).any():
                            print("Bad delta chi2 computed in {}".format(file))
                            continue
                        results.append(res)
                else:
                    res = np.loadtxt(f_path,delimiter=',')
                    if np.isnan(res).any() or np.isinf(res).any():
                        print("Bad delta chi2 computed in {}".format(f_path))
                        continue
                    results.append(res)
            results = np.array(results).flatten()
            print('saving asimov_deltachi2s.npy with {} entries'.format(results.shape[0]))
            np.save("asimov_deltachi2s",results)
            break
        else:
            print("Not a valid path")
            break

def MergeChi2s():
    while True:
        print("Type the path of the directory that contains the chi2 contours you want to merge")
        path = input()
        if len(path) == 0:
            print("Enter a valid path")
            break

        if os.path.isdir(path):
            data_files = [f for f in os.listdir(path) if "chi2_surface" in f]
            if len(data_files) == 0:
                print("No files found in {}".format(path))
                break
        else:
            print("directory does not exist")
            break

        # ----- sort file names to order PMNS parametesr ----- #
        start = '_m_'
        end = '_Ue4'
        m_names = [float(s[s.find(start)+len(start):s.rfind(end)]) for s in data_files]
        m_names = list(set(m_names))
        m_names.sort()

        start = 'Ue4_'
        end = '.dat.npy'
        e_names = [float(s[s.find(start)+len(start):s.rfind(end)]) for s in data_files]
        e_names = list(set(e_names))
        e_names.sort()

        datas = []
        asimovs = []

        for m in m_names:
            m_data   = []
            m_asimov = []
            for e in e_names:
                data = np.load(path+"chi2_surface_data_m_{}_Ue4_{}.dat.npy".format(m,e))
                asimov = np.load(path+"chi2_surface_pseudodata_m_{}_Ue4_{}.dat.npy".format(m,e))

                m_data.append(data)
                m_asimov.append(asimov)
            datas.append(m_data)
            asimovs.append(m_asimov)

        datas = np.array(datas)
        asimovs = np.array(asimovs)

        np.save("data_chi2s",datas)
        np.save("asimov_chi2s",asimovs)
        print("Done saving files")
        break

if __name__ in "__main__":
    print("Do you want to merge the chi2 contour files? (y/n)")
    ans = input().lower()
    if ans == 'y':
        MergeChi2s()

    print("Do you want to merge asimov delta chi2 files? (y/n)")
    ans = input().lower()
    if ans == 'y':
        MergeAsimovs()
