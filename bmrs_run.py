import os, sys, pickle
import matplotlib.pyplot as plt
import bmrs_utils

from bmrs_data import api_key__

# Options
saveData = 0
saveInt = 0
pltData = 0

# Choose data to pull - uncomment as required. See bau.ids
datatype = 'lolp' # LOLP and derated margins (DRMs)
datatype = 'ixtr' # interconnector flows
datatype = 'dswd' # derived system wide data
datatype = 'sysd' # system demand
datatype = 'flhh' # fuel hh, for wind on/offshore + solar
datatype = 'ftws' # forecast for wind on/offshore + solar
datatype = 'imbl' # day-ahead imbalance forecast

# load in bm_api_utils, bm_data:
bau = bmrs_utils.bm_api_utils(datatype,api_key__) 
bmd = bmrs_utils.bm_data()
if saveData or saveInt:
    from bmrs_data import saveDir

nSteps = bau.get_nSteps()
# nSteps = 2 # <--------- for debugging, uncomment this

for i in range(nSteps):
    print(i+1,'oo',nSteps)
    # get the api response
    r = bau.get_response(i)
    
    # dump the current response for manual inspection, if wanted
    if saveInt:
        fn = os.path.join(saveDir,'xmlDump.xml')
        with open(fn,'w') as file:
            print('XML saved to:\n\t',fn)
            file.write(r.text)
    
    # update the bm_data class
    bmd = bau.process_bm_xml(r,bmd,bau.ad)

bmd.process_bm_data(bau.ad)

if pltData:
    fig,ax = plt.subplots(figsize=(8,4))
    plt.plot_date(bmd.dtmsFull,bmd.dataOut,'.-')
    plt.legend(bmd.headings,loc='upper left')
    plt.xlabel('Date')
    plt.ylabel('Data out')
    plt.tight_layout()
    plt.show()

if saveData:
    fnOut = os.path.join(saveDir,'bmr_'+datatype+'.pkl')
    with open(fnOut,'wb') as file:
        print('\nSaving data to:\n\t',fnOut)
        pickle.dump({'data':bmd.dataOut,
                     'tDict':bmd.tDict,
                     'headings':bmd.headings
                     },file)