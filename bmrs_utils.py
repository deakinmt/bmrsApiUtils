import requests
from datetime import datetime, timedelta
import numpy as np
import xml.etree.ElementTree as ET
from functools import lru_cache
from dateutil import tz

# A few miscellaneous funcs...
def tset2stamps(t0,t1,dt):
    """Numpy timeseries starting at t0, ending t1, timestep dt."""
    return np.arange(t0,t1+dt,dt,dtype=object)

def tDict2stamps(tDict):
    """Parameter tDict has keys t0, t1, dt."""
    return tset2stamps(tDict['t0'],tDict['t1'],tDict['dt'])

@lru_cache(maxsize=16)
def loadSp2utcDict(yrA=None,yrB=None):
    """Creates a dict for converting (date,sp) tuples to UTC datetimes.
    
    yrA is the starting year, yrB is the end year (exclusive).
    """
    print('Loading sp2utc dict.')
    utc = tz.gettz('UTC')
    toZne = tz.gettz('London')
    
    if yrA is None:
        yrA = 2010
    if yrB is None:
        yrB = 2021
    i=0
    dts0 = [datetime(yrA,1,1)]
    dts = [datetime(yrA,1,1,tzinfo=utc)]
    while dts[-1]<datetime(yrB-1,12,31,23,30,tzinfo=utc):
        dts.append(dts[-1]+timedelta(0,1800))
        dts0.append(dts0[-1]+timedelta(0,1800)) # useful for later

    dtmsUk = [dt.astimezone(toZne) for dt in dts]
    dtsUk = [dtm.date() for dtm in dtmsUk]
    sps = 1
    spsUk = [1]
    dtm_1 = dtsUk[0].day
    for dtm in dtsUk[1:]:
        if dtm_1==dtm.day:
            sps+=1
        else:
            sps=1
        spsUk.append(sps)
        dtm_1 = dtm.day

    sp2utcDict = {}
    for dt,sps,dtUk in zip(dts0,spsUk,dtsUk):
        sp2utcDict[(dtUk,sps)] = dt
    return sp2utcDict



# Initialise session, and settlement period converter
sp2utcD = loadSp2utcDict()
requests.session()

class api_d:
    """Data for helping pull out data from the BMRS API.
    
    DOC for API:
    https://www.elexon.co.uk/guidance-note/bmrs-api-data-push-user-guide/
    
    """
    # with respect to the user guide (DOC for API) from 18 December 2019
    ids = {'lolp':'5.2.60',
           'ixtr':'5.2.18',
           'dswd':'5.2.51',
           'sysd':'5.2.46',
           'flhh':'5.1.23',
           'ftws':'5.1.21',
           'imbl':'5.2.38',
           }
    
    def __init__(self,datatype):
        # ---- Forming the URL:
        frmToFormat = {'lolp':['FromSettlementDate','ToSettlementDate'],
                        'ixtr':['FromDate','ToDate'],
                        'dswd':['FromSettlementDate','ToSettlementDate'],
                        'sysd':['FromDate','ToDate'],
                        'flhh':['SettlementDate'],
                        'ftws':['SettlementDate'],
                        'imbl':['FromDate','ToDate'],
                        }
        
        rpts = {'lolp':r'LOLPDRM/v1?',
                'ixtr':r'INTERFUELHH/v1?',
                'dswd':r'DERSYSDATA/v1?',
                'sysd':r'SYSDEM/v1?',
                'flhh':r'B1620/v1?',
                'ftws':r'B1440/v1?',
                'imbl':r'MELIMBALNGC/v1?',
                }
        
        # Start dates manually found from elexon website
        t0s = {'lolp':datetime(2015,11,5), # 23/3/2020
               'ixtr':datetime(2015,2,22), # 23/3/2020
               'dswd':datetime(2014,1,10), # 12/4/2020
               'sysd':datetime(2012,11,8), # 13/4/2020
               'flhh':datetime(2014,12,29), # 13/4/2020
               'ftws':datetime(2014,12,30), # 13/4/2020
               'imbl':datetime(2014,12,30), # 13/4/2020
               }
        
        dts = {'lolp':50,
               'ixtr':50,
               'dswd':40,
               'sysd':30,
               'flhh':1,
               'ftws':1,
               'imbl':50,
               }
        
        # ---- Extracting dates from returned data
        spFormat = {'lolp':'settlementPeriod',
                    'ixtr':'startTimeOfHalfHrPeriod',
                    'dswd':'settlementPeriod',
                    'sysd':'settlementPeriod',
                    'flhh':'settlementPeriod',
                    'ftws':'settlementPeriod',
                    'imbl':'settlementPeriod',
                    }
        
        sdFormat = {'lolp':'settlementDate',
                    'ixtr':'settlementDate',
                    'dswd':'settlementDate',
                    'sysd':'startTimeOfHalfHrPeriod',
                    'flhh':'settlementDate',
                    'ftws':'settlementDate',
                    'imbl':'settlementDate',
                    }
        
        # ---- Gettung the data from the XML that is returned
        containers = {'lolp':[],
                    'ixtr':[],
                    'dswd':[],
                    'sysd':{},
                    'flhh':{},
                    'ftws':{},
                    'imbl':[],
                    }
        
        recordKey = {'lolp':None,
                    'ixtr':None,
                    'dswd':None,
                    'sysd':'recordType',
                    'flhh':'powerSystemResourceType',
                    'ftws':'powerSystemResourceType',
                    'imbl':None,
                    }
        
        headsets = {'lolp':['drm12Forecast','lolp12Forecast',
                            'drm8HourForecast','lolp8HourForecast',
                            'drm4HourForecast','lolp4HourForecast',
                            'drm2HourForecast','lolp2HourForecast',
                            'drm1HourForecast','lolp1HourForecast'],
                    'ixtr':['int'+ctry+'Generation' \
                            for ctry in ['fr','irl','ned','ew','nem']],
                    'dswd':['systemSellPrice','systemBuyPrice'],
                    'sysd':['demand'],
                    'flhh':['quantity'],
                    'ftws':['quantity'],
                    'imbl':['margin','imbalanceValue'],
                    }
        
        dataClass = {'lolp':float,
                    'ixtr':int,
                    'dswd':float,
                    'sysd':float,
                    'flhh':float,
                    'ftws':float,
                    'imbl':float,
                    }
        
        
        self.datatype = datatype
        self.f2f = frmToFormat[datatype]
        self.spf = spFormat[datatype]
        self.sdf = sdFormat[datatype]
        self.rpt = rpts[datatype]
        self.t0 = t0s[datatype]
        self.hst = headsets[datatype]
        self.dcls = dataClass[datatype]
        self.dt = timedelta(dts[datatype])
        self.ctr = containers[datatype]
        self.rkey = recordKey[datatype]

class bm_data:
    """A class for concatenating sequential BMRS API data.
    
    """
    sDates = []
    sPrds = []
    data = []
    def process_bm_data(self,api_d):
        self.dtms = [sp2utcD[(self.sDates[i].date(),self.sPrds[i])]
                                        for i in range(len(self.sDates))]
        self.tDict = {'t0':min(self.dtms),
                      't1':max(self.dtms),
                      'dt':abs(self.dtms[1]-self.dtms[0])
                      }
        self.dtmsFull = tDict2stamps(self.tDict)
        
        if type(self.data[0]) is dict:
            # get unique headings
            allKeys = [list(row.keys()) for row in self.data]
            self.headings = list(set(
                                [key for keys in allKeys for key in keys]))
            
            # effectively vvvvvv, with error checking.
            # data = [[row[h] for h in headings] for row in self.data]
            data = []
            for row in self.data:
                data.append([])
                for h in self.headings:
                    try:
                        data[-1].append(row[h])
                    except KeyError:
                        data[-1].append(np.nan)
        else:
            data = self.data
            self.headings = api_d.hst
        
        # force all data to the correct point in time
        self.dataOut = np.nan*np.zeros((len(self.dtmsFull),len(data[0])))
        for d,v in zip(self.dtms,data):
            idx = int((d-self.tDict['t0'])/self.tDict['dt'])
            idxSel = np.where(np.isnan(self.dataOut[idx]))[0]
            self.dataOut[idx][idxSel] = [v[i] for i in idxSel]

class bm_api_utils():
    """A class for calling BMRS using requests, then processing the XML data.
    
    """
    def __init__(self,datatype,api_key):
        host = 'https://api.bmreports.com'
        port = ':443'
        self.key = api_key
        apikey = 'APIKey='+self.key
        
        self.ad = api_d(datatype)
        self.url0 = host + port + r'/BMRS/' + self.ad.rpt + apikey
        self.datatype = datatype
    
    def get_nSteps(self):
        return int(np.ceil((datetime.today()-self.ad.t0)/self.ad.dt))
    
    def get_url(self,i):
        fsd, tsd = self.get_time_url(i)
        srv = 'ServiceType=xml'
        url = '&'.join((self.url0,fsd,tsd,srv))
        return url
    
    def get_time_url(self,i):
        dA = self.ad.t0+(i*self.ad.dt)
        fsd = self.ad.f2f[0] + '={:%Y-%m-%d}'.format(dA)
        if len(self.ad.f2f)==2:
            dB = (self.ad.t0+((i+1)*self.ad.dt)) - timedelta(1)
            tsd = self.ad.f2f[1] + '={:%Y-%m-%d}'.format(dB)
        else:
            tsd = 'Period=*'
        return fsd, tsd
    
    def get_response(self,i):
        url = self.get_url(i)
        r = requests.get(url)
        r.raise_for_status
        return r
    
    @staticmethod
    def process_bm_xml(r,bmd,api_d):
        root = ET.XML(r.text)
        rl = root.find('responseBody').find('responseList').findall('item')
        
        dspPrev = [0,0]
        for item in rl:
            sdate0 = (item.find(api_d.sdf).text).split('-')
            sdate = datetime(*([int(x) for x in sdate0]))
            sp = int(item.find(api_d.spf).text)
            
            # sometimes there are multiple dates/sps:
            if [sdate,sp]!=dspPrev:
                bmd.sDates.append(sdate)
                bmd.data.append(api_d.ctr.copy())
                bmd.sPrds.append(sp)
                dspPrev = [bmd.sDates[-1],bmd.sPrds[-1]]
            
            for head in api_d.hst:
                if item.find('activeFlag').text=='Y':
                    if api_d.datatype in ['lolp','imbl']:
                        try:
                            val_ = api_d.dcls(item.find(head).text)
                        except AttributeError:
                            val_ = np.nan
                    else:
                        val = item.find(head)
                        if val.text is None:
                            val_ = np.nan
                        else:
                            val_ = api_d.dcls(val.text)
                else:
                    val_ = np.nan
                
                if type(api_d.ctr) is dict:
                    rcd = item.find(api_d.rkey).text
                    bmd.data[-1].update({rcd:val_})
                else:
                    bmd.data[-1].append(val_)
        
        return bmd