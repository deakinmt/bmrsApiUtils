import requests, csv
from datetime import datetime, timedelta
import numpy as np
import xml.etree.ElementTree as ET
from functools import lru_cache
from dateutil import tz

# Choose data to pull - uncomment as required. See bau.ids
opts_list = {
        'lolp':'LOLP and derated margins (DRMs)',
        'ixtr':'interconnector flows',
        'dswd':'derived system wide data',
        'sysd':'system demand',
        'flhh':'fuel hh, for wind on/offshore + solar, based on PSRtype',
        'flxn':"fuel hh, based on Elexon's list (a la ESPENI)",
        'ftws':'forecast for wind on/offshore + solar [day ahead]',
        'ftws_c':'forecast for wind on/offshore + solar [current]',
        'ftws_i':'forecast for wind on/offshore + solar [intraday]',
        'imbl':'day-ahead imbalance forecast',
        'rsys':'rolling system demand',
        'ttrs':'temperatures',
        'xchg':'gbp-euro exchange rates',
        'dsps':'Detailed system prices detsysprices',
        'gcpu':'Installed generating capacity per unit',
        }

# A few miscellaneous funcs...
def tset2stamps(t0,t1,dt):
    """Numpy timeseries starting at t0, ending t1, timestep dt."""
    return np.arange(t0,t1+dt,dt,dtype=object)

def tDict2stamps(tDict):
    """Parameter tDict has keys t0, t1, dt."""
    return tset2stamps(tDict['t0'],tDict['t1'],tDict['dt'])

# date to string functions
d2s = lambda d: d.isoformat()[:13].replace('-','').replace('T','')
s2d = lambda s: datetime(*[int(v) for v in [s[:4],s[4:6],s[6:8],s[8:]]])
m2s = lambda m: '-'.join([f'{ii:02d}' for ii in [m.year,m.month,m.day]])
def s2m(ss):
    # if np.isnan(ss):
    if ss is str:
        return datetime(*[int(v) for v in s.split('-')])
    else:
        return ss


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

def sp2timedelta(sp_):
    hrMnSc = np.array(sp_.split(':')).astype('int')
    nSeconds = int(sum(np.array([3600,60,1])*hrMnSc))
    return timedelta(0,nSeconds)

# Initialise session, and settlement period converter
sp2utcD = loadSp2utcDict(yrA=2000,yrB=2022)
utc2spD = {v:k for k,v in sp2utcD.items()}

requests.session()

class api_d:
    """Data for helping pull out data from the BMRS API.
    
    DOC for API:
    https://www.elexon.co.uk/guidance-note/bmrs-api-data-push-user-guide/
    
    """
    # with respect to the user guide (DOC for API) from 18 December 2019
    ids = {
           'lolp':'5.2.60',
           'ixtr':'5.2.18',
           'dswd':'5.2.51',
           'sysd':'5.2.46',
           'flhh':'5.1.23',
           'flxn':'5.2.17',
           'ftws':'5.1.21',
           'ftws_c':'5.1.21', # 'current' version of ftws 
           'ftws_i':'5.1.21', # 'intraday' version of ftws 
           'imbl':'5.2.38',
           'rsys':'5.2.12',
           'ttrs':'5.2.1',
           'xchg':'5.2.64',
           'dsps':'5.2.52',
           'gcpu':'5.1.19',
           }
    
    def __init__(self,datatype):
        # ---- Forming the URL:
        frmToFormat = {
                        'lolp':['FromSettlementDate','ToSettlementDate'],
                        'ixtr':['FromDate','ToDate'],
                        'dswd':['FromSettlementDate','ToSettlementDate'],
                        'sysd':['FromDate','ToDate'],
                        'flhh':['SettlementDate'],
                        'flxn':['FromDate','ToDate'],
                        'ftws':['SettlementDate'],
                        'ftws_c':['SettlementDate'],
                        'ftws_i':['SettlementDate'],
                        'imbl':['FromDate','ToDate'],
                        'rsys':['FromDateTime','ToDateTime'],
                        'ttrs':['FromDate','ToDate'],
                        'xchg':['SettlementDayFrom','SettlementDayTo'],
                        'dsps':['SettlementDate',],
                        'gcpu':['Year',],
                        }
        
        rpts = {
                'lolp':r'LOLPDRM/v1?',
                'ixtr':r'INTERFUELHH/v1?',
                'dswd':r'DERSYSDATA/v1?',
                'sysd':r'SYSDEM/v1?',
                'flhh':r'B1620/v1?',
                'flxn':r'FUELHH/v1?',
                'ftws':r'B1440/v1?',
                'ftws_c':r'B1440/v1?',
                'ftws_i':r'B1440/v1?',
                'imbl':r'MELIMBALNGC/v1?',
                'rsys':r'ROLSYSDEM/v1?',
                'ttrs':r'TEMP/v1?',
                'xchg':r'EURGBFXDATA/v1?',
                'dsps':r'DETSYSPRICES/v1?',
                'gcpu':r'B1420/v1?',
                }
        
        # Start dates manually found from elexon website
        t0s = {
               'lolp':datetime(2015,11,5), # 23/3/2020
               'ixtr':datetime(2015,2,22), # 23/3/2020
               'dswd':datetime(2014,1,10), # 12/4/2020
               'sysd':datetime(2012,11,8), # 13/4/2020
               'flhh':datetime(2014,12,29), # 13/4/2020
               'flxn':datetime(2015,2,22), # 19/3/2021
               'ftws':datetime(2014,12,30), # 13/4/2020
               'ftws_c':datetime(2018,12,12), # 23/3/2021
               'ftws_i':datetime(2018,12,12), # 23/3/2021
               'imbl':datetime(2014,12,30), # 13/4/2020
               'rsys':datetime(2015,7,14), # 10/6/2020
               'ttrs':datetime(2012,11,10), # 06/7/2020
               'xchg':datetime(2019,12,12), # 06/7/2020
               'dsps':datetime(2014,1,10), # 16/3/2021
               # 'dsps':datetime(2019,12,7), # 16/3/2021
               'gcpu':datetime(2000,1,1), # 18/3/2021
               }
        
        dts = {'lolp':[50],
               'ixtr':[50],
               'dswd':[40],
               'sysd':[30],
               'flhh':[1],
               'flxn':[50],
               'ftws':[1],
               'ftws_c':[1],
               'ftws_i':[1],
               'imbl':[50],
               'rsys':[5],
               'ttrs':[30],
               'xchg':[30],
               'dsps':[0,60*30,],
               'gcpu':None,
               }
        
        # ---- Extracting dates from returned data
        spFormat = {
                    'lolp':'settlementPeriod',
                    'ixtr':'startTimeOfHalfHrPeriod',
                    'dswd':'settlementPeriod',
                    'sysd':'settlementPeriod',
                    'flhh':'settlementPeriod',
                    'flxn':'settlementPeriod',
                    'ftws':'settlementPeriod',
                    'ftws_c':'settlementPeriod',
                    'ftws_i':'settlementPeriod',
                    'imbl':'settlementPeriod',
                    'rsys':'publishingPeriodCommencingTime',
                    'ttrs':None,
                    'xchg':None,
                    'dsps':'settlementPeriod',
                    'gcpu':None,
                    }
        
        spFormatFunc = {
                    'lolp':int,
                    'ixtr':int,
                    'dswd':int,
                    'sysd':int,
                    'flhh':int,
                    'flxn':int,
                    'ftws':int,
                    'ftws_c':int,
                    'ftws_i':int,
                    'imbl':int,
                    'rsys':sp2timedelta,
                    'ttrs':None,
                    'xchg':None,
                    'dsps':int,
                    'gcpu':None,
                    }
        
        sdFormat = {
                    'lolp':'settlementDate',
                    'ixtr':'settlementDate',
                    'dswd':'settlementDate',
                    'sysd':'startTimeOfHalfHrPeriod',
                    'flhh':'settlementDate',
                    'flxn':'startTimeOfHalfHrPeriod',
                    'ftws':'settlementDate',
                    'ftws_c':'settlementDate',
                    'ftws_i':'settlementDate',
                    'imbl':'settlementDate',
                    'rsys':'settDate',
                    'ttrs':'publishingPeriodCommencingTime',
                    'xchg':'settlementDay',
                    'dsps':'settlementDate',
                    'gcpu':'year',
                    }
        
        # ---- Gettung the data from the XML that is returned
        containers = {
                    'lolp':[],
                    'ixtr':[],
                    'dswd':[],
                    'sysd':{},
                    'flhh':{},
                    'flxn':[],
                    'ftws':{},
                    'ftws_c':{},
                    'ftws_i':{},
                    'imbl':[],
                    'rsys':[],
                    'ttrs':[],
                    'xchg':[],
                    'dsps':[],
                    'gcpu':[],
                    }
        
        recordKey = {
                    'lolp':None,
                    'ixtr':None,
                    'dswd':None,
                    'sysd':'recordType',
                    'flhh':'powerSystemResourceType',
                    'flxn':None,
                    'ftws':'powerSystemResourceType',
                    'ftws_c':'powerSystemResourceType',
                    'ftws_i':'powerSystemResourceType',
                    'imbl':None,
                    'rsys':[],
                    'ttrs':None,
                    'xchg':None,
                    'dsps':None,
                    'gcpu':None,
                    }
        
        headsets = {
                    'lolp':['drm12Forecast','lolp12Forecast',
                            'drm8HourForecast','lolp8HourForecast',
                            'drm4HourForecast','lolp4HourForecast',
                            'drm2HourForecast','lolp2HourForecast',
                            'drm1HourForecast','lolp1HourForecast'],
                    'ixtr':['int'+ctry+'Generation' \
                                for ctry in ['fr','irl','ned','ew','nem',
                                                        'elec','ifa2','nsl']],
                    'dswd':['systemSellPrice','systemBuyPrice'],
                    'sysd':['demand'],
                    'flhh':['quantity'],
                    'flxn':['ccgt','oil','coal','nuclear','wind','ps',
                            'npshyd','ocgt','other','intfr','intirl',
                            'intned','intew','biomass','intnem','intelec',
                            'intifa2','intnsl',],
                    'ftws':['quantity'],
                    'ftws_c':['quantity'],
                    'ftws_i':['quantity'],
                    'imbl':['margin','imbalanceValue'],
                    'rsys':['fuelTypeGeneration'],
                    'ttrs':['temperature',
                            'normalReferenceTemperature',
                            'lowReferenceTemperature',
                            'highReferenceTemperature',],
                    'xchg':['settlementExchangeRate'],
                    'dsps':['recordType','id',
                            'soFlag','storFlag','cadlFlag',
                            'offerPrice','bidPrice',
                            'offerVolume','bidVolume',
                            'acceptanceId',
                            ],
                    'gcpu':[
                            'timeSeriesID',
                            'powerSystemResourceType',
                            'registeredResourceEICCode',
                            'bMUnitID',
                            'nGCBMUnitID',
                            'registeredResourceName',
                            'activeFlag',
                            'implementationDate',
                            'decommissioningDate',
                            ],
                    }
        
        dataClass = {
                    'lolp':float,
                    'ixtr':int,
                    'dswd':float,
                    'sysd':float,
                    'flhh':float,
                    'flxn':float,
                    'ftws':float,
                    'ftws_c':float,
                    'ftws_i':float,
                    'imbl':float,
                    'rsys':float,
                    'ttrs':float,
                    'xchg':float,
                    'dsps':lambda x: x,
                    'gcpu':lambda x: x,
                    }
        
        
        self.datatype = datatype
        self.f2f = frmToFormat[datatype]
        self.spf = spFormat[datatype]
        self.spF = spFormatFunc[datatype]
        self.sdf = sdFormat[datatype]
        self.rpt = rpts[datatype]
        self.t0 = t0s[datatype]
        self.hst = headsets[datatype]
        self.dcls = dataClass[datatype]
        self.dt = None if dts[datatype] is None else timedelta(*dts[datatype])
        self.ctr = containers[datatype]
        self.rkey = recordKey[datatype]

class bm_data:
    """A class for concatenating sequential BMRS API data.
    
    """
    sDates = []
    sPrds = []
    data = []
    
    def process_bm_data(self,api_d):
        """Process the bm data ready for saving.
        
        """
        # First, if dsps, then remove NULL elements so int() can be taken
        if api_d.datatype=='dsps':
            data_ = []
            for row in self.data:
                data_.append([])
                nh = len(api_d.hst)
                idsel = api_d.hst.index('acceptanceId')
                for i,v_ in enumerate(row):
                    if i%nh==idsel:
                        v = -1 if v_=='NULL' else v_
                    else:
                        v = v_
                    data_[-1].append(v)
            
            self.data = data_
        
        # Get the timestamps
        if api_d.spF==int or api_d.spF is None:
            self.dtms = [sp2utcD[(self.sDates[i].date(),self.sPrds[i])]
                                        for i in range(len(self.sDates))]
        elif api_d.spF==sp2timedelta:
            self.dtms = [date+dt for date,dt in zip(self.sDates,self.sPrds)]
        
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
        
        
        if api_d.datatype=='dsps':
            self.dataOut = {t:self.l2t_dsps(r,len(self.headings)) 
                                for t,r in zip(self.dtmsFull,data,)}
        elif api_d.datatype=='gcpu':
            self.dataOut = {t.year:self.l2t_gcpu(r,len(self.headings)) 
                                for t,r in zip(self.dtmsFull,data,)}
        else:
            # force all data to the correct point in time
            self.dataOut = np.nan*np.zeros((len(self.dtmsFull),len(data[0])))
            for d,v in zip(self.dtms,data):
                idx = int((d-self.tDict['t0'])/self.tDict['dt'])
                idxSel = np.where(np.isnan(self.dataOut[idx]))[0]
                self.dataOut[idx][idxSel] = [v[i] for i in idxSel]
    
    @staticmethod
    def l2t_dsps(ll,nh):
        """Convert a single list of dsps to a formatted table."""
        r2r = lambda rr: [rr[0][0],rr[1],rr[2]=='T',rr[3]=='T',rr[4]=='T',
            float(rr[5]),float(rr[6]),float(rr[7]),float(rr[8]),int(rr[9])]
        tbl = []
        for ii in range(2*nh,len(ll),nh):
            tbl.append(r2r(ll[ii:ii+nh]))
        
        return tbl
    
    @staticmethod
    def l2t_gcpu(ll,nh,):
        """Convert the gcpu list into a table to output."""
        r2r = lambda rr: rr[:6] + [s2m(rr[6]),s2m(rr[7]),]
        tbl = []
        for ii in range(0,len(ll),nh):
            tbl.append(r2r(ll[ii:ii+nh]))
        
        return tbl

class bm_api_utils():
    """A class for calling BMRS using requests, then processing the XML data.
    
    """
    def __init__(self,datatype,api_key):
        host = 'https://api.bmreports.com'
        port = ':443'
        apikey = 'APIKey='+api_key
        
        self.ad = api_d(datatype)
        self.datatype = datatype
        
        self.url0 = host + port + r'/BMRS/' + self.ad.rpt + apikey
        
        # So far, it seems only the ftws options need an appended query.
        if self.datatype=='ftws_c':
            self.url0 = self.url0 + '&processType=Current'
        elif self.datatype=='ftws_i':
            self.url0 = self.url0 + '&processType=Intraday'
    
    def get_nSteps(self):
        if not self.ad.dt is None:
            return int(np.ceil((datetime.today()-self.ad.t0)/self.ad.dt))
        else:
            return datetime.today().year - self.ad.t0.year
    
    def get_url(self,i):
        fsd, tsd = self.get_time_url(i)
        srv = 'ServiceType=xml'
        url = '&'.join((self.url0,fsd,tsd,srv))
        return url
    
    def get_time_url(self,i):
        if self.ad.dt is None:
            # dt is None when using years instead of dates
            fsd=f'{self.ad.f2f[0]}={self.ad.t0.year+i}'
            return fsd,''
        else:
            dA = self.ad.t0+(i*self.ad.dt)
            
            if self.ad.datatype=='dsps':
                # dsps sets the settlement period INSTEAD of a 'to' date
                dA, sp  = utc2spD[dA]
                tsd = f'SettlementPeriod={sp}'
            
            # Usually a date is ok, for 'FromDateTime' extra info is required.
            if self.ad.f2f[0]=='FromDateTime':
                formatStr = lambda d_: '={:%Y-%m-%d %H:%M:%S}'.format(d_)
                d_dB = timedelta(0,5*60)
            else:
                formatStr = lambda d_: '={:%Y-%m-%d}'.format(d_)
                d_dB = timedelta(1)
            
            fsd = self.ad.f2f[0] + formatStr(dA)
            
            if len(self.ad.f2f)==2:
                dB = (self.ad.t0+((i+1)*self.ad.dt)) - d_dB
                tsd = self.ad.f2f[1] + formatStr(dB)
            elif self.ad.datatype=='dsps':
                pass
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
        
        if api_d.datatype=='gcpu':
            yr = int(root.find('responseMetadata').find(
                                            'queryString').text[-4:])
        
        dspPrev = [0,0]
        for item in rl:
            sdate0 = (item.find(api_d.sdf).text).split('-')
            sdate = datetime(*([int(x) for x in sdate0])) if len(sdate0)==3\
                                        else datetime(yr,1,1,)
            
            # Get the settlement period
            sp = 1 if api_d.spf is None else \
                                    api_d.spF( item.find(api_d.spf).text )
            
            # sometimes there are multiple dates/sps:
            if [sdate,sp]!=dspPrev:
                bmd.sDates.append(sdate)
                bmd.data.append(api_d.ctr.copy())
                bmd.sPrds.append(sp)
                dspPrev = [bmd.sDates[-1],bmd.sPrds[-1]]
            
            for head in api_d.hst:
                if item.find('activeFlag') is None or \
                                        item.find('activeFlag').text=='Y':
                    if api_d.datatype in ['lolp','imbl','dsps','gcpu',]:
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

def data2csv(fn,data,head=None):
    """Write list of lists 'data' to a csv.
    
    If head is passed, put as the first row.
    """
    if not head is None:
        data = [head] + data
    
    with open(fn, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(data)
