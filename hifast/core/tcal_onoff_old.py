"""
Created on Sun Jun 16 13:09:57 2019

@author: yixiancao
@author: YingJing
"""
import sys
from astropy.io import fits
import numpy as np
import functools
from astropy import constants as const
from scipy.io.idl import readsav

from .ugdopplerfast import ugdopplerfast

# f21 = 1420.405751
# z = -0.001001 # for m31
# f0 = f21 / (1+z)
#%%
class FastRawData(object):
    """ RAW data from FAST observations.
    
    Attributes
    ----------
    filename: str
        Name of Raw data file (FITS file)
    fileinfo: dict
        Infomation inferred from the filename
    obs: dict
        Basic observation informaiton from observation log (NOW hd0 header). 
    data: Numpy ndarray
        Data
    
    Methods
    ----------
    listInfo():
        Display the basic infomation about the data file.

    """
    
    
    def __init__(self, basename,start=1,stop=None, verbose=False):
        self.filenames = [basename+'%04d.fits'%i for i in range(start,stop+1)]
        self.verbose=verbose
        if self.verbose:
            print('files:')
            for filename in self.filenames:
                print(filename)
        self.hduls = [fits.open(filename) for filename in self.filenames]
#        self.hdul.info() 
        self.hd0s = [hdul[0].header for hdul in self.hduls]
        self.hd1s = [hdul[1].header for hdul in self.hduls]
#        self.hdul.readall()
        self.tdata = [hdul[1].data for hdul in self.hduls] # BinTable DATA
        #[hdul.close() for hdul in self.hduls]

    @property
    def data(self):
        data = np.vstack([tdata['DATA'] for tdata in self.tdata])
        return data
    @property
    def mjds(self):
        mjds = np.hstack([tdata['UTOBS'] for tdata in self.tdata])
        return mjds

#     @property
#     def vrad(self):
#         global f0
#         frad = (f0-self.freq)/f0 * const.c.value/1e3 # km/s
#         return frad
    
    
    @property
    def fileinfo(self):
        fitsname = dict(value = self.filenames[0].split('/')[-1], desc = "First FITS file name")  
        dummy = fitsname['value']
        tinfo = dummy.split('-')[0]
        tinfo = tinfo.split('_') 
        #project = dict(value = 'FAST_M31', desc = "Project name")
        source = dict(value = tinfo[0], desc = "Source name")
        ipos = dict(value  = int(tinfo[1][1:]), desc = "Position number")
        #iobs = dict(value = None, desc = "Observation number")
        track_mode = dict(value = tinfo[2], desc = "Observation mode") # e.g. Tracking, drifting
        
        finfo = dummy.split('-')[1]
        finfo = finfo.split('_')
        ibeam = dict(value = int(finfo[0][1:]), desc = "Beam number") 
        #ifile = dict(value = int(finfo[-1][0:4]), desc = "File number") #ith file for the observation
        
        filetype = dict(value = 'psr' if len(finfo) < 3 else 'spec' +  finfo[1],
                         desc = "Observation type")
            
#         fileinfo = dict(fitsname = fitsname, project=project, source = source, ipos = ipos, track_mode = track_mode,
#                         iobs =iobs, ibeam = ibeam, ifile = ifile, filetype = filetype)
        fileinfo = dict(fitsname = fitsname, source = source, ipos = ipos, track_mode = track_mode,
                         ibeam = ibeam, filetype = filetype)

        return fileinfo
    
    def listInfo(self, list_type = 'all'):
        finfo = self.fileinfo
        print("FAST RAW DATA FILE Summary")
        print("==================")
        for key, value in finfo.items():
            print (value['desc'], ':' , value['value'])

    #@property 
    def obs(self,chunk_num=0):
        """ Dictionary of observation parameters. 
            From HD0 or from observation log
        """
        hd0 = self.hd0s[chunk_num]
        time = {"t_obs": dict(value  = hd0['DATE'], 
                              desc = 'Observation time (starting to record data)') , 
                "t_samp": dict(value = None, 
                               desc = 'Sampling time interval (s)') 
                }        
        track = {"Mode": dict(value = self.fileinfo['track_mode'], 
                                    desc =  "Tracking mode (e.g. tracking, drifting, etc)"),
                "RA": dict(value = np.double(0.0), desc = "Tracking RA"), 
                "DEC": dict(value = np.double(0.0), desc = "Tracking DEC")}
        
        noise_diode = {"Mode": dict(value = "",  desc = "Noise diode mode (Modulate, ON, or OFF)"),
                      "Power": dict(value = "",  desc = "Noise diode power (High/Low)"), 
                      "Delay": dict(value = 0, desc =  "Noise Delay in unit of tunit"), 
                      "On": dict(value = 0,  desc = "Noise diode ON in unit of tunit"), 
                      "Off":  dict(value = 0,  desc = "Noise diode OFF in unit of tunit"), 
                      "tunit":  dict(value = 4e-9, desc =  "Units of delay/on/off time (s)")
                      } #
        receiver = {"Frontend": dict(value = '', desc =  "Receiver frontend (e.g. 19BEAM)"), 
                    "rfgain": dict(value = 0.0, desc =  "Receiver rfgain (dB)"),  
                    "dgain": dict(value = 0.0, desc = "Receiver dgain"), 
                    "Backend": dict(value = '', desc = "Receiver Backend ID (e.g. 'MB4K')") 
                    }     
        
        obs = dict(time = time,
               track = track,
               noise_diode = noise_diode,
               receiver = receiver)
        
        return obs
    
#     def listObs(self, list_type = 'all'):
# #        finfo = self.fileinfo
#         obs = self.obs

#         print("FAST RAW DATA Observation Summary")
#         print("==================")
#         print("File parth: " , self.filename)
        
#         for keygroup in obs.keys():
#             print("------------------")
#             print(keygroup.upper() + ' INFO')
#             for key, value in obs[keygroup].items():
#                 print (value['desc'], ':' , value['value'])
                
#     def get_baseline(self, spec, deg = 1, fwid = 15):
#         """ Get the baseline for a single spectra. 
        
#         """
#         selchan = np.argwhere((np.abs(self.freq-f0) < fwid) & (np.abs(self.freq-f0) > 0.1 * fwid))
#         selchan = selchan.flatten()
#         pfit = np.polyfit(self.freq[selchan], spec[selchan], deg) 
#         pfunc = np.poly1d(pfit) 
#         return pfunc

class tcal_onoff(FastRawData):
    """ FAST raw data from spec backend.  
    """
#     def __init__(self, filename):
#         super().__init__(filename)
#         self.filetype = self.fileinfo['filetype']['value']
#         if self.filetype == 'specN':
#             suffix = "_N"
#         elif self.filetype == 'specW':
#             suffix = "_W"
#         else:
#             suffix = ""
            
#         self.psrfile = self.filename.replace(suffix,'')
#         self.nspecfile = self.filename.replace(suffix,'_N')
#         self.wspecfile = self.filename.replace(suffix,'_W')

    @property
    def freq(self):
        nchan = self.tdata[0]['NCHAN'][0]
        freq0 = self.tdata[0]['FREQ'][0]
        chanwidth = self.tdata[0]['CHAN_BW'][0]
        freq=(np.arange(nchan))*chanwidth
        freq=freq+freq0
        return freq
             
#     @property
#     def obs(self):
#         # For M31 project observed with psr
#         psrfile = self.psrfile
#         obs = FastRawPsr(psrfile).obs
# #        obs['time']['t_obs'] = obs['time']['t_obs'].replace('Z','')
#         obs['time']['t_obs']['value'] = self.hd0['DATE'].replace('Z','')
#         obs["time"]["t_int"] = dict(value = self.tdata['EXPOSURE'][0], 
#                                     desc = 'Integration time (s)')
#         return obs
    
    def get_spectrum(self, nint = 0, pl = 'xx'):
        """" Get one spectrum for inspection. """
        ipl = 0 if pl == 'xx' else 1
        return self.data[nint, :, ipl]    
    
    def sep_onoff(self, n1, n2, freq_range=None, ave_cycle= True):
        """
        separate the noise on and off spec

        Parameters
        ----------
        n1 : int
            n1 times the Sampling Time
        n2 : int 
            n1 times the Sampling Time
        ave_cycle: bool
            if True, average the value of on and off in same cycle, respectively.
        Returns
        -------
        on and off 
            
        """
        if freq_range is not None:
            self.use= (self.freq< freq_range[1] ) & (self.freq> freq_range[0])
            data_t= self.data[:,self.use,:2]
        else:
            data_t= self.data[:,:,:2] ## on use XX,YY
        
        self.on, self.off, self.drop_num= _sep_onoff(data_t, n1, n2, ave_cycle)
        if self.drop_num >0:
            print("There are %d incomplete cycle sample data in the tail. Drop..."% self.drop_num)
            
    def sep_mjds_onoff(self, n1, n2, ave_cycle= True):
        """
        separate the noise on and off spec

        Parameters
        ----------
        n1 : int
            n1 times the Sampling Time
        n2 : int 
            n1 times the Sampling Time
        ave_cycle: bool
            if True, average the value of on and off in same cycle, respectively.
        Returns
        -------
        on and off 
            
        """
        data_t= self.mjds
        
        self.on_mjds, self.off_mjds, _= _sep_onoff(data_t, n1, n2, ave_cycle)
    
    def get_T(self,):
        tcal_dir= '../Tcal/'
        s_type= self.fileinfo['filetype']['value'][-1].lower()
        #if s_type!='n':
        #    s_type='w'
        s_type='w'
        nB= self.fileinfo['ibeam']['value']
        ## for very low z, narrow Tcal is sufficient for narrow spec
        tc_info= readsav(tcal_dir+'median_20190115.Tcal-results.HI_%s.high.sav'%s_type)['high_%s'%s_type][0]
        tc_freq= tc_info['freq']
        tc = np.mean(tc_info['M%02d_TC'%nB], axis=0)

        #ind   = (freqcal >= freqlims[0]) & (freqcal <= freqlims[1])
        tc_spec_freq = np.interp(self.freq,tc_freq,tc)

        on_m = np.mean(self.on[:,:,:2],axis=(0,2))
        off_m = np.mean(self.off[:,:,:2],axis=(0,2))
        T= off_m/((on_m-off_m)/tc_spec_freq)

        self.T=T
        self.tc_spec_freq= tc_spec_freq
        self.on_m= on_m
        self.off_m= off_m
        
     
    def get_vlsr(self, ra, dec):
        '''
        ra, dec: deg; array or scalar
           len(ra) ==1 or len(ra)==len(self.mjds)
        '''
        restfreq = np.array([1420.405751])
        c        =  2.99792458e5  # km/s
        velo     =  c*(restfreq[0]-self.freq)/restfreq[0] 
        ra,dec    =  np.asarray(ra),np.asarray(dec)
        
        if len(ra)==1 and len(dec)==1:
            mjd=self.mjds[10:11]
        else:
            mjd= self.mjds
        jd = mjd+2400000.5

        vlsrcor  = ugdopplerfast(ra,dec,jd)
        #print( vlsrcor,mjd[300])

        self.vlsr  = np.vstack([velo-ii for ii in vlsrcor])
        if len(vlsrcor)==1:
            self.vlsr= self.vlsr[0]
#     def calib_nint(self,nint = 0):
#         """ Scale T_cal to T_sys for pulsar spectrum. """
#         subint = nint * 20
#         nfile = (subint+1)//256 + 1
#         if nfile > 1:
#             self.psrfile  = self.psrfile.replace('0001',str(nfile).rjust(4,'0'))
#             self.psr = FastRawPsr(self.psrfile)
#         on, off = self.psr.get_avg_spectra(nint = nint)
        
#         t_cal = 12
#         tsys = t_cal * (on + off * 1.4)/2.4/(on-off)
                
#         return tsys 

def _sep_onoff(data_t, n1, n2, ave_cycle= True):
    n_on,n_off= n1,n2
    inds= np.arange(len(data_t))
    drop= len(inds)%(n_on+n_off)
    inds= inds[:len(inds)-drop]# drop data of incomplete cycle in the tail
    nth= inds%(n_on+n_off)
    #inds[nth<n_on]
    #inds[nth>=n_on]
    
    on= data_t[np.where(nth<n_on)[0]].astype('float64') #use float64 to avoid inaccurate when using np.mean
    off= data_t[np.where(nth>=n_on)[0]].astype('float64') #
    #average on and off in same cycle
    if ave_cycle:
        on= np.mean(on.reshape(tuple([-1,n_on]+list(data_t.shape[1:]))),axis=1)
        off= np.mean(off.reshape(tuple([-1,n_off]+list(data_t.shape[1:]))),axis=1)
    return on, off, drop
