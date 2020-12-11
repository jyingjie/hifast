#!/usr/bin/env python
# coding: utf-8

import os
import re
import sys
from glob import glob

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fname',
                        help='file name; hdf5 file with "mjd" filed or KY file(.xlsx)')
    parser.add_argument('-f', '--force', action='store_true',
                        help='overwriting file if out file exists')
    parser.add_argument('--kyfiles', nargs='*',
                        help='KY files, if not given, guessing from fname')
    parser.add_argument('--outdir',
                       help='output file directory; Default is same with input file if input hdf5 file, "./" if input .xlsx file.')
    parser.add_argument('--t_extra', type=float,
                       help='max allowed extrapolate time; unit: second; defalut is 3s.')
    parser.add_argument('--plot', action='store_true',
                       help='plot the ra dec in pdf image')
    parser.add_argument('-s','--strict', action='store_true',
                       help='(deprecate), cal radec of drift mode strictly from ky')
    parser.add_argument('-l','--loose', action='store_true',
                       help='using the endding of the recorded KY position for Drift mode')

    args = parser.parse_args()

    fname= args.fname

    outdir= args.outdir
    kyfiles= args.kyfiles
    plot= args.plot
    
    outpart = '-radec'
    if outdir is None: outdir = os.path.dirname(fname)
    fileout = os.path.join(outdir, '.'.join(os.path.basename(fname).split('.')[:-1]) + f'{outpart}.hdf5')
    if os.path.exists(fileout):
        if args.force:
            print(f"will overwrite the existing out file {fileout}")
        else:
            print(f"File exists {fileout}")
            print('exit... Using -f to overwrite it.')
            sys.exit()



import json
import numpy as np
import scipy.interpolate as interp
import pandas as pd
import xlrd
import h5py

from astropy import units as u
from astropy.time import Time
from astropy.utils import iers
#iers.Conf.iers_auto_url.set("https://datacenter.iers.org/data/9/finals2000A.all")
from .kypara2radec import kypara2radec
from .util import save_dict_hdf5
    
def guess_ky(fname, obs_mjd, ky_dir=None):
    """
    fname: string
    obs_mjd: scale, UTC
    """
    obs_data= (Time(obs_mjd,format='mjd')+8*u.hour).to_datetime().strftime('%Y_%m_%d')
    
    ky_dir_default=[os.path.expanduser("~")+'/KY/',
                    '/data/inspur_disk06/fast_data/KY/',
                    '/data31/KY/',]
    if ky_dir is None:
        for ky_dir in ky_dir_default:
            if os.path.exists(ky_dir):
                break
    if ky_dir is None:
        raise(ValueError('KY dir not find'))
    keys= re.match('^.*M[0-1][0-9]_[W,N,F]',os.path.basename(fname)).group(0)[:-6].split('_')
    kyck_files=[] #set variable first, in case len(keys)==1
    for i in range(1,len(keys)):
        kyck_files= glob(f"{ky_dir}/**/{'_'.join(keys[:-i])}*{obs_data}*.xlsx", recursive=True)
        if len(kyck_files)>=1:
            break
    if len(kyck_files)==0:
        raise(ValueError("can not find KY file similar to input file in %s,  please use --kyfiles to appoint KY file" % ky_dir))
    return kyck_files
    
def process_ky(kyck_files, mjds_src, tlim, ky_type=None):
    '''
    check KY files, load kc data
    '''
    for kyck_file in kyck_files:
        print(f"check {kyck_file}")
        this=True
        utcoffset = 8*u.hour
        #read ky data
        sheet_name='整控-馈源舱数据'
        book = xlrd.open_workbook(kyck_file)
        if sheet_name not in book.sheet_names():
            sheet_name=0 #if sheet_name is not exist, use first sheet
        kyck_data = pd.read_excel(book,sheet_name=sheet_name)
        if len(kyck_data)==0:
            print('empty sheet')
            continue
        kyck_data=kyck_data.to_dict('series')
        systime = np.array(kyck_data['SysTime'],dtype=str)
        systime = Time(systime) - utcoffset
        mjd = systime.mjd
        if mjds_src is None:
            # if only deal with ky file
            return kyck_data, mjd, kyck_file
        
        diff= mjds_src.min()-mjd.min()
        if -tlim< diff< 0:
            if ky_type!='drift':
                print(f"input mjd is earlier {abs(diff)*24*3600}s than the mjd in KY. Notice: ra and dec will be extrapolated.")
            else:
                print(f"input mjd is earlier {abs(diff)*24*3600}s than the mjd in KY.")
        if -tlim > diff:
            this=False
            print(f"input mjd is earlier {abs(diff)*24*3600}s than the mjd in KY. Abort.")
            #raise(ValueError(f"input mjd exceed {abs(diff)*24*3600}s to the mjd in KY,"))
        if ky_type=='drift':
            diff= mjd.max()-mjds_src.min()
            if -tlim< diff< 0:
                print(f"the begin of input mjd is later {abs(diff)*24*3600}s than the end of mjd in KY. ")
            if -tlim > diff:
                this=False
                print(f"the begin of input mjd is later {abs(diff)*24*3600}s than the end of mjd in KY. Abort.")
                #raise(ValueError(f"input mjd exceed {abs(diff)*24*3600}s to the mjd in KY,"))
        
        diff= mjd.max() - mjds_src.max()
        if -tlim< diff< 0:
            if ky_type!='drift':
                print(f"input mjd is later {abs(diff)*24*3600}s than the mjd in KY, Notice: ra and dec will be extrapolated")
            else:
                print(f"input mjd is later {abs(diff)*24*3600}s than the mjd in KY.")
        if -tlim > diff:
            #raise(ValueError(f"input mjd exceed {abs(diff)*24*3600}s to the mjd in KY,"))
            if ky_type!='drift':
                this=False
                print(f"input mjd is later {abs(diff)*24*3600}s than the mjd in KY. Abort.")

        if this:
            break

    if not this:
        print(f"input mjd are not in the mjd range of those KY files, try to use --kyfiles to appoint KY file")
        sys.exit()
    else:
        return kyck_data, mjd, kyck_file

def kydata2radec(kyck_data, mjd, nBs='All',ky_type=None):
    """
    calculating ra dec of ky from ky data
    """
    if nBs== 'All':
        nBs= range(1,20)
   
    multibeamAngle= kyck_data['SDP_AngleM']

    #实测中心波束相对中心的全局坐标
    globalCenterX = kyck_data['SDP_PhaPos_X']
    globalCenterY = kyck_data['SDP_PhaPos_Y']
    globalCenterZ = kyck_data['SDP_PhaPos_Z']
    #实测下平台的全局姿态角
    globalYaw = kyck_data['SDP_SwtDPose_Y']
    globalPitch = kyck_data['SDP_SwtDPose_P']
    globalRoll = kyck_data['SDP_SwtDPose_R']
    
    if ky_type=='drift':
        std= np.std(multibeamAngle[-10:])/np.pi*180
        if std > 3:
            raise(ValueError(f'multibeamAngle variation at last 10 points is {std} degree.'))
        full_fun= lambda x:np.full(len(mjd), np.mean(x[-5:]))
        multibeamAngle= full_fun(multibeamAngle)
        globalCenterX = full_fun(globalCenterX)
        globalCenterY = full_fun(globalCenterY)
        globalCenterZ = full_fun(globalCenterZ)
        globalYaw =     full_fun(globalYaw)
        globalPitch =   full_fun(globalPitch)
        globalRoll =    full_fun(globalRoll)
        
    kypara2radec_ufun= np.frompyfunc(kypara2radec,9,1)
    
    radec= {}
    radec['mjd']= mjd
    radec['angle']= multibeamAngle
    for nB in nBs:
        res= kypara2radec_ufun(mjd, multibeamAngle, nB, globalCenterX,  globalCenterY, globalCenterZ, globalYaw, globalPitch, globalRoll)
        res= np.vstack(res)
        res= res*180/np.pi # convert radians to degree
        radec['ra'+str(nB)]= res[:,0]
        radec['dec'+str(nB)]= res[:,1]
    return radec

def radec_src(radec_ky, mjds_src):
    """
    get ra dec of source through interpolating that of ky
    """
    mjds_ky= radec_ky['mjd']
    radec_src={}
    radec_src['mjd']= mjds_src
    print('interpolating...')
    radec_src['is_extrapo']= (mjds_src < np.nanmin(mjds_ky)) | (mjds_src > np.nanmax(mjds_ky))
    for key in radec_ky:
        if 'ra' in key or 'dec' in key or 'angle' in key:
            value = radec_ky[key]
            if 'ra' in key:
                value = _tight_ra(value)
            radec_src[key]= interp.interp1d(mjds_ky, value, kind='linear', fill_value ='extrapolate')(mjds_src)
            if 'ra' in key:
                radec_src[key][radec_src[key]<0] += 360 # PyAstronomy.pyasl don't support negative ra
    return radec_src

def main(mjds_src, tlim, fname_g=None, kyck_files=None,outdir=None, ky_type=None):
    # cache_dir
    cache_dir= os.path.expanduser('~/.cache/fast_python/kyradec/')
    os.makedirs(cache_dir, exist_ok=True)
    
    kyck_file=None
    if mjds_src is not None:
        cache_file= cache_dir + 'known_kyfile.json'
        obs_data= (Time(mjds_src[0],format='mjd')+8*u.hour).to_datetime().strftime('%Y_%m_%d')
        cache_key= os.path.basename(fname_g).split('-')[0] + '-' + obs_data
        #if cached file is not exits, init
        if not os.path.exists(cache_file):
                with open(cache_file, "w") as f:
                    json.dump({}, f)
        #try to use cache when kyck_files is not assigned
        if kyck_files is None:
            try:
                with open(cache_file,'r') as f:
                    kyck_file= json.load(f)[cache_key]
                print('get kyfile name from ', cache_file)
                print('using ', kyck_file)
                kyck_data=None 
            except KeyError:
                pass
    #if can't get kyck_file
    if kyck_file is None:
        if kyck_files is None and mjds_src is not None:
            obs_mjd= mjds_src[0]
            kyck_files= guess_ky(fname_g,obs_mjd)
        # test kyfile and load data
        kyck_data, mjds_ky, kyck_file= process_ky(kyck_files, mjds_src, tlim, ky_type=ky_type)
        #cache kyck_file
        if mjds_src is not None:
            with open(cache_file, "r+") as f:
                cache = json.load(f)
                cache.update({cache_key:kyck_file})
                f.seek(0)
                json.dump(cache, f, indent=4)

    # obtain radec of ky
    cache_fname= os.path.join(cache_dir,'.'.join(os.path.basename(kyck_file).split('.')[:-1])+'_cache.npy')
    if os.path.exists(cache_fname) and ky_type!='drift':
        print('using cached ra dec of KY: ',cache_fname)
        radec_ky= np.load(cache_fname,allow_pickle=True)[()]
    else:
        if kyck_data is None:
            kyck_data, mjds_ky, kyck_file= process_ky([kyck_file,], mjds_src, tlim, ky_type=ky_type)
        if ky_type=='drift':
            return kydata2radec(kyck_data, mjds_src, ky_type=ky_type)
        else:
            radec_ky= kydata2radec(kyck_data, mjds_ky, ky_type=ky_type)
            np.save(cache_fname,radec_ky)
    if mjds_src is None:
        return radec_ky
    else:
        return radec_src(radec_ky, mjds_src)

def _tight_ra(ra):
    ra_s= np.sort(ra)
    diff_s= np.diff(ra_s)
    ind_max= np.argmax(diff_s)
    if diff_s[ind_max] < ra_s[0]+360-ra_s[-1]:
        return ra
    else:
        ra= np.copy(ra)
        is_c= ra>=ra_s[ind_max+1]
        ra[is_c]= ra[is_c]-360
        return ra
    
def plot_radec(radec,outname=None,ax=None):
    if ax is None:
        from matplotlib import pyplot as plt
        plt.switch_backend('agg')
        fig,ax = plt.subplots(1,1,figsize=(8,8))
    ind_sort= np.argsort(radec['mjd'])
    for i in range(1,20):
        ax.plot(_tight_ra(radec[f'ra{i}'][ind_sort]),radec[f'dec{i}'][ind_sort])
    
    ax.set_xlabel('ra')
    ax.set_ylabel('dec')
    fig.tight_layout()
    if outname is not None:
        ax.set_title(outname)
        fig.tight_layout()
        fig.savefig(outname)
if __name__ == '__main__':
    
    if args.t_extra is None:
        tlim= 3/3600/24
    else:
        tlim= args.t_extra/3600/24
    
    # load the input mjds
    ky_type=None
    add_mjd=False
    if 'xls' in fname.split('.')[-1]:
        mjds_src= None
        kyfiles= [fname]
        if outdir is None:
            outdir= './'
    else:
        if 'npy' in fname.split('.')[-1]:
            mjds_src= np.load(fname,allow_pickle=True)[()]['mjd']
        elif 'hdf5' in fname.split('.')[-1]:
            with h5py.File(fname,'r') as f:
                mjds_src= f['mjd'][:]
        if 'arcdrift' in os.path.basename(fname) and args.loose:
            ky_type='drift'
        if outdir is None:
            outdir= os.path.dirname(fname)
        
        sep_match= re.findall(r'specs_T_[0-9]{4}_[0-9]{4}',os.path.basename(fname))
        if len(sep_match)>0:
            if 'specs_T_0001' not in sep_match[0]:
                add_mjd= True
                try:
                    fname1st=glob(re.sub(r'specs_T_[0-9]{4}_[0-9]{4}','specs_T_0001_[0-9][0-9][0-9][0-9]',fname))[0]
                except:
                    raise(OSError("can not find first spec file"))
                with h5py.File(fname1st,'r') as f:
                    mjd_0= f['mjd'][:][0]
                mjds_src= np.insert(mjds_src,0,mjd_0)
    radec= main(mjds_src, tlim, fname_g=fname, kyck_files= kyfiles, outdir=outdir, ky_type=ky_type)
    if add_mjd:
        len_mjd= len(radec['mjd'])
        for key in radec.keys():
            try:
                length= len(radec[key])
            except:
                continue
            if length==len_mjd:
                radec[key]=radec[key][1:]
    #saving
    ##record history
    from .util import rec_his
    header=rec_his(args=args)
    
    print('Saving...')
    save_dict_hdf5(fileout, radec, header=header)
    print(f"Saved to {fileout}")
    if plot:
        plot_radec(radec, fileout + '.pdf')