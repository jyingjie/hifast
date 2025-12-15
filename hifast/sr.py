###############
# Author
# Qingze Chen
###############

__all__ = ['parser', 'IO']

# Cell
from hifast.utils.io import *
import numpy as np
from spectral_cube import SpectralCube
from tqdm import tqdm
import numpy as np
#from scipy.ndimage import rotate
#nbdev_comment _all_ = ['parser']

# Internal Cell
sep_line = '##'+'#'*70+'##'
parser = ArgumentParser(prog=f"python -m hifast.{os.path.basename(sys.argv[0])[:-3]}", formatter_class = formatter_class, allow_abbrev=False,
                        description='Stray-radiation correction', )
add_common_argument(parser)
parser.add_argument('fpath',
                    help='input file path.')
parser.add_argument('--frange', type=float, nargs=2, default=[0, float('inf')],
                    help='Limit frequence range')
parser.add_argument('--no_radec', action='store_true', not_in_write_out_config_file=True,
                    help="if set, don't check or add ra dec")

parser.add_argument('--fitspath', required=True,
                   help='the path of the FITS file.')
# parser.add_argument('--pat', '--pattern_fpath', 
#                    help='path of the file stored pattern of 19-beams.')
parser.add_argument('--ra_range', type=float,nargs=2, default=[0, float('inf')],
                   help='ra range; unit: deg')
parser.add_argument('--dec_range', type=float,nargs=2, default=[0, float('inf')],
                   help='dec_range; unit: deg')




###### load fits cube ############
class in_fits(object):
    def __init__(self,f_file):
        hdu=fits.open(f_file)
        self.hdu=fits.open(f_file)
        i=self.hdu[0].header
        self.ARefP=i['CRPIX1']
        self.DRefP=i['CRPIX2']
        self.VRefP=i['CRPIX3']
        self.UnitA=i['CDELT1']
        self.UnitD=i['CDELT2']
        self.UnitV=i['CDELT3']
        self.StartA=i['CRVAL1']
        self.StartD=i['CRVAL2']
        self.StartV=i['CRVAL3']
        self.numV=i['NAXIS3']
        self.FitsDat=hdu[0].data
        self.CUnit3=i['CUNIT3']
###### load pattern and ita ############
class pattern_ita(object):
    def __init__(self, p_file=os.path.dirname(__file__) + '/core/data/pattern_eta.npz'):
        PI=np.load(p_file)
        self.ita=PI['eta']
        self.pattern=PI['pattern']

# Cell
class IO(BaseIO):
    ver = 'old'
    def _get_fpart(self,):
        return '-srmod'
    
    def _get_rotate_angle(self):
        args = self.args
        # use M01 M08
        m01_path = args.fpath.replace(f'M{self.nB:02d}', 'M01')
        m08_path = args.fpath.replace(f'M{self.nB:02d}', 'M08')

        from .ripple.util import Read_hdf5, Args
        h1 = Read_hdf5(Args(m01_path))
        h8 = Read_hdf5(Args(m08_path))

        from astropy.coordinates import SkyCoord
        from astropy import units as u
        coord1 = SkyCoord(h1.ra[0], h1.dec[0], unit = (u.deg, u.deg))
        coord8 = SkyCoord(h8.ra[0], h8.dec[0], unit = (u.deg, u.deg))
        
        # position_angle is the angle with N longitude
        rot_angle = 90 - coord1.position_angle(coord8).deg
        
        self.angle = None
        for a in [23.4, 53.4, 0. ]:
            if np.abs(rot_angle - a) < 1:
                self.angle = a
        
        if self.angle is None:
            raise ValueError("Only support multibeam rotation angle 23.4, 53.4 and 0.")
        else:
            print("Rotation angle:", self.angle)
         
    def rotate(self,img):
        args = self.args
        
        self._get_rotate_angle()
        
        angle = (90+self.angle)*np.pi/180  #23.4&53.4//53.4 is the true direction. 180 deg is for substraction
        # 读取库图片 Attention 转换 默认为int8 运算时可能会溢出
        # 设置新的图像大小
        h,w = img.shape[0],img.shape[1]
        newW = int(w*abs(np.cos(angle)) + h*abs(np.sin(angle)))+1
        newH = int(w*abs(np.sin(angle)) + h*abs(np.cos(angle)))+1
        # Attention dtype
        newimg1 = np.zeros((newH,newW,3))
        newimg2 =  np.zeros((newH,newW,3))
        newimg3 =  np.zeros((newH,newW,3))
        # 设置旋转矩阵 scr -> dex
        trans1 = np.array([[1,0,0],[0,-1,0],[-0.5*w,0.5*h,1]])
        trans1 = trans1.dot(np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]]))
        trans1 = trans1.dot(np.array([[1,0,0],[0,-1,0],[0.5*newW,0.5*newH,1]]))
        # des -> src
        trans2 = np.array([[1,0,0],[0,-1,0],[-0.5*newW,0.5*newH,1]])
        trans2 = trans2.dot(np.array([[np.cos(angle),np.sin(angle),0],[-np.sin(angle),np.cos(angle),0],[0,0,1]]))
        trans2 = trans2.dot(np.array([[1,0,0],[0,-1,0],[0.5*w,0.5*h,1]]))
        # 开始旋转
        for x in range(w):
            for y in range(h):
                newPos = np.array([x,y,1]).dot(trans1)
                newimg1[int(newPos[1])][int(newPos[0])] = img[y][x]

        for x in range(newW):
            for y in range(newH):
                srcPos = np.array([x,y,1]).dot(trans2)
                if srcPos[0] >= 0 and srcPos[0] < w and srcPos[1] >= 0 and srcPos[1] < h:
                    # 双线性内插
                    bix,biy = int(srcPos[0]),int(srcPos[1]) # 取左上角坐标
                    # 避免最后一行溢出
                    if bix < w-1 and biy < h-1:
                        # 沿 X 方向线性内插
                        rgbX1 = img[biy][bix] + (img[biy][bix+1] - img[biy][bix])*(srcPos[0]-bix)
                        rgbX2 = img[biy+1][bix] + (img[biy+1][bix+1] - img[biy+1][bix])*(srcPos[0]-bix)
                        # 沿 Y  方向内插
                        rgb = rgbX1 + (rgbX2-rgbX1)*(srcPos[1]-biy)
                        newimg3[y][x] = rgb
        return newimg3[:,:,2]
    
    def func(self, x, y):
        if x<=39 and x>=0 and y<=39 and y>=0:
            return self.kernel[int(x)][int(y)]
        else:
            return 0
        
    def rotatefunc(self, x, y):
        if x<=39 and x>=0 and y<=39 and y>=0:
            return self.rotate(self.kernel)[int(x)][int(y)]
        else:
            return 0
    #SPfunc=np.vectorize(rotatefunc)
    def func2(self, x, y):
        #b=get_nB(self.args.fpath)-1
        r1=self.rotate(self.kernel)
        r2=r1[7:47,7:47]
        if x<=39 and x>=0 and y<=39 and y>=0:
            return r2[int(x)][int(y)]
        else:
            return 0
    #SPfunc=np.vectorize(rotatefunc)
    @property
    def SPfunc(self,):
        return np.vectorize(self.func2)

################################
    def vel_loop(self,n_channel,de,al,OpCube,dec,Ta,p,mpix,sk,C_start,ita,r2):
        krange=13
        if mpix==10:
            krange=20
        for c in range(n_channel-1):
            # print('       ',c)
            cor=0
            for k in range(krange):#Dec DIRECTION
                if de[k]<0 or de[k]>39:
                    cor=cor
                else:
                    for l in range(int(krange)):#R.A DIRECTION
                        #print(al[l])
                        if al[l]<=0 or al[l]>39:
                            cor=cor
                        else:
                        #print(SPfunc(de[l],al[k]))   test
                            cor=cor+OpCube[c][k][l]*r2[al[l],de[k]]
                #multiple with 3 because of the pixel size
                    #SSprint(cor)
                    #multiple with 3 because of the pixel size
                    #cor=cor+OpCube[c][k][l]*newfunc(20-3*np.cos(dec[b,1270+p]*np.pi/180)*(13-k)-BiasA, 3*l-1-BiasD) #old version
            if cor<0 or Ta[p,C_start+c]<0:
                #pointscube[0][p][c]=Ta[p,C_start+c]+0.0
                self.pointscube[0][p][c]=(Ta[p,C_start+c]/ita-cor/ita*((1-ita)*mpix*mpix*np.cos(dec[p]*np.pi/180)/sk))
            else:  
                self.pointscube[0][p][c]=(Ta[p,C_start+c]/ita-cor/ita*((1-ita)*mpix*mpix*np.cos(dec[p]*np.pi/180)/sk))#normalize it 


    def cut_spec_num(self):
        args = self.args
        s2p = self.s2p[:]
        Ta0=self.s2p
        # merge polar
        if len(Ta0.shape) == 3:
            Ta0 = np.nanmean(Ta0, axis = -1)
            
        ra_range = args.ra_range
        dec_range = args.dec_range
        freq_range = args.frange

        ra0=self.ra
        dec0=self.dec
        vel=self.vel
        freq=self.freq

        okno = (ra0 > ra_range[0]) & (ra0 < ra_range[1]) & (dec0 > dec_range[0]) & (dec0 < dec_range[1])
        self.ra=ra0[okno]
        self.dec=dec0[okno]
        self.s2p=Ta0[okno, :]
        self.mjd = self.mjd[okno]
        
        return okno
                
#############################
    def gen_s2p_out(self,):
        args = self.args
        
        # cut spec in ra dec range
        okno = self.cut_spec_num()
        
        # load fits
#############################################################        
        fitspath=args.fitspath
        fits_cube=in_fits(fitspath)
        self.PI=pattern_ita()
        ARefP=fits_cube.ARefP
        DRefP=fits_cube.DRefP
        VRefP=fits_cube.VRefP
        UnitA=fits_cube.UnitA
        UnitD=fits_cube.UnitD
        UnitV=fits_cube.UnitV
        StartA=fits_cube.StartA
        StartD=fits_cube.StartD
        StartV=fits_cube.StartV
        numV=fits_cube.numV
        FitsDat=fits_cube.FitsDat
        hdu=fits.open(fitspath)
        SCcube = SpectralCube.read(fitspath).with_spectral_unit(u.km / u.s) 
        fitswcs=WCS(hdu[0].header)
#to function 

########################################
        
        
        #### parameters #########               
        P_start=0# the start of points
        n_point=len(self.ra) # the number of points
        
        #  spectral axis cut
        vel = self.vel
        fitsvel=SCcube.spectral_axis.value
        
        ch0=min(fitsvel[0],vel[0])
        ch1=max(fitsvel[-1],vel[-1])
        self.is_use_freq = is_use_vel = (vel >= ch1) & (vel <= ch0)
        print("Spectral range:", ch1, ch0)
        self.is_use_freq2 = is_use_vel
        velv = vel[is_use_vel]
        
        C_start=(np.arange(len(vel))[is_use_vel])[0]         # the start of channels
        n_channel=len(velv) 
        
########### interrupt when channel of velocity not consistent
        if np.abs(fitsvel[1]-fitsvel[0]-(velv[1]-velv[0]))>0.01:
            print("The channel of velocity not consistent")
            print("fits:"+str(fitsvel[1]-fitsvel[0]))
            print("hdf5:"+str((vel[1]-vel[0])))

            raise(ValueError("The channel of velocity not consistent"))
            
        #########################
        self.pointscube=np.zeros((1,n_point,n_channel))# beams,points,channels
        beam_start = time.perf_counter()
        
        b = self.nB - 1
        ita=self.PI.ita[b]
        self.kernel=kernel=self.PI.pattern[b]
        r1=self.rotate(kernel)
        sumk=np.sum(kernel)
        r2=r1[7:47,7:47]
        x = np.linspace(0,38,39) 
        y = np.linspace(0,38,39) 
        x, y = np.meshgrid(x, y)
     #   arcmin30=IO.SPfunc(x,y)          
    #    newfunc = interpolate.interp2d(x, y, arcmin30, kind='cubic')#newfunc为一个函数 


        if np.abs(UnitD)>0.024 and np.abs(UnitD)<0.026:
            print("pix size is 90 arcsec")
            mDel=7
            mAlf=7
            mpix=3
        elif np.abs(UnitD)>0.016 and np.abs(UnitD)<0.017:
            print("pix size is 60 arcsec") 
            mDel=10
            mAlf=10
            mpix=2
        else:
            raise(ValueError('size of pixels not support'))
        
        ra = self.ra
        dec = self.dec
        Ta = self.s2p
        
        iter_ = tqdm(range(n_point), desc='CPU 0: ', mininterval=3)
        for p in iter_:#p:point
    #        print('    ',p)
            NAlf=(ra[P_start+p]-StartA)/UnitA   #how many pixels from this point to ref point
            NDel=(dec[P_start+p]-StartD)/UnitD
            Alf=int(np.round(fitswcs.world_to_pixel_values(ra[P_start+p],dec[P_start+p], 0)[0]))      
            Del=int(np.round(fitswcs.world_to_pixel_values(ra[P_start+p],dec[P_start+p], 0)[1]))
            BiasA=NAlf+ARefP-Alf
            BiasD=NDel+DRefP-Del
            #print('fits start at',int(round((velv[0]-fitsvel[0])/UnitV)))
            #print('fits ends at',int(round((velv[-1]-fitsvel[-1])/UnitV)))
            Refcube=FitsDat[int(round((velv[0]-fitsvel[0])/UnitV)):int(round((velv[0]-fitsvel[0])/UnitV))+n_channel,Del-mDel:Del+mDel,Alf-mAlf:Alf+mAlf]
            #Refcube=FitsDat[int(round((velv[0]-fitsvel[0])/UnitV)):int(round((velv[0]-fitsvel[0])/UnitV))+n_channel,Del-mDel:Del+mDel,Alf-mAlf:Alf+mAlf] #Cube:[freq,dec,ra]
            OpCube0=copy.deepcopy(Refcube)
            OpCube=np.nan_to_num(OpCube0)
            #%%
            #print(np.shape(OpCube))
            #%%
            if(Del<mDel or Del>np.shape(FitsDat)[1]-mDel or Alf<mAlf or Alf>np.shape(FitsDat)[2]-mAlf):
                OpCube=np.zeros([n_channel,mDel*2,mAlf*2])   
            al=[int(20-mpix*(l-mAlf)) for l in range(mAlf*2+1)]
            de=[int(20-mpix*(k-mDel)) for k in range(mDel*2+1)] 

#            if p%200==0:
#                print ("point="+str(p)+"   time=t+"+str((np.round(time.perf_counter()-beam_start,4))))
            self.vel_loop(n_channel,de,al,OpCube,dec,Ta,p,mpix,sumk,C_start,ita,r2)
            
        self.s2p_out = self.pointscube.transpose(1, 2, 0)
    
    
    def __call__(self, save=True):
        self.gen_s2p_out()
        self.gen_dict_out()
        # save to hdf5 file
        if save:
            self.save()
          #

# Cell
if __name__ == '__main__':
    import time
    import numpy as np
    from glob import glob
    from scipy.stats import binned_statistic_2d
    from astropy.io import fits
    from astropy import wcs
    import h5py
    import matplotlib.pyplot as plt
    from scipy import interpolate
    from PIL import Image
    from scipy import ndimage
    import copy
    from mpl_toolkits.mplot3d import Axes3D 
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord  # High-level coordinates
    from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
    from astropy.coordinates import Angle, Latitude, Longitude  # Angles
    import astropy.units as u
    ##########################################
    args_ = parser.parse_args()
    # print(parser.format_help())
    # print("----------")
    print(parser.format_values())  # useful for logging where different settings came from
    print(args_)
    io = IO(args_)
    io()

