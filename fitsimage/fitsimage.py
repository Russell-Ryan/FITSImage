import numpy as np
from astropy.wcs import WCS



from astropy.io import fits
from astropy import nddata
from scipy import ndimage
#https://docs.scipy.org/doc/scipy-1.3.0/reference/


# https://docs.astropy.org/en/stable/nddata/index.html
class FITSImage(WCS):
    def __init__(self,*args):
        n=len(args)
        if n==1:
            if isinstance(args[0],FITSImage):
                self.loadData(args[0].image,args[0].header)
        elif n==2:
            if isinstance(args[0],str):
                self.loadFile(*args)
            elif isinstance(args[0],fits.hdu.hdulist.HDUList):
                self.loadHDU(*args)
            elif isinstance(args[0],np.ndarray):
                self.loadData(*args)
            else:
                raise NotImplementedError
        else:
            pass
            

    
    #def __init__(self,image=None,header=None):
    #    if image is not None and header is not None:
    #        self.loadData(image,header)

    def loadFile(self,filename,exten):
        ''' read a fits file by exten '''
        with fits.open(filename) as hdul:
            self.loadHDU(hdul,exten=exten)
            #image=hdul[exten].data
            #header=hdul[exten].header
        #self.load(image,header)
        
    def loadHDU(self,hdul,exten):
        image=hdul[exten].data
        header=hdul[exten].header
        self.loadData(image,header)
        
    def loadData(self,image,header):
        ''' load the data into the object '''
        self.image=image
        self.header=header
        WCS.__init__(self,self.header)

    def astrom(self):
        pass
    #https://reproject.readthedocs.io/en/stable/api/reproject.reproject_interp.html#reproject.reproject_interp


    def extract(self,x0,x1,y0,y1,**kwargs):
        ''' extract a sub array.  equivalent to hextract '''

        setLTV=lambda k,v: hdr[k]-v if k in hdr else -v
        
        # use astropy to cut
        shape=(y1-y0+1,x1-x0+1)
        position=((y0+y1+1)/2.,(x0+x1+1)/2.)
        img=nddata.extract_array(self.image,shape,position,**kwargs)

        # update the new header keywords
        hdr=self.header.copy()
        hdr['NAXIS1']=img.shape[1]
        hdr['NAXIS2']=img.shape[0]
        hdr['CRPIX1']-=x0
        hdr['CRPIX2']-=y0
        hdr['LTV1']=setLTV('LTV1',x0)
        hdr['LTV2']=setLTV('LTV2',y0)
        
        # add new keywords
        hdr['XMIN']=(x0,'lower x-bound')
        hdr['XMAX']=(x1,'upper x-bound')
        hdr['YMIN']=(y0,'lower y-bound')
        hdr['YMAX']=(y1,'upper y-bound')
        
        # update the history
        history='Extracted from region (x,y)=[{}:{},{}:{}]'.format(x0,x1,y0,y1)
        hdr.add_history(history)

        # create the output
        output=type(self)(img,hdr)

        return output

    def rebin(self,binfact):
        if np.isscalar(binfact):
            xbin=binfact
            ybin=binfact
        else:
            xbin=binfact[0]
            ybin=binfact[1]
        xratio=1./xbin
        yratio=1./ybin

        lamb=yratio/xratio
        pixratio=xratio*yratio
            
        
        # trim the image first (to deal with edge effects
        x=xbin*(self.shape[1]//xbin)-1
        y=ybin*(self.shape[0]//ybin)-1

        sub=self.extract(0,x,0,y)
        
        # block average the image
        img=nddata.block_reduce(sub.image,(xbin,ybin),func=np.average)
                
        # get the header
        hdr=sub.header.copy()

        # update the main the WCS values
        hdr['CRPIX1']=(hdr['CRPIX1']-1.)*xratio+1.
        hdr['CRPIX2']=(hdr['CRPIX2']-1.)*yratio+1.
        if 'CDELT1' in hdr:
            hdr['CDELT1']=hdr['CDELT1']/xratio
        if 'CDELT2' in hdr:
            hdr['CDELT2']=hdr['CDELT2']/yratio

    
        
        print("Need to sort out LTVs")
        #if 'LTV1' in hdr:
        #    hdr['LTV1']=hdr['LTV1']*xratio
        #if 'LTV2' in hdr:
        #    hdr['LTV2']=hdr['LTV2']*yratio
        hdr['CD1_1']=hdr['CD1_1']/xratio
        hdr['CD1_2']=hdr['CD1_2']/xratio
        hdr['CD2_1']=hdr['CD2_1']/yratio
        hdr['CD2_2']=hdr['CD2_2']/yratio
        
        # update the distortion
        if self.sip is not None:
            print("HELP, must update SIP")
        

        # update the BSCALE
        if not isinstance(img.dtype,np.unsignedinteger):
            if 'BSCALE' in hdr and hdr['BSCALE']!=1 and hdr['BSCALE']!=0:
                hdr['BSCALE']=hdr['BSCALE']/pixratio
            if 'BZERO' in hdr and hdr['BZERO']!=0:
                hdr['BZERO']=hdr['BZERO']/pixratio
                
        
        # update the history
        history='Block averaged image with factor {}'.format(binfact)
        hdr.add_history(history)


        # create the output
        output=type(self)(img,hdr)

        return output

        
    def rot(self):
        #https://docs.scipy.org/doc/scipy-0.16.0/reference/ndimage.html
        pass

    def rotate(self):
        pass

    def reverse(self):
        pass

    
    def writefits(self,filename,overwrite=True):
        hdul=fits.HDUList()
        hdul.append(self.ImageHDU())
        hdul.writeto(filename,overwrite=overwrite)
        hdul.close()
        
    
    def ImageHDU(self,**kwargs):
        return fits.ImageHDU(data=self.image,header=self.header,**kwargs)
    

    
    def __getitem__(self,k):
        if isinstance(k,str):
            return self.header[k]
        elif isinstance(k,tuple):
            return self.image[k[0],k[1]]
        else:
            raise NotImplementedError(type(k))

    def __setitem__(self,k,v):
        if isinstance(k,str):
            self.header[k]=v
        elif isinstance(k,tuple):
            self.image[k[0],k[1]]=v
        elif isinstance(k,np.ndarray):
            try:
                print(k.dtype)
                self.image[k]=v
            except:
                print('keyword is invalid numpy datatype')
        else:
            print('Invalid keyword type',type(k))
                

            
    def __sub__(self,val):
        return self.image-val

    def __rsub__(self,val):
        return val-self.image

    def __eq__(self,val):
        return self.image == val

    def __ne__(self,val):
        return self.image != val

    def __gt__(self,val):
        return self.image > val

    def __ge__(self,val):
        return self.image >= val
        
    def __lt__(self,val):
        return self.image < val
    
    def __le__(self,val):
        return self.image <= val

    def __contains__(self,k):
        return k in self.header
    
    @property
    def shape(self):
        return self.image.shape

    @property
    def dtype(self):
        return self.image.dtype
    

    # convenience functions to avoid that dumb 0
    def xy2ad(self,x,y):
        return self.all_pix2world(x,y,0)

    def ad2xy(self,a,d):
        return self.all_world2pix(a,d,0)

    def xy2xy(self,x,y,obj):
        a,d=self.all_pix2world(x,y,0)
        return self.all_world2pix(a,d,0)
      
            
        
if __name__=='__main__':

        
    filename='/Users/rryan/icoi3immq_flt.fits'
    exten=1
    img=FITSImage(filename,exten)
    
    sub=img.extract(490,530,490,510)
    sub.writefits('sub.fits')
    binned=sub.rebin(3)

    binned.writefits('python.fits')
    
    
    x=np.array([400,500])
    y=np.array([400,500])
    #x=np.array([500])
    #y=np.array([500])
    #x=500
    #y=500
    
