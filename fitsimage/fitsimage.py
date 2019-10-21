import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from astropy import nddata


# https://docs.astropy.org/en/stable/nddata/index.html
class FITSImage(WCS):
    def __init__(self,image=None,header=None):
        if image is not None and header is not None:
            self.load(image,header)

    def read(self,filename,exten=1):
        ''' read a fits file by exten '''
        with fits.open(filename) as hdul:
            image=hdul[exten].data
            header=hdul[exten].header
        self.load(image,header)
        
    
    def load(self,image,header):
        ''' load the data into the object '''
        self.image=image
        self.header=header
        WCS.__init__(self,self.header)

    def astrom(self):
        pass
        
    def extract(self,x0,x1,y0,y1,**kwargs):
        ''' extract a sub array.  equivalent to hextract '''
        
        # use astropy to cut
        shape=(y1-y0+1,x1-x0+1)
        position=((y0+y1+1)/2.,(x0+x1+1)/2.)
        img=nddata.extract_array(self.image,shape,position,**kwargs)

        # cut via direct image
        #x0=np.clip(x0,0,self.shape[1])
        #x1=np.clip(x1+1,0,self.shape[1])
        #y0=np.clip(y0,0,self.shape[0])
        #y1=np.clip(y1+1,0,self.shape[0])
        #img=self.image[y0:y1,x0:x1]

        # update the new header keywords
        hdr=self.header.copy()
        hdr['NAXIS1']=img.shape[1]+1
        hdr['NAXIS2']=img.shape[0]+1
        hdr['CRPIX1']-=x0
        hdr['CRPIX2']-=y0
        hdr['LTV1']=-x0
        hdr['LTV2']=-y0
        hdr.add_history('Extracted from region (x,y)=[{}:{},{}:{}]'.format(x0,x1,y0,y1))

        
        #FITSImage(img,hdr)     
        
        # output
        return type(self)(img,hdr)

    def rebin(self):
        pass

    def rot(self):
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
        #hdr=self.to_header(relax=True)
        #header=self.header.copy()
        #for h in hdr:
        #    header[h]=hdr[h]
        out=fits.ImageHDU(data=self.image,header=self.header,**kwargs)
                
        return out
    

    
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

        

            
    @property
    def shape(self):
        return self.image.shape

            
        
if __name__=='__main__':
    
    filename='/Users/rryan/icoi3immq_flt.fits'
    exten=1
    x=FITSImage()
    x.read(filename,exten=exten)
    sub=x.extract(490,530,490,510)
    print(sub.wcs.crpix,x.wcs.crpix)

    sub.writefits('t.fits')

        
