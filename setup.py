from setuptools import setup,find_packages


ld='A class to facilitate various operations for fits images with astrometric headers, modeled after the h* routines in IDL astrolib'



#https://docs.python.org/2/distutils/setupscript.html
setup(name='fitsimage',
      version='1.0',
      author='Russell Ryan',
      author_email='rryan@stsci.edu',
      keywords='image processing fits header astrometric',
      description='Processing images with astrometry',
      long_description=ld,
      maintainer='Russell Ryan',
      license='MIT',
      url='https://github.com/Russell-Ryan/fitsimage',
      platforms='posix',
      install_requires=['numpy','astropy','pkginfo'],
      classifiers=['Development Status :: 3 Alpha',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Astronomy',],
      packages=find_packages())
