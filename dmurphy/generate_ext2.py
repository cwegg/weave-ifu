import astropy.io.fits as pyfits
from numpy import random
import numpy
fname = 'WA_FITS_catalogue_180523-132206.fits'

fits = pyfits.open(fname)


data = fits[1].data
cols = ['TARGID','TEMPLATE','MAG','FILTERID','REDSHIFT','VELOCITY','FWHM','IFU_SPAXEL','DITHER_ID']
fmts = ['9A','50A','E','I','E','E','E','3A','20A']



template = []
_data = {}
_data['TEMPLATE'] = []
columns = []

for c,f in zip(cols,fmts):
    if c == 'TARGID':
        col_data = numpy.array([tid for tid in fits[1].data['TARGID']])
    elif c == 'TEMPLATE':
        col_data = numpy.array(['tmp_%d.%d_%s.dat'%(random.randint(50),random.randint(5),random.choice(['qe','xe','ps','nw','ae'])) for targ in fits[1].data['TARGID']])
    elif c == 'MAG':
        col_data = numpy.array([random.normal(19.0,0.65)  for targ in fits[1].data['TARGID']])
    elif c == 'FILTERID':
        col_data = numpy.array([1 for targ in fits[1].data['TARGID']],dtype=numpy.int)
    elif c == 'REDSHIFT':
        col_data = numpy.array([max(0.0,random.normal(1.2,1.0))  for targ in fits[1].data['TARGID']])
    elif c == 'VELOCITY':
        col_data = numpy.array([random.normal(0.0,150.0)  for targ in fits[1].data['TARGID']])
    elif c == 'FWHM':
        col_data = numpy.array([random.normal(1.0,0.2)  for targ in fits[1].data['TARGID']])
    elif c == 'IFU_SPAXEL':
        col_data = numpy.array([tid for tid in fits[1].data['IFU_SPAXEL']])
    elif c == 'DITHER_ID':
        col_data = numpy.array(['' for tid in fits[1].data['IFU_SPAXEL']])

        
    disp = f
    if 'A' in f:
        disp = 'A'+f.replace('A','')
    column = pyfits.Column(name=c,format=pyfits.column._ColumnFormat(f),disp=disp,array=col_data)
    columns.append(column)


print 'Creating new Coldef'
new_coldef = pyfits.ColDefs(columns)
print 'Creating new table'
#new_table = pyfits.BinTableHDU.from_columns(new_coldef,header=base[1].header)   ##ha - the origin of woe: providing the original header...
new_table = pyfits.BinTableHDU.from_columns(new_coldef)
new_table.update()
fits.append(new_table)

fits.writeto(fname.replace('.fit','ex2.fit'),overwrite=True)
