#!/usr/bin/env python
# This program search for correlations between 2 fits catalogues
# Author: Evandro M. Ribeiro (evandromartinezribeiro@gmail.com)

# modules for Python 2 and 3 compatibility
from __future__ import absolute_import, division, print_function
from builtins import input  # installation: pip install future

# modules for this script
import time     # time options library
import pprocess # paralel python library
import math     # math operations
import numpy as np # numerical recipes for arrays
import pyfits as fits # FITS manipylation module
#from astropy.io import fits    #FITS manipulation module
# transition to astropy in progress, problems with NULL and NaN values


#------------------------------------------------------------------------------
#This function make the search in cat2
def compare(start, nrows,rows, tol, ra1, ra2, de1, de2,rows1):
    found = 0
    result = np.array([])
    if start+nrows > rows:
        nrows = rows - start

    #to search the whole cat1 use: 'in range(rows1)'
    # or 'in range ([start line], [end line])'
    for i in range(rows1):
        for j in range(start,nrows+start):
            #Angular distance between 2 objects
            #General Formula
            degtorad = math.pi/180.0
            print(de1[i]*degtorad, math.radians(de1[i]))
            tmp1 = math.sin(math.radians(de1[i]))*math.sin(math.radians(de2[j]))
            tmp2 = math.cos(math.radians(de1[i]))*math.cos(math.radians(de2[j]))
            tmp3 = math.cos(math.radians(ra1[i] - ra2[j]))
            tmp = tmp1 + tmp2*tmp3
            dist1 = math.degrees(math.acos(tmp))*3600.0
            #More precise formula to small distances
            if (abs(dist1) < 600.0):
                x = (0.5*(de1[i]+de2[j]))
                tmp4 = (ra1[i]-ra2[j])*math.cos(math.degrees(x))
                dist2 = math.sqrt(tmp4**2 + (de1[i]-de2[j])**2)*3600.0
            elif (abs((180.0*3600.0)-dist1) < 600.0):
                x = (0.5*(de1[i]+de2[j]))
                tmp4 = (ra1[i]-ra2[j])*math.cos(math.degrees(x))
                dist2 = math.sqrt(tmp4**2 + (de1[i]-de2[j])**2)*3600.0
            else:
                dist2 = dist1
            #dist2 = dist1
            #print(dist1, dist2)
            if dist2 < tol:
                found = found+1
                #vector results keep: found, i,j and dist
                a = np.array(found)
                a = np.append(a,i)
                a = np.append(a,j)
                a = np.append(a,dist2)
                result = np.append(result, a)
        if found == 0:
            return 0
        else:
            return result

cat1name = str(input('Type the name of your smallest catalogue with the extension (.fits or .fit): '))
ra1col = int(input('Type the collumn with the Right Assencion in the catalogue: '))
de1col = int(input('Type the collumn with the Declination in the catalogue: '))
cat2name = str(input('Type the name of your bigest catalogue with the extension (.fits or .fit): '))
ra2col = int(input('Type the collumn with the Right Assencion in the catalogue (to 2XMMiDR3 is 16): '))
de2col = int(input('Type the collumn with the Declination in the catalogue (to 2XMMiDR3 is 17): '))

#--EDIT HERE---------------------------------------------------------------------------------------------
cat1 = fits.open(cat1name)              # The small catalogue is cat1
rows1 = cat1[1].header['NAXIS2']
print(rows1)
cat2 = fits.open(cat2name)              # The big catalogue is cat2
rows = cat2[1].header['NAXIS2']
ra1 = cat1[1].data.field(ra1col)        #the number of the column with the RA position on cat1 in the field argument (start from zero)
de1 = cat1[1].data.field(de1col)        #the number of the column with the Declination position on cat1 in the field argument (start from zero)
ra2 = cat2[1].data.field(ra2col)        #the number of the column with the RA coordinate on cat2 in the field argument (to 2XMMiDR3 is 16)
de2 = cat2[1].data.field(de2col)        #the number of the column with the Declination coordinate on cat2 in the field argument (to 2XMMiDR3 is 17)

#ATENTION !!!  the Follow block will define the output file format
#--------------------------------
t1 = cat1[1].columns    #wich columns of cat1 you want in your output file (use .columns[start:end] to define)
for i in xrange(len(t1)):
    t1[i].name += '1'
t1 = t1[0:2]
t2 = cat2[1].columns    #wich columns of cat2 you want in your output file (use .columns[start:end] to define)
for i in xrange(len(t2)):
    t2[i].name += '2'
t2 = t2[0:2]
t3 = fits.Column(name='Distance',format='E', unit='arcsec')     #Create a extra column with the distance between the found objects
t4 = fits.Column(name='ID',format='E', unit='number')           #Create a extra column with an ID number to the found objects
t=t1+t2+t3+t4

totalcol = len(t)
coldefs = fits.ColDefs(t)

tol = float(input('Choose your search tolerance in arcsec:'))       #your radii search's tolerance in arcsec
ncpus = int(input('How many processor will you use?:'))       #number of processor to run te search in paralel

#---------------------------------------------------------------------------------------------
#auxiliar variables and inicializations
nrows = 0
ignore = 0

#-----------------------------------------------------------------------------
# pprocess stuff
presults = pprocess.Map(limit=ncpus, reuse=1)
parallel_correlation = presults.manage(pprocess.MakeReusable(compare))
#---------------------------------------------------------------------------------------------
#divide the jobs acording to number of processors
while nrows < 1:
        if rows % ncpus == 0:
                nrows = rows//ncpus
        else:
                rows=rows-1
                ignore=ignore+1

inputs = range(0,rows+ignore,nrows)

print("\n Initializing the process with", ncpus, "processor(s)\n")
inputs = range(0,rows+ignore,nrows)
tic = time.time()
[parallel_correlation(inpute, nrows, rows+ignore, tol, ra1, ra2, de1, de2, rows1) for inpute in inputs]

results = []
aux=0
for presult in presults :
    if type(presult) != type(1):
        if len(presult) > 1:
            for pvalue in presult:
                results.append(pvalue)
        else:
            results.append(presult)

results = np.array(results)
print('time = {0}'.format(time.time() - tic))

#---------------------------------------------------------------------------------------------
#create and write the output.fits file
#ATENTION see below the correct fields to aux array to match the columns define in the startup
results = np.reshape(results, (-1,4))
if results.shape[0] == 0:
        print('\n\n Sorry we did not find corresponding objects\n\n')
        pass
else:
        new = fits.TableHDU.from_columns(coldefs, nrows=results.shape[0])
        #new = fits.new_table(coldefs, nrows=results.shape[0])
        for k in xrange(results.shape[0]):
                line1 = results[k][1]
                line2 = results[k][2]
                aux = np.array(cat1[1].data[line1][0:2])        #all columns from cat1 (else to limited columns, use[line1][start:end])
                aux = np.append(aux,cat2[1].data[line2][0:2])   #all columns from cat2 (else to limited columns, use[line1][start:end])
                aux = np.append(aux,results[k][3])              #Distance between the objects
                aux = np.append(aux,line1)                      #ID of the result
                a=[]
                for m in xrange(totalcol):
                        a.append(aux[m])
                new.data[k] = a
        new.writeto('output.fits')
        print('\n\n We found', results.shape[0],'corresponding objects.')
        print(' Please, see the output.fits file to results\n\n')

#--------------------------------------------------------------------------------------------
cat2.close()
cat1.close()
