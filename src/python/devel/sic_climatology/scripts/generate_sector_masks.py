import cdms2
import MV2 as MV
import string
import os,sys

from pcmdi_metrics.pcmdi.pmp_parser import *

#Factors
factor1=1.e-6   #model units are m^2, converting to km^2
factor2=1.e-2   #model units are %, converting to non-dimen
a=6371.009      #Earth radii in km
pi=22./7.
factor3=4.*pi*a*a       # Earth's surface area in km2
dc=0.15         #minimum ice concentration contour

pin = '/work/gleckler1/processed_data/cmip5clims-historical/sic/cmip5.MOD.historical.r1i1p1.mo.seaIce.OImon.sic.ver-1.1980-2005.SC.nc'

pins = string.replace(pin,'MOD','*')

lst = os.popen('ls ' + pins).readlines()

mods = []
mods_failed = []
for l in lst:
  mod  = string.split(l,'.')[1]
  if mod not in mods: mods.append(mod)

#w =sys.stdin.readline()

var = 'sic'
factor2 = 1

mods = ['ACCESS1-3']

for mod in mods:
 try:
   fc = string.replace(pin,'MOD',mod)
   f=cdms2.open(fc)
   sic=f(var, squeeze=1)
   sic_grid=sic.getGrid()
   lat=sic.getLatitude()
   lon=sic.getLongitude()
   sic=MV.multiply(sic,factor2)

   print 'CMIP5-native= ',MV.max(sic)
   if MV.rank(lat)==1:
         tmp2d=f(var,time=slice(0,1),squeeze=1)
         lats=MV.zeros(tmp2d.shape)
         for ii in range (0,len(lon)):
           lats[:,ii]=lat[:]
   else:
           lats=lat

   if MV.rank(lon)==1:
         tmp2d=f(var,time=slice(0,1),squeeze=1)
         lons=MV.zeros(tmp2d.shape)
         for ii in range (0,len(lat)):
           lons[ii,:]=lon[:]
   else:
           lons=lon

   f.close()


#######################################################
### areacello
   area_dir = '/work/cmip5/fx/fx/areacello/'
   alist = os.listdir(area_dir) # LIST OF ALL AREACELLO FILES

   for a in alist:
     if string.find(a,mod) != -1: 
      areaf = a 
      print mod,' ', a
#  w = sys.stdin.readline()

   g = cdms2.open(area_dir + areaf)

   try:
    area = g('areacello')
   except:
     area = g('areacella')
   area = MV.multiply(area,factor1)
   area = MV.multiply(area,factor1)

   g.close()

   land_mask=MV.zeros(area.shape)

# Reading the ocean/land grid cell masks (sftof/sftlf)
   mask_dir = '/work/cmip5/fx/fx/sftof/'
   mlist = os.listdir(mask_dir) # LIST OF ALL AREACELLO FILES

   for m in mlist:
     if string.find(m,mod) != -1:
      maskf = m
      print mod,' ', m

   s = cdms2.open(mask_dir + maskf)
   try:
    frac = s('sftof')
    if (mod != 'MIROC5' and mod!= 'GFDL-CM2p1' and mod!='GFDL-CM3' and mod!='GFDL-ESM2M'):
       frac = MV.multiply(frac,factor2)
       print mod,MV.max(frac)
    area = MV.multiply(area,frac)
    land_mask = MV.multiply(1,(1-frac))
   except:
      frac = s('sftlf')
   if (mod != 'MIROC5' or mod!= 'GFDL-CM2p1' or mod!='GFDL-CM3' or mod!='GFDL-ESM2M'):
     frac = MV.multiply(frac,factor2)
   area = MV.multiply(area,(1-frac))
   land_mask = MV.multiply(1,frac)
   s.close()

# Creating regional masks
# GFDL and bcc model grids have shift of 80
   if (mod[0:4] == 'GFDL' or mod[0:3] == 'bcc') :
           lons_a=MV.where(MV.less(lons,-180.),lons+360.,lons)
           lons_p=MV.where(MV.less(lons,0.),lons+360.,lons)
   else:
           lons_a=MV.where(MV.greater(lons,180.),lons-360.,lons)
           lons_p=lons

   print 'CMIP5'
   print 'lons_na= ',MV.min(lons_a),MV.max(lons_a)
   print 'lons_np= ',MV.min(lons_p),MV.max(lons_p)
#  mask_arctic=MV.zeros(area.shape)
#      mask_antarctic=MV.zeros(area.shape)
   mask_ca=MV.zeros(area.shape)
   mask_na=MV.zeros(area.shape)
   mask_np=MV.zeros(area.shape)
   mask_sa=MV.zeros(area.shape)
   mask_sp=MV.zeros(area.shape)
   mask_io=MV.zeros(area.shape)

###############################################################
#Arctic Regions
#Central Arctic
   lat_bound1=MV.logical_and(MV.greater(lats,80.),MV.less_equal(lats,90.))
   lat_bound2=MV.logical_and(MV.greater(lats,65.),MV.less_equal(lats,90.))
   lon_bound1=MV.logical_and(MV.greater(lons_a,-120.),MV.less_equal(lons_a,90.))
   lon_bound2=MV.logical_and(MV.greater(lons_p,90.),MV.less_equal(lons_p,240.))
   reg1_ca=MV.logical_and(lat_bound1,lon_bound1)
   reg2_ca=MV.logical_and(lat_bound2,lon_bound2)
   mask_ca=MV.where(MV.logical_or(reg1_ca,reg2_ca),1,0)
   mask_ca=MV.where(MV.equal(land_mask,0),0,mask_ca)           # 0 - Land

   g = cdms2.open('data/mask_' + mod + '_ca.nc','w+')
   mask_ca.id = 'mask'
   g.write(mask_ca)
#  land_mask.id = 'sftof'
#  g.write(land_mask)
   g.close()

#NA region
   lat_bound=MV.logical_and(MV.greater(lats,45.),MV.less_equal(lats,80.))
   lon_bound=MV.logical_and(MV.greater(lons_a,-120.),MV.less_equal(lons_a,90.))
   lat_bound3=MV.logical_and(MV.greater(lats,45.),MV.less_equal(lats,50.))
   lon_bound3=MV.logical_and(MV.greater(lons_a,30.),MV.less_equal(lons_a,60.))
   reg3=MV.logical_and(lat_bound3,lon_bound3)

   mask_na=MV.where(MV.logical_and(lat_bound,lon_bound),1,0)
   mask_na=MV.where(MV.equal(reg3,True),0,mask_na)   # Masking out the Black and Caspian Seas
   mask_na=MV.where(MV.equal(land_mask,True),0,mask_na)           # 0 - Land
   mask_na=MV.where(MV.equal(land_mask,0),0,mask_na)           # 0 - Land

   g = cdms2.open('data/mask_' + mod + '_na.nc','w+')
   mask_na.id = 'mask'
   g.write(mask_na)
   g.close()
   
#NP region
   lat_bound=MV.logical_and(MV.greater(lats,45.),MV.less_equal(lats,65.))
   lon_bound=MV.logical_and(MV.greater(lons_p,90.),MV.less_equal(lons_p,240.))
   mask_np=MV.where(MV.logical_and(lat_bound,lon_bound),1,0)
   mask_np=MV.where(MV.equal(land_mask,0),0,mask_np)           # 0 - Land

   g = cdms2.open('data/mask_' + mod + '_np.nc','w+')
   mask_np.id = 'mask'
   g.write(mask_np)
   g.close()


#Antarctic Regions
   lat_bound=MV.logical_and(MV.greater(lats,-90.),MV.less_equal(lats,-55.))

#SA region
   lon_bound=MV.logical_and(MV.greater(lons_a,-60.),MV.less_equal(lons_a,20.))
   mask_sa=MV.where(MV.logical_and(lat_bound,lon_bound),1,0)
   mask_sa=MV.where(MV.equal(land_mask,0),0,mask_sa)           # 0 - Land

   g = cdms2.open('data/mask_' + mod + '_sa.nc','w+')
   mask_sa.id = 'mask'
   g.write(mask_sa)
   g.close()

#SP region
   lon_bound=MV.logical_and(MV.greater(lons_p,90.),MV.less_equal(lons_p,300.))
   mask_sp=MV.where(MV.logical_and(lat_bound,lon_bound),1,0)
   mask_sp=MV.where(MV.equal(land_mask,0),0,mask_sp)           # 0 - Land

   g = cdms2.open('data/mask_' + mod + '_sp.nc','w+')
   mask_sp.id = 'mask'
   g.write(mask_sp)
   g.close()

#IO region
   lon_bound=MV.logical_and(MV.greater(lons_p,20.),MV.less_equal(lons_p,90.))
   mask_io=MV.where(MV.logical_and(lat_bound,lon_bound),1,0)
   mask_io=MV.where(MV.equal(land_mask,0),0,mask_io)           # 0 - Land

   g = cdms2.open('data/mask_' + mod + '_io.nc','w+')
   mask_io.id = 'mask'
   g.write(mask_io)
   g.close()

# Calculate the Total Sea Ice Concentration/Extent/Area
#       ice_area = MV.multiply(sic,1)                                #SIC
   area = MV.multiply(1,area)                              #SIE
   ice_area = MV.multiply(sic,area)                            #SIA
   ice_area = MV.where(MV.greater_equal(sic,0.15),ice_area,0.)                     #Masking out the sic<0.15

#    arctic=MV.logical_and(MV.greater_equal(lats,35.),MV.less(lats,87.2))         #SSM/I limited to 87.2N
   mask_arctic=MV.logical_and(MV.greater_equal(lats,45.),MV.less(lats,90.))         #Adding currently in SSM/I 100% in the area >87.2N
   mask_antarctic=MV.logical_and(MV.greater_equal(lats,-90.),MV.less(lats,-55.))

 except:
  'Failed for ', mod
  mods_failed.append(mod)

print 'failed for ', mods_failed

# Calculate the Sea Ice Covered Area


