#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 11:45:23 2017

@author: vinay
"""
import os
import numpy as np
import math
from scipy import ndimage
import gdal

def water_estimation(loc, data_list, W_para, b_para):
    
    lis = open(loc + data_list, 'r')
    lis = lis.readlines()
    
    for i in range(len(lis)):
        dirctory = lis[i].strip()
        location= dirctory
        print (location)
    
        # reading the Neaural network percepetron derived parametes
        W=np.load(W_para)
        b=np.load(b_para)
        # for import raster
        os.chdir(loc + dirctory)
        
        data = {}   
        print ("Data:", dirctory)                        
        fl = open(dirctory + "_MTL.txt", 'r')
        file_lines = fl.readlines()
        
        for line in file_lines:
           values = line.split(' = ')
           if len(values) != 2:
               continue
        
           data[values[0].strip()] = values[1].strip().strip('"')
    
        Ml_rad = float(data['RADIANCE_MULT_BAND_10'])
        Ad_rad = float(data['RADIANCE_ADD_BAND_10'])
        k1 = float(data['K1_CONSTANT_BAND_10'])
        k2 = float(data['K2_CONSTANT_BAND_10'])
        
        bands = []
        for i in range(2,11):
            if i == 8 or i ==9:
                continue
            print ("B"+str(i))
            
            band = gdal.Open(dirctory + '_B'+ str(i)+".TIF")
            band = band.ReadAsArray()
            print (band.shape)
    
            if i <=7:
                Ml= float(data['REFLECTANCE_MULT_BAND_'+ str(i)])
                Ad = float(data['REFLECTANCE_ADD_BAND_'+ str(i)])
                print Ml
                print Ad
                band = band * Ml + Ad
    
            else:
                band = np.where((band > 0), (k2 / (np.log(k1/(Ml_rad*band + Ad_rad) +1))
                                   - 273.15), 0)
                print ("Temprature band" + str(band.shape))
            bands.append(band)
            
        Bands = np.array(bands)       
        print ("final shape" + str(Bands.shape))    
        Blue =  Bands[0]   
        Green = Bands[1]
        Red = Bands[2]
        NIR = Bands[3]
        SWIR1 = Bands[4]
        SWIR2 = Bands[5]
    
    
        A = Green - NIR
        B = Blue - NIR
        C = Red - SWIR1
        D = SWIR1
        E = SWIR2
        
        print ("calculating non-water(w1) and water (w2)")
        
        w1 = W[0,0] * A + W[0,1] * B + W[0,2]  * C + W[0,3] * D + W[0,4] * E + b[0]
        w2 = W[1,0] * A + W[1,1] * B + W[1,2]  * C + W[1,3] * D + W[1,4] * E + b[1]
        wat_nonwat_soft = np.exp(w2)/(np.exp(w1)+np.exp(w2))
        test = np.where(w2>w1,1,0)  
         
        # clculate relative azimuth from the inclination of the scene
        # ax,ay, lower point values for x and y
        # bx,by, upper point values for x and y
    
        image = Blue
        for y in range(image.shape[0]):
           for x in range(image.shape[1]):
              if (image[y,x] > 0.0):
                 print image[y,x]
                 b = [y,x]
                 print (y,x)
                 break
            
           if (image[y,x] > 0.0):    
              break
        image_tr = np.transpose(image,(1,0))
        for x in range(image_tr.shape[0]):
           for y in range(image_tr.shape[1]):
               if (image_tr[x,y] > 0.0):
                 print image_tr[x,y]
                 a = [x,y]
                 print (x,y)
                 break
            
           if (image_tr[x,y] > 0.0):    
              break
          
        bx = b[1] * 30
        ax = a[0] * 30
        by = b[0] * 30
        ay = a[1] * 30
        print ("%s,%f,%s,%f,%s,%f,%s,%f" % ("bx:",bx,"by:",by,"ax:",ax,"ay:",ay))
        dx= bx - ax
        dy= ay - by
        print ("%s,%f,%s,%f" % ("dx:",dx,"dy:",dy))
        m= math.sqrt((dx**2)+(dy**2))
        r= math.atan2(dx/m, dy/m)
        print r
        d= (r*180)/(math.pi)
        print ("d",d)
        sun_az = float(data['SUN_AZIMUTH']) 
        rel_az_pv=(sun_az-d) + 90
        #convert to radian
        rel_az_pv = rel_az_pv * math.pi /180
        rel_az_nv=sun_az -(d + 90)
        #convert to radian
        rel_az_nv = rel_az_nv * math.pi /180
        #convert to radian
        sol_z = 90 - float(data['SUN_ELEVATION']) 
        sol_z= sol_z * math.pi/180
        print np.cos(sol_z)
        
        # Fmask to generate pixel-wise sensor zenith     
        os.system("gdal_merge.py -separate -of HFA -co COMPRESSED=YES -o ref.img L*_B[1-7,9].TIF")
        os.system("gdal_merge.py -separate -of HFA -co COMPRESSED=YES -o thermal.img L*_B1[0,1].TIF")
        os.system("fmask_usgsLandsatMakeAnglesImage.py -m *_MTL.txt -t ref.img -o angles.img")
    
        # import fmask senith angle to grass gis 
        fmask_angl = gdal.Open("angles.img")
        fmask_angl = fmask_angl.ReadAsArray()
        fmask_angl = fmask_angl[1]
        sen_z = fmask_angl * 3.14/180
    
        print (sen_z.shape[0],sen_z.shape[1]) 
        print sen_z.shape[0]
        
        #create posetive and negative patches
        imgR = ndimage.rotate(sen_z, d, reshape=False)
        del sen_z
        div  = imgR.shape[1]/2
        print ("div"+str(div))
        negative = np.arange((imgR.shape[0]) * div).reshape(imgR.shape[0],div)
        if (imgR.shape[1] % 2 != 0):
            posetive = np.arange((imgR.shape[0]) * (div+1)).reshape(imgR.shape[0],(div+1))
        else:
            posetive = np.arange((imgR.shape[0]) * div).reshape(imgR.shape[0],div)
        negative = negative.astype(np.float32)
        print ("negative: ", str(negative.shape))
        posetive = posetive.astype(np.float32)
        posetive = np.zeros_like(posetive)
        print ("posetive: ", str(posetive.shape))
        neg_y = negative.shape[1]
        #pos_y = posetive.shape[1]
        negative[:,:] =  (np.arccos(np.cos(sol_z)*np.cos(imgR[:,:neg_y])-
                                    np.sin(sol_z)* np.sin(imgR[:,:neg_y]) * np.cos(rel_az_nv)) * 180/np.pi)  
        posetive[:,:] =  (np.arccos(np.cos(sol_z)*np.cos(imgR[:,neg_y:])-
                                    np.sin(sol_z)* np.sin(imgR[:,neg_y:]) * np.cos(rel_az_pv)) * 180/np.pi)
    
        ps_ng =np.concatenate((negative,posetive), axis=1)
        del negative
        del posetive
                                    
        # rotate back to same angle
        rotate_bak = ndimage.rotate(ps_ng, -(d), reshape=False)
        sun_gl = rotate_bak
        wmask_glnt1 = np.where(sun_gl < 20, wat_nonwat_soft +(1/sun_gl), 0)
        wmask_glnt2 = np.where((sun_gl > 20) & (sun_gl < 35), wat_nonwat_soft +(1/(sun_gl*2)),0)
        wmask_glnt3 = np.where(sun_gl >=35, wat_nonwat_soft +(1/(sun_gl*3)),0)
        wmask_glnt = wmask_glnt1 + wmask_glnt2 + wmask_glnt3
        wmsk_glnt_crct = np.where(wmask_glnt > 0.5, 1, 0) 
        
        BT = Bands[6] 
        MNDWI = (Green-SWIR1)/(Green+SWIR1)
        NDWI= (Green-NIR)/ (Green + NIR)
        snow = np.where((wmsk_glnt_crct == 1) & (MNDWI > NDWI + 0.7) & (BT < 8), 1, 0)
        water_non_water = np.where(snow == 1, 0,wmsk_glnt_crct)
        # convert nan value region to zero
        water_non_water = np.where(Bands[1] != -0.1,water_non_water,0 ) 
        np.save(dirctory + 'water_non_water.npy', water_non_water)# water estimated saved as numpy file
        
        os.chdir(loc)
        
        try:
    	     os.remove(loc+dirctory+'/ref.img')
    	     os.remove(loc+dirctory+'/angles.img')
        except:
            pass
        return (water_non_water)
if __name__ == "__main__":

    loc=('Landsat_data/')
    data_list= 'list_.txt'
    W_para = 'W_best_initpretrained.npy'
    b_para = 'b_best_initpretrained.npy'
    water_non_water = water_estimation(loc,data_list,W_para,b_para)
