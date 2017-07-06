# -*- coding: utf-8 -*-

from lxml import etree as ET
import numpy as np
import os
from copy import deepcopy

from natsort import natsorted
from scipy import interpolate


'''
 Adrian Jacobo - June 2017
 Fix image names in TrackEm2 exported files by matching to the image names.
'''

dataset_name='N1.5'

# Path to the TrackEM2 traces
TEM2_path ='/Volumes/Eliot_Datasets/Original_Datasets/'+dataset_name+'/N1.5_TRAKEM/N1.5/'

# Path to the Reconstruct traces
img_path = TEM2_path #'/Volumes/Eliot_Datasets/Original_Datasets/'+dataset_name+'/'

# Output path
outpath ='/Volumes/Eliot_Datasets/Original_Traces/'+dataset_name+'/'+dataset_name+'_FromTEM2/'

# Find all the reconstruct files in a dataset
#TEM2_filenames=natsorted([i for i in os.listdir(TEM2_path) if (('newSeries export' in i) and ('.ser' not in i)) ])
#R_filenames=natsorted([i for i in os.listdir(R_path) if (('newSeries export' in i) and ('.ser' not in i)) ])

TEM2_filenames=natsorted([i for i in os.listdir(TEM2_path) if ((dataset_name in i) and ('.ser' not in i) and ('tif' not in i) ) ])
#img_filenames=natsorted([i for i in os.listdir(img_path) if ((dataset_name in i) and ('.tif'  in i))])
img_filenames=natsorted([i for i in os.listdir(img_path) if ('.tif'  in i)])

nfiles=len(TEM2_filenames)

for z in range(nfiles):

    TEM2_name=TEM2_filenames[z]

    # Read the TEM2 exported file for slice z and get the root of the xlm tree
    TEM2_tree = ET.parse(TEM2_path+TEM2_name)
    TEM2_root = TEM2_tree.getroot()


    TEM2_img_tag = TEM2_root.find('./Transform/Image')
    TEM2_img_tag.set('src',img_filenames[z]) #Replace the image name by the one in the list of image names
    TEM2_img_tag.set('brightness','0') #Replace the brightness by the one in the reconstruct file
    TEM2_img_tag.set('contrast','1') #Replace the contrast by the one in the reconstruct file

    #TEM2_section_index =  TEM2_root.attrib['index']
    #TEM2_root.set('index',str(int(TEM2_section_index)+1))

    N_tree = deepcopy(TEM2_tree)
    N_root = N_tree.getroot()

    for transform in N_root.findall('./Transform'):
        N_root.remove(transform)

    N_root.append(deepcopy(TEM2_img_tag.getparent()))

    for contour in TEM2_root.findall('./Transform/Contour'):
        name = contour.attrib['name']
        ct_output=''
        transform = contour.getparent()
        new_transform=ET.SubElement(N_root,'Transform',attrib=transform.attrib)
        new_contour=ET.SubElement(new_transform,'Contour',attrib=contour.attrib)
        lpoints=contour.attrib['points'].split(',')[0:-1] #take up to the second to last element to eliminate remining spaces
        lpoints = np.array([map(float,i.split()) for i in lpoints])

        # Find distances between points, if they are too small interpolate to smooth the curves
        diff = np.diff(lpoints,axis=0)
        dist = np.linalg.norm(diff,axis=1)
        # Interpolate curves and replace them into lpoints
        if (0 < np.min(dist)<0.005):
            # append the starting x,y coordinates
            x=lpoints[:,0]
            y=lpoints[:,1]
            print contour.attrib['name'],len(x)
            x = np.r_[x, x[0]]
            y = np.r_[y, y[0]]
            # fit splines to x=f(u) and y=g(u), treating both as periodic. also note that s=0
            # is needed in order to force the spline fit to pass through all the input points.
            tck, u = interpolate.splprep([x, y], s=0, per=True)
            n=len(lpoints)/2
            print ' interpolated:',len(lpoints),'-->',n
            # evaluate the spline fits for 1000 evenly spaced distance values
            xi, yi = interpolate.splev(np.linspace(0, 1, n), tck)
            lpoints=np.stack((xi,yi),axis=1)

        for i in range(lpoints.shape[0]):
            lpoints[i,0]=lpoints[i,0]
            lpoints[i,1]=lpoints[i,1]
            ct_output=ct_output+'%1.4f' % lpoints[i,0]+' '+'%1.4f' % lpoints[i,1]+','+'\n'
        new_contour.set('points',ct_output)

    filestring = ET.tostring(N_tree,pretty_print=True,
                      xml_declaration=True, encoding=TEM2_tree.docinfo.encoding,
                      doctype=TEM2_tree.docinfo.doctype)

    filestring = filestring.replace('><','> \n <').replace('&#10;','\n \t')

    outname= ''.join(TEM2_name.split(' export',2))
    print TEM2_name,'---->',outname,'. Linked: ',img_filenames[z]
    f=open(outpath+outname,'w')
    f.write(filestring)
    f.close()
