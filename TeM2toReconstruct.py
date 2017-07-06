# -*- coding: utf-8 -*-

from lxml import etree as ET
import numpy as np
import os
from copy import deepcopy

from natsort import natsorted

'''
 Adrian Jacobo - May 2017
 Fix image names in TrackEm2 exported files
'''

# Path to the TrackEM2 traces
TEM2_path ='/Volumes/Eliot_Datasets/Original_Datasets/DS1/20170530_TrakEM-export_final/'

# Path to the Reconstruct traces
R_path ='/Volumes/Eliot_Datasets/Original_Traces/DS1/DS1_v10/'

# Output path
outpath ='/Volumes/Eliot_Datasets/Original_Traces/DS1/DS1_FromTEM2/'

# Find all the reconstruct files in a dataset
TEM2_filenames=natsorted([i for i in os.listdir(TEM2_path) if (('FULL-REG' in i) and ('.ser' not in i)) ])
R_filenames=natsorted([i for i in os.listdir(R_path) if (('FULL-REG' in i) and ('.ser' not in i)) ])

nfiles=len(TEM2_filenames)

for z in range(nfiles):

    TEM2_name=TEM2_filenames[z]
    R_name=R_filenames[z]
    print TEM2_name,'---->',R_name
    # Read the Reconstruct trace file for slice z and get the root of the xlm tree
    R_tree = ET.parse(R_path+R_name)
    R_root = R_tree.getroot()

    # Get the name of the image
    R_img_tag = R_root.find('./Transform/Image')
    R_img_name = R_img_tag.attrib['src']

    R_domain=R_root.findall("./Transform/Contour[@name='domain1']")[0]
    R_domain_points = R_domain.attrib['points']

    # Read the TEM2 exported file for slice z and get the root of the xlm tree
    TEM2_tree = ET.parse(TEM2_path+TEM2_name)
    TEM2_root = TEM2_tree.getroot()

    TEM2_img_tag = TEM2_root.find('./Transform/Image')
    TEM2_img_tag.set('src',R_img_name) #Replace the image name by the one in the reconstruct file
    TEM2_img_tag.set('brightness','0') #Replace the brightness by the one in the reconstruct file
    TEM2_img_tag.set('contrast','1') #Replace the contrast by the one in the reconstruct file

    TEM2_section_index =  TEM2_root.attrib['index']
    TEM2_root.set('index',str(int(TEM2_section_index)+1))

    N_tree = deepcopy(TEM2_tree)
    N_root = N_tree.getroot()

    for transform in N_root.findall('./Transform'):
        N_root.remove(transform)

    N_root.append(TEM2_img_tag.getparent())
    for contour in TEM2_root.findall('./Transform/Contour'):
        transform = contour.getparent()
        new_transform=ET.SubElement(N_root,'Transform',attrib=transform.attrib)
        new_contour=ET.SubElement(new_transform,'Contour',attrib=contour.attrib)
        lpoints=contour.attrib['points'].split(',')[0:-1] #take up to the second to last element to eliminate remining spaces
        lpoints = np.array([map(float,i.split()) for i in lpoints])
        ct_output=''
        for i in range(lpoints.shape[0]):
            ct_output=ct_output+'%1.4f' % lpoints[i,0]+' '+'%1.4f' % lpoints[i,1]+','+'\n'
            new_contour.set('points',ct_output)

    filestring = ET.tostring(N_tree,pretty_print=True,
                      xml_declaration=True, encoding=TEM2_tree.docinfo.encoding,
                      doctype=TEM2_tree.docinfo.doctype)

    filestring = filestring.replace('><','> \n <').replace('&#10;','\n \t')
    f=open(outpath+R_name,'w')
    f.write(filestring)
    f.close()
