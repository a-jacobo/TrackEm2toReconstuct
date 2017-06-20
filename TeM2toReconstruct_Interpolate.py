# -*- coding: utf-8 -*-

from lxml import etree as ET
import numpy as np
import os
from copy import deepcopy

from natsort import natsorted
from scipy import interpolate
from scipy.stats import mode

''' 
 Adrian Jacobo - May 2017
 Fix image names in TrackEm2 exported files using the original Reconstruct file.
 Register traces back to images.
 Remove jaggedness from TrackEm2 traces using interpolation.
'''

# Path to the TrackEM2 traces
#TEM2_path ='/Volumes/Eliot_Datasets/TrakEM-Reconstruct Tests/Stair_Case/'
# Path to the Reconstruct traces
#R_path ='/Volumes/Eliot_Datasets/TrakEM-Reconstruct Tests/Stair_Case/'
# Output path
#outpath ='/Volumes/Eliot_Datasets/TrakEM-Reconstruct Tests/Stair_Case/Smoothed/'

dataset_name='T5.9'


# Path to the TrackEM2 traces
TEM2_path ='/Volumes/Eliot_Datasets/Original_Datasets/'+dataset_name+'/20170608_T5.9_TrakEM-export/'

# Path to the Reconstruct traces
R_path ='/Volumes/Eliot_Datasets/Original_Datasets/'+dataset_name+'/'

# Output path
outpath ='/Volumes/Eliot_Datasets/Original_Traces/'+dataset_name+'/'+dataset_name+'_FromTEM2/'

# Find all the reconstruct files in a dataset

dataset_name='T5-9'
TEM2_filenames=natsorted([i for i in os.listdir(TEM2_path) if ((dataset_name in i) and ('.ser' not in i) and ('tif' not in i)  ) ])
R_filenames=natsorted([i for i in os.listdir(R_path) if ((dataset_name in i) and ('.ser' not in i)
                                                         and ('.xml' not in i) and ('TrakEM' not in i) and ('tif' not in i)
                                                         and ('DM4'not in i)) ])

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

    N_root.append(deepcopy(TEM2_img_tag.getparent()))

    shifts=np.empty((0,2),dtype=np.float)
    #Fast way to compute shifts between traces of the TEM2 exported and original Reconstruct traces
    for contour in R_root.findall('./Transform/Contour'):
        name = contour.attrib['name']
        branch = TEM2_root.findall("./Transform/Contour[@name='"+name+"']")
        if len(branch)!=0:
            contour2 = branch[0]
            lpoints=contour.attrib['points'].split(',')[0:-1] #take up to the second to last element to eliminate remining spaces
            lpoints = np.array([map(float,i.split()) for i in lpoints])
            lpoints2=contour2.attrib['points'].split(',')[0:-1] #take up to the second to last element to eliminate remining spaces
            lpoints2 = np.array([map(float,i.split()) for i in lpoints2])
            shifts=np.vstack((shifts,np.mean(lpoints,axis=0)-np.mean(lpoints2,axis=0)))

    shift = mode(np.around(shifts,decimals=4),axis=0)[0][0]
    print 'Shifting dataset by ',shift

    for contour in TEM2_root.findall('./Transform/Contour'):
        name = contour.attrib['name']
        ct_output=''
        transform = contour.getparent()
        new_transform=ET.SubElement(N_root,'Transform',attrib=transform.attrib)
        new_contour=ET.SubElement(new_transform,'Contour',attrib=contour.attrib)
        if name !='domain1':
            lpoints=contour.attrib['points'].split(',')[0:-1] #take up to the second to last element to eliminate remining spaces
            lpoints = np.array([map(float,i.split()) for i in lpoints])
           
            # Find distances between points, if they are too small interpolate to smooth the curves
            #print contour.attrib['name']
            diff = np.diff(lpoints,axis=0)
            dist = np.linalg.norm(diff,axis=1)
            # Interpolate curves and replace them into lpoints
            if ((0 < np.min(dist)<0.005) and (lpoints.shape[0] > 10)):
                # append the starting x,y coordinates
                x=lpoints[:,0]
                y=lpoints[:,1]

                x = np.r_[x, x[0]]
                y = np.r_[y, y[0]]
                try:
                    # fit splines to x=f(u) and y=g(u), treating both as periodic. also note that s=0
                    # is needed in order to force the spline fit to pass through all the input points.
                    tck, u = interpolate.splprep([x, y], s=0, per=True)
                    n=len(x)/2
                    #print ' interpolated:',len(x),'-->',n
                    # evaluate the spline fits for 1000 evenly spaced distance values
                    xi, yi = interpolate.splev(np.linspace(0, 1, n), tck)
                    lpoints=np.stack((xi,yi),axis=1)
                except:
                    print contour.attrib['name'],len(x)
                    print 'Unable to interpolate'
                    print x
                    print i
        
            for i in range(lpoints.shape[0]):
                lpoints[i,0]=lpoints[i,0]+shift[0] #Shift new dataset
                lpoints[i,1]=lpoints[i,1]+shift[1]
                ct_output=ct_output+'%1.4f' % lpoints[i,0]+' '+'%1.4f' % lpoints[i,1]+','+'\n'
            new_contour.set('points',ct_output)
        else:
        # For domain1 get the boundary from the original file.
            R_domain=R_root.findall("./Transform/Contour[@name='domain1']")[0]
            new_contour.set('points',R_domain.attrib['points'])

    filestring = ET.tostring(N_tree,pretty_print=True,
                      xml_declaration=True, encoding=TEM2_tree.docinfo.encoding,
                      doctype=TEM2_tree.docinfo.doctype)

    filestring = filestring.replace('><','> \n <').replace('&#10;','\n \t')
    f=open(outpath+R_name,'w')
    f.write(filestring)
    f.close()
#    exit()


#    filestring = ET.tostring(TEM2_tree,pretty_print=True,
#                         xml_declaration=True, encoding=TEM2_tree.docinfo.encoding,
#                         doctype=TEM2_tree.docinfo.doctype)


#    for contour in TEM2_root.findall('./Transform/Contour'):
#        if contour.attrib['name']== 'domain1':
#                contour.set('points',R_domain_points)
#        else:
#                contour.set('mode','11')

#    filestring = ET.tostring(TEM2_tree,pretty_print=True,
#                xml_declaration=True, encoding=TEM2_tree.docinfo.encoding,
#                doctype=TEM2_tree.docinfo.doctype)

#    TEM2_tree.write(outpath+R_name,pretty_print=True,
#                    xml_declaration=True, encoding=TEM2_tree.docinfo.encoding,
#                    doctype=TEM2_tree.docinfo.doctype)

#    f=open(outpath+R_name,'w')
#    f.write(filestring)
#    f.close()
#        content = f.read()
#        f.seek(0, 0)
#        f.write(line.rstrip('\r\n') + '\n' + content)

