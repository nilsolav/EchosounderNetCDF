# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 11:57:07 2019

@author: sindrev


WORK IN PROGRESS

TODO: 
    1. Finish the erased regions
    2. Test function on data
    3. Write to netcdf structure

"""




def readLSSSwork(nc_raw_file,work_file):
    '''
    Description: 
        This script reads the LSSS work information and convert it to polygons
        that can be used to store the interpretation mask as a seperate netcdf 
        file
        
    Author: 
        Sindre Vatnehol
        Institute of Marine Research
        Bergen, Norway
        
    '''
    
    
    
    #Load som packages
    from netCDF4 import Dataset
    import numpy as np
    import xmltodict
    from datetime import datetime
    from scipy.spatial import Delaunay
    
    
    
    
    
    #Some other functions
    def alpha_shape(points, alpha, only_outer=True):
        """
        Compute the alpha shape (concave hull) of a set of points.
        :param points: np.array of shape (n,2) points.
        :param alpha: alpha value.
        :param only_outer: boolean value to specify if we keep only the outer border
        or also inner edges.
        :return: set of (i,j) pairs representing edges of the alpha-shape. (i,j) are
        the indices in the points array.
        """
        assert points.shape[0] > 3, "Need at least four points"
    
        def add_edge(edges, i, j):
            """
            Add an edge between the i-th and j-th points,
            if not in the list already
            """
            if (i, j) in edges or (j, i) in edges:
                # already added
                assert (j, i) in edges, "Can't go twice over same directed edge right?"
                if only_outer:
                    # if both neighboring triangles are in shape, it's not a boundary edge
                    edges.remove((j, i))
                return
            edges.add((i, j))
    
        tri = Delaunay(points)
        edges = set()
        # Loop over triangles:
        # ia, ib, ic = indices of corner points of the triangle
        for ia, ib, ic in tri.vertices:
            pa = points[ia]
            pb = points[ib]
            pc = points[ic]
            # Computing radius of triangle circumcircle
            # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
            a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
            b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
            c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
            s = (a + b + c) / 2.0
            area = np.sqrt(s * (s - a) * (s - b) * (s - c))
            circum_r = a * b * c / (4.0 * area)
            if circum_r < alpha:
                add_edge(edges, ia, ib)
                add_edge(edges, ib, ic)
                add_edge(edges, ic, ia)
        return edges
    
        
    def find_edges_with(i, edge_set):
        i_first = [j for (x,j) in edge_set if x==i]
        i_second = [j for (j,x) in edge_set if x==i]
        return i_first,i_second
    
    def stitch_boundaries(edges):
        edge_set = edges.copy()
        boundary_lst = []
        while len(edge_set) > 0:
            boundary = []
            edge0 = edge_set.pop()
            boundary.append(edge0)
            last_edge = edge0
            while len(edge_set) > 0:
                i,j = last_edge
                j_first, j_second = find_edges_with(j, edge_set)
                if j_first:
                    edge_set.remove((j, j_first[0]))
                    edge_with_j = (j, j_first[0])
                    boundary.append(edge_with_j)
                    last_edge = edge_with_j
                elif j_second:
                    edge_set.remove((j_second[0], j))
                    edge_with_j = (j, j_second[0])  # flip edge rep
                    boundary.append(edge_with_j)
                    last_edge = edge_with_j
    
                if edge0[0] == last_edge[1]:
                    break
    
            boundary_lst.append(boundary)
        return boundary_lst

    
    
    
    
    #Define a structure type
    class structtype(): 
        pass
    
    
    
    
    #For bookkeeping
    Out = structtype()
    Out.ExcludedRegion = dict()
    Out.SchoolRegion = dict()
    Out.LayerRegion = dict()
    
    
    
    #Get the information from the netcdf file of the acoustic data
    #What we need is information of the time and the depth 
    inn = Dataset(nc_raw_file)
    
    
    
    #Get the reference time and depth from acoustic file
    ref_time = inn.groups['Sonar'].groups['Beam_group2'].variables['ping_time'][:].data/1e9
    dept_ref = np.arange(len(inn.groups['Sonar'].groups['Beam_group2'].variables['Power'][0]))*(inn.groups['Sonar'].groups['Beam_group2'].variables['sample_interval'][0]*1490/2)
    sample_space = (inn.groups['Sonar'].groups['Beam_group2'].variables['sample_interval'][0]*1490/2)







    #Parse the LSSS work file
    with open(work_file) as fd:
        doc = xmltodict.parse(fd.read())
    
    
    
    
    
    
    
    
    
    '''
    #Start processing the excluded regions
    #For this, only the start time and number of pings are given. 
    '''
    i=0
    
    
    
    #If there is some excluded region
    if not doc['regionInterpretation']['exclusionRanges']== None: 
        
     
        if type(doc['regionInterpretation']['exclusionRanges']['timeRange']) == list: 
            for exclude_region in doc['regionInterpretation']['exclusionRanges']['timeRange'] :
                Out.ExcludedRegion[i]=structtype()
                    
                
                #Get start time and number of pings of the excluded region
                time_start = float(exclude_region['@start'])
                numberOfPings = int(exclude_region['@numberOfPings'])
                
                
                #Convert to number of seconds since 1601
                #This is because different time convention are used between lsss and the ices format
                second_since_epoch = (datetime.utcfromtimestamp(time_start)-datetime(1601, 1, 1, 0, 0, 0) ).total_seconds()
                
                
                #Grab ping number and depth 
                Ping = np.arange(np.argmin(abs(second_since_epoch-ref_time)),(np.argmin(abs(second_since_epoch-ref_time))+numberOfPings))
                Depth = np.hstack((np.zeros((len(Ping))),np.nanmax(dept_ref)*np.ones((len(Ping)))))
                
                
                #Do some modification, so it will look nice when plotting
                Ping = np.hstack((Ping,np.flip(Ping)))
                ping_=[]
                for ping in Ping: 
                    ping_ = np.hstack((ping_,ref_time[ping]))
                
                Out.ExcludedRegion[i].Time = ping_
                Out.ExcludedRegion[i].Ping = Ping
                Out.ExcludedRegion[i].Depth = Depth
                
            
                i+=1
                
        #If there is ony one region
        else: 
            exclude_region=doc['regionInterpretation']['exclusionRanges']['timeRange']
            
            Out.ExcludedRegion[i]=structtype()
                    
            print(exclude_region['@start'])
            #Get start time and number of pings of the excluded region
            time_start = float(exclude_region['@start'])
            numberOfPings = int(exclude_region['@numberOfPings'])
            
            
            #Convert to number of seconds since 1601
            #This is because different time convention are used between lsss and the ices format
            second_since_epoch = (datetime.utcfromtimestamp(time_start)-datetime(1601, 1, 1, 0, 0, 0) ).total_seconds()
            
            
            #Grab ping number and depth 
            Ping = np.arange(np.argmin(abs(second_since_epoch-ref_time)),(np.argmin(abs(second_since_epoch-ref_time))+numberOfPings))
            Depth = np.hstack((np.zeros((len(Ping))),np.nanmax(dept_ref)*np.ones((len(Ping)))))
            
            
            #Do some modification, so it will look nice when plotting
            Ping = np.hstack((Ping,np.flip(Ping)))
            ping_=[]
            for ping in Ping: 
                ping_ = np.hstack((ping_,ref_time[ping]))
            
            Out.ExcludedRegion[i].Time = ping_
            Out.ExcludedRegion[i].Ping = Ping
            Out.ExcludedRegion[i].Depth = Depth
        
        
        
        
        
        
        
        
        
        
    '''
    #Start processing school information
    '''
    
    #Test if there is schools in the work file
    if not doc['regionInterpretation']['schoolInterpretation']== None: 

        
        
        ik=0
        if type(doc['regionInterpretation']['schoolInterpretation']['schoolMaskRep']) == list: 
            for layer in doc['regionInterpretation']['schoolInterpretation']['schoolMaskRep']: 
                
                
                #For bookkeeping
                Out.SchoolRegion[ik]=structtype()
                Depth = []
                Ping = []
                Interpretation = dict() 
                
                
                #Get species fractin
                if (type(layer['speciesInterpretationRoot']['speciesInterpretationRep']))==list: 
                    iik=0
                    Interpretation=structtype()
                    Interpretation.Freq = []
                    Interpretation.Spec = []
                    Interpretation.Fraction = []
                    for freq in layer['speciesInterpretationRoot']['speciesInterpretationRep']: 
                        if (type(freq['species'])) == list: 
                            for spec in freq['species']:
                                Interpretation.Freq = np.hstack((Interpretation.Freq,freq['@frequency']))
                                Interpretation.Spec = np.hstack((Interpretation.Spec,spec['@ID']))
                                Interpretation.Fraction = np.hstack((Interpretation.Fraction,spec['@fraction']))
                        else: 
                            spec = freq['species']
                            Interpretation.Freq = np.hstack((Interpretation.Freq,freq['@frequency']))
                            Interpretation.Spec = np.hstack((Interpretation.Spec,spec['@ID']))
                            Interpretation.Fraction = np.hstack((Interpretation.Fraction,spec['@fraction']))
                            
                    iik+=1
                else: 
                    print(layer['speciesInterpretationRoot']['speciesInterpretationRep'])
                    asdf
            
                
                
                for ping in layer['pingMask']: 
                    depth_i = np.array(ping['#text'].split(' '),float)
                    
                    for i in np.arange(0,np.int(len(depth_i)/2)):
                        depth_ii = np.arange(depth_i[0+2*i],depth_i[1+2*i],sample_space)
                        Depth = np.hstack((Depth,depth_ii))
                        Ping = np.hstack((Ping,np.repeat(int(ping['@relativePingNumber']),len(depth_ii))))                  
        
                points = np.vstack([Ping, Depth]).T
                
                edges = alpha_shape(points, alpha=1, only_outer=True)
                edges = stitch_boundaries(edges)
                Ping = []
                Depth = []
          
                for i, j in edges[0]:
                    Ping = np.hstack((Ping,points[[i, j], 0][0]))
                    Depth = np.hstack((Depth,points[[i, j], 1][0]))
                    
                ping_=[]
                for ping in Ping: 
                    ping_ = np.hstack((ping_,ref_time[np.int(ping)]))
                
                Out.SchoolRegion[ik].Time = ping_
                Out.SchoolRegion[ik].Ping = Ping
                Out.SchoolRegion[ik].Depth = Depth
                Out.SchoolRegion[ik].Interpretation = Interpretation
                ik+=1
                
                        
        else: 
            layer = doc['regionInterpretation']['schoolInterpretation']['schoolMaskRep']
            Out.SchoolRegion[ik]=structtype()
            Depth = []
            Ping = []
            Interpretation = structtype()
                
            if (type(layer['speciesInterpretationRoot']['speciesInterpretationRep']))==list: 
                iik=0
                Interpretation.Freq = []
                Interpretation.Spec = []
                Interpretation.Fraction = []
                for freq in layer['speciesInterpretationRoot']['speciesInterpretationRep']: 
                    if (type(freq['species'])) == list: 
                        for spec in freq['species']:
                            
                            Interpretation.Freq = np.hstack((Interpretation.Freq,freq['@frequency']))
                            Interpretation.Spec = np.hstack((Interpretation.Spec,spec['@ID']))
                            Interpretation.Fraction = np.hstack((Interpretation.Fraction,spec['@fraction']))
                    else: 
                        spec = freq['species']
                        Interpretation.Freq = np.hstack((Interpretation.Freq,freq['@frequency']))
                        Interpretation.Spec = np.hstack((Interpretation.Spec,spec['@ID']))
                        Interpretation.Fraction = np.hstack((Interpretation.Fraction,spec['@fraction']))
                        
                iik+=1
            else: 
                print(layer['speciesInterpretationRoot']['speciesInterpretationRep'])
                asdf
                    
                    
                    
            for ping in layer['pingMask']: 
                depth_i = np.array(ping['#text'].split(' '),float)
                
                for i in np.arange(0,np.int(len(depth_i)/2)):
                    depth_ii = np.arange(depth_i[0+2*i],depth_i[1+2*i],sample_space)
                    Depth = np.hstack((Depth,depth_ii))
                    Ping = np.hstack((Ping,np.repeat(int(ping['@relativePingNumber']),len(depth_ii))))                  
    
            points = np.vstack([Ping, Depth]).T
            
            edges = alpha_shape(points, alpha=1, only_outer=True)
            edges = stitch_boundaries(edges)
            Ping = []
            Depth = []
      
            for i, j in edges[0]:
                Ping = np.hstack((Ping,points[[i, j], 0][0]))
                Depth = np.hstack((Depth,points[[i, j], 1][0]))
                
            ping_=[]
            for ping in Ping: 
                ping_ = np.hstack((ping_,ref_time[np.int(ping)]))
            
            Out.SchoolRegion[ik].Time = ping_
            Out.SchoolRegion[ik].Ping = Ping
            Out.SchoolRegion[ik].Depth = Depth
            Out.SchoolRegion[ik].Interpretation = Interpretation
        
        
        
        
        
        
    '''
    #Start processing layer information
    
    TODO: 
        Do a proper bugtesting o fthis one
    '''
    #Find catthegory
    i=0
    
#    num_layers = len(doc['regionInterpretation']['layerInterpretation']['layerDefinitions'])
    print(type(doc['regionInterpretation']['layerInterpretation']['layerDefinitions']['layer']))
    if type(doc['regionInterpretation']['layerInterpretation']['layerDefinitions']['layer'])==list:
        for layer in doc['regionInterpretation']['layerInterpretation']['layerDefinitions']['layer']: 
            Out.LayerRegion[i] = structtype()
            Depth = []
            Ping = []
            for bound in ((layer['boundaries'])): 
                for bound_i in (layer['boundaries'][bound]): 
                    
                    for bound_ in doc['regionInterpretation']['layerInterpretation']['boundaries'][bound]: 
                        if bound_['@id']==bound_i['@id']: 
                            try: 
                                depth = np.array(bound_['curveRep']['depths'].replace('\n',' ').split(' '),float)
                                ping = np.arange(np.int(bound_['curveRep']['pingRange']['@startOffset']),np.int(bound_['curveRep']['pingRange']['@startOffset'])+np.int(bound_['curveRep']['pingRange']['@numberOfPings']))
                                if bound_i['@isUpper'] == 'true':
                                    Depth=np.hstack((Depth,np.flip(depth)))
                                    Ping = np.hstack((Ping,np.flip(ping)))
                                else: 
                                    Depth=np.hstack((Depth,depth))
                                    Ping = np.hstack((Ping,ping))
                                    
                            except: 
                                d=1
                
            Ping = np.hstack((Ping,Ping[0]))
            Depth = np.hstack((Depth,Depth[0]))
            
            
            ping_=[]
            for ping in Ping: 
                ping_ = np.hstack((ping_,ref_time[np.int(ping)]))
            
            Out.LayerRegion[i].Time = ping_
            Out.LayerRegion[i].Ping = Ping
            Out.LayerRegion[i].Depth = Depth
            i=i+1
    else: 
        
        Out.LayerRegion[i] = structtype()
        Depth = []
        Ping = []
        layer = doc['regionInterpretation']['layerInterpretation']['layerDefinitions']['layer']
#        print(layer)
##        
        bound = doc['regionInterpretation']['layerInterpretation']['layerDefinitions']['layer']['boundaries']['curveBoundary']
        
        for bound_i in bound: 
            for bound_ in doc['regionInterpretation']['layerInterpretation']['boundaries']['curveBoundary']: 
                if bound_['@id']==bound_i['@id']: 
                    depth = np.array(bound_['curveRep']['depths'].replace('\n',' ').split(' '),float)
                    ping = np.arange(np.int(bound_['curveRep']['pingRange']['@startOffset']),np.int(bound_['curveRep']['pingRange']['@startOffset'])+np.int(bound_['curveRep']['pingRange']['@numberOfPings']))
                    if bound_i['@isUpper'] == 'true':
                        Depth=np.hstack((Depth,np.flip(depth)))
                        Ping = np.hstack((Ping,np.flip(ping)))
                    else: 
                        Depth=np.hstack((Depth,depth))
                        Ping = np.hstack((Ping,ping))
                        
            
        Ping = np.hstack((Ping,Ping[0]))
        Depth = np.hstack((Depth,Depth[0]))
        
        print(Depth)
        ping_=[]
        for ping in Ping: 
            ping_ = np.hstack((ping_,ref_time[np.int(ping)]))
        
        Out.LayerRegion[i].Time = ping_
        Out.LayerRegion[i].Ping = Ping
        Out.LayerRegion[i].Depth = Depth
#        
#        
        
        
        
        
        
    '''
    #Start processing erased region information
    
    This is a TODO
    '''
        
    try: 
        doc['regionInterpretation']['masking']['mask']
        run=True
    except: 
        run=False
            
    if run: 
        for mask in doc['regionInterpretation']['masking']['mask']: 
            if mask['@channelID'] == '2': 
                Ping= []
                Depth = []
                
                for ping in mask['ping']: 
                    tick = 0
                    for depth in np.asarray(ping['#text'].split(' '),float): 
                        if depth <0: 
                            depth = 0
                        elif depth >500: 
                            depth = 500
                            
                        Ping = np.hstack((Ping,int(ping['@pingOffset'])))
                        if tick == 0: 
                            Depth = np.hstack((Depth,depth))
                            tick +=1
                        else: 
                            Depth = np.hstack((Depth,depth+Depth[-1]))
                            
                            
        
    return(Out)
        
    
    