
def readLSSSwork(nc_raw_file,work_file):
    '''
    Description: 
        This script reads the LSSS work information and convert it to polygons
        that can be used to store the interpretation mask as a seperate netcdf 
        file
    '''
    
    
    
    #Load som packages
    from netCDF4 import Dataset
    import numpy as np
    import xmltodict
    from datetime import datetime
    from scipy.spatial import Delaunay
    
    
    
    
    
    #Define a structure type
    class structtype(): 
        pass
    
    
    
    
    #Define the output structure
    Out = structtype()
    Out.ExcludedRegion = dict()
    Out.SchoolRegion = dict()
    Out.LayerRegion = dict()
    Out.ErasedRegion = dict()
    
    
    

    #Get the information from the netcdf file of the acoustic data
    #What we need is information of the time and the depth 
    inn = Dataset(nc_raw_file)
    
    
    
    #Get the reference time and depth from acoustic file
    #TODO: 
    #    investigate in how to make this function independent of the raw data
    ref_time = inn.groups['Sonar'].groups['Beam_group2'].variables['ping_time'][:].data/1e9
    dept_ref = np.arange(len(inn.groups['Sonar'].groups['Beam_group2'].variables['Power'][0]))*(inn.groups['Sonar'].groups['Beam_group2'].variables['sample_interval'][0]*1490/2)
#    sample_space = (inn.groups['Sonar'].groups['Beam_group2'].variables['sample_interval'][0]*1490/2)







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
#                Depth = np.hstack((np.zeros((len(Ping))),np.nanmax(dept_ref)*np.ones((len(Ping)))))
                
                
#                print(len(Depth))
#                print(len(Ping))
                
                #Do some modification, so it will look nice when plotting
#                Ping = np.hstack((Ping,np.flip(Ping)))
                ping_=[]
                Depth = []
                for ping in Ping: 
                    ping_ = np.hstack((ping_,ref_time[ping]))
                    Depth.append(np.array([-1000000,1000000]))

                Depth = np.array(Depth)

                #Write to the output structure
                Out.ExcludedRegion[i].Time = ping_
                Out.ExcludedRegion[i].Ping = Ping
                Out.ExcludedRegion[i].Depth = Depth
                
            
                i+=1
      



          
        #If there is ony one region
        else: 
            exclude_region=doc['regionInterpretation']['exclusionRanges']['timeRange']
            
            Out.ExcludedRegion[i]=structtype()
                    
            #Get start time and number of pings of the excluded region
            time_start = float(exclude_region['@start'])
            numberOfPings = int(exclude_region['@numberOfPings'])
            
            
            #Convert to number of seconds since 1601
            #This is because different time convention are used between lsss and the ices format
            second_since_epoch = (datetime.utcfromtimestamp(time_start)-datetime(1601, 1, 1, 0, 0, 0) ).total_seconds()
            
            
            #Grab ping number and depth 
            Ping = np.arange(np.argmin(abs(second_since_epoch-ref_time)),(np.argmin(abs(second_since_epoch-ref_time))+numberOfPings))
#            Depth = np.hstack((np.zeros((len(Ping))),np.nanmax(dept_ref)*np.ones((len(Ping)))))
            
            
            
            #Do some modification, so it will look nice when plotting
#            Ping = np.hstack((Ping,np.flip(Ping)))
            ping_=[]
            Depth = []
            for ping in Ping: 
                ping_ = np.hstack((ping_,ref_time[ping]))
                Depth.append(np.array([-1000000,1000000]))

            Depth = np.array(Depth)
            
            Out.ExcludedRegion[i].Time = ping_
            Out.ExcludedRegion[i].Ping = Ping
            Out.ExcludedRegion[i].Depth = Depth
        
        
        
        
        
        
        
        
    '''
    #Start processing layer information
    
    TODO: 
        Do a proper bugtesting o fthis one
    '''
    #Find catthegory
    i=0
    
#    num_layers = len(doc['regionInterpretation']['layerInterpretation']['layerDefinitions'])
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
        
        ping_=[]
        for ping in Ping: 
            ping_ = np.hstack((ping_,ref_time[np.int(ping)]))
        
        Out.LayerRegion[i].Time = ping_
        Out.LayerRegion[i].Ping = Ping
        Out.LayerRegion[i].Depth = Depth
#        
#        
        
        
        
        
        
    '''
    Start processing erased region information
    
    
    This saves mask information (start-stop)
    
    
    '''
        
    
    #Check if information is within the work file
    try: 
        doc['regionInterpretation']['masking']['mask']
        run=True
    except: 
        run=False
            
        
    #If it is, continiue
    #This works for both one and t
    if run: 
        
        
        i = 0
        #Run throguh each channel
        for mask in doc['regionInterpretation']['masking']['mask']: 
            
            #Define a structure
            Out.ErasedRegion[i]=structtype()
            
            
            #For bookkeeping
            Ping= []
            Depth = []
            
            
            
            #Store data
            for ping in mask['ping']: 
                Ping = np.hstack((Ping,int(ping['@pingOffset'])))
                depth = np.cumsum(np.asarray(ping['#text'].split(' '),float))
                Depth.append(depth) # 
                
            
            #Grab time
            Depth = np.asarray(Depth)
            ping_=[]
            for ping in Ping: 
                ping_ = np.hstack((ping_,ref_time[np.int(ping)]))    
            
            
            #Store data into structure
            Out.ErasedRegion[i].Time = ping_
            Out.ErasedRegion[i].Ping = Ping
            Out.ErasedRegion[i].Depth = Depth 
            Out.ErasedRegion[i].Channel = mask['@channelID'] 

            i=i+1

        
                            
            



    #School mask
    
    layer = doc['regionInterpretation']['schoolInterpretation']['schoolMaskRep']

    try: 
        layer[0]
        run=True
    except KeyError:
        run = False
          
    
    if run == False: 
        ii = 0
        Out.SchoolRegion[ii]=structtype()
        Ping= []
        Depth = []
    
        
        for ping in layer['pingMask']: 
            Ping = np.hstack((Ping,int(ping['@relativePingNumber'])))
            depth = (np.asarray(ping['#text'].split(' '),float))
            Depth.append(depth) # 
        ping_=[]
        for ping in Ping: 
            ping_ = np.hstack((ping_,ref_time[np.int(ping)]))    
            
        Out.SchoolRegion[ii].Time = ping_
        Out.SchoolRegion[ii].Ping = Ping
        Out.SchoolRegion[ii].Depth = Depth 
            
    else: 
        
        for ii in range(len(layer)): 
            Out.SchoolRegion[ii]=structtype()
            
            Ping= []
            Depth = []
        
            for ping in layer[ii]['pingMask']: 
                Ping = np.hstack((Ping,int(ping['@relativePingNumber'])))
                depth = (np.asarray(ping['#text'].split(' '),float))
                Depth.append(depth) # 
            ping_=[]
            for ping in Ping: 
                ping_ = np.hstack((ping_,ref_time[np.int(ping)]))    
                
        
    
            Out.SchoolRegion[ii].Time = ping_
            Out.SchoolRegion[ii].Ping = Ping
            Out.SchoolRegion[ii].Depth = Depth 
    return(Out)
        
    
    
