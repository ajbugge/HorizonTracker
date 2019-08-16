
'''Author: Aina Juell Bugge. aina.juell.bugge@gmail.com'''
'''Code for the paper: "Automatic extraction of dislocated horizons from 3D seismic data using non-local trace matching", Aina Juell Bugge, Jan Erik Lie, Andreas Kjelsrud Evensen, Jan Inge Faleide, and Stuart Clark, Geophysics, 2019.'''

import scipy.io
import numpy as np
import time
from tslearn import metrics
from scipy.signal import argrelextrema

def CreateGrid (seismic_data, grid_step):
    xls=[]
    ils=[]
    for i in range (0, seismic_data.shape[1], grid_step):
        xls.append(i)
    for i in range (0, seismic_data.shape[2], grid_step):
        ils.append(i)
    grid=[]
    for ii in xls:
        for jj in ils:
            seismictrace=seismic_data[:, ii, jj]
            nnz_total=np.count_nonzero(seismictrace)
            
            if seismictrace.any() and nnz_total>0:
                trace=[ii, jj]
                grid.append(trace)   
    print('Number of traces in grid:', len(grid))
    values_to_sort_from=[]
    for a in grid:
        z=a[0]+a[1]
        values_to_sort_from.append(z)
    grid = [x for _,x in sorted(zip(values_to_sort_from,grid))]
    return grid


#grid_step2=5
#horizontal_warp_distance =100

import sys
import time
from IPython.display import clear_output

def DynamicTimeWarping (seismic_data, grid, window_width):
    path_info=[]
    start = time.time()
    
    for i in range (0, len(grid), 1):  
        reference_trace=seismic_data[:, grid[i][0], grid[i][1]]
        reference_trace=np.trim_zeros(reference_trace, 'b')
        block_size=0

        for ii in range (0, len(grid), 1): 
            if abs(grid[i][0]-grid[ii][0])<=window_width and abs(grid[i][1]-grid[ii][1])<=window_width:
                block_size=block_size+1
                matched_trace=seismic_data[:, grid[ii][0], grid[ii][1]]
                matched_trace=np.trim_zeros(matched_trace, 'b')
                path, sim= metrics.dtw_path(reference_trace, matched_trace)   #print(grid[i][0], grid[i][1], grid[ii][0], grid[ii][1], sim)
                
                #if sim < similarity_score:
                PI=[grid[i][0], grid[i][1], grid[ii][0], grid[ii][1], path]
                path_info.append(PI)
        stop = time.time()
        t=(stop-start)/60
        
        clear_output()
        print('reference trace:', grid[i][0], grid[i][1],'(',i,'out of:',len(grid),')',' time in minutes: ', t, 'Sub-block size:', block_size)
        
        
    return path_info


# events can be peaks, troughs or peaks and troughs


import time
from scipy.signal import argrelextrema
def ExtractHorizon (seismic_data, grid, path_info, events, event_spacing_parameter):
    temporary_horizons=[]
    start = time.time()
    it=0
    
    start_horizon_tracker=[] #can be a list of several starting points. (e.g. [0, 0], [20, 20], [40, 20] if grid step=20)
    starting_point=[0, 0]  
    start_horizon_tracker.append(starting_point)

    for start_trace in start_horizon_tracker:
        reference_trace=seismic_data[:, start_trace[0], start_trace[1]]
        reference_trace=np.trim_zeros(reference_trace, 'b')
        
        if events == 'peaks':
            amplitudes= argrelextrema(reference_trace, np.greater, order=event_spacing_parameter)
            amplitudes=amplitudes[0]
        if events == 'troughs':
            amplitudes= argrelextrema(reference_trace, np.less, order=event_spacing_parameter)
            amplitudes=amplitudes[0]
            
        if events == 'both':
            peaks=argrelextrema(reference_trace, np.greater, order=event_spacing_parameter)
            troughs=argrelextrema(reference_trace, np.less, order=event_spacing_parameter)
            amplitudes=np.insert(peaks[0], [0], troughs[0])
            amplitudes=np.sort(amplitudes)
            
        #print("number of horizons to track:", len(amplitudes), "starting trace:", start_trace)

        for event in amplitudes:
            reference_trace=seismic_data[:, start_trace[0], start_trace[1]]
            reference_trace=np.trim_zeros(reference_trace, 'b')       
            reflector=[]

            for k in range (0, len(path_info)):
                t1 = path_info[k][0]
                t2 = path_info[k][1]
                if t1 == start_trace[0] and t2 == start_trace[1]:
                    matched_trace=seismic_data[:, t1, t2]
                    matched_trace=np.trim_zeros(matched_trace)
                    path= path_info[k][4]
                    for ll in range (0, len(path)):
                        if path[ll][0]==event:
                            point=[path[ll][1], path_info[k][2], path_info[k][3]]
                            reflector.append(point)

            for i in range (0, len(grid)):  
                reference_trace=seismic_data[:, grid[i][0], grid[i][1]]
                tempevent=[]
                for k in reflector:
                    t1 =grid[i][0]
                    t2 = grid[i][1]
                    if t1 ==k[1] and t2 ==k[2]:
                        tempevent.append(k[0])
                if tempevent:
                    counts = np.bincount(tempevent)
                    event=np.argmax(counts)
                    for k in range (0, len(path_info)):
                        t1 = path_info[k][0]
                        t2 = path_info[k][1]
                        if t1 == grid[i][0]and t2 == grid[i][1]:
                            matched_trace=seismic_data[:, t1, t2]
                            path= path_info[k][4]
                            for ll in range (0, len(path)):
                                if path[ll][0]==event:
                                    point=[path[ll][1], path_info[k][2], path_info[k][3]]
                                    reflector.append(point)
            temporary_horizons.append(reflector)

            stop = time.time()
            t=(stop-start)/60
            clear_output()
            print("status:", 'horizon number', it, 'of', len(amplitudes))
            print('time in minutes: ', t)
            it=it+1
            
    return temporary_horizons




def HorizonAccuracyFilter (seismic_data, input_horizons, grid, path_info,print_output, hitrate_per_trace_percent=0, hardcoded_minimum_hitrate=0):
    # hitrate_per_trace_percent is calculated from how many times that trace is revisited. 
    # hardcoded_minimum_hitrate is a hardcoded integer, such as 1 or 2'''
    
    output_horizons=[]
    output_horizons_binary=[]
    
    if hitrate_per_trace_percent ==0 and hardcoded_minimum_hitrate == 0:
        print('No filtercriterions given. User has to set one or both of hitrate_per_trace_percent andhardcoded_minimum_hitrate')

    else:
        if hitrate_per_trace_percent > 0:
            hitrate_per_trace_percent=hitrate_per_trace_percent/100
             # only if hitrate_per_trace_percent is non None
            trace_iteration=[]
            for i in range (0, len(grid)):  
                reference_trace=seismic_data[:, grid[i][0], grid[i][1]]
                reference_trace=np.trim_zeros(reference_trace, 'b')       
                reflector=[]
                count=0
                for k in range (0, len(path_info)):
                    t1 = path_info[k][2]
                    t2 = path_info[k][3]
                    if t1 == grid[i][0] and t2 == grid[i][1]:
                        count=count+1
                trace_iteration.append(count)  
        print('filter criterion:', hitrate_per_trace_percent, '%. and hardcoded hitrate in general:', hardcoded_minimum_hitrate)

        
        it=0
        for reflector in input_horizons:
            new_reflector=[]
            new_binary_reflector=[]
            count_hits=[]
            for i in range (0, len(grid)):
                hits_per_trace=[]
                for ii in range (0, len(reflector)):
                    if grid[i][0] == reflector[ii][1] and grid[i][1] == reflector[ii][2]:
                        hit=reflector[ii][0]
                        hits_per_trace.append(hit) #count all hits       
                accepted_reflector_points=[]
                hits_per_trace=np.array(hits_per_trace)
                unique, counts = np.unique(hits_per_trace, return_counts=True)
                counts_per_hit= np.asarray((unique, counts)).T

                if hitrate_per_trace_percent > 0 and hardcoded_minimum_hitrate == 0:
                    MINIMIM_HIT_NUMBER=trace_iteration[i]*hitrate_per_trace_percent 
                elif hitrate_per_trace_percent == 0 and hardcoded_minimum_hitrate > 0:
                    MINIMIM_HIT_NUMBER=hardcoded_minimum_hitrate
                elif hitrate_per_trace_percent > 0 and hardcoded_minimum_hitrate > 0:
                    percent_to_value=trace_iteration[i]*hitrate_per_trace_percent 
                    criterion = np.add(hardcoded_minimum_hitrate, percent_to_value)
                    MINIMIM_HIT_NUMBER=np.max(criterion)
                for cph in counts_per_hit:
                    if cph[1] > MINIMIM_HIT_NUMBER:
                        count_hits.append(cph[1])
                        for _ in range (cph[1]):
                            accepted_reflector_points.append(cph[0])

                if accepted_reflector_points:
                    unique, counts = np.unique(accepted_reflector_points, return_counts=True)

                    if len(unique) == 1:
                        event=unique[0]
                        binary_trace=[event, grid[i][0], grid[i][1]]
                        new_binary_reflector.append(binary_trace)
                        for _ in range (counts[0]):
                            new_trace= [event, grid[i][0], grid[i][1]]
                            new_reflector.append(new_trace)  

                    if len (unique) > 1:
                        idx=np.argmax(counts) # get the most abundant
                        event=unique[idx]
                        binary_trace=[event, grid[i][0], grid[i][1]]
                        new_binary_reflector.append(binary_trace)
                        for _ in range (counts[idx]):
                            new_trace= [event, grid[i][0], grid[i][1]]
                            new_reflector.append(new_trace) 
                            
                        # ALSO INCLUDE CLOSE POINTS WITH HIGH HITRATE
                        #for kk in range(0, len(unique)):
                            #multiple_event=unique[kk]
                            #if np.abs(multiple_event-event) > 0 and np.abs(multiple_event-event) < 5: # and close points
                                #binary_trace=[multiple_event, grid[i][0], grid[i][1]]
                                #new_binary_reflector.append(binary_trace)
                                #for _ in range (counts[kk]):
                                    #new_trace= [multiple_event, grid[i][0], grid[i][1]]
                                    #new_reflector.append(new_trace) 
                                    
            if print_output==True:
                print('horizon',it, 'reflector points:', len(new_binary_reflector))
                
            output_horizons_binary.append(new_binary_reflector)
            output_horizons.append(new_reflector)
            it=it+1
    return output_horizons_binary, output_horizons




from matplotlib import pyplot as plt
from skimage.morphology import square, disk, dilation

def ShowHorizon (seismic_data, horizons, horizon_number, view, line_number, selem):
    correlated_cube= np.zeros(seismic_data.shape, dtype=np.float32)
    if horizon_number == 'all':
        num=0
        correlated_cube= np.zeros(seismic_data.shape, dtype=np.float32)
        for i in range (0, len(horizons)):
            reflector=horizons[i]
            num=num+1
            for p in reflector:
                x=p[0]
                y=p[1]
                z=p[2]
                correlated_cube[x,y,z]=num
    else:
        reflector=horizons[horizon_number]
        for p in reflector:
                x=p[0]
                y=p[1]
                z=p[2]
                correlated_cube[x,y,z]=correlated_cube[x,y,z]+1 #+1 FOR Å FÅ EN IDE OM HITRATE

    fig, ax = plt.subplots()
    if view == 'inline':
        A=correlated_cube[:,line_number,:]
        B= seismic_data[:,line_number,:].copy()
    if view == 'crossline':
        A=correlated_cube[:,:,line_number]
        B= seismic_data[:,:,line_number].copy()
        
        
    A = dilation(A, selem)
    masked_data = np.ma.masked_where(A < 1, A)
    plt.imshow(B, 'gray', aspect='auto')
    plt.imshow(masked_data, 'jet', alpha=0.9, aspect='auto')
    
    if horizon_number == 'all':
        cbar = plt.colorbar(ax=None)
        cbar.set_label('Horizon number')
    else:
        cbar = plt.colorbar(ax=None)
        cbar.set_label('Autotracking accuracy (Hit rate)', rotation=270)
        cbar.ax.set_yticklabels([]) #no labels on colorbar
    
    plt.xlabel('intersecting lines')
    plt.title('Autotracked horizons superimposed onto seismic')
    plt.ylabel('Samples')
    plt.tight_layout()
    plt.show()
    
    return


import time
from IPython.display import clear_output

def regrid_horizons_DTW(seismic_data, grid, grid_step, filtered_horizons_binary, dense_grid, window_size):
    horizontal_warp_distance=window_size*grid_step
    dense_path_info=[]
    start = time.time()

    for i in range (0, len(grid)):  
        reference_trace=seismic_data[:, grid[i][0], grid[i][1]]
        reference_trace=np.trim_zeros(reference_trace, 'b')
        block_size=0

        for iii in range (0, len(dense_grid)):  
            if abs(grid[i][0]-dense_grid[iii][0])<=horizontal_warp_distance and abs(grid[i][1]-dense_grid[iii][1])<=horizontal_warp_distance:
                block_size=block_size+1
                matched_trace=seismic_data[:, dense_grid[iii][0], dense_grid[iii][1]]
                matched_trace=np.trim_zeros(matched_trace, 'b')
                path, sim= metrics.dtw_path(reference_trace, matched_trace)
                PI=[grid[i][0], grid[i][1], dense_grid[iii][0], dense_grid[iii][1], path]
                dense_path_info.append(PI)

        stop = time.time()
        t=(stop-start)/60
        clear_output()
        print('(',i,'out of:',len(grid),')',' time in minutes: ', t, 'Sub-block size:', block_size)


    regridded_horizons=[]
    it=0
    for reflector in filtered_horizons_binary:
        new_horizon=[]
        it=it+1
        for i in range (0, len(reflector)):
            events=[]
            reference_trace=seismic_data[:, reflector[i][1], reflector[i][2]]
            reference_trace=np.trim_zeros(reference_trace, 'b')
            event=reflector[i][0]
            for k in range (0, len(dense_path_info)):
                if dense_path_info[k][0] == reflector[i][1] and dense_path_info[k][1] == reflector[i][2]:
                    path= dense_path_info[k][4]
                    for ll in range (0, len(path)):

                        if path[ll][0]==event:
                            point=[path[ll][1], dense_path_info[k][2], dense_path_info[k][3]]
                            event_dense_path=path[ll][1]
                            new_horizon.append(point)
        regridded_horizons.append(new_horizon)
        
    return regridded_horizons, dense_path_info




from scipy import interpolate


def interpolateHorizons2D (seismic_data,unwrapped3D, horizons, dense_grid_step, phase_criterion, vertical_criterion, horizontal_criterion):
    interpolated_horizons=[]
    for horizon_number in range (0, len(horizons)):
        reflector=horizons[horizon_number]
        tracked_cube= np.zeros(seismic_data.shape, dtype=np.float32)
        interpolated_horizon=[]

        phase_criterion_interpolated=phase_criterion/dense_grid_step #criterion for interpolated neighbors
        for XL in range (0, seismic_data.shape[1]):
            for p in reflector:
                y=p[1]
                if y == XL:
                    x=p[0]
                    z=p[2]
                    tracked_cube[x,y,z]=1
            A=tracked_cube[:,XL,:]
            a=np.where(A>0)
            data = np.array(a)

            sorted_samples = [x for _,x in sorted(zip(data[1],data[0]))]
            sorted_inlines = np.sort(data[1])

            # for two neighbor points in XL direction
            for start in range (0, len(data[0])):
                x= sorted_inlines[start:start+2] #IL number
                y = sorted_samples[start:start+2]  #time sample

                if len(x) > 1 and np.abs(y[0]-y[1])<=vertical_criterion and np.abs(x[0]-x[1])<=horizontal_criterion:
                    p1_unwrapped=unwrapped3D[y[0], XL, x[0]]
                    p2_unwrapped=unwrapped3D[y[1], XL, x[1]]

                    if np.abs(p1_unwrapped-p2_unwrapped) <= phase_criterion:

                        f = interpolate.interp1d(x, y)
                        xnew = np.arange(np.min(x), np.max(x), 1) # range between two neighboring points
                        ynew = f(xnew)   # use interpolation function returned by `interp1d`
                        mask = ~np.isnan(xnew)
                        xnew = xnew[mask]
                        ynew = ynew[mask]
                        mask = ~np.isnan(ynew)
                        xnew = xnew[mask]
                        ynew = ynew[mask]
                        phase_values=[]
                        interpolated_points=[]
                        for p in range (0, len(xnew)):
                            int_x=np.int(xnew[p])
                            int_y=np.int(ynew[p])
                            phase_value=unwrapped3D[int_y, XL, int_x] 
                            phase_values.append(phase_value)  #print('SAMPLE:', int_x, 'inline:', int_y, 'phase:',  phase_value)
                            if len(phase_values) == 1:
                                ppp=[int_y, int_x]
                                interpolated_points.append(ppp)
                            # diff mellom hvert eneste nabopoints
                            if len(phase_values) > 1:
                                last_Two_phase_values=phase_values[-2:]

                                if np.abs(last_Two_phase_values[0]-last_Two_phase_values[1]) <= phase_criterion_interpolated:
                                    ppp=[int_y, int_x]
                                    interpolated_points.append(ppp)
                        for p in interpolated_points:
                            ppp=[np.int(p[0]), XL, np.int(p[1])]
                            interpolated_horizon.append(ppp)

        reflector=interpolated_horizon.copy()
        tracked_cube= np.zeros(seismic_data.shape, dtype=np.float32)

        for IL in range (0, seismic_data.shape[2]):
            for p in reflector:
                y=p[2]
                if y == IL:
                    x=p[0]
                    z=p[1]
                    tracked_cube[x,z,y]=1 #'''+1 FOR Å FÅ EN IDE OM HITRATE'''
            A=tracked_cube[:,:,IL]
            a=np.where(A>0)
            data = np.array(a)
            sorted_samples = [x for _,x in sorted(zip(data[1],data[0]))]
            sorted_inlines = np.sort(data[1])

            for start in range (0, len(data[0])):
                x= sorted_inlines[start:start+2]
                y = sorted_samples[start:start+2]
                if len(x) > 1 and np.abs(y[0]-y[1])<=vertical_criterion and np.abs(x[0]-x[1])<=horizontal_criterion:
                    p1_unwrapped=unwrapped3D[y[0], XL, x[0]]
                    p2_unwrapped=unwrapped3D[y[1], XL, x[1]]
                    if np.abs(p1_unwrapped-p2_unwrapped) <= phase_criterion:
                        f = interpolate.interp1d(x, y)
                        xnew = np.arange(np.min(x), np.max(x), 1) # range between two neighboring points
                        ynew = f(xnew)   # use interpolation function returned by `interp1d`
                        mask = ~np.isnan(xnew)
                        xnew = xnew[mask]
                        ynew = ynew[mask]
                        mask = ~np.isnan(ynew)
                        xnew = xnew[mask]
                        ynew = ynew[mask]

                        phase_values=[]
                        interpolated_points=[]
                        for p in range (0, len(xnew)):
                            int_x=np.int(xnew[p])
                            int_y=np.int(ynew[p])
                            phase_value=unwrapped3D[int_y, int_x, IL] 
                            phase_values.append(phase_value)  #print('SAMPLE:', int_x, 'inline:', int_y, 'phase:',  phase_value)
                            if len(phase_values) == 1:
                                ppp=[int_y, int_x]
                                interpolated_points.append(ppp)

                            # for every pair of neighbors of the interpolated points
                            if len(phase_values) > 1:
                                last_Two_phase_values=phase_values[-2:]
                                if np.abs(last_Two_phase_values[0]-last_Two_phase_values[1]) <= phase_criterion_interpolated:
                                    ppp=[int_y, int_x]
                                    interpolated_points.append(ppp)

                        for p in interpolated_points:
                            ppp=[np.int(p[0]), np.int(p[1]), IL]
                            interpolated_horizon.append(ppp)
        interpolated_horizons.append(interpolated_horizon)
    return interpolated_horizons




def interpolateHorizons3D (unwrapped3D, horizons, dense_grid_step, phase_criterion, vertical_criterion, close_neighbors):
    interpolated_horizons=[]
    from scipy import interpolate
    for horizon_number in range (0, len(horizons)):
        reflector=horizons[horizon_number]
        interpolated_points=[]
        for jj in range (0, len(reflector), 2):  # hvis denne er 1 blir noen doble?? ikke sikkert det går med 2

            p=reflector[jj] # for each point on the reflector
            p1_unwrapped=unwrapped3D[p[0], p[1], p[2]]
            points_for_interpolation=[]
            for p2 in reflector:
                if np.abs(p[1]-p2[1]) <= dense_grid_step and np.abs(p[2]-p2[2]) <= dense_grid_step and np.abs(p[0]-p2[0]) <= vertical_criterion:
                    points_for_interpolation.append(p2)

            if len(points_for_interpolation) >= close_neighbors:
                x = [item[1] for item in points_for_interpolation]
                y = [item[2] for item in points_for_interpolation]
                z = [item[0] for item in points_for_interpolation]

                sorted_x = np.unique(np.sort(x))
                sorted_y = np.unique(np.sort(y))
                sorted_z=[]
                for i in sorted_x:
                    for ii in sorted_y:
                        for pp in points_for_interpolation:
                            if i==pp[1] and ii==pp[2]:
                                sorted_z.append(pp[0])

                f = interpolate.interp2d(x, y, z)
                xnew = np.arange(np.min(x), np.max(x)+1, 1)
                ynew = np.arange(np.min(y), np.max(y)+1, 1)
                znew = f(xnew, ynew)
                znew=np.round(znew).astype(int)

                for i in range (0, len(xnew)):
                    for ii in range (0, len(ynew)):
                        iii=znew[ii][i]

                        if iii <= unwrapped3D.shape[0] and iii <= np.max(z) and iii>= np.min(z) :
                            newpoint=[iii, xnew[i], ynew[ii]]
                            pnew_unwrapped=unwrapped3D[iii, xnew[i], ynew[ii]]
                            if np.abs(p1_unwrapped-pnew_unwrapped) <= phase_criterion:
                                interpolated_points.append(newpoint)

        interpolated_horizons.append(interpolated_points)
    return interpolated_horizons

