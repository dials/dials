def filter_reflections_by_roi_volume(region_of_interest, percent):
    """Filter the reflections by roi volume.
    
    Calculate the volume of each reflection, filter out the 1.0-percent largest
    reflection volumes and return an array containing True/False if the volume
    if valid.
    
    Args:
        region_of_interest The array of rois
        percent The percent of volumes to keep (0.0 < percent <= 1.0)
    
    Returns:
        A boolean array, True/False is the roi valid
    
    """
    from dials.array_family import flex
    from heapq import nlargest

    # Check given percentage
    if percent <= 0 or percent > 1.0:
        raise ValueError
        
    # A Calculate the volume of each region of interest
    calculate_roi_volume = lambda roi: ((roi[1] - roi[0]) * 
                                        (roi[3] - roi[2]) * 
                                        (roi[5] - roi[4]))
    volume = map(calculate_roi_volume, region_of_interest)
    
    # Calculate the volume limit below which 99% of reflections are
    n_reflections = len(volume)
    volume_limit = nlargest(int((1.0 - percent) * n_reflections), volume)[-1]
    
    # Create an array which is true if reflection volume is below the limit
    result = flex.bool(n_reflections)
    for i, v in enumerate(volume):
        result[i] = v < volume_limit
    
    return result
