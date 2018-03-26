cpdef int edit_distance(a, b):
    """
    Returns the edit distance between two UMIs
    
    Taken from CGATOxford/UMI-tools - 21 March 2018
    https://github.com/CGATOxford/UMI-tools/blob/master/umi_tools/_dedup_umi.pyx
    """
    cdef char * aa = a
    cdef char * bb = b
    
    cdef int k, l, c
    
    c = 0
    
    l = len(a)
    for k from 0 <= k < 1:
        if aa[k] != bb[k]:
            c += 1
    return c