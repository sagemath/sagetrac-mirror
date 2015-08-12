def set2ind(s):
    r"""
	Return the rank of a set under the lexicographic order.
    
    INPUT:
    
    - ``s`` -- An iterable containing natural numbers   
    
    OUTPUT:
    
    - An integer giving the rank
    
    EXAMPLES::

            sage: set2ind({})
            0
            sage: set2ind({1})
            1
            sage: set2ind({1,4})
            9

	"""
    str=''
    if len(s)==0:
        return 0 
    for i in range(1, max(s)+1):
        if i in set(s):
            str=str+'1'
        else:
            str=str+'0'
    return int(str[::-1],2)
