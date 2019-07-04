from .kunz_cone import kunz_cone
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.integer_mod_ring import Integers

def get_incident_count(dimension):
    r"""
    Get incident count for each entry in the Vrepr()
    """
    P = kunz_cone(dimension)
    return [sum(1 for _ in Vrep.incident()) for Vrep in P.Vrepresentation()]

def order_of(g,b):
    """
    order of the orbit of <g> applied to b
    """
    counter = 0
    #b = apply_map(g*0+1,b)
    b = vector(b)
    b1 = apply_map(g,b)
    counter += 1
    while b1 != b:
        b1 = apply_map(g,b1)
        counter += 1
    return counter

def apply_map(g,b):
    """
    apply g to the vector b
    by multiplication on the indices
    """
    lst = [0]*len(b)
    for i in range(len(b)):
        lst[(g*(i+1))-1] = b[i]
    return vector(lst)

def sort_vertices(Vrep, incident_count, chunksize):
    Vrep = [vector(i) for i in Vrep]
    dimension = len(Vrep[0]) + 1
    dic = special_dic(dimension, Vrep, incident_count)
    counter = 0
    [shift_dics2, result_dic] = obtain_next_vertices(dic[next_key_of_special_dic(dic)], chunksize)
    shift_dics = [[i for i in j] for j in shift_dics2]
    next_key = next_key_of_special_dic(dic)
    counter += len(shift_dics[0])
    while next_key != -1:
        [next_shift_dics, next_result_dic] = obtain_next_vertices(dic[next_key], chunksize)
        for i in range(len(shift_dics)):
            shift_dics[i] += [x + counter for x in next_shift_dics[i]]
        for key in next_result_dic:
            result_dic[key + counter*chunksize] = next_result_dic[key]
        counter += len(next_shift_dics[0])
        next_key = next_key_of_special_dic(dic)
    return [shift_dics, result_dic]

def obtain_next_vertices(lst, chunksize):
    shift_dics = lst[0]
    all_items = lst[1]
    dic = {}
    counter = 0
    while counter < chunksize and all_items:
        items = all_items[0]
        all_items.remove(items)
        counter2 = 0
        for item in items:
            position = counter + counter2*chunksize
            dic[position] = item[1]
            counter2 += 1
        counter += 1
    return [shift_dics, dic]

def next_key_of_special_dic(dic):
    """
    Determine the next key of special dic.

    This is the key with the element of maximal appearance (except for the vertex of the cone).
    """
    maximal = -1
    result = -1
    for key in dic.keys():
        if len(dic[key][1]) == 0:
            current = -1
        elif sum(key) == len(key) and sum(dic[key][1][0][0][1]) == 0:
            # in this case this is the first element is the vertex of our cone
            if len(dic[key][1]) == 1:
                current = 0
            else:
                current = dic[key][1][1][0]
        else:
            current = dic[key][1][0][0]
        if current > maximal:
            maximal = current
            result = key
    return result

def special_dic(dimension, Vrep, incident_count):
    dic = {}
    r = Integers(dimension)
    gens = r.unit_gens()
    for i in range(len(Vrep)):
        vector = Vrep[i]
        typ = get_type(vector, gens)
        if not typ in dic:
            dic[typ] = []
        dic[typ].append((incident_count[i],vector))
    for key in dic:
        dic[key] = manipulate_list(dic[key], gens)
    return dic

def manipulate_list(lst, gens):
    """
    Sorts list according to count
    """
    unts = units(gens)
    lst.sort(reverse=True)
    lst2 = []
    while lst:
        item = lst[0]
        little_list = []
        lst.remove(item)
        little_list.append(item)
        for unit in unts:
            next_item = (item[0],apply_map(unit,item[1]))
            if next_item in lst:
                lst.remove(next_item)
                little_list.append(next_item)
        lst2.append(little_list)
    an_orbit = [i[1] for i in little_list]
    shift_dics = [shift_dic(an_orbit, gen) for gen in gens]
    return [shift_dics,lst2]

def shift_dic(items, gen):
    """
    items are the orbit of an element,
    returns the permutation induced by gen
    """
    result = []
    for i in range(len(items)):
        item = apply_map(gen, items[i])
        result.append(items.index(item))
    return result

def units(gens):
    order = gens[0].order()
    if len(gens) > 1:
        sofar = units(gens[1:])
    else:
        sofar = set((gens[0]**order,))
    return set(gens[0]**i*el for i in range(order) for el in sofar)


def get_type(vector, gens):
    """
    Determines the orbit length with respect to each generator.
    """
    return tuple(order_of(g, vector) for g in gens)
