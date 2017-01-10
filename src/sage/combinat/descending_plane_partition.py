r"""
Descending plane partitions

AUTHORS:

- Colton Keller and Jessica Striker (2017): Initial version
"""
#*****************************************************************************
#       Copyright (C) 2017 Colton Keller, Jessica Striker <jessicapalencia@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
# python3

from copy import copy

class DescendingPlanePartition():
    r"""
    A descending plane partition. 
    
    A descending plane partition (DPP) is a strict shifted plane partition with 
    the additional restrictions that: the number of parts in each row is less 
    than the greatest part in that row, and the number of parts in each row is 
    greater than or equal to the largest part in the next row.
    
    These were introduced in [MiRoRu]_.

    REFERENCES:

    .. [MiRoRu] \W. H. Mills, David P Robbins, Howard Rumsey Jr.,
       *Alternating sign matrices and descending plane partitions*,
       Journal of Combinatorial Theory, Series A,
       Volume 34, Issue 3, May 1983, Pages 340--359.
       http://www.sciencedirect.com/science/article/pii/0097316583900687
    """
    def __init__(self, DPP):
        """
        Initialize ``self``.
        
        EXAMPLES::
            
            sage: DPP = DescendingPlanePartition([[5,4,4],[3,2]])
            sage: DPP #not sure if this is correct
            A descending plane partition: 
            544
             32
        """
        if DescendingPlanePartitions.is_dpp(DPP):
            self.DPP = DPP
        else:
            raise ValueError("Invalid descending plane partition")
        #self.DPP = DPP
        
    def __repr__(self):
        r"""
        Return a representation of ``self``.
        
        EXAMPLES::
            
            sage: DPP = DescendingPlanePartition([[5,4,4],[3,2]])
            sage: DPP #not sure if this is correct
            A descending plane partition: 
            544
             32
        """
        return "A descending plane partition: \n" + self._pretty_string()
    #self.dpp= DPP
    
    def _pretty_string(self):
        r"""
        Return a formatted string representing ``self``.
        
        EXAMPLES::
            
            sage: DPP = DescendingPlanePartition([[5,4,4],[3,2]])
            sage: DPP._pretty_string()
            '544\n 32\n'
        
        """
        rstr = ""
        for i in range(len(self.DPP)):
            for j in range(i):
                rstr += " "
            for element in self.DPP[i]:
                rstr += str(element)
            rstr += "\n"
        return rstr
    
    def pp(self):
        r"""
        Return a formatted string representing ``self``.
       
        """
        print self._pretty_string()
    
    def _is_columnar_descending(self, dpp):
        """
        Return True if each column of the descending plane partition is strictly 
        decreasing and False otherwise.

        """
        #len == 2 len dpp[-2]>len dpp[-1]
        if len(dpp[-1])>=len(dpp[-2]):
            return False
        for i in range(len(dpp[-1])):
            if not dpp[-1][i] < dpp[-2][i+1]:
                return False
        return True

    def _increment_row(self, dpprow):
        """
        Return the next descending plane partition row given a starting 
        descending plane partition row.

        """
        if type(dpprow) is type([]):
            if len(dpprow)< (dpprow[0] - 1):
                dpprow.append(1)
                return dpprow
            else:
                for i in range(len(dpprow)):
                    if dpprow[i] == 1:
                        dpprow[i] = dpprow[i]+1
                        return dpprow[:i+1]
                if dpprow[-1] < dpprow[0]:
                    for i in range(len(dpprow)):
                        if dpprow[i]==dpprow[-1]:
                            dpprow[i] = dpprow[i] + 1
                            return dpprow[:i+1]
                return [dpprow[0] + 1]
    
    def getNextDPP(self):
        """
        Return the next descending plane partition in the iterator's sequence. 
        
        EXAMPLES::
        
            sage: DPP = DescendingPlanePartition([[5,4,4],[3,3],[2]])
            sage: DPP.getNextDPP()
            [[5, 4, 4, 1]]
        
        """
        return self._increment_2_rows(self.DPP)
    
    def _increment_2_rows(self,DPP):
        """
        Increment the descending plane partition by recursively traversing the 
        rows of the descending plane partition.

        """
        dpp = copy(DPP)
        if len(dpp)>2:
            subdpp = self._increment_2_rows(dpp[1:])
            dpp = dpp[:1]
            dpp.extend(subdpp)
        elif len(dpp) ==1:
            if len(dpp[0]) == 0:
                dpp = [[2]]
                return dpp
            if len(dpp[0])>1 and dpp[0][1]>2:
                dpp.append([2])
                return dpp
            else:
                dpp[0] = self._increment_row(dpp[0])
                return dpp
        else:
            if len(dpp[-1])>1 and dpp[-1][1] > 2:
                dpp.append([2])
                return dpp
            dpp[-1] = self._increment_row(dpp[-1])
        while len(dpp[-2]) >= dpp[-1][0]:
            if not self._is_columnar_descending(dpp):
                dpp[-1] = self._increment_row(dpp[-1])
            else:
                return dpp
        if (len(dpp[-2]) < dpp[-1][0]):
            dpp = dpp[:-1]
            dpp[-1] = self._increment_row(dpp[-1])
        return dpp
    
    def increment(self):
        
        """
        Return the next descending plane partition in the iterator's sequence.
        
        EXAMPLES::
            
            sage: DPP = DescendingPlanePartition([[5,4,4],[3,3],[2]])
            sage: DPP.increment()
            sage: DPP
            A descending plane partition:
            5441
            
        """
        
        self.DPP = self.getNextDPP()
        
    def sum(self):
        
        """
        Return the sum of ``self``.
        
        EXAMPLES:
            sage: DPP = DescendingPlanePartition([[5,4,4],[3,3],[2]])
            sage: DPP.sum()
            21
        
        """
        
        sum = 0
        for row in self.DPP:
            for entry in row:
                sum = sum + entry
        return sum

    
    def special_parts(self):
        
        """
        Return the number of special parts in ``self``.

        A special part of a descending plane partition is a part, $p$, of the 
        descending plane partition such that $p$ is less than its position from 
        the left of start of the row.
        
        EXAMPLES::

            sage: DPP = DescendingPlanePartition([[5,4,4],[3,3],[2]])
            sage: DPP.special_parts()
            0
            sage: DPP = DescendingPlanePartition([[5,2,2,2]])
            sage: DPP.special_parts()
            2
        
        """
        
        sparts = 0
        for i in range(len(self.DPP)):
            for j in range(len(self.DPP[i])):
                if self.DPP[i][j] <= j:
                    sparts = sparts +1
        return sparts

    def regular_parts(self):
        
        """
        Return the number of regular parts in ``self``.

        A regular part of a descending plane partition is a part, $p$, of the 
        descending plane partition such that $p$ is greater than or equal to its 
        position from the left of start of the row.
        
        EXAMPLES::

            sage: DPP = DescendingPlanePartition([[5,4,4],[3,3],[2]])
            sage: DPP.regular_parts()
            6
            sage: DPP = DescendingPlanePartition([[5,2,2,2]])
            sage: DPP.regular_parts()
            2
        
        """
        
        rparts = 0
        for i in range(len(self.DPP)):
            for j in range(len(self.DPP[i])):
                if self.DPP[i][j] > j:
                    rparts = rparts +1
        return rparts
    
    def __getitem__(self, rowIndex):
        
        """
        Return the row with index rowIndex of the descending plane partition.
        
        EXAMPLES::
        
            sage: DPP = DescendingPlanePartition([[5,4,4],[3,3],[2]])
            sage: DPP[1]
            [3, 3]
        """
        
        return self.DPP[rowIndex]
    
    def __copy__(self):
        
        """
        Return a copy of the descending plane partition.

        EXAMPLES::

            sage: DPP = DescendingPlanePartition([[5,4,4],[3,2]])
            sage: DPP.__copy__()
            A descending plane partition: 
            544
             32
            <BLANKLINE>
        """
  
        return DescendingPlanePartition(self.DPP)
    
    def is_Catalan(self):
        
        """
        Return True if this DPP is a Catalan DPP.
        
        A Catalan DPP is a row $a_{i,1} \dots a_{i,{\lambda}}$ of a descending 
        plane partition such that each entry satisfies $a_{i,j} \leq a_{i,1} - j + 1$
        
        EXAMPLES::

            sage: DPP = DescendingPlanePartition([[5,4,2]])
            sage: DPP.is_Catalan()
            True
            sage: DPP = DescendingPlanePartition([[5,4,4],[3,2]])
            sage: DPP.is_Catalan()
            False
        
        """
        
        if len(self.DPP) == 1:
            first = self.DPP[0][0]
            flag = True
            for i in range(len(self.DPP[0])):
                if self.DPP[0][i] > first - i:
                    flag = False
                    break
            return flag
        else:
            return False
        
    def get_catalan_dpp_path(self):
        
        """
        Return the Catalan DPP path that corresponds to ``self`` if this is a 
        valid Catalan descending plane partition.
        
        A Catalan DPP path is a word composed of $1$'s and $-1$'s such that the 
        partial sum of each word is always greater than or equal to $0$ when 
        taken from left to right, and the complete sum is greater than $0$ if 
        the word is non-empty.
        
        EXAMPLES::

            sage: dpp = DescendingPlanePartition([[3,2]])
            sage: dpp.get_catalan_dpp_path()
            [1, -1, 1]
        
        """
        
        if self.is_Catalan():
            dpprow = self.DPP[0]
            path = []
            for i in range(len(dpprow)-1):
                dif = dpprow[i]-dpprow[i+1]
                for j in range(dif):
                    path.append(1)
                path.append(-1)
            for k in range(dpprow[-1]-1):
                path.append(1)
            return path
        
    @staticmethod
    def catalan_dpp_path_to_dpp(path):
        
        """
        Return the corresponding Catalan DPP from the given Catalan DPP Path.
        
        EXAMPLES::

            sage: DescendingPlanePartition.catalan_dpp_path_to_dpp([1,-1,1])
            A descending plane partition:
            32
        """
        
        row = []
        descents = 0
        for entry in path:
            if entry == 1:
                descents += 1
        n = descents + 1
        row.append(n)
        for i in range(len(path)):
            if path[i] == 1:
                n -= 1
            else:
                row.append(n)
        return DescendingPlanePartition([row])
    
    def _latex_(self):
        
        """
        Return the latex representation of ``self``.

        """
        
        xaxis = "210"
        yaxis = "-30"
        zaxis = "90"
        def topside(x,y,z):
            return "\\fill[fill=white, draw=black,shift={(" + xaxis + ":" + str(x) + ")},shift={(" + yaxis + ":" + str(y) + ")},  \nshift={(" + zaxis + ":" + str(z) + ")}] (0,0) -- (30:1) -- (0,1) --(150:1)--(0,0);\n"
        def leftside(x,y,z):
            return "\\fill[fill=lightgray, draw=black,shift={("+xaxis+":"+ str(x) + ")},shift={("+yaxis +":"+str(y)+")},  \nshift={("+zaxis+":"+str(z)+")}] (0,0) -- (0,-1) -- (210:1) --(150:1)--(0,0);\n"
        def rightside(x,y,z):
            return "\\fill[fill=darkgray, draw=black,shift={("+xaxis+":"+str(x) + ")},shift={("+yaxis+":"+str(y)+")},  \nshift={(" +zaxis+ ":" + str(z) + ")}] (0,0) -- (30:1) -- (-30:1) --(0,-1)--(0,0);\n"
        def cube(x,y,z):
            return topside(x,y,z) + leftside(x,y,z) + rightside(x,y,z)
        myC = copy(self.DPP)
        for i in range(len(myC)):
            for j in range(i):
                myC[i].insert(0,0)
        x = -1
        latstr = ""
        for row in myC:
            x += 1
            y = -1
            for element in row:
                y += 1
                z = -1
                if element != 0:
                    for k in range(element):
                        z += 1
                        latstr+=cube(x,y,z)
        return "\\begin{tikzpicture}\n" + latstr + "\\end{tikzpicture}\n"
    
class DescendingPlanePartitions():
    
    def __init__(self, n = None, is_Catalan = False):
        
        """
        Initialize ``self`` to represent either regular or Catalan DPPs with 
        largest part less than or equal to `n`.
        
        EXAMPLES::

            sage: DPPs = DescendingPlanePartitions(3, True)
            sage: DPPs
            Catalan Descending Plane Partitions with largest part at most 3
        """
        
        if n is None:
            #infinite order
            self.n = n
        elif n >=0:
            self.n = n
        else:
            raise ValueError("Descending plane partitions must have a positive greatest part")
        self.catalan = is_Catalan
    
    @staticmethod
    def is_dpp(dpp):
        
        """
        Return True if dpp is a valid descending plane partition.
    
        EXAMPLES::
        
        sage: DescendingPlanePartitions.is_dpp([[4,4],[2]])
        True
        """
        
        if type(dpp) is type([]):
            #check for acceptable parameter
            for i in range(len(dpp)):
                #loop over all sublists (each row in dpp)
                if type(dpp[i]) is type([]):
                    #check for acceptable parameter
                    for j in range(len(dpp[i])):
                        #loop over each entry in the sublist (each column in dpp)
                        if j != 0:
                            #if the row has more than one entry,
                            #check for descending rows
                            if dpp[i][j]>dpp[i][j-1]:
                                return False
                        else:
                            #on the first entry of each row, check row length less than first entry
                            if not len(dpp[i]) < dpp[i][0]:
                                return False
                            if i != 0:
                                if len(dpp[i-1]) < dpp[i][0]:
                                    return False

                        if i != 0:
                            #if there are more than one rows
                            if len(dpp[i-1]) > j+1:
                                #if there is an entry directly above current entry
                                if not dpp[i-1][j+1] > dpp[i][j]:
                                    #strictly descending vertically check
                                    return False
                            else:
                                #no entry directly above, strictly descending vertically fails
                                return False
                else:
                    #unacceptable parameter
                    return False

        else:
            if type(dpp) is type(DescendingPlanePartition):
                return True
            else:
                return False
        return True #no errors found
    
    @staticmethod
    def is_Catalan_dpp(dpp):
        """
        Return True if a DPP is a Catalan DPP.
        
        A Catalan DPP is a row $a_{i,1} \dots a_{i,{\lambda}}$ of a descending 
        plane partition such that each entry satisfies $a_{i,j} \leq a_{i,1} - j + 1$
        
        EXAMPLES::

            sage: DescendingPlanePartitions.is_Catalan_dpp([[4,3]])
            True
            sage: DescendingPlanePartitions.is_Catalan_dpp([[3,3]])
            False
        
        """
        if DescendingPlanePartitions.is_dpp(dpp):
            DPP = dpp
            if not isinstance(dpp,DescendingPlanePartition):
                DPP = DescendingPlanePartition(dpp)
            return DPP.is_Catalan()
        else:
            return False
                
    
#    def __contains__(self, DPP):
#        
#        """
#        Return True if a DPP is a part of this DescendingPlanePartitions object.
#        
#        EXAMPLES::
#            
#            sage: DPPs = DescendingPlanePartitions(3)
#            sage: DPP = DescendingPlanePartition([[3,3]])
#            sage: DPP in DPPs
#            True
#        
#        """
#        
#        if DescendingPlanePartitions.is_dpp(DPP):
#            if self.n is None:
#                return True
#            elif len(DPP[0]) ==0:
#                return True
#            else:
#                if DPP[0][0] <= self.n:
#                    return True
#                else:
#                    return False
#        else:
#            return False
    
    def __iter__(self):
        """
        Iterate over ``self``.

        """
        DPP = DescendingPlanePartition([[]])
        yield copy(DPP)
        if self.n is None:
            while(True):
                DPP.increment()
                yield copy(DPP)
        else:
            if self.n < 2:
                return
            while(True):
                DPP.increment()
                if(DPP[0][0]<= self.n):
                    yield copy(DPP)
                else:
                    return
    
    def __repr__(self):
        """
        Return the string representing this object.

        """
        if self.n == None:
            if self.catalan:
                return "The class of all Catalan Descending Plane partitions"
            else:
                return "The class of all Descending Plane Partitions"
        else:
            if self.catalan:
                return "Catalan Descending Plane Partitions with largest part at most {}".format(int(self.n))
            else:
                return "Descending Plane Partitions with largest part at most {}".format(int(self.n))
