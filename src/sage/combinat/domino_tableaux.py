
from copy import copy

def zeros(n):
    """returns a list of n zeros"""
    if n == 0:
        return []
    else:
        return zeros(n - 1) + [0]

def reverse(L):
    """returns the reverse of a list L"""
    if L == []:
        return []
    else:
        return reverse(L[1:]) + [L[0]]    
    
    
def Remove_Element(e,L):
    """Removes the element e from list L"""
    if L == []:
        return []
    if L[0] == e:
        return Remove_Element(e,L[1:])
    else:
        return [L[0]] + Remove_Element(e,L[1:])

def mergesort(L,M):
    """merge sorts the lists L,M"""
    if L == []:
        return M
    if M == []:
        return L
    if L[0] <= M[0]:
        return [L[0]] + mergesort(L[1:],M)
    else:
        return [M[0]] + mergesort(L,M[1:])
    
def sort(L):
    """sorts a list L"""
    n = len(L)
    if L == []:
        return []
    if n == 1:
        return L
    else:
        m = n//2
        return mergesort(sort(L[:m]), sort(L[m:]))
        

def Remove_Repeats(L):
    """Removes the repeats that show up in list L and sorts it"""
    M = sort(L)
    if M == []:
        return M
    else:
        ans = [M[0]]
        for j in range(1,len(M)):
            if not(M[j] == M[j - 1]):
                ans = ans + [M[j]]
        return ans 
    

def Tab_Elements(Tab):
    """returns the elements that show up in tableaux Tab"""
    List = [] # Intitialize list
    for i in range(0,len(Tab)):
        for j in range(0,len(Tab[i])):
            List = List + [Tab[i][j]]
    return Remove_Repeats(List)

def Cell_Present(Tab,Cell):
    """returns True if the given cell shows up in Tab"""
    [a,b] = Cell
    if a < len(Tab) and b < len(Tab[a]):
        return True
    else:
        return False

def Tab_Entry(Tab,Cell):
    """returns the entry that is in given cell"""
    [a,b] = Cell
    if a < len(Tab) and b < len(Tab[a]):
        return Tab[a][b]
    else:
        print "Entry not found at", Cell, "in", Tab
    
def Transpose_Part(lam):
    """returns the transpose of a partition lam, where the parts are listed in decreasing order"""
    if lam == []:
        return []
    else:
        ans = []
        m = len(lam)   # Fix length as m
        for k in range(1,lam[0] + 1):
            i = 0
            while i < m and lam[i] >= k:
                i = i + 1  # Increase index by 1 as long as it is value is bigger than k
            ans = ans + [i]
        return ans

def shape(Tab):
    """returns the partition shape associated to tableaux Tab"""
    return map(len,Tab)

def Col_heights(Tab):
    """returns the list of column heights to tableaux Tab"""
    return Transpose_Part(shape(Tab))

def Col_height(Tab,j):
    """returns the height of column j in tableaux Tab"""
    lam = shape(Tab)
    m = len(lam)   # Fix length as m
    i = 0
    while i < m and lam[i] >= j + 1:
        i = i + 1  # Increse index by 1 as long as it is value is bigger than j
    return i
    
def Tab_Cells(Tab,k):
    """return the cells in tableaux Tab that have entry k"""
    ans = []
    height = len(Tab)
    for i in range(0,height):
        for j in range(0,len(Tab[i])):
            if Tab[i][j] == k:  # Checks whether each entry is equal to k
                ans = ans + [[i,j]]  # If so, add it to the list
    return ans

def Transpose(Tab):
    """returns the transpose of Tab, swapping rows and columns"""
    Part = shape(Tab)
    new_shape = Transpose_Part(Part)
    n = len(new_shape)      # Record number of rows in transpose
    new_Tab = zeros(n)
    for i in range(0,n):
        new_Tab[i] = zeros(new_shape[i])  # Fill tableuax with zeros in shape of transpose
    for i in range(0,n):
        for j in range(0,new_shape[i]):
            new_Tab[i][j] = Tab[j][i]
    return new_Tab

def Restrict_Tab(Tab,k):
    """restrict the tableaux to entries that are at most k"""
    new_height = 0
    while new_height < len(Tab) and Tab[new_height][0] <= k:
        new_height = new_height + 1   # Find the new height of tab, that is how many entries in first col are <= k
    new_Tab = zeros(new_height)  # Create list of 0s with right length to put tab in
    for i in range(0,new_height):   # For each row i
        row_length = 0
        m = len(Tab[i])
        while row_length < m and Tab[i][row_length] <= k:
            row_length = row_length + 1    # Row length counts the number of entris in each row that are <= k
        new_Tab[i] = Tab[i][:row_length]    # Replace old row by new row
    return new_Tab

def Check_Rows_Increase(Tab):
    """returns True if the Rows weakly increase and False otherwise"""
    height = len(Tab)
    for i in range(0,height):
        for j in range(0,len(Tab[i]) - 1):
            if Tab[i][j] > Tab[i][j + 1]:  # Checks whether the rows weakly increase or not
                return False  # If they don't, return false
    return True
    
def Check_Cols_Increase(Tab):
    """returns True if the Cows weakly increase and False otherwise"""
    return Check_Rows_Increase(Transpose(Tab))
    
def Check_Tableaux(Tab):
    """returns True if the rows and columns increase"""
    return Check_Rows_Increase(Tab) and Check_Cols_Increase(Tab)
    
def Check_Domino_Tab(Tab):
    """sees whether tableaux Tab is a domino tableaux with entry list L or not"""
    L = Tab_Elements(Tab)  # Let L list the elements of the tableaux T
    if Check_Tableaux(Tab) == False:
        return False     # if we are not a tableaux, it is false
    for k in L:
        Cells = Tab_Cells(Tab,k)
        if not(len(Cells) == 2):
            return False
        if len(Cells) == 2:
            A = Cells[0]
            B = Cells[1]     # Finds where each entry occurs in the tab
            Diff = [B[0] - A[0],B[1] - A[1]]   # Find the difference between the Cells
            if not(Diff == [0,1]) and not(Diff == [1,0]):
                return False     # If the entries are not adjacent return false
    return True
    
def Add_Domino(Tab,C1,C2,k):
    """Adds a new domino with entry k in cells C1 and C2, making a larger tableaux, assumes this is possible"""
    new_Tab = copy(Tab)  # Intitialize new Tableaux, not attached to old one
    i = C1[0]
    if C2[0] == i:   # if the cells are in the same row
        if i < len(Tab):  # If row i is still right of the tableaux ...
            j = len(Tab[i])  # j is the corresponding column after this row
            if (C1[1] == j) and (C2[1] == j + 1):   # check to see if they appear at the end of the row
                if i == 0 or len(Tab[i - 1]) >= (j + 1):    # Checks to see this doesn't violate partition shape
                    new_Tab[i] = Tab[i] + [k,k]    # If they do, add horiz domino to the end of this tableaux
                    return new_Tab
        if i == len(Tab): # If the row appears at the end
             if C1[1] == 0 and C2[1] == 1:   # check to see if they appear at the end of the row
                if i == 0 or len(Tab[i - 1]) > 1:    # Checks to see this doesn't violate partition shape
                    new_Tab = Tab + [[k,k]]    # If it doesn't, add domino to the end of this tableaux
                    return new_Tab

    if C2[0] == i + 1:   # if the cells are in adjacent rows
        if i < len(Tab) - 1:  # If row i is still right of the tableaux...
            j = len(Tab[i])   # j is the corresponding column after this row
            # print "j =", j
            if len(Tab[i + 1]) == j:  # These adjacent rows must be the same length
                if (C1[1] == j) and (C2[1] == j) :   # check to see if they appear at the end of the row
                    if i == 0 or len(Tab[i - 1]) >= (j + 1):    # Checks to see this doesn't violate partition shape
                        new_Tab[i] = Tab[i] + [k]
                        new_Tab[i + 1] = Tab[i + 1] + [k]   # If it doesn't, add vert domino to the end of this tableaux
                        return new_Tab
        if i == len(Tab): # If the row appears at the end
            if C1[1] == 0 and C2[1] == 0:   # check to see if they appear at the end of the row
                # print "Adding to First Column"
                new_Tab = Tab + [[k]] + [[k]]    # If it doesn't, add vertical domino below tableaux
                return new_Tab
    print "Unable to add this domino at" 
    print C1
    print C2
    return
     
    
def Add_Domino_2(Tab,C1,C2,k):
    """Adds a new domino with entry k in cells C1 and C2 if it can - But also does adjustments so that it adds back in a close place,
    as in the generalized RSK algorithm"""
    i = C1[0]
    if C2[0] == i:   # if the cells are in the same row
        if i < len(Tab):  # If row i is still right of the tableaux ...
            j = len(Tab[i])  # j is the corresponding column after this row
            # print "j =", j
            if (C1[1] == j-1) and (C2[1] == j):     # If it hits one suqare at the end of the row
                return Add_Domino(Tab,[i,j],[i + 1,j],k)   # Rotate to a veritcal domino at the end
            if (C1[1] == j - 2) and (C2[1] == j - 1):    # If it hits last two squares at the end of the row
                # print "Row Drop"
                if i < len(Tab) - 1:
                    m = len(Tab[i + 1])
                if i == len(Tab) - 1:
                    m = 0
                return Add_Domino(Tab,[i+1,m],[i+1,m+1],k)   # Bump the domino down to next row
        return Add_Domino(Tab,C1,C2,k)
            
    if C2[0] == (i + 1):   # if the cells are in the same row
        if i < len(Tab):  # If row i is still right of the tableaux ...
            j = C1[1]  # j is the corresponding column thid domino appears in, this is different from above
            # print "j =", j
            r = Col_height(Tab,j)  # r is column height of j-th col 
            # print "r =", r
            s = Col_height(Tab,j + 1)   # s is column height of (j + 1)-th col
            # print "s =", s
            if C1 == [r-1,j] and C2 == [r,j] :     # If it hits one square at the end of the col
                # print "Turn horiz"
                return Add_Domino(Tab,[r,j],[r,j + 1],k)   # Rotate to a horiz domino at the end
            if C1 == [r-2,j] and C2 == [r-1,j]:    # If it hits last two squares at the end of the col
                # print "Column Drop"
                return Add_Domino(Tab,[s,j + 1],[s+1,j + 1],k)   # Bump the vertical domino to the right
        return Add_Domino(Tab,C1,C2,k)
        
def Add_Domino_3(Tab,C,k):
    """Adds a new domino with entry k in cells C = [C1,C2] if it can - But also does adjustments so that it adds back in a close place,
    as in the generalized RSK algorithm""" 
    C1 = C[0]
    C2 = C[1]
    return Add_Domino_2(Tab,C1,C2,k)
        
def Alpha_Bump(Tab,k,sign):
    """Bumps k into tab T, to get a new Tableaux with k into it"""
    if sign == 1:  # In this case we do horizontal insertion and row bumping
        Elts = Tab_Elements(Tab)   # Creates a list of the elements of tableaux Tab, sorted
        New_Elts = mergesort(Elts,[k])   # Add k to the list of elements - we can mergesort here 
        new_Tab = []
        for i in New_Elts:
            if i == k:
                if new_Tab == []:
                    C1 = [0,0]
                    C2 = [0,1]
                else:
                    C1 = [0,len(new_Tab[0])]
                    C2 = [0,len(new_Tab[0]) + 1]    # If we are inserting k, we add it to the first row
            else:
                Cells = Tab_Cells(Tab,i)
                C1 = Cells[0]
                C2 = Cells[1]                # We try to add i back to where it was, but make adjustments if necessary
            new_Tab = Add_Domino_2(new_Tab,C1,C2,i)     # We add back the cells this way in order
            # print new_Tab
            
    if sign == -1:  # In this case we do vertical insertion and column bumping
        Elts = Tab_Elements(Tab)   # Creates a list of the elements of tableaux Tab, sorted
        New_Elts = mergesort(Elts,[k])   # Add k to the list of elements - we can mergesort here 
        new_Tab = []
        for i in New_Elts:
            if i == k:
                if new_Tab == []:
                    C1 = [0,0]
                    C2 = [1,0]
                else:
                    C1 = [len(new_Tab),0]
                    C2 = [len(new_Tab) + 1,0]    # If we are inserting k, we add it to the first row
            else:
                Cells = Tab_Cells(Tab,i)
                C1 = Cells[0]
                C2 = Cells[1]                # We try to add i back to where it was, but make adjustments if necessary
            new_Tab = Add_Domino_2(new_Tab,C1,C2,i)     # We add back the cells this way in order
            # print new_Tab
    return new_Tab

def Tab_Diff(Tab1,Tab2):
    """returns the difference between the two tableaux Tab1,Tab2,and Tab2 is bigger than Tab1"""
    lam1 = shape(Tab1)
    lam2 = shape(Tab2)   # Record shapes for tableaux
    m = len(lam1)
    n = len(lam2)   # list lengths
    Diff = []   # Intitialize List of differences
    lam1 = lam1 + zeros(n - m)  # Pad 
    for j in range(0,n):
        if not(lam1[j] == lam2[j]):
            Extra_Cols = range(lam1[j], lam2[j])  # add all entries in row j in Tab2 but not Tab2
            Diff = Diff + map(lambda x: [j,x], Extra_Cols)
    return Diff
    
            
def Domino_RSK(w):
    """Performs Domino RSK on the Word w, assumes that each element has an element and a sign in the form [a,1] or [a,-1]"""
    L = []  # Initialize left (insertion) Tableaux
    R = []  # Initialize right (recording) Tableaux
    n = len(w)
    for j in range(0,n):
        k = w[j][0]
        sign = w[j][1]
        New_L = Alpha_Bump(L,k,sign)   # Bump in the nest element in the insertion tableaux
        # print New_L
        Diff = Tab_Diff(L,New_L)
        R = Add_Domino_3(R,Diff,j+1)    # Record what cells got added in the recording tableaux
        # print R
        L = New_L
    return [L,R]


def Remove_Cell(Tab,Cell):
    """removes the Cell Cell from the tableaux"""
    [i,j] = Cell
    if Tab == []:
        print "Cannot remove from empty tableaux"
        return
    if i < len(Tab):
        if j == len(Tab[i]) - 1:     # if j is at the end of the row
            if j == 0:      # If row is last column
                if i == len(Tab) - 1:   # Check to see we are actually in the last row
                    return Tab[:i]   # If so, remove last row
            else:
                return Tab[:i] + [Tab[i][:j]] + Tab[(i + 1):]  # Else, we just remove last element from row
    print "Cannot remove Cell", Cell, "from", Tab
    return

def Remove_Domino(Tab,C1,C2):
    """removes the Domino C1,C2 from the tableaux"""
    Middle_Tab = Remove_Cell(Tab,C2)   # First, remove C2
    return Remove_Cell(Middle_Tab,C1)  # Then remove C1

def Change_Entry(Tab,Cell,k):
    """changes the given Cell in Tab to have entry k, includes adding a new entry"""
    [i,j] = Cell
    if i < len(Tab):
        if j <= len(Tab[i]):
            return Tab[:i] + [Tab[i][:j] + [k] + Tab[i][(j + 1):]] + Tab[(i + 1):]     # We first return rows above it
                            # Then, change j-th place in ith row, then return rest of the row
    if i == len(Tab):
        if j == 0:
            return Tab + [[k]]
    print "Cannot change Cell", Cell, "in", Tab
    return
                              
def New_Cells(Tab,C1,C2,D1,D2):
    """Take in insertion cells C1,C2, and avoidance cells D1,D2 - as in inverse to RSK, and returns new insertion and avoidance cells"""
    # new_Tab = Tab
    i = C1[0]
    j = C1[1]    # Store [i,j] as position of first cell in
    if C2 == D2:
        if C2[0] == i:  # Here domino we are adding is horizontal
            if C1 == D1:   # And domino we need to avoid is in same place
                if i == 0:
                    # print "Found domino to bump out"
                    return [1]  # In this case, we don't add anything back and just return Tab itself in a list
                                    # We add a 1 to indicate this came from a row bump
                else:
                    # print "Bumping up row domino"
                    r = len(Tab[i - 1])  # Find length of previous row
                    Cells = [[i-1,r-2],[i-1,r-1],[i-1,r-2],[i-1,r-1]]   # Inserion cells, avoidance cells are both bumped up to previous row
            if D1[0] == i - 1 and D1[1] == j + 1:
                Cells = [[i - 1,j],[i,j],[i - 1,j],[i - 1,j +1]]  # Insertion cells rotate up, Avoidance Cells rotate left
        if C2[0] == i + 1:  # Here domino we are adding is horizontal
            if C1 == D1:   # And domino we need to avoid is in same place
                if j == 0:    # We are first column
                    # print "Found domino to bump out"
                    return [-1]  # In this case, we don't add anything back and just return Tab itself in a list
                                     # We add a -1 to indicate this came from a row bump
                else:
                    # print "Bumping left column domino"
                    r = Col_height(Tab,j-1)  # Find length of previous column
                    Cells = [[r - 2,j-1],[r-1,j-1],[r-2,j-1],[r-1,j-1]]
            if D1[0] == i + 1 and D1[1] == j - 1:
                Cells = [[i,j - 1],[i,j],[i,j-1],[i+1,j - 1]]    # Insertion cells rotate left, Avoidance Cells rotate up
        return Cells  # Here we return the new cells we want
    
    if C1 == D1:
        if C2[0] == i:  # Here domino we are adding is horizontal
            if D2[0] == i + 1 and D2[1] == j:   # And domino we need to avoid is in same place
                # print "Hello Row Turn"
                if i == 0:
                    # print "Found row domino to bump out"
                    return [1]  # In this case, we don't add anything back and just return Tab itself in a list
                                    # We add a 1 to indicate this came from a row bump
                else:
                    # print "Bumping up row domino"
                    r = len(Tab[i - 1])  # Find length of previous row
                    Cells = [[i-1,r-2],[i-1,r-1]]
        if C2[0] == i + 1:  # Here domino we are adding is horizontal
            if D2[0] == i and D2[1] == j + 1:   # And domino we need to avoid is in same place
                if j == 0:    # We are first column
                    # print "Found column domino to bump out"
                    return [-1]  # In this case, we don't add anything back and just return Tab itself in a list
                                     # We add a -1 to indicate this came from a row bump
                else:
                    # print "Bumping left column domino"
                    r = Col_height(Tab,j-1)  # Find length of previous column
                    Cells = [[r-2,j-1],[r-1,j-1]]
        return Cells

    else:   # Otherwise, we just C1, C2 alone, and return them
        return [C1,C2,D1,D2]   # This indiciates our insertion and avoidance cells do not intersect, so we leave them alone
        

def Beta_Unbump(Tab,D1,D2):
    """Returns the tableaux obtained by bumping out cells D1,D2 - assumes they appear on the boundary, and also returns the value of the 
    domino bumped out"""
    Elts = Tab_Elements(Tab)   # Creates a list of the elements of tableaux Tab, sorted
    m = len(Elts)        # Store m as number of elements
    new_Tab = Remove_Domino(Tab,D1,D2) # Initialize new tableaux which we will be updating - we will have to remove 2 spots mentioned
    E1 = D1                            # But we still have information we need
    E2 = D2   # Store these cells as new values which we will update, they represent where bump is, where we must avoid
    # print "E1,E2 =", E1,E2
    for j in reverse(Elts):   # Loop over Tab Elements in decreasing order
        # print "Hello",j
        # print Tab
        [C1,C2] = Tab_Cells(Tab,j)       # List the Cells where j originally was
        # print "C1,C2 =", C1,C2
        # print "Tab Restriction =", Restrict_Tab(Tab,j)
        M = New_Cells(Restrict_Tab(Tab,j),C1,C2,E1,E2)   # Try to add them back to restricted tab, but avoid cells E1,E2
        # print M
        if len(M) == 0:    # This is where there is no interference from cells we need to avoid
            new_Tab = Change_Entry(new_Tab,C1,j)
            new_Tab = Change_Entry(new_Tab,C2,j)    # Here we just change the entries at C1,C2, but do not update E1,E2
        if len(M) == 1:       # In this case we found the element that just came in
            # print "Found the number that just came in"
            return (new_Tab,j,M[0])    # And we return the new_Tab together with the value and the sign
        if len(M) == 4:   # In this case we bump up or right, but nor out
            [I1,I2,E1,E2] = M    # Update the Cells which we must avoid
            # print "E1,E2 =", E1,E2      # And the cells which we insert into
            new_Tab = Change_Entry(new_Tab,I1,j)
            new_Tab = Change_Entry(new_Tab,I2,j)   # Overwriting Entries is bad and transfers to original Tab, but 
            # print "new_Tab is", new_Tab, "at", j                  # we made it so change entry does not do that
    print "Unbumping failed"
    return Remove_Domino(new_Tab,D1,D2)  # If bumping fails somehow, just remove the domino

def Domino_RSK_Inverse(L,R):
    """First, checks all necessary conditions, and then returns the unique word that is sent to P,Q under RSK"""
    if Check_Domino_Tab(L) == True and Check_Domino_Tab(R) == True and shape(L) == shape(R):
        #  Makes sure all requisite conditions are met for inverse to be defined
        return Domino_RSK_Inverse_2(L,R)
        
def Domino_RSK_Inverse_2(L,R):
    """returns the unique word that is sent to P,Q under RSK, assuming requisite conditions met"""
    Elts = Tab_Elements(R)  # Lists the tab elements of R
    new_L = copy(L)   # Create a new copy of L to do stuff with
    ans = []
    for k in reverse(Elts):  # Loops over elements of L in reverse order
        [D1,D2] = Tab_Cells(R,k)  # Stores the cells in which k shows up in R as a domino
        [new_L,x,sign] = Beta_Unbump(new_L,D1,D2)   # Un-bump this domino from L leaving behind a new Tab with x bumped out with sign
                                                    # indicating whether it came from a row or column
        ans = [[x,sign]] + ans   # Add this element to the front of the list
    return ans

def grid_label(Type,Cell):
    """returns the grid label of the Cell using the given Type, where Type is some cell labelled with x"""
    [i,j] = Cell
    [r,s] = Type
    x_diff = i - r
    y_diff = j - s
    if x_diff % 2 == 0:
        if y_diff % 2 == 0:
            return 'X'
        if y_diff % 2 == 1:
            return 'Y'
    if x_diff % 2 == 1:
        if y_diff %2 == 0:
            return 'Z'
        if y_diff % 2 == 1:
            return 'W'
        
def grid_label_B(Cell):
    """returns the grid label of the Cell starting from labelling [0,0] with X"""
    return grid_label([0,0],Cell)
    

def Filled_Corners(Type,Tab):
    """returns the cells of the tableaux that are filled corners"""
    Part = shape(Tab)
    height = len(Part)
    ans = []
    for i in range(0,height):       # Look at the end of each row, and ...
        if grid_label(Type,[i, Part[i] - 1]) == 'X':
            if i < (height - 1) and Part[i] > Part[i + 1]:   # Check if this is actually a corner in grid sense,
                            # And that it does not lie above anything
                ans = ans + [[i, Part[i] - 1]]         # If so, add it to list of Filled Corners
            if i == height - 1:     # The case where we are in the last row, we don't check below 
                ans = ans + [[i,Part[i] - 1]]          
    return ans

def Move_Cell(Tab,Start):
    """Returns the new cells for the given domino off Start position as if moving through a cycle"""
    new_Tab = Tab
    [i,j] = Start    # Let i,j denote the starting cell and k be the value in that cell
    k = Tab[i][j]
    Cells = Tab_Cells(Tab,k)    # Let Cells denot the cells for domino k
    [C1,C2]= Cells
    if C1 == [i,j] and C2 == [i,j + 1]:
        Pivot_Cell = [i-1,j + 2]
        Case = "A"
    if C1 == [i - 1,j] and C2 == [i,j]:
        Pivot_Cell = [i - 2,j + 1]
        Case = "A"
    if C1 == [i,j] and C2 == [i + 1,j]:
        Pivot_Cell = [i + 2, j - 1]
        Case = "B"
    if C1 == [i,j - 1] and C2 == [i,j]:     # Based on where this domino shows up, we know what pivot cell has to be, and 
        Pivot_Cell = [i + 1,j - 2]        # We split this into two cases depending on relative posntion of pivot cell to domino
        Case = "B"
    [a,b] = Pivot_Cell
    if a == -1:
        return [0,b]
    if b == -1:
        return [a,0]   # Take care of edge cases, literally
    if Case == "A":         # Here the pivot cell is in top right of the relevant 3 by 3 square
        if Cell_Present(Tab,Pivot_Cell) == False:  # In this case, we go left of the pivot cell so as to stay a tableaux
            return [a,b - 1]
        if k < Tab[a][b]:     # Left of the pivot cell so as to remain standard
            return [a,b - 1]
        if k > Tab[a][b]:     # Below the pivot cell so as to remain standard
            return [a + 1,b]
    if Case == "B":       # Here the pivot cell is in bottom left of the relevant 3 by 3 square
        if Cell_Present(Tab,Pivot_Cell) == False:   # In this case, we go above of the pivot cell so as to stay a tableaux
            return [a - 1,b]
        if k < Tab[a][b]:    # Above the pivot cell so as to remain standard
            return [a - 1,b]
        if k > Tab[a][b]:    # Right of the pivot cell so as to remain standard
            return [a,b + 1]


def Move_Through_Open_Cycle(Tab,Start):
    """Returns the new Tab one gets after moving the given cycle starting at the filled corner"""
    Active_Cell = Start    # The Active Cell keeps track of where we are and what entry we are changing
    new_Tab = Remove_Cell(Tab,Start)   # We begin a new Tab by removing Start which gets rid of Filled Corner
    while Cell_Present(Tab,Active_Cell) == True:   # Note that we reference original Tab for all info we need
        new_Cell = Move_Cell(Tab,Active_Cell)      # We keep moving through a cycle until we get out of the tableaux
        k = Tab_Entry(Tab,Active_Cell)   # Store k as the entry we are adjusting
        new_Tab = Change_Entry(new_Tab,new_Cell,k)
        Active_Cell = new_Cell
        if Active_Cell == [0,0]:  # W  terminate immediately if we reach the upper left corner of the Tab
            return new_Tab
    return new_Tab

def Specialize(Type,Tab):
    """returns the special tableaux that is equivalent to Tab through open cycles"""
    new_Tab = Tab
    for Start in Filled_Corners(Type,Tab):
        new_Tab = Move_Through_Open_Cycle(new_Tab,Start)   # Simply move through all of the open cycles in the Tab,
                            # Noting that they do not interfer with open another
    return new_Tab

def Specialize_B(Tab):
    """returns the special tableaux that is equivalent to Tab through open cycles, assuming that Type is [0,0]"""
    return Specialize([0,0],Tab)

def Left_Specialize(Type,w):
    """Simply specializes the left Tableaux under Domino RSK of word w"""
    L = Domino_RSK(w)
    return Specialize(Type,L)

def Left_Specialize_B(w):
    """Simply specializes the left Tableaux under Domino RSK of word w"""
    L = Domino_RSK(w)
    return Specialize_B(L)