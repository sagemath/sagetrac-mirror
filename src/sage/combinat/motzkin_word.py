from sage.structure.list_clone import ClonableArray

class MotzkinWord(ClonableArray):
    def __init__(self, parent, word):
        for i in range(1,len(word)+1):
            if sum(word[1:i])<0:
                raise ValueError("A Motzkinword is not allowed to go beneath the starting level, error at position {}".format(i))
            if word[i] notin {1,-1,0}:
                raise ValueError("Only steps in {1,-1,0} are allowed, however you used {}.".format(word[i]))
        ClonableArray.__init__(self, parent,word)
    def __repr__(self):
        return "MotzkinWord {}".format(self.word)

