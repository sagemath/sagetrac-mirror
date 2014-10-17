class BlockCodeTools(object):
	#Parent class
	#It's a toolbox with useful functions, used in many different codes
    #Some values are
    #   F field of code's elements
    #   n length of the code
    #   k dimension of the code
    #   d minimum distance of the code


    def __contains__(self, x):
        return self.iscodeword(x)

    def __eq__(self, other):
        return self.generator_matrix() == other.generator_matrix()

    def generator_matrix(self):
        return "TODO gen mat"

    def parity_check_matrix(self):
        Pc=self.generator_matrix().right_kernel()
        return Pc.basis_matrix() 

    def encode(self, word):
        return vector(word)*self.generator_matrix()
    
    def dimension(self):
        return self.generator_matrix().rank()

    def random_codeword(self):
        k = self.dimension()
        return random_vector(self.F, k)*self_generator_matrix()

    def minimum_distance(self):
        return "TODO min dist"

    def syndrome(self, r):
        return self.parity_check_matrix()*r.column()
    
