class modsym_dist_aws(modsym_dist):
	"""this special class arises as specializes of families of OMSs.  The only difference is the filtration is slightly modified to fit in families."""

	def normalize(self):
		"""normalizes all values of self via the AWS filtration"""
		assert self.valuation()>=0, "moments are not integral"
		v=[]
		for j in range(0,len(self.data)):
			v=v+[self.data[j].normalize_aws()]
		return modsym_dist_aws(self.level,v,self.manin)

