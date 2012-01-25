"""
Qhull tests

TESTS:
    sage: import sage.tests.benchmark
    sage: testsys = [[exp(t/10)*sin(t),exp(t/10)*cos(t)] for t in srange(0,6*pi,.25)]
    sage: input_sys = '2\n'+str(len(testsys))+'\n'+join([str(q[0]) + ' ' + str(q[1]) for q in testsys],sep='\n')
    sage: from sage.misc.misc import tmp_filename
    sage: tname = tmp_filename()
    sage: f = file(tname,'w')
    sage: f.write(input_sys)
    sage: f.close()
    sage: import subprocess as S
    sage: out = S.check_output('cat %s | qconvex i'%tname, shell=True)
    sage: print(out)
          22
          73 72 
          75 74 
          74 73 
          57 56 
          56 55 
          54 75 
          55 54 
          70 69 
          67 66 
          64 63 
          71 70 
          72 71 
          68 67 
          69 68 
          65 64 
          66 65 
          58 57 
          59 58 
          61 60 
          60 59 
          62 61 
          63 62
"""

