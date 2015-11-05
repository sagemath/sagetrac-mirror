Test combining various modifiers::

    sage: sys.maxint  # long time, abs tol 0.001
    2147483646.999            # 32-bit
    9223372036854775806.999   # 64-bit

Testing test limits::

    sage: 1 # long time: 10s
    1
    sage: 2 # long time: 20s
    2
    sage: 3 # long time: 1m
    3
    sage: 4 # long time: 4h
    4
