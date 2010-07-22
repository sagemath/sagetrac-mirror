.. _chapter-help:

获取帮助
============

Sage拥有强大的内置文档，只需要输入函数或者常数的名字，
再加个问号即可：

.. skip

::

    sage: tan?
    Type:        <class 'sage.calculus.calculus.Function_tan'>
    Definition:  tan( [noargspec] )
    Docstring: 
    
        The tangent function
    
        EXAMPLES:
            sage: tan(pi)
            0
            sage: tan(3.1415)
            -0.0000926535900581913
            sage: tan(3.1415/4)
            0.999953674278156
            sage: tan(pi/4)
            1
            sage: tan(1/2)
            tan(1/2)
            sage: RR(tan(1/2))
            0.546302489843790
    sage: log2?
    Type:        <class 'sage.functions.constants.Log2'>
    Definition:  log2( [noargspec] )
    Docstring: 
    
        The natural logarithm of the real number 2.
        
        EXAMPLES:
            sage: log2
            log2
            sage: float(log2)
            0.69314718055994529
            sage: RR(log2)
            0.693147180559945
            sage: R = RealField(200); R
            Real Field with 200 bits of precision
            sage: R(log2)
            0.69314718055994530941723212145817656807550013436025525412068
            sage: l = (1-log2)/(1+log2); l
            (1 - log(2))/(log(2) + 1)
            sage: R(l)
            0.18123221829928249948761381864650311423330609774776013488056
            sage: maxima(log2)
            log(2)
            sage: maxima(log2).float()
            .6931471805599453
            sage: gp(log2)
            0.6931471805599453094172321215             # 32-bit
            0.69314718055994530941723212145817656807   # 64-bit
    sage: sudoku?
    File:        sage/local/lib/python2.5/site-packages/sage/games/sudoku.py
    Type:        <type 'function'>
    Definition:  sudoku(A)
    Docstring: 
    
        Solve the 9x9 Sudoku puzzle defined by the matrix A.
    
        EXAMPLE:
            sage: A = matrix(ZZ,9,[5,0,0, 0,8,0, 0,4,9, 0,0,0, 5,0,0,
        0,3,0, 0,6,7, 3,0,0, 0,0,1, 1,5,0, 0,0,0, 0,0,0, 0,0,0, 2,0,8, 0,0,0,
        0,0,0, 0,0,0, 0,1,8, 7,0,0, 0,0,4, 1,5,0,   0,3,0, 0,0,2,
        0,0,0, 4,9,0, 0,5,0, 0,0,3])
            sage: A
            [5 0 0 0 8 0 0 4 9]
            [0 0 0 5 0 0 0 3 0]
            [0 6 7 3 0 0 0 0 1]
            [1 5 0 0 0 0 0 0 0]
            [0 0 0 2 0 8 0 0 0]
            [0 0 0 0 0 0 0 1 8]
            [7 0 0 0 0 4 1 5 0]
            [0 3 0 0 0 2 0 0 0]
            [4 9 0 0 5 0 0 0 3]
            sage: sudoku(A)
            [5 1 3 6 8 7 2 4 9]
            [8 4 9 5 2 1 6 3 7]
            [2 6 7 3 4 9 5 8 1]
            [1 5 8 4 6 3 9 7 2]
            [9 7 4 2 1 8 3 6 5]
            [3 2 6 7 9 5 4 1 8]
            [7 8 2 9 3 4 1 5 6]
            [6 3 5 1 7 2 8 9 4]
            [4 9 1 8 5 6 7 2 3]

Sage还提供了“Tab补全”功能：输入函数名的前面几个字母，
然后按tab键。比如，如果你输入 ``ta`` 再按 ``TAB``,
Sage就会列出 ``tachyon, tan, tanh, taylor``.
这是一个查找Sage函数或者其他结构的好方法。


.. _section-functions:

函数，缩进和计数
====================================

在Sage中定义一个新的函数，需要使用 ``def`` 命令，
并且在变量列表后跟一个冒号。比如：

::

    sage: def is_even(n):
    ...       return n%2 == 0
    ...
    sage: is_even(2)
    True
    sage: is_even(3)
    False

注：根据你所阅读的本教程的版本的不同，在这个例子中，
你可能会看到第二行有三个点"``...``" 。不要输入它们，
它们只是强调一下代码是缩进的。不管是什么情况，在程序块的最后，
按 [Return/Enter] 插入一个空行以结束函数的定义。

你没有指定输入参数的类型。你可以指定多个输入，
每个参数都可以带一个可选的默认值。比如下面的函数中，
如果不指定 ``divisor`` 的话，默认取 ``divisor=2``. 

::

    sage: def is_divisible_by(number, divisor=2):
    ...       return number%divisor == 0
    sage: is_divisible_by(6,2)
    True
    sage: is_divisible_by(6)
    True
    sage: is_divisible_by(6, 5)
    False

在调用函数时，你还可以明确的指定一个或多个参数的值。
如果你明确指定参数的值，参数可以以任何顺序出现。

.. link

::

    sage: is_divisible_by(6, divisor=5)
    False
    sage: is_divisible_by(divisor=2, number=6)
    True

与其他很多语言不同，Python中的程序块不用花括号或者begin，end来标记，
而是用精确的缩进来标记。比如下面的代码有一个语法错误, ``return``
语句与它上面的语句缩进的不完全一致。

.. skip

::

    sage: def even(n):
    ...       v = []
    ...       for i in range(3,n):
    ...           if i % 2 == 0:
    ...               v.append(i)
    ...      return v
    Syntax Error:
           return v

修正缩进格数之后，函数就对了：

::

    sage: def even(n):
    ...       v = []
    ...       for i in range(3,n):
    ...           if i % 2 == 0:
    ...               v.append(i)
    ...       return v
    sage: even(10)
    [4, 6, 8]

多数情况下，一行结束后会开始一个新行，这时行尾不需要分号。
但是如果要将多个语句放在同一行，就要用分号隔开：

::

    sage: a = 5; b = a + 3; c = b^2; c
    64

如果你要将一行代码分开放在多行，要在行尾使用反斜杠：

::

    sage: 2 + \
    ...      3
    5

在Sage中，通过遍历一个范围内的整数进行计数。
比如下面代码中的第一行相当于C++或者Java中的 ``for(i=0; i<3; i++)``:

::

    sage: for i in range(3):
    ...       print i
    0
    1
    2

下面的第一行相当于 ``for(i=2;i<5;i++)``.

::

    sage: for i in range(2,5):
    ...       print i
    2
    3
    4

第三个参数控制步长，下面的第一行相当于
``for(i=1;i<6;i+=2)``.

::

    sage: for i in range(1,6,2):
    ...       print i
    1
    3
    5

可能你经常需要将Sage中的计算结果以漂亮的表格形式输出，
一个简单的方法是使用格式化字符串。下面，我们计算数的平方和立方，
并建立一个有三列的表格，每一列都是6个字符宽。

::

    sage: for i in range(5):
    ...       print '%6s %6s %6s'%(i, i^2, i^3)
         0      0      0
         1      1      1
         2      4      8
         3      9     27
         4     16     64

Sage中最最基本的数据结构是list，跟字面意思一样，list就是任意对象的列表。
比如我们刚才用到的 ``range`` 命令就产生一个list：

::

    sage: range(2,10)
    [2, 3, 4, 5, 6, 7, 8, 9]

下面是更复杂的list：

::

    sage: v = [1, "hello", 2/3, sin(x^3)]
    sage: v
    [1, 'hello', 2/3, sin(x^3)]

象其他很多语言一样，list的下标以0开始计数。

.. link

::

    sage: v[0]
    1
    sage: v[3]
    sin(x^3)

使用 ``len(v)`` 得到 ``v`` 的长度，使用 ``v.append(obj)``
向 ``v`` 的末尾添加新的对象，使用 ``del v[i]``
删除 ``v`` 的第 :math:`i` 个元素：

.. link

::

    sage: len(v)
    4
    sage: v.append(1.5)
    sage: v
    [1, 'hello', 2/3, sin(x^3), 1.50000000000000]
    sage: del v[1]
    sage: v
    [1, 2/3, sin(x^3), 1.50000000000000]

另一个重要的数据结构是dictionary（或associative array）。
用法和list类似，但它几乎可以使用所有的对象进行索引（指标必须是固定的）：

::

    sage: d = {'hi':-2,  3/8:pi,   e:pi}
    sage: d['hi']
    -2
    sage: d[e]
    pi

你可以使用“类”定义新的数据结构。将数学对象用类进行封装是一个强大的技术，
可以帮你简化和组织Sage程序。下面我们定义一个类来表示不超过 *n*
的正偶数列表，它由内置类型 ``list`` 继承而来。

::

    sage: class Evens(list):
    ...       def __init__(self, n):
    ...           self.n = n
    ...           list.__init__(self, range(2, n+1, 2))
    ...       def __repr__(self):
    ...           return "Even positive numbers up to n."

在建立对象时，调用 ``__init__`` 方法进行初始化；
``__repr__`` 方法打印对象。我们在 ``__init__``
方法的第二行调用list的constructor方法。
我们可以象下面一样建立 ``Evens`` 类的一个对象：
（译注：原文中使用的变量名为 ``e``, 考虑到 ``e`` 是内置的常数，
因此换成了 ``ee`` ）

.. link

::

    sage: ee = Evens(10)
    sage: ee
    Even positive numbers up to n.

注意 ``ee`` 使用我们定义的 ``__repr__`` 方法进行输出。
要使用 ``list`` 的函数才能查看隐含的数据列表：

.. link

::

    sage: list(ee)
    [2, 4, 6, 8, 10]

我们还可以访问 ``n`` 属性或者将 ``ee`` 当做list。

.. link

::

    sage: ee.n
    10
    sage: ee[2]
    6
