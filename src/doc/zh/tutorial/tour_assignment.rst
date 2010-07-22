
赋值，等式和算术运算
====================================

除了少数例外，Sage使用了Python语言，
因此大多数关于Python的入门书籍将帮助你学习Sage。

Sage使用 ``=`` 进行赋值，使用 ``==``, ``<=``, ``>=``, ``<`` 和 ``>`` 
进行比较：

::

    sage: a = 5
    sage: a
    5
    sage: 2 == 2
    True
    sage: 2 == 3
    False
    sage: 2 < 3
    True
    sage: a == 5
    True

Sage提供了基本的数学运算：

::

    sage: 2**3    #  ** 的意思是指数
    8
    sage: 2^3     #  ^ 和 ** 一个意思（这点不象Python）
    8
    sage: 10 % 3  #  对于整数参数，% 的意思是取模，即求余数
    1
    sage: 10/4
    5/2
    sage: 10//4   #  对于整数参数，// 返回整数商
    2
    sage: 4 * (10 // 4) + 10 % 4 == 10
    True
    sage: 3^2*4 + 2%5
    38

象 ``3^2*4 + 2%5`` 这样的表达式的计算结果取决于运算的顺序。
计算顺序由“运算符优先级表”指定，参见 :ref:`section-precedence` 。

Sage还提供了许多常见的数学函数，这里是几个例子：

::

    sage: sqrt(3.4)
    1.84390889145858 
    sage: sin(5.135)
    -0.912021158525540 
    sage: sin(pi/3)
    1/2*sqrt(3)

象最后一个例子那样，有些数学表达式返回“精确”的值，而不是近似的数值结果。
要得到一个近似的数值解，使用函数 ``n`` 或者方法 ``n``
（两者的全名都是 ``numerical_approx``, 并且函数 ``N`` 和 ``n`` 是一样的）。
它们都有可选参数 ``prec`` 和 ``digits`` ，前者指定结果的二进制位数，
即bit数，后者指定结果的十进制位数。默认精度是53 bit。

::

    sage: exp(2)
    e^2
    sage: n(exp(2))
    7.38905609893065
    sage: sqrt(pi).numerical_approx()
    1.77245385090552
    sage: sin(10).n(digits=5)
    -0.54402
    sage: N(sin(10),digits=10)
    -0.5440211109 
    sage: numerical_approx(pi, prec=200)
    3.1415926535897932384626433832795028841971693993751058209749

Python是动态类型的，因此每一个赋给变量的值都有一个类型，
但是在给定作用域内，一个给定的变量可以接受任何Python类型的值。

::

    sage: a = 5   # a是整数
    sage: type(a)
    <type 'sage.rings.integer.Integer'>
    sage: a = 5/3  # 现在a是有理数
    sage: type(a)
    <type 'sage.rings.rational.Rational'>
    sage: a = 'hello'  # 现在a是字符串
    sage: type(a)
    <type 'str'>

C语言就很不一样了，它是静态类型的。一个被声明为整数的变量，
在它的作用域内只能接受整数值。

Python中一个潜在的容易混淆的地方是，以0开头的整数是八进制数，
也就是以8为基的数。

::

    sage: 011
    9
    sage: 8 + 1
    9
    sage: n = 011
    sage: n.str(8)   # 将n以8进制字符串形式输出
    '11'

这与C语言的规定是一致的。
