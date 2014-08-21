************
简介
************

阅读本教程的全部内容最多只需要3、4个小时。
你可以阅读HTML或PDF版本，或者在Sage notebook中点击 ``Help`` 
，再点击 ``Tutorial`` ，边阅读边使用 Sage。

虽然Sage主要是用Python实现的，但是不懂Python也可以阅读本教程。
如果你想要学习一下Python（一种非常有趣的语言），
网上有很多关于Python的优秀资源，比如 [PyT]_ 和 [Dive]_ 。
如果你只是想快速的尝试一下Sage，阅读本教程就对了。比如：

::

    sage: 2 + 2
    4
    sage: factor(-2007)
    -1 * 3^2 * 223
    
    sage: A = matrix(4,4, range(16)); A
    [ 0  1  2  3]
    [ 4  5  6  7]
    [ 8  9 10 11]
    [12 13 14 15]
    
    sage: factor(A.charpoly())
    x^2 * (x^2 - 30*x - 80)
    
    sage: m = matrix(ZZ,2, range(4))
    sage: m[0,0] = m[0,0] - 3
    sage: m
    [-3  1]
    [ 2  3]
    
    sage: E = EllipticCurve([1,2,3,4,5]); 
    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 
    over Rational Field
    sage: E.anlist(10)
    [0, 1, 1, 0, -1, -3, 0, -1, -3, -3, -3]
    sage: E.rank()
    1
    
    sage: k = 1/(sqrt(3)*I + 3/4 + sqrt(73)*5/9); k
    36/(20*sqrt(73) + 36*I*sqrt(3) + 27)
    sage: N(k)
    0.165495678130644 - 0.0521492082074256*I
    sage: N(k,30)      # 30 "bits"
    0.16549568 - 0.052149208*I
    sage: latex(k)
    \frac{36}{20 \, \sqrt{73} + 36 i \, \sqrt{3} + 27}

安装
====

如果你没有安装Sage，只是想试几个命令，可以使用在线的Sage notebook：
http://www.sagenb.org 。

要在自己的电脑上安装Sage，请参考Sage主页 [Sage]_ 上的Sage安装指南。
这里我们只强调两点：


#. Sage的安装包是“内置电池”的。也就是说，虽然Sage用到了Python，
   IPython，PARI，GAP，Singular，Maxima，NTL，GMP等等一些软件，
   但是你不需要单独安装这些软件，因为它们已经包含在Sage的发行版里了。
   然而，要使用Sage的一些特定的功能，比如Macaulay或者KASH，
   你必须安装相关的Sage可选包或者已经单独安装了这些软件。Macaulay和KASH
   都是Sage的扩展包（输入 ``sage -optional`` 可以得到可选扩展包列表，
   或者在Sage网站上浏览“下载”页）。

#. 安装编译好的二进制版本（可在Sage网站找到）可能比安装源码版更容易、
   更快。只需要解压缩之后运行 ``sage`` 即可。


使用Sage的方法
==============

使用Sage的方法有好几种。


-  **Notebook图形界面：** 参见参考手册中关于Notebook的章节，
   以及下面的 :ref:`section-notebook` ；

-  **交互命令行：** 参见 :ref:`chapter-interactive_shell` ；

-  **程序：** 在Sage中编写解释型或者编译型的程序（参见
   :ref:`section-loadattach` 和 :ref:`section-compile` ）；

-  **脚本：** 在独立的Python脚本中调用Sage库文件（参见
   :ref:`section-standalone` ）。


Sage的长期目标
==============

-  **实用：** Sage的预期用户是学数学的学生（从高中生到研究生）、
   教师以及研究人员。我们的目标是在代数、几何、数论、微积分、
   数值计算等领域提供可用于探索和尝试的软件。
   Sage使得进行与数学对象有关的交互实验变得容易。

-  **高效：** 越快越好。Sage使用高度优化的成熟软件，如GMP，
   PARI，GAP和NTL。这样，Sage的某些运算非常快。

-  **免费、开源：** 源代码必须可以自由的获取，并且有较好的可读性，
   这样用户才能真正了解系统是如何运行的，并且更容易进行扩展。
   就像数学家要深入理解一个定理的话，就要仔细地阅读定理的证明，
   最起码要浏览一下。搞计算的人应该可以通过阅读源码来了解计算是如何进行的。
   如果你在论文中使用Sage进行计算，你可以确保读者能够免费得到Sage及其源码。
   并且你可以打包或者重新发布你所使用的Sage版本。

-  **易于编译：** Linux，OS X和Windows的用户应该很容易使用源代码编译Sage。
   这为用户修改系统提供了便利。

-  **协作：** 为其他计算机代数系统提供健壮的接口，包括PARI，GAP，
   Singular，Maxima，KASH，Magma，Maple和Mathematica。
   Sage希望统一并扩展现有的数学软件。

-  **文档完善：** 教程，编程指南，参考手册和基本指南要包含大量的例子，
   以及对数学背景的讨论。

-  **可扩展：** 可以定义新的数据类型或者从内置的类型中继承，
   可以使用其他语言编写的代码。

-  **用户友好：** 给定对象所提供的功能应该是清晰易懂的，
   文档和源码应该易于查看。用户支持要达到比较高的水平。

.. [Dive] Dive into Python, Freely available online at 
          http://diveintopython.org

.. [PyT] The Python Tutorial, http://www.python.org/

.. [Sage] Sage, http://www.sagemath.org
