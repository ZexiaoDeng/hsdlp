# 齐次自对偶嵌入算法



## 测试算例

测试算例下面给出几个算例进行验证

$$
\begin{equation}
    \begin{array}{crrrrrrc}
    \textrm{maximize} & x_1 & - 2x_2 & + x_3 & & &\\
    \textrm{subject to} & x_1 & + x_2 & - 2x_3 & + x_4 & = & 0, \\
    & 2x_1 & - x_2 & + 4x_3 & & \leq & 0, \\
    & -x_1 & + 2x_2 & - 4x_3 &  & \leq & 0, \\
    & x_j & \geq 0, & j = 1,& \cdots, & 4. &
    \end{array}
\end{equation}
$$



上式LP问题的最优解为$ \left(x_1, x_2, x_3, x_4\right) = \left(0, 12, 5, 8\right) $，目标函数值为$ f_{min} = -19 $。

## 参考文献

[1] Robert J. Vanderbei. Linear Programming Foundations and Extensions. Springer New York Heidelberg

Dordrecht London, Princeton, New Jersey, USA, 4nd edition, 2014.
