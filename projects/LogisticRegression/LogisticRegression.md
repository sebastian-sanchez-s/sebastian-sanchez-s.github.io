Suppose we want to label an object in any of $N$ categories.
The relevant information of the object is summarized in
a vector with $d$ entries.
Then, one can state the problem as:
find a function $f: R^d \to R^N$ such that the $i$-entry of $f(X)$ is
the probability to label $X$ in the $i$ category.

In order to find such a function, we need samples $(X^k, Y^k)$,
$k=1,\ldots, K$, to learn from.
For now, the following conditions arise:

1. $f(X) \ge 0$ and $f(X)\cdot 1 = 1$.
2. $f(X^k) = p^k$ should peak only in the coordinate $\textup{argmax}_{i} Y^{k}_i$.
This is handled by the information theoretic entropy, acting on
$p=f(X)$; 
maximizing this formula means that we want all probabilities to
be as close to zero unless it is necessary.
$$
    H(p) = \sum_{i=1}^{N} p_i \log(1/p_i).
$$

3. The probability $p^k$ should match the data. Pure interpolation,
$p^k_i = Y^k_i$, does not make much sense since samples are
random.
Instead, we should ask for something more structural, such 
as matching the expected value of each coordinate $i$
between the prediction and the empirical data.
$$
    \sum_{k} p^k_i X^k_i = \sum_{k} Y^k_i X^k_i
    \quad
    (i=1:N).
$$

This is enough to propose an optimization problem:

$$
\left\{
\begin{aligned}
    \max \quad &\sum_{k=1}^{K} \sum_{i=1}^{N} p^k_i \log(1/p^k_i)
    \\
    \textup{s.t.}\quad
    & p^k_i \ge 0 &&\quad (i=1:N, k=1:K)\\
    & \sum_{i} p^k_i = 1 &&\quad (k=1:K)\\
    & \sum_{k} p^k_i X^k_i = \sum_{k} Y^k_i X^k_i &&\quad (i=1:N)
\end{aligned}
\right.
$$

The Lagrangian reads
$$
    L(p, a, b) = 
    -\sum_{k=1}^{K} \sum_{i=1}^{N} p^k_i \log(p^k_i)
    +
    \sum_{k=1}^{K} a_k (\sum_{i=1}^{N} p^k_i - 1)
    +
    \sum_{i=1}^{N} b_i (\sum_{k=1}^{N} p^k_i X^k_i - Y^k_i X^k_i),
$$
and the first order conditions
$$
    0 
    = \partial_{p^k_i} L
    = -\log(p^k_i) - 1 + a_k + b_i X^k_i
$$
implies that
$$
    p^k_i = \exp(b_i X^K_i + a_k - 1).
$$
Next,
$$
    1 
    = \sum_{i=1}^{N} p^k_i 
    = \sum_{i=1}^{N} \exp(b_i X^k_i + a_k - 1)
    = e^{a_k-1} \sum_{i=1}^{N} \exp(b_i X^k_i).
$$
Hence,
$$
    p^k_i
    = \exp(b_i X^k_i + a_k - 1)
    = \frac{\exp(b_i X^k_i)}{\sum_{j=1}^{N} \exp(b_j X^k_j)}.
$$

This gives us a nice ansatz for $f$. All we need to do is set
$$
    f_b(X)
    = \frac{e^{b_i X_i}}{\sum_{j=1}^N e^{b_j X_j}}
    = \frac{1}{1 + \sum_{j\ne i} e^{b_j X_j}}
$$
and find the coefficients that match $p^k_i$.
