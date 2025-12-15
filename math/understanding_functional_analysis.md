# Understanding the Adjoint Operator #

Let $X$ and $Y$ be normed spaces, and $T: X\to Y$ a bounded linear transformation.
The canonical problem of interest lies in the equation
$$
    Ax = y,
$$
where $y$ is a datum and $x$ is unknown.
Of course, we would like $A^{-1}$ to exist so that $x = A^{-1} y$.
Instead of trying to solve the _strong_ problem $Ax = y$,
we look for solutions of the _weak_ problem
$$
    \ell(Ax) = \ell(y)
    \qquad \ell\in Y^{*}.
$$
It's convenient to think of $X^{*}$ and $Y^{*}$ as questions or
linear statistics.
In this sense, weak solutions are ''essentially'' the same
as strong solution, since they respond the same to the same questions.

So, how do we recover $x$? Firstly, we need to change the questions
in $Y^{*}$ to question in $X^{*}$.
The adjoint $A^{*}\colon Y^{*} \to X^{*}$ is defined to do exactly that
by the following equation
$$
    \ell(Ax) 
    = (A^{*}\ell)(x)
    \qquad \forall \ell\in Y^{*}.
$$
If $A^{*}$ happen to be invertible, we get that
