def TDMAsolver(a, b, c, d):
    """
    Solve a tridiagonal system of equations AX = D using the Thomas algorithm.

    Parameters:
    a : list of floats
        Sub-diagonal elements of the matrix (a[0] is unused -> Must add a 0 at the begining of the list).
    b : list of floats
        Diagonal elements of the matrix.
    c : list of floats
        Super-diagonal elements of the matrix (c[-1] is unused -> Must add a 0 at the end of the list).
    d : list of floats
        Right-hand side of the equation.

    Returns:
    x : list of floats
        Solution to the system.
    """
    n = len(d)

    # Make copies of the input arrays to avoid modifying them
    cp = c[:]
    dp = d[:]

    # Forward elimination
    for i in range(1, n):
        w = a[i] / b[i - 1]
        b[i] -= w * cp[i - 1]
        dp[i] -= w * dp[i - 1]

    # Back substitution
    x = [0] * n
    x[-1] = dp[-1] / b[-1]
    for i in range(n - 2, -1, -1):
        x[i] = (dp[i] - cp[i] * x[i + 1]) / b[i]

    return x
