import datetime as dtm

def TDMAsolver(a, b, c, d):

    """Solve a tridiagonal system of equations AX = D using the Thomas algorithm.

    Parameters:
    a : list of floats
        Sub-diagonal elements of the matrix (a[0] is unused -> Must add a 0 at the beginning of the list).
    b : list of floats
        Diagonal elements of the matrix.
    c : list of floats
        Super-diagonal elements of the matrix (c[-1] is not used in calculations, but the list should have the same length as b and d).
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

def runtime_program (beginning_date_and_time):

    """Calculates and prints program runtime based on start time.

    Parameters
    ----------
    beginning_date_and_time : datetime.datetime
    
    Returns
    -------
    None
    """

    ending_date_and_time=dtm.datetime.now()
    runtime = ending_date_and_time - beginning_date_and_time

    #type(runtime) =  datetime.timedelta, it should be convert in a string to be easily manipulated and printed
    string_runtime = str(runtime)

    hour = string_runtime.split(":")[0]
    minute = string_runtime.split(":")[1]
    second = str(round(float(string_runtime.split(":")[2]),2))

    # These different cases are made to avoid printing 00h00min00,4s for the shortest program (just for the visual aspect)
    if hour == "0":

        if minute == "00":

            print("Runtime : "+second+"s")

        else:

            print("Runtime : "+minute+"min"+second+"s")
    
    else:

        print("Runtime : "+hour+"h"+minute+"min"+second+"s")
    
    return None

def norm_L2 (matrix):

    "Calculate the L2 norm of a matrix"

    Nx = len(matrix)
    Ny = len(matrix[0])

    norm = 0

    for i in range(len(matrix)):
        for j in range(len(matrix[0])):

            norm += matrix[i][j]*matrix[i][j]
    
    norm = norm/(Nx*Ny)

    return norm