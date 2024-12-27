
def derivative(x, y):
    """Computes the derivative of a function from x and y values.

    Parameters:
        x (array): xvalues.
        y (array): yvalues.

    Returns:
        derivative (array): the derivative of the function from x to y
    """
    if len(x) != len(y):
        raise ValueError("x and y shall have the same length")

    return (y[1:] - y[:-1]) / (x[1:] - x[:-1])

