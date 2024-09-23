import subprocess
import sys
import numpy as np
import cmath
import matplotlib.pyplot as plt


# Function to install a package if it is not present
def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])


# Attempt to import required packages, and if they are not installed, install them
def ensure_dependencies():
    try:
        import numpy
    except ImportError:
        print("NumPy is not installed. Installing...")
        install('numpy')
        import numpy

    try:
        import cmath
    except ImportError:
        raise ImportError("cmath is part of the Python standard library and should be available.")

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib is not installed. Installing...")
        install('matplotlib')
        import matplotlib.pyplot as plt

    try:
        global mplcursors
        import mplcursors
    except ImportError:
        print("mplcursors is not installed. Installing...")
        install('mplcursors')
        import mplcursors


# Call to ensure all modules are installed
ensure_dependencies()


# Function to calculate the parabola using the Muller method
def P(x, a, b, c):
    return a * x ** 2 + b * x + c


# Function to perform interpolation using the Muller method
def muller_interpolation(f, x1, x2, x3, epsilon=1e-7):
    i = 0
    first_iteration = True

    print("n\txn\t\tf(xn)")
    print(f"1\t{x1}\t\t{f(x1)}")
    print(f"2\t{x2}\t\t{f(x2)}")
    print(f"3\t{x3}\t\t{f(x3)}")

    while abs(f(x3)) > epsilon:
        q = (x3 - x2) / (x2 - x1)
        a = q * f(x3) - q * (1 + q) * f(x2) + q ** 2 * f(x1)
        b = (2 * q + 1) * f(x3) - (1 + q) ** 2 * f(x2) + q ** 2 * f(x1)
        c = (1 + q) * f(x3)

        if first_iteration:
            plt.figure()
            y = np.linspace(min(x1, x2, x3) - 1, max(x1, x2, x3) + 1, 100)
            plt.plot(y, P(y, a, b, c), label='Interpolated Parabola')
            plt.scatter([x1, x2, x3], [f(x1), f(x2), f(x3)], color='red', label='Data Points')
            plt.title("Parabola - First Iteration")
            plt.legend()
            mplcursors.cursor(hover=True)  # Adding interactivity
            plt.show()
            first_iteration = False

        r = x3 - (x3 - x2) * ((2 * c) / (b + cmath.sqrt(b ** 2 - 4 * a * c)))
        s = x3 - (x3 - x2) * ((2 * c) / (b - cmath.sqrt(b ** 2 - 4 * a * c)))

        xplus = r if abs(f(r)) < abs(f(s)) else s
        xplus = xplus.real if xplus.imag == 0 else xplus

        print(f"{i + 4}\t{round(xplus, 5)}\t\t{round(f(xplus), 5)}")

        x1, x2, x3 = x2, x3, xplus
        i += 1

    print(f"\nFunction obtained after interpolation: P(x) = {a} * x^2 + {b} * x + {c}")
    return lambda x: a * x ** 2 + b * x + c  # Returning the constructed parabolic function


# Function to perform Richardson extrapolation for the derivative
def richardson(f, x, n, h):
    d = np.zeros((n + 1, n + 1), float)

    for i in range(n + 1):
        d[i, 0] = 0.5 * (f(x + h) - f(x - h)) / h
        powerOf4 = 1
        for j in range(1, i + 1):
            powerOf4 = 4 * powerOf4
            d[i, j] = d[i, j - 1] + (d[i, j - 1] - d[i - 1, j - 1]) / (powerOf4 - 1)
        h = 0.5 * h

    return d[n, n]  # Most accurate value for the derivative


# Main function for user selection
def main():
    print("Select an option:")
    print("1: Calculate the derivative of a function at a point using Richardson")
    print("2: Perform interpolation using Muller and then calculate the derivative at a point using Richardson")

    choice = input("Enter your choice (1 or 2): ")

    if choice == "1":
        func_str = input("Enter the function (e.g., x**3): ")
        f = eval(f"lambda x: {func_str}")
        x = float(input("Enter the point at which to calculate the derivative: "))
        h = float(input("Enter the initial step size (h): "))
        n = int(input("Enter the number of extrapolation steps: "))
        result = richardson(f, x, n, h)
        print(f"The estimated derivative at {x} is: {result}")

    elif choice == "2":
        x1 = float(input("Enter the point x1: "))
        x2 = float(input("Enter the point x2: "))
        x3 = float(input("Enter the point x3: "))

        # Build the parabolic function using Muller interpolation
        interpolated_func = muller_interpolation(lambda x: x ** 2, x1, x2, x3)

        print("Interpolation performed using Muller, now calculating the derivative of the constructed function.")

        x = float(input("Enter the point at which to calculate the derivative: "))
        h = float(input("Enter the initial step size (h) for calculating the derivative: "))
        n = int(input("Enter the number of extrapolation steps: "))
        result = richardson(interpolated_func, x, n, h)
        print(f"The estimated derivative at {x} is: {result}")

    else:
        print("Invalid choice.")


# Run the program
if __name__ == "__main__":
    main()
