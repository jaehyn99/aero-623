import numpy as np
import matplotlib.pyplot as plt
import sys

def main():
    if (len(sys.argv) < 2):
        print('Pass at least one argument: Norm')
    norm = np.loadtxt(sys.argv[1])
    iter = np.arange(norm.size)
    plt.figure()
    plt.semilogy(iter, norm)
    plt.xlabel("Iteration number")
    plt.ylabel("L1 residual norm")
    plt.show()

if __name__ == "__main__":
    main()
