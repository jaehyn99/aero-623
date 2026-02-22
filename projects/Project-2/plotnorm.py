import numpy as np
import matplotlib.pyplot as plt
import sys

def main():
    if (len(sys.argv) < 2):
        print('Pass at least one argument: Norm')
    
    argc = len(sys.argv)
    plt.figure()
    legend = [] if argc <= 2 else sys.argv[1:]
    print(legend)
    for i in range(1, argc):
        norm = np.loadtxt(sys.argv[i])
        iter = np.arange(norm.size)
        plt.semilogy(iter, norm)
    if argc > 2:
        plt.legend(["RK2", "RK3"])
    plt.xlabel("Iteration number")
    plt.ylabel("L1 residual norm")
    plt.show()

if __name__ == "__main__":
    main()
