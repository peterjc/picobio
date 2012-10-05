import numpy as np
from matplotlib import pyplot as plt

filename = "test.cov"
png_filename = "test.png"

def load(filename):
    h = open(filename)
    line = h.readline()
    assert line.startswith(">")
    name = line[1:].split(None,1)[0]
    values = []
    while line:
        line = h.readline()
        if not line or line[0] == ">":
            break
        values.append([float(v) for v in line.rstrip("\n").split("\t")])
    h.close()
    return name, np.array(values, np.float)

def stack(name, values, filename):
    x = range(values.shape[1])
    print values.shape
    print values.sum(axis=1)
    #Based on example here http://stackoverflow.com/questions/2225995/how-can-i-create-stacked-line-graph-with-matplotlib
    y_stack = np.cumsum(values, axis=0)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.fill_between(x, 0, y_stack[0,:], facecolor="#CC6666", alpha=.7)
    ax1.fill_between(x, y_stack[0,:], y_stack[1,:], facecolor="#1DACD6", alpha=.7)
    ax1.fill_between(x, y_stack[1,:], y_stack[2,:], facecolor="#6E5160")
    plt.show()
    plt.savefig(filename)

name, values = load(filename)
stack(name, values, png_filename)

