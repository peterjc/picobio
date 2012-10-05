import numpy as np
from matplotlib import pyplot as plt

filename = "test.cov"
png_filename = "test.png"

def load(filename):
    h = open(filename)
    line = h.readline()
    assert line.startswith(">")
    while line and line[0] == ">":
        name = line[1:].split(None,1)[0]
        values = []
        while line:
            line = h.readline()
            if not line or line[0] == ">":
                break
            values.append([float(v) for v in line.rstrip("\n").split("\t")])
        yield name, np.array(values, np.float)
    h.close()

def stack(data, filename):
    fig = plt.figure()
    total = len(data)
    fig = plt.figure()
    for i, (name, values) in enumerate(data):
        x = range(values.shape[1])
        print i, name, values.shape
        print values.sum(axis=1)
        y_stack = np.cumsum(values, axis=0)
        ax1 = fig.add_subplot(total, 1, i+1)
        ax1.fill_between(x, 0, y_stack[0,:], facecolor="#CC6666", alpha=.7)
        ax1.fill_between(x, y_stack[0,:], y_stack[1,:], facecolor="#1DACD6", alpha=.7)
        ax1.fill_between(x, y_stack[1,:], y_stack[2,:], facecolor="#6E5160")
    plt.show()
    plt.savefig(filename)

data = list(load(filename))
stack(data, png_filename)
