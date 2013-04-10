
def bhthreshold(histo):
    import numpy

    start = 0
    end = len(histo)
    mid = int((start + end) / 2.0)

    wl = numpy.sum(histo[start:mid])
    wr = numpy.sum(histo[mid:end])

    while (start <= end):
        print start, end
        if wr > wl:
            end -= 1
            wr -= histo[end]
            if ((start + end) / 2) < mid:
                wr += histo[mid]
                mid -= 1
                wl -= histo[mid]

        else:
            start += 1
            wl -= histo[start]
            if ((start + end) / 2) > mid:
                wl += histo[mid + 1]
                wr -= histo[mid + 1]
                mid += 1

    return mid

def maximum_deviation(histo):
    import numpy
    from math import sqrt

    x0 = numpy.argmax(histo) + 0.5
    y0 = numpy.max(histo)
    x1 = len(histo) - 0.5
    y1 = histo[len(histo)-1]
    m = (y1 - y0) / (x1 - x0)
    c = y0 - m * x0

    dmax = 0
    imax = -1
    for i, y in enumerate(histo):
        x = i + 0.5
        d = abs(m * x - y + c) / sqrt(m**2 + 1)
        if d > dmax:
            dmax = d
            imax = i

    return imax


def normal_threshold(histo, nsdev=10):
    import numpy
    mean = numpy.mean(histo)
    sdev = numpy.std(histo)
    return mean + nsdev * sdev


def percentage_threshold(histo, percent=0.9999):

    cumulative = histo.copy()
    for i, h in enumerate(histo):
        if i > 0:
            cumulative[i] += cumulative[i-1]

    for i, h in enumerate(cumulative):
        if h > percent:
            return i

    raise ValueError("Bad value")
