

def func(a):
    return a * a

from multiprocessing import Pool
pool = Pool(processes=4)



print pool.map(func, [1, 2, 3, 4, 1, 2, 3, 4])
