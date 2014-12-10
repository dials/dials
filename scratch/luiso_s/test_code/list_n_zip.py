def getKey(item):
  return item[3]


if( __name__ == "__main__" ):
  a= [
      ["h", "i", "j", "k", "l", "m"],
      ["u", "v", "w", "x", "y", "z"],
      ["a", "b", "c", "d", "e", "f"],
      [ 1 ,  5 ,  3 ,  4 ,  2 ,  6 ]
                                    ]

  print "a =", a
  print "zip(a) =", zip(a)
  print "zip(*a) =", zip(*a)
  '''
  print "tuple(zip(*a)) =", tuple(zip(*a))
  print "list(tuple(zip(*a))) =", list(tuple(zip(*a)))
  '''
  print "___________________________________________"
  print "              sorted(zip(*a)) =", sorted(zip(*a))
  print "sorted(zip(*a), key = getKey) =", sorted(zip(*a), key = getKey)


  print "___________________________________________"
  '''
  lst = [4,2,4,2,5,8,5,2,54,7,9]
  print "lst =", lst
  print "sorted(lst) =", sorted(lst)
  print "tuple(lst) =", tuple(lst)
  '''
