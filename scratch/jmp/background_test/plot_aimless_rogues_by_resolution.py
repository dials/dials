

if __name__ == '__main__':

  import sys
  filename = sys.argv[1]

  with open(filename) as infile:

    H = []
    K = []
    L = []
    I = []
    S = []
    for line in infile.readlines():
      tokens = line.split()
      if len(tokens) in [19, 20, 21]:
        try:
          H.append(int(tokens[0]))
          K.append(int(tokens[1]))
          L.append(int(tokens[2]))
          I.append(int(tokens[7]))
          S.append(int(tokens[8]))
        except Exception, e:
          pass
    assert len(H) == len(K)
    assert len(H) == len(L)
    assert len(H) == len(I)
    assert len(H) == len(S)

    with open("rogues.txt", "w") as outfile:
      print "Writing to file: rogues.txt"
      print "Columns: H K L I SigI"
      for h, k, l, i, s in zip(H, K, L, I, S):
        outfile.write("%d %d %d %f %f\n" % (h,k,l,i,s))

