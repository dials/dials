def run():
  from cctbx import sgtbx

  awesomeness_index = []
  for i in range(230):
    sg = sgtbx.space_group_info(number=i+1).group()
    awesomeness_index.append(len(sg)/(i+1))

  from dials.array_family import flex
  perm = flex.sort_permutation(flex.double(awesomeness_index), reverse=True)
  for rank, i in enumerate(perm):
    sgi = sgtbx.space_group_info(number=i+1)
    print "%3i %.2f" %(rank+1, awesomeness_index[i]), sgi

  from matplotlib import pyplot
  pyplot.scatter(range(1, 231), awesomeness_index)
  pyplot.xlabel('Space group number')
  pyplot.ylabel('Space group awesomeness index')
  pyplot.show()


if __name__ == '__main__':
  run()
