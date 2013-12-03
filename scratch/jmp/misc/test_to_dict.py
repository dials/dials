



from dxtbx.model import Detector
d = Detector()
p1 = d.add_panel()
p2 = d.add_panel()
p3 = d.add_panel()
p4 = d.add_panel()


root = d.hierarchy()
g = root.add_group()
g.add_panel(d[0])
g.add_panel(d[1])
root.add_panel(d[2])
root.add_panel(d[3])
print d.to_dict()
