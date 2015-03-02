import numpy as np

blah = [1,2,3]
bluh = [11,12,13]
dico = {"A": blah, "B":bluh}
print([x[0] for x in dico.values()])
