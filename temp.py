from main import NumericalProcedures
import matplotlib.pyplot as plt

procedures = NumericalProcedures(-5, 0, 1000, 2)
'''
err = procedures.error()
plt.plot(procedures.x, err)
plt.show()
'''
print(procedures)
