#import sys
import subprocess

print('##############################################################################')
print(' Couette ')
print('##############################################################################')
output = subprocess.check_output(['python','couette/compare_couette_solutions.py'])
print(output)



print('##############################################################################')
print(' Poiseuille')
print('##############################################################################')
output = subprocess.check_output(['python','poiseuille/compare_poiseuille_solutions.py'])
print(output)

