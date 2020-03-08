import os 
os.chdir('/home/fangyuan/projects/hic/codes')

for f in os.listdir('../data/cool_balanced'): 
    os.system(f'cooler balance {f}') 
