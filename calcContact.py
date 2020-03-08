import os 
meta = 't47d'
os.chdir(f'/home/fangyuan/projects/hic/data/cool_unbalanced/{meta}')

# os.system('mkdir /home/fangyuan/projects/hic/data/contact/{meta}')
for f in os.listdir():
    if f.endswith('.cool') and 'mapq' in f:
        sp = f.split('.')[0]
        os.system(f'cooler dump -t pixels --header --join {f} > /home/fangyuan/projects/hic/data/contact/{meta}/{sp}.txt')