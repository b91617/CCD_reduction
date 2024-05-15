#%%
import os
from pathlib import Path

DIR = '/Users/b91617/Library/CloudStorage/OneDrive-cycu.org.tw/LOT/RR Lyra/LOT20220723/calibrated_data/reduced_lights'
files = os.listdir(DIR)

p = Path(DIR)

# %%
for file in files:
    print(file, end='\n')

#%%
for file in files:
    with open("list.txt", "w") as f:
        print(file, file=f)
