#%%
import astropy.io.fits as ap
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

fits_image_filename = ''
hdul = ap.open(fits_image_filename)
data = hdul[0].data

plt.figure()
plt.imshow(data, cmap='gray')
plt.colorbar()
# %%
