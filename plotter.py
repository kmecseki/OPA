import matplotlib.pyplot as plt
import pandas as pd


df = pd.read_csv('output/Spec_sig.dat', sep='\t', header=None, names=['wl', 'int'])

# Plot the data
plt.plot(df['wl'], df['int'])
plt.xlabel('Wavelength')
plt.ylabel('Intensity')
plt.title('Spectrum')
plt.show()

