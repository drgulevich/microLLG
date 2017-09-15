from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

interpolation='nearest'
#interpolation='gaussian'
#interpolation=None

### ---  Display --- 
rc('font', family='serif')
rc('font', size='14')

#cmap=plt.get_cmap('bwr')
#cmap=plt.get_cmap('seismic_r')
cmap=plt.get_cmap('coolwarm_r')

### --- Parameters ---

fig, ax = subplots(figsize=(6,6))
#ax.set_title('(a)')

frame=-1
mx=magdata[frame,1:Nx+1,1:Ny+1,0]
my=magdata[frame,1:Nx+1,1:Ny+1,1]
mz=magdata[frame,1:Nx+1,1:Ny+1,2]

im=ax.imshow(mz.T,interpolation=interpolation, cmap = cmap, origin='lower',vmin=-1,vmax=1,zorder=1)
width=0.0025
scale=1.0

Q = ax.quiver(mx.T,my.T,pivot='mid',zorder=2,width=width, scale=scale, scale_units='x', headwidth=6, headlength=8)

#props = dict(boxstyle='round',facecolor='white', alpha=0.7)
#ax.text(0.06, 0.94, '(d)', fontsize=15, verticalalignment='top', bbox=props,transform=ax.transAxes)
#ax.text(0.1, 0.9, '(a)', fontsize=20, verticalalignment='top', transform=ax.transAxes, color='w')

fig.colorbar(im, label=r'$m_z$',orientation='vertical')

#ax.set_xticks([])
#ax.set_yticks([])
show()
