#!/usr/bin/env python
###
### Animation example: http://matplotlib.org/examples/animation/dynamic_image.html
###
from pylab import *
import matplotlib.animation as animation

interpolation='nearest'
#interpolation='gaussian'
#interpolation=None

#cmap=plt.get_cmap('bwr')
#cmap=plt.get_cmap('seismic_r')
cmap=plt.get_cmap('coolwarm_r')

### --- Parameters ---

#fig, ax = subplots(figsize=(6,6))
fig, ax = subplots(figsize=(12,12))
#subplots_adjust(left=0, right=1, bottom=0, top=1)

mx=magdata[0,1:Nx+1,1:Ny+1,0]
my=magdata[0,1:Nx+1,1:Ny+1,1]
mz=magdata[0,1:Nx+1,1:Ny+1,2]

im=ax.imshow(mz.T,interpolation=interpolation, cmap = cmap, origin='lower',vmin=-1,vmax=1,zorder=1)
#im=ax.imshow(mz.T,interpolation=interpolation, cmap = cmap, origin='lower',extent=[1,Nx,1,Ny],vmin=-1,vmax=1,zorder=1)
#im=imshow(magdata[0,:,:,2].T,interpolation=interpolation, cmap = cmap, origin='lower',vmin=-1,vmax=1)

#width=0.0016
#scale=1

#width=0.0012
#scale=0.8

width=0.0015
scale=1.1

#X, Y = meshgrid(np.arange(1,Nx+1),np.arange(1,Ny+1))
#Q = ax.quiver(X, Y, mx.T,my.T,pivot='mid',zorder=2,width=width, scale=scale, scale_units='x')

Q = ax.quiver(mx.T,my.T,pivot='mid',zorder=2,width=width, scale=scale, scale_units='x')

#Q = ax.quiver(X, Y, U, V, pivot='mid', color='r', units='inches')

mt = text(.5, .5, 't=%.2f' % 0., fontsize=15)
#mt = text(1.5, 1.5, 't=%.2f' % 0., fontsize=15)

#time_text = text(.5, .5, '', fontsize=15)

def init():
    return updatefig(0)
	
def updatefig(frame):
    data=magdata[frame,1:Nx+1,1:Ny+1,2].T
    im.set_array(data)
    Q.set_UVC(magdata[frame,1:Nx+1,1:Ny+1,0].T, magdata[frame,1:Nx+1,1:Ny+1,1].T)
    mt.set_text('t=%.2f' % (frame*countout*dt))
    return im,Q,mt,

def animate_as_gif(frame):
    return updatefig(frame)

#export = True
export = False

if(export==True):
    anim = animation.FuncAnimation(fig, animate_as_gif, np.arange(0, Nframes), init_func=init, interval=100, blit=True, repeat=False)
    anim.save('animation.gif', writer='imagemagick')
else:
    #anim = animation.FuncAnimation(fig, updatefig, np.arange(1, Nframes), init_func=init, interval=500, blit=True, repeat=False)
    anim = animation.FuncAnimation(fig, updatefig, np.arange(0, Nframes), init_func=init, interval=100, blit=True, repeat=False)   
#    anim = animation.FuncAnimation(fig, updatefig, np.arange(0, Nframes,100), init_func=init, interval=100, blit=True, repeat=False)   

fig.tight_layout()
gca().set_aspect('equal', adjustable='box')
#axis('off')

show()
