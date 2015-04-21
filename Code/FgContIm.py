from astropy import log
log.setLevel('ERROR')
import os, glob
import aplpy       # http://aplpy.readthedocs.org/en/stable/howto_subplot.html
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

path = '../FinalData/'
Plotpath = '../Figure'
label = dict.fromkeys(['VLA', 'cont', 'lin'])
for k in label.iterkeys():
    files = glob.glob(path+'*'+k+'*.fits')
    label[k] = files

# set this ratio for 3 panels (~ square each)
fig = plt.figure(figsize=(9, 3))


# Change contours of yellow VLA 9GHz, take out CARMA zoomed out panel, zoom out VLA panel to inlcude noise.

def get_RA_dec(f):
    """
    Get RA, dec from fits

    Input:
    f: str
        filename with .fits extension
    """
    import atpy
    t = atpy.Table(f)
    return t.ra, t.dec


def sigma_contour_array(sigma):
    """
    return list of sigma multiples, starting 4, 8, 16, 32..
     -- > change to higher sigma, and multiply by ~ 2 for VLA 9 GHz
    """
    arr = range(11)
    arr.remove(0)
    arr.remove(1)
    arr = [2**i * sigma for i in arr]
    arr.append(list(np.array(arr) * -1)[0])
    print arr
    return arr

def sigma_contour_CARMA(sigma):
    """
    return list of sigma multiples,
    """
    arr = range(-3, 11)
    arr.remove(0)
    arr.remove(1)
    arr.remove(-1)
    arr = [i * np.sqrt(2) * sigma for i in arr]
    return arr


def standard_plot_setup(sp, size):
    sp.set_frame_color('black')
    sp.frame.set_linewidth(1)
    sp.set_system_latex(True)
    sp.recenter(x=ra_center, y=dec_center, radius=size)   # width = blah, height = blah     # Scaling and panning

    sp.tick_labels.set_font(size='medium', weight='bold')      # size='10'
    sp.tick_labels.set_style('colons')
    sp.tick_labels.set_xformat('hh:mm:ss')
    sp.tick_labels.set_yformat('dd:mm:ss')
    sp.ticks.set_color('black')
    sp.ticks.set_linewidth(2)
    sp.ticks.set_length(10)
    sp.ticks.set_minor_frequency(4)
    # sp.ticks.set_xspacing(45*15/3600.)        # deg

    sp.axis_labels.set_font(size='medium', weight='bold')      # (size='12')
    sp.axis_labels.set_xtext('Right Ascension (J2000)')
    sp.axis_labels.set_ytext('Declination (J2000)')
    sp.axis_labels.set_xpad(5)
    sp.axis_labels.set_ypad(-15)
    sp.axis_labels.set_xposition('bottom')


def inset_frmt(sp):
    """
    Hide the ticks, tick labels and axis labels (e.g. RA & Dec.)
    """
    sp.ticks.hide()
    sp.tick_labels.hide()
    sp.axis_labels.hide()


def setup_beam(sp, loc='bottom left', c='black', hatch=None, a=0.8, lw=3):
    """
    Auto fetch beam info from header
    Inputs:
    --------
    loc: str
    hatch: str
        one of /,|,-,+,x,o,O,.,*
    """
    sp.show_beam(corner=loc, color=c, fill=True, edgecolor='grey', facecolor='black', linestyle='solid', linewidth=lw, frame=True, alpha=a)
    if hatch != None:
        sp.beam.set_hatch(hatch)


# def setup_beam2(sp, loc='top right', major, minor, angle, c='black', hatch=None, a=0.8, lw=3):
#     """
#     Manual input beam inf

#     major: float
#         deg
#     minor: float
#         deg
#     PA: float
#         deg
#     """

#     sp.add_beam(0,0,0)
#     sp.beam.set_frame(True)
#     sp.beam.show(major, major, angle)
#     sp.beam.set(corner=loc, color=c, fill=True, edgecolor='grey', facecolor='black', linestyle='solid', linewidth=lw, frame=True, alpha=a)
#     sp.beam.set_major(major)  # degrees
#     sp.beam.set_minor(major)  # degrees
#     sp.beam.set_angle(angle)  # degrees
#     sp.beam.set_corner(loc)
#     sp.beam.set_alpha(a)
#     sp.beam.set_color('white')
#     sp.beam.set_edgecolor('white')
#     sp.beam.set_facecolor('gray')
#     sp.beam.set_linestyle('solid')
#     sp.beam.set_linewidth(1)


def markers_cross(sp, ra, dec, layer=None, ec='yellow', fc='none', mk='+', s=500, a=1.0, lw=2):
    """
    Inputs:
    -------
    ra: float
    dec: float
    layer: str
    """
    sp.show_markers(ra, dec, layer=layer, edgecolor=ec, facecolor=fc, marker=mk, s=s, alpha=a, linewidth=lw)


def put_label(sp, x, y, text, layer, c='yellow', s='x-large', w='bold'):
    """
    Inputs:
    -------
    x: float
    y: float
    text: str
    layer: str
    """
    sp.add_label(x, y, text, relative=True, color=c, size=s, layer=layer, weight=w)


def frmt_colorbar(sp):
    """
    format colorbar
    """
    sp.colorbar.set_axis_label_font(size=15)
    sp.colorbar.set_width(0.1)
    sp.colorbar.set_axis_label_pad(15)
    sp.colorbar.set_axis_label_text('$Flux Density [mJy]$')


def setup_scalebar(sp, lg, label, unit=None, loc='bottom right', a=0.9, c='white', lw=3):
    """
    unit: str
        unit for label
    lg: float
        length [degrees]
    color: str
        white, lime, Red, DogderBlue, yellow, black, etc
    label: str
    """
    sp.show_scalebar(lg, color=c, corner=loc, alpha=a, linestyle='solid', linewidth=lw)
    sp.scalebar.set_font(size='large', weight='bold')
    # if unit == None:
    #     unit = r'$kpc / {\prime\prime}$'
    sp.scalebar.set_label(label)


########################################
# user define area
########################################
sigma_9ghz = 6.41094e-05
sigma_cont = 0.829972e-3       # 0.38466e-3    # 0.829972e-3
sigma_line = 1.57804

vla_min = -0.000270099
vla_max = 0.00960869
max_cont = 5.347e-3
min_cont = -1.94e-3
max_line = 9.609e-3
min_line = -0.2700985e-3

VLA_ra = 1.448437499980E+02
VLA_dec = 8.325666666846E+01
sizep = 0.0025
sizePcont_noise = 0.006
sizePLine = 0.0057
ra_center = 144.84933
dec_center = 83.257175

# H_0 = 70, Omega_m = 0.3, flat Universe
z_radio = 0.685
scale_radio = 7.083     # kpc/"
z_SMG = 2.221
scale_SMG = 8.252       # kpc/"

# # positions for crosses, centering on radio core
ra_cross, dec_cross = ra_center, dec_center
row_a = 0.05
width = 0.32
x0 = 0.05
dy = 1.0


########################################
# intialize base figure
########################################
fvla = aplpy.FITSFigure(label['VLA'][0], \
        figure=fig, subplot=[x0,row_a,width,dy])
# fvla.show_grayscale(stretch='log', vmin=vla_min, vmax=vla_max, vmid=-0.001)

fcont = aplpy.FITSFigure(label['cont'][0], \
        figure=fig, subplot=[x0+width,row_a,width,dy])
fcont.show_colorscale(cmap=mpl.cm.jet, stretch='log', vmin=min_cont, vmax=max_cont, vmid=-0.5)

flin = aplpy.FITSFigure(label['lin'][0], \
        figure=fig, subplot=[x0+width*2,row_a,width,dy])
flin.show_colorscale(cmap=mpl.cm.jet)#, vmin=min_line, vmax=max_line, vmid=min_line-0.01, stretch='log')

########################################
# Contours
########################################
fvla.show_contour(label['VLA'][0], colors="lime", levels=sigma_contour_array(sigma_9ghz), linewidths=2, layer='fg')
fvla.show_contour(label['cont'][0], colors='red', levels=sigma_contour_CARMA(sigma_cont), linewidths=2, layer='fg_cont')
fcont.show_contour(label['cont'][0], colors='white', levels=sigma_contour_CARMA(sigma_cont), linewidths=1.5, layer='fg')
flin.show_contour(label['lin'][0], colors="white", levels=sigma_contour_CARMA(sigma_line), linewidths=2, layer='mol')


########################################
# beam
########################################
setup_beam(fvla)
bmaj = fcont._header['BMAJ']
bmin = fcont._header['BMIN']
bpa = fcont._header['BPA']
#setup_beam(fvla, loc='top right', major=bmaj, minor=bmin, angle=bpa)

setup_beam(fcont)
setup_beam(flin)


########################################
# scale bar
########################################
lg = 1./3600

# --> change to ~ 20kpc
setup_scalebar(fvla, lg, str(scale_radio))
setup_scalebar(fcont, lg, str(scale_radio))
setup_scalebar(flin, lg, str(scale_SMG))


########################################
# axes
########################################
standard_plot_setup(fvla, sizep)
standard_plot_setup(fcont, sizePcont_noise)
standard_plot_setup(flin, sizePLine)
# fcont.tick_labels.hide()
fcont.axis_labels.hide()
# flin.tick_labels.hide()
flin.axis_labels.hide()


########################################
# markers
########################################
markers_cross(fvla, ra_cross, dec_cross, layer='marker_set_1')
markers_cross(fcont, ra_cross, dec_cross, layer='marker_set_1')
markers_cross(flin, ra_cross, dec_cross, layer='marker_set_1')


########################################
# Labels
########################################
put_label(fvla, 0.31, 0.95, 'VLA 9GHz', 'titleBand')
put_label(fvla, 0.2625, 0.9, '3C220.3', 'titleObj')
put_label(fcont, 0.31, 0.95, 'CARMA 104GHz', 'titleBand')
put_label(fcont, 0.2625, 0.9, '3C220.3', 'titleObj')
put_label(flin, 0.31, 0.95, 'CARMA CO(3-2)', 'titleBand')
put_label(flin, 0.2625, 0.9, 'SMM J0939+8315', 'titleObj')

labsize = 'xx-large'
labc = 'limegreen'
put_label(fvla, 0.80, 0.85, '(a)', 'ref', c=labc, s=labsize)
put_label(fcont, 0.80, 0.85, '(b)', 'ref', c=labc, s=labsize)
put_label(flin, 0.80, 0.85, '(c)', 'ref', c=labc, s=labsize)

########################################
# Colorbar
########################################
axisflin = fig.add_axes([0.92,0.19,0.02,0.68])
normflin = mpl.colors.Normalize(vmin=min_line, vmax=max_line)
cbflin = mpl.colorbar.ColorbarBase(axisflin, cmap=mpl.cm.jet, norm=normflin, orientation='vertical')
cbflin.set_label('mJy')

fig.canvas.draw()
plt.show()


# if save==True:
#     os.system('rm -rf ImagePanels.png ImagePanels.eps')
#     fig.savefig("ImagePanels.eps", dpi=600)
#     fig.savefig("ImagePanels.png", dpi=300)

import sys; sys.exit()