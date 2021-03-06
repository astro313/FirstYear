from astropy import log
log.setLevel('ERROR')
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
from FgContIm import *


font = {'family': 'Arial',
        'weight': 'normal'}
mpl.rc('font', **font)

path = '../FinalData/'
Plotpath = '../Figure/'
label = dict.fromkeys(['VLA', 'cont', 'lin', 'SMA'])
for k in label.iterkeys():
    files = glob.glob(path+'*'+k+'*.fits')
    label[k] = files

print label

fig = plt.figure(3, figsize=(12, 5))
fig_line = plt.figure(2, figsize=(12, 5))


########################################
# user define area
########################################
sigma_9ghz          = 6.41094e-05
sigma_cont          = 5e-4
sigma_line          = 1.03 # 1.365  # chan 65,85     #     1.0  # chan72,90
sigma_SMA           = 0.000836951

vla_min             = -0.000270099
vla_max             = 0.00960869
max_cont            = 4.796309E-03
min_cont            = -0.02 # -2.038660E-03
max_line            = 8.93562
min_line            = -8.836704
SMA_min             = -0.00210908
SMA_max             = 0.00634805

VLA_ra              = 1.448437499980E+02
VLA_dec             = 8.325666666846E+01
sizePcont_noise     = 0.020
sizePLine           = 0.02 #0.012
sizevla             = 0.0029
ra_center           = 144.84933
dec_center          = 83.257175

# WMAP9, flat Universe
z_radio             = 0.685
scale_radio         = 7.169     # kpc/"
z_SMG               = 2.221
scale_SMG           = 8.406       # kpc/"

# # positions for crosses, centering on radio core
ra_cross, dec_cross = ra_center, dec_center
row_a               = 0.10
width               = 0.35
x_gap               = 0.05
x0                  = 0.10
dy                  = 0.90


# positions for crosses, centering on SMA peaks
RA1  = 144.845625      # 09:39:22.950
Dec1 = 83.2575472222222  # +83:15.27.17


RA2  = 144.850916666  # 09:39:24.220
Dec2 = 83.2575047   # +83:15:27.017


RA3  = 144.853454166      # 09:39:24.829
Dec3 = 83.257044         # +83:15.25.36



########################################
# intialize base figure
########################################

fcont = aplpy.FITSFigure(label['cont'][0], \
        figure=fig, subplot=[x0,row_a,width,dy])
fcont.show_colorscale(cmap=mpl.cm.jet, stretch='log', vmin=min_cont, vmax=max_cont, vmid=-10000000)   # vmid= min_cont-0.0001 ~ orange; min_cont-0.0000001 ~ red

fvla = aplpy.FITSFigure(label['VLA'][0], \
        figure=fig, subplot=[x0+width+2*x_gap, row_a, width, dy])
fvla.show_grayscale(stretch='log', vmin=vla_min, vmax=vla_max, vmid=-0.001)


flin = aplpy.FITSFigure(label['lin'][0], \
        figure=fig_line, subplot=[x0,row_a,width,dy])
flin.show_colorscale(cmap=mpl.cm.jet, vmin=min_line+1, vmax=max_line-1, vmid=-115, stretch='log')

fSMA = aplpy.FITSFigure(label['SMA'][0], figure=fig_line, subplot=[x0+width+2*x_gap, row_a, width, dy])
fSMA.show_grayscale(stretch='log', vmin=SMA_min, vmax=SMA_max, vmid=-0.1)

########################################
# Contours
########################################
fcont.show_contour(label['cont'][0], colors='white', alpha=1, levels=sigma_contour_cont(sigma_cont), linewidths=1.8)#, layer='fg')
fvla.show_contour(label['VLA'][0], colors="lime", levels=sigma_contour_array(sigma_9ghz), linewidths=2)#, layer='fg')
fvla.show_contour(label['cont'][0], colors='red', levels=sigma_contour_cont(sigma_cont), linewidths=2, layer='fg_cont')
flin.show_contour(label['lin'][0], colors="white", alpha=1, levels=sigma_contour_CARMA(sigma_line), linewidths=1.8, layer='mol')
fSMA.show_contour(label['SMA'][0], colors="lime", levels=sigma_contour_CARMA(sigma_SMA), linewidths=1.8, layer='bf_cont')
fSMA.show_contour(label['lin'][0], colors="red", levels=sigma_contour_CARMA(sigma_line), linewidths=2, layer='mol')

########################################
# beam
########################################
setup_beam(fvla)
bmaj_cont = fcont._header['BMAJ']
bmin_cont = fcont._header['BMIN']
bpa_cont = fcont._header['BPA']
# setup_beam2(fvla, bmaj_cont, bmin_cont, bpa_cont, 1)
setup_beam(fcont)
setup_beam(flin)
setup_beam(fSMA)
bmaj_line = flin._header['BMAJ']
bmin_line = flin._header['BMIN']
bpa_line = flin._header['BPA']
# setup_beam2(fSMA, bmaj_line, bmin_line, bpa_line, 1)


########################################
# scale bar
########################################
lg_1arcsec = 1./3600
lg_20kpc_fg = lg_1arcsec * 20./scale_radio
lg_20kpc_bg = lg_1arcsec * 20./scale_SMG
setup_scalebar(fvla, lg_20kpc_fg, str('20kpc'))
setup_scalebar(fcont, lg_20kpc_fg, str('20kpc'), c='black')
setup_scalebar(flin, lg_20kpc_bg, str('20kpc'), c='black')
setup_scalebar(fSMA, lg_20kpc_bg, str('20kpc'))

########################################
# axes
########################################
standard_plot_setup(fcont, ra_center, dec_center, sizePcont_noise) #, tc='white')
standard_plot_setup(fvla, ra_center, dec_center, sizevla, tc='white')
standard_plot_setup(flin, ra_center, dec_center, sizePLine)
standard_plot_setup(fSMA, ra_center, dec_center, sizevla, tc='white')
# fcont.tick_labels.hide()
fvla.axis_labels.hide()
# flin.tick_labels.hide()
fSMA.axis_labels.hide()


########################################
# markers
########################################
markers_cross(fvla, ra_cross, dec_cross)# , layer='marker_set_1')
markers_cross(fcont, ra_cross, dec_cross, ec='black') #), layer='marker_set_1')
markers_cross(flin, ra_cross, dec_cross, ec='black')#, layer='marker_set_1', s=500)
markers_cross(fSMA, ra_cross, dec_cross)#, layer='marker_set_1')

########################################
# markers to satisfy referee
########################################
markers_cross(fvla, RA1, Dec1, ec='cyan')
markers_cross(fvla, RA2, Dec2, ec='cyan')
markers_cross(fvla, RA3, Dec3, ec='cyan')


########################################
# Labels
########################################
put_label(fvla, 0.40, 0.93, 'VLA 9GHz, CARMA 104GHz', 'titleBand')
put_label(fvla, 0.17, 0.87, '3C220.3', 'titleObj')
put_label(fcont, 0.31, 0.90, 'CARMA 104GHz', 'titleBand', c='black')#, w=
#'heavy')
put_label(fcont, 0.2625, 0.85, '3C220.3', 'titleObj', c='black')#, w='extra bold')
put_label(flin, 0.31, 0.93, 'CARMA CO(3-2)', 'titleBand', c='black')
put_label(flin, 0.31, 0.87, 'SMM J0939+8315', 'titleObj', c='black')
put_label(fSMA, 0.40, 0.935, 'SMA 1 mm, CARMA CO(3-2)', 'titleBand')
put_label(fSMA, 0.2625, 0.85, 'SMM J0939+8315', 'titleObj')

labsize = 'xx-large'
labc = 'white'

# put_label(fvla, 0.80, 0.925, '(a)', 'ref', c=labc, s=labsize)
# put_label(fcont, 0.80, 0.925, '(b)', 'ref', c=labc, s=labsize)
# put_label(flin, 0.80, 0.925, '(c)', 'ref', c=labc, s=labsize)
# put_label(fSMA, 0.80, 0.925, '(d)', 'ref', c=labc, s=labsize)
########################################
# Colorbar
########################################
# axisflin = fig.add_axes([0.92,0.19,0.02,0.68])
# normflin = mpl.colors.Normalize(vmin=min_line, vmax=max_line)
# cbflin = mpl.colorbar.ColorbarBase(axisflin, cmap=mpl.cm.jet, norm=normflin, orientation='vertical')
# cbflin.set_label('mJy')
# fig.canvas.draw()
# fig_line.canvas.draw()


if __name__ == '__main__':
    """
    run script.py True
    sys.argv[1] determines whether to save all the plots or not
    """
    import sys, os
    if len(sys.argv) < 2:
        errmsg = "Invalid number of arguments: {0:d}\n  run script.py Save_True"
        raise IndexError(errmsg.format(len(sys.argv)))
    saveFig = True if sys.argv[1].lower() == 'true' else False
    if saveFig == True:
        outfile = 'ContPanel'
        outfile_line = 'LinePanel'
        os.system('rm -rf ' + outfile + '.png' + ' ' + outfile + '.eps' + ' ' +outfile + '.pdf')
        os.system('rm -rf ' + outfile_line + '.png' + ' ' + outfile_line + '.eps' + ' ' + outfile_line + '.pdf')
        fig_line.savefig('../Figure/'+outfile_line+'.pdf')
        fig_line.savefig('../Figure/'+outfile_line+'.eps')
        fig_line.savefig('../Figure/'+outfile_line+'.png', dpi=300)
        fig.savefig('../Figure/'+outfile+'.pdf')
        fig.savefig('../Figure/'+outfile+'.eps')
        fig.savefig('../Figure/'+outfile+'.png', dpi=300)
    else:
        plt.show()

