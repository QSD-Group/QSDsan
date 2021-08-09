#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as LSC
from colorpalette import Color, Palette

__all__ = ('hex2rgb', 'rgb2hex', 'test_colormaps', 'palettes', 'colormaps')


# %%

# =============================================================================
# Util functions
# =============================================================================

def hex2rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))


def rgb2hex(r, g, b):
    return f'#{r:02x}{g:02x}{b:02x}'


def test_colormaps(colormaps=[]):
    '''
    Generate a gradient plot to test the colormap from [1]_.

    References
    ----------
    .. [1] `Creating Colormaps in Matplotlib <https://matplotlib.org/stable/tutorials/colors/colormap-manipulation.html>`_
    '''
    np.random.seed(19680801)
    data = np.random.randn(30, 30)
    n = len(colormaps)
    fig, axs = plt.subplots(1, n, figsize=(n * 2 + 2, 3),
                            constrained_layout=True, squeeze=False)
    for [ax, cmap] in zip(axs.flat, colormaps):
        psm = ax.pcolormesh(data, cmap=cmap, rasterized=True, vmin=-4, vmax=4)
        fig.colorbar(psm, ax=ax)
    plt.show()


# %%

# =============================================================================
# Group colors, note that bg colors are to be used as the background
# (lighter than fg colors)
# =============================================================================

palettes = {} # color palettes
RGBn = {} # normalized (i.e., 0-1) RGBn values
colormaps = {} # colormaps for gradient coloring

# Guest group colors
Guest = Palette()
Guest.purple = Color('purple', (162, 128, 185))
Guest.blue = Color('blue', (96, 193, 207))
Guest.green = Color('green', (121, 191, 130))
Guest.yellow = Color('yellow', (243, 195, 84))
Guest.orange = Color('orange', (249, 143, 96))
Guest.red = Color('red', (237, 88, 111))
Guest.gray = Color('gray', (144, 145, 142))

RGBn['Guest'] = [
    Guest.red.RGBn,
    Guest.purple.RGBn,
    Guest.blue.RGBn,
    Guest.green.RGBn,
    Guest.yellow.RGBn
    ]

palettes['Guest'] = Guest
colormaps['Guest'] = LSC.from_list('Guest', RGBn['Guest'])
# # If want to anchor color at a specific point
# colormaps['Guest'] = LSC.from_list('Guest', list(zip([0, 0.2, 0.4, 0.6, 1], RGBn['Guest'])))


# CABBI colors
CABBI = Palette()
CABBI.feedstock = Color('feedstock', (45, 130, 63))
CABBI.conversion = Color('conversion', (243, 195, 84))
CABBI.sustainability = Color('sustainability', (21, 145, 118))

CABBI.green_light = Color('green_light', fg=(142, 173, 62))
CABBI.green_dark = Color('green_dark', fg=(0, 127, 61), bg=(59, 164, 89))

CABBI.yellow = Color('yellow', fg=(252, 184, 19), bg=(255, 221, 80))

CABBI.blue_light = Color('blue_light', fg=(0, 169, 150), bg=(178, 224, 229))
CABBI.blue_dark = Color('blue_dark', fg=(26, 132, 118), bg=(130, 207, 208))

CABBI.brown_light = Color('brown_light', fg=(152, 135, 110), bg=(225, 222, 218))
CABBI.brown_dark = Color('brown_dark', fg=(64, 58, 72), bg=(152, 135, 110))

RGBn['CABBI'] = [
    CABBI.brown_dark.RGBn,
    CABBI.green_dark.RGBn,
    CABBI.green_light.RGBn,
    CABBI.yellow.RGBn_bg,
    CABBI.conversion.RGBn,
    CABBI.yellow.RGBn
    ]

palettes['CABBI'] = CABBI
colormaps['CABBI'] = LSC.from_list('CABBI', RGBn['CABBI'])


# test_colormaps([i for i in colormaps.values()])