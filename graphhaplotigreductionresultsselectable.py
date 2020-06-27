#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""
import pandas as pd
from bokeh.io import output_notebook, show
from bokeh.models import ColumnDataSource, HoverTool, WheelZoomTool, Select
from bokeh.plotting import figure, output_file, show, curdoc
from bokeh.palettes import Spectral5
from bokeh.layouts import row, column
#this is for jupyter
#output_notebook()

SIZES = list(range(6, 22, 3))
COLORS = Spectral5
N_SIZES=len(SIZES)
N_COLORS=len(COLORS)

mydf = pd.read_csv('haplotigreduction.csv')
mydf.columns = ['asmsize', 'n50size', 'cbusco', 'sbusco', 'dbusco', 'fbusco', 'mbusco', 'nbuscos', 'mincontig', 'ovlrange', 'pct']
columns = sorted(mydf.columns)
output_file('haplotigreductiongraph.html')

def create_figure():
    xs = mydf[x.value].values
    ys = mydf[y.value].values
    x_title = x.value.title()
    y_title = y.value.title()

    kw = dict()
    kw['x_range'] = sorted(set(map(str, xs)))
    kw['y_range'] = sorted(set(map(str, ys)))
    kw['title'] = '%s vs %s' % (x_title, y_title)
    TOOLTIPS = [('Index','$index'),('(sbusco,dbusco)','($x{0.0000},$y{0.0000})'),('asmSize','@asmsize'),('n50','@n50size'),('min contig','@mincontig'),('Overlap Range','@ovlrange'),('% Identity','@pct{0.0000000}'),('mbusco','@mbusco{0.0000}')]
    #p = figure(title='Single vs Duplicate BUSCOS', x_axis_label='Single BUSCO"s', y_axis_label='Duplicate BUSCOS"s', tools="pan,wheel_zoom,box_zoom,reset,hover", tooltips=TOOLTIPS)
    p = figure(plot_height=600, plot_width=800, tools="pan,wheel_zoom,box_zoom,reset,hover", tooltips=TOOLTIPS, **kw)
    p.xaxis.axis_label = x_title
    p.yaxis.axis_label = y_title
    p.xaxis.major_label_orientation = pd.np.pi / 4
    sz = 9
    if size.value != 'None':
        if len(set(mydf[size.value])) > N_SIZES:
            groups = pd.qcut(mydf[size.value].values, N_SIZES, duplicates='drop')
        else:
            groups = pd.Categorical(mydf[size.value])
        sz = [SIZES[xx] for xx in groups.codes]

    c = '#31AADE'
    if color.value != 'None':
        if len(set(mydf[color.value])) > N_SIZES:
            groups = pd.qcut(mydf[color.value].values, N_COLORS, duplicates='drop')
        else:
            groups = pd.Categorical(mydf[color.value])
        c = [COLOR[xx] for xx in groups.codes]

    #p.circle(x=xs, y=ys, color=c, size=sz, line_color='white', alpha=0.6, hover_color='white', hover_alpha=0.5)
    p.scatter(x=xs, y=ys, color=c, size=sz, line_color='white', alpha=0.6, hover_color='white', hover_alpha=0.5)
    return p

def update(attr, old, new):
    layout.children[1] = create_figure()

x = Select(title='X-Axis', value='sbusco', options=columns)
x.on_change('value', update)
y = Select(title='Y-Axis', value='dbusco', options=columns)
y.on_change('value', update)

size = Select(title='Size', value='None', options=['None'] + continuous)
size.on_change('value', update)

color = Select(title='Color', value='None', options=['None'] + continuous)
color.on_change('value', update)

controls = column([x, y, color, size], width=200)
layout = row(controls, create_figure())

curdoc().add_root(layout)
curdoc().title = "Crossfilter"
