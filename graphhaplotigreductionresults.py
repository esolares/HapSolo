#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""
import pandas as pd
from bokeh.io import output_notebook, show
from bokeh.models import ColumnDataSource, HoverTool, WheelZoomTool
from bokeh.plotting import figure, output_file, show
output_notebook()

mydf = pd.read_csv('haplotigreduction.csv')
mydf.columns = ['asmsize', 'n50size', 'cbusco', 'sbusco', 'dbusco', 'fbusco', 'mbusco', 'nbuscos', 'mincontig', 'ovlrange', 'pct']
output_file('haplotigreductiongraph.html')
TOOLTIPS = [('Index','$index'),('(sbusco,dbusco)','($x{0.0000},$y{0.0000})'),('asmSize','@asmsize'),('n50','@n50size'),('min contig','@mincontig'),('Overlap Range','@ovlrange'),('% Identity','@pct{0.0000000}'),('mbusco','@mbusco{0.0000}')]
p = figure(title='Single vs Duplicate BUSCOS', x_axis_label='Single BUSCO"s', y_axis_label='Duplicate BUSCOS"s', tools="pan,wheel_zoom,box_zoom,reset,hover", tooltips=TOOLTIPS)
p.scatter('sbusco', 'dbusco', source = mydf)
top5asms = mydf.sort_values(['mbusco','dbusco','sbusco','asmsize'],ascending=[True,True,False,False]).head(5)
print(top5asms)
show(p)
