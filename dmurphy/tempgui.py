"""
Dr Kenneth Duncan - April 2018
"""
import numpy as np
import ipywidgets as widgets
from IPython.display import display
from collections import OrderedDict
from six.moves import urllib
from math import cos,radians

obstemp = "http://casu.ast.cam.ac.uk/~dmurphy/opr3/swg/resources/obstemp.dat"
progtemp = "http://casu.ast.cam.ac.uk/~dmurphy/opr3/swg/resources/progtemp.dat"

style = {'description_width': 'initial'}

Scodes = []
for letter in range(97, 121):
    Scodes.append(chr(letter).upper())

try:
    obstemp_data = urllib.request.urlopen(obstemp).read().decode('utf-8').split('\n')
except:
    raise SystemExit('Could not get data from URL %s'%(obstemp))
    a = 1
    #raise a warning here

try:
    progtemp_data = urllib.request.urlopen(progtemp).read().decode('utf-8').split('\n')
except:
    raise SystemExit('Could not get data from URL %s'%(progtemp))
    a = 1
    #raise a warning here



dicts = {}
dicts['seeing_max'] = {}
dicts['transparency_min'] = {}
dicts['elevation_min'] = {}
dicts['moondist_min'] = {}
dicts['skybright_max'] = {}

ordering = {}
for key in dicts.keys():
    ordering[key] = {}

dlookup = {}
dlookup[0] = 'seeing_max'
dlookup[1] = 'transparency_min'
dlookup[2] = 'elevation_min'
dlookup[3] = 'moondist_min'
dlookup[4] = 'skybright_max'

for line in obstemp_data:
    if '=' in line:
        _line = line.split('\t')[1].replace('\n','').split('#')[0]
        key = _line.split(':')[0]
        for i in np.arange(len(key)):
            if key[i] != '_':
                dicts[dlookup[i]][key[i]] = _line.split('=')[1].split()[0]
#                ordering[dlookup[i]][float(_line.split('=')[1].split()[0])] = key[i]
                ordering[dlookup[i]][_line.split('=')[1].split()[0]] = key[i]
_dicts = {}
for d in dicts.keys():
    _dicts[d] = OrderedDict()
    _vals = ordering[d].keys()
    _vals.sort()
    #highest -> lowest
    _vals.reverse()
    for v in _vals:
        _dicts[d][ordering[d][v]] = dicts[d][ordering[d][v]]

dicts = _dicts
        
seeing_dict = OrderedDict(zip(dicts['seeing_max'].values(),
                              dicts['seeing_max'].keys()))

transparency_dict = OrderedDict(zip(dicts['transparency_min'].values(),
                                    dicts['transparency_min'].keys()))

#print dicts['transparency_min']

elevation_dict = OrderedDict(zip(dicts['elevation_min'].values(),
                                 dicts['elevation_min'].keys()))

moon_dict = OrderedDict(zip(dicts['moondist_min'].values(),
                            dicts['moondist_min'].keys()))


skybright_dict = OrderedDict(zip(dicts['skybright_max'].values(),
                                 dicts['skybright_max'].keys()))

airmasses = list(map(lambda x: '{:1.1f}'.format(1.0/(cos(radians(90-float(x))))),
                 elevation_dict.keys()))
airmass_dict = OrderedDict(zip(dicts['elevation_min'].keys(), airmasses))





seeing = widgets.FloatSlider(
            value=0.9,
            min=0.7,
            max=3.0,
            step=0.1,
            description='Seeing Max:',
            disabled=False,
            continuous_update=True,
            orientation='horizontal',
            readout=True,
            readout_format='.1f',
    style = style, layout=widgets.Layout(width='30%'))

transparency = widgets.RadioButtons(
            options=transparency_dict,
            value='A',
            description='Transparency Min:',
            disabled=False, style=style)

elevation = widgets.RadioButtons(
            options=elevation_dict,
            value='A',
            description='Elevation Min:',
            disabled=False, style = style)

moon = widgets.RadioButtons(
            options=moon_dict,
            value='E',
            description='Moon Dist. Min:',
            disabled=False, style = style)

skybright = widgets.RadioButtons(
            options=skybright_dict,
            value='B',
            description='Sky Brightness Max:',
            disabled=False, style = style)

plate = widgets.ToggleButtons(
    options=['PLATE_A', 'PLATE_B','LIFU'],
    description='Plate:',
    disabled=False,
    button_style='', # 'success', 'info', 'warning', 'danger' or ''
    value='PLATE_A'
)

hour_angle_range = widgets.FloatRangeSlider(
    value=[-2.0, 2.0],
    min=-5.0,
    max=5.0,
    step=0.1,
    description='HA range:',
    disabled=False,
    continuous_update=False,
    orientation='horizontal',
    readout=True,
    readout_format='.1f',
)

def print_obs(seeing, transparency, elevation, moon, skybright):
    global obstemp
    obstemp = seeing_dict[str(seeing)] + transparency + elevation + moon + skybright
    print('OBSTEMP: {0}'.format(obstemp))

def print_airmass(elevation):
    print('Max. Airmass: {0}'.format(airmass_dict[elevation]))

def getobstemp(seeing, transparency, elevation, moon, skybright):
    obstemp = seeing_dict[seeing] + transparency_dict[transparency] + \
              elevation_dict[elevation] + moon_dict[moon] + \
              skybright_dict[skybright]
    return obstemp


out = widgets.interactive_output(print_obs,
                                 {'seeing': seeing, 'transparency': transparency,
                                  'elevation': elevation, 'moon': moon, 'skybright':skybright})

airmass_out = widgets.interactive_output(print_airmass, {'elevation': elevation})

printouts = widgets.VBox([airmass_out, out])

row1 = widgets.HBox([seeing, transparency, elevation])
row2 = widgets.HBox([moon, skybright, printouts])

obs_controls = widgets.VBox([row1, row2])

vals = [("MOS","low","VPH1", "VPH1"),
        ("MOS","high","VPH2", "VPH2"),
        ("MOS","high","VPH2", "VPH3"),
        ("LIFU", "low", "VPH1", "VPH1"),
        ("LIFU", "high", "VPH2", "VPH2"),
        ("LIFU", "high", "VPH2", "VPH3"),
        ("mIFU", "low",  "VPH1", "VPH1"),
        ("mIFU", "high", "VPH2", "VPH2"),
        ("mIFU", "high", "VPH2", "VPH3")]

junk = """
for line in progtemp_data:
    if '=' in line:
        _line = line.split('\t')[1].replace('\n','').split('#')[0]
        key = _line.split(':')[0]
        for i in xrange(len(key)):
            if key[i] != '_':
                dicts[dlookup[i]][key[i]] = _line.split('=')[1].split()[0]

"""

N = widgets.RadioButtons(
            options=OrderedDict((('Type: {0:<10s} Resolution: {1} Red Arm: {2} Blue Arm: {3}'.format(*vals[0]), '1'),
                                 ('Type: {0:<10s} Resolution: {1} Red Arm: {2} Blue Arm: {3}'.format(*vals[1]), '2'),
                                 ('Type: {0:<10s} Resolution: {1} Red Arm: {2} Blue Arm: {3}'.format(*vals[2]), '3'),
                                 ('Type: {0:<10s} Resolution: {1} Red Arm: {2} Blue Arm: {3}'.format(*vals[3]), '4'),
                                 ('Type: {0:<10s} Resolution: {1} Red Arm: {2} Blue Arm: {3}'.format(*vals[4]), '5'),
                                 ('Type: {0:<10s} Resolution: {1} Red Arm: {2} Blue Arm: {3}'.format(*vals[5]), '6'),
                                 ('Type: {0:<10s} Resolution: {1} Red Arm: {2} Blue Arm: {3}'.format(*vals[6]), '7'),
                                 ('Type: {0:<10s} Resolution: {1} Red Arm: {2} Blue Arm: {3}'.format(*vals[7]), '8'),
                                 ('Type: {0:<10s} Resolution: {1} Red Arm: {2} Blue Arm: {3}'.format(*vals[8]), '9')
                                )),
            value='1',
            description='Instrument Configuration:',
            disabled=False, style = style, layout=widgets.Layout(width='100%'))

orb_vals = []
orb_vals.append((1,1620,"both"))
orb_vals.append((2,720,"both"))
orb_vals.append((3,420,"both"))

orb_vals.append((1,3420,"both"))
orb_vals.append((2,1620,"both"))
orb_vals.append((3,1020,"both"))
orb_vals.append((4,720,"both"))
orb_vals.append((6,420,"both"))

orb_vals.append((1,7020,"both"))
orb_vals.append((2,3420,"both"))
orb_vals.append((4,1620,"both"))
orb_vals.append((6,1020,"both"))
orb_vals.append((8,720,"both"))
orb_vals.append((12,420,"both"))

orb_vals = {}
orb_vals['000'] = (1,1620,"both")
orb_vals['044'] = (2,720,"both")
orb_vals['066'] = (3,420,"both")

orb_vals['100'] = (1,3420,"both")
orb_vals['122'] = (2,1620,"both")
orb_vals['133'] = (3,1020,"both")
orb_vals['144'] = (4,720,"both")
orb_vals['166'] = (6,420,"both")

orb_vals['200'] = (1,7020,"both")
orb_vals['211'] = (2,3420,"both")
orb_vals['222'] = (4,1620,"both")
orb_vals['233'] = (6,1020,"both")
orb_vals['244'] = (8,720,"both")
orb_vals['266'] = (12,420,"both")


OBlenkenobi = widgets.Dropdown(options=['30m','1hr','2hr'],layout=widgets.Layout(flex="1 1 auto"),value='1hr',description='OB length')
#Use the widget, Luke

ORB_selector = {}
ORB_selector['30m'] = OrderedDict((
                ("{0} Exposures Exp.Time={1} Arm={2}".format(*orb_vals['000']), "000"),
                ("{0} Exposures Exp.Time={1} Arm={2}".format(*orb_vals['044']), "044"),
                ("{0} Exposures Exp.Time={1} Arm={2}".format(*orb_vals['066']), "066")))

ORB_selector['1hr'] = OrderedDict((
                ("{0} Exposures Exp.Time={1} Arm={2}".format(*orb_vals['100']), "100"),
                ("{0} Exposures Exp.Time={1} Arm={2}".format(*orb_vals['122']), "122"),
                ("{0} Exposures Exp.Time={1} Arm={2}".format(*orb_vals['133']), "133"),
                ("{0} Exposures Exp.Time={1} Arm={2}".format(*orb_vals['144']), "144"),
		("{0} Exposures Exp.Time={1} Arm={2}".format(*orb_vals['166']), "166")))

ORB_selector['2hr'] = OrderedDict((
                ("{0} Exposures Exp.Time={1} Arm={2}".format(*orb_vals['200']), "200"),
                ("{0} Exposures Exp.Time={1} Arm={2}".format(*orb_vals['211']), "211"),
                ("{0} Exposures Exp.Time={1} Arm={2}".format(*orb_vals['222']), "222"),
                ("{0} Exposures Exp.Time={1} Arm={2}".format(*orb_vals['233']), "233"),
                ("{0} Exposures Exp.Time={1} Arm={2}".format(*orb_vals['244']), "244"),
		("{0} Exposures Exp.Time={1} Arm={2}".format(*orb_vals['266']), "266")))

init = OBlenkenobi.value

ORB_inter = widgets.RadioButtons(options=ORB_selector[init],value='133',
                                       description='Exposures:',
                                       style=style, layout=widgets.Layout(width='100%'))

ORB = ORB_inter


def set_orb(key):
    ORB_inter.options = ORB_selector[key]

def print_orb(key):
    #this is pretty much a dummy function
    return

interlink2 = widgets.interactive(print_orb, key=ORB_inter, layout=widgets.Layout(flex="1 1 auto"))
interlink1 = widgets.interactive(set_orb, key=OBlenkenobi, layout=widgets.Layout(flex="2 1 auto"))

I = widgets.RadioButtons(
            options=OrderedDict((
                ('1x', '1'), ('2x', '2'), ('3x', '3'),
                ('4x', '4'))),
            value='1',
            description='Spectral Binning:',
            disabled=False, style = style, layout=widgets.Layout(width='30%'))


I = widgets.RadioButtons(
            options=OrderedDict((
                ('1x', '1'), ('2x', '2'), ('3x', '3'),
                ('4x', '4'))),
            value='1',
            description='Spectral Binning:',
            disabled=False, style = style, layout=widgets.Layout(flex="0 1 auto"))



X = widgets.IntSlider(
            value=0,
            min=0,
            max=7,
            step=1,
            description='Contextual OB Actions:',
            disabled=False,
            continuous_update=True,
            orientation='horizontal',
            readout=True,
            readout_format='d',
    style = style, layout=widgets.Layout(width='30%'))

grouped = widgets.Checkbox(
    value=False,
    description='Chain observations:',
    disabled=False
)

def printprog(N, ORB, I, X, grouped):
    global progtemp
    global blue_vph
    global red_vph
    global resolution
    global obstype
    if X > 1:
        progtemp = N + ORB + I + '.' + str(X)
        if grouped:
            progtemp += '+'
    else:
        progtemp = N + ORB + I
    print('PROGTEMP: {0}'.format(progtemp))
    #print orb_vals,ORB
    config_index = int(N)-1
    obstype, resolution, red_vph, blue_vph = vals[config_index]


out = widgets.interactive_output(printprog,
                                 {'N': N, 'ORB': ORB_inter,
                                  'I': I, 'X': X, 'grouped':grouped})

row1 = widgets.HBox([N])

row2 = widgets.HBox(layout=widgets.Layout(width='100%',display='inline-flex',flex_flow='row wrap', justify_content='space-between'))
elem = []
elem.append(interlink1)
elem.append(interlink2)
elem.append(I)
row2.children=[i for i in elem]

row3 = widgets.HBox([X, grouped, out])

prog_controls = widgets.VBox([row1, row2, row3],layout=widgets.Layout(width='100%'))


row4 = widgets.HBox([plate, hour_angle_range])
misc_controls = widgets.VBox([row4])
