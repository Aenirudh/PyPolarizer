""" POLARIZATION ELLIPSE MODE """

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider, Button

from py_pol import degrees, np
from py_pol.stokes import Stokes
from py_pol.mueller import Mueller
from py_pol.drawings import draw_poincare_sphere, draw_on_poincare

import pdb

Stokes_points2 = []

""" CLASSES FOR WAVEPLATES, IN/OUT SIGNAL VECTORS, AND VISUALIZING """
class WAVEPLATES():
	def __init__(self, radii, turns, rotation, matrix_name, quarter=False):
		self.name = matrix_name
		self.radii = radii
		self.turns = turns
		self.plate_angle = rotation
		self.Mueller = Mueller(matrix_name)
		if (quarter): self.Mueller.quarter_waveplate(angle=self.plate_angle*degrees)
		else: self.Mueller.half_waveplate(angle=self.plate_angle*degrees)

	def ROTATE_PLATE(self, val):
		if (self.name=='QWP_A'):
			self.Mueller.quarter_waveplate(angle=val*degrees)

		if (self.name=='HWP'):
			self.Mueller.half_waveplate(angle=val*degrees)

		if (self.name=='QWP_B'):
			self.Mueller.quarter_waveplate(angle=val*degrees)

		TWEAK_SIGNAL(False, True, True)

class STOKES():
	def __init__(self, s_name):
		self.a = 1
		self.b = 0
		self.ang = 0
		self.ph = 0
		self.sVECTOR = Stokes(s_name)
		#self.sVECTOR.linear_light(angle=45*degrees)
		self.sVECTOR.elliptical_light(a=self.a, b=self.b, phase=self.ph*degrees, angle=self.ang*degrees, pol_degree=1)

	def E_FIELD_PARAMETERS(self):
		E0x, E0y, E0_unpol = self.sVECTOR.parameters.amplitudes()
		delay = self.sVECTOR.parameters.delay()

		angles = np.linspace(0, 360 * degrees, 256)
		Ex = E0x * np.cos(angles)
		Ey = E0y * np.cos(angles + delay)

		E_unpolarized_x = E0_unpol * np.cos(angles)
		E_unpolarized_y = E0_unpol * np.sin(angles)

		size_unpol = np.sqrt(E_unpolarized_x**2 + E_unpolarized_y**2).max()

		return Ex, Ey, E_unpolarized_x, E_unpolarized_y, size_unpol, delay

	def UPDATE_ELLIPSE_a(self, val):
		self.a = val
		self.sVECTOR.elliptical_light(a=self.a, b=self.b, phase=self.ph*degrees, angle=0*degrees, pol_degree=1)
		TWEAK_SIGNAL(True, True, True)

	def UPDATE_ELLIPSE_b(self, val):
		self.b = val
		self.sVECTOR.elliptical_light(a=self.a, b=self.b, phase=self.ph*degrees, angle=0*degrees, pol_degree=1)
		TWEAK_SIGNAL(True, True, True)

	def UPDATE_ELLIPSE_ang(self, val):
		self.ang = val
		self.sVECTOR.elliptical_light(a=self.a, b=self.b, phase=self.ph*degrees, angle=self.ang*degrees, pol_degree=1)
		TWEAK_SIGNAL(True, True, True)

	def UPDATE_ELLIPSE_p(self, val):
		self.ph = val
		self.sVECTOR.elliptical_light(a=self.a, b=self.b, phase=self.ph*degrees, angle=0*degrees, pol_degree=1)
		TWEAK_SIGNAL(True, True, True)

class VISUALIZE():
	def __init__(self, fig, plot_num, title):
		self.ax = fig.add_subplot(plot_num)
		self.title = title

	def ELLIPSIFY(self, SIGNAL):
		Ex, Ey, E_unpolarized_x, E_unpolarized_y, size_unpol, delta = SIGNAL.E_FIELD_PARAMETERS()
		self.ax.clear()
		phase_shift = round(np.degrees(delta), 2)
		label = 'PHASE SHIFT: ' + str(phase_shift) + ' deg'
		self.ax.plot(Ex, Ey, 'k', lw=2, label=label)

		self.ax.grid(True)
		self.ax.axis('equal')
		self.ax.axis('square')
		self.ax.set_title(self.title)
		self.ax.set_xlabel('$E_x$', fontsize=22)
		self.ax.set_ylabel('$E_y$', fontsize=22)
		self.ax.set_xlim(-1, 1)
		self.ax.set_ylim(-1, 1)
		self.ax.legend()

def TWEAK_SIGNAL(TWEAK_INPUT=False, TWEAK_OUTPUT=False, UPDATE_PLOT=False):
	if (TWEAK_INPUT):
		global INPUT
		global VIS_IN
		VIS_IN.ELLIPSIFY(INPUT)

	if (TWEAK_OUTPUT):
		global OUTPUT
		global QWP_A
		global QWP_B
		global HWP
		global Stokes_points2

		""" FULL TRANSFORMATION """
		OUTPUT.sVECTOR = QWP_B.Mueller*HWP.Mueller*QWP_A.Mueller*INPUT.sVECTOR
		
		""" IDEAL CASES """
		#OUTPUT.sVECTOR = QWP_A.Mueller*INPUT.sVECTOR
		#OUTPUT.sVECTOR = HWP.Mueller*INPUT.sVECTOR
		#OUTPUT.sVECTOR = HWP.Mueller*QWP_A.Mueller*INPUT.sVECTOR
		
		Stokes_points2.append(OUTPUT.sVECTOR)

		if (UPDATE_PLOT):
			global VIS_OUT
			VIS_OUT.ELLIPSIFY(OUTPUT)

def CREATE_POINCARE_IN(event):
	global INPUT
	INPUT.sVECTOR.draw_poincare()
	plt.show()

def CREATE_POINCARE_OUT(event):
	global OUTPUT
	OUTPUT.sVECTOR.draw_poincare()
	plt.show()

def CREATE_POINCARE_TRANS(event):
	global fig
	global Stokes_points2
	ax2, fig2=draw_poincare_sphere(stokes_points=Stokes_points2 , kind='line', color='b', label='rotating quarter wave', filename='poincare3.png')
	Stokes_points2.clear()
	plt.show()

"""

STARTING POLARIZATION CONTROLLER SIMULATION THE WAVEPLATES, STOKES VECTOR, MUELLER MATRICES ARE INITIALIZED BELOW.

"""

""" SETUP FIGURE OBJECT """
fig = plt.figure()
plt.subplots_adjust(left=0.2, bottom=None, right=None, top=1.2, wspace=0.5, hspace=None)

""" WAVEPLATES """
QWP_A = WAVEPLATES(2.1, 1, 0, 'QWP_A', True)
QWP_B = WAVEPLATES(2.1, 1, 0, 'QWP_B', True)
HWP = WAVEPLATES(2.1, 2, 0, 'HWP')

""" IN/OUT STOKES """
INPUT = STOKES('Source')
OUTPUT = STOKES('Output')
TWEAK_SIGNAL(False, True, False)

""" GENERATING POLARIZATION ELLIPSES FOR IN/OUT """
VIS_IN = VISUALIZE(fig, 121, 'INPUT SOP')
VIS_IN.ELLIPSIFY(INPUT)

VIS_OUT = VISUALIZE(fig, 122, 'OUTPUT SOP')
VIS_OUT.ELLIPSIFY(OUTPUT)

""" SLIDERS """
SCOLORa = 'darkgoldenrod'
SCOLORb = 'goldenrod'

""" Input Ellipse Parameters """
SLIDER_IN_a_POS = plt.axes([0.2, 0.24, 0.65, 0.02], facecolor=SCOLORa)
SLIDER_IN_a = Slider(SLIDER_IN_a_POS, 'Input Amp x', -1, 1, valinit=1, valstep=0.1)

SLIDER_IN_b_POS = plt.axes([0.2, 0.20, 0.65, 0.02], facecolor=SCOLORa)
SLIDER_IN_b = Slider(SLIDER_IN_b_POS, 'Input Amp y', -1, 1, valinit=0, valstep=0.1)

SLIDER_IN_p_POS = plt.axes([0.2, 0.16, 0.65, 0.02], facecolor=SCOLORa)
SLIDER_IN_p = Slider(SLIDER_IN_p_POS, 'Input Ellipse Phase', -180, 180, valinit=0, valstep=5)

SLIDER_IN_ang_POS = plt.axes([0.2, 0.12, 0.65, 0.02], facecolor=SCOLORa)
SLIDER_IN_ang = Slider(SLIDER_IN_ang_POS, 'Input Ellipse Angle', -180, 180, valinit=0, valstep=5)
#SLIDER_IN_p = Slider(SLIDER_IN_p_POS, 'Input Ellipse Phase', 0, 180, valinit=0)

""" Waveplate Parameters """
SLIDER_A_POS = plt.axes([0.2, 0.08, 0.65, 0.02], facecolor=SCOLORb)
SLIDER_QWP_A = Slider(SLIDER_A_POS, 'QWP A Angle', -180, 180, valinit=0, valstep=5)

SLIDER_B_POS = plt.axes([0.2, 0.04, 0.65, 0.02], facecolor=SCOLORb)
SLIDER_HWP = Slider(SLIDER_B_POS, 'HWP Angle', -180, 180, valinit=0, valstep=0.5)

SLIDER_C_POS = plt.axes([0.2, 0.00, 0.65, 0.02], facecolor=SCOLORb)
SLIDER_QWP_B = Slider(SLIDER_C_POS, 'QWP B Angle', -180, 180, valinit=0, valstep=1)

""" Poincare Button """
p_axeso = plt.axes([0.8, 0.28, 0.1, 0.04])
POINMEo = Button(p_axeso, 'SPHERE OUT', color='lightgoldenrodyellow', hovercolor='0.975')

p_axesi = plt.axes([0.1, 0.28, 0.1, 0.04])
POINMEi = Button(p_axesi, 'SPHERE IN', color='lightgoldenrodyellow', hovercolor='0.975')

""" Event Check """
SLIDER_IN_a.on_changed(INPUT.UPDATE_ELLIPSE_a)
SLIDER_IN_b.on_changed(INPUT.UPDATE_ELLIPSE_b)
SLIDER_IN_ang.on_changed(INPUT.UPDATE_ELLIPSE_ang)
SLIDER_IN_p.on_changed(INPUT.UPDATE_ELLIPSE_p)

SLIDER_QWP_A.on_changed(QWP_A.ROTATE_PLATE)
SLIDER_HWP.on_changed(HWP.ROTATE_PLATE)
SLIDER_QWP_B.on_changed(QWP_B.ROTATE_PLATE)

""" POINCARE SPOT """
POINMEo.on_clicked(CREATE_POINCARE_OUT)
POINMEi.on_clicked(CREATE_POINCARE_IN)

""" POINCARE TRANSFORMATION HISTORY """
#POINMEo.on_clicked(CREATE_POINCARE_TRANS)

plt.show()
