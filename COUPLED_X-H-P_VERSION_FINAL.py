#Total coupling simulation
import math;import numpy as np;from scipy.integrate import odeint;import matplotlib.pyplot as plt;from numpy import *;import numpy, scipy.optimize
from scipy.signal import argrelextrema
import matplotlib.lines as lines
from matplotlib.lines import Line2D 
row = 1030;g = 9.81;H = 1.0;h = 25;L = 35.6
#rob = 748#962.57
rob = 962.57
lb = 18.6
lmda = 2.885*L

le = L-lb;De = 2.0;lc = lb-De/2;lf = lb-De;B = 9.91;Ae = De*B;Ac=B*(lb-lf)
d1 = 1.56;d2 = d1+De;ld = (d1+d2)/2#2.56

#########Surging
AB = B*(le+lf)
mB = rob*AB*d1
Vs = B*L*d1
ax = row*Vs
#############Heaving
Vz = B*L*d1
az = 0.7*row*Vz
bzr = 0.25*row*Vz*math.sqrt(g/L)
############Pitching
#Iy=5.28e7#
Iy=mB *(L**2+d1**2)/12
#HB=math.sqrt(12*5.28e7/mB-L**2)
#7.48e-10
print(Vz,mB,Iy)

################Linear Theory
wT = math.sqrt(lmda*(2*math.pi)/(g*math.tanh(2*math.pi*h/lmda)))
k = 2*math.pi/lmda
U = 0.0
omgae = 2*math.pi/wT+U*k;

Crx = (1+d1/d2)*(0.1667*math.pi/omgae-0.5)*row*Vs*math.sqrt(g/L)

print(wT,lmda,omgae)
####################################################################################
# solve model with ODEINT

def model_coupling(w,t):
    eta0  = H/2*math.cos(omgae*t)  
    deta0 = -H*omgae/2*math.sin(omgae*t)
    eta1  = H/2*math.cos(omgae*t-k*le)  
    deta1 = -H*omgae/2*math.sin(omgae*t-k*le) 
    eta2  = H/2*math.cos(omgae*t+k*lb)
    deta2 = -H*omgae/2*math.sin(omgae*t+k*lb) 
    eta3  = H/2*math.cos(omgae*t+k*lf)
    deta3 = -H*omgae/2*math.sin(omgae*t+k*lf)

    SHYP = (math.sinh(k*h-k*d1)-math.sinh(k*h-k*d2))/math.cosh(k*h)
    Ve = g*H/(2*omgae*De)*SHYP*math.cos(omgae*t-k*le)
    Pe = row*g*H/(2*k*De)*(SHYP+k*De)*math.cos(omgae*t-k*le)
    #########Surging
    SHYP1 = (math.sinh(k*eta1-k*h)-math.sinh(k*h-k*d1))/math.cosh(k*h)
    SHYP2 = (math.sinh(k*eta2+k*h)-math.sinh(k*h-k*d2))/math.cosh(k*h)
    Fwx = row*g*B*(eta1/k*SHYP1-eta1**2/2-eta2/k*SHYP2+eta2**2/2)
    C1x = 1/(mB+ax)*(-Crx-2*row*Ae*Ve)
    CFx = 1/(mB+ax)*(row*Ae) 
    Sx =  1/(mB+ax)*(row*Ae*Ve**2+Fwx)
    ##########Heaving
    Fw=row*g*H*B/(2*k)*( math.cosh(k*h-k*d2)/math.cosh(k*h)*( math.sin(omgae*t+k*lb)-math.sin(omgae*t-k*le) )
                         +math.sin(omgae*t+k*lf)-math.sin(omgae*t-k*le) )
    C1z= 1/(mB+az)*(-bzr)
    C2z= 1/(mB+az)*(-row*g*AB)
    C3z= 1/(mB+az)*(2*U*az)
    C4z= 1/(mB+az)*(-1)*(0.5*row*g*B*(lf**2-le**2)-U*bzr)
    C5z= 1/(mB+az)*(2*row*Ae*Ve)
    CFz= 1/(mB+az)*(-row*Ae)
    Sz=  1/(mB+az)*(-row*Ae*Ve**2+Fw)
    ###########Pitching
    Mw=row*g*B/k**2*( math.cosh(k*h-k*d2)/math.cosh(k*h)*( eta2+k*lb/omgae*deta2-eta1-k*le/omgae*deta1 )
                      +eta3-k*lf/omgae*deta3-eta1-k*le/omgae*deta1 )   
    C1thta= 1/(Iy+az*L**2/12)*(-bzr*L**2/12) 
    C2thta= 1/(Iy+az*L**2/12)*(-row*g*(le**3+lb**3)/3) 
    C3thta= 1/(Iy+az*L**2/12)*(-row*g*B*(lb**2-le**2)/3) 
    C4thta= 1/(Iy+az*L**2/12)*(lc-ld)*(2*row*Ae*Ve) 
    CFthta= 1/(Iy+az*L**2/12)*(lc-ld)*(-row*Ae) 
    Sthta=  1/(Iy+az*L**2/12)*((lc-ld)*(-row*Ae*Ve**2)+Mw)
    
    x = w[0]
    dx = w[1]
    z = w[2]
    dz = w[3]
    thta = w[4]
    dthta = w[5]
    
    wdot = [[],[],[],[],[],[]]
    wdot[0] = dx
    wdot[1] = C1x*dx + CFx*dx**2 + Sx
    wdot[2] = dz
    wdot[3] = C1z*dz+ C2z*z+ C3z*dthta+ C4z*thta+ C5z*dx+ CFz*dx**2 + Sz
    wdot[4] = dthta
    wdot[5] = C1thta*dthta+ C2thta*thta+ C3thta*z+ C4thta*dx+ CFthta*dx**2 + Sthta

    return wdot

t = np.linspace(0,500,5000)
y = odeint(model_coupling,[0,0,0,0,0,0],t)

#----------------------------------------------------------------------

t0=3950
ft=10.2
#_pc=pc[t0:t0+int(wT*ft)]

_t=t[t0:t0+int(wT*ft)]

etagc  = H/2*np.cos(omgae*t)
_etagc = etagc[t0:t0+int(wT*ft)]

detagc = -H*omgae/2*np.sin(omgae*t)
_detagc = detagc[t0:t0+int(wT*ft)]

etabow  = H/2*np.cos(omgae*t+k*lb)
_etabow = etabow[t0:t0+int(wT*ft)]

detabow = -H*omgae/2*np.sin(omgae*t+k*lb)
_detabow = detabow[t0:t0+int(wT*ft)]

etastr  = H/2*np.cos(omgae*t-k*le)
_etastr = etastr[t0:t0+int(wT*ft)]

detastr = -H*omgae/2*np.sin(omgae*t-k*le)
_detastr = detastr[t0:t0+int(wT*ft)]

c = lmda/wT
#----------------------------------------------------------------------------------
mlb = argrelextrema(_etabow, np.greater)[0] #array of indexes of the locals maxima
ylb = [_etabow[i] for i in mlb] #array of max values
ylb=ylb[0]
mle = argrelextrema(_etastr, np.greater)[0] #array of indexes of the locals maxima
yle = [_etastr[i] for i in mle] #array of max values
yle=yle[0]
mlgc = argrelextrema(_etagc, np.greater)[0] #array of indexes of the locals maxima
ylgc = [_etagc[i] for i in mlgc] #array of max values
ylgc=ylgc[0]
lbe=c*(_t[mle]-_t[mlb])
lbgc=c*(_t[mlgc]-_t[mlb])
#----------------------------------------------------------------------
fig = plt.figure()
plt.title('Instantenous free-surface elevation ($\eta$) from the still-water level at center of gravity (CG), bow and stern \n'
          + r'(Wave speed, c=%.2f m/s)'%c,fontsize=14)
ax = fig.add_subplot(1,1,1)                                                      

ax.plot(_t,_etagc,'k-', label=r"$\eta_{GC}, m/s$")
plt.annotate(r'$\eta_{ _{GC}}$', xy=(_t[-1], _etagc[-1]), xytext=(30,10),
         textcoords='offset points',fontsize=16, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.3',color='k'))

ax.plot(_t,_etabow,'k--', label=r"$\eta_{BOW}, m/s$")
plt.annotate(r'$\eta_{ _{BOW}}$', xy=(_t[0], _etabow[0]), xytext=(-35,-20), 
                  textcoords='offset points',fontsize=16, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='k'))

ax.plot(_t,_etastr,'k-.', label=r"$\eta_{STR}, m/s$")
plt.annotate(r'$\eta_{ _{STR}}$', xy=(_t[-1], _etastr[-1]), xytext=(35,-15), 
         textcoords='offset points',fontsize=16, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.3',color='k'))

ax.plot(_t[mlb], ylb, 'rs')
ax.plot(_t[mlgc], ylgc, 'ks')
ax.plot(_t[mle], yle, 'bs')

ax.axvline(_t[mlb], linestyle='--', color='0.75')
ax.axvline(_t[mlgc], linestyle='--', color='0.75')
ax.axvline(_t[mle], linestyle='--', color='0.75')

x1 = [_t[mlgc],_t[mlb]]
y1 = [ylgc,ylb]
line = Line2D(x1, y1)
ax.add_line(line)
plt.annotate(r'$L_{CG-BOW}=%.2f m$'%lbgc, xy=(mean(x1), mean(y1)), xytext=(-75,-40),
         textcoords='offset points',fontsize=11, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.2',color='k'))
x2 = [_t[mle],_t[mlb]]
y2 = [yle+0.03,ylb+0.03]
line = Line2D(x2, y2)
ax.add_line(line)
plt.annotate(r'$L_{BOW-STR}=%.2f m$'%lbe, xy=(mean(x2), mean(y2)), xytext=(0,-50),
         textcoords='offset points',fontsize=11, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2',color='k'))
plt.xlabel('Time,s')
plt.xlim(min(_t),max(_t))
ax.set_ylabel(r"$Velocities,\;\; (\frac{m}{s}) \;\;and\;\;  Elevations, \;\;(m)$",fontsize=11, color='k')
plt.grid()
plt.show()
#----------------------------------------------------------------------
#c_z = (diff(sign(diff(_z))) < 0).nonzero()[0] + 1 # local max
#------------------------------------------------------------------------
t0=3910
ft=10.2
_t=t[t0:t0+int(wT*ft)]

etastr  = H/2*np.cos(omgae*t-k*le)
_etastr = etastr[t0:t0+int(wT*ft)]

detastr = -H*omgae/2*np.sin(omgae*t-k*le)
_detastr = detastr[t0:t0+int(wT*ft)]

SHYP = (math.sinh(k*h-k*d1)-math.sinh(k*h-k*d2))/math.cosh(k*h)
_Ve = np.empty_like(_t)
_Ve = g*H/(2*omgae*De)*SHYP*np.cos(omgae*_t-k*le)
_Pe = row*g*H/(2*k*De)*(SHYP+k*De)*np.cos(omgae*_t-k*le)

fig = plt.figure()
plt.title('Wave-induced dynamic pressure ($P_e$), water (particle) velocity at the mouth ($V_e$) and\n'
          + r'free-surface elevation at stern ($\eta_{STR}$)',fontsize=14)
ax = fig.add_subplot(1,1,1)

ax.plot(_t,_etastr,'k-', label=r"$\eta_{STR}, m$")
plt.annotate(r'$\eta_{ _{STR}}$', xy=(_t[0], _etastr[0]), xytext=(-30,5),
         textcoords='offset points',fontsize=16, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.95',color='k'))

ax.plot(_t,_Ve,'k--', label=r"$V_e, m/s$")
plt.annotate(r'$V_e$', xy=(_t[0], _Ve[0]), xytext=(-20,-25),
         textcoords='offset points',fontsize=16, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.95',color='k'))

ax2 = ax.twinx()
ax2.plot(_t,_Pe,'k.', label=r"$P_e, kPa$")
plt.annotate(r'$P_e$', xy=(_t[-1], _Pe[-1]), xytext=(20,30),
         textcoords='offset points',fontsize=16, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='k'))

plt.xlabel('Time,s')
plt.xlim(min(_t),max(_t))
ax.set_ylabel(r"$Velocities,\;\; (\frac{m}{s}) \;\;and\;\;  Elevations, \;\;(m)$",fontsize=11, color='k')
ax2.set_yticks([-10000,-5000,0,5000,10000])
ax2.set_yticklabels([-10.,-5.,0,5.,10.])
ax2.set_ylabel(r"Wave-induced dynamic pressure, Pa",fontsize=11, color='k')
ax.grid()
plt.xlabel('Time,s')
plt.show()
#----------------------------------------------------------------------
_ts=t[t0:t0+int(10*wT*ft)]
_Ves = g*H/(2*omgae*De)*SHYP*np.cos(omgae*_ts-k*le)
x=y[:,0]

vx=y[:,1]

_vxs=vx[t0:t0+int(10*wT*ft)]

_xs=x[t0:t0+int(10*wT*ft)]

_vzcs=_Ves-_vxs
#----------------------------------------------------------------------
f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(_ts,_xs,'k-', label=r"$x, m/s$")
plt.annotate(r'$x$', xy=(_ts[-1], _xs[-1]), xytext=(25,0),
         textcoords='offset points',fontsize=16, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.95',color='k'))
axarr[0].set_title('Surging$\;\;$ motion, $x$',fontsize=14)
axarr[0].set_ylabel(r"$Surging\;\; motion\;\;(m)$",fontsize=11, color='k')

axarr[0].grid()

axarr[1].plot(_ts,_vxs,'k-', label=r"$\dot(x), m/s$")
plt.annotate(r'$\dot{x}$', xy=(_ts[0], _vxs[0]), xytext=(-45,-10),
         textcoords='offset points',fontsize=16, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.5',color='k'))

axarr[1].plot(_ts,_Ves,'k--', label=r"$V_e, m/s$")
plt.annotate(r'$V_e$', xy=(_ts[0], _Ves[0]), xytext=(-45,10),
         textcoords='offset points',fontsize=16, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.5',color='k'))

axarr[1].plot(_ts,_vzcs,'k-.', label=r"$(V_e-\dot{x}), m/s$")
plt.annotate(r'$(V_e-\dot{x})$', xy=(_ts[0], _vzcs[0]), xytext=(-60,15),
         textcoords='offset points',fontsize=14, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='k'))

axarr[1].set_title('Surging velocity ($\dot{x}$) with wave-induced water velocity ($V_e$) at the mouth and\n'
                   + r' absolute velocity of the internal free-surface, $\dot{z}_c=(V_e-\dot{x})$')
axarr[1].grid()
ax2 = axarr[1].twinx()
ax2.set_ylabel(r"$Velocities,\;\; (\frac{m}{s}) $",fontsize=11, color='k')
ax2.set_yticks([])
plt.xlabel('Time,s')
plt.xlim(min(_ts),max(_ts))
plt.show()
#----------------------------------------------------------------------
_vx=vx[t0:t0+int(wT*ft)]

_vzc=_Ve-_vx
_vzc_mn=mean(_vzc)
a = np.array(_t)
a.fill(_vzc_mn)
_vzc=_vzc-a

vz=y[:,3]
_vz=vz[t0:t0+int(wT*ft)]

vthta=y[:,5]
_vthta=vthta[t0:t0+int(wT*ft)]

_Vc= (_Ve -_vx)*Ae/Ac - _vz - _vthta*lc

def fit_sin(tt, yy):
    '''Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase",
        "offset", "freq", "period" and "fitfunc"'''
    tt = numpy.array(tt)
    yy = numpy.array(yy)
    ff = numpy.fft.fftfreq(len(tt), (tt[1]-tt[0]))   # assume uniform spacing
    Fyy = abs(numpy.fft.fft(yy))
    guess_freq = abs(ff[numpy.argmax(Fyy[1:])+1])   # excluding the zero frequency "peak", which is related to offset
    guess_amp = numpy.std(yy) * 2.**0.5
    guess_offset = numpy.mean(yy)
    guess = numpy.array([guess_amp, 2.*numpy.pi*guess_freq, 0., guess_offset])
    def sinfunc(t, A, w, p, c):  return A * numpy.sin(w*t + p) + c
    popt, pcov = scipy.optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    A, w, p, c = popt
    f = w/(2.*numpy.pi)
    fitfunc = lambda t: A * numpy.sin(w*t + p) + c
    fit_para[0]=A
    fit_para[1]=w
    fit_para[2]=p
    fit_para[3]=c

    return {"amp": A, "omega": w, "phase": p, "offset": c, "freq": f, "period": 1./f,
            "fitfunc": fitfunc, "maxcov": numpy.max(pcov), "rawres": (guess,popt,pcov),
            "fit_para": fit_para}

N, amp, omega, phase, offset= int(ft*10), 0., 0., 0.0,0.
tt = _t
yy = _vzc

fit_para=[0.0,0.0,0.0,0.0]
res = fit_sin(tt, yy)
print( "Amplitude=%(amp)s, Angular freq.=%(omega)s, phase=%(phase)s, offset=%(offset)s, Max. Cov.=%(maxcov)s" % res )
print(fit_para)
Amp=fit_para[0]
Angfreq=fit_para[1]
phase=fit_para[2]
offset=fit_para[3]

vzc_fit= Amp*sin(Angfreq*_t+phase)+offset*1.0

_vzowc=(_Ve-_vx) - _vz - _vthta*lc

#----------------------------------------------------------------------
fig = plt.figure()
plt.title('Pitching velocity, ($\dot{\Theta}$), heaving velocity, ($\dot{z}$),\n'
                   + r' absolute ($\dot{z}_c=V_e-\dot{x}$) and relative ($V_c$) velocities of the internal free-surface',
          fontsize=14)
ax = fig.add_subplot(1,1,1)

ax.plot(_t,_vz,'k-', label=r"$\dot(z), m/s$")
plt.annotate(r'$\dot{z}$', xy=(_t[0], _vz[0]), xytext=(-45,-10),
         textcoords='offset points',fontsize=16, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.5',color='k'))

ax.plot(_t,_vthta,'k--', label=r"$\dot{\theta}, rd/s$")
plt.annotate(r'$\dot{\theta}$', xy=(_t[0], _vthta[0]), xytext=(-45,0),
         textcoords='offset points',fontsize=16, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.5',color='k'))

ax.plot(_t,vzc_fit,'k-.', label=r"$(V_e-\dot{x}), m/s$") 
plt.annotate(r"$\dot{z}_{c}= %.2f*sin(%.2f*t%.2f)$"
             %(Amp,Angfreq,phase), xy=(_t[35], vzc_fit[35]), xytext=(120,-10),
         textcoords='offset points',fontsize=11, ha='center', va='bottom',color='k',
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='k'))

ax2 = ax.twinx()
ax2.plot(_t,_Vc,'k.-', label=r"$V_c, m/s$")
plt.annotate(r'$V_c$', xy=(_t[-1], _Vc[-1]), xytext=(40,-30),rotation=90,
         textcoords='offset points',fontsize=16, ha='center', va='bottom',color='k',
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.7',color='k'))
ax2.set_ylabel(r"$=(V_e-\dot{x})\;\frac{A_e}{A_c}-\dot{z}-lc\;\dot{\theta} \;\; [\frac{m}{s}] $",fontsize=14, color='k')
plt.xlabel('Time,s')
plt.xlim(min(_t),max(_t))
plt.xlabel('Time,s')
ax.grid()
plt.show()
#----------------------------------------------------------------------
t1=3903
_t1=t[t1:t1+int(wT*ft)]
z=y[:,2]
_z=z[t1:t1+int(wT*ft)]
thta=y[:,4]
_thta=thta[t1:t1+int(wT*ft)]
zc_fit=- Amp/Angfreq*cos(Angfreq*_t1+phase)+offset*0.0
_zc=zc_fit

b_z = (diff(sign(diff(_z))) > 0).nonzero()[0] + 1 # local min
c_z = (diff(sign(diff(_z))) < 0).nonzero()[0] + 1 # local max

b_thta = (diff(sign(diff(_thta))) > 0).nonzero()[0] + 1 # local min
c_thta = (diff(sign(diff(_thta))) < 0).nonzero()[0] + 1 # local max

b_zc = (diff(sign(diff(_zc))) > 0).nonzero()[0] + 1 # local min
c_zc = (diff(sign(diff(_zc))) < 0).nonzero()[0] + 1 # local max


_zowc = _zc - _z - lc*_thta#pcMc=BME/Dc*OWC_Rel_Disp

b_zowc = (diff(sign(diff(_zowc))) > 0).nonzero()[0] + 1 # local min
c_zowc = (diff(sign(diff(_zowc))) < 0).nonzero()[0] + 1 # local max

#_thtadeg=_thta*180/pi
AmpAngfreq=-Amp/Angfreq
phasecos=phase+math.pi/2

fig = plt.figure()
plt.title('Instanteneous heaving, ($z$), pitching effect, ($lc\;\Theta$) and oscillating water column, ($OWC$)\n'
                   + r'(Wave speed, c=%.2f m/s)'%c,fontsize=14)
ax = fig.add_subplot(1,1,1)


ax.plot(_t1,_z,'k-', label=r"$z,\;\;m$")
plt.annotate(r'$z$', xy=(_t1[0], _z[0]), xytext=(-35,-30),
         textcoords='offset points',fontsize=14, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.5',color='k'))

ax.plot(_t1,_zc,'k--', label=r"$z_c,\;\;m$")
plt.annotate(r"$z_c= %.2f*sin(%.2f*t%.2f)$"%(AmpAngfreq,Angfreq,phasecos),
             xy=(_t1[20], _zc[20]), xytext=(100,-10),
         textcoords='offset points',fontsize=11, ha='center', va='bottom',color='k',
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.5',color='k'))

ax2 = ax.twinx()
ax2.plot(_t1,lc*_thta,'k-.', label=r"$lc\;\theta,\;\;m$")
plt.annotate(r'$lc\;\Theta$', xy=(_t1[-1], lc*_thta[-1]), xytext=(35,10),rotation=0,
         textcoords='offset points',fontsize=16, ha='center', va='bottom',color='k',
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.35',color='k'))

ax2.plot(_t1,_zowc,'k.-', label=r"$OWC$") 
plt.annotate(r"$OWC$", xy=(_t1[-1], _zowc[-1]), xytext=(40,30),rotation=90,
         textcoords='offset points',fontsize=12, ha='center', va='bottom',color='k',
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.5',color='k'))
ax2.set_ylabel(r"$=z_c-z-lc\;\theta \;\; [m]$",fontsize=14, color='k')
                      
ax.legend(loc='best')
plt.xlabel('Time,s')
plt.xlim(min(_t1),max(_t1))
plt.grid()
plt.show()
#---------------------------------------------------------------------------

norm_z= np.empty_like(_z)
norm_thta= np.empty_like(_thta)
norm_zc= np.empty_like(_zc)
norm_zowc= np.empty_like(_zowc)

i = 0

while i < len(_t1):
    
    if _zc[i]<0:
        norm_zc[i] = _zc[i]/abs(_zc[b_zc])
    else:
        norm_zc[i] = _zc[i]/_zc[c_zc]
    
    if _z[i]<0:
        norm_z[i] = _z[i]/abs(_z[b_z])
    else:
        norm_z[i] = _z[i]/_z[c_z]
        
    if _thta[i]<0:
        norm_thta[i] = _thta[i]/abs(_thta[b_thta])
    else:
        norm_thta[i] = _thta[i]/_thta[c_thta]
        
    if _zowc[i]<0:
        norm_zowc[i] = _zowc[i]/abs(_zowc[b_zowc])
    else:
        norm_zowc[i] = _zowc[i]/_zowc[c_zowc]  
    i+= 1

#----------------------------------------------------------------------
mlb = argrelextrema(norm_thta, np.greater)[0] #array of indexes of the locals maxima
ylb = [norm_thta[i] for i in mlb] #array of max values
ylb=ylb[0]

mlgc = argrelextrema(norm_z, np.greater)[0] #array of indexes of the locals maxima
ylgc = [norm_z[i] for i in mlgc] #array of max values
ylgc=ylgc[0]

lbgc=c*(_t1[mlb]-_t1[mlgc])

mlc = argrelextrema(_zc, np.greater)[0] #array of indexes of the locals maxima
ylc = [_zc[i] for i in mlb] #array of max values
ylc=ylc[0]

lcgc=c*(_t1[mlc]-_t1[mlgc])

fig = plt.figure()
plt.title('Normalized: pitching, ($\Theta/\Theta_{max}$), heaving, ($z/z_{max}$), absolute and relative ($z_c,\;OWC$) internal free-surface displacement\n'
                   + r'(Wave speed, c=%.2f m/s)'%c,fontsize=14)
ax = fig.add_subplot(1,1,1)

ax.plot(_t1,norm_z,'k-', label=r"$\frac{z}{z_{max}}$")
plt.annotate(r'$\frac{z}{z_{max}}$', xy=(_t1[0], norm_z[0]), xytext=(-45,-10),
         textcoords='offset points',fontsize=14, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='k'))

ax.plot(_t1,norm_thta,'k:', label=r"${\Theta}{\Theta_{max}}$")
plt.annotate(r'$\frac{\Theta}{\Theta_{max}}$', xy=(_t1[0], norm_thta[0]), xytext=(-45,-20),
         textcoords='offset points',fontsize=14, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='k'))

ax.plot(_t1,norm_zc,'k--', label=r"$\frac{z_c}{z_c{_{max}}}$")
plt.annotate(r'$\frac{z_c}{z_c{_{max}}}$', xy=(_t1[0], norm_zc[0]), xytext=(-45,20),
         textcoords='offset points',fontsize=14, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.5',color='k'))

ax.plot(_t1,norm_zowc,'k.-', label=r"$(\frac{OWC}{OWC_{max}}$") 
plt.annotate(r"$\frac{OWC}{OWC_{max}}$", xy=(_t1[0], norm_zowc[0]), xytext=(-55,20),
         textcoords='offset points',fontsize=12, ha='center', va='bottom',color='k',
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='k'))
#ax.legend(loc='best')
ax.plot(_t1[mlb], ylb, 'rs')
ax.plot(_t1[mlgc], ylgc, 'ks')
ax.plot(_t1[mlc], ylc, 'gs')

ax.axvline(_t1[mlb], linestyle='--', color='0.75')
ax.axvline(_t1[mlgc], linestyle='--', color='0.75')
ax.axvline(_t1[mlc], linestyle='--', color='0.75')

x1 = [_t1[mlgc],_t1[mlb]]
y1 = [0.8,0.8]
line = Line2D(x1, y1)
ax.add_line(line)
plt.annotate(r'$\Delta L_(\theta_{max}-z_{max})=%.2f m$'%lbgc, xy=(mean(x1), mean(y1)), xytext=(0,-50),
         textcoords='offset points',fontsize=11, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2',color='k'))

x2 = [_t1[mlgc],_t1[mlc]]
y2 = [0.9,0.9]
line = Line2D(x2, y2)
ax.add_line(line)
plt.annotate(r'$\Delta L_(z_{c_{max}}-z_{max})=%.2f m$'%lcgc, xy=(mean(x2), mean(y2)), xytext=(0,-50),
         textcoords='offset points',fontsize=11, ha='center', va='bottom',color='k',
         bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.3),
         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2',color='k'))

x_ticks = np.arange(t1/10.0, t1/10.0*wT, wT/2.0)
x_labels = [0,r"$\pi/\omega_e$",r"$2\pi/\omega_e$"]                       

ax.set_xticks(x_ticks)                                                           
ax.set_xticklabels(x_labels)
ax.set_yticks([-1.0,-0.5,0,0.5,1.0])
ax.set_yticklabels([-1.0,-0.5,0,0.5,1.0])
plt.xlim(min(_t1),max(_t1))
plt.ylim(-1,1)
plt.xlabel("time (s)")
plt.grid()
plt.show()
#----------------------------------------------------------------------

##fig = plt.figure()
##plt.title('Instantenous of normalized heaving, pitching, OWC responses\n' + r'$(\omega_e$=%f rd/s)'%omgae)
##ax = fig.add_subplot(1,1,1)                                                      
##ax.plot(_t1,norm_z, 'k-', label=r"$z/z_{max}$")                                                          
##ax.plot(_t1,norm_thta,'k:', label=r"$\theta/\theta_{max}$")
##ax.plot(_t1,norm_zc,'k--', label=r"$zc/zc_{max}$")
##ax.plot(_t1,norm_zowc,'k.:', label=r"$zc/zc_{max}$")
##
##ax.legend(loc='best')
##
##x_ticks = np.arange(t1/10.0, t1/10.0*wT, wT/2.0)
##x_labels = [0,r"$\pi/\omega_e$",r"$2\pi/\omega_e$"]                       
##
##ax.set_xticks(x_ticks)                                                           
##ax.set_xticklabels(x_labels)
##ax.set_yticks([-1.0,-0.5,0,0.5,1.0])
##ax.set_yticklabels([-1.0,-0.5,0,0.5,1.0])
##plt.xlim(min(_t1),max(_t1))
##plt.ylim(-1,1)
##plt.xlabel("time (s)")
##plt.grid()
##plt.show()

   
##Ac=B*(lb-lf)
##Dc=4
##BME=141919
##pcMc=BME*_zowc/Dc
##_vzc=vzc_fit
##PMc=abs(pcMc*Ac*_vzowc)
##
##fig = plt.figure()
##plt.title('Pressure in the air chamber and converted power response\n' + r'$(\omega_e$=%f rd/s)'%omgae)
##ax1 = fig.add_subplot(1,1,1)                                                      
##ax1.plot(_t,pcMc/1e3,'b:')
##ax2 = ax1.twinx()
##ax2.plot(_t,PMc/1e3,'g--')
##ax1.set_ylabel(r"$pc_{Mc}, kPa$", color='b')
##ax2.set_ylabel(r"$P_{Mc}, kw$", color='g')
##ax1.legend(loc='best')
##plt.xlabel('Time,s')
##plt.xlim(min(_t),max(_t))
##plt.grid()
##plt.show()
##
