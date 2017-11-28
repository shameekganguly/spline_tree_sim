import numpy as np
import matplotlib.pyplot as plt

def dsinc(x):
	try:
		ret_arr = np.array([(y*np.cos(y) - np.sin(y))/(y**2) if abs(y) > 1e-5 else 0 for y in x])
	except TypeError:
		ret_arr = (x*np.cos(x) - np.sin(x))/(x**2) if abs(x) > 1e-5 else 0
	return ret_arr

def sinc(x):
	try:
		ret_arr = np.array([np.sin(y)/y if abs(y) > 1e-5 else 0 for y in x])
	except TypeError:
		ret_arr = np.sin(x)/x if abs(x) > 1e-5 else 0
	return ret_arr

def dpsi_dalp(psi, alp, beta):
	return -np.cos(beta)*np.sin(alp)/np.cos(psi)

def dpsi_dbeta(psi, alp, beta):
	return -np.sin(beta)*np.cos(alp)/np.cos(psi)

def dgam_dalp(psi, alp, beta):
	return -2.0*dpsi_dalp(psi, alp, beta)

def dgam_dbeta(psi, alp, beta):
	return -2.0*dpsi_dbeta(psi, alp, beta)

def px_dalp(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	return (s**2)/l*dsinc(gam*s/l)*dgam_dalp(psi, alp, beta)

def py_dalp(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	tb = np.tan(beta)
	sa = np.sin(alp)
	ca = np.cos(alp)
	ret_val = s*s/2.0/l*((sinc(s*gam/2.0/l)**2) + gam*s/l*sinc(s*gam/2.0/l)*dsinc(s*gam/2.0/l))*tb/np.sqrt(tb**2 + sa**2)*dgam_dalp(psi, alp, beta)
	ret_val -= s*s/2.0/l*gam*(sinc(s*gam/2.0/l)**2)*tb*sa*ca/np.power(tb**2 + sa**2, 1.5)
	return ret_val

def pz_dalp(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	tb = np.tan(beta)
	sa = np.sin(alp)
	ca = np.cos(alp)
	ret_val = -s*s/2.0/l*((sinc(s*gam/2.0/l)**2) + gam*s/l*sinc(s*gam/2.0/l)*dsinc(s*gam/2.0/l))*sa/np.sqrt(tb**2 + sa**2)*dgam_dalp(psi, alp, beta)
	ret_val -= s*s/2.0/l*gam*(sinc(s*gam/2.0/l)**2)*(tb**2)*ca/np.power(tb**2 + sa**2, 1.5)
	return ret_val

def px_dbeta(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	return (s**2)/l*dsinc(gam*s/l)*dgam_dbeta(psi, alp, beta)

def py_dbeta(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	tb = np.tan(beta)
	sa = np.sin(alp)
	ca = np.cos(alp)
	cb = np.cos(beta)
	ret_val = s*s/2.0/l*((sinc(s*gam/2.0/l)**2) + gam*s/l*sinc(s*gam/2.0/l)*dsinc(s*gam/2.0/l))*tb/np.sqrt(tb**2 + sa**2)*dgam_dbeta(psi, alp, beta)
	ret_val += s*s/2.0/l*gam*(sinc(s*gam/2.0/l)**2)*(sa**2)/(cb**2)/np.power(tb**2 + sa**2, 1.5)
	return ret_val

def pz_dbeta(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	tb = np.tan(beta)
	sa = np.sin(alp)
	ca = np.cos(alp)
	cb = np.cos(beta)
	ret_val = -s*s/2.0/l*((sinc(s*gam/2.0/l)**2) + gam*s/l*sinc(s*gam/2.0/l)*dsinc(s*gam/2.0/l))*sa/np.sqrt(tb**2 + sa**2)*dgam_dbeta(psi, alp, beta)
	ret_val += s*s/2.0/l*gam*(sinc(s*gam/2.0/l)**2)*(sa*tb)/(cb**2)/np.power(tb**2 + sa**2, 1.5)
	return ret_val

def Jv_s(s, alp, beta):
	return np.array(
		[[px_dalp(s, alp, beta), px_dbeta(s, alp, beta)],
		[py_dalp(s, alp, beta), py_dbeta(s, alp, beta)],
		[pz_dalp(s, alp, beta), pz_dbeta(s, alp, beta)],
		]
	)

def Rxx_dalpha(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	tb = np.tan(beta)
	sa = np.sin(alp)
	ca = np.cos(alp)
	cb = np.cos(beta)
	ret_val = -s/l*np.sin(s*gam/l)*dgam_dalp(psi, alp, beta)
	return ret_val

def Rxy_dalpha(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	tb = np.tan(beta)
	sa = np.sin(alp)
	ca = np.cos(alp)
	cb = np.cos(beta)
	ret_val = s/l*np.cos(s*gam/l)*tb/np.sqrt(tb**2 + sa**2)*dgam_dalp(psi, alp, beta)
	ret_val -= np.sin(s*gam/l)*(ca*sa*tb)/np.power(tb**2 + sa**2, 1.5)
	return ret_val

def Rxz_dalpha(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	tb = np.tan(beta)
	sa = np.sin(alp)
	ca = np.cos(alp)
	cb = np.cos(beta)
	ret_val = -s/l*np.cos(s*gam/l)*sa/np.sqrt(tb**2 + sa**2)*dgam_dalp(psi, alp, beta)
	ret_val -= np.sin(s*gam/l)*ca*(tb**2)/np.power(tb**2 + sa**2, 1.5)
	return ret_val

def Ryx_dalpha(s, alp, beta):
	return -Rxy_dalpha(s, alp, beta)

def Rzx_dalpha(s, alp, beta):
	return -Rxz_dalpha(s, alp, beta)

def Ryy_dalpha(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	tb = np.tan(beta)
	sa = np.sin(alp)
	ca = np.cos(alp)
	cb = np.cos(beta)
	ret_val = 2.0*ca*sa*(1 - np.cos(s*gam/l))/(tb**2 + sa**2)
	ret_val -= s/l*np.sin(s*gam/l)*dgam_dalp(psi, alp, beta)
	ret_val *= tb**2/(tb**2 + sa**2)
	return ret_val

def Rzz_dalpha(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	tb = np.tan(beta)
	sa = np.sin(alp)
	ca = np.cos(alp)
	cb = np.cos(beta)
	ret_val = 2.0*ca*(tb**2)/sa*(np.cos(s*gam/l) - 1)/(tb**2 + sa**2)
	ret_val -= s/l*np.sin(s*gam/l)*dgam_dalp(psi, alp, beta)
	ret_val *= sa**2/(tb**2 + sa**2)
	return ret_val

def Ryz_dalpha(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	tb = np.tan(beta)
	sa = np.sin(alp)
	ca = np.cos(alp)
	cb = np.cos(beta)
	ret_val = s/l*np.sin(s*gam/l)*sa*tb/(tb**2 + sa**2)*dgam_dalp(psi, alp, beta)
	ret_val += (1 - np.cos(s*gam/l))*ca*tb*(tb**2 - sa**2)/np.power(tb**2 + sa**2, 2)
	return ret_val

def Rzy_dalpha(s, alp, beta):
	return Ryz_dalpha(s, alp, beta)

def Rxx_dbeta(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	tb = np.tan(beta)
	sa = np.sin(alp)
	ca = np.cos(alp)
	cb = np.cos(beta)
	ret_val = -s/l*np.sin(s*gam/l)*dgam_dbeta(psi, alp, beta)
	return ret_val

def Rxy_dbeta(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	tb = np.tan(beta)
	sa = np.sin(alp)
	ca = np.cos(alp)
	cb = np.cos(beta)
	ret_val = -s/l*np.cos(s*gam/l)*tb/np.sqrt(tb**2 + sa**2)*dgam_dbeta(psi, alp, beta)
	ret_val -= np.sin(s*gam/l)*(sa**2)/(cb**2)/np.power(tb**2 + sa**2, 1.5)
	return ret_val

def Ryx_dbeta(s, alp, beta):
	return -Rxy_dbeta(s, alp, beta)

def Rxz_dbeta(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	tb = np.tan(beta)
	sa = np.sin(alp)
	ca = np.cos(alp)
	cb = np.cos(beta)
	ret_val = s/l*np.cos(s*gam/l)*sa/np.sqrt(tb**2 + sa**2)*dgam_dbeta(psi, alp, beta)
	ret_val -= np.sin(s*gam/l)*sa*tb/(cb**2)/np.power(tb**2 + sa**2, 1.5)
	return ret_val

def Rzx_dbeta(s, alp, beta):
	return -Rxz_dbeta(s, alp, beta)

def Ryy_dbeta(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	tb = np.tan(beta)
	sa = np.sin(alp)
	ca = np.cos(alp)
	sb = np.sin(beta)
	cb = np.cos(beta)
	ret_val = 2.0*(sa**2)/(cb*sb)*(np.cos(s*gam/l)-1)/(tb**2 + sa**2)
	ret_val -= s/l*np.sin(s*gam/l)*dgam_dbeta(psi, alp, beta)
	ret_val *= tb**2/(tb**2 + sa**2)
	return ret_val

def Rzz_dbeta(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	tb = np.tan(beta)
	sa = np.sin(alp)
	ca = np.cos(alp)
	sb = np.sin(beta)
	cb = np.cos(beta)
	ret_val = 2.0*tb/(cb**2)*(1-np.cos(s*gam/l))/(tb**2 + sa**2)
	ret_val -= s/l*np.sin(s*gam/l)*dgam_dbeta(psi, alp, beta)
	ret_val *= sa**2/(tb**2 + sa**2)
	return ret_val

def Ryz_dbeta(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	tb = np.tan(beta)
	sa = np.sin(alp)
	ca = np.cos(alp)
	cb = np.cos(beta)
	ret_val = s/l*np.sin(s*gam/l)*sa*tb/(tb**2 + sa**2)*dgam_dbeta(psi, alp, beta)
	ret_val += (1 - np.cos(s*gam/l))*sa/(cb**2)*(sa**2 - tb**2)/np.power(tb**2 + sa**2, 2)
	return ret_val

def Rzy_dbeta(s, alp, beta):
	return Ryz_dbeta(s, alp, beta)

def lp_alpha(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	ret_val = s*s/2.0/l*((sinc(s*gam/2.0/l)**2) + s*gam/l*sinc(s*gam/2.0/l)*dsinc(s*gam/2.0/l))*dgam_dalp(psi, alp, beta)
	ret_val *= s*s/2.0/l*gam*(sinc(s*gam/2.0/l)**2)
	return ret_val

def lp_beta(s, alp, beta):
	l = 1.0
	psi = np.arcsin(np.cos(alp)*np.cos(beta))
	gam = np.pi - 2.0*psi
	ret_val = s*s/2.0/l*((sinc(s*gam/2.0/l)**2) + s*gam/l*sinc(s*gam/2.0/l)*dsinc(s*gam/2.0/l))*dgam_dbeta(psi, alp, beta)
	ret_val *= s*s/2.0/l*gam*(sinc(s*gam/2.0/l)**2)
	return ret_val

# test
# alp = np.pi/10.0
# beta = np.pi/10.0
# s = 0.732
# print Jv_s(s, alp, beta)

alp = np.linspace(-np.pi*2.0/3.0, np.pi*2.0/3.0, 20)
beta = np.linspace(-np.pi*2.0/3.0, np.pi*2.0/3.0, 20)
grad = np.zeros((20,20))
for i in range(20):
	for j in range(20):
		grad[i,j] = lp_beta(0.5, alp[i], beta[j])

plt.imshow(grad, cmap='hot', interpolation='nearest')
plt.colorbar()
plt.show()
