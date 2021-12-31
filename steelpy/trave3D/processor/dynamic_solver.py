# 
# Copyright (c) 2009-2021  steelpy
#
# Python stdlib imports
#
# package imports
#import numpy as np
#from scipy.linalg import eig
from steelpy.trave3D.processor.operations import zeros
from steelpy.trave3D.processor.jacobi import jacobi
#
#
#
# --------------------
# eigenvalue solver
# --------------------
#
def eigen(neqmas, ibandm, ivib, neq, iband, x, mtx):
    """
    Solve eigenvalue problem using jacobi rotations
    """
    # stiffness
    #a = np.zeros((neq, iband), dtype = np.float64, order = 'F')
    a = zeros(neq, neq)
    # mass
    #b = np.zeros((neq, iband), dtype = np.float64, order = 'F')
    b = zeros(neq, neq)
    #
    rtol = 1e-12
    iwidth = neq
    # x  = rdstff(iwidth, iband)
    # reassign to full form  
    for i in range(neq):
        jmax = min(iband, neq - i)
        for j in range(jmax):
            a[i][i + j] = x[i][j]
            a[i + j][i] = x[i][j]
    #
    #ilog = open('stadyn.log','a') 
    # 
    if ivib == 2 :
        x = mtx  # rdmass(neqmas, iband)
        # reassign to [b] matrix 
        if ibandm == 1 :
            for i in range(neq):
                b[i][0] = x[i][0]
        else :
            for i in range(neq):
                jmax = min(ibandm, neq - i)
                
                for j in range(jmax):
                    b[i][i + j] = x[i][j]
                    b[i + j][i] = x[i][j]
        #
        #ilog.write("@@ eign:  reloaded [k] [m]  ok \n")
        print("@@ eign:  reloaded [k] [m]  ok")
    elif ivib == 1 :
        iwidth =  neq
        x = mtx # rdgeom(iwidth, iband)
        # reassign to [b] 
        for i in range(neq):
            jmax = min(iband, neq - i)
            for j in range(jmax):
                b[i][i + j] = x[i][j]
                b[i + j][i] = x[i][j]
        #ilog.write("@@ eign: raloaded [k] [g] ok  ok \n")
        print("@@ eign: raloaded [k] [g] ok")
    #
    ibandm = neqmas
    nsmax = 15
    _eigv, _x = jacobi(a, b, neq, rtol, nsmax, ibandm, x, neqmas)
    _eigv, _x = eigsrt(_eigv, _x, neq)
    #
    #eigv, v = eig(a, b)
    #print ("{:}".format(eigv))
    #print ("====")
    #print ("{:}".format(v))
    #v1 = v.T
    #x = np.dot(a,v) - eigv*np.dot(b,v)
    #x = x.real
    #eigv = eigv.real
    #print ("====")
    #z = np.dot(v.T, np.dot(b, v))
    #z = np.dot(np.dot(v.T, b), v)
    #print ("{}".format(z))
    #print ("====")
    #x = np.array([v[:,i]/np.sqrt(abs(z[i,i])) for i in range(neq)])
    #x = x.T
    #print("{}".format(x))
    #print  np.divide(v,b)
    #print [v[:,i]/math.sqrt(abs(b[i,i])) for i in range(neq)]
    
    #print"----------"
    #print eigv
    
    #print np.linalg.eigvals(a)
    #print np.linalg.eigvals(b)
    
    # eigv = jacobi(a, b, eigv, d, neq, rtol, nsmax, ibandm, x, neqmas)
    
    print("@@ No of sweeps".format(nsmax))
    #ilog.write("@@ No of sweeps".format(nsmax))
    #
    #eigv, x = eigsrt(eigv, x)
    #print"----------"
    #print eigv
    
    # rewind isnp
    #isnp = open('stadyn.snp','w')
    #for i in range(neq):
    #    isnp.write("{: 1.6e} ".format(float(_eigv[i])))
    #    for j in range(neq):
    #        isnp.write("{: 1.6e} ".format(float(_x[j][i])))
    #    isnp.write(" \n") 
    # L84:
    #isnp.close()
    #ilog.close()
    #
    return _x, _eigv
#
def eigsrt(eigv, x, neq):
    """
    Sort the eigenvalues in ascending order
    """
    #
    for i in range(neq-1):
        k = i
        p = eigv[i]
        # search for lowest value
        for j in range(i+1, neq):
            if abs(eigv[j]) < abs(p) :
                k = j
                p = eigv[j]
        # L11:
        # re-arrange vectors 
        if k != i :
            eigv[k] = eigv[i]
            eigv[i] = p
            
            for j in range(neq):
                p = x[j][i]
                x[j][i] = x[j][k]
                x[j][k] = p
            # L12:
    # L13:
    #
    return eigv, x
#
#
# --------------------           
#  transient 
# --------------------                      
#
def trnsient(stf, mass, load, disp, vel, acc, fmag, olddis,
             wk, damp, jbc, loadin, maxnode):
    """
    Transient  analysis by newmark time integration
    """
    # parameter( sigma=0.5e0, alpha=0.25e0)
    # Dim loadin(maxnode*3) as double
    ildin = maxnode*3
    #
    # get things # ready for _time integration loop
    # open report.txt for app# end as #1
    print (" ")
    print (" --> ")
    # read (*,*) deltat, npt, iprcnt
    if npt > ildin :
        print("@@ _time steps npt > load size ", npt," > ",ildin)
        print("@@ _time steps npt > load size ", npt," > ",ildin)
        npt = ildin - 1
    # end if
    print(deltat, npt, iprcnt," ::dt no pts print")
    startt = 0
    # endt = real(npt) * deltat
    #
    # set integration constants for newmark method 
    a0 = 1.0/(alpha*deltat*deltat)
    a1 = sigma/(alpha*deltat)
    a2 = 1.0/(alpha*deltat)
    a3 = 1.0/(alpha*2.0) - 1.0
    a4 = sigma/alpha - 1.0
    a5 = (sigma/alpha - 2.0) *0.5*deltat
    a6 = (1.0 - sigma) *deltat
    a7 = sigma*deltat
    #
    # read stiffness, mass & load 
    iwidth = iband
    rdstff(stf, iwidth)
    iwidth = ibandm
    rdmass(mass,iwidth)
    rdload(fmag)
    print("@@ reloaded [k] [m] {p} ok")
    print("@@ reloaded [k] [m] {p} ok")
    #
    print("@@ dammath.ping coeffs: ", dampkk, dampmm)
    if dampkk > 0.0 or dampmm > 0.0 : 
        idamp = 1
    else : 
        idamp = 0
    #
    if idamp == 1 :
        for i in range(neq):
            damp[i] = dampkk * stf[i][1] + dampmm * mass[i][1] + dampcc
    
    # form effective stiffness matrix
    if ibandm == 1 :
        for i in range(neq):
            stf[i][1] += a0 * mass[i][1]
    else :
        for i in range(neq):
            for j in range(iband):
                stf[i][j] += a0 * mass[i][j]
    
    if idamp == 1 :
        for i in range(neq):
            stf[i][1] += a1 * damp[i]
    #
    # decompose effective stiffness matrix 
    ier1 = 0
    udu(stf, neq, iband, ier1)
    if ier1 == 0 :
        print("error: zero diagonal term")
        #close #1:exit def
        sys.exit("error: zero diagonal term")
    #
    print("@@ udu: mults & divs",ier1)
    print("@@ udu: mults & divs",ier1)
    #
    # input load history from a file and interpolate
    getload(npt, deltat, loadin, ildin, ilog)
    #
    print(" ")
    print("choose nodal output: ")
    nout = 0
    # L27:
    nout += 1
    print("type: node no |  dof  |  rate  | <<0 0 0 for No end>>")
    print(" --> ")
    # read(*,*) inode, idof, irate
    jdof = (inode-1) *6 + idof
    inod[nout][1] = jbc[jdof]
    inod[nout][2] = irate
    print(inode, idof, irate," :: node dof rate")
    print("@@ eqn no ",inod[nout][1])
    if inode != 0 :  
        print (" goto L27")
    nout -= 1
    #
    # start _time integration loop
    print("@@ beginning transient analysis")
    print("@@ beginning transient analysis")
    #
    # initialize tines, displacement, velocity, accln
    print("@@ initial disps. vels and accels set=0")
    for i in range(neq):
        disp[i] = 0.0e0
        vel[i]  = 0.0e0
        acc[i]  = 0.0e0
        # read(*,*) acc(i)
    #L620:
    
    forc = loadin[1]
    nodeout(_time,forc,disp,vel,acc,neq,nout,inod,idyn)
    kount = 0
    print(" ")
    # big _time loop
    for itime in range(1, npt):
        if mod((itime - 1) ,10) == 0 : 
            print(itime-1)
        _time = real(itime-1) * deltat
        kount += 1
        #
        # save displacements
        for i in range(neq):
            olddis[i] = disp[i]
        #
        # fmag says where the load is applied
        # done this way in case distributed load applied
        for i in range(neq):
            load[i] = loadin[itime] * fmag[i]
        #
        # form effective load vector
        if idamp == 1 :
            for i in range(neq):
                atermc = a1*disp[i] + a4*vel[i] + a5*acc[i]
                load[i] += atermc * damp[i]
        #
        # the effective acceleration
        for i in range(neq):
            aterm   = a0*disp[i] + a2*vel[i] + a3*acc[i]
            disp[i] = aterm
        #
        if ibandm == 1 :
            for i in range(neq):
                load[i] += disp[i] * mass[i,1]
        else :
            abband(mass, disp, load, neq, ibandm)
        # 
        # solve for new displacements : udu al# ready obtained  
        #                             : do back - defstitution
        ier1 = 0
        bak(stf,load,neq,iband,wk,ier1)
        #
        # obtain new velocities, accelerations
        for i in range(neq):
            oldvel = vel[i]
            oldacc = acc[i]
            disp[i] = wk[i]
            acc[i] = a0*(disp[i] - olddis[i]) - a2*oldvel - a3*oldacc
            vel[i] = oldvel + a6*oldacc + a7*acc(i)
        #
        # print out nodal results
        if iprcnt == kount :
            kount = 0
            forc = loadin[itime]
            nodeout(_time, forc, disp, vel, acc, neq,nout, inod, idyn)
    #
    # L30:
    # next itime
    # bottom of _time integration loop   
    #
    # 122   format(1x,8(g12.5, 1x) )        
    # 123   format(1x, i3, 1x, 6(g12.5,1x) )
    #
    # close #1
    # end def
#
def abband(matrix, vecin, vecout, neq, iband):
    """
    Multiplies [ banded ]{vector} = {vector}
    """
    for i in range(neq):
        jlim = max(0, i - iband + 1)
        for j in range(jlim, i):
            val = vecin[j]
            vecout[i] += val * matrix[j][i - j + 1]
        #
        jlim = min(iband, neq - i + 1)
        for j in range(1, jlim):
            val = vecin[i + j - 1]
            vecout[i] += val * matrix[i][j]
    #
    return vecout
#
def getload(npt, dt, loadin, nmax, ilog):
    """
    Gets applied load history by interpolation from input
    """
    iload = 24
    #
    # read data file
    # open report.txt for app# end as #1
    print(" ")
    # read(*,"(1a40) ") fy12
    print(fy12," ::filename")
    #
    # open(unit=iload, file = fy12)
    # rewind iload
    #
    # adjust _times and interpolate 
    #
    _time = 0.0
    # read(iload,*,# end=240) t1, f1
    tzero = t1
    t1 = t1 - tzero
    loadin[0] = f1
    #
    k = 0
    _time += dt
    for i in range(1, 4100):
        t0 = t1
        f0 = f1
        # read(iload,*,# end=220) t1, f1
        t1 -= tzero
        # L214:
        if k == npt : 
            continue  # print " goto L220"
        if t1  >  _time :
            xk4 = (f1 - f0) / (t1 - t0) * (_time - t0) + f0
            k += 1
            _time += dt
            loadin[k] = xk4
            print (" goto L214")
    #L210:
    #
    #L220:
    print("@@ # of force points # read ", i-1)
    print("@@ # of force points # read ", i-1)
    #     pad reminder with zeroes
    for i in range(k+1, npt):
        loadin[i] = 0.0
    #L230:
    #
    # close (iload)
    # close #1:exit def
    #
    # L240:
    print("@@ no data in load file")
#
# 
#
