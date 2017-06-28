# A script to map the Drude topology file into the MPID 
# Jing Huang (jing.huang.ff@gmail.com)
# June, 2017
# some data structures for handling rtf IO are based on the source code of openmm python API

import math, sys, copy

def numeq(a,b):
    if math.fabs(a-b) < 1E-9:
        return True
    else:
        return False

def numeqloose(a,b):
    if math.fabs(a-b) < 1E-4:
        return True
    else:
        return False

def dipole(q,c,scale):
    u = [0.0,0.0,0.0]
    nm=len(q)
    for i in range(3):
        for m in range(nm):
            u[i] += q[m]*c[m][i]
        u[i] = u[i]*scale
    return u

def quadrupole(q,c,scale):
    Q = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
    nm=len(q)
    for i in range(3):
        for j in range(3):
            for k in range(nm):
                d2=0.0
                if ( i == j ):
                    d2 = (c[k][0]*c[k][0]) + (c[k][1]*c[k][1]) + (c[k][2]*c[k][2])
                Q[i][j] = Q[i][j] + q[k] * ( 3. * c[k][i]*c[k][j] - d2 )
            Q[i][j] = 0.5*Q[i][j]*scale
    return Q

def octupole(q,c,scale):
    O = [[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]],[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]],[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]]
    nm=len(q)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for m in range(nm):
                    d2 = (c[m][0]*c[m][0]) + (c[m][1]*c[m][1]) + (c[m][2]*c[m][2])
                    d3 = 0.0
                    if (i == j):
                        d3 += c[m][k]
                    if (j == k):
                        d3 += c[m][i]
                    if (i == k):
                        d3 += c[m][j]
                    O[i][j][k] = O[i][j][k] + q[m] * (15*c[m][i]*c[m][j]*c[m][k]-3*d2*d3)
                O[i][j][k] = O[i][j][k]*scale/6
    return O

def hexadecapole(q,c,scale):
    H=[[[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]],[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]],[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]] for i in range(3)]
    nm=len(q)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(nm):
                        d2 = (c[m][0]*c[m][0]) + (c[m][1]*c[m][1]) + (c[m][2]*c[m][2])
                        H[i][j][k][l] = H[i][j][k][l] + q[m] * (105*c[m][i]*c[m][j]*c[m][k]*c[m][l]  \
                        -15*d2*(c[m][i]*c[m][j]*delta(k,l)+c[m][i]*c[m][k]*delta(j,l)+c[m][i]*c[m][l]*delta(j,k)+c[m][j]*c[m][k]*delta(i,l)+c[m][j]*c[m][l]*delta(i,k)+c[m][k]*c[m][l]*delta(i,j)) \
                        +3*d2*d2*(delta(i,j)*delta(k,l)+delta(i,k)*delta(j,l)+delta(i,l)*delta(j,k)) )
                    H[i][j][k][l] = H[i][j][k][l]*scale/24
    return H

def stripcomment0(f):
    # remove everything in an array after !
    f1='!'
    if f1 in f:
        i=f.index(f1)
        return f[:i]
    else:
        return f

def stripcomment(f):
    # remove everything in an array after !
    i=len(f)
    for f1 in f:
        if f1.startswith('!'):
            i=f.index(f1)
            break
    return f[:i]

def setupswm6():
    l = 0.9572 #OH
    theta = 104.520 #HOH
    d = 0.247 #O-M
    a = 0.315 #O-LP
    phi = 101.098 #LP-O-LP

    mass=[15.9994, 1.008, 1.008, 0.0, 0.0, 0.0] # O, H1, H2, OM, LP1, LP2
    q=[0.2880, 0.5307, 0.5307, -1.1334, -0.1080, -0.1080]
    #coor=[[0.0,0.0,0.0],
    #      [0.0,l*math.sin(math.radians(theta/2)),l*math.cos(math.radians(theta/2))],
    #      [0.0,-1*l*math.sin(math.radians(theta/2)),l*math.cos(math.radians(theta/2))],
    #      [0.0,0.0,d],
    #      [a*math.sin(math.radians(phi/2)),0.0,-1*a*math.cos(math.radians(phi/2))],
    #      [-1*a*math.sin(math.radians(phi/2)),0.0,-1*a*math.cos(math.radians(phi/2))]]
    coor=[[0.0,0.0,0.0],
          [l*math.sin(math.radians(theta/2)),0.0,l*math.cos(math.radians(theta/2))],
          [-1*l*math.sin(math.radians(theta/2)),0.0,l*math.cos(math.radians(theta/2))],
          [0.0,0.0,d],
          [0.0,a*math.sin(math.radians(phi/2)),-1*a*math.cos(math.radians(phi/2))],
          [0.0,-1*a*math.sin(math.radians(phi/2)),-1*a*math.cos(math.radians(phi/2))]]
    return q, coor

def centerofmass(mass,c):
    com=[]
    for i in range(3):
        a = 0.0
        for j in range(6):
            a += mass[j]*c[j][i]
        com.append(a/sum(mass))
    return com

def lppos(dist,angle,dihe):
    return [dist*math.sin(math.radians(angle))*math.cos(math.radians(dihe)), dist*math.sin(math.radians(angle))*math.sin(math.radians(dihe)), dist*math.cos(math.radians(angle))]

def comp_aniso(a):
    # calculate the anisotropy, just consider a few situations
    a11=a.anisotropy[4]
    a22=a.anisotropy[5]
    a33=3.0-a11-a22
    a11 = a.alpha * a11
    a22 = a.alpha * a22
    a33 = a.alpha * a33

    nlp=len(a.lonepairs)
    if nlp == 2:
        if set([a.lonepairs[0][0], a.lonepairs[1][0]]) == set(a.anisotropy[2:4]):
            pos1=lppos(a.lonepairs[0][6],a.lonepairs[0][7],a.lonepairs[0][8])
            pos2=lppos(a.lonepairs[1][6],a.lonepairs[1][7],a.lonepairs[1][8])
            if numeq(pos1[1],pos2[1]) and numeq(pos1[2],pos2[2]) and a.lonepairs[0][1]==a.lonepairs[1][1] and a.lonepairs[0][2]==a.lonepairs[1][2] and a.lonepairs[0][3]==a.lonepairs[1][3] and a.lonepairs[0][5]==a.lonepairs[1][5]:
                if a.lonepairs[0][5].lower().startswith('rela'):
                    str0=' ZX '+a.lonepairs[0][1]+' '+a.lonepairs[0][2]+' '+a.lonepairs[0][3]
                elif a.lonepairs[0][5].lower().startswith('bise'):
                    str0=' BISECT '+a.lonepairs[0][1]+' '+a.lonepairs[0][2]+' '+a.lonepairs[0][3]
                else:
                    print 'unsupport lonepair type', a.anisotropy, a.lonepairs
            elif numeq(pos1[0],pos2[0]) and numeq(pos1[2],pos2[2]) and a.lonepairs[0][1]==a.lonepairs[1][1] and a.lonepairs[0][2]==a.lonepairs[1][2] and a.lonepairs[0][3]==a.lonepairs[1][3] and a.lonepairs[0][5]==a.lonepairs[1][5]:
                if a.lonepairs[0][5].lower().startswith('rela'):
                    str0=' ZX '+a.lonepairs[0][1]+' '+a.lonepairs[0][2]+' '+a.lonepairs[0][3]
                    # in this case it's ZY so have to exchange a11 and a22 for the general printansio to work
                    aexch=a11
                    a11=a22
                    a22=aexch
                else:
                    print 'unsupport lonepair type', a.anisotropy, a.lonepairs
            else:
                print 'unsupport anisotropic polarizability', a.anisotropy, a.lonepairs
        elif a.lonepairs[0][0] not in a.anisotropy[1:4] and a.lonepairs[1][0] not in a.anisotropy[1:4]:  # in the case no lonepair is used in defining anisotropy
            str0=' ZBISECT '+a.anisotropy[1]+' '+a.anisotropy[2]+' '+a.anisotropy[3]
    elif nlp==1:
        pos1=lppos(a.lonepairs[0][6],a.lonepairs[0][7],a.lonepairs[0][8])
        # here often angle and dihe are set to 179.99 but meant to be 180.0, thus numeqloose
        if numeqloose(pos1[0],0.0) and numeqloose(pos1[1],0.0) and a.lonepairs[0][0] == a.anisotropy[1] and a.lonepairs[0][5]:
            str0=' BISECT '+a.lonepairs[0][1]+' '+a.lonepairs[0][2]+' '+a.lonepairs[0][3]
        else:
            print 'unsupport lonepair type', a.anisotropy, a.lonepairs
    elif nlp==0:
        str0=' ZBISECT '+a.anisotropy[1]+' '+a.anisotropy[2]+' '+a.anisotropy[3]
    elif nlp==3:

# TODO: do sanity check here
        #here only support the third lonepair is a dummy one

        # ANISOTROPY O1 LPX LP1A LP1B A11 0.8889  A22 1.2222 LONEPAIR bisector LPX  O1 C1 C4 distance 0.10 angle   0.0 dihe   0.0 LONEPAIR bisector LP1A O1 C1 C4 distance 0.35 angle 110.0 dihe  90.0 LONEPAIR bisector LP1B O1 C1 C4 distance 0.35 angle 110.0 dihe 270.0 0.618
        if numeq(a.lonepairs[2][4], 0.0):
            pos1=lppos(a.lonepairs[0][6],a.lonepairs[0][7],a.lonepairs[0][8])
            pos2=lppos(a.lonepairs[1][6],a.lonepairs[1][7],a.lonepairs[1][8])
            pos3=lppos(a.lonepairs[2][6],a.lonepairs[2][7],a.lonepairs[2][8])
            print a.anisotropy, a.lonepairs
            print pos1, pos2, pos3
            str0='  BISECT    '+a.lonepairs[0][1]+' '+a.lonepairs[0][2]+' '+a.lonepairs[0][3]

    else:
        print 'unsupport anisotropic polarizability:', nlp, 'lonepairs'

    a.printaniso=str0+' '+str(a33)+' '+str(a22)+' '+str(a11)



def lp_to_mpole(a):
    # generate up to octupole
    # unit conversion factors
    debye2eA=0.20819434
    bohr2A=0.52917721092

    nlp=len(a.lonepairs)
    if nlp == 2:
        q, coor, str0 = multipole2(a)
    elif nlp == 1:
        q, coor, str0 = multipole1(a)
    elif nlp == 3:
        if numeq(a.lonepairs[2][4], 0.0):
            # the third lonepair only for anisotropy
            q, coor, str0 = multipole2(a)
        else:
	    q, coor, str0 = multipole3(a)
            print a.lonepairs

    #convert A to Bohr so that the results are all in a.u. 
    for ci in coor:
        for j in range(3):
            ci[j] = ci[j]/bohr2A

    o1=dipole(q,coor,1.0)
    o2=quadrupole(q,coor,1.0)
    o3=octupole(q,coor,1.0)

    #update the charge of parent atom
    a.charge = sum(q)

    #print out dipole to octupole 
    str1='%7.5f %7.5f %7.5f -' % (o1[0],o1[1],o1[2])
    str2='%7.5f %7.5f %7.5f %7.5f %7.5f %7.5f -' % (o2[0][0], o2[0][1],o2[1][1],o2[0][2],o2[1][2],o2[2][2]) #xx, xy, yy, xz, yz, zz
    str3='%7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f ' % (o3[0][0][0],o3[0][0][1],o3[0][1][1],o3[1][1][1],o3[0][0][2],o3[0][1][2],o3[1][1][2],o3[0][2][2],o3[1][2][2],o3[2][2][2]) #xxx, xxy, xyy, yyy, xxz, xyz, yyz, xzz, yzz, zzz

    a.printmpole=str0+'\n'+str1+'\n'+str2+'\n'+str3

def multipole1(a):
    if a.lonepairs[0][5].lower().startswith('bise'):
        q=[a.charge, a.lonepairs[0][4]]
        coor=[[0.0,0.0,0.0]]
        coor.append(lppos(a.lonepairs[0][6],a.lonepairs[0][7],a.lonepairs[0][8]))
        octstr='OPOLE '+a.name+'  BISECT  '+a.lonepairs[0][1]+' '+a.lonepairs[0][2]+' '+a.lonepairs[0][3]+' -'
    else: 
        print "lonepair types/combinations not supported", a.lonepairs

    return q, coor, octstr

def multipole2(a):
    # subroutine with two relative or bisector lonepairs
    #sanity check
    if (a.lonepairs[0][5].lower().startswith('rela') and a.lonepairs[1][5].lower().startswith('rela') and a.lonepairs[0][1]==a.lonepairs[1][1] and a.lonepairs[0][2]==a.lonepairs[1][2] and a.lonepairs[0][3]==a.lonepairs[1][3]):
        octstr='OPOLE '+a.name+'  ZX  '+a.lonepairs[0][1]+' '+a.lonepairs[0][2]+' '+a.lonepairs[0][3]+' -'
    elif (a.lonepairs[0][5].lower().startswith('bise') and a.lonepairs[1][5].lower().startswith('bise') and a.lonepairs[0][1]==a.lonepairs[1][1] and a.lonepairs[0][2]==a.lonepairs[1][2] and a.lonepairs[0][3]==a.lonepairs[1][3]):
        octstr='OPOLE '+a.name+'  BISECT  '+a.lonepairs[0][1]+' '+a.lonepairs[0][2]+' '+a.lonepairs[0][3]+' -'
    elif (a.lonepairs[0][5].lower().startswith('bise') and a.lonepairs[1][5].lower().startswith('bise') and a.lonepairs[0][1]==a.lonepairs[1][1] and a.lonepairs[0][2]==a.lonepairs[1][3] and a.lonepairs[0][3]==a.lonepairs[1][2]): # the case of reverting two bisector atoms
    # reorder their position and map phi into 180 - phi
        a.lonepairs[1][2]=a.lonepairs[0][2]
        a.lonepairs[1][3]=a.lonepairs[0][3]
        a.lonepairs[1][8] = 180.0 - a.lonepairs[1][8]
        octstr='OPOLE '+a.name+'  BISECT  '+a.lonepairs[0][1]+' '+a.lonepairs[0][2]+' '+a.lonepairs[0][3]+' -'
    else:
        print "lonepair types/combinations not supported!", a.lonepairs #sys.exit("lonepair types/combinations not supported!")
    q=[a.charge, a.lonepairs[0][4], a.lonepairs[1][4]]
    coor=[[0.0,0.0,0.0]]
    coor.append(lppos(a.lonepairs[0][6],a.lonepairs[0][7],a.lonepairs[0][8]))
    coor.append(lppos(a.lonepairs[1][6],a.lonepairs[1][7],a.lonepairs[1][8]))

    return q, coor, octstr

def multipole3(a):
    # simple sanity check
    if (a.lonepairs[0][5].lower().startswith('bise') and a.lonepairs[1][5].lower().startswith('bise') and a.lonepairs[2][5].lower().startswith('bise')  and a.lonepairs[0][1]==a.lonepairs[1][1] and a.lonepairs[0][2]==a.lonepairs[1][2] and a.lonepairs[0][3]==a.lonepairs[1][3] and a.lonepairs[0][1]==a.lonepairs[2][1] and a.lonepairs[0][2]==a.lonepairs[2][2] and a.lonepairs[0][3]==a.lonepairs[2][3]):
        octstr='OPOLE '+a.name+'  BISECT  '+a.lonepairs[0][1]+' '+a.lonepairs[0][2]+' '+a.lonepairs[0][3]+' -'
    q=[a.charge, a.lonepairs[0][4], a.lonepairs[1][4], a.lonepairs[2][4]]
    coor=[[0.0,0.0,0.0,0.0]]
    coor.append(lppos(a.lonepairs[0][6],a.lonepairs[0][7],a.lonepairs[0][8]))
    coor.append(lppos(a.lonepairs[1][6],a.lonepairs[1][7],a.lonepairs[1][8]))
    coor.append(lppos(a.lonepairs[2][6],a.lonepairs[2][7],a.lonepairs[2][8]))

    return q, coor, octstr

def getFieldPairs(fields):
    pairs = []
    for i in range(len(fields)/2):
        pairs.append((fields[2*i], fields[2*i+1]))
    return pairs

class Residue(object):
    def __init__(self, name, charge):
        self.name = name
        self.charge = charge
        self.atoms = []
        #self.atomMap = {}
        self.deletions = []
        self.bonds = []
        self.externalBonds = ['N', 'C']
        if name == 'GLY':
            self.patches = ['GNTE', 'CTEG']
        elif name == 'PRO':
            self.patches = ['PROP', 'CTEP']
        else:
            self.patches = ['NTER', 'CTER']
        self.lonepairs = []
        self.group = []
        self.ic = []
        self.patch = []
        self.impr = []
        self.cmap = []
        self.donor = []

    def addAtom(self, atom):
        self.atoms.append(atom)

    def deleteAtom(self, atom):
        # update the group definition
        lpindex = self.atoms.index(atom)
        ngroup=len(self.group)
        for i in range(ngroup):
            if self.group[i] > lpindex:
                self.group[i] = self.group[i] - 1
        self.atoms.remove(atom)
        self.lonepairs.append(atom)

    def hasLonepair(self, aname):
        for a in self.lonepairs:
            if aname == a.name:
                return True
        return False
        
class Atom(object):
    def __init__(self, fields):
        self.name = fields[1]
        self.atomClass = fields[2]
        self.charge = float(fields[3])
        self.polarizable = False
        self.anisotropic = False
        self.haslonepair = False
        self.lonepairs = []
        for param, value in getFieldPairs(fields[4:]):
            if param == 'ALPHA':
                self.polarizable = True
                self.alpha = -1*float(value) # need to be float in case of anisotropic
            elif param == 'THOLE':
                self.thole = value
            elif param == 'TYPE':
                self.drudeType = value
        if self.polarizable:
            if 'thole' not in dir(self):
                self.thole = 1.3
        #self.type = len(atomTypes)
        #atomTypes.append(self)

    def addLonepair(self, fields):
        self.haslonepair = True
        self.lonepairs.append(fields)
        self.printmpole = None

    def addAnisotropy(self, fields):
        self.anisotropic = True
        self.anisotropy = fields
        self.printaniso = None


if __name__ == "__main__":
    residue = None
    residues = []
    section = 0


    inputfile=sys.argv[1]
    print "* convert Drude topology file "+inputfile+" to MPID."
    for line in open(inputfile):
        fields = line.split()
        if line[0:4]== 'RESI':
            residue = Residue(fields[1], fields[2])
            residues.append(residue)
            section = 1
        elif line[0:4]== 'PRES': # mpole currently doesn't support patch
            section = 2
        elif line[0:35]=='AUTOGENERATE ANGLES DIHEDRALS DRUDE':
            print 'AUTOGENERATE ANGLES DIHEDRALS MPOLE'
        else:
            if section == 0:
                print line,
            elif section == 1: # process one residue
                if len(fields) == 0:
                    continue
                if fields[0] == 'ATOM':
                    afield = stripcomment(fields)
                    residue.addAtom(Atom(afield))
                elif fields[0] == 'GROUP':
                    residue.group.append(len(residue.atoms))
                elif fields[0] == 'BOND':
                    afield = stripcomment(fields[1:])
                    for name1, name2 in getFieldPairs(afield):
                        residue.bonds.append((name1, name2))
                elif fields[0] == 'ANISOTROPY':
                    params = dict(getFieldPairs(fields[5:]))
                    for a in residue.atoms:
                        if a.name == fields[1]:
                           lpparent=a
                    anisoinfo = fields[1:5]+[float(params['A11']), float(params['A22'])]
                    lpparent.addAnisotropy(anisoinfo)
                elif fields[0] == 'LONEPAIR':
                    params = dict(getFieldPairs(fields[6:]))
                    for a in residue.atoms:
                        if a.name == fields[2]:
                           lp=a
                        elif a.name == fields[3]:
                           lpparent=a
                    lpinfo = fields[2:6]+[lp.charge, fields[1], float(params['distance']), float(params['angle']), float(params['dihe'])]
                    residue.deleteAtom(lp)
                    lpparent.addLonepair(lpinfo)
                elif fields[0] == 'IC':
                    residue.ic.append(line)
                elif fields[0].lower() == 'patch':
                    residue.patch.append(line)
                elif fields[0] == 'IMPR':
                    residue.impr.append(line)
                elif fields[0] == 'CMAP':
                    residue.cmap.append(line)
                elif fields[0].startswith('DONO') or fields[0].startswith('ACCE'):
                    residue.donor.append(line)


            else:
                continue

    for r in residues:
        # remove bonds that have lonepairs, note that there are situations as +N, -C
        for b in list(r.bonds): # make bonds a new list to loop over
            if (r.hasLonepair(b[0]) or r.hasLonepair(b[1])):
                r.bonds.remove(b)
        for a in r.atoms:   
            if a.haslonepair:
            # map atoms with lonepairs into octpole
                lp_to_mpole(a)
            if a.anisotropic:
            # compute the anisotropy
                comp_aniso(a)

    for r in residues:
        print 'RESI ', r.name, r.charge

        #print out the atom with group information
        aindex=0
        for a in r.atoms:
            if aindex in r.group:
                print 'GROUP'
            print 'ATOM ', a.name, a.atomClass, a.charge
            aindex += 1
        print

        #print out polarizability
        for a in r.atoms:
            if a.polarizable:
                if a.anisotropic:
                    print 'APOLARIZE', a.name, a.printaniso, a.thole
                else:
                    print 'POLARIZE', a.name, a.alpha, a.thole
        print

        #print out the bond
        for b in r.bonds:
            print 'BOND ', b[0], b[1]
        print

        #print out the IMPR and CMAP part
        for info in r.impr:
            print info,
        for info in r.cmap:
            print info,
        print

        #print out multipole
        for a in r.atoms:
            if a.haslonepair:
            # map atoms with lonepairs into octpole
                print a.printmpole
        print

        # print out IC table and patch information
        for info in r.donor:
            print info,
        #for info in r.ic: # ignore IC table as some IC such as the one in SWM4 water, has lonepairs in it
        #    print info,
        for info in r.patch:
            print info,
        print

    print 'end'


