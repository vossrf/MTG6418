# Discrete Math routines to simplify diagram display in Sagemath/Python
# Prof. Richard F. Voss rvoss@fau.edu
import matplotlib.pyplot as plt  
from matplotlib.patches import Polygon
from math import pi,sin,cos,sqrt

#fig = None   # globals for matplotlib
#ax  = None
# fig, ax = plt.subplots(1,1)

fig = None; ax = None   # defaults

def setFigAx(f=None,a=None): # set globals
    global ax,fig
    if f==None or a==None: fig,ax = plt.subplots(1,1)
    else:
        fig = f
        ax = a    

Nopts = {'fill':True,       # default options for circular nodes
         'facecolor':'lightyellow',
         'edgecolor':'black',
         'linewidth':2,
         'zorder':1}    
Topts = {'fontsize':12,       # default options for text 
         'color':'black',
         'backgroundcolor': 'none',
         'fontweight':'normal',   # could be bold
         'ha':'center',
         'va':'center',
         'zorder':3}
Aopts = {'color':'darkgreen',   # default options for arrows
         'width':0.01,   
         'head_width':0.1,
         'length_includes_head':True,
         'overhang':0.3,
         'zorder':2}
A2opts = {'color':'darkblue',   # default options for 2-headed arrows
         'width':0.03,   
         'head_width':0.15,
         'length_includes_head':True,
         'overhang':0.3,
         'zorder':2}
Bopts = {'fill':True,       # default options for boxes
         'facecolor':(.85,.9,1),
         'edgecolor':'black',
         'linewidth':1,
         'zorder':0}  
         
def updateDict(old,update):
    "update values of existing keys in old dictionary to new values, return copy of updated"
    new = {}           # new dictionary combining old and updates
    for key in old:    # only update existing keys
        if key in update:
            new[key] = update[key]
        else:
            new[key] = old[key]
    return new        
         
def Rarrow(xy1,xy2,r=0.3,**kwargs):
    "draw arrow from xy1 to xy2, endpoints reduced by distance r"
    opts = updateDict(Aopts,kwargs)  # non-default options
    d = ((xy1[0]-xy2[0])**2+(xy1[1]-xy2[1])**2)**0.5  # distance between endpoints
    if d<2*r: return
    f = r/d
    x1 = xy1[0]+f*(xy2[0]-xy1[0])
    y1 = xy1[1]+f*(xy2[1]-xy1[1])
    x2 = xy2[0]+f*(xy1[0]-xy2[0])
    y2 = xy2[1]+f*(xy1[1]-xy2[1])
    ax.arrow(x1,y1,x2-x1,y2-y1,**opts)    
    
def R2arrow(xy1,xy2,r=0.3,**kwargs):
    "draw double headed arrow from xy1 to xy2, endpoints reduced by distance r"
    opts = updateDict(A2opts,kwargs)  # non-default options
    d = ((xy1[0]-xy2[0])**2+(xy1[1]-xy2[1])**2)**0.5  # distance between endpoints
    if d<2*r: return
    f = r/d
    x1 = xy1[0]+f*(xy2[0]-xy1[0])
    y1 = xy1[1]+f*(xy2[1]-xy1[1])
    x2 = xy2[0]+f*(xy1[0]-xy2[0])
    y2 = xy2[1]+f*(xy1[1]-xy2[1])
    ax.arrow(x1,y1,x2-x1,y2-y1,**opts)     
    ax.arrow(x2,y2,x1-x2,y1-y2,**opts)   
    opts.update(head_length=0,color='lightyellow',width=.001,head_width=0)
    ax.arrow(x2,y2,x1-x2,y1-y2,**opts)
    
def Node(xy,r,**kwargs):
    "add circular node to graphics, return position tuple"
    opts = updateDict(Nopts,kwargs)  # non-default options    
    c = plt.Circle(xy,r,**opts)
    ax.add_artist(c)

def Text(xy,t,**kwargs):
    opts = updateDict(Topts,kwargs)  # non-default options 
    ax.text(xy[0],xy[1],t,**opts)

def Box(xy1,xy2,**kwargs):
    opts = updateDict(Bopts,kwargs)  # non-default options
    x1 = min(xy1[0],xy2[0])
    x2 = max(xy1[0],xy2[0])
    y1 = min(xy1[1],xy2[1])
    y2 = max(xy1[1],xy2[1])    
    poly = Polygon(((x1,y1),(x2,y1),(x2,y2),(x1,y2)),**opts)
    ax.add_patch(poly)

def DrawRelation(lab,S1,r1,S2,r2=[],S3=[],f=None,a=None):
    "create diagram of relation r1 (and optionally r2) as lists of arrow ends, lab list of left, arrow, right"
    global ax,fig
    setFigAx(f,a)
    n1 = len(S1); n2 = len(S2); n3 = len(S3)
    nm = max(n1,n2,n3)
    Box((0.4,nm-n1+0.4),(1.6,nm+0.6))   # background boxes
    Box((3.4,nm-n2+0.4),(4.6,nm+0.6))
    if n3>0: Box((6.4,nm-n3+0.4),(7.6,nm+0.6)) # optional 3rd
    r = 0.30   # size of nodes    
    y = nm                  # y coordinate
    XY1 = {}                # dictionary of node positions
    for e in S1:            # generate nodes, save positions
        XY = (1,y)          # coordinates
        y -= 1
        XY1[e] = XY
        Node(XY,r)
        Text(XY,str(e))     # graphics
    y = nm                  # y coordinate        
    XY2 = {}                # dictionary of node positions
    for e in S2:            # generate nodes, save positions
        XY = (4,y)          # coordinates
        y -= 1
        XY2[e] = XY
        Node(XY,r)
        Text(XY,str(e))     # graphics    
    y = nm                  # y coordinate        
    XY3 = {}                # dictionary of node positions
    for e in S3:            # generate nodes, save positions
        XY = (7,y)          # coordinates
        y -= 1
        XY3[e] = XY
        Node(XY,r)
        Text(XY,str(e))     # graphics         
    for e1,e2 in r1:        # draw the r1 arrows   
        Rarrow(XY1[e1],XY2[e2],r)
    for e2,e3 in r2:        # draw the r2 arrows   
        Rarrow(XY2[e2],XY3[e3],r)        
    nm += 1    
    
    ax.set_ylim(0.3,nm+0.4)
    Rarrow((1,nm),(4,nm),r,color='black')            
    Text((1,nm),lab[0])
    Text((4,nm),lab[2])
    Text((2.5,nm-0.1),lab[1],va='top')
    if n3>0:              # label 2nd relation
        Rarrow((4,nm),(7,nm),r,color='black')            
        Text((7,nm),lab[4])
        Text((5.5,nm-0.1),lab[3],va='top') 
        ax.set_xlim(0.3,7.7)    
    else: ax.set_xlim(0.3,4.7)
    ax.set_axis_off()
    ax.set_aspect('equal')   
    plt.tight_layout()    
    return ax,fig
    
def DiGraph(connects,prt=0,dir=1,rotate=0.0,nodemult=1.0,order=[],nodeXY=[],f=None,a=None,**kwargs):
    "create DiGraph based on connects between nodes, use circular arrangement of nodes"
    global ax,fig
    setFigAx(f,a)
    
    nodes = {}    # empty dictionary with center positions
    edges = {}    # edge dictionary with dir 0,1,2
    for e in sorted(connects):    # dictionary of nodes, add positions later
        n1,n2 = e
        if n1 not in nodes: nodes[n1] = 0
        if n2!=None and n2 not in nodes: nodes[n2] = 0
        if e in edges:
            print(" DiGraph error, duplicate edge", e)
            continue
        elif (n2,n1) in edges:  # reverse already present
            edges[(n2,n1)] = 2  # make it a 2-way connection
        elif n1==n2:
            edges[e] = 0        # self-connection
        else:
            edges[e] = 1        # 1-way connection
            
    if len(order)>0:            # specified list for order of nodes
        for n in nodes:
            if n not in order:
                print("node",n,"not in order, no DiGraph possible")
                return 
        if len(nodeXY)==len(order):     # exact placement of nodes
            xbar = ybar = 0
            for k in range(len(order)): # find center of notes 
                xbar += nodeXY[k][0]
                ybar += nodeXY[k][1]
            xbar /= len(nodeXY)    
            ybar /= len(nodeXY)    
            for k in range(len(order)): # save exact placement of notes 
                xy = (nodeXY[k][0]-xbar,nodeXY[k][1]-ybar)
                nodes[order[k]] = xy
    else:
        order = sorted(nodes.keys())  # default order on circle
    if prt: print('order:',order)    
    Nn = len(nodes)
    
    if len(nodeXY)<1:        # default arrange nodes in circle based on order
        dTheta = 2*pi/Nn     # angle between nodes
        if dir<0: dTheta *= -1    # reverse direction around circle
        R = 1.5/sqrt(sin(dTheta)**2+(1-cos(dTheta))**2)  # radius for min separation
        theta = rotate       # starting angular position
        r = 0.2*nodemult     # node size
        for n in order:      # position each node on circle
            if prt: print('dir %2f theta %0.2f dTheta %0.2f' % (dir,theta,dTheta))
            xy = R*sin(theta),R*cos(theta)
            nodes[n] = xy
            theta += dTheta
    if prt:                # print nodes and edges         
        print(Nn,'nodes')                
        for n in order:
            print("%s: at %5.2f %5.2f" % (n, nodes[n][0],nodes[n][1]))
        print(len(edges),'edges')    
        for e in sorted(edges):
            print(e,edges[e])       
    r = 0.2*nodemult         # node size
    xmin = xmax = nodes[order[0]][0]
    ymin = ymax = nodes[order[0]][1]
    for n in order:          # draw nodes and find bounding box
            xy = nodes[n]
            xmin = min(xmin,xy[0])
            xmax = max(xmax,xy[0])
            ymin = min(ymin,xy[1])
            ymax = max(ymax,xy[1])
            Node(xy,r,**kwargs)       # draw node circle
            Text(xy,str(n),**kwargs)  # text label            
    dr = 2.2*r    
    ax.set_xlim(xmin-dr,xmax+dr)
    ax.set_ylim(ymin-dr,ymax+dr)
    
# draw the directed edges
    for e in edges: 
        n1,n2 = e
        if n1==None or n2==None: continue   # no real connection
        if n1==n2:     # self connection as circle
            x,y = nodes[n1]      # original position of node
            R = sqrt(x*x+y*y)
            dR = (R+0.8*r)/R
            xy = x*dR,y*dR
            ax.add_artist(plt.Circle(xy,1.2*r,fill=False,ec='darkgreen',linewidth=2,zorder=0,**kwargs))
            continue
        dir = edges[e]
        if dir==2:     # double headed
            R2arrow(nodes[n1],nodes[n2],r,**kwargs)
        elif dir==1:   # forward
            Rarrow(nodes[n1],nodes[n2],r,**kwargs)
        else:          # reverse
            Rarrow(nodes[n2],nodes[n1],r,**kwargs)
    ax.set_axis_off()
    ax.set_aspect('equal')        
    plt.tight_layout()            
    return 
  
def Hasse(r,prt=0,**kwargs):
    "create Hasse diagram from list of connections r"
    out = {}    # empty dictionary of connections leaving each node
    rk = []                       # list of connections we will keep
    for n1,n2 in r:               # assign pair connections to out
        if n2==None or n1==n2: continue  # drop loops
        rk.append((n1,n2))    
        if n1 in out:          
            out[n1].add(n2)       # add another outgoing node
        else:                     # initial outgoing node
            out[n1]=set((n2,))  
        if n2 not in out:
            out[n2]=set()         # make sure all nodes included
    for n1 in out:                # all possible starting nodes
        for n2 in out[n1]:        # all possible 1-step destinations
            for n3 in out[n2]:    # all 3-step destinations after n1,n2
                if n3 in out[n1]: # remove direct connectio n1,n3
                    e13 = n1,n3
                    if e13 in rk: 
                        rk.remove(e13) # remove   
    if prt: print('rk =',rk)
    nodes = {}         # empty dictionary of node in,out connections
    for n1,n2 in rk:   # compile note inputs and outputs
        if n1 not in nodes:
            nodes[n1]=[[],[n2],-1]  # initial output node and rank
        else:
            nodes[n1][1].append(n2) # another output
        if n2 not in nodes: 
            nodes[n2]=[[n1],[],-1]  # initial input node and rank
        else:
            nodes[n2][0].append(n1) # another input            
    bot = []         # initial list of nodes with no inputs
    for n in nodes:
        if len(nodes[n][0])==0:
            bot.append(n)
            nodes[n][2]=0           # bot nodes have rank 0
    if prt: print('bot =',bot) 
    rerank = 1
    while rerank>0:      # work up network step by step
        rerank = 0
        for n in nodes:                 # renumber output connections
            rank = nodes[n][2]+1        # current rank+1
            if rank==0: continue        # unranked
            for n2 in nodes[n][1]:      # rerank outputs    
                if rank>nodes[n2][2]:
                    nodes[n2][2]=rank   # 1 step beyond previous
                    rerank += 1
    if prt: print('nodes =',nodes)
    ranks = {}        # dictionary of nodes at each rank
    for n in nodes:
        rank = nodes[n][2]              # rank for this node
        if rank in ranks:
            ranks[rank].append(n)
        else:
            ranks[rank] = [n]
    if prt: print('ranks =',ranks)
    order = []         # empty list for nodes
    XY = []            # empty list for ndde positions
    for rank in sorted(ranks):
        nR = sorted(ranks[rank]) # list of nodes with this rank
        x = (1-len(nR))/2        # left most position
        for n in nR:             # all nodes at this rank
            order.append(n)      # add node to list
            XY.append((x,rank))  # position
            x += 1
    if prt: 
        for k in range(len(order)):
            print("%2d: %s %4.1f %3.1f"%(k,order[k],XY[k][0],XY[k][1]))
    DiGraph(rk,order=order,nodeXY=XY,**kwargs)
    
def isRefl(r):  
    "decide if the pair connections on set A are reflexive"
    nodes = set()      # empty set to collect all nodes
    for n1,n2 in r:    # unpack each pair into 2 nodes
        nodes.add(n1)  # add to set of all nodes
        if n2!=None: nodes.add(n2)
    for n in nodes:         # check every node
        if (n,n) not in r:  # not connected to itself?
            return False    # yes, we're done
    return True  

def isASym(r):  
    "decide if the pair connections r are antisymmetric"
    for n1,n2 in r:            # test all connections 
        if n1==n2 or n2==None: continue    # skip loops and solitary nodes
        if (n2,n1) in r:       # reverse found in r
            return False
    return True    # all connections checked, no 2-way 
 
def isSym(r):  
    "decide if the pair connections r are symmetric"
    for n1,n2 in r:            # test all connections 
        if n1==n2 or n2==None: continue    # skip loops and solitary nodes
        if (n2,n1) not in r:  # no reverse found in r
            return False
    return True    # all connections checked, all 2-way 
    
def isTran(r):  
    "decide if every 2-step connection also has a 1-step shortcut"
    out = {}    # empty dictionary of connections leaving each node
    for n1,n2 in r:               # assign pair connections to out
        if n2==None: continue  
        if n1 in out:          
            out[n1].add(n2)       # add new outgoing node
        else:                     # initial outgoing node
            out[n1] = set((n2,))
        if n2 not in out:         # make sure n2 in out dictionary
            out[n2] = set()
    for n1 in out:                # all possible starting nodes
        for n2 in out[n1]:        # all possible 1-step destinations
            for n3 in out[n2]:    # all 3-step destinations
                if n3 not in out[n1]: # no direct connection from n1
                    return False
    return True                   # no 2-step without 1-step    
    
def isEquiv(r):
    " decide if relation r is an Equivalence Relation"
    if isRefl(r) and isSym(r) and isTran(r): return True
    return False   
        
def isPO(r):
    " decide if relation r is a Partial Ordering"
    if isRefl(r) and isASym(r) and isTran(r): return True
    return False
def isEquiv(r):
    " decide if relation r is an Equivalence Relation"
    if isRefl(r) and isSym(r) and isTran(r): return True
    return False