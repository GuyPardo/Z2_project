# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 11:49:08 2025

@author: owner
"""
import importlib
import qutip as qp
import matplotlib.pyplot as plt

import numpy as np
import lattice_utils
importlib.reload(lattice_utils)
from lattice_utils import *
from qutip_utils import *

J = 0
m = 0
lam_B = 1
lam_E = 0
lat =square_lattice(2, [3,2], pbc = True)
# l0 = lat.links[0]
# l1 = lat.links[1]
# lat.plot()

Nx = lat.size[0]
Ny = lat.size[1]
rangex = range(Nx)
rangey = range(Ny)
descriptor = f"${lat.size[0]}\\times{lat.size[1]}$ sites, $J={J}$, $m={m}$, $\\lambda_B={lam_B}$, $\\lambda_E={lam_E}$"

# add l and u boundary links:
# for i in range(1,Nx-1):
#     lat.sites.append(site([i,-1]))
#     lat.links.append(link(lat.get_site([i,-1]),lat.get_site([i,0])))
#     lat.sites.append(site([i,Ny]))
#     lat.links.append(link(lat.get_site([i,Ny-1]),lat.get_site([i,Ny])))
n = lat.num_links() # number of qubits
lat.plot(); plt.show()
#%%
####### exact Hamiltonian ######
Id = qp.identity([2 for _ in range(n)])
sig = [qp.sigmax(), qp.sigmay(), qp.sigmaz()] 

class OpVisualizer:
    def __init__(self):
        self.enabled = True
        self.reset()

    def reset(self):
        self.ops = []  # list of dicts

    def add_1q(self, link, label):
        if self.enabled:
            self.ops.append({
                "type": "1q",
                "links": (link,),
                "label": label
            })

    def add_2q(self, link1, link2, label):
        if self.enabled:
            self.ops.append({
                "type": "2q",
                "links": (link1, link2),
                "label": label
            })


visualizer = OpVisualizer()

from matplotlib.patches import FancyArrowPatch

def draw_arrow(ax, start, end, color="black", lw=2):
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        linewidth=lw,
        color=color,
        mutation_scale=15,
        shrinkA=0,
        shrinkB=0,
    )
    ax.add_patch(arrow)

def link_midpoint(link):
    return 0.5 * (link.sites[0].coordinates + link.sites[1].coordinates)

def plot_lattice_with_ops(lat, ops, title=None):
    fig, ax = plt.subplots(figsize=(4, 4))
    lat.plot()
    ax.set_aspect("equal")

    if title is not None:
        ax.set_title(title)

    for entry in ops:
        kind = entry["type"]
        links = entry["links"]
        label = entry["label"]

        # ---------- 1-qubit operator ----------
        if kind == "1q":
            lgt = links[0]
            if lat.link_exists(lgt[0], lgt[1]):
                link = lat.get_link(lgt[0], lgt[1])
                pos = link_midpoint(link)
                ax.text(
                    pos[0],
                    pos[1],
                    label,
                    ha="center",
                    va="center",
                    fontsize=10,
                    bbox=dict(facecolor="white", edgecolor="black", alpha=0.8),
                )

        # ---------- 2-qubit operator (arrow) ----------
        elif kind == "2q":
            lgt1, lgt2 = links
            if lat.link_exists(lgt1[0], lgt1[1]) and lat.link_exists(lgt2[0], lgt2[1]):
                link1 = lat.get_link(lgt1[0], lgt1[1])
                link2 = lat.get_link(lgt2[0], lgt2[1])

                start = link_midpoint(link1)
                end = link_midpoint(link2)

                draw_arrow(ax, start, end)

                mid = 0.5 * (start + end)
                ax.text(
                    mid[0],
                    mid[1],
                    label,
                    fontsize=9,
                    ha="center",
                    va="center",
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.6),
                )

    plt.show()

def op(qobj, linkargs, *, _label=None):
    if lat.link_exists(*linkargs):
        idx = lat.link_index(*linkargs)
        operator = tensor_op(qobj, qp.identity(2), idx, n)

        if _label is not None:
            visualizer.add_1q(
                LGT_notation(linkargs[0], linkargs[1]),
                _label
            )
    else:
        operator = Id

    return operator



# def op(qobj, linkargs):
#     if lat.link_exists(*linkargs):
#         idx = lat.link_index(*linkargs)
#         operator = tensor_op(qobj, qp.identity(2), idx, n)
#     # if the given link coordinates are not in the lattice, return Id
#     else:
#        operator =  Id
#     return operator



def X(*linkargs):
    return op(qp.sigmax(), linkargs)
def Y(*linkargs):
    return op(qp.sigmay(), linkargs)
def Z(*linkargs):
    return op(qp.sigmaz(), linkargs)



def HB(lam_B):
    H = qp.qzero([2 for _ in range(n)])
    # plaquette is labeld by nbottom left site
    for i in range(Nx):
        for j in range(Ny):
            # plaq_exists  = lat.link_exists([i,j],0)*lat.link_exists([i,j],1)*lat.link_exists([i+1,j],1)*lat.link_exists([i,j+1],0)
            if lat.plaquette_exists([i,j]):
                H = H + Y([i,j],0) * Y([i+1,j],1) * X([i,j],1) * X([i,j+1],0) * Z([i-1,j+1],0) * Z([i,j+1],1)
    return -lam_B*H
        

def Hm(m):
    H = qp.qzero([2 for _ in range(n)])
    for i in rangex:
        for j in rangey:
            H = H + Z([i,j],0) * Z([i,j],1) * Z([i-1,j],0) * Z([i,j-1],1)
    return -m/2*H

def HE(lam_E):
    # TODO if we add ancilliary links then this is wrong
    H = qp.qzero([2 for _ in range(n)])
    for l in lat.links:
        H = H + tensor_op(qp.sigmaz(), qp.identity(2), lat.links.index(l), n)
    return lam_E*H

def Hint(J):
    H = qp.qzero([2 for _ in range(n)])
    
    # bulk
    for i in range(-1,Nx+1):
        for j in range(-1,Ny+1):
            if lat.link_exists([i,j], 0):
                # horizontal 2-body
                # H = H + Y([i,j],0) * Z([i+1,j-1],1)
                # horizontal 6-body
                H = H + Y([i,j],0) * Z([i,j-1],1) * Z([i-1,j],0) * Z([i,j],1) * Z([i+1,j],1) * Z([i+1,j],0)
            if lat.link_exists([i, j], 1):
                # vertical two-body
                H = H + Y([i,j],1) * Z([i,j],0)
                # vertical 6-body
                # H = H + Y([i,j],1) * Z([i,j-1],1) * Z([i-1,j],0) * Z([i-1,j+1],0) * Z([i,j+1],1) * Z([i,j+1],0)

    return J/2*H



####### Trotterization ###########

# def U(a,b,l1,l2):
#     # if one of the two links is not in the lattice, return Id
#     if lat.link_exists(*l1) and lat.link_exists(*l2):
#         return 1/2*(Id +  op(sig[a],l2) + op(sig[b], l1)  - op(sig[a], l2)*op(sig[b], l1))
#     else:
#         return Id


def U(a, b, l1, l2):
    # use opposite convention of Tomer:
    # U(a,b,l1,l2) is represented by an arrow from l1 to l2
    # to fit Tomer's graphs, shift the operators
    if lat.link_exists(*l1) and lat.link_exists(*l2):

        # ---- visualization (single record!) ----
        visualizer.add_2q(
            LGT_notation(l1[0], l1[1]),  # order matters!
            LGT_notation(l2[0], l2[1]),
            f"{['X','Y','Z'][a]}{['X','Y','Z'][b]}"
        )

        return 1/2 * (
            Id
            + op(sig[a], l1)
            + op(sig[b], l2)
            - op(sig[a], l1) * op(sig[b], l2)
        )
    else:
        return Id


# def ZiX(l1,l2):
#     # if one of the two links is not in the lattice, return Id
#     if lat.link_exists(*l1) and lat.link_exists(*l2):
#         return 1/2*(Id +  op(sig[2],l2) + 1j*op(sig[0], l1)  - 1j*op(sig[2], l2)*op(sig[0], l1))
#     else:
#         return Id


def ZiX(l1, l2):
    if lat.link_exists(*l1) and lat.link_exists(*l2):

        visualizer.add_2q(
            LGT_notation(l1[0], l1[1]),
            LGT_notation(l2[0], l2[1]),
            "ZiX"
        )

        return 1/2 * (
            Id
            + op(sig[2], l1)
            + 1j * op(sig[0], l2)
            - 1j * op(sig[2], l1) * op(sig[0], l2)
        )
    else:
        return Id

# def miX(l1,l2):
#     # if one of the two links is not in the lattice, return Id
#     if lat.link_exists(*l1) and lat.link_exists(*l2):
#         return 1/2*(Id - 1j* op(sig[0],l2) + op(sig[2], l1)  + 1j*op(sig[0], l2)*op(sig[2], l1))
#     else:
#         return Id

def miX(l1, l2):
    if lat.link_exists(*l1) and lat.link_exists(*l2):

        visualizer.add_2q(
            LGT_notation(l1[0], l1[1]),
            LGT_notation(l2[0], l2[1]),
            "miX"
        )

        return 1/2 * (
            Id
            - 1j * op(sig[0], l1)
            + op(sig[2], l2)
            + 1j * op(sig[0], l1) * op(sig[2], l2)
        )
    else:
        return Id


# def RX(theta, link):
#     return op((-1j*theta/2*qp.sigmax()).expm(), link)
# def RY(theta, link):
#     return op((-1j*theta/2*qp.sigmay()).expm(), link)
# def RZ(theta, link):
#     return op((-1j*theta/2*qp.sigmaz()).expm(), link)

def RX(theta, link):
    return op(
        (-1j * theta / 2 * qp.sigmax()).expm(),
        link,
        _label="RX"
    )

def RY(theta, link):
    return op(
        (-1j * theta / 2 * qp.sigmay()).expm(),
        link,
        _label="RY"
    )

def RZ(theta, link):
    return op(
        (-1j * theta / 2 * qp.sigmaz()).expm(),
        link,
        _label="RZ"
    )


def rotate90(link):
    N = max(Nx,Ny)
    origin = [(N-1)/2,(N-1)/2]
    site_vector = np.array(link[0]) - np.array( origin)
    new_site_vector = np.array([site_vector[1], - site_vector[0]])
    new_site = origin + new_site_vector
    if link[1] == 0:
        new_site = new_site - np.array([0,1])
        new_direction = 1
    elif link[1] == 1:
        new_direction = 0
    return (new_site, new_direction)

def flipx(link):
    N = max(Nx,Ny)
    origin = [(N-1)/2,(N-1)/2]
    # origin = [(Nx-1)/2,(Ny-1)/2]
    site_vector = np.array(link[0]) - np.array(origin) 
    direction=link[1]
    new_site_vector = np.array([-site_vector[0], site_vector[1]]) + np.array(origin) 
    if link[1]==0:
        new_site_vector = new_site_vector - np.array([1,0])
    return (new_site_vector,direction)

def shifty(link):
    return (np.array(link[0] + np.array([0,1])), link[1])
    
        

def hor2ver(link):
    return flipx(rotate90(link))

#%%
def trotter_step_E_J(dt):
    trotter = [Id for _ in range(30)]
    trotter_all = Id
    shift = 15

    ##1
    for i in rangex:
        for j in rangey:
            if i % 2:  # odd collum
                trotter[0] = trotter[0] * RZ(2 * dt * lam_E, ([i, j], 0))
                trotter[0] = trotter[0] * RZ(2 * dt * lam_E, ([i, j], 1))
                # trotter[shift] = trotter[shift]*RZ(2 * dt * lam_E, hor2ver(([i,j],0)))
                # trotter[shift] = trotter[shift]*RZ(2 * dt * lam_E, hor2ver(([i,j],1)))
            else:
                trotter[0] = trotter[0] * U(2, 2, ([i, j], 0), ([i, j], 1))
                # trotter[shift] = trotter[shift] * U(2,2,hor2ver(([i,j],0)), hor2ver(([i,j],1)))
                # if j==Ny-1:
                #     trotter[0] * trotter[0] * U(2,2, ([i-1,j],0), ([i,j-1],1)) * U(2,1, ([i,j],0), ([i,j-1],1))

    ##2
    for i in rangex:
        for j in rangey:
            if i % 2:  # odd collum
                trotter[1] = trotter[1] * U(2, 1, ([i, j], 0), ([i, j], 1))
                # trotter[1+shift] = trotter[1+shift] * U(2,1,hor2ver(([i,j],0)), hor2ver(([i,j],1)))
            else:  # even collum
                trotter[1] = trotter[1] * RY(dt * J, ([i, j], 1))
                # trotter[1+shift] = trotter[1+shift]*RY( dt * J, hor2ver(([i,j],1)))
    ##3
    for i in rangex:
        for j in rangey:
            if i % 2:  # odd collum
                trotter[2] = trotter[2] * U(2, 2, ([i - 1, j], 0), ([i, j], 1))
                # trotter[2+shift] = trotter[2+shift] * U(2,2,hor2ver(([i-1,j],0)), hor2ver(([i,j],1)))
            else:  # even collum
                trotter[2] = trotter[2] * U(2, 1, ([i - 1, j], 0), ([i, j - 1], 1))
                # trotter[2+shift] = trotter[2+shift] * U(2,1,hor2ver(([i-1,j],0)), hor2ver(([i,j-1],1)))

    ##4
    for i in rangex:
        for j in rangey:
            if i % 2:  # odd collum
                # if j==Ny-1:
                #     trotter[3] = trotter[3] * U(2,1,([i,j],0), ([i,j-1],1))
                pass
            else:  # even collum
                trotter[3] = trotter[3] * U(2, 2, ([i, j - 1], 1), ([i, j], 0))
                # trotter[3+shift] = trotter[3+shift] * U(2,2,hor2ver(([i,j-1],1)), hor2ver(([i,j],0)))

    ##5
    for i in rangex:
        for j in rangey:
            if i % 2:  # odd collum
                pass
            else:  # even collum
                trotter[4] = trotter[4] * RY(dt * J, ([i, j], 0))
                # trotter[4+shift] = trotter[4+shift] * RY(dt * J, hor2ver(([i,j],0)))

    ##6
    trotter[5] = trotter[3].dag()
    # trotter[5+shift] = trotter[3+shift].dag()

    ##7
    trotter[6] = trotter[2].dag()
    # trotter[6+shift] = trotter[2+shift].dag()

    ##8
    # here there isi a bit of uncertainty
    for i in rangex:
        for j in rangey:
            if i % 2:  # odd collum
                # pass
                trotter[7] = trotter[7] * ZiX(([i, j], 0), ([i, j], 1)).dag()
                # trotter[7] = trotter[7]*U(2,2, ([i+1,j],1),([i,j+1],1))
                # trotter[7+shift] = trotter[7+shift]*ZiX(hor2ver(([i,j],0)),hor2ver(([i,j],1))).dag()
            else:  # even collum
                # pass
                trotter[7] = trotter[7] * ZiX(([i, j], 0), ([i, j], 1))
                # trotter[7+shift] = trotter[7+shift]*ZiX(hor2ver(([i,j],0)),hor2ver(([i,j],1)))

    ##9
    for i in rangex:
        for j in rangey:
            if i % 2:  # odd collum
                trotter[8] = trotter[8] * U(2, 1, ([i - 1, j], 0), ([i, j - 1], 1))
                # trotter[8+shift] = trotter[8+shift] * U(2,1, hor2ver(([i-1,j],0)),hor2ver(([i,j-1],1)))

            else:  # even collum
                trotter[8] = trotter[8] * U(2, 2, ([i - 1, j], 0), ([i, j], 1))
                # trotter[8+shift] = trotter[8+shift] * U(2,2, hor2ver(([i-1,j],0)),hor2ver(([i,j],1)))

    ##10
    for i in rangex:
        for j in rangey:
            if i % 2:  # odd collum
                trotter[9] = trotter[9] * U(2, 2, ([i, j - 1], 1), ([i, j], 0))
                # trotter[9+shift] = trotter[9+shift] * U(2,2,hor2ver(([i,j-1],1)), hor2ver(([i,j],0)))
            else:  # even collum
                pass
    ##11
    for i in rangex:
        for j in rangey:
            if i % 2:  # odd collum
                trotter[10] = trotter[10] * RY(dt * J, ([i, j], 0))
                # trotter[10+shift] = trotter[10+shift] * RY(dt * J, hor2ver(([i,j],0)))
            else:  # even collum
                pass

    ##12
    trotter[11] = trotter[9].dag()
    # trotter[11+shift] = trotter[9+shift].dag()

    ##13
    trotter[12] = trotter[8].dag()
    # trotter[12+shift] = trotter[8+shift].dag()

    ##14
    for i in rangex:
        for j in rangey:
            if i % 2:  # odd collum
                trotter[13] = trotter[13] * RY(dt * J, ([i, j], 1))
                # trotter[13+shift] = trotter[13+shift] * RY(dt * J, hor2ver(([i,j],1)))
            else:  # even collum
                trotter[13] = trotter[13] * U(2, 1, ([i, j], 0), ([i, j], 1))
                # trotter[13+shift] = trotter[13+shift] * U(2,1,hor2ver(([i,j],0)), hor2ver(([i,j],1)))

    ##15
    for i in rangex:
        for j in rangey:
            if i % 2:  # odd collum
                trotter[14] = trotter[14] * U(2, 2, ([i, j], 0), ([i, j], 1))
                # trotter[14+shift] = trotter[14+shift] * U(2,2,hor2ver(([i,j],0)), hor2ver(([i,j],1)))
            else:  # even collum
                trotter[14] = trotter[14] * RZ(2 * dt * lam_E, ([i, j], 0))
                trotter[14] = trotter[14] * RZ(2 * dt * lam_E, ([i, j], 1))

    for op in trotter:
        trotter_all =  op*trotter_all

    return trotter_all


def trotter_step_B_m(dt, plot_ops=False):
    trotter = [Id for _ in range(22)]
    trotter_all = Id
    diagonal = lambda i, j: [(j - i) % 4 == k for k in range(4)]
    ##1
    visualizer.reset()
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[2]:
                trotter[0] = trotter[0] * U(2, 1, ([i, j], 0), ([i, j], 1))
                trotter[0] = trotter[0] * U(2, 1, ([i, j - 1], 1), ([i - 1, j], 0))
    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 1")

    ##2
    visualizer.reset()
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[1]:  # 1/4 diagonal
                trotter[1] = trotter[1] * U(2, 1, ([i, j], 1), ([i, j], 0))
                trotter[1] = trotter[1] * U(2, 1, ([i - 1, j], 0), ([i, j - 1], 1))
    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 2")

    ##3
    visualizer.reset()
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[1]:
                trotter[2] = trotter[2] * U(2, 1, ([i, j - 1], 1), ([i, j], 0))
            if diagonal(i, j)[2]:
                trotter[2] = trotter[2] * U(2, 1, ([i - 1, j], 0), ([i, j], 1))
    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 3")
        visualizer.reset()
    ##4

    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[1]:
                trotter[3] = trotter[3] * RZ(-dt * m, ([i,j],0))
            if diagonal(i, j)[2]:
                trotter[3] = trotter[3] * RZ(-dt * m, ([i,j],1))
    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 4")
        visualizer.reset()

    ##5
    # trotter[4] = trotter[2].dag()
    visualizer.reset()
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[1]:
                trotter[4] = trotter[4] * U(2, 1, ([i, j - 1], 1), ([i, j], 0))
            if diagonal(i, j)[2]:
                trotter[4] = trotter[4] * U(2, 1, ([i - 1, j], 0), ([i, j], 1))
    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 5")
        visualizer.reset()
    ##6
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[0]:
                trotter[5] = trotter[5] * U(2,1,([i,j],1),([i,j],0))
                trotter[5] = trotter[5] * U(2, 1, ([i-1, j], 0), ([i, j-1], 1))
                # trotter[5] = trotter[5] * U(1,2,([i,j],1),([i,j],0))
                # trotter[5] = trotter[5] * U(1, 2, ([i-1, j], 0), ([i, j-1], 1))
            if diagonal(i, j)[2]:
                trotter[5] = trotter[5] * miX(([i,j],1), ([i,j],0))
                trotter[5] = trotter[5] * miX(([i-1, j], 0), ([i, j-1], 1))
    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 6")
        visualizer.reset()


    ##7
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[1]:
                trotter[6] = trotter[6] * U(0, 2, ([i,j-1],1),([i,j],0))
            if diagonal(i, j)[2]:
                trotter[6] = trotter[6] * U(0, 2, ([i,j-1],1),([i,j],0))

    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 7")
        visualizer.reset()
    ##8
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[1]:
                trotter[7] = trotter[7] * RX(-2 * dt * lam_B, ([i,j], 0))
                trotter[7] = trotter[7] * RX(-2 * dt * lam_B, ([i-1, j], 0))
                # trotter[7] = trotter[7] * X(*([i,j], 0))
                # trotter[7] = trotter[7] * X(*([i-1, j], 0))
    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 8")
        visualizer.reset()
    ##9
    # trotter[8] = trotter[6].dag()
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[1]:
                trotter[8] = trotter[8] * U(0, 2, ([i,j-1],1),([i,j],0))
            if diagonal(i, j)[2]:
                trotter[8] = trotter[8] * U(0, 2, ([i,j-1],1),([i,j],0))

    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 9")
        visualizer.reset()
    ##10a
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[0]:
                #TODO check order of multiplication
                # trotter[9] = trotter[9] * U(1 ,2, ([i, j], 1), ([i, j], 0)) * U(2, 1, ([i, j], 1), ([i, j], 0))
                # trotter[9] = trotter[9] * U(1, 2 ,([i-1, j], 0), ([i, j-1], 1)) * U(2,1, ([i-1, j], 0), ([i, j-1], 1))
                trotter[9] = trotter[9] * U(2, 1, ([i, j], 1), ([i, j], 0))
                trotter[9] = trotter[9] * U(2, 1, ([i - 1, j], 0), ([i, j - 1], 1))
            if diagonal(i, j)[2]:
                trotter[9] = trotter[9] * U(2, 2, ([i, j], 1), ([i, j], 0))
                trotter[9] = trotter[9] * U(2, 2, ([i-1, j], 0), ([i, j-1], 1))
    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 10a")
        visualizer.reset()

    #10b
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[0]:
                #TODO check order of multiplication
                trotter[10] = trotter[10] * U(1, 2, ([i, j], 1), ([i, j], 0))
                trotter[10] = trotter[10] * U(1, 2, ([i - 1, j], 0), ([i, j - 1], 1))

    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 10b")
        visualizer.reset()
    ##11
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[1] or diagonal(i, j)[3]:
                trotter[11] = trotter[11] * U(2, 1, ([i, j], 1), ([i, j], 0))
                trotter[11] = trotter[11] * U(2, 1, ([i-1, j], 0), ([i, j-1], 1))


    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 11")
        visualizer.reset()

    ##12
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[0]:
                trotter[12] = trotter[12] * U(2, 1, ([i-1, j], 0), ([i, j], 1))
            if diagonal(i, j)[3]:
                trotter[12] = trotter[12] * U(2, 1, ([i, j-1], 1), ([i, j], 0))
    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 12")
        visualizer.reset()

    ##13
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[0]:
                trotter[13] = trotter[13] * RZ(-dt * m, ([i,j],1))
            if diagonal(i, j)[3]:
                trotter[13] = trotter[13] * RZ(-dt * m, ([i, j], 0))
    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 13")
        visualizer.reset()
    ##14
    # trotter[13] = trotter[11].dag()
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[0]:
                trotter[14] = trotter[14] * U(2, 1, ([i-1, j], 0), ([i, j], 1))
            if diagonal(i, j)[3]:
                trotter[14] = trotter[14] * U(2, 1, ([i, j-1], 1), ([i, j], 0))
    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 14")
        visualizer.reset()

    ##15a
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[0]:
                # trotter[14] = trotter[14] * miX(([i, j], 1), ([i, j], 0))
                # trotter[14] = trotter[14] * miX(([i-1, j], 0), ([i, j-1], 1))
                trotter[15] = trotter[15] * U(1,2,([i, j], 1), ([i, j], 0))
                trotter[15] = trotter[15] * U(1,2,([i-1, j], 0), ([i, j-1], 1))
            if diagonal(i, j)[2]:
                trotter[15]= trotter[15] * U(2,1,([i, j], 1), ([i, j], 0))
                trotter[15] = trotter[15] * U(2, 1, ([i-1, j], 0), ([i, j-1], 1))

    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 15a")
        visualizer.reset()


    ##15b
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[0]:
                # trotter[14] = trotter[14] * miX(([i, j], 1), ([i, j], 0))
                # trotter[14] = trotter[14] * miX(([i-1, j], 0), ([i, j-1], 1))
                trotter[16] = trotter[16] * U(2,2,([i, j], 1), ([i, j], 0))
                trotter[16] = trotter[16] * U(2,2,([i-1, j], 0), ([i, j-1], 1))

    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 15b")
        visualizer.reset()
    ##16
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[0]:
                trotter[17] = trotter[17] * U(0,2,([i, j-1], 1), ([i, j], 0))
            if diagonal(i, j)[3]:
                trotter[17] = trotter[17] * U(0,2,([i, j-1], 1), ([i, j], 0))
    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 16")
        visualizer.reset()
    ##17
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[3]:
                if lat.plaquette_exists([i, j-1]):
                    trotter[18] = trotter[16] * RX(-2 * dt * lam_B, ([i, j], 0))
                if lat.plaquette_exists([i-1, j - 1]):
                    trotter[18] = trotter[16] * RX(-2 * dt * lam_B, ([i - 1, j], 0))

    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 17")
        visualizer.reset()
    ##18
    # trotter[17] =  trotter[15].dag()
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[0]:
                trotter[19] = trotter[19] * U(0,2,([i, j-1], 1), ([i, j], 0))
            if diagonal(i, j)[3]:
                trotter[19] = trotter[19] * U(0,2,([i, j-1], 1), ([i, j], 0))
    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 18")
        visualizer.reset()

    ##19
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[0]:
                trotter[20] = trotter[20] * U(2,2,([i, j], 1), ([i, j], 0))
                trotter[20] = trotter[20] * U(2, 2, ([i-1, j], 0), ([i, j-1], 1))
            if diagonal(i, j)[2]:
                trotter[20] = trotter[20] * U(2,1,([i, j], 1), ([i, j], 0))
                trotter[20] = trotter[20] * U(2, 1, ([i-1, j], 0), ([i, j-1], 1))
    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 19")
        visualizer.reset()

    ##20
    for i in rangex:
        for j in rangey:
            if diagonal(i, j)[3]:
                trotter[21] = trotter[21] * U(2,1,([i, j], 1), ([i, j], 0))
                trotter[21] = trotter[21] * U(2, 1, ([i-1, j], 0), ([i, j-1], 1))
    if plot_ops:
        plot_lattice_with_ops(lat, visualizer.ops, title="Step 20")
        visualizer.reset()
    # combine:
    for op in trotter:
        trotter_all = op*trotter_all

    return trotter_all



trotter_step_B_m(1,plot_ops=True)

H = HE(lam_E) + Hint(J) + HB(lam_B) + Hm(m)
#%% test time evolution
from tqdm import tqdm

T = 1
# print(H)
prop = (-1j*T*H).expm()
err_sum = []
n_trotter_range = np.array(range(50,1000, 20))
# n_trotter_range = np.array([200])
for n_trotter in tqdm(n_trotter_range):
    dt = T/n_trotter
    tr = trotter_step_B_m(dt)*trotter_step_E_J(dt)
    # prop_tr = Id
    prop_tr = tr**n_trotter
    # for _ in range(n_trotter):
    #     prop_tr = tr*prop_tr
    err = np.abs(prop.full() - prop_tr.full())
    err_sum.append(np.sum(err)/np.sum(np.abs(prop.full())))
#%%

plt.figure()
plt.imshow(np.real(prop_tr.full()))
plt.colorbar()
plt.title("Trotter")
plt.show()

plt.figure()
plt.imshow(np.real(prop.full()))
plt.title("exact")
plt.colorbar()
plt.show()

plt.figure()
plt.imshow(np.abs(err))
plt.title("err")
plt.colorbar()
plt.show()


plt.figure()
# plt.semilogy(n_trotter_range,err_sum)
plt.plot(n_trotter_range,err_sum)
plt.xlabel("# of trotter steps")
plt.title(descriptor)
plt.ylabel("time evolution error")
plt.show()

#%% 0 helper
#
# from matplotlib.patches import FancyArrowPatch
#
# def draw_arrow(ax, start, end, color="black", lw=2, arrowstyle="-|>"):
#     arrow = FancyArrowPatch(
#         start,
#         end,
#         arrowstyle=arrowstyle,
#         linewidth=lw,
#         color=color,
#         shrinkA=0,
#         shrinkB=0,
#         mutation_scale=15,
#     )
#     ax.add_patch(arrow)
# def link_midpoint(link):
#     return 0.5 * (link.sites[0].coordinates + link.sites[1].coordinates)
#
#
# def annotate_link(ax, link, label, color="black"):
#     pos = link_midpoint(link)
#     ax.text(
#         pos[0],
#         pos[1],
#         label,
#         ha="center",
#         va="center",
#         fontsize=10,
#         color=color,
#         bbox=dict(facecolor="white", edgecolor=color, alpha=0.8),
#     )

#%% 1 plotting func
# def plot_step(lat, ops, title):
#     fig, ax = plt.subplots(figsize=(4, 4))
#     lat.plot()
#     ax.set_aspect("equal")
#     ax.set_title(title)
#
#     for links, label, color in ops:
#         if not isinstance(links, (list, tuple)):
#             links = [links]
#         for lgt in links:
#             if lat.link_exists(lgt[0], lgt[1]):
#                 link = lat.get_link(lgt[0], lgt[1])
#                 annotate_link(ax, link, label, color)
#
#     plt.show()

#
# def plot_lattice_with_ops(lat, ops, title=None):
#     """
#     ops entries:
#       (LGT_notation, label, color)                    -> 1-qubit
#       ((LGT_notation, LGT_notation), label, color)   -> 2-qubit arrow
#     """
#     fig, ax = plt.subplots(figsize=(4, 4))
#     lat.plot()
#     ax.set_aspect("equal")
#
#     if title is not None:
#         ax.set_title(title)
#
#     for item in ops:
#         links, label, color = item
#
#         # ---- 1-qubit operator ----
#         if isinstance(links, LGT_notation):
#             if lat.link_exists(links[0], links[1]):
#                 link = lat.get_link(links[0], links[1])
#                 pos = link_midpoint(link)
#                 ax.text(
#                     pos[0],
#                     pos[1],
#                     label,
#                     color=color,
#                     ha="center",
#                     va="center",
#                     fontsize=10,
#                     bbox=dict(facecolor="white", edgecolor=color, alpha=0.8),
#                 )
#
#         # ---- 2-qubit operator (arrow) ----
#         elif isinstance(links, tuple) and len(links) == 2:
#             a, b = links
#             if lat.link_exists(a[0], a[1]) and lat.link_exists(b[0], b[1]):
#                 link_a = lat.get_link(a[0], a[1])
#                 link_b = lat.get_link(b[0], b[1])
#
#                 start = link_midpoint(link_a)
#                 end   = link_midpoint(link_b)
#
#                 draw_arrow(ax, start, end, color=color)
#
#                 # operator label near the arrow midpoint
#                 mid = 0.5 * (start + end)
#                 ax.text(
#                     mid[0],
#                     mid[1],
#                     label,
#                     color=color,
#                     fontsize=9,
#                     ha="center",
#                     va="center",
#                     bbox=dict(facecolor="white", edgecolor="none", alpha=0.6),
#                 )
#
#         else:
#             raise TypeError("Invalid operator format")
#
#     plt.show()

# #%% 2 lattice setip
# lat = square_lattice(2, 5)
# rangex = range(lat.size[0])
# rangey = range(lat.size[1])
#
# diagonal = lambda i, j: [(j - i) % 4 == k for k in range(4)]

#%%   001

# ops = []
#
# for i in rangex:
#     for j in rangey:
#         if diagonal(i, j)[2]:
#             # ZY diagonal arrow
#             ops.append((
#                 (
#                     LGT_notation([i, j], 1),   # Z on this link
#                     LGT_notation([i, j], 0),   # Y on this link
#                 ),
#                 "ZY",
#                 "black"
#             ))
#
#             ops.append((
#                 (
#                     LGT_notation([i, j-1], 1),
#                     LGT_notation([i-1, j], 0),
#                 ),
#                 "ZY",
#                 "black"
#             ))
#
# plot_lattice_with_ops(lat, ops, title="Step 1")

# ops = []
#
# for i in rangex:
#     for j in rangey:
#         if diagonal(i, j)[2]:
#             ops.append((
#                 (
#                     LGT_notation([i, j], 1),
#                     LGT_notation([i, j], 0),
#                 ),
#                 "ZY",
#                 "black"
#             ))
#             ops.append((
#                 (
#                     LGT_notation([i, j-1], 1),
#                     LGT_notation([i-1, j], 0),
#                 ),
#                 "ZY",
#                 "black"
#             ))
#
# plot_step(lat, ops, "Step 1")
#
