# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 14:21:43 2025

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
from collections.abc import Iterable
import itertools
    


class LGT_notation(tuple):
    # a class for LGT representation of links as a tuple with site coordinates and link direction (int)
    def __new__(cls, coordinates, direction):
        if type(direction) == int and (isinstance(coordinates, Iterable) ) :
            coordinates = np.array(coordinates)
            return super().__new__(cls, [coordinates, direction])
        else:
            raise ValueError
            
    def __eq__(self, other):
        if not isinstance(other, LGT_notation):
           # don't attempt to compare against unrelated types
           return NotImplemented

        return (self[0] == other[0]).all() and self[1] == other[1]
            
    
class site:
    coordinates:np.array  # should be a numpy array 
    out_links = None
    in_links = None
    def __init__(self, coordinates):
        self.coordinates = np.array(coordinates)
        self.out_links = []
        self.in_links = []
    def __repr__(self):
        return f"<site :{self.coordinates}>"
        
    def dimension(self):
        return len(self.coordinates)
    
    def links(self):
        return [*self.in_links, *self.out_links]

def check_dimension(*sites):
    #returns an error if sites are not all the same dimension
        for s in sites:
            if not s.dimension() == sites[0].dimension():
                raise ValueError("Sites must have the same dimsnsion")
                
def vector_2_direction_index(vec):
    #TODO this is very ugly
    if len(vec)==1:
        return 0
    if len(vec) == 2:
        if (vec==np.array([1,0])).all():
            return 0
        elif (vec == np.array([0,1])).all():
            return 1
        else:
            Warning("vec must be a standard unit vector, returning nan")
            return np.nan
    if len(vec) == 3:
        if (vec==np.array([1,0,0])).all():
            return 0
        elif (vec == np.array([0,1,0])).all():
            return 1
        elif (vec == np.array([0,0,1])).all():
            return 2
        else:
            Warning("vec must be a standard unit vector, returning nan")
            return np.nan
    raise RuntimeError(
        "Direction must not be inferred from vectors when using PBC"
    )


class link:
    def __init__(self, site_from, site_to, direction):
        check_dimension(site_from, site_to)
        self.sites = (site_from, site_to)
        self._direction = direction

        site_from.out_links.append(self)
        site_to.in_links.append(self)

    def direction(self):
        return self._direction

    def vector(self):
        # purely geometric, may be non-unit for PBC
        return self.sites[1].coordinates - self.sites[0].coordinates

    def LGT_notation(self):
        return LGT_notation(self.sites[0].coordinates, self._direction)

    
 
class lattice:
    sites = None
    links = None
    
    def __init__(self, *sites):
        # input verification:
        check_dimension(*sites)
        self.sites = list(sites) #TODO maybe this can be a tuple. I don't know
        self.links = []
        missing_sites=[]
        for site in self.sites:
            for link in site.links():
                
                if link not in self.links:
                    self.links.append(link)
                #check that link does not connect to a site outside the lattice        
                for link_site in link.sites:
                    if not any (s==link_site for s in self.sites):
                        missing_sites.append(link_site)
                
        #inform user and propose to add them
        if len(missing_sites)>0:
            user_input = input(f"{len(missing_sites)} sites have links to sites that are not in the lattice. do you want to add them? (yes/no): ")
            if user_input.lower() in ["yes", "y"]:
               self.sites = [*self.sites, *missing_sites]
               print(f"added {len(missing_sites)} sites")
            else:
                print("lattice will contain links to sites that are not in the lattice. This will probably cause bugs")
        
                        
    def get_connectivity_mat(self):
        pass #TODO
    
    def get_site(self, coordinates, n=1):
        # return lattice site from coordinates. optionally return also the index in the sites list.
        # in the case that there is more than one wit the same coordinates, get the first n ( in a list, if n>1). 
        relevant_sites = [site for site in self.sites if (site.coordinates==coordinates).all()]
        n = min(len(relevant_sites),n)
        if len(relevant_sites)==0:
            raise ValueError("no site with the given coordinates")
            
        return relevant_sites[0:n] if n>1 else relevant_sites[0] 
    

    
    def get_link(self, *args, n=1):
        # return lattice link from two site objects,  or from LGT notation: tuple format: (coordinates:nparray, direction:int))
        if type(args[0])==site and type(args[1])==site :
            relevant_links = [link for link in self.links if (link.sites==(args[0], args[1]))]
        elif (type(args[0])==list or type(args[0])==np.ndarray ) and type(args[1]) == int:
            # args[ = np.array(args[0])        
            relevant_links = [link for link in self.links if link.LGT_notation() == LGT_notation(args[0], args[1])]
        
        n = min(len(relevant_links),n)
        return relevant_links[0:n] if n>1 else relevant_links[0] 
    
    def link_index(self, *args, n=1):
        if isinstance(args[0], link):
            l=args[0]
        else:
            l = self.get_link(*args,n)
        return self.links.index(l)
    
    def dimension(self):
        return self.sites[0].dimension()

    def plot(self, show_pbc=True):
        dim = self.dimension()
        size = np.array(self.size)

        # plot sites
        coords = np.array([s.coordinates for s in self.sites])
        if dim == 2:
            plt.scatter(coords[:, 0], coords[:, 1], zorder=3)

        segments = []

        for l in self.links:
            r0 = l.sites[0].coordinates.astype(float)
            r1 = l.sites[1].coordinates.astype(float)
            d = r1 - r0

            wrapped = False
            wrap_dir = None

            # detect wrapping
            for mu in range(dim):
                if abs(d[mu]) > 1:
                    wrapped = True
                    wrap_dir = mu
                    break

            if not wrapped or not show_pbc:
                segments.append([tuple(r0), tuple(r1)])
            else:
                # split wrapped link
                mu = wrap_dir
                sign = np.sign(d[mu])

                # first segment: r0 → boundary
                r_mid1 = r0.copy()
                r_mid1[mu] = size[mu] - 0.5 if sign > 0 else -0.5

                # second segment: opposite boundary → r1
                r_mid2 = r1.copy()
                r_mid2[mu] = -0.5 if sign > 0 else size[mu] - 0.5

                segments.append([tuple(r0), tuple(r_mid1)])
                segments.append([tuple(r_mid2), tuple(r1)])

        lc = mc.LineCollection(segments, linewidths=2)
        plt.gca().add_collection(lc)

        plt.gca().set_aspect('equal')
        plt.xlim(-.5, size[0]-.5)
        plt.ylim(-.5, size[1]-.5)
        plt.grid(True)

    def num_links(self):
        return len(self.links)
    
    def num_sites(self):
        return len(self.sites)


    

class square_lattice(lattice):
    size = None
    # boundary_sites = None
    # boundary_links = None
    
    
    def __init__(self,dim, size, pbc=False):
        if not isinstance(size, Iterable):
            size = [size for i in range(dim)]
        
        ranges = [range(size[n]) for n in range(dim)]
        
        self.sites =[]
        self.size=size
        for coordinates in itertools.product(*ranges):
            self.sites.append(site(coordinates))
            
        self.links = []
        unit_vecs = np.eye(dim)

        for s in self.sites:
            for mu in range(dim):
                step = unit_vecs[mu]
                target_coords = s.coordinates + step

                if pbc:
                    target_coords = target_coords % self.size

                try:
                    next_site = self.get_site(target_coords)
                    self.links.append(link(s, next_site, direction=mu))
                except ValueError:
                    pass

    def link_exists(self, *link_args):

       try:
           self.get_link(*link_args)
           return True
       except (IndexError):
           return False
    
    def plaquette_exists(self, s):
        if isinstance(s, site):
            coordinates = s.coordinates

        elif isinstance(s, Iterable):
            coordinates = s
        i = coordinates[0]
        j = coordinates[1]
        return self.link_exists([i, j], 0) * self.link_exists([i, j], 1) * self.link_exists([i + 1, j],1) * self.link_exists([i, j + 1], 0)


    
            
if __name__ =="__main__":
    lat2 = square_lattice(2,[3,4], pbc=False)
    lat2.plot()
    plt.show()