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
            


class link:
    sites = None # a directional list of size 2
    def __init__(self, *sites):
        check_dimension(*sites)
        self.sites = sites
        self.sites[0].out_links.append(self)
        self.sites[1].in_links.append(self)
    def __repr__(self):
        return f"<link : {self.sites}>"
    
    def vector(self):
        return self.sites[1].coordinates - self.sites[0].coordinates
    
    def dimension(self):
        return self.sites[0].dimension()
        
    def direction(self):
        ## works only for standard square lattice unit directions
        try:
            return vector_2_direction_index(self.vector())
        except Warning:
            Warning("direction index is only defined for standard unit vectors ( positive directions of a square lattice)")
                        
    def LGT_notation(self):
        return LGT_notation(self.sites[0].coordinates, self.direction())
    
 
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
    
    def plot(self):
        sites_coordinates = []
        for i in range(self.dimension()):
            sites_coordinates.append([site.coordinates[i] for site in self.sites])    
        
        sites_coordinates = np.array(sites_coordinates)
        links_coordinates = [ [tuple(link.sites[0].coordinates),tuple(link.sites[1].coordinates)] for link in self.links]
        #TODO: 3D
        if self.dimension()==2:
            xx,yy =     np.meshgrid(sites_coordinates[0],sites_coordinates[1])
            plt.scatter(xx,yy)
            
        plt.gca().add_collection(mc.LineCollection(links_coordinates, linewidths=2))
    
    def num_links(self):
        return len(self.links)
    
    def num_sites(self):
        return len(self.sites)


    

class square_lattice(lattice):
    size = None
    # boundary_sites = None
    # boundary_links = None
    
    
    def __init__(self,dim, size):
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
            for i in range(dim):
                try:
                    next_site = self.get_site(s.coordinates+unit_vecs[i])
                    self.links.append(link(s, next_site))
                except ValueError:
                    # do skip the cases where next site does not exist
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
    lat2 = square_lattice(2,10)
    lat2.plot()