#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 13:55:04 2022

@author: lara
"""

import numpy as np
r3 = np.sqrt(3)
import matplotlib.pyplot as plt
import json



directory = '/Users/lara/Documents/SelfAssembly2/Lattice/'



n_particles_max = 6

    

def rotation_matrix_axis_angle(x,y,z, cos):
    """(x,y,z) must be unitary"""
    sin = np.sqrt(1-cos**2)
    R = [[cos+x**2*(1-cos),    x*y*(1-cos) - z*sin,   x*z*(1-cos) + y*sin  ],
         
         [x*y*(1-cos)+z*sin,   cos + y**2*(1-cos),    y*z*(1-cos) - x*sin ],
         
         [z*x*(1-cos)-y*sin,    z*y*(1-cos)+x*sin,     cos + z**2*(1-cos) ]]
    
    R= np.array(R)    
    return(R)#/np.linalg.norm(R))
def get_axis_cos(u,v):
    cos = u.dot(v)
    axis = np.cross(u,v)
    return(axis, cos)
    
def rotation_matrix(u,v):
    """u and v are unitary vectors
    v = R.u
    """
    
    axis, cos = get_axis_cos(u,v)
    norm = np.linalg.norm(axis)
    if norm >1e-5 :
        x,y,z = axis/np.linalg.norm(axis)
    else :
        x, y, z= 0,0,0
    return(rotation_matrix_axis_angle(x, y, z, cos))
    


class LatticeParticle:   
    
    def __init__(self, lattice_name, generators, edges, neighbors ,basis,face_type
                 , rotation_axis, rotation_angle):
        """
        orientation = 0 means no particle
        edge = 0 means structrue of one single particle and not two
        Therefore, orientation and structures are indexed from one 
        for two particles scaffold
        """
        
        self.lattice = lattice_name
        self.generators = np.array(generators)
        self.edges = np.array(edges, dtype=float) 
        
        self.n_edges = len(edges)
        self.neighbors = np.array(neighbors, dtype=float) 
        self.n_neighbors = len(edges)*2
        self.neighbors_cartesian = np.zeros((self.n_neighbors,3))
        
        self.n_elem_rotations = len(rotation_angle)
        self.elementary_rotations = np.zeros((self.n_elem_rotations, 3, 3))
        for i in range (self.n_elem_rotations):
            cos = np.cos(rotation_angle[i])
            x,y,z = rotation_axis[i]
            self.elementary_rotations[i] = rotation_matrix_axis_angle(x, y, z, cos)
            
        
        for i in range (self.n_neighbors):
            self.neighbors_cartesian[i] = self.lattice_to_cartesian(self.neighbors[i])
            norm = np.linalg.norm(self.neighbors_cartesian[i])
            self.neighbors_cartesian[i] = self.neighbors_cartesian[i]/norm
            
            
            
        
        """The reference edges are the contact direction for which we will define 
        the structures"""
        self.n_face_type = len(np.unique(face_type))        
        self.face_type = np.array(face_type)
        self.reference_edges_id_left = []
        self.reference_edges_id_right = []
        for i in range (self.n_face_type):
            edge_ref = np.where(self.face_type==np.unique(face_type)[i])[0][0] +1
            self.reference_edges_id_left.append( edge_ref )
            self.reference_edges_id_right.append( self.find_opposite_edge_id(1))
            
        
        
        """Find all the orientations from elementary rotations"""
        self.orientations = np.array([np.arange(0, self.n_neighbors, 1, dtype=int)])
        self.generate_all_orientations()
        self.n_orientations = len(self.orientations)
        
        """Find all structures for one type of particle"""
        # self.map_exchange_left_right = {} #o1, o2, edge is equivalent to o2', o1', -edge
        self.enumerate_one_particle_structures()
        self.enumerate_two_particle_structures_empty_full()
        self.enumerate_two_particles_structures_full_full()
        
        
        print(self.lattice,': Successfully built,', 
              self.n_neighbors, 'neighbors,', 
              self.n_orientations, 'orientations' )
        
        """Several types of particle"""
        
        self.all_structures = {}
        self.Nstructures = {} 
        self.detail_Nstructures = {} #how much structures of type 1particle, 
        # 2particle empty-full, and 2 particles full-full
        for n_particles in range (1,n_particles_max+1):
            self.create_complete_structure_map(n_particles)
            print(n_particles, 'particle type ->', self.Nstructures[n_particles])
        
      
   
    
    
    def isIdentity(self, R):
        if np.all(np.abs(R-np.eye(3))<1e-5):
            return(True)
        else:
            return(False)
        
    def orientationExists(self, anOrientation):
        searchZeros = np.abs(self.orientations -anOrientation)
        indices = np.where((searchZeros==0).all(axis=1))[0]
        if len(indices)>0:
            return(True)
        else:
            return(False)
        
    def generate_all_orientations(self):
        self.matrices = [np.eye(3)]
        n_orientations_found = 1
        index_to_rotate = 0
        while index_to_rotate<n_orientations_found:
            for i in range(self.n_elem_rotations):
                # print(i)
                R = self.elementary_rotations[i]
                matrix = R
                
                previous_matrix_applied = self.matrices[index_to_rotate] #rotation matrix applied so far
                while not self.isIdentity(matrix):
                    new_orientation = self.rotate(self.orientations[index_to_rotate], matrix)
                    
                    if not self.orientationExists(new_orientation):
                        # we found a new state
                        self.orientations = np.concatenate((self.orientations, np.array([new_orientation])), axis=0)
                        self.matrices.append(matrix.dot(previous_matrix_applied))
                    matrix = matrix.dot(R)
                n_orientations_found = len(self.orientations)
                index_to_rotate+=1
                
            
        
    
    def predict_n_structures(self, n_particles):
        n=1 #empty-empty
        n+= self.n_neighbors*n_particles #empty-full
        
        for t in range(self.n_face_type) :
            n_this_face = len(np.where(self.face_type==t)[0])*2 * n_particles
            n += n_this_face*(n_this_face+1)/2
        return(int(n))
            
    def lattice_to_cartesian(self, coordinate):
        d =  len(self.generators)
        x,y,z = 0,0,0
        for i in range(d):
            dx, dy, dz = coordinate[i] * self.generators[i]
            x+=dx
            y+=dy
            z+=dz
        return(x,y,z)
    
        
    
    def plot_orientations(self, rotations_to_plot=None, reindex=True):
        colors = ['blue', 'red', 'green', 'orange', 'pink', 'purple', 'cyan', 'grey', 'black', 'magenta', 'gold', 'darkblue']
        
        if rotations_to_plot is None :
            # rotations_to_plot = range(self.n_rotations)
            
            rotations_to_plot = range(self.n_orientations)
        
        
            
        f = plt.figure(figsize=(7,7))
        # f, axes = plt.subplots((1, self.n_rotations))
        n1 = int(np.sqrt(len(rotations_to_plot)))
        n2 = len(rotations_to_plot)//n1 
        print(n1,n2)
        
        index_plot = 1
        for rot in rotations_to_plot :
            
            axes = f.add_subplot(n1,n2,index_plot, projection='3d')
            index_plot+=1
            for k in range(self.n_neighbors) :
                x, y, z = self.lattice_to_cartesian(self.neighbors[self.orientations[rot]][k])
                norm = np.sqrt(x**2+y**2+z**2)
                x,y,z=x/norm, y/norm, z/norm
                axes.plot([0,x], [0,y], [0,z], marker='o', color=colors[k])
                if reindex :
                    axes.text(x,y,z, str(k+1), color=colors[k])
                else :
                    axes.text(x,y,z, str(k), color=colors[k])
            if reindex :
                axes.set_title('Orientation '+str(rot+1), fontsize=8)
            else :
                axes.set_title('Orientation '+str(rot), fontsize=8)
        f.suptitle(self.lattice)
        
        
        
    def find_rotation_id(self, an_orientation):
            sames = np.where((an_orientation==self.orientations).all( axis=1))[0]  
            return(sames)
        
    def get_particle_state(self, particle, orientation):
        if orientation==0:
            return(0)
        else :
            return(particle*self.n_orientations+orientation)
    
    def get_particle_id_orientation(self, particle_state):
        if particle_state ==0 :
            particle_id=0
            orientation=0
        else :
            particle_id = (particle_state-1)//self.n_orientations
            orientation = (particle_state-1)%self.n_orientations+1
        return(particle_id, orientation)
        
    def check_rotations_iteration(self, printInfo=True):
        """Iterate rotations until you find identity"""
        nmax = 20
        all_good = True
        forgotten_rotations = []
        seed_forgotten, n_iterations = [], []
        
        
            
        for i in range(self.n_rotations):
            """Check if repeating n times goes back to identity"""
            a_rotation = self.rotations[i]
            # print(a_rotation)
            my_seq = self.rotations[0].copy()
            back = False
            n=0
            while not back and n < nmax:
                my_seq = my_seq[a_rotation]
                back = np.all(my_seq ==self.rotations[0] )
                n+=1
                sames = self.find_rotation_id( my_seq)
                if len(sames)==0 :
                    forgotten_rotations.append(my_seq)
                    seed_forgotten.append(i)
                    n_iterations.append(n)
            """Check if is unique in the roations list"""
            sames = self.find_rotation_id( a_rotation)
            unique= False
            if len(sames)==1:
                unique=True
            if printInfo :
                print('rotation is identity after '+str(i)+ ' '+str(back) +' after '+str(n)+' and is unique '+str(unique))
            all_good = all_good&unique&back
            
        if len(forgotten_rotations)!=0 :
            print('Forgotten rotations')
            print(forgotten_rotations)
            print(seed_forgotten)
            print(n_iterations)
            all_good = False
        return(all_good)
    

    

    def find_face_id(self, orientation, edge_id):
        """If particle is in given orientation, 
        what face is facing the edge number edge_id
        If the orientation is 0, then the face is -1"""
        if orientation==0 :
           face_id = -1
        else :
            
            face_id = np.where(self.orientations[orientation-1] == edge_id-1)[0][0]+1
            
        return(face_id)
    
    def find_edge_id(self, edge_vector, latticeSpace=True):
        """Returns the index of a given edge, in lattice coordinate"""
        # print(edge_vector)
        if latticeSpace :
            searchZeros = np.abs(self.neighbors -edge_vector)
        else :
            searchZeros = np.abs(self.neighbors_cartesian -edge_vector)
        try :
            edge_id = np.where((searchZeros<1e-4).all(axis=1))[0][0]+1
        except :
            print(edge_vector, "not found")
        return(edge_id)
        
    
    
    def rotate(self, current_orientation, matrix):
        new_orientation = np.zeros(self.n_neighbors, dtype=int)
        for face in range(self.n_neighbors):
            current_edge_id = current_orientation[face]
            vector_edge = self.neighbors_cartesian[current_edge_id]
            new_vector_edge = matrix.dot(vector_edge)
            new_edge_id = self.find_edge_id(new_vector_edge, latticeSpace=False)-1
            # print(current_edge_id,vector_edge,new_vector_edge,  new_edge_id)
            new_orientation[face] = new_edge_id
        return(new_orientation)
    
    def rotate_id(self, current_orientation_id, matrix):
        current_orientation = self.orientations[current_orientation_id-1]
        new_orientation = self.rotate(current_orientation, matrix)
        new_orientation_id = self.find_rotation_id(new_orientation)[0]+1
        return(new_orientation_id)

    
    # def find_rotation_matrix_from_edges(self, edeg_id_init, edge_id_final) :
    #     vector_edge_0 = self.neighbors_cartesian[edeg_id_init-1]
    #     vector_edge_1 = self.neighbors_cartesian[edge_id_final-1]
    #     matrix = rotation_matrix(vector_edge_0, vector_edge_1)
    #     return(matrix)
    

    
    def find_equivalent_orientation_1particle(self, orientation, face):
        """What is the orientation where the chosen face is at edge=1"""
        
        edge_ref_left = 0
        orientation_reference = 0
        
        """We look for an orientation where the contact face is at the edge of reference"""
        for i in range(len(self.reference_edges_id_left)):
            reference_edge = self.reference_edges_id_left[i]
            possible_new_orientation_particle = np.where(self.orientations[:,face-1]==reference_edge-1)[0]
            if len(possible_new_orientation_particle)>0:
                edge_ref_left = reference_edge
                orientation_reference = possible_new_orientation_particle[0]+1
        return(edge_ref_left, orientation_reference)
        
    
    def find_rotation_matrix(self, orientation_init, orientation_final):
        matrix_init_to_1 = np.linalg.inv(self.matrices[orientation_init-1])
        matrix_1_to_final = self.matrices[orientation_final-1]
        matrix_init_to_final= matrix_1_to_final.dot(matrix_init_to_1)
        return(matrix_init_to_final)
    
    def find_equivalent_orientations(self, orientation1, orientation2, edge):
        """Given the edge that is the direction of the contact, and the current 
        orientations of the particle, what are the equivalent orientations if you 
        have the contact in direction edge=1"""
        
        face1 = self.find_face_id(orientation1, edge)
        opposite_edge = self.find_opposite_edge_id(edge)
        face2 = self.find_face_id(orientation2, opposite_edge)
        
        """First particle on the left of the reference contact-> 1L (1 is Left)"""
        edge_ref_left_1L, new_o1_1L = self.find_equivalent_orientation_1particle(orientation1, face1)
        matrix_1L = self.find_rotation_matrix (orientation1,new_o1_1L )
        new_o2_1L = self.rotate_id(orientation2, matrix_1L)
       
        
        """Second particle on the right of the reference contact -> 1R (1 is Right)"""
        edge_ref_left_1R, new_o1_1R = self.find_equivalent_orientation_1particle(orientation2, face2)
        matrix_1R = self.find_rotation_matrix (orientation2,new_o1_1R )
        new_o2_1R = self.rotate_id(orientation1, matrix_1R)
        
        
        """Choose between two options : the one where new_o1 is the smallest"""
        if new_o1_1L <= new_o1_1R :
            newo1, newo2, edge_ref = new_o1_1L, new_o2_1L, edge_ref_left_1L 
            first_is_left = True
        else :
            
            newo1, newo2, edge_ref = new_o1_1R, new_o2_1R, edge_ref_left_1R
            first_is_left = False
            
        is_dimer = (new_o1_1L==new_o1_1R) and (new_o2_1L==new_o2_1R)
        # if not (new_o1_1L, new_o2_1L, edge_ref_left_1L) in self.map_exchange_left_right.keys():            
        #     self.map_exchange_left_right[(new_o1_1L, new_o2_1L, edge_ref_left_1L)] = (new_o1_1R, new_o2_1R, edge_ref_left_1R)
        #     self.map_exchange_left_right[(new_o1_1R, new_o2_1R, edge_ref_left_1R)] = (new_o1_1L, new_o2_1L, edge_ref_left_1L)
        
        return(newo1, newo2, edge_ref, first_is_left,is_dimer)
    
    # def find_equivalent_empty_full_structure(self, orientation, edge):
    #     """What is the equivalent o2 and edge such that (o1, 0, edge) <=> (0, o2, -edge)"""
    #     opposite_edge = self.find_opposite_edge_id(edge)
    #     face = self.find_face_id(orientation, edge)
        
    #     new_orientation = np.where(self.orientations[:,face-1]==opposite_edge-1)[0][0] +1
    #     return(new_orientation, opposite_edge)
        
    def enumerate_two_particles_structures_full_full(self):
        """for 2 particles with orientation>1, maps it to the reference structure"""
        self.reference_structures_ff = []
        self.structure_map_2p_ff = {}
        
        self.inversions = {}
        self.dimers = {}
        
        
        for edge in range(1,self.n_neighbors+1):
            for o1 in range(1, self.n_orientations+1):
                for o2 in range(1, self.n_orientations+1):
                    new_o1, new_o2, new_edge, first_is_left, is_dimer = self.find_equivalent_orientations(o1, o2, edge)
                    if (new_o1, new_o2, new_edge) not in self.reference_structures_ff :
                        s = len(self.reference_structures_ff)
                        self.reference_structures_ff.append((new_o1, new_o2, new_edge))
                    else :
                        s = self.reference_structures_ff.index((new_o1, new_o2, new_edge))
                    self.structure_map_2p_ff[(o1, o2, edge)] = s
                    self.inversions[(o1, o2, edge)] = (not first_is_left)
                    self.dimers[(o1,o2,edge)] = is_dimer
        
        self.Nstructures_2p_ff = len(self.reference_structures_ff)
        self.verify_one_map(self.structure_map_2p_ff,'structure_map_2p_ff, 1type')
            
    def enumerate_one_particle_structures(self):
        """We want to leave the possibility that a single orientation has an energy 
        (if there is a magnetic field for instance)
        To encode this as the rest, we use the following convention
        (orientation, orientation, edge=0, structure)
        There is (n_orientations+1) extra structures that are added 
        at the beginnning of the structure list"""
        
        self.structure_map_1p = {}     
  
        for o in range(self.n_orientations+1):
            self.structure_map_1p[(o,o,0)] = o
            
        self.Nstructures_1p = len(self.structure_map_1p)
        self.verify_one_map(self.structure_map_1p,'structure_map_1p, 1type')
            
            
        
    def enumerate_two_particle_structures_empty_full(self):
        """We look at the empty-full pairs"""
        
            
        self.reference_structures_ef = []
        self.structure_map_2p_ef = {}
        
        """Empty-Empty"""
        for edge in range(1, self.n_neighbors+1):
            self.structure_map_2p_ef[(0,0,edge)]=0
        self.reference_structures_ef.append((0,0,1))
        
        """Full site on the left of the reference contact"""
        for edge in range(1, self.n_neighbors+1):
            for o1 in range (1, self.n_orientations+1):
                face = self.find_face_id(o1, edge)
                new_edge, new_o1 = self.find_equivalent_orientation_1particle(o1, face)
                
                if (new_o1, 0, new_edge) not in self.reference_structures_ef :
                    s = len(self.reference_structures_ef)
                    self.reference_structures_ef.append((new_o1, 0, new_edge))
                else :
                    s = self.reference_structures_ef.index((new_o1, 0, new_edge))
                self.structure_map_2p_ef[(o1, 0, edge)] = s
                    
        """Full site on the right of the reference contact"""
        for edge in range(1, self.n_neighbors+1):
            for o2 in range (1, self.n_orientations+1):
                opposite_edge = self.find_opposite_edge_id(edge)
                face = self.find_face_id(o2, opposite_edge)
                new_edge, new_o1 = self.find_equivalent_orientation_1particle(o1, face)
                s = self.reference_structures_ef.index((new_o1, 0, new_edge))
                self.structure_map_2p_ef[(0, o2, edge)] = s
                    
                
        # o2, opposite_edge = self.find_equivalent_empty_full_structure(o1, edge)
        # self.structure_map_2p_ef[(0, o2, edge)] = s
                
        self.Nstructures_2p_ef = len(self.reference_structures_ef)       
        self.verify_one_map(self.structure_map_2p_ef,'structure_map_2p_ef, 1type')
        
        
        
    def find_opposite_edge_id(self, edge_id):
        """If edge_id is the id of the vector pointing towards the particle
        what is the id of the same vector pointing outwards the particle"""
        
        
        if edge_id ==0 :
            return(0)
        else :
            edge = self.neighbors[edge_id-1]
            opposite_edge = -edge
            searchZeros = np.abs(self.neighbors -opposite_edge)
            opposite_edge_id = np.where((searchZeros==0).all(axis=1))[0][0]
            opposite_edge_id +=1
            return(opposite_edge_id)
            # return((edge_id+self.n_neighbors//2)%self.n_neighbors)
        
        
    # def faces_in_contact(self, orientation1, orientation2, edge_id):
    #     """Edge is the factor starting at particle1 and ending at particle2"""
    #     face1 = self.find_face_id(orientation1, edge_id)
    #     face2 = self.find_face_id(orientation2, self.find_opposite_edge_id(edge_id))
    #     return(face1, face2)
    
    
    
    
    def verify_one_map(self, aMap, map_name):
        """A few structures were redundant (the dimeric one for different particles) 
        For ex, (p1, o1)(p2,-o1)= s1 and (p2, o1)(p1, -o1)=s2, but now s2=s1 and
        there is no structrue for s2"""
        
        
        all_items = np.array(list(aMap.values()))
        unique_items = np.unique(all_items)
        n = len(unique_items)
        if np.all(unique_items == np.arange(0, n, 1)):
            print(self.lattice, 'map',map_name,'ok')
        if np.any(unique_items != np.arange(0, n, 1)) :
            print(self.lattice, 'map',map_name,'needs index shift')
            newMapping = {}
            for s in range(n):
                """Find the shift of the indices"""
                old_s = unique_items[s]
                newMapping[old_s] = s
            for key in aMap.keys():
                """Apply it to the structures"""
                old_s = aMap[key]
                aMap[key] = newMapping[old_s]
                
                
            
            # print(len(unique_items), 'unique items and goes from', 
            #       np.min(unique_items), 'to', np.max(unique_items))
   
    def verify_maps(self) :
        i=0
        for aMap in [self.structure_map_1p, self.structure_map_2p_ef, self.structure_map_2p_ff]:
            self.verify_one_map(aMap, str(i))
            i+=1
    
    
    
    
    
    
    """
    !
    !
    !
    !
    Several types of particles 
    !
    !
    !
    """
    def get_particle_pair_id(self, p1, p2, n_particles):
        """A pair of particle is associated with an id between 0 
        and n_particles*(n_particles+1)/2"""
        """
        p1/p2 0  1  2  3
        0     0  1  2  3
        1     4  5  6  7
        2     8  9 10 11
        etc
        """
        # return( np.ravel_multi_index((p1,p2),(n_particles,n_particles)) )
        return(p1*n_particles + p2)
    
    def find_structure_2p_ntypes(self, p1, p2, o1, o2, edge, n_particles):
        """Gets the corresponding structure for several particles"""
        s = self.structure_map_2p_ff[(o1, o2, edge)]
        inversion = self.inversions[(o1, o2, edge)]
        dimer = self.dimers[(o1,o2,edge)]
        if not inversion :
            particles_id = self.get_particle_pair_id(p1,p2,n_particles)
        else :
            particles_id = self.get_particle_pair_id(p2,p1,n_particles)
            
        if dimer :
            particles_id = self.get_particle_pair_id(min(p1,p2),max(p1,p2),n_particles)
        s_tot = (self.Nstructures_2p_ff)*particles_id + s
        return(s_tot)

    def all_two_particles_structures_full_full(self, n_particles):
        """There are several types of particles now"""
        self.structure_map_2p_ff_several = {}
        # self.inversions_several = {}
        for p1 in range(n_particles):
            for p2 in range (n_particles):
                for o1 in range (1, self.n_orientations+1):
                    for o2 in range (1, self.n_orientations+1):
                        for edge in range(1, self.n_neighbors+1):
                            state1 = self.get_particle_state(p1, o1)
                            state2 = self.get_particle_state(p2, o2)
                            s_tot = self.find_structure_2p_ntypes(p1, p2, o1, o2, edge,n_particles)
                            self.structure_map_2p_ff_several[(state1, state2, edge)] = s_tot
                            # self.inversions_several[(state1, state2, edge)] = inversion
        self.verify_one_map(self.structure_map_2p_ff_several,'structure_map_2p_ff, n='+str(n_particles))
        
        
    def find_structure_2p_ntypes_ef(self, particle, o1,  o2,  edge):
        """Gets the corresponding structure for several particles, 
        either o1 or o2 has to be zero """
        s = self.structure_map_2p_ef[(o1, o2, edge)]
        if o1!=0 or o2!=0 :
            s_tot = particle*(self.Nstructures_2p_ef-1) + s
        else :
            s_tot = 0
        return(s_tot)
    
    def all_two_particles_structures_empty_full(self, n_particles):
        """There are several types of particles now"""
        self.structure_map_2p_ef_several = {}
        # for edge in range(1, self.n_neighbors+1):
        #     self.structure_map_2p_ef_several[0,0,edge] = 0
            
        for p1 in range(n_particles):
            for (o1,o2, edge) in self.structure_map_2p_ef.keys():
                state1 = self.get_particle_state(p1, o1)
                state2 = self.get_particle_state(p1, o2)
                s_tot = self.find_structure_2p_ntypes_ef(p1, o1, o2, edge)
                self.structure_map_2p_ef_several[(state1, state2, edge)] = s_tot
        self.verify_one_map(self.structure_map_2p_ef_several,'structure_map_2p_ef, n='+str(n_particles))
        
            
    
    
    def find_structure_1p_ntypes(self, particle, orientation):
        s = self.structure_map_1p[(orientation, orientation, 0)]
        if orientation==0 :
            s_tot = 0
        else :
            s_tot = particle*(self.Nstructures_1p-1) + s
        return(s_tot)
                    
    
    def all_one_particle_structures(self, n_particles):
        self.structure_map_1p_several = {}
        for p in range(n_particles):
            for o in range(0, self.n_orientations+1):
                state = self.get_particle_state(p, o)
                s_tot = self.find_structure_1p_ntypes(p,o)
                self.structure_map_1p_several[(state, state, 0)] = s_tot
        self.verify_one_map(self.structure_map_1p_several,'structure_map_1p, n='+str(n_particles))
          
    
    
    
    
    def create_complete_structure_map(self, n_particles):
        """Create a list of (orientation1, orientation2, edge, structure)"""
        
        
        self.all_one_particle_structures(n_particles)
        self.all_two_particles_structures_empty_full(n_particles)     
        self.all_two_particles_structures_full_full(n_particles)
        
        self.all_structures[n_particles] = {}
        
        N1 = len(np.unique(list(self.structure_map_1p_several.values())))
        N2ef = len(np.unique(list(self.structure_map_2p_ef_several.values())))
        N2ff = len(np.unique(list(self.structure_map_2p_ff_several.values())))
        
        self.detail_Nstructures[n_particles] = [ N1, N2ef, N2ff]
        # print('n_particles', n_particles)
        # print( 
        #       '1pStructures', N1, 
        #       '2pStructure_empty-full', N2ef,
        #       '2pStructure_full-full', N2ff,
        #       'total', N1+N2ef+N2ff)
        print('(',N1, N2ef, N2ff,')')
        
        add = 0
        for (s1, s2, edge) in self.structure_map_1p_several.keys():
            s = self.structure_map_1p_several[(s1,s2,edge)]+add
            self.all_structures[n_particles][(s1,s2,edge)] =  s
        
        add = N1
        for (s1, s2, edge) in self.structure_map_2p_ef_several.keys():
            s = self.structure_map_2p_ef_several[(s1,s2,edge)]+add
            self.all_structures[n_particles][(s1,s2,edge)] = s
            
        add = N1+N2ef
        for (s1, s2, edge) in self.structure_map_2p_ff_several.keys():
            s = self.structure_map_2p_ff_several[(s1,s2,edge)]+add
            self.all_structures[n_particles][(s1,s2,edge)] = s
            
            
            
        self.verify_one_map(self.all_structures[n_particles], '-> total map')
        self.Nstructures[n_particles] = len(np.unique(list(self.all_structures[n_particles].values())))
        

    
    def save_lattice_file(self, directory):
        
        #CHECK IF NOT EXIST
        
        myDict = {}
        myDict['Generators'] = self.generators.tolist()
        myDict['Edges']=self.edges.tolist()
        myDict['N_Edges']=self.n_edges
        
        
        myDict['Neighbors']=self.neighbors.tolist()
        myDict['N_Neighbors']=self.n_neighbors
        
        myDict['N_Orientations']=self.n_orientations
        myDict['Orientations']=self.orientations.tolist()
        
        myDict['Faces_type']=self.face_type.tolist()
        for n_particles in range(1, n_particles_max+1):
            myDict['N_Structures_'+str(n_particles)]=self.Nstructures[n_particles]
            myDict['detail_N_Structures_'+str(n_particles)] = self.detail_Nstructures[n_particles]
        jsonString = json.dumps(myDict)
        
        filename = directory+self.lattice+".json"
        try : 
            open(filename, "r")
            fileExists = True
        except :
            fileExists = False
        
        if not fileExists :
            jsonFile = open(filename, "w")
            jsonFile.write(jsonString)
            jsonFile.close()
        else :
            print(filename+' already exists, we will not erase it')
        
    def save_structure_file(self, directory, n_particles):
        
        #CHECK IF NOT EXIST
        # self.create_structure_map(n_particles, print_check=False)
        filename = directory+self.lattice+'_structures_'+str(n_particles)+'particles.txt'
        
        
        try : 
            open(filename, "r")
            fileExists = True
        except :
            fileExists = False
        
        if not fileExists :
            file = open(filename, 'w')
            string =""
            for key in list(self.all_structures[n_particles].keys()):
                s1, s2, edge = key
                structure = self.all_structures[n_particles][key]
                string+=str(s1)+" "+str(s2)+" "+str(edge)+" "+str(structure)
                string+="\n"
            string+="#END"
            file.write(string)
            file.close()   
        else :
            print(filename+' already exists, we will not erase it')
        

myLattices = []

# r3 = np.sqrt(3)
# r22 = np.sqrt(2)

name = 'FaceCenteredCubic'
"""Generators : from real space cartesian coordinate to lattice coordinate"""
generators = [(1/2, 1/2, 0), (-1/2, 1/2, 0), (0, 1/2, 1/2)]
"""Edges : translation in the lattice space"""
edges = [(1,0,0), (0,1,0), (0,0,1), (0, -1, 1), (-1, 0, 1), (-1, -1, 1) ]
basis = [0,1,2]
neighbors = [np.array(e) for e in edges]+[-np.array(e) for e in edges]
face_type = [0,0,0,0,0,0]


face_center = [1,0,0]
# edge_center = 
rotation_axis = [[1,0,0] ]
rotation_angle = [72/180*np.pi]

rotations = [#groups (0,1,6,7), (3,2,4,5), (11,9,8,10)
            [0,1,2,3,4,5,6,7,8,9,10,11],
            [7,0,3,5,2,4,1,6,9,11,8,10], 
            [6,7,5,4,3,2,0,1,11,10,9,8], 
            [1,6,4,2,5,3,7,0,10,8,11,9],
             
              # #groups(3,10,9,4),  (0,11,1,2), (7,8,6,5)
            [11,2 ,0 ,10,3 ,7 ,5 ,8 ,6 ,4 ,9 ,1 ],
            [1 ,0 ,11,9 ,10,8 ,7 ,6 ,5 ,3 ,4 ,2 ],
            [2 ,11,1 ,4 ,9 ,6 ,8 ,5 ,7 ,10,3 ,0 ],
             
              # #groups (2,11,8,5), (1,9,6,4), (0,10,7,3)
            [10,9,11,0,1,2,4,3,5,6,7,8],
            [7,6,8,10,9,11,1,0,2,4,3,5],
            [3,4,5,7,6,8,9,10,11,1,0,2],
             
            #groups (1,2,3,7,8,9), (4,5,6), (11,0,10)
            [11,9,1,2,5,6,4,3,7,8,0,10],
            [10,2,3,7,6,4,5,8,9,1,11,0], 
            [11,3,7,8,5,6,4,9,1,2,0,10],
            [0,7,8,9,4,5,6,1,2,3,10,11],
            [10,8,9,1,6,4,5,2,3,7,11,0],
            
            #groups (0,3,5,6,9,11) , (2,4,1), (10,7,8)]
            [3,2,4,5,1,6,9,8,10,11,7,0],
            [11,4,1,0,2,3,5,10,7,6,8,9],
            [9,2,4,11,1,0,3,8,10,5,7,6],
            [6,1,2,9,4,11,0,7,8,3,10,5],
            [5,4,1,6,2,9,11,10,7,0,8,3],
            
            
            #groups (0,2,4,6,8,10) , (3,7,5), (11,1,9)]
            [2,9,4,5,6,7,8,3,10,11,0,1],
            [10, 11,  0,  7,  2,  3,  4,  5,  6,  1,  8,  9],
            [ 8,  9, 10,  5,  0,  7,  2,  3,  4, 11,  6,  1],
            [ 6,  1,  8,  3, 10,  5,  0,  7,  2,  9,  4, 11],
            [ 4, 11,  6,  7,  8,  3, 10,  5,  0,  1,  2,  9]
            
    ]

# FCC = LatticeParticle(name, generators, edges, neighbors, basis, face_type, rotation_axis, rotation_angle)
# myLattices.append(FCC)

name = 'Triangular'
generators = [(1,0, 0), (1/2, r3/2, 0)]
edges = [(1,0,0), (0,1,0), (-1, 1, 0)]
neighbors = [np.array(e) for e in edges]+[-np.array(e) for e in edges]
basis = [0,1]
rotation_axis = [[0,0,1]]
rotation_angle = [np.pi/3]

rotations = [[0,1,2,3,4,5],
              [5,0,1,2,3,4],
              [4,5,0,1,2,3],
              [3,4,5,0,1,2],
              [2,3,4,5,0,1],
              [1,2,3,4,5,0]
              ]

face_type = [0,0,0]
Triangular = LatticeParticle(name, generators, edges, neighbors, basis, face_type, rotation_axis, rotation_angle)
myLattices.append(Triangular)


name = 'Hexagonal'
generators = [(1,0, 0), (1/2, r3/2, 0), (0, 0, 1)]
edges = [(1,0,0), (0,1,0), (-1, 1, 0), (0,0,1)]
neighbors = [np.array(e) for e in edges]+[-np.array(e) for e in edges]
basis = [0,1,3]
face_type = [0,0,0,1]

rotation_axis = [[0,0,1], [1,0,0]]
rotation_angle = [np.pi/3, np.pi]


rotations = [[0,1,2,3,4,5,6,7],
              [6,0,1,3,2,4,5,7],
              [5,6,0,3,1,2,4,7],
              [4,5,6,3,0,1,2,7],
              [2,4,5,3,6,0,1,7],
              [1,2,4,3,5,6,0,7],
             
              [6,5,4,7,2,1,0,3],
              [0,6,5,7,4,2,1,3],
              [1,0,6,7,5,4,2,3],
              [2,1,0,7,6,5,4,3],
              [4,2,1,7,0,6,5,3],
              [5,4,2,7,1,0,6,3]]

Hexagonal = LatticeParticle(name, generators, edges, neighbors, basis, face_type, rotation_axis, rotation_angle)
myLattices.append(Hexagonal)


name = 'Square'
generators = [(1,0,0), (0,1,0)]
edges = [(1,0,0), (0,1,0)]
neighbors = [np.array(e) for e in edges]+[-np.array(e) for e in edges]
basis = [0,1]
rotation_axis = [[0,0,1]]
rotation_angle = [np.pi/2]

face_type = [0,0]
rotations = [[0,1,2,3],
             [3,0,1,2],
             [2,3,0,1],
             [1,2,3,0] ]
Square = LatticeParticle(name, generators, edges, neighbors, basis, face_type, rotation_axis, rotation_angle)
myLattices.append(Square)

name = 'Cubic'
generators = [(1,0,0), (0,1,0), (0,0,1)]
edges = [(1,0,0), (0,1,0), (0,0,1)]
neighbors = [np.array(e) for e in edges]+[-np.array(e) for e in edges]


rotation_axis = [[0,0,1], [1,0,0]]
rotation_angle = [np.pi/2, np.pi/2]

basis = [0,1,2]

rotations = [#(0,1,3,4)
            [0,1,2,3,4,5],
            [0,5,1,3,2,4],
            [0,4,5,3,1,2],
            [0,2,4,3,5,1],
            
            [1,0,5,4,3,2],
            [1,5,3,4,2,0],
            [1,3,2,4,0,5],  
            [1,2,0,4,5,3],
             
            [2,1,3,5,4,0],
            [2,3,4,5,0,1],
            [2,4,0,5,1,3],
            [2,0,1,5,3,4],
            
            [3,4,2,0,1,5],
            [3,2,1,0,5,4],
            [3,1,5,0,4,2],
            [3,5,4,0,2,1],
            
            [4,0,2,1,3,5],
            [4,2,3,1,5,0],
            [4,3,5,1,0,2],
            [4,5,0,1,2,3],
            
            [5,1,0,2,4,3],
            [5,3,1,2,0,4],
            [5,0,4,2,3,1],
            [5,4,3,2,1,0],
            
            
            
            ]


face_type = [0,0,0]
Cubic = LatticeParticle(name, generators, edges, neighbors, basis, face_type, rotation_axis, rotation_angle)
myLattices.append(Cubic)


name = 'BodyCenteredCubic'
generators = [(1,0,0), (0,1,0), (1/2,1/2,1/2)]

edges = [(0,0,1), (-1,0,1), (-1,-1,1), (0,-1,1)]
neighbors = [np.array(e) for e in edges]+[-np.array(e) for e in edges]
basis = [0,1,3]
face_type = [0,0,0,0]

rotations = [#(0,1,2,3), (4,5,6,7)
            [0,1,2,3,4,5,6,7],
            [3,0,1,2,7,4,5,6],
            [2,3,0,1,6,7,4,5],
            [1,2,3,0,5,6,7,4],
            
            #(0,3,5,6), (1,2,4,7)
            [3,2,4,5,7,6,0,1],
            [5,4,7,6,1,0,3,2],
            [6,7,1,0,2,3,5,4],
             
              #(0,1,7,6), (2,4,5,3)
            [1,7,4,2,5,3,0,6],
            [7,6,5,4,3,2,1,0],
            [6,0,3,5,2,4,7,1]]
# BCC = LatticeParticle(name, generators, edges, neighbors, basis,face_type, rotation_axis, rotation_angle)
# myLattices.append(BCC)
        

for lattice in myLattices:
    lattice.save_lattice_file(directory)
    for n_particles in range (1, n_particles_max+1):
        lattice.save_structure_file(directory, n_particles)
    

    
    

    
    
    