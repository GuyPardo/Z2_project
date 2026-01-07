# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 12:01:25 2025

@author: owner
"""

import qutip as qp

def tensor_op(obj, default_obj, index, N):
    obj_list = []
    for i in range(index):
        obj_list.append(default_obj)
    obj_list.append(obj)
    for i in range(index+1, N):
        obj_list.append(default_obj)
    return qp.tensor(obj_list)

def tensor_op_sim(objs, default_obj, indices, N):
    obj_list = [default_obj for _ in range(N)]
    
    for pos, ind in enumerate(indices):
        obj_list[ind] = objs[pos]
    
    return qp.tensor(obj_list)
    