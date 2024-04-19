#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 10:40:26 2024

@author: aellien
"""
import ray
@ray.remote
def f(x):
    return x * x

if __name__ == "__main__":
    
    
    ray.init()
    futures = [f.remote(i) for i in range(1000)]
    print(ray.get(futures)) # [0, 1, 4, 9]
    ray.close()
    
    
    