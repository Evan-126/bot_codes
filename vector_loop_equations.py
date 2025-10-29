## this programs primary purpose is so i never have to write out the vector loop equations again
## it defines a function that takes in the lengths of the 4 bars and theta2 (angle between ground and input link)
## and returns theta3 and theta4 (angles between ground and coupler and output link respectively)
## additionally, it currently plots mechanical advantage vs theta2 for a given 4 bar linkage

import math as m
import matplotlib.pyplot as plt

def theta_3_and_4(l1, l2, l3, l4, theta2, state="closed"):
    K1 = l1 / l2
    K2 = l1 / l4
    K3 = (l2**2 - l3**2 + l4**2 + l1**2) / (2 * l2 * l4)
    K4 = l1 / l3
    K5 = (l4**2 - l1**2 - l2**2 - l3**2) / (2 * l2 * l3)

    q = m.cos(theta2)
    p = m.sin(theta2)

    # Match slides notation exactly
    A = q - K1 - K2*q + K3
    B = -2*p
    C = K1 - (K2+1)*q + K3
    D = q - K1 + K4*q + K5
    E = -2*p
    F = K1 + (K4 - 1)*q + K5


    # print(f"A={A}, B={B}, C={C}, D={D}, E={E}, F={F}") # For debugging
    # print(f"B^2 - 4AC = {B**2 - 4*A*C}")
    # print(f"E^2 - 4DF = {E**2 - 4*D*F}")
    # Avoid domain errors in sqrt for impossible configurations
    try:
        if state == "open":
            theta4 = 2 * m.atan2((-B + m.sqrt(B**2 - 4*A*C)), 2*A)
            theta3 = 2 * m.atan2((-E + m.sqrt(E**2 - 4*D*F)), 2*D)
        else:
            theta4 = 2 * m.atan2((-B - m.sqrt(B**2 - 4*A*C)), 2*A)
            theta3 = 2 * m.atan2((-E - m.sqrt(E**2 - 4*D*F)), 2*D)
    except ValueError as err:
        print("Math domain error occurred (check physical feasibility of linkage and inputs).")
        return None, None

    return theta3, theta4
