import os
import csv
import openpyxl
from tkinter import *
import tkinter.filedialog



def prompt():
    print('Please choose an option from the list:')
    print('1. Find the binomial distribution for a single peptide')
    print('2. Find the binomial distribution for all peptides within an experiment')
    print('Please enter 1 or 0')
    selection = input('Enter your choice: ').strip()
    while check_selection(selection) == False:
        print('\n')
        print('This is an invalid input.')
        print('Please enter 1 or 0')
        selection = input('Enter your choice: ').strip()
    return selection

def check_selection(selection):
    if selection == '1' or selection == '2':
        return True
    else:
        return False
    
def find_folder():
    Tk().withdraw()
    print("Please select a directory.")
    folder_name = tkinter.filedialog.askdirectory(initialdir = os.getcwd(), title= "Please select a directory that contains your files")    
    print("You choice is %s" % folder_name)
    return folder_name


def find_mass(peptide):
    #to_do: look up associated mass for each amino acid in the peptide, then find sum of these 
    #this association can be found within VB code ilchenko
    return totMass

'''
this is going to need more explanation. I need to know what Isotope/Isotope calcuation
is called. 
'''
def calc_abundance():
    #this can be found on the web
    #these abundances are given
    X2, X13, X14, X16, X17, X18, X32, X33, X34 = 0.00015, 0.011074,0.99635,0.99759, 0.00037, 0.00204, 0.95,0.0076,0.0422
    X1, X12, X15, = 1-X2, 1-X13, 1-X14
    Y2, Y13, Y15, Y17, Y18, Y33, Y34 = X2/X1, X13/X12, X15/X14, X17/X16, X18/X16, X33/X32, X34/X32
    D2 = X1
    D1 = 1- D2
    YD2 = D2/D1

    abundance = [Y2, Y13, Y15, Y17, Y18, Y33, Y34, YD2]
    return abundance


def find_total_elements(peptide):
    peplist = list(peptide)
    # 'key' : Sa,Ca,Na,Oa,Ha,Da(hellenstein)
    Sp,Cp,Np,Op,Hp,Dp = 0,0,0,0,0,0
    AA_dict = {
        'G': [0,2,1,1,3,2.06],
        'A': [0,3,1,1,7,4],
        'S': [0,3,1,2,5,2.61],
        'P': [0,5,1,1,7,2.59],
        'V': [0,5,1,1,9,0.42],
        'T': [0,4,1,7,0.2],
        'C': [1,3,1,1,5,1.62],
        'c': [1,5,2,2,6,1.62],
        'L': [0,6,1,1,11,0.6],
        'l': [0,6,1,1,11,1], #this may be either lowercase L or I
        'N': [0,4,2,2,6,1.89],
        'O': [0,5,2,1,10,1],
        'D': [0,4,1,3,5,1.89],
        'Q': [0,5,2,2,8,3.95],
        'K': [0,6,2,1,12,0.54],
        'E': [0,5,1,3,7,3.95],
        'M': [1,5,1,1,9,1.12],
        'm': [1,5,1,2,9,1.12],
        'H': [0,9,1,1,9,0.32],
        'F': [0,6,4,1,12,3.43],
        'R': [0,6,4,1,12,3.43],
        'Y': [0,9,1,2,9,0.42],
        'W': [0,11,2,1,10,0.08]
        }
    for i in peplist:
        Sp += AA_dict[i][0]
        Cp += AA_dict[i][1]
        Np += AA_dict[i][2]
        Op += AA_dict[i][3]
        Hp += AA_dict[i][4]
        Dp += AA_dict[i][5]
    elements = [Sp,Cp,Np,Op,Hp,Dp]    
    return elements

def binomial_distribution(peptide):
    abd = calc_abundance()
    tot_elements = find_total_elements(peptide)
    return 

def single_distribution(): #this is for debugging purposes
    peptide = input('Please enter peptide: ').strip()
    #BWE = input('Please enter BWE as a decimal: ').strip()
    binomial_distribution(peptide)
    return

def main():

    selection = prompt()
    if selection == '1':
        single_distribution()
    elif selection == '2':
        full_distribution()      
    
    
    


main()
print('Done')
