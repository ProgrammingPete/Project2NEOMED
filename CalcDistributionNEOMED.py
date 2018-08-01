import os
import csv
import openpyxl
from tkinter import *
import tkinter.filedialog



def prompt():
    print('Please choose an option from the list:')
    print('1. Find the binomial distribution for a single peptide [DEBUGGING PURPOSES]')
    print('2. Find the binomial distribution for all peptides within an experiment')
    print('Please enter 1 or 2')
    selection = input('Enter your choice: ').strip()
    while check_selection(selection) == False:
        print('\n')
        print('This is an invalid input.')
        print('Please enter 1 or 0')
        selection = input('Enter your choice: ').strip()
        
    Ne = input('Enter the number of Isotopomers you want to use, up to two: ').strip()    
    while check_selection(Ne) == False:
        print('Enter 1 or 2')
        Ne = input('Enter the number of Isotopomers you want to use, up to three: ').strip()
    BWE = input('Please enter BWE as a decimal: ').strip()    
    return selection, BWE, Ne

def check_selection(selection):
    if selection == '1' or selection == '2':
        return True
    else:
        return False
    
def find_folder():
    window = Tk()    
    print("Please select a directory.")
    folder_name = tkinter.filedialog.askdirectory(initialdir = os.getcwd(), title= "Please select a directory that contains your files")    
    print("Your choice is %s" % folder_name)
    window.withdraw()
    return folder_name

def read_file(folder_name):
    needed_data = []
    for file_name in os.listdir(folder_name): #this creates a list from the directory (list of specimens) 
        file_name_split = file_name.split('.')        
        if len(file_name_split) < 2:
            pass
        elif file_name_split[1] == 'Quant':
            
            needed_data.append(file_name_split[0])
            #extract peptide names
            file_path = os.path.sep.join([folder_name,file_name])
            with open(file_path) as file:
                hitpeptides = False
                for row in file:
                    row = row.split(',')
                    if hitpeptides == True:
                        needed_data.append(row[0])                        
                    if row[0] == 'Peptide':
                        hitpeptides = True
                        continue                       
            needed_data.append('\n')  # this is to specify end of specimen
        else:
            pass
    print('Done reading')
    return needed_data        


def find_mass(peptide):
    """
    This is not actually needed within this program,
    it is only to check if all amino acids were included into the protein 
    """
    return totMass


def calc_abundance(BWE = 0):
    """
    % abundances
    """
    BWE = float(BWE)
    X2, X13, X14, X16, X17, X18, X32, X33, X34 = 0.00015, 0.011074,0.99635,0.99759, 0.00037, 0.00204, 0.95,0.0076,0.0422
    X1, X12, X15, = 1-X2, 1-X13, 1-X14
    Y2, Y13, Y15, Y17, Y18, Y33, Y34 = X2/X1, X13/X12, X15/X14, X17/X16, X18/X16, X33/X32, X34/X32
    if(BWE ==0):
        D2 = X2        
    else:
        D2 = BWE
    D1 = 1- D2
    YD2 = D2/D1
    abundance = [Y2, Y13, Y15, Y17 , Y33,YD2, Y18, Y34]
    return abundance

def calc_totLabel(M10, M20, M30):
    """
    this calculates the total labeling given M10,M20,M30....
    1- sum(m0 to m2)/ sum(m1 to m2)

    Ne is not needed due to the fact that the excluded distributions will be 0 by default (defined in binomial_distribution() arguments)
    """
    msumD = 1 + M10+ M20 + M30
    totalLabel = 1 - ((M10+ M20+ M30)/ msumD)
    return totalLabel
    
def calc_plateau(Y0, P, Ne, elements, BWE):
    N = elements[5]
    hp = elements[0] - N  #hp defined by total H minus Dp
    hp0 = elements[0]               #save original value of hp0
    elements[0] = hp                #this is then passed into the elements  list as such
    abd = calc_abundance(BWE)
    #abd[0] is YD0 in Visual Basic, and abd[5] is YD2 
    YD0 = abd[0]
    YD2 = abd[5]
    #calculate beta
    beta = 0
    for i in range(5):
        beta += elements[i]*abd[i]
    #calculate gamma
    gamma = 0
    for i in range(5):
        gamma += elements[i]*(elements[i]-1)*((abd[i])**2)
    #calculate omega
    omega = 0
    for i in range(4):
        if(i == 1):
            continue
        else:
            omega += elements[1]*abd[1]*elements[i]*abd[i]
    #calculate sigma 
    sigma = 0
    sigma= elements[3]*abd[6] + elements[4]*abd[7] + elements[3]*abd[3]*elements[4]*abd[4] + omega
    
    #plat1 is the denominator to the equation
    if(Ne == '1'):
        plat1 = N*(YD0-YD2)
    elif(Ne == '2'):
        plat1 = N*(YD2 - YD0)*(1+beta+((N-1)/2)*(YD2-YD0))
    elif(N == '3'):
        plat1 = N*(YD2-YD0)*(1+(gamma/2)+sigma+((N-1)/2)*(YD2-YD0)*(1+beta)+(N-1)*(N-2)*(YD2**2 + YD2*YD0)/6)
    elif(N == '4'):
        pass
    
    Plat = 1-1/(1/(1-Y0)+plat1)
    elements[0] = hp0 # restore back to original Hp 
    return Plat

def find_total_elements(peptide):
    """
        'key' : Sa,Ca,Na,Oa,Ha,Da(hellenstein). this is theoretical values given the amino acid code
    """
    peplist = list(peptide) #this creates an array of the characters
    Sp,Cp,Np,Op,Hp,Dp = 0,0,0,0,0,0
    AA_dict = {
        'G': [0,2,1,1,3,2.06],
        'A': [0,3,1,1,7,4],
        'S': [0,3,1,2,5,2.61],
        'P': [0,5,1,1,7,2.59],
        'V': [0,5,1,1,9,0.42],
        'T': [0,4,1,2,7,0.2],
        'C': [1,3,1,1,5,1.62],
        'c': [1,5,2,2,6,1.62],
        'L': [0,6,1,1,11,0.6],
        'I': [0,6,1,1,11,1],
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
        try:
            Sp += AA_dict[i][0]
            Cp += AA_dict[i][1]
            Np += AA_dict[i][2]
            Op += AA_dict[i][3]
            Hp += AA_dict[i][4]
            Dp += AA_dict[i][5]
        except KeyError:
            print("Key not found, skipping:", i)
            print("Continuing")
            pass
    elements = [Hp,Cp,Np,Op,Sp, Dp]   
    return elements

def binomial_distribution(peptide,Ne, elements, BWE = 0, M10 = 0, M20 = 0, M30 = 0, M40 = 0):
    """ Comment:
        This deals with a single peptide. If there are more than one,
        this will loop. 
    """
    
    abd = calc_abundance(BWE)
    m202, m203 = 0, 0
    for i in range(6):
        M10 += abd[i]*elements[i]

    if(Ne == '1'):
        print(M10)
        return M10,M20,M30
    elif(Ne == '2'):
        for i in range(6):
            M20 += (elements[i]*(elements[i]-1)*(abd[i]**2))/2
        for i in range(6):
            m202 += abd[i]*elements[i]
            for j in range(1+i,6):
                m203 += abd[j]*elements[j]
            m202 = m202*m203
        M20 += m202 + elements[3]*abd[6]+elements[4]*abd[7]
        return M10,M20,M30
    elif(Ne == '3'):
        #TO DO ....
        for i in range(6):
            M20 += (elements[i]*(elements[i]-1)*(elements[i]-2)*(abd[i]**3))/6
        return M10,M20,M30
    
def single_peptide(peptide,Ne, BWE):
    print(peptide)
    elements = find_total_elements(peptide)                               #find total elements in a peptide
    M10, M20, M30 = binomial_distribution(peptide,Ne, elements)           #calculate distribution at start
    print("M10 at Start: ", M10)
    print("M20 at Start: ", M20)
    TL0 = calc_totLabel(M10 , M20, M30)                                   #calcualte TL at start
    M10, M20, M30 = binomial_distribution(peptide ,Ne, elements,BWE)      #calculate distribution at Plateau BWE.. elements become overwritten
    print("M10 at Plateau (BWE): ", M10)
    print("M20 at Plateau (BWE): ", M20)
    TL2 = calc_totLabel(M10, M20, M30)                                    #calculate TL at plateau
    theoryPlat = calc_plateau(TL0, TL2, Ne, elements,BWE)                 #calculate theory plateau... elements become overwritten

    #find expirement data
    #using this data, find the plateau
    return

def single_distribution(BWE,Ne): #this is for debugging purposes
    peptide = input('Please enter peptide: ').strip()    
    single_peptide(peptide, Ne, BWE)    
    return

def full_distribution(BWE, Ne):
    #use a queue for iteration 
    from collections import deque
    peptide = ''
    folder_name = find_folder()
    needed_data = read_file(folder_name)
    Queue = deque(needed_data)

    print('Finding values for:', Queue.popleft())
    while(len(Queue) != 0 ):
        peptide = Queue.popleft()
        if (len(Queue) == 1):
            continue
        elif(peptide == '\n'):
            print('Finding values for:', Queue.popleft())
            continue
        peptide.strip()
        #single_peptide(peptide, Ne, BWE)       
        
    return

def main():
    selection, BWE, Ne  = prompt()
    if selection == '1':
        single_distribution(BWE, Ne)
    elif selection == '2':
        full_distribution(BWE, Ne)      
    


main()
print('Done')
