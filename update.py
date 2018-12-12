import numpy as np
from biopandas.mol2 import PandasMol2
import pandas as pd




def splitter(str):
        a=str.split()
        return a
##this function permitt to separate all the char in the string removing tab in t
##here, this way i can obtain a new string with alla the terms i need

def readmol2():
	fn=open(input(),"r")
	mol2=[x for x in fn.readlines()]
	fn.close()
	return (mol2)

##this is the function that take in input the path of the file we want to read and store it in fn-variable.
##then i use the readlines() module to read separately all the line of this file one by one, and i do this x times= number of file's lines.at the end i obtain a mol2 array where each line is stored as a string.
#bond_row.append()


ls=[]   #ls is the array where i'm going to store all the information splitted.then i have to convert them by columns

mol=readmol2()
for i in range(len(mol)):
	b=splitter(mol[i]) #b i s a general name var, for each line i create ls
	ls.append(b)

## now i have to look for the section we need so i have to determine the index position where it start. 

ind1=mol.index('@<TRIPOS>ATOM\n')
ind2=mol.index('@<TRIPOS>BOND\n')
ind3=mol.index('@<TRIPOS>HEADTAIL\n')
ind4=mol.index('@<TRIPOS>RESIDUECONNECT\n')
ind5=mol.index('@<TRIPOS>SUBSTRUCTURE\n')
## i do this for all sections (2 of them are still missed),then I print the interested section in a separated array(name_section)
##sections are limitated by strings with @TRIPOS,so i just have to print the string immediatly after the beginning af the section till the beginning of the next one(ind1+1:ind2), with ind2 not considered. 



atom_section=ls[ind1+1:ind2]
bond_section=ls[ind2+1:ind5]
headtail_section=ls[ind3+1:ind4]
residueconnect_section=ls[ind4+1:]


##now i want to create different array that can store all the same type of data so i create different arrays for atom_id , atom_name ecc..
## i know that for each array i have in ls[], in position 0 I have atom_id, in 1 i have the name and so on. so it is easy to separate and put them in different arrays, that im going to use to create df.
##ATOM SECTION
atom_id=[]
atom_name=[]
x_coor=[]
y_coor=[]
z_coor=[]
atom_type=[]
subst_id=[]
subst_name=[]
charge=[]



for i in range(len(atom_section)):
	
	atom_id.append(int(atom_section[i][0]))
	atom_name.append(atom_section[i][1])
	x_coor.append(float(atom_section[i][2]))
	y_coor.append(float(atom_section[i][3]))
	z_coor.append(float(atom_section[i][4]))
	atom_type.append(atom_section[i][5])
	subst_id.append(int(atom_section[i][6]))
	subst_name.append(atom_section[i][7])
	charge.append(float(atom_section[i][8]))



##BOND_SECTION
bond_id=[]
origin_atom_id=[]
target_atom_id=[]
bond_type=[]

for i in range(len(bond_section)):
	
	bond_id.append(int(bond_section[i][0]))
	origin_atom_id.append(int(bond_section[i][1]))
	target_atom_id.append(int(bond_section[i][2]))
	bond_type.append(bond_section[i][3])

##HEADTAIL_SECTION

Head=[]
Tail=[]
for i in range(len(atom_section)):
	if i in range(len(headtail_section)):
		a=headtail_section[0][i]
		b=headtail_section[1][i]
		Head.append(a)
		Tail.append(b)
	else:
		Head.append(None)
		Tail.append(None)
##BONDING ROW
######################################################
bond_row_head=[]
bond_row_tail=[]
add_row=[None for i in range(len(atom_id))]

ind=True
ind0=True
for i in range(len(atom_id)):
	if atom_name[i]==Head[0]:
		ind0=atom_id[i]

	elif atom_name[i]==Tail[0]:
		ind=atom_id[i]

	else:
		ind=False
		ind0=False


for j in range(len(bond_id)):

	if origin_atom_id[j]==ind0:
		bond_row_head.append(target_atom_id[j])

	if target_atom_id[j]==ind0:
		bond_row_head.append(origin_atom_id[j])

	if origin_atom_id[j]==ind:
		bond_row_tail.append(target_atom_id[j])

	if target_atom_id[j]==ind:
		bond_row_tail.append(origin_atom_id[j])
	else:
		continue
if ind0==True:
	add_row[ind0]=bond_row_head
if ind==True:
	add_row[ind]=bond_row_tail


################################################
list_keys=['atom_id','atom_name','x_coor','y_coor','z_coor','atom_type','subst_id','subst_name','charge','head','tail','add row']
list_values=[atom_id,atom_name,x_coor,y_coor,z_coor,atom_type,subst_id,subst_name,charge,Head,Tail,add_row]

zipped=list(zip(list_keys,list_values))     #i found a better way to create dataframe, this way it looks more ordered
data=dict(zipped)
nodf=pd.DataFrame(data)
df=nodf.set_index('atom_id')     # i have changed the first atomcolumn with atom_id column because they were redundant
print(df)

list_keys_bond=['Bond_id','Origin_atom_id','Target_atom_id','Bond_type']
list_values_bond=[bond_id,origin_atom_id,target_atom_id,bond_type]
zipped_bond=list(zip(list_keys_bond,list_values_bond))
data_Bond=dict(zipped_bond)            #same thing for bond_section.
nodf_bond=pd.DataFrame(data_Bond)	##I know it is not necessary to print bond section 
df_Bond=nodf_bond.set_index('Bond_id')##but i  did it for checK
print(df_Bond)


#print(pmol.df.head(10))
#print(pmol.df.describe())
#print(pmol.df.info())
#print(pmol.df['atom_type'] == 'H')
#print(pmol.df[pmol.df['atom_type'] == 'H'])
#print(pmol.df[pmol.df['atom_name']=='C1']['atom_id'])
#aname=pmol.df['atom_name']
#xyz = pmol.df[['x', 'y', 'z']].values

##CLASSES

class CapFragment():
	def __init__(self,df):
		self.df = df
	def change_first_name(self, new_name):
		self.names[0] = new_name
	def shift(self, vector):
		self.coords = self.coords + vector                         



class LinkFragment():
	def __init__(self,df):
		self.df = df
	def change_first_name(self, new_name):
		self.names[0] = new_name
	def shift(self, vector):
		self.coords = self.coords + vector
