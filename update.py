import math
import numpy as np
#from biopandas.mol2 import PandasMol2
import pandas as pd
import sys
import os
import copy
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

#Thiol(Fragments.split(',')[0])

def splitter(str):
        a=str.split()
        return a
##this function permitt to separate all the char in the string removing tab in t
##here, this way i can obtain a new string with alla the terms i need
def readmol2(filename):
	fn=open(filename,"r")
	mol2=[x for x in fn.readlines()]
	fn.close()
	return (mol2)

	##this is the function that take in input the path of the file we want to read and store it in fn-variable.
	##then i use the readlines() module to read separately all the line of this file one by one, and i do this x times= number of file's lines.at the end i obtain a mol2 array where each line is stored as a strin.
	#bond_row.append()

def reader(filename):
	ls=[]   #ls is the array where i'm going to store all the information splitted.then i have to convert them by columns

	mol=readmol2(filename)
	for i in range(len(mol)):
		b=splitter(mol[i]) #b i s a general name var, for each line i create ls
		ls.append(b)

	## now i have to look for the section we need so i have to determine the index position where it start.

	ind1=mol.index('@<TRIPOS>ATOM\n')
	ind2=mol.index('@<TRIPOS>BOND\n')
	ind3=mol.index('@<TRIPOS>HEADTAIL\n')
	ind4=mol.index('@<TRIPOS>RESIDUECONNECT\n')
	ind5=mol.index('@<TRIPOS>SUBSTRUCTURE\n')
	## i do this for all sections (2 of them are still missed),then I prin t the interested section in a separated array(name_section)
	##sections are limitated by strings with @TRIPOS,so i just have to prin t the string immediatly after the beginning af the section till the beginning of the next one(ind1+1:ind2), with ind2 not considered.



	atom_section=ls[ind1+1:ind2]
	bond_section=ls[ind2+1:ind5]
	headtail_section=ls[ind3+1:ind4]
	residueconnect_section=ls[ind4+1:]

	np_headtail_section=np.array(headtail_section)
	np_atom_section=np.array(atom_section)
	##now i want to create different list that can store all the same type of data so i create different arrays for atom_id , atom_name ecc..
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
	#i know that the first column of each row is the atom name so faor all the lines
	#of atom section i have to append to atom_id list the first element
	#i do the same things for all the other columns


	atom_id= np.array(np_atom_section[:,0],dtype='int')
	atom_name= np.array(np_atom_section[:,1],dtype='str')
	x_coor= np.array(np_atom_section[:,2],dtype='float')
	y_coor= np.array(np_atom_section[:,3],dtype='float')
	z_coor= np.array(np_atom_section[:,4],dtype='float')
	atom_type= np.array(np_atom_section[:,5],dtype='str')
	subst_id= np.array(np_atom_section[:,6],dtype='int')
	subst_name= np.array(np_atom_section[:,7],dtype='str')
	charge= np.array(np_atom_section[:,8],dtype='float')


	##BOND_SECTION
	np_atom_name=np.array(atom_name)
	np_atom_id=np.array(atom_id)

	bond_id=[]
	origin_atom_id=[]
	target_atom_id=[]
	bond_type=[]

	for i in range(len(bond_section)):

		bond_id.append(bond_section[i][0])
		origin_atom_id.append(bond_section[i][1])
		target_atom_id.append(bond_section[i][2])
		bond_type.append(bond_section[i][3])

	##HEADTAIL_SECTION
	np_origin_atom_id=np.array(origin_atom_id, dtype = 'int')
	np_target_atom_id=np.array(target_atom_id, dtype = 'int')
	np_bond_id=np.array(bond_id, dtype='int')

	Head=[]
	Tail=[]
	np_head=np.zeros(len(atom_id),dtype='bool')
	np_tail=np.zeros(len(atom_id),dtype='bool')
	#if np_head()

	headtail_section

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
	#########################################################################
	bond_row_head=[]
	bond_row_tail=[]
	np_row=np.zeros(len(atom_id),dtype='bool')
	#np_row.fill(None)


	#add_row=[[]]*len(atom_id)
	#add_row=list(np_row)
	#ind=True

	#idnd0=True
	#for i in range(len(atom_id)):
	#	if atom_name[i]==Head[0]
	#		ind0=atom_id[i]
	#	else:
	#		ind0==False
	#	if atom_name[i]==Tail[0]:
	#		ind=atom_id[i]
	#	else:
	#		ind==False

	#ind0=np_atom_id[np.where(np_atom_name==Head[0])[0]]

	theres_head = np_atom_name == Head[0]
	if np.any(theres_head):
		ind0=int(np_atom_id[np.where(theres_head)[0]])

		s=np.where(ind0 == np_origin_atom_id)[0]
		np_bond_row_head=(np_target_atom_id[s])
		s=np.where(ind0 == np_target_atom_id)[0]
		np_bond_row_head=np.append(np_bond_row_head,np_origin_atom_id[s])

	theres_tail = np_atom_name == Tail[0]

	if np.any(theres_tail):
		ind=int(np_atom_id[np.where(theres_tail)[0]])
		s=np.where(ind == np_origin_atom_id)[0]
		np_bond_row_tail=(np_target_atom_id[s])
		s=np.where(ind == np_target_atom_id)[0]
		np_bond_row_tail=np.append(np_bond_row_tail,np_origin_atom_id[s])
	else:
		ind=0
		np_bond_row_tail=[]

	Head = np.zeros(len(np_atom_id), dtype = 'bool')
	Tail = np.zeros(len(np_atom_id), dtype = 'bool')
	Head[theres_head] = True
	Tail[theres_tail] = True

	#a_name= headtail_section[0]

	#in this section i want to add an other column which is seet to None for all the lines exept for the tail and the head column. in that rows i have to put a list of atom_id bonded to them.
	#Indeed i created 2 list: one for atoms bonded to the tail and one for the head

	#schedule:fragments can have head or tail or both, then i have to check in what case we are. the procedure is the same for both situations.
	#first i have to look at the first row of Tail where there is the atom_name of the tail part. thenm ii have to go to the Bond_section because in the second and in the third columns we have the origin atom_id and which atom is bonded to.
	#so i need the atom_id and i take it in the first part with ind and ind0
	#then i look which atom, the tail atom is bonded to and append all in the lists created.
	##############################################################
	list_keys=['atom_id','atom_name','x_coor','y_coor','z_coor','atom_type','subst_id','subst_name','charge']
	list_values=[atom_id,atom_name,x_coor,y_coor,z_coor,atom_type,subst_id,subst_name,charge]
	list_xyz=['x_coor','y_coor','z_coor']
	list_coor=[x_coor,y_coor,z_coor]
	#list_ht=['a_name','subst_id']
	#list_ht_values=[headtail_section[0],headtail_section[1]]


	zipp=list(zip(list_xyz,list_coor))
	zipped=list(zip(list_keys,list_values))

	#ht=dict(zi)
	xyz_coor=dict(zipp)
	data=dict(zipped)

	df=pd.DataFrame(data)
	df=df.set_index('atom_id')     # i have changed the first atomcolumn with atom_id column because they were redundant

	list_keys_bond=['Bond_id','Origin_atom_id','Target_atom_id','Bond_type']
	list_values_bond=[np_bond_id,np_origin_atom_id,np_target_atom_id,bond_type]

	zipped_bond=list(zip(list_keys_bond,list_values_bond))
	data_Bond=dict(zipped_bond)

	#nodf_bond=pd.DataFrame(data_Bond)
	nodf_Bond=pd.DataFrame(data_Bond)
	df_Bond=nodf_Bond.set_index('Bond_id')

	#df_ht= pd.DataFrame(headtail_section)
	#df_ht_new= df_ht.set_index('0')

	if ind!=0:
		type= 'Linking Fragment'
	if ind==0:
		type= 'Capping Fragment'




	return df,ind0,np_bond_row_head,ind,np_bond_row_tail,type,df_Bond

'''i have to return all that values that are the ones i need to represent all the attirbutes of my fragment'''

def center(df):
	Baricentro=np.array(np.mean(df[['x_coor','y_coor','z_coor']].values,axis=0),dtype='float')
	#i made the baricentro to shift the fragment in the origin, without axis=0 i have average of all the terms but i need a mean of x,y,z separately)

	xyz_matrix=df[['x_coor','y_coor','z_coor']].values

	subtr=np.subtract(xyz_matrix,Baricentro)
	df[['x_coor','y_coor','z_coor']]=subtr
	#subtr is the df with the coor shifted in the origine
	return df

def rotate(df):

	pca=PCA(n_components=1)
	pca.fit(df[['x_coor','y_coor','z_coor']])
	teta = math.asin(pca.components_[0][-1])
	phi = -math.acos(pca.components_[0][0]/(math.cos(teta)))
	if pca.components_[0][1] < 0:
		phi = -phi

	rotation_matrix_z=np.array([[math.cos(phi),-math.sin(phi),0],[math.sin(phi), math.cos(phi), 0],[0,0,1]])
	rotation_matrix_y=np.array([[math.cos(teta),0,math.sin(teta)],[0,1,0],[-math.sin(teta),0,math.cos(teta)]])

	rotation_yz= rotation_matrix_y.dot(rotation_matrix_z)
	product=rotation_yz.dot(df[['x_coor','y_coor','z_coor']].T).T

	df[['x_coor','y_coor','z_coor']]=product
	#product is the df with the xyz coor that has been rotated in the x direction
	return df

def orientation(df, ind, ind0):

	if ind!= 0 and df.loc[ind0].x_coor <df.loc[ind].x_coor:
		#print("Bandersnatch")
		phi=math.pi
		rotation_matrix_z=np.array([[math.cos(phi),-math.sin(phi),0],[math.sin(phi), math.cos(phi), 0],[0,0,1]])
		rotation= rotation_matrix_z.dot(df[['x_coor','y_coor','z_coor']].T)
		df[['x_coor','y_coor','z_coor']]=rotation.T

	if ind==0 and df.loc[ind0].x_coor <0:
		#print("Bandersnatch")
		phi=math.pi
		rotation_matrix_z=np.array([[math.cos(phi),-math.sin(phi),0],[math.sin(phi), math.cos(phi), 0],[0,0,1]])
		rotation= rotation_matrix_z.dot(df[['x_coor','y_coor','z_coor']].T)
		df[['x_coor','y_coor','z_coor']]=rotation.T
	return df
		#this method check if the head is in the right side, that is the negative one.
		#if the head is positive we neew to rotate the fragment with an angle of pi

def tetrahedron(a):

	v = a.df.loc[a.np_bond_row_tail][['x_coor','y_coor','z_coor']].values
	a1=np.linalg.norm(v,axis=1)
	v_n=np.divide(v,a1)

	v_t=a.df.loc[a.atom_id_tail][['x_coor','y_coor','z_coor']].values
	w=np.array(v-v_t, dtype='float')

	w1=np.linalg.norm(w,axis=1)

	w=np.divide(w,w1)

	B=np.array(np.mean(w, axis=0),dtype='float')

	B_n=np.array(np.linalg.norm(B), dtype='float')
	p= v_t-(1.54*B)/B_n

	return p

class Fragment:
	def __init__( self,filename): #i have to initialize the class so i need self and the filename,than all what i read in the file from the function reader is stored in the attributes and is referred to the OBJECT Fragment
		self.df, self.atom_id_head, self.np_bond_row_head,self.atom_id_tail,self.np_bond_row_tail,self.type,self.df_Bond=reader(filename)
		self.df = center(self.df)
		#show_pca(self.df)
		self.df= rotate(self.df)
		self.df= orientation(self.df,self.atom_id_tail, self.atom_id_head)
		#show_pca(self.df)

	def shift(self, vector):
		self.coords = self.coords - vector


class Thiol:
	def __init__(self,frag):
		self.df= frag.df
		self.df_Bond= frag.df_Bond
		self.atom_id_head, self.np_bond_row_head,self.atom_id_tail,self.np_bond_row_tail=frag.atom_id_head,frag.np_bond_row_head, frag.atom_id_tail,frag.np_bond_row_tail
	def polymerize(self,frag2):
		if frag2.type== 'Linking Fragment':
			new_tail_id=len(self.df)+frag2.atom_id_tail
			new_np_bond_row_tail= len(self.df)+frag2.np_bond_row_tail

		else:
			new_tail_id=0
			new_np_bond_row_tail=[]

		self.df,self.df_Bond= function(self,frag2)
		self.atom_id_tail= new_tail_id
		self.np_bond_row_tail=new_np_bond_row_tail
		self.df = orientation(self.df, self.atom_id_tail, self.atom_id_head)
		#show_pca(self.df)

		#self.atom_id_tail= len(self.df)+frag2.atom_id_tail
		#self.np_bond_row_tail= len(self.df)+frag2.np_bond_row_tail




#we need to decide what is the head and what is the tail
#we have gold core and the S of thiol
#we always have a head
#Fragments can have different xyz coor so the problem is that i want the same referring system
#the idea is to put for each fragment the baricentro in (x,y,z)=(0,0,0) and the axis of the molecule on the x direction

#f=Fragment('/Users/matteo/Desktop/LIG2.mol2')
#f is the Fragment created from the file indicated. The file must be a string so be careful


def show_pca(df):
	xyz = df[['x_coor', 'y_coor', 'z_coor']].values
	pca = PCA(n_components = 1)
	pca.fit(xyz)
	p = pca.components_[0]
	#print(p)
	space = np.linspace(0,10,500)
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.set_xlim((-9, 9))
	ax.set_ylim((-9, 9))
	ax.set_zlim((-9, 9))
	ax.set_xlabel("X", fontsize = 22)
	ax.set_ylabel("Y", fontsize = 22)
	ax.set_zlabel("Z", fontsize = 22)
	ax.scatter(df.x_coor, df.y_coor, df.z_coor, c='b')
	ax.plot(space*p[0], space*p[1], space*p[2], c='r')
	#plt.show()

def function(a,b):
	b_copy = [b].copy()
	b_copy = b_copy[0]
	 #i have to made a copy of the b df because otherwise this funcion would modify the real df coming from the reader function
	p=tetrahedron(a)#p is the point of the tetrahedron of fragment a starting from the tail

	#v_head_b= #is the position of the head of the b fragment

	w=p-b.df.loc[b.atom_id_head][['x_coor','y_coor','z_coor']].values # to shift the head of b in the tetrahedron point i need w
	#v_b= #but i have to shift all the fragment not only the head
	b_copy.df[['x_coor','y_coor','z_coor']]=w + b.df[['x_coor','y_coor','z_coor']].values

	#b_copy.df[['subst_id']]=a.df[['subst_id'][-1]]+1

	#last_row[['Origin_atom_id']].values[0][0]=a.df_Bond.loc[a.atom_id_tail]
	b_copy.df.index=b_copy.df.index+a.df.index[-1]
	b_copy.df_Bond.index=b_copy.df_Bond.index+a.df_Bond.index[-1]

	new_origin_atom_id=b.df_Bond[['Origin_atom_id']]+len(a.df)
	new_target_atom_id=len(a.df)+b.df_Bond[['Target_atom_id']]

	b_copy.df_Bond[['Origin_atom_id']]= new_origin_atom_id
	b_copy.df_Bond[['Target_atom_id']]= new_target_atom_id

	b_copy.df[['atom_name']] = renamer(a,b_copy)

	df_new=a.df.append(b_copy.df)

	df_Bond_new= a.df_Bond.append(b_copy.df_Bond)


	list_columns=['Bond_id','Origin_atom_id','Target_atom_id','Bond_type']
	list_columnsval=[[df_Bond_new.index[-1]+1],[a.atom_id_tail],[b.atom_id_head+len(a.df)],[1]]
	zi=list(zip(list_columns,list_columnsval))
	new_row=dict(zi)
	last_row=pd.DataFrame(new_row)
	last_row=last_row.set_index("Bond_id")
	df_Bond_new=df_Bond_new.append(pd.DataFrame(last_row))

	df_new= center(df_new)
	#show_pca(df_new)
	df_new= rotate(df_new)

	print(df_new)

	return df_new,df_Bond_new

def renamer(a,b_copy):
	q=[]
	t=[]
	count_C=0
	count_h=0
	count_ox=0
	list_C=list(a.df.atom_name.values)
	list_b_copy=list(b_copy.df.atom_name)

	for i in range(len(list_C)):
		y= list(list_C[i])
		t.append(y)
		if t[i][0]=='C' or t[i][0]=='c':
			count_C+=1
		if t[i][0]=='H' or t[i][0]=='h':
			count_h+=1
		if t[i][0]=='O' or t[i][0]=='o':
			count_ox+=1

	for i in range(len(b_copy.df.atom_name)):
		k= list(list_b_copy[i])
		t.append(k)
		if t[i][0]=='C' or t[i][0]=='c':
			t[i][1]= count_C+1
			count_C+=1
		if t[i][0]=='H' or t[i][0]=='h':
			t[i][1]= count_h+1
			count_h+=1
		if t[i][0]=='O' or t[i][0]=='o':
			t[i][1]= count_ox+1
			count_ox+=1
		q.append(t[i][0]+str(t[i][1]))


	return q





Fragments=list(sys.argv[1].split(','))
times=np.array(sys.argv[2].split(','),dtype='int')
gna=os.path.realpath(Fragments[0])
m=Fragment(gna)
t=Thiol(m)
times[0] = times[0] - 1

if sum(times)>5:
	print('\n Error: you have insered too many fragments!\n')
if times[-1]!= 1 and len(times)!=1:
	print('\n You can\'t polymerize more than 1 capping fragment!\n')
else:
	for i in range(len(Fragments)):
		gna=os.path.realpath(Fragments[i])

		for j in range(times[i]):
			m=Fragment(gna)
			t.polymerize(m)

#pol= reader(Fragments.split(',')[0])
#def rename(df_new):

def writer(thiol):
	fp= open('Thiol_prova.mol2','w')
	fp.write("@<TRIPOS>MOLECULE \n")
	fp.write("THIOL\n")
	fp.write("{0:6}{1:6}{2:6}\n".format(len(t.df),len(t.df_Bond),t.df[['subst_id']].loc[t.atom_id_head].values[0]))
	fp.write("SMALL\n")
	fp.write("USER_CHARGES\n")
	fp.write("@<TRIPOS>ATOM \n")

	for ndx,atom in t.df.iterrows():
		fp.write("{0:>4} {1:<4} {2:>10.4f} {3:>10.4f} {4:>10.4f} {5:>2} {6:>3} {7:>5} {8:>7.4f}\n".format( ndx, atom[['atom_name']].values[0],
		float(atom[['x_coor']].values[0]),
		float(atom[['y_coor']].values[0]),
		float(atom[['z_coor']].values[0]),
		atom[['atom_type']].values[0],
		atom[['subst_id']].values[0],
		atom[['subst_name']].values[0],
		float(atom[['charge']].values[0])))
	fp.write("@<TRIPOS>BONDS\n")

	for ndx,bond in t.df_Bond.iterrows():
		fp.write("{0:>5} {1:>6} {2:>6} {3:>2}\n".format(ndx,bond[['Origin_atom_id']].values[0],bond[['Target_atom_id']].values[0],bond[['Bond_type']].values[0]))
	fp.write("@<TRIPOS>SUBSTRUCTURE\n")
	fp.write("{0:>4} {1:>4} {2:>4} \n".format(t.df[['subst_id']].loc[t.atom_id_head].values[0],t.df[['subst_name']].loc[t.atom_id_head].values[0],1))

	fp.write("@<TRIPOS>HEADTAIL\n")
	#fp.write(t.df[['atom_name']].loc[t.atom_id_head].values[0],t.df[['subst_id']].loc[t.atom_id_head].values[0],"\n")
	fp.write("{0:1} {1:1} \n".format(t.df[['atom_name']].loc[t.atom_id_head].values[0],t.df[['subst_id']].loc[t.atom_id_head].values[0]))
	fp.write("{0:1} {1:1} \n".format(0,0))
	fp.write("@<TRIPOS>RESIDUECONNECT\n")
	fp.write("{0:1} {1:1} 0 0 0 0 ".format(t.df[['atom_name']].loc[t.atom_id_head].values[0],t.df[['subst_id']].loc[t.atom_id_head].values[0],0,0,0,0,0))
	fp.close()

file_text=writer(t)

'''
i dont need a linking class because i put all in one class called Fragment:
here we have head and Tail or just a head or just a Tail. ind and ind0 doesn't exist always and this is why i had to initialize that values to ind=0
the same for the list that represent all the atoms_id bonded to the tail/head group. i initialized it like an empty list=[]
'''
