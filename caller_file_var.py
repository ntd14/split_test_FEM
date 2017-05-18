
import os as os
import math as math
import csv as csv
    
wd_path = "/home/nick/git/split_test_FEM/"
os.chdir(wd_path)
mesh_list = os.listdir("./final_xml/")
num_of_samples = 1000

for ii in range(0,len(mesh_list)):
    mesh_list[ii] = ("./final_xml/"+mesh_list[ii])

samp_list = os.listdir("./samples_var/")
num_exist_samps = len(samp_list)

seed_num = num_exist_samps;
th = math.pi/5

for ii in range(num_exist_samps, num_of_samples):
    seed_num = seed_num + 1
    list_out = [[] for kk in range(len(mesh_list))]
    for jj in range(0,len(mesh_list)):
        if mesh_list[jj].startswith("./final_xml/centre_cut.xml"):
            mesh_path = wd_path+"final_xml/centre_cut.xml"
            theta = 0
            cutp = 1
            cutn = -1

        elif mesh_list[jj].startswith("./final_xml/qtest.xml"):
            mesh_path = wd_path+"final_xml/qtest.xml"
            theta = math.pi/2
            cutp = 1
            cutn = -1
            
        elif mesh_list[jj].startswith("./final_xml/cut"):
            if mesh_list[jj]=="./final_xml/cut_1.xml":
                mesh_path = wd_path+"final_xml/cut_1.xml"
                theta = 0
                cutp = 5
                cutn = 3
            elif mesh_list[jj]=="./final_xml/cut_2.xml":
                mesh_path = wd_path+"final_xml/cut_2.xml"
                theta = 0
                cutp = 9
                cutn = 7
            elif mesh_list[jj]=="./final_xml/cut_3.xml":
                mesh_path = wd_path+"final_xml/cut_3.xml"
                theta = 0
                cutp = -3
                cutn = -5
            elif mesh_list[jj]=="./final_xml/cut_4.xml":
                mesh_path = wd_path+"final_xml/cut_4.xml"
                theta = 0
                cutp = -7
                cutn = -9
            else:
                print("cut* file with undefined properties")
                
        elif mesh_list[jj].startswith("./final_xml/rot"):
            if mesh_list[jj]=="./final_xml/rot_1.xml":
                mesh_path = wd_path+"final_xml/rot_1.xml"
                theta = -1*th
                cutp = 1
                cutn = -1  
            elif mesh_list[jj]=="./final_xml/rot_2.xml":
                mesh_path = wd_path+"final_xml/rot_2.xml"
                theta = -2*th
                cutp = 1
                cutn = -1 
            elif mesh_list[jj]=="./final_xml/rot_3.xml":
                mesh_path = wd_path+"final_xml/rot_3.xml"
                theta = -3*th
                cutp = 1
                cutn = -1   
            elif mesh_list[jj]=="./final_xml/rot_4.xml":
                mesh_path = wd_path+"final_xml/rot_4.xml"
                theta = -4*th
                cutp = 1
                cutn = -1  
            else:
                print("rot* file with undefined properties")     
                           
        else:
            print("strange file present check dir")
            
        os.system("python ./split_test_FEM_var.py "+str(mesh_path)+' '+str(theta)+' '+str(cutp)+' '+str(cutn)+' '+str(seed_num))
        with open('tfile_var.csv', 'r') as fh:
            mes = fh.read()
            fh.close()
        list_out[jj] = [str(mesh_path), str(theta), str(cutp), str(cutn), str(seed_num), str(mes)]
        
    with open('./samples_var/sample_'+str(ii)+'.csv', 'w+') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerows(list_out)


        
