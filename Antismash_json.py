#coding=utf-8
import json
import os
import pandas as pd
import numpy as np
import re
import glob
from pathlib import Path
import sys
import subprocess


path=sys.argv[1]
outpath=sys.argv[2]


def extract_missing_files(zip_path, extract_to):
    # 获取目标文件夹中已存在的文件
    command = f'unzip -u {zip_path} -d {extract_to}'
    subprocess.run(command, shell=True, check=True)

# def my_generator(path):
#     print(path)
#     for subfolder in os.listdir(path):
#         subfolder_path = os.path.join(path, subfolder)
#         if os.path.isdir(subfolder_path):
#             json_files = glob.glob(os.path.join(subfolder_path, '*.json'))
#             if len(json_files) < 1:
#                 print(subfolder_path)
#                 zip_path = glob.glob(os.path.join(subfolder_path, '*.zip'))
#                 print(zip_path)
#                 if os.path.exists(zip_path[0]):
#                     extract_missing_files(zip_path[0], subfolder_path)
#                     json_files = glob.glob(os.path.join(subfolder_path, '*.json'))
#             for json_file in json_files:
#                 yield json_file

# 提取缺失文件的函数
# def extract_missing_files(zip_path, extract_to):
#     # 使用 zipfile 解压文件
#     with zipfile.ZipFile(zip_path, 'r') as zip_ref:
#         zip_ref.extractall(extract_to)

# 改进的 my_generator 函数
def my_generator(path):
    print(f"Processing path: {path}")
    for entry in os.scandir(path):
        if entry.is_dir():  # 仅处理目录
            subfolder_path = entry.path
            json_files = glob.glob(os.path.join(subfolder_path, '*.json'))

            if len(json_files) < 1:  # 如果没有找到 .json 文件
                print(f"No JSON files in {subfolder_path}, looking for .zip files.")
                zip_path = glob.glob(os.path.join(subfolder_path, '*.zip'))
                if zip_path:  # 如果找到了 .zip 文件
                    print(f"Found zip file: {zip_path[0]}, extracting.")
                    extract_missing_files(zip_path[0], subfolder_path)
                    # 再次检查是否提取了 .json 文件
                    json_files = glob.glob(os.path.join(subfolder_path, '*.json'))
            
            # 返回每个 .json 文件
            for json_file in json_files:
                yield json_file

n=0
ALL_df = pd.DataFrame()

for file_name in my_generator(path):
    print(file_name)
    n=n+1
    r=0
    #name = file_name.strip(".json")
    name_dir = os.path.dirname(file_name)
    name='10.8.5.43:8888'+name_dir+'/index.html'
    Pos_Charge=['Arg',"His","Lys","Orn","Dab","Dap"]
    try: 
        file = open(file_name, encoding = 'utf-8')
        lines = [line.strip('\n') for line in file.readlines()]
        file.close()
        string = "".join(lines)
        root = json.loads(string)
        records = root['records']
        Areas_df = pd.DataFrame()
        
        for record in records:
            record_id = record["id"]
            #从features中提取BGC长度、domain方向信息
            features = record["features"]
            areas = record["areas"]
            modules = record['modules']
            region_data = []
            ID_data=[]
            areas_data = []
            Length=''
            i=0
            r=r+1
            start=0
            end=0
            Condensation = '0'
            CAL_Domain = 0
            if 'source' in record['annotations']:
                annotations=record['annotations']['source']
            else:
                annotations=record['description'].strip(record_id).split(',')[0]
            for feature in features:
                if feature['type'] == 'aSDomain' and feature['qualifiers']['aSDomain'][0] == 'Condensation' and 'Condensation_Starter' in feature['qualifiers']['domain_id'][0]:
                    Condensation = '1'
                if feature['type'] == 'aSDomain' and feature['qualifiers']['aSDomain'][0] == 'CAL_domain':
                    CAL_Domain = CAL_Domain + 1
                if feature['type'] == 'aSDomain' and (feature['qualifiers']['aSDomain'][0] == 'AMP-binding' or feature['qualifiers']['aSDomain'][0] == 'Thioesterase'or feature['qualifiers']['aSDomain'][0] == 'CAL_domain'):
                    ID = feature['qualifiers']['domain_id'][0]
                    if 'join' in feature['location']:
                        ID_Start = int(min(re.findall("\d+\.?\d*", feature['location'])))
                        ID_End = int(max(re.findall("\d+\.?\d*", feature['location'])))
                        if feature['location'][-3] == "+":
                            Strand = 1
                        else:
                            Strand = -1
                    else:
                        ID_Start = int(feature['location'].strip("(+)(-)").strip('[]').split(":")[0].strip('<').strip('>'))
                        ID_End = int(feature['location'].strip("(+)(-)").strip('[]').split(":")[1].strip('<').strip('>'))
                        if feature['location'][-2] == "+":
                            Strand = 1
                        else:
                            Strand = -1
                    ID_data.append([ID, Strand, ID_Start, ID_End])
            ID_df = pd.DataFrame(ID_data, columns = ["ID","Strand", "ID_Start", "ID_End"])

            for area in areas:
                i=i+1
                j=0
                Te=0
                z=0
                string=','
                name1 = f"{name}#r{r}c{i}"
                if 'NRPS' in area['products'] and len(area['products']) > 1:
                    Product = string.join(str(item) for item in area['products'])                    
                #if 'NRPS' in area['products']:
                    ID=''
                    line_add = ''
                    line_sco = ''
                    CAL_add = ''
                    CAL_Score_add = ''
                    A_Domain = ''
                    Score = ''
                    areas_ID =  record_id + "." + str(i)
                    start_ls=[]
                    end_ls=[]
                    candid_start_end=[]
                    candid_kind=[]

                    for t in area['protoclusters']:
                        if area['protoclusters'][str(t)]['product'] == 'NRPS':
                            start = area['protoclusters'][str(t)]['start']
                            end = area['protoclusters'][str(t)]['end']
                            NRP_protoclu=t
                    start_end=str(start)+'_'+str(end)
                    #根据位置判断NRPS对应的是否为single
                    # for  candid in area['candidates']:
                    #     candid_start=candid['start']
                    #     candid_end=candid['end']
                    #     candid_start_end.append(str(candid_start)+'_'+str(candid_end))
                    #     candid_kind.append(candid['kind'])
                    # if start_end not in candid_start_end:
                    #     continue
                    # elif candid_kind[candid_start_end.index(start_end)]!='single':
                    #     continue
                    ID = ID_df[(ID_df['ID_Start'] > start) & (ID_df['ID_End'] < end)]['ID'].tolist()
                    location = str(start) + '-' + str(end)
                    list=[]
                    if (len(modules['antismash.modules.nrps_pks']["region_predictions"])==0):
                        continue
                    Antismash_Prediction = modules['antismash.modules.nrps_pks']["region_predictions"][str(i)][0]['polymer']
                    knownclusterblast_ID=(name_dir+"/svg/knownclusterblast_r"+str(r)+"c"+str(i)+"_all.svg")
                    known_ID=os.popen('grep '+ " BGC.*class " + knownclusterblast_ID + '|'+ "sed " + " -n " + " 1p " + "|"+"sed "+"s/.*label=//g"+"|"+"sed "+"s/description=.*//g")
                    known_ID=known_ID.readline().strip('\n').strip('" ')
                    known_per=os.popen('grep '+ " BGC.*class " + knownclusterblast_ID + '|'+ "sed " + " -n " + " 1p " + "|"+"sed "+"s/.*\(//g"+ "|" + "sed "+"s/of.*//g")
                    known_per=known_per.readline().strip('\n')


                    for id in range(len(ID)):
                        data = ID[id]
                        if  "Thioesterase" in data:
                            Te=1
                            if  ID_df[ID_df['ID'] == data]["Strand"].iloc[0] == 1 :
                                line_add =  line_add + ','
                                line_sco =  line_sco + ','
                                CAL_add = CAL_add + ','
                                CAL_Score_add = CAL_Score_add + ','
                            if  ID_df[ID_df['ID'] == data]["Strand"].iloc[0] == -1 :
                                line_add =  ',' + line_add
                                line_sco =  ',' + line_sco
                                CAL_add = "," + CAL_add
                                CAL_Score_add = "," + CAL_Score_add
                        else:
                            j = j + 1
                            if 'nrpys' in modules['antismash.modules.nrps_pks']["domain_predictions"][data] and 'stachelhaus_matches' in modules['antismash.modules.nrps_pks']["domain_predictions"][data]['nrpys']:
                                if modules['antismash.modules.nrps_pks']["domain_predictions"][data]['nrpys']['stachelhaus_matches']!=[] and len(modules['antismash.modules.nrps_pks']["domain_predictions"][data]['nrpys']['stachelhaus_matches']) == 1 and 'substrates' in modules['antismash.modules.nrps_pks']["domain_predictions"][data]['nrpys']['stachelhaus_matches'][0]:
                                    if len(modules['antismash.modules.nrps_pks']["domain_predictions"][data]['nrpys']['stachelhaus_matches'][0]['substrates']) > 1:
                                        A_Domain=''
                                        for y in modules['antismash.modules.nrps_pks']["domain_predictions"][data]['nrpys']['stachelhaus_matches'][0]['substrates']:
                                            A_Domain_part1 = y['short']
                                            A_Domain = A_Domain + "(" + A_Domain_part1 + ")"
                                    else:
                                        A_Domain = modules['antismash.modules.nrps_pks']["domain_predictions"][data]['nrpys']['stachelhaus_matches'][0]['substrates'][0]['short']
                                    Score = str(round(modules['antismash.modules.nrps_pks']["domain_predictions"][data]['nrpys']['stachelhaus_matches'][0]['aa10_score']*100)) + "%"
                                elif modules['antismash.modules.nrps_pks']["domain_predictions"][data]['nrpys']['stachelhaus_matches']!=[] and len(modules['antismash.modules.nrps_pks']["domain_predictions"][data]['nrpys']['stachelhaus_matches']) > 1:
                                    A_Domain = ''
                                    for x in modules['antismash.modules.nrps_pks']["domain_predictions"][data]['nrpys']['stachelhaus_matches']:
                                        if len(x['substrates']) > 1:
                                            A_Domain_part2 = ''
                                            for m in x['substrates']:
                                                A_Domain_part3 = m['short']
                                                A_Domain_part2 = A_Domain_part2 + "(" + A_Domain_part3 + ")"
                                        else:
                                            A_Domain_part2 = x['substrates'][0]['short']
                                        Score_part = str(round(x['aa10_score']*100)) + "%"
                                        A_Domain = A_Domain + "[" + A_Domain_part2 + "]"
                                    Score = str(Score) + "[" + Score_part + "]"
                                else:
                                    A_Domain = '*'
                                    Score = 0
                            if'minowa_cal' in modules['antismash.modules.nrps_pks']["domain_predictions"][data]:
                                print(modules['antismash.modules.nrps_pks']["domain_predictions"][data])
                            if 'minowa_cal' in modules['antismash.modules.nrps_pks']["domain_predictions"][data] and 'predictions' in modules['antismash.modules.nrps_pks']["domain_predictions"][data]['minowa_cal'] and modules['antismash.modules.nrps_pks']["domain_predictions"][data]['minowa_cal']['predictions']!=[]:
                                CAL = modules['antismash.modules.nrps_pks']["domain_predictions"][data]['minowa_cal']['predictions'][0][0]
                                CAL_Score = modules['antismash.modules.nrps_pks']["domain_predictions"][data]['minowa_cal']['predictions'][0][1]
                                if  ID_df[ID_df['ID'] == data]["Strand"].iloc[0] == 1 :
                                    CAL_add = CAL_add + '|' + CAL
                                    CAL_Score_add = CAL_Score_add + '|' + str(CAL_Score)
                                else:
                                    CAL_add = CAL + '|' + CAL_add
                                    CAL_Score_add = str(CAL_Score) + '|' + CAL_Score_add
                                print(CAL_add)
                            if any(substring in A_Domain for substring in Pos_Charge):
                                z = z+1
                            if  ID_df[ID_df['ID'] == data]["Strand"].iloc[0] == 1 :
                                line_add =  line_add + '|' + A_Domain
                                line_sco =   line_sco + '|' + str(Score)
                            if  ID_df[ID_df['ID'] == data]["Strand"].iloc[0] == -1 :
                                line_add =  A_Domain + '|' + line_add
                                line_sco =   str(Score) + '|' + line_sco
                    areas_data.append([Product,name1, areas_ID, j, z, Condensation, Te,CAL_Domain, location, Antismash_Prediction, line_add, line_sco,CAL_add, CAL_Score_add, known_ID, known_per, annotations])
                    print(areas_data)
            areas_df = pd.DataFrame(areas_data, columns = ["Product","Name","Region_ID","A_Domain_Num","Postive_Charge_Num","Condensation_Start","Thioesterase","CAL_Domain","Location", "Antismash_Products","A_Doamin", "A_Doamin_Scores","CAL_Domain","CAL_Domain_Scores", "Similarity_ID", "Similarity", "Source"])
            Areas_df = Areas_df._append(areas_df)
    except json.decoder.JSONDecodeError as e:
        print("json.decoder.JSONDecodeError:", e)
    ALL_df =  ALL_df._append(Areas_df)
    #print(ALL_df)
    print(n)
ALL_df.to_csv(outpath, sep=',', index=False, header=True)
